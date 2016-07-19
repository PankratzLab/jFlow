package org.genvisis.cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.CrossValidation;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.RegressionModel;
import org.genvisis.stats.StatsCrossTabs;
import org.genvisis.stats.Stepwise;
import org.genvisis.stats.LeastSquares.LS_TYPE;
import org.genvisis.stats.StatsCrossTabs.STAT_TYPE;
import org.genvisis.stats.StatsCrossTabs.StatsCrossTabRank;
import org.genvisis.stats.StatsCrossTabs.VALUE_TYPE;
import org.genvisis.stats.Stepwise.StepWiseSummary;

/**
 * <p>
 * Class to compute residuals wrt PCs and summarize median values, currently aimed at Mitochondrial copy number, but is extensible to other data
 * 
 *
 */
public class PrincipalComponentsResiduals implements Cloneable, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final int NUM_PC_SVD_OVERIDE = 160;
	private static final String[] MT_REPORT = { "DNA", "FID", "IID", "Sex", "median_MT_LRR_raw", "median_MT_LRR_PC_residuals", "median_MT_LRR_PC_residuals_inverseTransformed" };
	public static final String[] MT_REPORT_EXT = { ".report.txt" };
	public static final String[] MT_REPORT_MARKERS_USED = { ".MitoMarkers.MarkersUsed.txt", ".MitoMarkers.RawValues.txt" };
	private static final String[] MT_RESIDUAL_CROSS_VALIDATED_REPORT = { "Time Completed", "Time to complete(seconds)", "PC", "Cross-validation Average SSerr", "Cross-validation Average R-squared", "Average Standard Error of Betas", "Full model R-squared", "Full model SSerr" };

	private String markersToAssessFile, output, residOutput, pcFile;
	private String[] markersToAssess, samplesToReport, allProjSamples, pcTitles;
	private double[] assesmentData, residuals, invTResiduals;
	private double[][] fullData, pcBasis;
	private byte[][] abGenotypesAfterFilters;
	private float gcThreshold;
	/**
	 * if true, does not do a sample lookup for loading
	 */
	protected boolean useIID = false;
	protected Logger log;
	protected Project proj;
	protected int numComponents, totalNumComponents;
	private int numSamples;// numComponents are loaded, totalNumComponents are present in the pc file
	protected Hashtable<String, Integer> samplesInPc;// stores the index of a sample in the pc file
	protected boolean[] samplesToUse;// corresponds to the samples in the PC
	private boolean printFull, homozygousOnly, recomputeLRR, sortedByProject;
	private GcAdjustorParameters params;

	/**
	 * @param proj
	 *            current project
	 * @param pcFile
	 *            pcs to compute residuals from
	 * @param markersToAssessFile
	 *            markers to compute median values for, and regress out pcs
	 * @param numComponents
	 *            number of components to regress out
	 * @param printFull
	 *            print all data associated with the markers to assess.
	 * @param gcThreshold
	 *            threshold to use a marker
	 * @param homozygousOnly
	 *            only use homozygous calls
	 * @param recomputeLRR
	 *            recompute Log R Ratios for median markers
	 * @param output
	 * @param log
	 */
	public PrincipalComponentsResiduals(Project proj, String pcFile, String markersToAssessFile, int numComponents, boolean printFull, float gcThreshold, boolean homozygousOnly, boolean recomputeLRR, String output) {
		super();
		this.numComponents = numComponents;
		this.proj = proj;
		this.log = proj.getLog();
		this.numSamples = 0;
		this.samplesInPc = new Hashtable<String, Integer>();
		this.markersToAssessFile = markersToAssessFile;
		this.printFull = printFull;
		this.gcThreshold = gcThreshold;
		this.homozygousOnly = homozygousOnly;
		this.recomputeLRR = recomputeLRR;
		this.output = output;
		this.pcFile = pcFile;
		this.useIID = false;
		loadPcFile(pcFile, useIID);
		parseSamplesToUse();
		this.sortedByProject = determineSortedByProject();
	}

	public GcAdjustorParameters getParams() {
		return params;
	}

	public void setParams(GcAdjustorParameters params) {
		this.params = params;
	}

	public String[] getPcTitles() {
		return pcTitles;
	}

	/**
	 * Compute the median data in project order, and then set back to pcfile order, just in case
	 */
	public void computeAssessmentDataMedians() {
		this.markersToAssess = PrincipalComponentsCompute.sortByProjectMarkers(proj, HashVec.loadFileToStringArray(proj.PROJECT_DIRECTORY.getValue() + markersToAssessFile, false, new int[] { 0 }, true));
		getData();
		double[] projectOrderMedians = getLRRMedian();
		if (printFull) {
			printFull();
		}
		setAssesmentDataSortByPCs(projectOrderMedians);
	}

	public void setHomozygousOnly(boolean homozygousOnly) {
		this.homozygousOnly = homozygousOnly;
	}

	/**
	 * Regress out the Pcs from the assessment data, return the model's R-squared value for checking
	 */
	public double computeResiduals() {
		// TODO, add svdRegression option
		// RegressionModel model = (RegressionModel) new LeastSquares(assesmentData, prepPcs(pcBasis));
		// LS_TYPE lType = LS_TYPE.REGULAR;
		// if (numComponents > NUM_PC_SVD_OVERIDE) {
		// log.reportTimeInfo("Number of components " + numComponents + " greater than " + NUM_PC_SVD_OVERIDE + ", switching to " + LS_TYPE.QR_DECOMP + " decomp regression");
		// lType = LS_TYPE.QR_DECOMP;
		// }
		RegressionModel model = (RegressionModel) new LeastSquares(assesmentData, prepPcs(pcBasis), null, false, true, LS_TYPE.REGULAR);// auto switch in reg model
		double R2 = Double.NaN;
		if (!model.analysisFailed()) {
			this.residuals = model.getResiduals();
			R2 = model.getRsquare();
			log.report("" + model.getRsquare());
		} else {
			log.reportError("Error - the regression model has failed and residuals could not be computed");
			this.residuals = null;
		}
		return R2;
	}

	/**
	 * Inverse normalize the residuals
	 */
	public void computeInverseTransformResiduals() {
		if (residuals == null) {
			log.reportError("Error - must compute residuals first");
			System.exit(1);
		} else {
			this.invTResiduals = Array.inverseNormalize(residuals);
		}
	}

	public void setResidOutput(String residOutput) {
		this.residOutput = residOutput;
	}

	public void setResiduals(double[] residuals) {
		this.residuals = residuals;
	}

	/**
	 * Summarize the medians, residuals from the PCs, and inverse transformed values in the same order as the input pc file
	 */
	public void summarize(String output) {
		this.residOutput = ext.rootOf(output) + MT_REPORT_EXT[0];
		SampleData sampleData = proj.getSampleData(0, false);
		try {
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + residOutput)) {
				Files.backup(residOutput, proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + residOutput));
			writer.print(Array.toStr(MT_REPORT));
			for (int i = 0; i < numComponents; i++) {
				writer.print("\t" + PrincipalComponentsCompute.PC_STRING + (i + 1));
			}
			writer.println();
			int modelCount = 0;
			// samplesToReport is in the same sample order as the pc file
			for (int i = 0; i < samplesToReport.length; i++) {
				String[] samp = sampleData.lookup(samplesToReport[i]);
				writer.print(samp[0] + "\t" + samp[1]);
				if (sampleData.getSexForIndividual(samp[0]) == -1) {
					writer.print("\tNA");
				} else {
					writer.print("\t" + sampleData.getSexForIndividual(samp[0]));
				}
				if (Double.isNaN(assesmentData[i])) {
					writer.print("\tNA\tNA\tNA");
				} else {
					writer.print("\t" + assesmentData[i] + "\t" + (residuals == null ? "NA" : residuals[modelCount]) + "\t" + (invTResiduals == null ? "NA" : invTResiduals[modelCount]));
					modelCount++;
				}
				for (int k = 0; k < numComponents; k++) {
					writer.print("\t" + pcBasis[k][i]);
				}
				writer.println();
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + residOutput + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + residOutput + "\"");
			log.reportException(ioe);
		}
	}

	/**
	 * The residual file that the results were summarized to
	 */
	public String getResidOutput() {
		return residOutput;
	}

	/**
	 * Here we set the AssesmentData back to the same order as samples represented in the pc file
	 */
	private void setAssesmentDataSortByPCs(double[] projectOrderMedians) {
		this.assesmentData = new double[projectOrderMedians.length];
		int sampleIndex = 0;
		for (int i = 0; i < allProjSamples.length; i++) {
			// System.out.println(samplesToUse[i]);
			// System.out.println(allProjSamples[i]);
			// System.out.println(projectOrderMedians[sampleIndex]);
			// System.out.println(assesmentData.length);
			if (samplesToUse[i]) {
				assesmentData[samplesInPc.get(allProjSamples[i])] = projectOrderMedians[sampleIndex];
				sampleIndex++;
			}
		}
	}

	/**
	 * Since the samples in the principal component file may not represent the entire project, we create boolean[] samplesToUse to track the neccesary samples Also, the PC file may not represent the samples in the same order as the project so we need later step to ensure we report in the same order as the pc file.
	 */
	private void parseSamplesToUse() {
		this.allProjSamples = proj.getSampleList().getSamples();
		this.samplesToUse = new boolean[allProjSamples.length];
		int used = 0;
		for (int i = 0; i < allProjSamples.length; i++) {
			if (samplesInPc.containsKey(allProjSamples[i])) {
				samplesToUse[i] = true;
				used++;
			} else {
				samplesToUse[i] = false;
			}
		}
		if (used != samplesToReport.length) {
			log.reportError("Error - could not find all Samples from PCs in data, only found " + used);
			System.exit(1);
		}
	}

	/**
	 * Load each of the markers to assess (LRR, and AbGenotypesAfterFilters); We load genotypes to allow filtering by homozygous only, gc threshold, etc...
	 * <p>
	 * Data is loaded to assessmentData organized assessmentData[marker0][data for samples in PCs], ditto for abGenotypesAfterFilters
	 * <p>
	 * Note that the order of assessmentData[marker0] does not necessarily reflect the same order as the samples in the pc file
	 * 
	 */
	private void getData() {
		MDL mdl = new MDL(proj, proj.getMarkerSet(), markersToAssess, 2, 100);
		// MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markersToAssess);
		int count = numUsed(samplesToUse);
		this.fullData = new double[markersToAssess.length][count];
		this.abGenotypesAfterFilters = new byte[markersToAssess.length][count];
		float[] lrrs;
		byte[] abGenos;
		ClusterFilterCollection cluster;

		if (Files.exists(proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME))) {
			cluster = ClusterFilterCollection.load(proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME), false);
		} else {
			cluster = new ClusterFilterCollection();
			log.report("Info - did not find the cluster filter file " + proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME) + "; using original genotypes");
		}
		Hashtable<String, Integer> projectIndices = proj.getMarkerIndices();
		if (params != null && recomputeLRR) {
			proj.getLog().reportTimeError("recompute lrr was flagged AND gc correction parameters were passed to data load of median markers");
			return;
		}
		if (params != null) {
			proj.getLog().reportTimeInfo("Will be performing GC correction of median markers");
			if (params.getCentroids() != null) {
				proj.getLog().reportTimeInfo("Will be adjusting data for current centroids");
			}
		}
		// TODO, arg
		int index = 0;
		while (mdl.hasNext()) {
			// for (int index = 0;index < markersToAssess.length; index++) {
			MarkerData markerData = mdl.next();
			if (recomputeLRR) {
				// TODO, I think it is important to cluster on everyone at this point, since we may get inaccurate applications of the clusters
				// i.e Penncnv demands all three genotype clusters or none, we allow 1,2,and 3.
				// this could pose a problem if a genotype for a missing cluster is used
				lrrs = markerData.getRecomputedLRR_BAF(null, null, false, 1, gcThreshold, cluster, true, true, log)[1];
			} else {
				lrrs = markerData.getLRRs();
			}
			if (params != null) {
				lrrs = markerData.getGCCorrectedLRRBAF(params, projectIndices.get(markerData.getMarkerName()), proj.getLog())[1];
			}
			abGenos = markerData.getAbGenotypesAfterFilters(cluster, markersToAssess[index], gcThreshold, log);
			int sampleIndex = 0;
			for (int k = 0; k < samplesToUse.length; k++) {
				if (samplesToUse[k]) {
					fullData[index][sampleIndex] = lrrs[k];
					abGenotypesAfterFilters[index][sampleIndex] = abGenos[k];
					sampleIndex++;
				}
			}
			// markerDataLoader.releaseIndex(index);
			index++;
		}

		// markerDataLoader.reportWaitTimes();
	}

	/**
	 * Compute the median value (after filtering) across the markers to assess for each sample represented in the principal component file Note, medians are computed in a project order manner. A further step may be required if the Pc file is not in the same sample order as the project.
	 */
	public double[] getLRRMedian() {
		double[] medians = new double[fullData[0].length];
		// for sample
		for (int i = 0; i < fullData[0].length; i++) {
			// for marker
			ArrayList<Double> sampLRR = new ArrayList<Double>();
			for (int k = 0; k < fullData.length; k++) {
				// test for null in case we dropped data, test for individual NaN, test for homozygous, test for missing if gcThreshold greater than 0
				if (useMarker(i, k)) {
					sampLRR.add(fullData[k][i]);
				}
			}
			if (sampLRR.size() > 0) {
				medians[i] = Array.median(Array.toDoubleArray(sampLRR));
			} else {
				medians[i] = Double.NaN;
			}
		}
		return medians;
	}

	/**
	 * The main filter for whether to use a marker for a particular sample. The data must not be null, must not be NaN, the gcThreshold must be 0, or, the call must be gte 0 (since we loaded after filters and gcthreshold was taken into account there)
	 * <p>
	 * If the homozygous only option is flagged, the call must be 0 or 2.
	 * 
	 * @param sampIndex
	 *            the current sample
	 * @param markerIndex
	 *            the current marker
	 * @return
	 */
	private boolean useMarker(int sampIndex, int markerIndex) {
		boolean good = false;
		if (fullData[markerIndex] != null && !Double.isNaN(fullData[markerIndex][sampIndex]) && (gcThreshold == 0 || abGenotypesAfterFilters[markerIndex][sampIndex] >= 0) && (abGenotypesAfterFilters[markerIndex][sampIndex] == (byte) 0 || abGenotypesAfterFilters[markerIndex][sampIndex] == (byte) 2 || !homozygousOnly)) {
			good = true;
		}
		return good;
	}

	/**
	 * Prints two files, one corresponding to the raw data, and another boolean file. We do not sort the assessmentData by the order of principal components in the pcfile, so the reporting here may be in a different order.
	 */
	private void printFull() {
		String[] files = { output + MT_REPORT_MARKERS_USED[0], output + MT_REPORT_MARKERS_USED[1] };
		PrintWriter[] writers = getNWriters(proj, files);
		String[] projOrderedSubset = getUsedSubset(allProjSamples, samplesToUse, log);
		if (projOrderedSubset == null) {
			log.reportError("Error - could not print full data, please remove flag or contact jlanej@gmail.com");
			return;
		} else {
			writers[0].println(MT_REPORT[0] + "\t" + Array.toStr(markersToAssess));
			writers[1].println(MT_REPORT[0] + "\t" + Array.toStr(markersToAssess));
			for (int i = 0; i < fullData[0].length; i++) {
				writers[0].print(projOrderedSubset[i]);
				writers[1].print(projOrderedSubset[i]);
				// for marker
				for (int k = 0; k < fullData.length; k++) {
					if (useMarker(i, k)) {
						writers[0].print("\tTRUE");
						writers[1].print("\t" + fullData[k][i]);
					} else {
						writers[0].print("\tFalse");
						writers[1].print("\t.");
					}
				}
				writers[0].println();
				writers[1].println();
			}
			writers[0].close();
			writers[1].close();
		}
	}

	/**
	 * 
	 * @param samples
	 *            all samples in project
	 * @param samplesToUse
	 *            boolean[] representing samples in PCs
	 * @param log
	 * @return the extracted subset of samples represented in the PCs
	 */
	private static String[] getUsedSubset(String[] samples, boolean[] samplesToUse, Logger log) {
		ArrayList<String> used = new ArrayList<String>();
		if (samples.length != samplesToUse.length) {
			log.reportError("Error - mismatched number of samples when extracting samples used, this should not happen");
			return null;
		}
		for (int i = 0; i < samplesToUse.length; i++) {
			if (samplesToUse[i]) {
				used.add(samples[i]);
			}
		}
		return used.toArray(new String[used.size()]);
	}

	private static PrintWriter[] getNWriters(Project proj, String[] fileOuts) {
		PrintWriter[] writers = new PrintWriter[fileOuts.length];
		Logger log = proj.getLog();

		for (int i = 0; i < fileOuts.length; i++) {
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + fileOuts[i])) {
				Files.backup(fileOuts[i], proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			try {
				writers[i] = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + fileOuts[i]));
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + proj.PROJECT_DIRECTORY.getValue() + fileOuts[i] + "\" could not be written to (it's probably open)");
				log.reportException(fnfe);
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + proj.PROJECT_DIRECTORY.getValue() + fileOuts[i] + "\"");
				log.reportException(ioe);
			}
		}
		return writers;
	}

	/**
	 * Loads a principal component file to pcBasis[][] organized as pcBasis[Basis0][Basis for samples]. As the components are loaded, the corresponding sample (DNA) is tracked by the hashtable samplesInPc SamplestoReport maintains the order of samples represented by the pcFile
	 * 
	 * @param pcFile
	 */
	private void loadPcFile(String pcFile, boolean useIID) {
		String pcFilefull;
		if (Files.exists(pcFile) || proj == null) {
			pcFilefull = pcFile;
		} else {
			pcFilefull = proj.PROJECT_DIRECTORY.getValue() + pcFile;
		}
		log.reportTimeInfo("loading principal components from " + pcFilefull);

		SampleData sampleData = null;
		if (!useIID) {
			sampleData = proj.getSampleData(0, false);
		} else {
			log.reportTimeWarning("Using the IID as supplied in the pc file to load samples");
		}
		ArrayList<String> pcSamps = new ArrayList<String>();
		int sampIndex = 0;
		try {
			BufferedReader reader = Files.getReader(pcFilefull, proj == null ? false : proj.JAR_STATUS.getValue(), true, false);
			String[] line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equals("FID") || !line[1].equals("IID")) {
				log.reportError("Error - different format than expected; first column should be FID and second column should be IID, followed by PCs");
				return;
			}
			if ((line.length - 2) < numComponents) {
				log.reportError("Warning - cannot use " + numComponents + " components when only " + (line.length - 2) + " are provided, only loading " + (line.length - 2));
				numComponents = line.length - 2;
			}
			this.pcTitles = Array.subArray(line, 2, numComponents + 2);
			this.totalNumComponents = line.length - 2;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (useIID || sampleData.lookup(line[0] + "\t" + line[1]) != null) {
					samplesInPc.put(useIID ? line[1] : sampleData.lookup(line[0] + "\t" + line[1])[0], numSamples);
					pcSamps.add(useIID ? line[1] : line[0] + "\t" + line[1]);
					numSamples++;
				} else {
					log.reportError("Error - could not find " + line[0] + "\t" + line[1] + " in sample Data");
					return;
				}
			}
			reader.close();
			this.pcBasis = new double[numComponents][numSamples];
			reader = Files.getReader(pcFilefull, proj == null ? false : proj.JAR_STATUS.getValue(), true, false);
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				for (int i = 2; i < numComponents + 2; i++) {
					pcBasis[i - 2][sampIndex] = Double.parseDouble(line[i]);
				}
				sampIndex++;
			}
			reader.close();
			samplesToReport = pcSamps.toArray(new String[pcSamps.size()]);
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + pcFilefull + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + pcFilefull + "\"");
		}
	}

	/**
	 * Used to transpose from PC dominant, to sample dominant double[][] This is necesary for the regression models
	 */
	private static double[][] prepPcs(double[][] pcsBasis) {
		// sample pc values;
		double[][] prepped = new double[pcsBasis[0].length][pcsBasis.length];
		// 20 PC
		for (int i = 0; i < pcsBasis.length; i++) {
			// PC for sample
			for (int j = 0; j < pcsBasis[i].length; j++) {
				prepped[j][i] = pcsBasis[i][j];
			}

		}
		return prepped;
	}

	private static int numUsed(boolean[] samplesToUse) {
		return Array.booleanArraySum(samplesToUse);
	}

	//
	// The following methods are geared toward the cross validation of principal components
	//

	/**
	 * So we can track which pc File went in
	 */
	public String getOutput() {
		return output;
	}

	public int getNumComponents() {
		return numComponents;
	}

	public void setNumComponents(int numComponents) {
		this.numComponents = numComponents;
	}

	public int getNumSamples() {
		return numSamples;
	}

	public int getTotalNumComponents() {
		return totalNumComponents;
	}

	public double[][] getPcBasis() {
		return pcBasis;
	}

	public double[] getMedians() {
		return assesmentData;
	}

	/**
	 * Holds the index of the sample in the pc file, useful for sorting against the project
	 */
	public Hashtable<String, Integer> getSamplesInPc() {
		return samplesInPc;
	}

	public void setMedians(double[] medians) {
		this.assesmentData = medians;
	}

	public void setPcBasis(double[][] pcBasis) {
		this.pcBasis = pcBasis;
	}

	public String getPcFile() {
		return pcFile;
	}

	public PrincipalComponentsResiduals clone() {
		try {
			final PrincipalComponentsResiduals result = (PrincipalComponentsResiduals) super.clone();
			return result;
		} catch (final CloneNotSupportedException e) {
			log.reportError("Error - could not clone the residuals");
			log.reportException(e);
			return null;
		}
	}

	/**
	 * We trim the double[][] holding the basis vectors to a new number of components, this is so we do not have to keep loading the same pc file
	 */
	public static double[][] trimPcBasis(int numComponents, double[][] pcBasis, Logger log) {
		double[][] trimmed = new double[0][];
		if (numComponents > pcBasis.length) {
			log.reportError("Error - we cannot trim the basis vectors to " + numComponents + ", only " + pcBasis.length + " remaining");
			return trimmed;
		} else if (numComponents != pcBasis.length) {
			trimmed = new double[numComponents][];
			for (int i = 0; i < trimmed.length; i++) {
				trimmed[i] = pcBasis[i];
			}
		} else {
			// log.report("Retaining all components");
			trimmed = pcBasis;
		}
		return trimmed;
	}

	/**
	 * @param kFolds
	 *            the number of chunks to use for training/validation
	 * @param numComponents
	 *            the number of components to include, starting from PC 1
	 * @return the results of each training/validation combination
	 */
	public CrossValidation[] crossValidate(int kFolds, int numComponents, LS_TYPE lType) {
		return CrossValidation.kFoldCrossValidate(assesmentData, prepPcs(trimPcBasis(numComponents, pcBasis, log)), kFolds, false, lType, log);
	}

	/**
	 * Cross validate a pc file
	 * 
	 * @param kFolds
	 *            the number of chunks to use for training/validation
	 * @param numComponentsIter
	 *            the number of components to include, starting from PC 1
	 * @param numThreads
	 *            number of threads
	 * @param tmpOutput
	 *            report to a temporary file, if null, this will be skipped
	 * 
	 * @param tmpOutput
	 *            val_pcs If another {@linkPrincipalComponentsResiduals} is provided, this file will become the validation set
	 * @return the results of each training/validation combination
	 */
	public CrossValidation[][] crossValidate(int kFolds, int[] numComponentsIter, int numThreads, LS_TYPE lType, String tmpOutput, PrincipalComponentsResiduals val_pcs) {
		if (tmpOutput != null && Files.exists(tmpOutput)) {
			new File(tmpOutput).delete();
		}
		CrossValidation[][] crossValidations = new CrossValidation[numComponentsIter.length][];
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);// da pool of threads
		ArrayList<Future<CrossValidation[]>> tmpResults = new ArrayList<Future<CrossValidation[]>>();// stores the future CrossValidation[] that will be actualized once the thread has finished
		for (int i = 0; i < numComponentsIter.length; i++) {// need to submit the jobs first
			WorkerPCThread worker = new WorkerPCThread(assesmentData, prepPcs(trimPcBasis(Math.min(numComponentsIter[i], numComponents), pcBasis, log)), kFolds, tmpOutput, lType, val_pcs, log);
			tmpResults.add(executor.submit(worker));// tracks the future object
		}
		for (int i = 0; i < numComponentsIter.length; i++) {
			try {
				crossValidations[i] = tmpResults.get(i).get();// get is only applied after the job has finished
			} catch (InterruptedException e) {
				log.reportException(e);
			} catch (ExecutionException e) {
				log.reportException(e);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(7, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			return crossValidations;
		}
		return crossValidations;
	}

	/**
	 * WorkerThreads which process the validations at each PC requested
	 */

	private static class WorkerPCThread implements Callable<CrossValidation[]> {
		private PrincipalComponentsResiduals val_pcs;// can be null
		private double fullModelR2, fullModelSSerr;
		private double[] deps;
		private double[][] indeps;
		private int kFolds;
		private LS_TYPE lType;
		private String tmpOutput;
		private Logger log;

		public WorkerPCThread(double[] deps, double[][] indeps, int kFolds, String tmpOutput, LS_TYPE lType, PrincipalComponentsResiduals val_pcs, Logger log) {
			super();
			this.deps = deps;
			this.indeps = indeps;
			this.kFolds = kFolds;
			this.tmpOutput = tmpOutput;
			this.fullModelR2 = Double.NaN;
			this.fullModelSSerr = Double.NaN;
			this.lType = lType;
			this.val_pcs = val_pcs;
			this.log = log;
		}

		@Override
		public CrossValidation[] call() {// acts like run
			log.report(ext.getTime() + " Info - Starting validations for PC" + indeps[0].length + " on thread" + Thread.currentThread().getName());
			long time = System.currentTimeMillis();
			CrossValidation[] crossValidation;
			if (val_pcs == null) {
				crossValidation = CrossValidation.kFoldCrossValidate(deps, indeps, kFolds, true, lType, log);
			} else {
				// if we cannot trim to a particular number of components, the crossvalidation will fail and will be ignored
				double[][] val_basis = prepPcs(trimPcBasis(Math.min(indeps[0].length, val_pcs.getNumComponents()), val_pcs.getPcBasis(), log));
				crossValidation = CrossValidation.kFoldCrossValidateOutSample(deps, indeps, val_pcs.getMedians(), val_basis, kFolds, true, lType, log);
			}
			log.report(ext.getTime() + " Info - finished validations for PC" + indeps[0].length + " and took " + ext.getTimeElapsed(time));
			char S = 'S';
			computeFullModel();
			writeToTmpFile(crossValidation, ext.getTimeSince(time, S) + "", log);
			return crossValidation;
		}

		private void computeFullModel() {
			if (deps.length > indeps[0].length + 1) {
				RegressionModel model = (RegressionModel) new LeastSquares(deps, indeps, null, false, true, lType);
				if (!model.analysisFailed()) {
					fullModelR2 = model.getRsquare();
					fullModelSSerr = computeFullModelSSerr(model.getResiduals());
				}
			}
		}

		private static double computeFullModelSSerr(double[] residuals) {
			double SSerr = 0;
			for (int i = 0; i < residuals.length; i++) {
				SSerr += Math.pow(residuals[i], 2);
			}
			return SSerr;
		}

		/**
		 * We write to a temporary file so we can monitor progress
		 */
		private void writeToTmpFile(CrossValidation[] crossValidation, String time, Logger log) {
			if (tmpOutput != null) {
				if (!Files.exists(tmpOutput)) {
					Files.write(Array.toStr(MT_RESIDUAL_CROSS_VALIDATED_REPORT), tmpOutput);
				}
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(tmpOutput, true));
					writer.println(ext.getTime() + "\t" + time + "\t" + indeps[0].length + "\t" + CrossValidation.getEstimateError(crossValidation) + "\t" + CrossValidation.getAverageR2(crossValidation) + "\t" + CrossValidation.getAverageSEbetas(crossValidation) + "\t" + fullModelR2 + "\t" + fullModelSSerr);
					writer.close();
				} catch (Exception e) {
					log.reportError("Error writing to " + tmpOutput);
					log.reportException(e);
				}
			}
		}
	}

	// the following methods aid in correcting XY intensities
	/**
	 * Constructor mainly used by {@link PrincipalComponentsIntensity} Essentially handles the loading of the PCs
	 *
	 */

	public PrincipalComponentsResiduals(PrincipalComponentsResiduals principalComponentsResiduals) {
		this.numComponents = principalComponentsResiduals.getNumComponents();
		this.pcBasis = principalComponentsResiduals.getPcBasis();
		this.samplesToUse = principalComponentsResiduals.getSamplesToUse();
		this.proj = principalComponentsResiduals.getProj();
		this.samplesInPc = principalComponentsResiduals.getSamplesInPc();
		this.allProjSamples = principalComponentsResiduals.getAllProjSamples();
		this.sortedByProject = principalComponentsResiduals.isSortedByProject();
		this.log = proj.getLog();

	}

	public String[] getAllProjSamples() {
		return allProjSamples;
	}

	public Project getProj() {
		return proj;
	}

	/**
	 * @return boolean array representing all samples in the project that have a pc in this object
	 */
	public boolean[] getSamplesToUse() {
		return samplesToUse;
	}

	/**
	 * transposes from PC dominant to sample dominant (to be used for regressions etc...)
	 */
	public double[][] getPreppedPCs() {
		return prepPcs(pcBasis);
	}

	/**
	 * Transposes from PC dominant to sample dominant and trims the proper number of PCs, optionally can be sorted by project, if needed See {@link PrincipalComponentsResiduals#isSortedByProject()} and {@link PrincipalComponentsResiduals#determineSortedByProject()} to see if the file is sorted already
	 * 
	 * @param numComponents
	 * @param sortByProject
	 * @return
	 */
	public double[][] getTrimmedPreppedPCs(int numComponents, boolean sortByProject) {
		double[][] trimmedPrepped = prepPcs(trimPcBasis(numComponents, pcBasis, log));
		if (sortByProject) {
			return trimmedPrepped;
		} else {
			String[] samples = proj.getSamples();
			double[][] tmp = new double[samples.length][numComponents];
			for (int i = 0; i < tmp.length; i++) {
				if (!samplesInPc.containsKey(samples[i])) {
					Arrays.fill(tmp[i], Double.NaN);
				} else {
					tmp[i] = trimmedPrepped[samplesInPc.get(samples[i])];
				}
			}
			return tmp;
		}
	}

	public boolean isSortedByProject() {
		return sortedByProject;
	}

	public CrossValidation getCorrectedDataAt(double[] data, boolean[] samplesTobuildModel, int numComponentsForModel, LS_TYPE lType, String title, boolean verbose) {
		return getCorrectedDataAt(data, null, samplesTobuildModel, numComponentsForModel, lType, title, verbose);
	}

	/**
	 * @param data
	 *            the data to be corrected with the PCs {@link PrincipalComponentsResiduals } object, currently the data must correspond in length and order to all samples in the project. Set to NaN to mask
	 * @param extraIndeps
	 *            this can be null, but if not, this data will be included in the regression model with the PCs
	 * @param samplesTobuildModel
	 *            if not null, must correpsond in length and order to project samples, these individual will be used to build the regression model
	 * @param numComponentsForModel
	 *            number of components to use, must be >0
	 * @param svdRegression
	 *            if the number of components is large, use an svd based regression
	 * @param title
	 *            a title string for error reporting (marker name, phenotype, etc) * <br>
	 *            NOTE: the pc object does not have to be in project order, or contain all samples, but it should
	 * 
	 * @return a computed {@link CrossValidation}
	 * 
	 */
	public CrossValidation getCorrectedDataAt(double[] data, double[][] extraIndeps, boolean[] samplesTobuildModel, int numComponentsForModel, LS_TYPE lType, String title, boolean verbose) {
		int numSamples = proj.getSamples().length;
		boolean go = true;
		CrossValidation cval;
		if (numComponentsForModel <= 0) {
			if (extraIndeps == null) {
				proj.getLog().reportError("Error - number of components specified must be greater than 0");
				go = false;
			} else {
				proj.getLog().reportTimeWarning(numComponentsForModel + " components were specified, so only using the extra independen variables provided");
			}

		}
		if (data == null || willFailNAN(data, numComponentsForModel)) {
			int numNonNaN = data == null ? 0 : Array.removeNaN(data).length;
			proj.getLog().reportError("Error - there are not enough samples with non NAN (n=" + numNonNaN + ") data for " + title + " using " + numComponentsForModel + " " + (numComponentsForModel == 1 ? "principal component " : "principal components") + " to run a regression");
			go = false;
		}
		if (data.length != numSamples) {
			proj.getLog().reportError("Error - data points must be provided for every sample in the project, only found " + data.length);
			proj.getLog().reportError("      - data for unwanted samples can be set to NaN and masked ");
			go = false;
		}
		if (samplesTobuildModel != null && samplesTobuildModel.length != numSamples) {
			proj.getLog().reportError("Error - boolean definitions of samples to build the regression model must be provided for every sample in the project, only found  " + samplesTobuildModel.length);
			go = false;
		}

		if (numComponentsForModel > numComponents) {
			proj.getLog().reportError("Error - too many components (" + numComponentsForModel + ") were specified for the model, only have " + numComponents + " available");
			go = false;
		}

		if (extraIndeps != null && extraIndeps.length != samplesInPc.size()) {
			proj.getLog().reportError("The size of the independent variable array did not match the number of samples in the pc data");
			go = false;
		}
		if (!sortedByProject) {// this could happen if the extrapolated pcs are not used...encourage the full pc file since it is matched to the project
			proj.getLog().report("Warning - detected that the PC file is not perfectly matched to the project");
			if (numSamples != samplesInPc.size()) {
				proj.getLog().report("Warning - missing individuals (" + (numSamples - samplesInPc.size()) + ") will be set to NaN across all PCs");
			}
			proj.getLog().report("Warning - an extra step will be taken to allign the principal components to individuals in the project");
		}
		if (!go) {
			cval = new CrossValidation(new double[0], new double[0][0], new double[0], new double[0][0], true, lType, proj.getLog());
			cval.setAnalysisFailed(true);
		} else {
			// if (lType == LS_TYPE.SVD) {
			// lType = numComponentsForModel > NUM_PC_SVD_OVERIDE ? LS_TYPE.SVD : LS_TYPE.REGULAR;
			// if (lType == LS_TYPE.REGULAR) {
			// log.reportTimeWarning("Over-riding SVD method since " + numComponentsForModel + " < " + NUM_PC_SVD_OVERIDE);
			// }
			// }

			double[] train_deps = (samplesTobuildModel == null ? data : Array.subArray(data, samplesTobuildModel));
			double[][] train_indeps = numComponentsForModel > 0 ? getTrimmedPreppedIndepsProjectPCsFor(samplesTobuildModel, extraIndeps, numComponentsForModel, log) : Array.subArray(extraIndeps, samplesTobuildModel);
			double[][] val_indeps = numComponentsForModel > 0 ? getTrimmedPreppedIndepsProjectPCsFor(null, extraIndeps, numComponentsForModel, log) : extraIndeps;

			cval = new CrossValidation(train_deps, train_indeps, data, val_indeps, verbose, lType, proj.getLog());
			cval.train();
			cval.computePredictedValues();
			cval.computeResiduals();
		}
		return cval;
	}

	/**
	 */
	public static CrossValidation getCorrectedDataAt(PrincipalComponentsResiduals principalComponentsResiduals, float[] data, boolean[] samplesTobuildModel, int numComponentsForModel, LS_TYPE lType, String title) {
		return principalComponentsResiduals.getCorrectedDataAt(Array.toDoubleArray(data), samplesTobuildModel, numComponentsForModel, lType, title, true);
	}

	/**
	 * @param toExtract
	 * @param additionalData
	 *            tags this extra data onto the pcs, does not do any error checking
	 * @param numComponents
	 * @return
	 */
	private double[][] getTrimmedPreppedIndepsProjectPCsFor(boolean[] toExtract, double[][] additionalData, int numComponents, Logger log) {
		double[][] preppedPcs = getTrimmedPreppedProjectPCsFor(toExtract, numComponents);
		if (additionalData != null) {
			double[][] tmp = new double[preppedPcs.length][numComponents + additionalData[0].length];
			double[][] tmpAdd = toExtract == null ? additionalData : Array.subArray(additionalData, toExtract);
			if (preppedPcs.length != tmpAdd.length) {
				log.reportError("Mismatched array addition for additional data");
				return null;
			}
			for (int i = 0; i < preppedPcs.length; i++) {
				tmp[i] = Array.concatDubs(preppedPcs[i], tmpAdd[i]);
			}
			preppedPcs = tmp;
		}

		return preppedPcs;

	}

	/**
	 * Used to extract a subset of individuals from the PCs,
	 * 
	 * @param toExtract
	 *            must be the same length as the samples in the project
	 */
	public double[][] getTrimmedPreppedProjectPCsFor(boolean[] toExtract, int numComponents) {
		double[][] prepped;
		prepped = getTrimmedPreppedPCs(numComponents, sortedByProject);
		if (toExtract == null) {
			return prepped;
		} else {
			double[][] preppedFor = new double[Array.booleanArraySum(toExtract)][];
			int index = 0;
			for (int i = 0; i < toExtract.length; i++) {
				if (toExtract[i]) {
					preppedFor[index] = prepped[i];
					index++;
				}
			}
			return preppedFor;
		}
	}

	/**
	 * @return whether all the samples for a project are present, and whether they are in the same order
	 */
	private boolean determineSortedByProject() {
		if (allProjSamples.length != Array.booleanArraySum(samplesToUse)) {
			return false;
		} else {
			for (int i = 0; i < allProjSamples.length; i++) {
				if (!samplesInPc.containsKey(allProjSamples[i]) || samplesInPc.get(allProjSamples[i]) != i) {
					return false;
				}
			}
			return true;
		}
	}

	private static boolean willFailNAN(double[] data, float numComponents) {
		int count = 0;
		for (int i = 0; i < data.length; i++) {
			if (!Double.isNaN(data[i])) {
				count++;
			}
			if (count > numComponents) {
				return false;
			}
		}
		return true;
	}

	/**
	 * @param data
	 *            Currently, the data must be the same length as all samples in the project. To remove individuals from the regression, set values to NaN
	 * @param numComponents
	 *            number of components to include in the regression model
	 * @param title
	 *            optional to report if model has failed (could be a specific marker or trait, etc...) OLD OLD
	 * @return
	 */
	public double[] getCorrectdedDataAt(double[] data, int numComponents, String title) {
		String finalTitle = (title == null ? " data" : title);
		double[] correctedData = new double[data.length];
		if (data == null || willFailNAN(data, numComponents)) {
			proj.getLog().reportError("Error - there are not enough samples with non NAN data for " + title + " using " + numComponents + " " + (numComponents == 1 ? "principal component " : "principal components") + " to run a regression, skipping");
			return null;
		}
		if (data.length != samplesToUse.length) {// samplesToUse in project order;
			proj.getLog().reportError("Error - array of samples to correct and data must be the same length");
			return null;
		} else {
			setAssesmentDataSortByPCs(Array.subArray(data, samplesToUse));// sorts according to the order of samples in the pcFile
			RegressionModel model = (RegressionModel) new LeastSquares(assesmentData, getTrimmedPreppedPCs(numComponents, false));
			if (!model.analysisFailed()) {
				double[] residuals = model.getResiduals();
				int indexResid = 0;
				for (int i = 0; i < correctedData.length; i++) {
					if (samplesToUse[i] && !Double.isNaN(data[i])) {
						correctedData[i] = residuals[indexResid];
						indexResid++;
					} else {
						correctedData[i] = Double.NaN;
					}
				}
			} else {
				Arrays.fill(correctedData, Float.NaN);
				proj.getLog().report("Warning - could not correct" + finalTitle + ", regression model has failed...setting all values to NaN");
			}
		}
		return correctedData;
	}

	public void setMarkersToAssessFile(String markersToAssessFile) {
		this.markersToAssessFile = markersToAssessFile;
	}

	public void setRecomputeLRR(boolean recomputeLRR) {
		this.recomputeLRR = recomputeLRR;
	}

	public void setFullData(double[][] fullData) {
		this.fullData = fullData;
	}

	public double[][] getFullData() {
		return fullData;
	}

	public void setAbGenotypesAfterFilters(byte[][] abGenotypesAfterFilters) {
		this.abGenotypesAfterFilters = abGenotypesAfterFilters;
	}

	/**
	 * @param PC
	 *            one - based pc basis to extract
	 * @return
	 */
	public double[] getBasisAt(int PC) {
		double[] basis = null;
		if (PC <= 0) {
			log.reportTimeError("Requested PC must be greater than 0 (one -based extraction");
			return basis;
		}
		if (PC > numComponents) {
			log.reportTimeError("Requested PC must be less than or equal to the total number of components (" + numComponents + ")");
			return basis;
		} else {
			return pcBasis[PC - 1];
		}
	}

	private boolean verifyDataSampleSize(double[] data, Logger log) {
		boolean matched = true;
		if (data.length != samplesInPc.size()) {
			matched = false;
			log.reportTimeError("Input data (" + data.length + ") does not match the number of samples in the pc file (" + samplesInPc.size() + ")");
			log.reportTimeError("Consider masking ");

		}
		return matched;
	}

	/**
	 * @param data
	 *            organized as data[variable][dataForAllSamples]
	 * @param extraIndeps
	 *            organized as extraIndeps[sample][indepsForSample]
	 * @param samplesForRanking
	 *            stats will be calculated with these samples only
	 * @param titles
	 *            must have same length as data, essentially variable titles
	 * @param statType
	 *            see {@link STAT_TYPE}
	 * @param rankType
	 *            see {@link VALUE_TYPE}
	 * @param log
	 * @return the {@link StatsCrossTabRank } array for all data requested
	 * 
	 *         NOTE: the data is only ranked against individual PCs, not
	 */
	public StatsCrossTabRank[] getStatRanksFor(double[][] data, double[][] extraIndeps, boolean[] samplesForRanking, String[] titles, STAT_TYPE statType, VALUE_TYPE rankType, boolean stepwise, Logger log) {
		StatsCrossTabRank[] sRanks = new StatsCrossTabRank[data.length];
		for (int i = 0; i < sRanks.length; i++) {
			sRanks[i] = getStatRankFor(data[i], extraIndeps, samplesForRanking, titles[i], statType, rankType, stepwise, 1, log);
		}
		return sRanks;
	}

	public StatsCrossTabRank getStatRankFor(final double[] data, final double[][] extraIndeps, final boolean[] samplesForRanking, final String title, final STAT_TYPE statType, final VALUE_TYPE rankType, boolean stepwise, int numThreads, Logger log) {
		String[] allTitles = Array.concatAll(new String[] { title }, pcTitles);
		double[][] tmp = extraIndeps;
		if (verifyDataSampleSize(data, log)) {

			double[][] basis = getPcBasis();
			double[][] toRank = new double[basis.length + 1][];
			toRank[0] = samplesForRanking == null ? data : Array.subArray(data, samplesForRanking);
			for (int i = 0; i < basis.length; i++) {
				toRank[i + 1] = samplesForRanking == null ? basis[i] : Array.subArray(basis[i], samplesForRanking);
			}
			if (extraIndeps != null && samplesForRanking != null) {
				tmp = Array.subArray(extraIndeps, samplesForRanking);
			}
			if (stepwise) {
				double[][] basisToStep = new double[basis.length][];
				for (int i = 0; i < basis.length; i++) {
					basisToStep[i] = samplesForRanking == null ? basis[i] : Array.subArray(basis[i], samplesForRanking);
				}
				Stepwise stepwiseIt = new Stepwise(samplesForRanking == null ? data : Array.subArray(data, samplesForRanking), samplesForRanking == null ? basis : Array.subArray(getPreppedPCs(), samplesForRanking), NUM_PC_SVD_OVERIDE, true, numThreads);
				stepwiseIt.setVarNames(getPcTitles());
				StepWiseSummary stepWiseSummary = stepwiseIt.getStepWiseSummary(NUM_PC_SVD_OVERIDE, numThreads);

				if (stepWiseSummary != null) {
					StatsCrossTabRank statsCrossTabRank = new StatsCrossTabRank(title, stepWiseSummary.getOrderOfOriginal(), stepWiseSummary.getSigs(), stepWiseSummary.getStats(), Array.subArray(getPcTitles(), stepWiseSummary.getOrderOfOriginal()));
					return statsCrossTabRank;
				} else {
					log.reportTimeWarning("Stepwise regression did not find any signifcant variables");
					StatsCrossTabRank statsCrossTabRank = new StatsCrossTabRank(title, new int[] {}, new double[] {}, new double[] {}, new String[] {});
					return statsCrossTabRank;

				}

			} else {
				boolean[] mask = Array.booleanArray(toRank.length, false);
				mask[0] = true;
				StatsCrossTabs sCrossTabs = new StatsCrossTabs(toRank, tmp, mask, allTitles, statType, true, log);
				sCrossTabs.computeTable();
				return sCrossTabs.getInOrder(0, rankType, log);
			}
		}
		return null;
	}

	/**
	 * Constructor for non-project associated (such as {@link PrincipalComponentsPlink}) principal component manipulations
	 */
	public PrincipalComponentsResiduals(String pcFile, int numComponents, Logger log) {
		super();
		this.numComponents = numComponents;
		this.log = log;
		this.numSamples = 0;
		this.samplesInPc = new Hashtable<String, Integer>();
		this.pcFile = pcFile;
		this.useIID = true;
		loadPcFile(pcFile, useIID);
		this.sortedByProject = true;
	}

	public String[] getSamplesToReport() {
		return samplesToReport;
	}

	/**
	 * Constructor mainly used by {@link PrincipalComponentsIterator}
	 *
	 */

	public PrincipalComponentsResiduals(PrincipalComponentsResiduals principalComponentsResiduals, double[][] basis, double[] assesmentData) {
		this.numComponents = basis.length;
		this.pcBasis = basis;
		this.samplesToUse = principalComponentsResiduals.getSamplesToUse();
		this.proj = principalComponentsResiduals.getProj();
		this.samplesInPc = principalComponentsResiduals.getSamplesInPc();
		this.allProjSamples = principalComponentsResiduals.getAllProjSamples();
		this.sortedByProject = principalComponentsResiduals.isSortedByProject();
		this.assesmentData = assesmentData;
		this.log = proj.getLog();

	}

	public static class PrincipalComponentsIterator implements Iterator<PrincipalComponentsResiduals>, Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private PrincipalComponentsResiduals pcResids;
		private int[] order;
		private int index;

		public PrincipalComponentsIterator(PrincipalComponentsResiduals pcResids, int[] order) {
			super();
			this.pcResids = pcResids;
			this.order = order;
			this.index = 0;// one based extraction, but the first iter is 0 components
		}

		@Override
		public boolean hasNext() {
			return index < (order == null ? pcResids.getNumComponents() + 1 : order.length + 1); //+1 to account for PC0

		}

		public int[] getOrder() {
			return order;
		}

		public PrincipalComponentsResiduals getPcResids() {
			return pcResids;
		}

		@Override
		public PrincipalComponentsResiduals next() {
			double[][] newBasis = new double[index][];
			if (index > 0) {
				for (int i = 0; i < index; i++) {
					newBasis[i] = pcResids.getBasisAt(order == null ? i + 1 : order[i]);
				}
			}
			PrincipalComponentsResiduals newPcResiduals = new PrincipalComponentsResiduals(pcResids, newBasis, pcResids.getMedians());

			index++;
			return newPcResiduals;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}
	}

	/**
	 * Reformats the pcBasis so that the pcs are in strict project order, with missing samples set to Nan.<br>
	 * Only valid for a subset of the samples, i.e. PCs cannot contain non-project samples<br>
	 */
	public void fillInMissing() {
		String[] samples = proj.getSamples();
		if (ext.containsAll(samples, samplesInPc.keySet().toArray(new String[samplesInPc.size()]))) {
			log.reportTimeInfo("All project samples detected in pc file");
		} else {
			Arrays.fill(samplesToUse, true);
			log.reportTimeInfo("Reformatting pc basis");

			double[][] tmpBasis = new double[numComponents][samples.length];
			Hashtable<String, Integer> tmpSampsInPC = new Hashtable<String, Integer>();
			int[] sampleMap = new int[samples.length];
			Arrays.fill(sampleMap, -1);
			for (int i = 0; i < samples.length; i++) {
				if (ext.indexOfStr(samples[i], samples) < 0) {
					throw new IllegalStateException("All sample must come from the project " + proj.getPropertyFilename() + " , did not find " + samples[i] + " in the project");
				} else if (samplesInPc.containsKey(samples[i])) {
					sampleMap[i] = samplesInPc.get(samples[i]);
					tmpSampsInPC.put(samples[i], i);
				} else {
					// samplesToUse[i] = false;
					tmpSampsInPC.put(samples[i], i);
				}
			}
			log.reportTimeInfo("Filling in missing samples with NaN");

			for (int i = 0; i < pcBasis.length; i++) {

				for (int j = 0; j < samples.length; j++) {
					if (sampleMap[j] >= 0) {
						tmpBasis[i][j] = pcBasis[i][sampleMap[j]];
					} else {
						tmpBasis[i][j] = Double.NaN;
					}
				}
			}
			log.reportTimeInfo(Array.booleanArraySum(samplesToUse) + " project samples had PC values, " + (samples.length - Array.booleanArraySum(samplesToUse)) + " samples were set to NaN");

			samplesInPc = tmpSampsInPC;
			samplesToReport = samples;
			allProjSamples = samples;
			pcBasis = tmpBasis;
			sortedByProject = true;
			log.report(samplesInPc.size() + ": hash samples ; " + samplesToReport.length + " report Samples");
			log.report(allProjSamples.length + ": project samples ; " + pcBasis[0].length + " pc Samples");
		}
	}
}
// int trainIndex = 0;
// for (int i = 0; i < proj.getSamples().length; i++) {
// if (samplesTobuildModel[i]) {
// // System.out.println("TRAIN\t"+proj.getSamples()[trainIndex] + "\t" + data[trainIndex] + "\t" + Array.toStr(train_indeps[trainIndex]));
// if (train_deps[trainIndex] != data[i] || !Arrays.equals(train_indeps[trainIndex], val_indeps[i])) {
// System.out.println("DSFDSF\t" + proj.getSamples()[trainIndex] + "\t" + proj.getSamples()[i]);
// System.out.println("TRAIN\t" + proj.getSamples()[trainIndex] + "\t" + data[trainIndex] + "\t" + Array.toStr(train_indeps[trainIndex]));
// }
// trainIndex++;
//
// }
// // System.out.println("REG\t"+proj.getSamples()[i] + "\t" + data[i] + "\t" + Array.toStr(train_indeps[i]));
//
// }

// for (int i = 0; i < proj.getSamples().length; i++) {
// System.out.println(proj.getSamples()[i] + "\t" + data[i]);
// }
