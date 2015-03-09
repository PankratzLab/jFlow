package cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import stats.CrossValidation;
import stats.LeastSquares;
import stats.RegressionModel;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * <p>
 * Class to compute residuals wrt PCs and summarize median values, currently aimed at Mitochondrial copy number, but is extensible to other data
 * 
 *
 */
public class PrincipalComponentsResiduals implements Cloneable {
	private static final String[] MT_REPORT = { "DNA", "FID", "IID", "Sex", "median_MT_LRR_raw", "median_MT_LRR_PC_residuals", "median_MT_LRR_PC_residuals_inverseTransformed" };
	private static final String[] MT_REPORT_EXT = { ".report.txt" };
	private static final String[] MT_REPORT_MARKERS_USED = { ".MedianMarkers.MarkersUsed.txt", ".MedianMarkers.RawValues.txt" };
	private static final String[] MT_RESIDUAL_CROSS_VALIDATED_REPORT = { "Time Completed", "Time to complete(seconds)", "PC", "Cross-validation Average SSerr", "Cross-validation Average R-squared", "Average Standard Error of Betas", "Full model R-squared", "Full model SSerr" };

	private String markersToAssessFile, output, residOutput, pcFile;
	private String[] markersToAssess, samplesToReport, allProjSamples;
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

	/**
	 * Compute the median data in project order, and then set back to pcfile order, just in case
	 */
	public void computeAssessmentDataMedians() {
		this.markersToAssess = PrincipalComponentsCompute.sortByProjectMarkers(proj, HashVec.loadFileToStringArray(proj.getProjectDir() + markersToAssessFile, false, new int[] { 0 }, true));
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
		RegressionModel model = (RegressionModel) new LeastSquares(assesmentData, prepPcs(pcBasis));
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

	/**
	 * Summarize the medians, residuals from the PCs, and inverse transformed values in the same order as the input pc file
	 */
	public void summarize(String output) {
		this.residOutput = ext.rootOf(output) + MT_REPORT_EXT[0];
		SampleData sampleData = proj.getSampleData(0, false);
		try {
			if (Files.exists(proj.getProjectDir() + residOutput)) {
				Files.backup(residOutput, proj.getProjectDir(), proj.getProjectDir() + proj.getProperty(Project.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(proj.getProjectDir() + residOutput));
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
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + residOutput + "\"");
			log.reportException(ioe);
			System.exit(2);
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
	 * Note that the order of assessmentData[marker0] does not neccesarily reflect the same order as the samples in the pc file
	 * 
	 */
	private void getData() {
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markersToAssess);
		int count = numUsed(samplesToUse);
		this.fullData = new double[markersToAssess.length][count];
		this.abGenotypesAfterFilters = new byte[markersToAssess.length][count];
		float[] lrrs;
		byte[] abGenos;
		ClusterFilterCollection cluster;

		if (Files.exists(proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME))) {
			cluster = ClusterFilterCollection.load(proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME), false);
		} else {
			cluster = new ClusterFilterCollection();
			log.report("Info - did not find the cluster filter file " + proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME) + "; using original genotypes");
		}

		for (int i = 0; i < markersToAssess.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			if (recomputeLRR) {
				// TODO, I think it is important to cluster on everyone at this point, since we may get inaccurate applications of the clusters
				// i.e Penncnv demands all three genotype clusters or none, we allow 1,2,and 3.
				// this could pose a problem if a genotype for a missing cluster is used
				lrrs = markerData.getRecomputedLRR_BAF(null, null, false, 1, gcThreshold, cluster, true, true, log)[1];
			} else {
				lrrs = markerData.getLRRs();
			}
			abGenos = markerData.getAbGenotypesAfterFilters(cluster, markersToAssess[i], gcThreshold);
			int sampleIndex = 0;
			for (int k = 0; k < samplesToUse.length; k++) {
				if (samplesToUse[k]) {
					fullData[i][sampleIndex] = lrrs[k];
					abGenotypesAfterFilters[i][sampleIndex] = abGenos[k];
					sampleIndex++;
				}
			}
			markerDataLoader.releaseIndex(i);
		}

		markerDataLoader.reportWaitTimes();
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
			if (Files.exists(proj.getProjectDir() + fileOuts[i])) {
				Files.backup(fileOuts[i], proj.getProjectDir(), proj.getProjectDir() + proj.getProperty(Project.BACKUP_DIRECTORY));
			}
			try {
				writers[i] = new PrintWriter(new FileWriter(proj.getProjectDir() + fileOuts[i]));
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + proj.getProjectDir() + fileOuts[i] + "\" could not be written to (it's probably open)");
				log.reportException(fnfe);
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + proj.getProjectDir() + fileOuts[i] + "\"");
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
			pcFilefull = proj.getProjectDir() + pcFile;
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
			BufferedReader reader = Files.getReader(pcFilefull, proj == null ? false : proj.getJarStatus(), true, false);
			String[] line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equals("FID") || !line[1].equals("IID")) {
				log.reportError("Error - different format than expected; first column should be FID and second column should be IID, followed by PCs");
				return;
			}
			if ((line.length - 2) < numComponents) {
				log.reportError("Warning - cannot use " + numComponents + " components when only " + (line.length - 2) + " are provided, only loading " + (line.length - 2));
				numComponents = line.length - 2;
			}
			this.totalNumComponents = line.length - 2;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (useIID || sampleData.lookup(line[0] + "\t" + line[1]) != null) {
					samplesInPc.put(useIID ? line[1] : sampleData.lookup(line[0] + "\t" + line[1])[0], numSamples);
					pcSamps.add(useIID ? line[1] : line[0] + "\t" + line[1]);
					numSamples++;
				} else {
					log.reportError("Error - could not find " + line[0] + "\t" + line[1] + " in sample Data");
					System.exit(1);
				}
			}
			reader.close();
			this.pcBasis = new double[numComponents][numSamples];
			reader = Files.getReader(pcFilefull, proj == null ? false : proj.getJarStatus(), true, false);
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
	public CrossValidation[] crossValidate(int kFolds, int numComponents, boolean svdRegression) {
		return CrossValidation.kFoldCrossValidate(assesmentData, prepPcs(trimPcBasis(numComponents, pcBasis, log)), kFolds, false, svdRegression, log);
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
	public CrossValidation[][] crossValidate(int kFolds, int[] numComponentsIter, int numThreads, boolean svdRegression, String tmpOutput, PrincipalComponentsResiduals val_pcs) {
		if (tmpOutput != null && Files.exists(tmpOutput)) {
			new File(tmpOutput).delete();
		}
		CrossValidation[][] crossValidations = new CrossValidation[numComponentsIter.length][];
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);// da pool of threads
		ArrayList<Future<CrossValidation[]>> tmpResults = new ArrayList<Future<CrossValidation[]>>();// stores the future CrossValidation[] that will be actualized once the thread has finished
		for (int i = 0; i < numComponentsIter.length; i++) {// need to submit the jobs first
			WorkerPCThread worker = new WorkerPCThread(assesmentData, prepPcs(trimPcBasis(Math.min(numComponentsIter[i], numComponents), pcBasis, log)), kFolds, tmpOutput, svdRegression, val_pcs, log);
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
		private boolean svdRegression;
		private String tmpOutput;
		private Logger log;

		public WorkerPCThread(double[] deps, double[][] indeps, int kFolds, String tmpOutput, boolean svdRegression, PrincipalComponentsResiduals val_pcs, Logger log) {
			super();
			this.deps = deps;
			this.indeps = indeps;
			this.kFolds = kFolds;
			this.tmpOutput = tmpOutput;
			this.fullModelR2 = Double.NaN;
			this.fullModelSSerr = Double.NaN;
			this.svdRegression = svdRegression;
			this.val_pcs = val_pcs;
			this.log = log;
		}

		@Override
		public CrossValidation[] call() {// acts like run
			log.report(ext.getTime() + " Info - Starting validations for PC" + indeps[0].length + " on thread" + Thread.currentThread().getName());
			long time = System.currentTimeMillis();
			CrossValidation[] crossValidation;
			if (val_pcs == null) {
				crossValidation = CrossValidation.kFoldCrossValidate(deps, indeps, kFolds, true, svdRegression, log);
			} else {
				// if we cannot trim to a particular number of components, the crossvalidation will fail and will be ignored
				double[][] val_basis = prepPcs(trimPcBasis(Math.min(indeps[0].length, val_pcs.getNumComponents()), val_pcs.getPcBasis(), log));
				crossValidation = CrossValidation.kFoldCrossValidateOutSample(deps, indeps, val_pcs.getMedians(), val_basis, kFolds, true, svdRegression, log);
			}
			log.report(ext.getTime() + " Info - finished validations for PC" + indeps[0].length + " and took " + ext.getTimeElapsed(time));
			char S = 'S';
			computeFullModel();
			writeToTmpFile(crossValidation, ext.getTimeSince(time, S) + "");
			return crossValidation;
		}

		private void computeFullModel() {
			if (deps.length > indeps[0].length + 1) {
				RegressionModel model = (RegressionModel) new LeastSquares(deps, indeps, null, false, true, svdRegression);
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
		private void writeToTmpFile(CrossValidation[] crossValidation, String time) {
			if (tmpOutput != null) {
				if (!Files.exists(tmpOutput)) {
					Files.write(Array.toStr(MT_RESIDUAL_CROSS_VALIDATED_REPORT), tmpOutput);
				}
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(tmpOutput, true));
					writer.println(ext.getTime() + "\t" + time + "\t" + indeps[0].length + "\t" + CrossValidation.getEstimateError(crossValidation) + "\t" + CrossValidation.getAverageR2(crossValidation) + "\t" + CrossValidation.getAverageSEbetas(crossValidation) + "\t" + fullModelR2 + "\t" + fullModelSSerr);
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + tmpOutput);
					System.err.println(e);
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

	/**
	 * @param data
	 *            the data to be corrected with the PCs {@link PrincipalComponentsResiduals } object, currently the data must correspond in length and order to all samples in the project. Set to NaN to mask
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
	public CrossValidation getCorrectedDataAt(double[] data, boolean[] samplesTobuildModel, int numComponentsForModel, boolean svdRegression, String title, boolean verbose) {
		int numSamples = proj.getSamples().length;
		boolean go = true;
		CrossValidation cval;
		if (numComponentsForModel <= 0) {
			proj.getLog().reportError("Error - number of components specified must be greater than 0");
			go = false;
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
		if (!sortedByProject) {// this could happen if the extrapolated pcs are not used...encourage the full pc file since it is matched to the project
			proj.getLog().report("Warning - detected that the PC file is not perfectly matched to the project");
			if (numSamples != samplesInPc.size()) {
				proj.getLog().report("Warning - missing individuals (" + (numSamples - samplesInPc.size()) + ") will be set to NaN across all PCs");
			}
			proj.getLog().report("Warning - an extra step will be taken to allign the principal components to individuals in the project");
		}
		if (!go) {
			cval = new CrossValidation(new double[0], new double[0][0], new double[0], new double[0][0], true, svdRegression, proj.getLog());
			cval.setAnalysisFailed(true);
		} else {
			double[] train_deps = (samplesTobuildModel == null ? data : Array.subArray(data, samplesTobuildModel));
			double[][] train_indeps = getTrimmedPreppedProjectPCsFor(samplesTobuildModel, numComponentsForModel);
			double[][] val_indeps = getTrimmedPreppedProjectPCsFor(null, numComponentsForModel);
			cval = new CrossValidation(train_deps, train_indeps, data, val_indeps, verbose, svdRegression, proj.getLog());
			cval.train();
			cval.computePredictedValues();
			cval.computeResiduals();
		}
		return cval;
	}

	/**
	 */
	public static CrossValidation getCorrectedDataAt(PrincipalComponentsResiduals principalComponentsResiduals, float[] data, boolean[] samplesTobuildModel, int numComponentsForModel, boolean svdRegression, String title) {
		return principalComponentsResiduals.getCorrectedDataAt(Array.toDoubleArray(data), samplesTobuildModel, numComponentsForModel, svdRegression, title, true);
	}

	/**
	 * Used to extract a subset of individuals from the PCs,
	 * 
	 * @param toExtract
	 *            must be the same length as the samples in the project
	 */
	private double[][] getTrimmedPreppedProjectPCsFor(boolean[] toExtract, int numComponents) {
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

}
