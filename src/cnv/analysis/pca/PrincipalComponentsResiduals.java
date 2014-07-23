package cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import stats.LeastSquares;
import stats.RegressionModel;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;

/**
 * <p>
 * Class to compute residuals wrt PCs and summarize median values, currently aimed at Mitochondrial copy number, but is extensible to other data
 * 
 *
 */
public class PrincipalComponentsResiduals {
	private static final String[] MT_REPORT = { "DNA", "FID", "IID", "Sex", "median_MT_LRR_raw", "median_MT_LRR_PC_residuals", "median_MT_LRR_PC_residuals_inverseTransformed" };
	private static final String[] MT_REPORT_EXT = { ".report.txt" };
	private static final String[] MT_REPORT_MARKERS_USED = { ".MedianMarkers.MarkersUsed.txt", ".MedianMarkers.RawValues.txt" };

	private String markersToAssessFile, output, residOutput;
	private String[] markersToAssess, samplesToReport, allProjSamples;
	private double[] medians, residuals, invTResiduals;
	private double[][] assessmentData, pcBasis;
	private byte[][] abGenotypesAfterFilters;
	private float gcThreshold;
	private Logger log;
	private Project proj;
	private int numSamples, numComponents;
	private Hashtable<String, Integer> samplesInPc;
	private boolean[] samplesToUse;
	private boolean printFull, homozygousOnly;

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
	 * @param output
	 * @param log
	 */
	public PrincipalComponentsResiduals(Project proj, String pcFile, String markersToAssessFile, int numComponents, boolean printFull, float gcThreshold, boolean homozygousOnly, String output, Logger log) {
		super();
		this.numComponents = numComponents;
		this.log = log;
		this.proj = proj;
		this.numSamples = 0;
		this.samplesInPc = new Hashtable<String, Integer>();
		this.markersToAssessFile = markersToAssessFile;
		this.printFull = printFull;
		this.gcThreshold = gcThreshold;
		this.homozygousOnly = homozygousOnly;
		this.output = output;
		loadPcFile(pcFile);
		parseSamplesToUse();
	}

	/**
	 * Compute the median data in project order, and then set back to pcfile order, just in case
	 */
	public void computeAssesmentDataMedians() {
		this.markersToAssess = PrincipalComponentsCompute.sortByProjectMarkers(proj, PrincipalComponentsCompute.loadToArray(proj.getProjectDir() + markersToAssessFile, log));
		getData();
		double[] projectOrderMedians = getLRRMedian();
		if (printFull) {
			printFull();
		}
		setMedianSortByPCs(projectOrderMedians);
	}

	/**
	 * Regress out the Pcs from the medians
	 */
	public void computeResiduals() {
		RegressionModel model = (RegressionModel) new LeastSquares(medians, prepPcs(pcBasis));
		if (!model.analysisFailed()) {
			this.residuals = model.getResiduals();
			log.report("" + model.getRsquare());
		} else {
			log.reportError("Error - the regression model has failed and residuals could not be computed");
			this.residuals = null;
		}
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
	 * Here we set the medians back to the same order as samples represented in the pc file
	 */
	private void setMedianSortByPCs(double[] projectOrderMedians) {
		this.medians = new double[projectOrderMedians.length];
		int sampleIndex = 0;
		for (int i = 0; i < allProjSamples.length; i++) {
			if (samplesToUse[i]) {
				medians[samplesInPc.get(allProjSamples[i])] = projectOrderMedians[sampleIndex];
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
	 * The residual file that the results were summarized to
	 */
	public String getResidOutput() {
		return residOutput;
	}

	/**
	 * Load each of the markers to assess (LRR, and AbGenotypesAfterFilters); We load genotypes to allow filtering by homozygous only, gc threshold, etc... Data is loaded to assessmentData organized assessmentData[marker0][data for samples in PCs], ditto for abGenotypesAfterFilters Note that the order of assessmentData[marker0] does not neccesarily reflect the same order as the samples in the pc file
	 */
	public void getData() {
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markersToAssess, new Logger());
		int count = numUsed(samplesToUse);
		this.assessmentData = new double[markersToAssess.length][count];
		this.abGenotypesAfterFilters = new byte[markersToAssess.length][count];
		float[] lrrs;
		byte[] abGenos;
		ClusterFilterCollection cluster;

		if (Files.exists(proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME))) {
			cluster = ClusterFilterCollection.load(proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME), false);
		} else {
			cluster = new ClusterFilterCollection();
			log.report("Info - did not find the cluster filter file " + proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME) + " using original genotypes");
		}

		for (int i = 0; i < markersToAssess.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			lrrs = markerData.getLRRs();
			abGenos = markerData.getAbGenotypesAfterFilters(cluster, markersToAssess[i], gcThreshold);
			int sampleIndex = 0;
			for (int k = 0; k < samplesToUse.length; k++) {
				if (samplesToUse[k]) {
					assessmentData[i][sampleIndex] = lrrs[k];
					abGenotypesAfterFilters[i][sampleIndex] = abGenos[k];
					sampleIndex++;
				}
			}
			markerDataLoader.releaseIndex(i);
		}
	}

	/**
	 * Compute the median value (after filtering) across the markers to assess for each sample represented in the principal component file Note, medians are computed in a project order manner. A further step may be required if the Pc file is not in the same sample order as the project.
	 */
	private double[] getLRRMedian() {
		double[] medians = new double[assessmentData[0].length];
		// for sample
		for (int i = 0; i < assessmentData[0].length; i++) {
			// for marker
			ArrayList<Double> sampLRR = new ArrayList<Double>();
			for (int k = 0; k < assessmentData.length; k++) {
				// test for null in case we dropped data, test for individual NaN, test for homozygous, test for missing if gcThreshold greater than 0
				if (useMarker(i, k)) {
					sampLRR.add(assessmentData[k][i]);
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
		if (assessmentData[markerIndex] != null && !Double.isNaN(assessmentData[markerIndex][sampIndex]) && (gcThreshold == 0 || abGenotypesAfterFilters[markerIndex][sampIndex] >= 0) && (abGenotypesAfterFilters[markerIndex][sampIndex] == (byte) 0 || abGenotypesAfterFilters[markerIndex][sampIndex] == (byte) 2 || !homozygousOnly)) {
			good = true;
		}
		return good;
	}

	/**
	 * Prints two files, one corresponding to the raw data, and another boolean file. We do not sort the assessmentData by the order of principal components in the pcfile, so the reporting here may be in a different order.
	 */
	private void printFull() {
		String[] files = { output + MT_REPORT_MARKERS_USED[0], output + MT_REPORT_MARKERS_USED[1] };
		PrintWriter[] writers = getNWriters(proj, files, log);
		String[] projOrderedSubset = getUsedSubset(allProjSamples, samplesToUse, log);
		if (projOrderedSubset == null) {
			log.reportError("Error - could not print full data, please remove flag or contact jlanej@gmail.com");
			return;
		} else {
			writers[0].println(MT_REPORT[0] + "\t" + Array.toStr(markersToAssess));
			writers[1].println(MT_REPORT[0] + "\t" + Array.toStr(markersToAssess));
			for (int i = 0; i < assessmentData[0].length; i++) {
				writers[0].print(projOrderedSubset[i]);
				writers[1].print(projOrderedSubset[i]);
				// for marker
				for (int k = 0; k < assessmentData.length; k++) {
					if (useMarker(i, k)) {
						writers[0].print("\tTRUE");
						writers[1].print("\t" + assessmentData[k][i]);
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
		if (samples.length != samples.length) {
			log.reportError("Error - mismatched number of samples when extracting samples used ,this should not happen");
			return null;
		}
		for (int i = 0; i < samplesToUse.length; i++) {
			if (samplesToUse[i]) {
				used.add(samples[i]);
			}
		}
		return used.toArray(new String[used.size()]);
	}

	private static PrintWriter[] getNWriters(Project proj, String[] fileOuts, Logger log) {
		PrintWriter[] writers = new PrintWriter[fileOuts.length];
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
				if (Double.isNaN(medians[i])) {
					writer.print("\tNA\tNA\tNA");
				} else {
					writer.print("\t" + medians[i] + "\t" + (residuals == null ? "NA" : residuals[modelCount]) + "\t" + (invTResiduals == null ? "NA" : invTResiduals[modelCount]));
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
	 * Loads a principal component file to pcBasis[][] organized as pcBasis[Basis0][Basis for samples]. As the components are loaded, the corresponding sample (DNA) is tracked by the hashtable samplesInPc SamplestoReport maintains the order of samples represented by the pcFile
	 * 
	 * @param pcFile
	 */
	public void loadPcFile(String pcFile) {
		String pcFilefull = proj.getProjectDir() + pcFile;
		SampleData sampleData = proj.getSampleData(0, false);
		ArrayList<String> pcSamps = new ArrayList<String>();
		int sampIndex = 0;
		try {
			BufferedReader reader = Files.getReader(pcFilefull, proj.getJarStatus(), true, false);
			String[] line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equals("FID") || !line[1].equals("IID")) {
				log.reportError("Error - different format than expected; first column should be FID and second column should be IID, followed by PCs");
				return;
			}
			if ((line.length - 2) < numComponents) {
				log.reportError("Error - cannot use " + numComponents + " components when only " + (line.length - 2) + " are provided");
				return;
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (sampleData.lookup(line[0] + "\t" + line[1]) != null) {
					samplesInPc.put(sampleData.lookup(line[0] + "\t" + line[1])[0], numSamples);
					pcSamps.add(line[0] + "\t" + line[1]);
					numSamples++;
				} else {
					log.reportError("Error - could not find " + line[0] + "\t" + line[1] + " in sample Data");
					System.exit(1);
				}
			}
			reader.close();
			this.pcBasis = new double[numComponents][numSamples];
			reader = Files.getReader(pcFilefull, proj.getJarStatus(), true, false);
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
		int count = 0;
		for (int i = 0; i < samplesToUse.length; i++) {
			if (samplesToUse[i]) {
				count++;
			}
		}
		return count;
	}

}
