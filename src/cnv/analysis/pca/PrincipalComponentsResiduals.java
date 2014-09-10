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
	private double[][] assessmentData, pcBasis;
	private byte[][] abGenotypesAfterFilters;
	private float gcThreshold;
	protected Logger log;
	protected Project proj;
	protected int numComponents, totalNumComponents;
	private int numSamples;// numComponents are loaded, totalNumComponents are present in the pc file
	protected Hashtable<String, Integer> samplesInPc;// stores the index of a sample in the pc file
	protected boolean[] samplesToUse;// corresponds to the samples in the PC
	private boolean printFull, homozygousOnly, recomputeLRR;

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
		loadPcFile(pcFile);
		parseSamplesToUse();
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
	 * The residual file that the results were summarized to
	 */
	public String getResidOutput() {
		return residOutput;
	}

	/**
	 * Load each of the markers to assess (LRR, and AbGenotypesAfterFilters); We load genotypes to allow filtering by homozygous only, gc threshold, etc...
	 * <p>
	 * Data is loaded to assessmentData organized assessmentData[marker0][data for samples in PCs], ditto for abGenotypesAfterFilters
	 * <p>
	 * Note that the order of assessmentData[marker0] does not neccesarily reflect the same order as the samples in the pc file
	 * 
	 */
	public void getData() {
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markersToAssess);
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
		PrintWriter[] writers = getNWriters(proj, files);
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
	 * Loads a principal component file to pcBasis[][] organized as pcBasis[Basis0][Basis for samples]. As the components are loaded, the corresponding sample (DNA) is tracked by the hashtable samplesInPc SamplestoReport maintains the order of samples represented by the pcFile
	 * 
	 * @param pcFile
	 */
	public void loadPcFile(String pcFile) {
		String pcFilefull;
		if (Files.exists(pcFile)) {
			pcFilefull = pcFile;
		} else {
			pcFilefull = proj.getProjectDir() + pcFile;
		}
		log.report("Info - loading principal components from " + pcFilefull);
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
				log.reportError("Warning - cannot use " + numComponents + " components when only " + (line.length - 2) + " are provided, only loading " + (line.length - 2));
				numComponents = line.length - 2;
			}
			this.totalNumComponents = line.length - 2;
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

	// The following methods are geared toward the cross validation of principal components

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
		this.log = proj.getLog();

	}

	public String[] getAllProjSamples() {
		return allProjSamples;
	}

	public Project getProj() {
		return proj;
	}

	public boolean[] getSamplesToUse() {
		return samplesToUse;
	}

	public double[][] getPreppedPCs() {
		return prepPcs(pcBasis);
	}

	public double[][] getTrimmedPreppedPCs(int numComponents) {
		return prepPcs(trimPcBasis(numComponents, pcBasis, log));
	}

	/**
	 * @param data
	 *            Currently, the data must be the same length as all samples in the project. To remove individuals from the regression, set values to NaN
	 * @param numComponents
	 *            number of components to include in the regression model
	 * @param title
	 *            optional to report if model has failed (could be a specific marker or trait, etc...)
	 * @return
	 */
	public double[] getCorrectedDataAt(double[] data, int numComponents, String title) {
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
			setAssesmentDataSortByPCs(Array.subArray(data, samplesToUse));// sorts the lrrs according to the order in the pcFile
			RegressionModel model = (RegressionModel) new LeastSquares(assesmentData, getTrimmedPreppedPCs(numComponents));
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

}
