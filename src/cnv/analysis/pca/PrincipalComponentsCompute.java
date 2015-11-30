package cnv.analysis.pca;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.manage.MDL;
import cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import cnv.var.SampleData;
import common.Aliases;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import ejml.CommonOps;
import ejml.DecompositionFactory;
import ejml.DenseMatrix64F;
import ejml.SingularOps;
import ejml.SingularValueDecomposition;

/**
 * <p>
 * Computes basis (principal components) of a dataset , adapted from ejml
 * 
 * @author Peter Abeles , John Lane
 */
public class PrincipalComponentsCompute {
	public static final String[] OUTPUT_EXT = { ".PCs.txt", ".PCs.MarkerReport.txt", ".PCs.MarkerLoadings.txt", ".PCs.SingularValues.txt", ".Summary.txt" };
	public static final String PC_STRING = "PC";
	public static final String[] SAMPLE = { "sample", "FID", "IID" };
	public static final String MARKER = "markerName";
	private static final String[] MARKER_REPORT_SMALL = { "Marker", "Used for PCA (Contains at least one value that is not NaN)" };
	private static final String SV_STRING = "Singular Value";
	// principle component subspace is stored in the rows
	private DenseMatrix64F V_t;

	// how many principle components are used
	private int numComponents;
	// where the data is stored
	private DenseMatrix64F A = new DenseMatrix64F(1, 1);
	private int sampleIndex;
	// mean values of each element across all the samples, if center is flagged
	private double mean[];
	private double[] singularValues;
	private boolean[] samplesToUse;
	private String singularValuesFile;
	private String pcFile;
	private String markerLoadingFile;

	public String getSingularValuesFile() {
		return singularValuesFile;
	}

	public void setSingularValuesFile(String singularValuesFile) {
		this.singularValuesFile = singularValuesFile;
	}

	public String getPcFile() {
		return pcFile;
	}

	public void setPcFile(String pcFile) {
		this.pcFile = pcFile;
	}

	public String getMarkerLoadingFile() {
		return markerLoadingFile;
	}

	public void setMarkerLoadingFile(String markerLoadingFile) {
		this.markerLoadingFile = markerLoadingFile;
	}

	public PrincipalComponentsCompute() {
		this.A = new DenseMatrix64F(1, 1);
	}

	public boolean[] getSamplesToUse() {
		return samplesToUse;
	}

	public double getSingularValueAt(int index) {
		return singularValues[index];
	}

	public double[] getSingularValues() {
		return singularValues;
	}

	public void setup(int numSamples, int sampleSize) {
		mean = new double[sampleSize];
		A.reshape(numSamples, sampleSize, false);
		sampleIndex = 0;
		numComponents = -1;
	}

	public void addSample(double[] sampleData) {
		if (A.getNumCols() != sampleData.length)
			throw new IllegalArgumentException("Unexpected sample size");
		if (sampleIndex >= A.getNumRows())
			throw new IllegalArgumentException("Too many samples");

		for (int i = 0; i < sampleData.length; i++) {
			A.set(sampleIndex, i, sampleData[i]);
		}
		sampleIndex++;
	}

	public void computeBasis(int numComponents, boolean center) {
		if (numComponents > A.getNumCols())
			throw new IllegalArgumentException("More components requested than the data's length. Have " + A.getNumCols() + " available and " + numComponents + " were requested");
		if (sampleIndex != A.getNumRows())
			throw new IllegalArgumentException("Not all the data has been added");
		if (numComponents > sampleIndex)
			throw new IllegalArgumentException("More data needed to compute the desired number of components");

		this.numComponents = numComponents;
		// compute the mean of all the samples if center is flagged
		if (center) {
			for (int i = 0; i < A.getNumRows(); i++) {
				for (int j = 0; j < mean.length; j++) {
					mean[j] += A.get(i, j);
				}
			}
			for (int j = 0; j < mean.length; j++) {
				mean[j] /= A.getNumRows();
			}

			// subtract the mean from the original data
			for (int i = 0; i < A.getNumRows(); i++) {
				for (int j = 0; j < mean.length; j++) {
					A.set(i, j, A.get(i, j) - mean[j]);
				}
			}
		}

		// Compute SVD and save time by not computing U,
		// if loadings are wanted, we compute them later/on the fly to save memory
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(A.numRows, A.numCols, false, true, false);
		if (!svd.decompose(A))
			throw new RuntimeException("SVD failed");
		V_t = svd.getV(null, true);
		DenseMatrix64F W = svd.getW(null);
		// Singular values are in an arbitrary order initially
		SingularOps.descendingOrder(null, false, W, V_t, true);

		// strip off unneeded components and find the basis
		V_t.reshape(numComponents, mean.length, true);

		this.singularValues = Array.subArray(getDiagonal(W), 0, numComponents);// W is the sorted (descending order on the diagonal) matrix of singular values

	}

	public double[] getDiagonal(DenseMatrix64F W) {
		int numSingular = Math.min(W.numRows, W.numCols);
		double[] singularValues = new double[numSingular];
		for (int i = 0; i < numSingular; i++) {
			singularValues[i] = W.get(i, i);
		}
		return singularValues;
	}

	/**
	 * Extract a particular basis vector (principal component)
	 * 
	 * @param which
	 *            basis to extract
	 * @return
	 */
	public double[] getBasisVector(int which) {
		if (which < 0 || which >= numComponents)
			throw new IllegalArgumentException("Invalid component");

		DenseMatrix64F v = new DenseMatrix64F(1, A.numCols);
		CommonOps.extract(V_t, which, which + 1, 0, A.numCols, v, 0, 0);

		return v.data;
	}

	// main access point for project based analysis
	/**
	 * @param proj
	 * @param excludeSamples
	 *            exclude samples as defined in sample data
	 * @param numComponents
	 *            number of components to compute
	 * @param printFullData
	 *            print the full data that is loaded (use only for testing/to verify)
	 * @param center
	 *            center the marker about its mean
	 * @param reportMarkerLoadings
	 *            report loadings
	 * @param reportSingularValues
	 *            report singular values
	 * @param imputeMeanForNaN
	 *            impute the marker to a sample when it has a NaN value
	 * @param recomputeLRR
	 *            recompute Log R Ratios on the fly
	 * @param useFile
	 *            an optional file of samples to use
	 * @param output
	 *            the base name
	 * @param log
	 * @return the PrincipalComponentsCompute object with PCs computed
	 */
	public static PrincipalComponentsCompute getPrincipalComponents(Project proj, boolean excludeSamples, int numComponents, boolean printFullData, boolean center, boolean reportMarkerLoadings, boolean reportSingularValues, boolean imputeMeanForNaN, boolean recomputeLRR, String useFile, String output) {
		Logger log = proj.getLog();
		PrincipalComponentsCompute pcs = populateWithExistingFiles(proj, output, log);
		if (pcs != null) {
			return pcs; // we found all three necessary files, and so we just use them
		}

		boolean[] samplesToUse = getSamples(proj, excludeSamples, useFile);
		String[] markers = getMarkers(proj);
		int numSamples = Array.booleanArraySum(samplesToUse);

		if (numComponents > numSamples) {
			log.reportError("Error - cannot request more principal components (n=" + numComponents + ") than there are valid samples (n=" + numSamples + ")");
			return null;
		}

		if (numComponents > markers.length) {
			log.reportError("Error - cannot request more principal components (n=" + numComponents + ") than there are markers (n=" + markers.length + ")");
			return null;
		}

		// deals with NaN on the fly
		double[][] dataToUse = getData(proj, markers, samplesToUse, printFullData, imputeMeanForNaN, true, recomputeLRR, output);
		pcs = getPrincipalComponents(numComponents, center, dataToUse, true, log);
		double[][] pcsBasis = getPCs(pcs, numComponents, true, log);
		reportPCs(proj, pcs, numComponents, output, samplesToUse, pcsBasis);
		if (reportMarkerLoadings) {
			reportLoadings(proj, pcs, dataToUse, pcsBasis, markers, output);
		}
		if (reportSingularValues) {
			reportSingularValues(proj, pcs, output);
		}

		return pcs;
	}

	/**
	 * @param numComponents
	 *            number of components to compute
	 * @param center
	 *            center about the mean
	 * @param dataToUse
	 *            a double[][], organized as dataToUse[marker0][dataforMarkerAcrossSamples] (marker dominant, not sample dominant)
	 * @param log
	 * @return
	 */
	public static PrincipalComponentsCompute getPrincipalComponents(int numComponents, boolean center, double[][] dataToUse, boolean verbose, Logger log) {
		// count valid input (is not null )
		int initpc = getUseCount(dataToUse);
		PrincipalComponentsCompute pcs = new PrincipalComponentsCompute();
		int numSamples = 0;
		if (initpc > 0) {
			for (int i = 0; i < dataToUse.length; i++) {
				if (dataToUse[i] != null) {
					// initialize underlying arrays
					pcs.setup(initpc, dataToUse[i].length);
					numSamples = dataToUse[i].length;
				}
			}
			if (verbose) {
				log.report(ext.getTime() + " Using " + initpc + (initpc > 1 ? " inputs" : " input") + " with valid (not NaN) data for SVD across " + numSamples + " samples");
			}
			// adding on a per marker basis
			for (int i = 0; i < dataToUse.length; i++) {
				if (dataToUse[i] != null) {
					try {
						pcs.addSample(dataToUse[i]);
					} catch (IllegalArgumentException iae) {
						log.reportError("Error - could not add all data for SVD");
						log.reportException(iae);
					} catch (ArrayIndexOutOfBoundsException aioe) {
						log.reportError("Error - matrix for SVD ran out of space, the maximum number of markers can be " + (int) (Integer.MAX_VALUE / numSamples));
						log.reportException(aioe);
					}
				}
			}
			if (verbose) {
				log.report(ext.getTime() + " Added all valid input data to matrix, computing SVD");
			}
			try {
				pcs.computeBasis(numComponents, center);
			} catch (IllegalArgumentException iae) {
				log.reportError("Error - the number of components must be less than the number of markers used AND less than the number of samples used");
				log.reportError("Error - Please select a smaller number of principal components to compute, or add more markers or samples");
				log.reportException(iae);
			}
		} else {
			log.reportError("Error - no valid data was found");
		}
		return pcs;
	}

	/**
	 * @param proj
	 * @param excludeSamples
	 *            determine whether to exclude samples as defined in sample data
	 * @param useFile
	 *            if null, we get project samples, based on excludeSamples status
	 * @param log
	 * @return boolean[] of all samples in the project
	 */
	public static boolean[] getSamples(Project proj, boolean excludeSamples, String useFile) {
		if (useFile == null) {
			return getSamples(proj, excludeSamples);
		}
		if (!Files.exists(proj.PROJECT_DIRECTORY.getValue() + useFile)) {
			proj.getLog().reportError("Error - could not find the file of samples to use " + proj.PROJECT_DIRECTORY.getValue() + useFile);
			return null;
		} else {
			return getSamplesFromFile(proj, proj.PROJECT_DIRECTORY.getValue() + useFile);
		}
	}

	/**
	 * If the pcFile, singular value file, and marker loading file all exist, we return a {@link PrincipalComponentsCompute} with these files populated.
	 * <p>
	 * These three files are generally all that are needed. If they do not all exist, we return {@link null} and continue with the computation TODO, perhaps check the files for proper number of components and identical samples
	 */
	private static PrincipalComponentsCompute populateWithExistingFiles(Project proj, String output, Logger log) {
		PrincipalComponentsCompute pcs;
		String pcFile, singularValuesFile, markerLoadingFile;
		pcFile = ext.rootOf(output) + OUTPUT_EXT[0];
		singularValuesFile = ext.rootOf(output) + OUTPUT_EXT[3];
		markerLoadingFile = ext.rootOf(output) + OUTPUT_EXT[2];
		if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + pcFile) && Files.exists(proj.PROJECT_DIRECTORY.getValue() + singularValuesFile) && Files.exists(proj.PROJECT_DIRECTORY.getValue() + markerLoadingFile)) {
			log.report("Detected that the following principal component files already exist:\n" + proj.PROJECT_DIRECTORY.getValue() + pcFile + "\n" + proj.PROJECT_DIRECTORY.getValue() + singularValuesFile + "\n" + proj.PROJECT_DIRECTORY.getValue() + markerLoadingFile + "\n");
			log.report("Skipping the principal component computation and using these files instead.");
			log.report("If this is incorrect (using a different number of components, new samples, etc...),  please remove or change the name of the files listed above.\n Alternatively, specify a new analysis name");
			pcs = new PrincipalComponentsCompute();
			pcs.setPcFile(pcFile);
			pcs.setSingularValuesFile(singularValuesFile);
			pcs.setMarkerLoadingFile(markerLoadingFile);

		} else {
			pcs = null;
		}
		return pcs;
	}

	private static boolean[] getSamples(Project proj, boolean excludeSamples) {
		SampleList sampleList = proj.getSampleList();
		String[] samples = sampleList.getSamples();
		boolean[] samplesToUse = new boolean[samples.length];
		Logger log = proj.getLog();

		if (!proj.getSampleData(1, false).hasExcludedIndividuals() && excludeSamples) {
			log.reportError("Error - cannot exclude individuals for PCA , no factor named 'Exclude/CLASS=Exclude' in Sample Data");
			return null;
		} else if (excludeSamples) {
			SampleData sampleData = proj.getSampleData(1, false);
			int use = 0;
			for (int i = 0; i < samples.length; i++) {
				samplesToUse[i] = !sampleData.individualShouldBeExcluded(samples[i]);
				if (samplesToUse[i]) {
					use++;
				}
			}
			log.report("Computing PCs using " + use + " samples");
		} else {
			log.report("Computing PCs using all " + samples.length + " samples");
			Arrays.fill(samplesToUse, true);
		}

		return samplesToUse;
	}

	// can return null
	private static boolean[] getSamplesFromFile(Project proj, String useFile) {
		SampleData sampleData = proj.getSampleData(0, false);
		Logger log = proj.getLog();
		String[] samplesToUseFromFile = HashVec.loadFileToStringArray(useFile, false, false, new int[] { 0 }, true, true, "\t");
		//previous method causes issues with spaces in sample names
		// HashVec.loadFileToStringArray(useFile, false, new int[] { 0 }, true);
		String[] projSamples = proj.getSampleList().getSamples();
		boolean[] samplesToUse = new boolean[projSamples.length];
		Hashtable<String, Boolean> track = new Hashtable<String, Boolean>();
		int used = 0;

		for (int i = 0; i < samplesToUseFromFile.length; i++) {
			if (sampleData.lookup(samplesToUseFromFile[i]) != null) {
				track.put(sampleData.lookup(samplesToUseFromFile[i])[0], true);
			} else {
				log.reportError("Error -could not find sample " + samplesToUseFromFile[i] + " in sample data ");
				return null;
			}
		}
		for (int i = 0; i < projSamples.length; i++) {
			if (track.containsKey(projSamples[i])) {
				samplesToUse[i] = true;
				used++;
			} else {
				samplesToUse[i] = false;
			}
		}
		if (used != samplesToUseFromFile.length) {
			log.reportError("Error - " + used + " " + (used > 1 ? "samples were " : "sample was ") + " not found in the sample data file " + proj.getProperty(proj.SAMPLE_DATA_FILENAME));
			return null;
		}
		log.report(ext.getTime() + " Using the " + used + " samples in the project that passed QC " + (useFile != null ? " and that were also in the useFile" : ""));

		return samplesToUse;
	}

	// this version only supports target markers, and it takes them all
	private static String[] getMarkers(Project proj) {
		String[] markers = proj.getTargetMarkers();
		return sortByProjectMarkers(proj, markers);
	}

	/**
	 * This is necessary to allow marker data loader to load markers in order. Otherwise it can easily reach the read ahead limit, and parent thread cannot release
	 * 
	 * @param proj
	 * @param markers
	 * @return
	 */
	public static String[] sortByProjectMarkers(Project proj, String[] markers) {
		String[] projectMarkers, sorted;
		HashSet<String> tracker;
		int index;

		projectMarkers = proj.getMarkerNames();
		tracker = new HashSet<String>();

		for (int i = 0; i < markers.length; i++) {
			if (ext.indexOfStr(markers[i], Aliases.MARKER_NAMES) == -1) {
				tracker.add(markers[i]);
			}
		}

		index = 0;
		sorted = new String[tracker.size()];
		for (int i = 0; i < projectMarkers.length; i++) {
			if (tracker.contains(projectMarkers[i])) {
				sorted[index] = projectMarkers[i];
				tracker.remove(projectMarkers[i]);
				index++;
			}
		}

		return sorted;
	}

	/**
	 * Set up a matrix for the number of markers and samples
	 * 
	 * @param numMarkers
	 * @param samplesToUse
	 * @return
	 */
	private static double[][] getAppropriateArray(int numMarkers, boolean[] samplesToUse) {
		int use = 0;
		for (int i = 0; i < samplesToUse.length; i++) {
			if (samplesToUse[i]) {
				use++;
			}
		}
		return new double[numMarkers][use];
	}

	/**
	 * @param proj
	 *            current project
	 * @param markers
	 *            markers to load
	 * @param samplesToUse
	 *            samples to use
	 * @param printFullData
	 *            print all data loaded
	 * @param imputeMeanForNaN
	 *            if a sample has NaN, impute it with the mean of the current marker
	 * @param recomputeLRR
	 *            recompute Log R Ratios on the fly
	 * @param dealWithNaN
	 * @param output
	 * @param log
	 * @return
	 */
	public static double[][] getData(Project proj, String[] markers, boolean[] samplesToUse, boolean printFullData, boolean imputeMeanForNaN, boolean dealWithNaN, boolean recomputeLRR, String output) {
		double[][] dataToUse = getAppropriateArray(markers.length, samplesToUse);
		MDL mdl = new MDL(proj, markers, 2, 1);
		String gcCorrections = proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue()[0];
		proj.getLog().reportTimeInfo("Using gc paramter file "+gcCorrections);
		//TODO, arg
		GcAdjustorParameters parameters = GcAdjustorParameters.readSerial(gcCorrections, proj.getLog());
		// MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		boolean[] markerUsed = new boolean[markers.length];
		Logger log = proj.getLog();

		Arrays.fill(markerUsed, true);
		// for (int i = 0; i < markers.length; i++) {
		int i = 0;
		Hashtable<String , Integer> projectIndices =proj.getMarkerIndices();
		while (mdl.hasNext()) {
			float[][] data = new float[1][dataToUse[0].length];
			if (i % 1000 == 0) {
				float usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
				float freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
				float maxMemory = Runtime.getRuntime().maxMemory();
				log.report(ext.getTime() + "\tData loaded = " + Math.round(((double) i / (double) markers.length * 100.0)) + "%\tFree memory: " + Math.round(((double) freeMemory / (double) maxMemory * 100.0)) + "%");
			}
			// MarkerData markerData = markerDataLoader.requestMarkerData(i);
			MarkerData markerData = mdl.next();
			if (recomputeLRR) {
				// Warning - assuming only autosomal, not providing sex. Not caring about gc Threshold or call rate either
				// TODO Recompute only from samplesToUse??
				data[0] = markerData.getRecomputedLRR_BAF(null, null, false, 1, 0, null, true, true, log)[1];
			} else {
				data[0] = markerData.getLRRs();
			}
			if (parameters != null) {
				if (parameters.getCentroids() == null) {
					proj.getLog().reportTimeError("NO centroids available in load");
					return null;
				}
				data[0] = markerData.getGCCorrectedLRRBAF(parameters, projectIndices.get(markerData.getMarkerName()), proj.getLog())[1];
			}
			// please hard code dealWithNaN as it determines whether later data checks must be done;
			if (dealWithNaN && hasNAN(data)) {
				markerUsed[i] = false;
				if (imputeMeanForNaN) {
					data[0] = imputeMeanForNaN(markers[i], data[0], samplesToUse, log);
					markerUsed[i] = true;
				} else {
					data[0] = null;
				}
			}
			int sampleIndex = 0;
			// only use if have some data
			if (data[0] != null) {
				for (int k = 0; k < samplesToUse.length; k++) {
					if (samplesToUse[k]) {
						dataToUse[i][sampleIndex] = data[0][k];
						sampleIndex++;
					}
				}
			} else {
				// data is not valid,handle nulls when intitializing PCs
				dataToUse[i] = null;
				markerUsed[i] = false;
			}
			// markerDataLoader.releaseIndex(i);
			i++;
		}

		// markerDataLoader.reportWaitTimes();

		reportMarkersUsed(proj, markers, markerUsed, dataToUse, printFullData, samplesToUse, output);
		return dataToUse;
	}

	public static float[] imputeMeanForNaN(String markerName, float[] data, boolean[] samplesToUse, Logger log) {
		float sum = 0;
		int count = 0;
		float mean = 0;
		ArrayList<Integer> toImpute = new ArrayList<Integer>();
		for (int i = 0; i < data.length; i++) {
			if (Float.isNaN(data[i])) {
				toImpute.add(i);
			} else if (samplesToUse[i]) {
				sum += data[i];
				count += 1;
			}
		}
		if (count > 0) {
			mean = sum / count;
			for (int i = 0; i < toImpute.size(); i++) {
				data[toImpute.get(i)] = mean;
			}
		} else {
			log.reportError("Error - could not impute mean values for marker " + markerName + " (all values were NaN), skipping this marker...");
			data = null;
		}
		return data;
	}

	private static boolean hasNAN(float[][] data) {
		for (int i = 0; i < data.length; i++) {
			if (data[i] != null) {
				for (int j = 0; j < data[i].length; j++) {
					if (Float.isNaN(data[i][j])) {
						return true;
					}
				}
			}
		}
		return false;
	}

	private static int getUseCount(double[][] dataToUse) {
		int use = 0;
		for (int i = 0; i < dataToUse.length; i++) {
			if (dataToUse[i] != null) {
				use++;
			}
		}
		return use;
	}

	private static void reportMarkersUsed(Project proj, String[] markers, boolean[] markerUsed, double[][] dataToUse, boolean printFullData, boolean[] samplesToUse, String output) {
		String markersUsedForPCA = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(output) + OUTPUT_EXT[1];
		Logger log = proj.getLog();

		try {
			if (Files.exists(markersUsedForPCA)) {
				Files.backup(ext.rootOf(output) + OUTPUT_EXT[1], proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(markersUsedForPCA));
			if (!printFullData) {
				writer.println(Array.toStr(MARKER_REPORT_SMALL));
				for (int i = 0; i < markers.length; i++) {
					writer.println(markers[i] + "\t" + markerUsed[i]);
				}
			} else {
				printFullData(proj, writer, markers, markerUsed, dataToUse, samplesToUse);
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + markersUsedForPCA + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + markersUsedForPCA + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
	}

	private static void printFullData(Project proj, PrintWriter writer, String[] markers, boolean[] markerUsed, double[][] dataToUse, boolean[] samplesToUse) {
		String[] samples = proj.getSamples();
		writer.print(SAMPLE[0]);
		printSampleHeader(writer, samplesToUse, samples);
		for (int i = 0; i < markers.length; i++) {
			if (dataToUse[i] != null) {
				writer.print(markers[i]);
				printData(writer, dataToUse, i);
				writer.println();
			}
		}
	}

	private static void printData(PrintWriter writer, double[][] dataToUse, int markerIndex) {
		for (int j = 0; j < dataToUse[markerIndex].length; j++) {
			writer.print("\t" + dataToUse[markerIndex][j]);
		}
	}

	private static void printSampleHeader(PrintWriter writer, boolean[] samplesToUse, String[] samples) {
		for (int i = 0; i < samples.length; i++) {
			if (samplesToUse[i]) {
				writer.print("\t" + samples[i]);
			}
		}
		writer.println();
	}

	public static void reportSingularValues(Project proj, PrincipalComponentsCompute pcs, String output) {
		Logger log = proj.getLog();

		output = ext.rootOf(output) + OUTPUT_EXT[3];
		pcs.setSingularValuesFile(output);
		try {
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + output)) {
				Files.backup(output, proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output));
			double[] singularValues = pcs.getSingularValues();
			writer.println(PC_STRING + "\t" + SV_STRING);
			for (int i = 0; i < singularValues.length; i++) {
				writer.println(PC_STRING + (i + 1) + "\t" + singularValues[i]);
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + output + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + output + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
	}

	// sum of marker intensity * basis per sample, divided by singular values -computed on the fly to save memory.
	/**
	 * Report the marker loadings for the corresponding components, instead of storing loadings in memory, we track the related file
	 */
	private static void reportLoadings(Project proj, PrincipalComponentsCompute pcs, double[][] dataToUse, double[][] pcsBasis, String[] markers, String output) {
		Logger log = proj.getLog();

		output = ext.rootOf(output) + OUTPUT_EXT[2];
		pcs.setMarkerLoadingFile(output);
		try {
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + output)) {
				Files.backup(output, proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output));
			writer.print(MARKER);
			for (int i = 0; i < pcsBasis.length; i++) {
				writer.print("\t" + PC_STRING + (i + 1));
			}
			writer.println();
			for (int i = 0; i < dataToUse.length; i++) {
				if (dataToUse[i] != null) {
					writer.print(markers[i]);
					for (int j = 0; j < pcsBasis.length; j++) {
						double singularValue = pcs.getSingularValueAt(j);
						writer.print("\t" + getMarkerLoading(singularValue, dataToUse[i], pcsBasis[j]));
					}
					writer.println();
				}
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + output + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + output + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
	}

	/**
	 * @param singularValue
	 *            the corresponding singular value for the component
	 * @param data
	 *            the intensity data (across marker)
	 * @param basis
	 *            the basis for the component (across samples)
	 * @return sum of marker intensities * basis per sample divided by singular value
	 */
	private static double getMarkerLoading(double singularValue, double[] data, double[] basis) {
		double sum = 0;
		for (int i = 0; i < data.length; i++) {
			sum += data[i] * basis[i];
		}
		return sum / singularValue;
	}

	/**
	 * Extract a given number of principal components to a double[][] organized as double[pc0][pcs for samples]
	 * 
	 * @param pcs
	 *            PrincipalComponentsCompute object with computed basis vectors
	 * @param numComponents
	 *            number of components to extract
	 * @param log
	 * @return extracted components
	 */
	public static double[][] getPCs(PrincipalComponentsCompute pcs, int numComponents, boolean verbose, Logger log) {
		double[][] pcsBasis = new double[numComponents][];
		if (verbose) {
			log.report(ext.getTime() + " Extracting " + numComponents + " principle components");
		}
		for (int i = 0; i < numComponents; i++) {
			try {
				pcsBasis[i] = pcs.getBasisVector(i);
			} catch (IllegalArgumentException iae) {
				log.reportError("Error - retrieving component " + (i + 1));
				log.reportException(iae);
				System.exit(1);
			}
		}
		return pcsBasis;
	}

	private static double[][] reportPCs(Project proj, PrincipalComponentsCompute pcs, int numComponents, String output, boolean[] samplesToUse, double[][] pcsBasis) {
		SampleData sampleData = proj.getSampleData(0, false);
		Logger log = proj.getLog();

		output = ext.rootOf(output) + OUTPUT_EXT[0];
		pcs.setPcFile(output);
		try {
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + output)) {
				Files.backup(output, proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			log.report(ext.getTime() + " Free memory: " + Math.round(((double) Runtime.getRuntime().freeMemory() / (double) Runtime.getRuntime().totalMemory() * 100.0)) + "%");
			log.report(ext.getTime() + " Reporting top " + numComponents + " PCs");
			PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output));
			String[] samples = proj.getSampleList().getSamples();
			writer.print(SAMPLE[1] + "\t" + SAMPLE[2]);
			for (int i = 0; i < numComponents; i++) {
				writer.print("\t" + PC_STRING + (i + 1));
			}
			writer.println();
			int sampleIndex = 0;
			for (int i = 0; i < samples.length; i++) {
				if (samplesToUse[i]) {
					String[] samp = sampleData.lookup(samples[i]);
					writer.print(samp[1]);
					for (int k = 0; k < numComponents; k++) {
						writer.print("\t" + pcsBasis[k][sampleIndex]);
					}
					sampleIndex++;
					writer.println();
				}
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + output + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + output + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
		return pcsBasis;
	}
}
