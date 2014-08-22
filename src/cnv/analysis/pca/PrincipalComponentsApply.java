package cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

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
 * Class to apply Loadings of a PCA to samples given as input, currently only applies LRR loadings (setup for a single float[] for each marker)
 * 
 * 
 */
public class PrincipalComponentsApply {
	private static final String[] SINGULAR_HEADER = { "PC", "Singular Value" };
	private static final String[] LOADING_FIRST = { "markerName" };
	private Project proj;
	private MarkerLoadings[] markerLoadings;
	private SingularValues singularValues;
	private int numComponents;
	private boolean[] samplesToUse;
	private boolean imputeMeanForNaN, recomputeLRR;
	private String[] markers;
	private String extrapolatedPCsFile;
	// sample,PC#
	private double[][] extrapolatedPCs;
	private Logger log;

	/**
	 * @param proj
	 *            current project
	 * @param numComponents
	 *            number of components to apply
	 * @param singularFile
	 *            file of singular values
	 * @param markerLoadingFile
	 *            file of marker loadings
	 * @param samplesToUse
	 *            boolean array of samples to use
	 * @param imputeMeanForNaN
	 *            impute the mean of the marker for a sample value with NaN
	 * @param log
	 */
	public PrincipalComponentsApply(Project proj, int numComponents, String singularFile, String markerLoadingFile, boolean[] samplesToUse, boolean imputeMeanForNaN, boolean recomputeLRR) {
		super();
		this.proj = proj;
		this.log = proj.getLog();
		this.singularValues = new SingularValues(proj.getProjectDir() + singularFile, numComponents, log);
		this.numComponents = numComponents;
		this.samplesToUse = samplesToUse;
		this.markerLoadings = MarkerLoadings.getLoadings(proj.getProjectDir() + markerLoadingFile, numComponents, log);
		this.imputeMeanForNaN = imputeMeanForNaN;
		this.recomputeLRR = recomputeLRR;
		getMarkers();
		initExtPCS();
	}

	/**
	 * @return the created extrapolated Pc file, so we can pass the name to other functions
	 */
	public String getExtrapolatedPCsFile() {
		return extrapolatedPCsFile;
	}

	/**
	 * Load the necessary markers and apply each marker's loadings on the fly
	 */
	public void applyLoadings() {
		if (proj.getSampleList().getSamples().length != samplesToUse.length) {
			log.reportError("Error - the boolean array of samples to use does not equal the length of the samples in the project, exiting");
			return;
		} else {
			MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
			for (int i = 0; i < markers.length; i++) {
				if (i % 1000 == 0) {
					float usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
					float freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
					float maxMemory = Runtime.getRuntime().maxMemory();
					log.report(ext.getTime() + "\tData loaded = " + Math.round(((double) i / (double) markers.length * 100.0)) + "%\tFree memory: " + Math.round(((double) freeMemory / (double) maxMemory * 100.0)) + "%");
				}
				MarkerData markerData = markerDataLoader.requestMarkerData(i);
				float[] lrrs;
				if (recomputeLRR) {
					lrrs = markerData.getRecomputedLRR_BAF(null, null, false, 1, 0, null, true, true, log)[1];
				} else {
					lrrs = markerData.getLRRs();
				}
				if (!hasNAN(lrrs)) {
					applyMarkerLoading(lrrs, i);
				} else if (imputeMeanForNaN) {
					lrrs = PrincipalComponentsCompute.imputeMeanForNaN(markers[i], lrrs, samplesToUse, log);
					applyMarkerLoading(lrrs, i);
				} else {
					log.reportError("Warning - marker " + markers[i] + " contained a NaN value in the extrapolated dataset, skipping it for extrapolation");
				}
				markerDataLoader.releaseIndex(i);
			}
			applySingularValues();
		}
	}

	/**
	 * Summarize the newly built extrapolated principal components, in the standard Genvisis PC format of FID\tIID\tPC1\tPC2...
	 */
	public void reportExtropolatedPCs(String output) {
		SampleData sampleData = proj.getSampleData(0, false);
		this.extrapolatedPCsFile = output;
		try {
			if (Files.exists(proj.getProjectDir() + output)) {
				Files.backup(output, proj.getProjectDir(), proj.getProjectDir() + proj.getProperty(Project.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(proj.getProjectDir() + output));
			String[] samples = proj.getSampleList().getSamples();
			writer.print("FID\tIID");
			for (int i = 0; i < numComponents; i++) {
				writer.print("\t" + PrincipalComponentsCompute.PC_STRING + (i + 1));
			}
			writer.println();
			int sampleIndex = 0;
			for (int i = 0; i < samples.length; i++) {
				if (samplesToUse[i]) {
					String[] samp = sampleData.lookup(samples[i]);
					writer.print(samp[1]);
					for (int k = 0; k < numComponents; k++) {
						writer.print("\t" + extrapolatedPCs[sampleIndex][k]);
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
			log.reportError("Error writing to file \"" + output + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
	}

	/**
	 * builds the extrapolated PCs for each sample/component combination by summing the loadings of marker (at each component) multiplied by the intensity value (for the marker at that sample)
	 * 
	 * @param lrrs
	 *            input data for all samples
	 * @param markerIndex
	 *            used to match the current marker data with the marker loading
	 */
	private void applyMarkerLoading(float[] lrrs, int markerIndex) {
		int sampleIndex = 0;
		for (int i = 0; i < lrrs.length; i++) {
			if (samplesToUse[i]) {
				for (int j = 0; j < numComponents; j++) {
					extrapolatedPCs[sampleIndex][j] += markerLoadings[markerIndex].getLoadings()[j] * lrrs[i];
				}
				sampleIndex++;
			}
		}
	}

	/**
	 * After the loadings have been summed by applyMarkerLoading, the final step to build the extrapolated PCs is to divide each summation by the singular value for that component
	 */
	private void applySingularValues() {
		int sampleIndex = 0;
		// for each sample
		for (int i = 0; i < extrapolatedPCs.length; i++) {
			// for each pc
			for (int j = 0; j < extrapolatedPCs[i].length; j++) {
				extrapolatedPCs[sampleIndex][j] /= singularValues.getComponentLoading(j);
			}
			sampleIndex++;
		}
	}

	/**
	 * Properly formats extrapolatedPCs to an array organized as extrapolatedPCs[sample0][components for sample0]
	 */
	private void initExtPCS() {
		int numSamples = 0;
		for (int i = 0; i < samplesToUse.length; i++) {
			if (samplesToUse[i]) {
				numSamples++;
			}
		}
		this.extrapolatedPCs = new double[numSamples][numComponents];
	}

	private static boolean hasNAN(float[] lrrs) {
		// should not happen, but who knows
		if (lrrs == null) {
			return true;
		}
		for (int i = 0; i < lrrs.length; i++) {
			if (Float.isNaN(lrrs[i])) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Extract the marker names for all marker loadings. The data loaded from these markers will be used to build the extrapolated PCs. We assume that the markers in the marker loadings file were laid down in project order...If the markerDataLoader hangs, this is a good spot to sort
	 */
	private void getMarkers() {
		ArrayList<String> al = new ArrayList<String>();
		for (int i = 0; i < markerLoadings.length; i++) {
			al.add(markerLoadings[i].getMarker());
		}
		this.markers = al.toArray(new String[al.size()]);
	}

	/**
	 * Helper class to facilitate easy handling of the marker loadings
	 * 
	 */
	public static class MarkerLoadings {
		private String marker;
		private double[] loadings;

		public MarkerLoadings(String marker, double[] loadings) {
			super();
			this.marker = marker;
			this.loadings = loadings;
		}

		public double[] getLoadings() {
			return loadings;
		}

		public String getMarker() {
			return marker;
		}

		public static MarkerLoadings[] getLoadings(String markerLoadingFile, int numComponents, Logger log) {
			ArrayList<MarkerLoadings> ml = new ArrayList<MarkerLoadings>();
			try {
				BufferedReader reader = Files.getReader(markerLoadingFile, false, true, false);
				String[] line = reader.readLine().trim().split("\t");
				int[] indices = ext.indexFactors(LOADING_FIRST, line, true, true);
				if (indices == null || indices[0] != 0) {
					log.reportError("Error - Marker Loading file  must have " + Array.toStr(LOADING_FIRST) + " in the first column");
					System.exit(1);
				}
				if ((line.length - 1) < numComponents) {
					log.reportError("Error - Can only apply " + (line.length - 1) + " marker loadings (as provided in " + markerLoadingFile + ", not enough for " + numComponents + " components");
					log.reportError("Please apply a smaller number of components, or provide a file with more loadings");
					System.exit(1);
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t");
					String marker = line[indices[0]];
					double[] loadings = new double[numComponents];
					for (int i = 1; i <= numComponents; i++) {
						try {
							loadings[i - 1] = Double.parseDouble(line[i]);
						} catch (NumberFormatException nfe) {
							log.reportError("Error - could not parse marking loading value " + line[i] + " to a double");
							System.exit(1);
						}
					}
					ml.add(new MarkerLoadings(marker, loadings));
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + markerLoadingFile + "\" not found in current directory");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + markerLoadingFile + "\"");
			}
			return ml.toArray(new MarkerLoadings[ml.size()]);
		}
	}

	/**
	 * Helper class to retrieve and store singular values
	 * 
	 */
	public static class SingularValues {
		private String singularFile;
		private double[] singularValues;
		private int numComponents;
		private Logger log;

		public SingularValues(String singularFile, int numComponents, Logger log) {
			this.singularFile = singularFile;
			this.numComponents = numComponents;
			this.log = log;
			this.singularValues = new double[numComponents];
			loadSingularValues();

		}

		public double getComponentLoading(int component) {
			if (component < singularValues.length) {
				return singularValues[component];
			} else {
				log.reportError("Error - do not have singular value for this component " + component);
				return Double.NaN;
			}
		}

		private void loadSingularValues() {
			int numSingular = 0;
			try {
				BufferedReader reader = Files.getReader(singularFile, false, true, false);
				String[] line = reader.readLine().trim().split("\t");
				int[] indices = ext.indexFactors(line, SINGULAR_HEADER, true, false);
				if (indices == null) {
					log.reportError("Error - singular value file must have header " + Array.toStr(SINGULAR_HEADER));
					return;
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t");
					try {
						if (numSingular < numComponents) {
							singularValues[numSingular] = Double.parseDouble(line[indices[1]]);
							numSingular++;
						}
					} catch (NumberFormatException nfe) {
						log.reportError("Error - could not parse singular value " + line[indices[1]]);
						System.exit(1);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + singularFile + "\" not found in current directory");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + singularFile + "\"");
			}
			if (numSingular != singularValues.length) {
				log.reportError("Error - not enough singular values were found in " + singularFile + " for " + numComponents + " components");
				log.reportError("Please select a smaller number of components, or provide a file with more singular values");
			}
		}
	}
}
