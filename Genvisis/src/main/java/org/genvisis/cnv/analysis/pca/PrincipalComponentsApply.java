package org.genvisis.cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * <p>
 * Class to apply Loadings of a PCA to samples given as input, currently only applies LRR loadings
 * (setup for a single float[] for each marker)
 *
 *
 */
public class PrincipalComponentsApply {
	private static final String[] SINGULAR_HEADER = {"PC", "Singular Value"};
	private static final String[] LOADING_FIRST = {"markerName"};
	private final Project proj;
	private final MarkerLoadings[] markerLoadings;
	private final SingularValues singularValues;
	private final int numComponents;
	private final boolean[] samplesToUse;
	private final boolean imputeMeanForNaN, recomputeLRR;
	private String[] markers;
	private String extrapolatedPCsFile;
	// sample,PC#
	private double[][] extrapolatedPCs;
	private GcAdjustorParameters params;
	private final Logger log;

	/**
	 * @param proj current project
	 * @param numComponents number of components to apply
	 * @param singularFile file of singular values
	 * @param markerLoadingFile file of marker loadings
	 * @param samplesToUse boolean array of samples to use
	 * @param imputeMeanForNaN impute the mean of the marker for a sample value with NaN
	 * @param recomputeLRR recompute Log R Ratios on the fly
	 * @param log
	 */
	public PrincipalComponentsApply(Project proj, int numComponents, String singularFile,
																	String markerLoadingFile, boolean[] samplesToUse,
																	boolean imputeMeanForNaN, boolean recomputeLRR) {
		super();
		this.proj = proj;
		log = proj.getLog();
		singularValues = new SingularValues(proj.PROJECT_DIRECTORY.getValue()	+ singularFile,
																				numComponents, log);
		this.numComponents = numComponents;
		this.samplesToUse = samplesToUse;
		markerLoadings = MarkerLoadings.getLoadings(proj.PROJECT_DIRECTORY.getValue()
																								+ markerLoadingFile, numComponents, log);
		this.imputeMeanForNaN = imputeMeanForNaN;
		this.recomputeLRR = recomputeLRR;
		getMarkers();
		initExtPCS();
	}

	public GcAdjustorParameters getParams() {
		return params;
	}

	public void setParams(GcAdjustorParameters params) {
		this.params = params;
	}

	/**
	 * @return the created extrapolated Pc file, so we can pass the name to other functions
	 */
	public String getExtrapolatedPCsFile() {
		return extrapolatedPCsFile;
	}

	public void setExtrapolatedPCsFile(String extrapolatedPCsFile) {
		this.extrapolatedPCsFile = extrapolatedPCsFile;
	}

	/**
	 * @param output the filename to use as output
	 * @param warn report a message geared towards skipping the computation
	 * @return true if output exists
	 */
	public boolean outputExists(String output, boolean warn) {
		boolean exists = false;
		if (Files.exists(output)) {
			if (warn) {
				proj.getLog().report(
															"Detected that the following extrapolated principal component file already exists:\n"
															+ output + "\n");
				proj.getLog()
						.report("Skipping the principal component extrapolation and using this file instead");
				proj.getLog()
						.report("If this is incorrect (using a different number of components, new samples, etc...),  please remove or change the name of the file listed above.\n Alternatively, specify a new analysis name");
			}
			exists = true;
		}
		return exists;
	}

	/**
	 * Load the necessary markers and apply each marker's loadings on the fly
	 */
	public void applyLoadings() {
		if (proj.getSampleList().getSamples().length != samplesToUse.length) {
			log.reportError("Error - the boolean array of samples to use does not equal the length of the samples in the project, exiting");
			return;
		} else {

			if (params != null && recomputeLRR) {
				proj.getLog()
						.reportTimeError("recompute lrr was flagged AND gc correction parameters were passed to data load when applying PCs");
				return;
			}
			if (params != null) {
				proj.getLog().reportTimeInfo("Will be performing GC correction of input for applying PCs");
				if (params.getCentroids() != null) {
					proj.getLog().reportTimeInfo("Will be adjusting data for current centroids");
				}
			}
			Hashtable<String, Integer> projectIndices = proj.getMarkerIndices();

			MDL mdl = new MDL(proj, proj.getMarkerSet(), markers, 2, 100);
			// MarkerDataLoader markerDataLoader =
			// MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
			int index = 0;
			while (mdl.hasNext()) {
				// for (int index = 0; index < markers.length; index++) {
				if (index % 1000 == 0) {
					float usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
					float freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
					float maxMemory = Runtime.getRuntime().maxMemory();
					log.report(ext.getTime()	+ "\tData loaded = "
											+ Math.round(((double) index / (double) markers.length * 100.0))
											+ "%\tFree memory: "
											+ Math.round(((double) freeMemory / (double) maxMemory * 100.0)) + "%");
				}
				MarkerData markerData = mdl.next();
				float[] lrrs;
				if (recomputeLRR) {
					lrrs = markerData.getRecomputedLRR_BAF(null, null, false, 1, 0, null, true, true, log)[1];
				} else {
					lrrs = markerData.getLRRs();
				}

				if (params != null) {

					lrrs = markerData.getGCCorrectedLRRBAF(	params,
																									projectIndices.get(markerData.getMarkerName()),
																									proj.getLog())[1];
				}

				if (!hasNAN(lrrs)) {
					applyMarkerLoading(lrrs, index);
				} else if (imputeMeanForNaN) {
					lrrs = PrincipalComponentsCompute.imputeMeanForNaN(	markers[index], lrrs, samplesToUse,
																															log);
					applyMarkerLoading(lrrs, index);
				} else {
					log.reportError("Warning - marker "	+ markers[index]
													+ " contained a NaN value in the extrapolated dataset, skipping it for extrapolation");
				}
				index++;
				// markerDataLoader.releaseIndex(index);
			}
			// markerDataLoader.reportWaitTimes();
			applySingularValues();
		}
	}

	/**
	 * Summarize the newly built extrapolated principal components, in the standard Genvisis PC format
	 * of FID\tIID\tPC1\tPC2...
	 */
	public void reportExtropolatedPCs(String output) {
		SampleData sampleData = proj.getSampleData(0, false);
		extrapolatedPCsFile = output;
		try {
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + output)) {
				Files.backup(	output, proj.PROJECT_DIRECTORY.getValue(),
											proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}
			PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
																													+ output));
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
			log.reportError("Error: file \""	+ output
											+ "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error writing to file \"" + output + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
	}

	/**
	 * builds the extrapolated PCs for each sample/component combination by summing the loadings of
	 * marker (at each component) multiplied by the intensity value (for the marker at that sample)
	 *
	 * @param lrrs input data for all samples
	 * @param markerIndex used to match the current marker data with the marker loading
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
	 * After the loadings have been summed by applyMarkerLoading, the final step to build the
	 * extrapolated PCs is to divide each summation by the singular value for that component
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
	 * Properly formats extrapolatedPCs to an array organized as extrapolatedPCs[sample0][components
	 * for sample0]
	 */
	private void initExtPCS() {
		int numSamples = 0;
		for (boolean element : samplesToUse) {
			if (element) {
				numSamples++;
			}
		}
		extrapolatedPCs = new double[numSamples][numComponents];
	}

	private static boolean hasNAN(float[] lrrs) {
		// should not happen, but who knows
		if (lrrs == null) {
			return true;
		}
		for (float lrr : lrrs) {
			if (Float.isNaN(lrr)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Extract the marker names for all marker loadings. The data loaded from these markers will be
	 * used to build the extrapolated PCs. We assume that the markers in the marker loadings file were
	 * laid down in project order...If the markerDataLoader hangs, this is a good spot to sort
	 */
	private void getMarkers() {
		ArrayList<String> al = new ArrayList<String>();
		for (MarkerLoadings markerLoading : markerLoadings) {
			al.add(markerLoading.getMarker());
		}
		markers = al.toArray(new String[al.size()]);
	}

	/**
	 * Helper class to facilitate easy handling of the marker loadings
	 *
	 */
	public static class MarkerLoadings {
		private final String marker;
		private final double[] loadings;

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

		public static MarkerLoadings[] getLoadings(	String markerLoadingFile, int numComponents,
																								Logger log) {
			ArrayList<MarkerLoadings> ml = new ArrayList<MarkerLoadings>();
			try {
				BufferedReader reader = Files.getReader(markerLoadingFile, false, true, false);
				String[] line = reader.readLine().trim().split("\t");
				int[] indices = ext.indexFactors(LOADING_FIRST, line, true, true);
				if (indices == null || indices[0] != 0) {
					log.reportError("Error - Marker Loading file  must have "	+ Array.toStr(LOADING_FIRST)
													+ " in the first column");
					System.exit(1);
				}
				if ((line.length - 1) < numComponents) {
					log.reportError("Error - Can only apply "	+ (line.length - 1)
													+ " marker loadings (as provided in " + markerLoadingFile
													+ ", not enough for " + numComponents + " components");
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
							log.reportError("Error - could not parse marking loading value "	+ line[i]
															+ " to a double");
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
		private final String singularFile;
		private final double[] singularValues;
		private final int numComponents;
		private final Logger log;

		public SingularValues(String singularFile, int numComponents, Logger log) {
			this.singularFile = singularFile;
			this.numComponents = numComponents;
			this.log = log;
			singularValues = new double[numComponents];
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
					log.reportError("Error - singular value file must have header "
													+ Array.toStr(SINGULAR_HEADER));
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
				log.reportError("Error - not enough singular values were found in "	+ singularFile + " for "
												+ numComponents + " components");
				log.reportError("Please select a smaller number of components, or provide a file with more singular values");
			}
		}
	}
}
