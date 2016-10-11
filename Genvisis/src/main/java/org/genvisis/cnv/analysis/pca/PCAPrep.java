package org.genvisis.cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

/**
 * Methods for preparing a Project for PCA
 *
 */
public class PCAPrep {
	static final String MARKERS_TO_QC_FILE = "markers_to_QC.txt";
	static final String MARKERS_FOR_ABCALLRATE = "markers_ABCallRate.txt";


	public static int prepPCA(Project proj, int numThreads, String outputBase, boolean markerQC,
														double markerCallRateFilter, String useFile, SampleList sampleList,
														Logger log) {
		int[] counts;
		String markersForABCallRate;
		String markersForEverythingElse;
		// if marker QC is not flagged, sample qc is based on all target markers by default
		String markersToQC = proj.PROJECT_DIRECTORY.getValue() + outputBase + "_" + MARKERS_TO_QC_FILE;
		String markersABCallrate = proj.PROJECT_DIRECTORY.getValue()	+ outputBase + "_"
																+ MARKERS_FOR_ABCALLRATE;
		String baseLineMarkers = proj.PROJECT_DIRECTORY.getValue()	+ outputBase
															+ "_baselineMarkers.txt";
		String[] auto = proj.getAutosomalMarkers();
		ArrayList<String> tmp = new ArrayList<String>();
		for (int i = 0; i < auto.length; i++) {
			if (!proj.getArrayType().isCNOnly(auto[i])) {
				tmp.add(auto[i]);
			}
		}
		Files.writeIterable(tmp, baseLineMarkers);
		if (markerQC) {
			String markerQCFile = outputBase + "_markerQC.txt";
			proj.MARKER_METRICS_FILENAME.setValue(markerQCFile);
			qcMarkers(proj, baseLineMarkers, markersToQC, markersABCallrate, markerCallRateFilter,
								numThreads);
			markersForABCallRate = markersABCallrate;
			if (!Files.exists(markersForABCallRate)) {
				log.reportError("Error - markerQC was flagged but the file "	+ markersABCallrate
												+ " could not be found");
				return 1;
			}
		} else {
			markersForABCallRate = baseLineMarkers;
			writeMarkersToQC(proj, baseLineMarkers, markersToQC);
		}
		markersForEverythingElse = markersToQC;
		String qcFile = outputBase + "_lrr_sd.txt";
		proj.SAMPLE_QC_FILENAME.setValue(qcFile);

		counts = org.genvisis.cnv.qc.LrrSd.filterSamples(	proj, outputBase, markersForABCallRate,
																											markersForEverythingElse, numThreads, useFile,
																											false);
		if (counts == null || counts[1] != sampleList.getSamples().length) {
			if (counts == null || counts[1] == 0 && Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
				log.reportError("Error - was unable to parse QC file "	+ proj.SAMPLE_QC_FILENAME.getValue()
												+ ", backing up this file to "
												+ proj.BACKUP_DIRECTORY.getValue(false, false)
												+ " and re-starting sample qc");
				Files.backup(	ext.removeDirectoryInfo(proj.SAMPLE_QC_FILENAME.getValue()),
											proj.PROJECT_DIRECTORY.getValue(),
											proj.BACKUP_DIRECTORY.getValue(true, false), true);
			}
			counts = org.genvisis.cnv.qc.LrrSd.filterSamples(	proj, outputBase, markersForABCallRate,
																												markersForEverythingElse, numThreads,
																												useFile, false);
			if (counts == null || counts[1] != sampleList.getSamples().length) {
				if (counts == null) {
					log.reportError("Error - could not parse QC file ("	+ proj.SAMPLE_QC_FILENAME.getValue()
													+ ")");
				} else {
					log.reportError("Error - different number of samples (n="	+ counts[1]
													+ ") listed in the QC file (" + proj.SAMPLE_QC_FILENAME.getValue()
													+ ") compared to the number of samples in the project (n="
													+ sampleList.getSamples().length + ")");
					log.reportError("      - delete the QC file ("	+ proj.SAMPLE_QC_FILENAME.getValue()
													+ ") to regenerate it with the correct number of samples");
				}
				log.reportError("aborting...");
				return 2;
			}
		}
		if (counts == null || counts[0] == 0) {// no samples passed threshold, null case shouldn't
																						// happen but we will test anyway
			return 2;// message handled already
		}
		return 42;
	}

	/**
	 * Currently un-neccesary, but it is set up in case we want to QC the median markers at the same
	 * time
	 */
	private static void writeMarkersToQC(	Project proj, String targetMarkersFile,
																				String markersToQCFile) {
		String[] markers = null;
		if (targetMarkersFile == null) {
			markers = proj.getMarkerNames();
		} else {
			String[] markersToQC = {targetMarkersFile};
			markers = setMarkersToQC(proj, markersToQC);
		}
		// TODO remove CNVi probe markers (here?)
		Files.writeArray(markers, markersToQCFile);
	}

	/**
	 * Similar to above, currently un-neccesary, but it is set up in case we want to QC the median
	 * markers at the same time
	 */
	private static String[] setMarkersToQC(Project proj, String[] files) {
		ArrayList<String> markersToUse = new ArrayList<String>();
		for (String file : files) {
			String[] markers = HashVec.loadFileToStringArray(file, false, new int[] {0}, false);
			for (String marker : markers) {
				markersToUse.add(marker);
			}
		}
		return PrincipalComponentsCompute.sortByProjectMarkers(	proj,
																														markersToUse.toArray(new String[markersToUse.size()]));
	}

	/**
	 * Write a filtered list of markers to use for ab callRate in sample QC
	 */
	private static boolean filterMarkerMetricsFile(	Project proj, double markerCallRateFilter,
																									String markersABCallrate) {
		ArrayList<String> abMarkersToUse = new ArrayList<String>();
		BufferedReader reader;
		Logger log = proj.getLog();

		try {
			// reader = Files.getReader(proj.getFilename(proj.MARKER_METRICS_FILENAME), false, true,
			// false);
			reader = Files.getReader(proj.MARKER_METRICS_FILENAME.getValue(), false, true, false);
			String[] header = reader.readLine().trim().split("\t");
			int abIndex = ext.indexOfStr(MarkerMetrics.FULL_QC_HEADER[2], header);
			if (abIndex == -1) {
				// log.reportError("Error - the necessary marker metrics header " +
				// MarkerMetrics.FULL_QC_HEADER[2] + " was not found in the marker metrics file" +
				// proj.getFilename(proj.MARKER_METRICS_FILENAME));
				log.reportError("Error - the necessary marker metrics header "
													+ MarkerMetrics.FULL_QC_HEADER[2]
												+ " was not found in the marker metrics file"
												+ proj.MARKER_METRICS_FILENAME.getValue());
				return false;
			} else {
				String[] metrics;
				while (reader.ready()) {
					// proj.getArrayType().isCNOnly(markerName)
					metrics = reader.readLine().trim().split("\t");
					try {
						double callRate = Double.parseDouble(metrics[abIndex]);
						if (callRate >= markerCallRateFilter) {
							abMarkersToUse.add(metrics[0]);
						}
					} catch (NumberFormatException nfe) {
						log.report("Warning - found an invalid number "	+ metrics[abIndex] + " for marker"
												+ metrics[0] + " skipping this marker");
					}
				}
			}
			if (abMarkersToUse.size() == 0) {
				log.reportError("Error - no markers passed the callRate threshold. Please consider lowering threshold, or ensure that markers can have call rates (not cnv only probes)");
				return false;
			} else {
				log.report("Sample call rate will be computed with " + abMarkersToUse.size() + " markers");
				Files.writeArray(	abMarkersToUse.toArray(new String[abMarkersToUse.size()]),
													markersABCallrate);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			// log.reportError("Error: file \"" + proj.getFilename(proj.MARKER_METRICS_FILENAME) + "\" not
			// found in current directory");
			log.reportError("Error: file \""	+ proj.MARKER_METRICS_FILENAME.getValue()
											+ "\" not found in current directory");
		} catch (IOException ioe) {
			// log.reportError("Error reading file \"" + proj.getFilename(proj.MARKER_METRICS_FILENAME) +
			// "\"");
			log.reportError("Error reading file \"" + proj.MARKER_METRICS_FILENAME.getValue() + "\"");
		}

		return true;
	}

	/**
	 * @param proj
	 * @param targetMarkersFile
	 * @param markersToQCFile
	 * @param markersABCallrate
	 * @param markerCallRateFilter
	 * @param numthreads
	 */
	private static void qcMarkers(Project proj, String targetMarkersFile, String markersToQCFile,
																String markersABCallrate, double markerCallRateFilter,
																int numthreads) {
		Logger log;
		String markerMetricsFilename;

		log = proj.getLog();
		markerMetricsFilename = proj.MARKER_METRICS_FILENAME.getValue(true, false);
		// skip if marker qc file exists
		if (Files.exists(markerMetricsFilename)	&& new File(markerMetricsFilename).length() > 0
				&& Files.exists(markersToQCFile)
				&& Files.countLines(markerMetricsFilename, 1) >= Files.countLines(markersToQCFile, 0)) {
			log.report("Marker QC file "	+ proj.MARKER_METRICS_FILENAME.getValue(true, false)
									+ " exists");
			log.report("Skipping Marker QC computation for the analysis, filtering on existing file");
		} else {
			log.report("Computing marker QC for "
									+ (targetMarkersFile == null	? "all markers in project."
																								: "markers in " + targetMarkersFile));
			writeMarkersToQC(proj, targetMarkersFile, markersToQCFile);
			boolean[] samplesToExclude = new boolean[proj.getSamples().length];
			Arrays.fill(samplesToExclude, false);

			PSF.checkInterrupted();
			MarkerMetrics.fullQC(	proj, samplesToExclude, ext.removeDirectoryInfo(markersToQCFile), false,
														numthreads);

			// MarkerMetrics.fullQC(proj, samplesToExclude, null, false, numthreads);
		}

		PSF.checkInterrupted();
		filterMarkerMetricsFile(proj, markerCallRateFilter, markersABCallrate);
	}
}
