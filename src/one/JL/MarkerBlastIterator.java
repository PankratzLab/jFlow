package one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import seq.analysis.Blast;
import seq.analysis.Blast.BlastResults;
import stats.Histogram.DynamicAveragingHistogram;
import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.ExtProjectDataParser.Builder;
import cnv.qc.MarkerBlast;
import cnv.qc.MarkerMetrics;
import cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import cnv.qc.MarkerBlast.MarkerBlastResult;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.Sort;
import common.ext;
import filesys.Segment;

public class MarkerBlastIterator {

	private static final String[] PLOT_FILE_HEADER = new String[] { "MarkerName", "PercentAppropriateMatch", "EvalAppropriateMatch", "AvgCrossHybPercentMatch", "AvgCrossHybEval", "AvgCrossHybLength", "NumCrossHyb" };
	private static final String[] X_COLUMNS_PLOT = new String[] { "AvgCrossHybPercentMatch", "AvgCrossHybEval", "AvgCrossHybLength", "NumCrossHyb" };
	private static final int[][] QC_GROUPINGS = new int[][] { { 1 }, { 2, 3, 4 }, { 5, 6 }, { 7, 8, 9 }, { 10, 11, 12 }, { 23 } };
	private static final String[] QC_TITLES = new String[] { "CallRate", "MeanClusterTheta", "DiffTheta", "SDClusterTheta", "MeanClusterR", "LRR_SD" };

	private enum BLAST_METRICS {
		CROSS_HYBE_PERCENT_MATCH(0, 1, 2), CROSS_HYBE_EVAL(0, 1, 2), CROSS_HYBE_LENGTH(10, 50, 0);

		private double minVal;
		private double maxVal;
		private int sigFigs;

		private BLAST_METRICS(double minVal, double maxVal, int sigFigs) {
			this.minVal = minVal;
			this.maxVal = maxVal;
			this.sigFigs = sigFigs;
		}

		public double getMinVal() {
			return minVal;
		}

		public void setMinVal(double minVal) {
			this.minVal = minVal;
		}

		public double getMaxVal() {
			return maxVal;
		}

		public void setMaxVal(double maxVal) {
			this.maxVal = maxVal;
		}

		public int getSigFigs() {
			return sigFigs;
		}

		public void setSigFigs(int sigFigs) {
			this.sigFigs = sigFigs;
		}

	}

	public static void blastIter(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, int blastWordSize, int numThreads, boolean reportToTmp, String[] otherMarkersToExtract) throws FileNotFoundException {
		Logger log = proj.getLog();
		final Hashtable<String, Integer> indices = proj.getMarkerIndices();
		final int[] blastWordSizes = new int[] { blastWordSize };
		final int[] reportWordSizes = new int[] { 0 };
		final byte[] chrs = proj.getMarkerSet().getChrs();
		final int[] pos = proj.getMarkerSet().getPositions();
		final double appropriateMatchPercent = 100;
		final double mafFilter = 0.0;
		final double minCrossHybLength = 0;
		final double maxExpect = 20;
		final int plotFontSize = 6;
		final double yMin = 0.00001;
		int oneHitWonderDef = proj.ARRAY_TYPE.getValue() == ARRAY.ILLUMINA ? 20 : 20;

		final int crossHybBufferDistance = proj.getArrayType().getProbeLength();
		if (!Files.exists(proj.MARKER_METRICS_FILENAME.getValue())) {
			proj.getLog().reportTimeInfo("Creating file " + proj.MARKER_METRICS_FILENAME.getValue());
			MarkerMetrics.fullQC(proj, null, null, numThreads);
		}

		MarkerBlastResult[] results = new MarkerBlastResult[blastWordSizes.length * reportWordSizes.length];
		int index = 0;
		for (int i = 0; i < blastWordSizes.length; i++) {
			for (int j = 0; j < reportWordSizes.length; j++) {
				results[index] = MarkerBlast.blastEm(proj, fileSeq, type, blastWordSizes[i], reportWordSizes[j], numThreads, reportToTmp);
				index++;
			}
		}

		for (int i = 0; i < results.length; i++) {
			Hashtable<String, Integer> trackResults = new Hashtable<String, Integer>();
			ArrayList<MarkerIterationSummary> summaries = new ArrayList<MarkerIterationSummary>();
			ArrayList<String> markersOffTargetPerfectMatch = new ArrayList<String>();
			for (int j = 0; j < results[i].getTmpFiles().length; j++) {
				proj.getLog().reportTimeInfo("Processing " + results[i].getTmpFiles()[j] + " blast results");

				try {
					BufferedReader reader = Files.getAppropriateReader(results[i].getTmpFiles()[j]);
					int numEntries = 0;
					while (reader.ready()) {
						String[] line = reader.readLine().trim().split("\t");
						if ((line.length == Blast.BLAST_HEADER.length-1||line.length == Blast.BLAST_HEADER.length) && !line[0].startsWith(Blast.BLAST_HEADER[0])) {
							numEntries++;
							if (numEntries % 100000 == 0) {
								proj.getLog().reportTimeInfo("Processed " + numEntries + " blast results");
								// reader.close();
								// break;
							}
							BlastResults blastResults = new BlastResults(line, log);
							String marker = blastResults.getQueryID();
							if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
								if (marker.endsWith("_A") || marker.endsWith("_B")) {
									marker = marker.substring(0, marker.length() - 2);
								} else {
									proj.getLog().reportTimeError("Query id did not end with _A or _B which is required for an AFFY array");
								}

							}
							boolean hasMarker = trackResults.containsKey(marker);
							int summaryindex = -1;
							if (hasMarker) {
								summaryindex = trackResults.get(marker);
							} else {
								summaryindex = summaries.size();
								trackResults.put(marker, summaryindex);
								summaries.add(new MarkerIterationSummary(marker));
							}

							Segment blastSeg = blastResults.getSegment();

							int markerIndex = indices.get(marker);
							Segment markerSeg = new Segment(chrs[markerIndex], pos[markerIndex] - 1, pos[markerIndex] + 1);
							if (markerSeg.getChr() > 0 && blastSeg.getChr() > 0) {
								boolean straight = blastResults.getMismatches() == 0 && blastResults.getGapOpens() == 0;
								boolean match = blastResults.getPercentIdentity() == appropriateMatchPercent && results[i].getSequenceSize() == blastResults.getAlignmentLength();
								if (match) {
									match = blastResults.getEvalue() < maxExpect;
								}
								if (blastSeg.overlaps(markerSeg) && straight && match) {
									summaries.get(summaryindex).setHasAppropriateMatch(true);
									summaries.get(summaryindex).setEvalueAppropriateMatch(blastResults.getEvalue());
									summaries.get(summaryindex).setPercentAppropriateMatch(blastResults.getPercentIdentity());
									// blastSeg.getChr()!=markerSeg.getChr()&&
								} else if (straight && blastResults.getEvalue() < maxExpect && blastResults.getAlignmentLength() > minCrossHybLength && !blastSeg.overlaps(markerSeg.getBufferedSegment(crossHybBufferDistance))) {// we do not want to count probe A/B if probe A is perfect and probe B is not but is in the right spot
									if (blastResults.getAlignmentLength() == proj.getArrayType().getProbeLength()) {
										if (!markersOffTargetPerfectMatch.contains(marker)) {
											markersOffTargetPerfectMatch.add(marker + "\t" + markerSeg.getUCSClocation() + blastSeg.getUCSClocation());
										}
									}

									summaries.get(summaryindex).getCrossHybEvalue().add(blastResults.getEvalue());
									summaries.get(summaryindex).getCrossHybPercentMatch().add(blastResults.getPercentIdentity());
									summaries.get(summaryindex).getCrossHybLength().add((double) blastResults.getAlignmentLength());
								}
							}
						}
					}

					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \"" + results[i].getTmpFiles()[j] + "\" not found in current directory");
				} catch (IOException ioe) {
					log.reportError("Error reading file \"" + results[i].getTmpFiles()[j] + "\"");
					return;
				}

			}
			Files.writeList(markersOffTargetPerfectMatch.toArray(new String[markersOffTargetPerfectMatch.size()]), results[i].getOutput() + ".OFF_target_PM.TXT");
			Builder builder = new ExtProjectDataParser.Builder();
			builder.separator("\t");
			builder.sampleBased(false);
			builder.requireAll(true);
			builder.dataKeyColumnName(MarkerMetrics.FULL_QC_HEADER[0]);
			ExtProjectDataParser parser = builder.build(proj, proj.MARKER_METRICS_FILENAME.getValue());
			parser.loadData();

			DynamicAveragingHistogram[][] dHistograms = new DynamicAveragingHistogram[BLAST_METRICS.values().length][parser.getNumericData().length];
			for (int k = 0; k < BLAST_METRICS.values().length; k++) {
				BLAST_METRICS metric = BLAST_METRICS.values()[k];
				switch (metric) {
				case CROSS_HYBE_EVAL:
					metric.setMaxVal(maxExpect);
					break;
				case CROSS_HYBE_LENGTH:
					metric.setMaxVal(proj.getArrayType().getProbeLength());
					break;
				case CROSS_HYBE_PERCENT_MATCH:
					break;
				default:
					proj.getLog().reportTimeError("Invalid Metric " + metric);
					break;
				}
				for (int l = 0; l < parser.getNumericData().length; l++) {
					dHistograms[k][l] = new DynamicAveragingHistogram(metric.getMinVal(), metric.getMaxVal(), metric.getSigFigs());
				}

			}
			int numMarkers = 0;
			ArrayList<String> oneHitters = new ArrayList<String>();
			ArrayList<String> notOneHitters = new ArrayList<String>();
			ArrayList<Double> notOneHittersMaxCrossHybe = new ArrayList<Double>();
			ArrayList<String> noAppropriateMatch = new ArrayList<String>();
			ArrayList<Double> noAppropriateMatchMaxCrossHybe = new ArrayList<Double>();
			MarkerSet markerSet = proj.getMarkerSet();
			for (int k = 0; k < summaries.size(); k++) {
				numMarkers++;
				int markerIndex = indices.get(summaries.get(k).getMarkerName());
				if (summaries.get(k).isHasAppropriateMatch()) {
					if (summaries.get(k).getCrossHybLength().size() == 0 || Array.max(Array.toDoubleArray(summaries.get(k).getCrossHybLength())) < oneHitWonderDef) {
						oneHitters.add(summaries.get(k).getMarkerName());
					} else {
						notOneHitters.add(summaries.get(k).getMarkerName() + "\t" + markerSet.getChrs()[markerIndex]);
						notOneHittersMaxCrossHybe.add(Array.max(Array.toDoubleArray(summaries.get(k).getCrossHybLength())));
					}
					if (parser.getNumericDataForTitle("MAF")[markerIndex] > mafFilter) {
						if (summaries.get(k).getCrossHybLength().size() > 0) {
							for (int l = 0; l < parser.getNumericData().length; l++) {

								for (int m = 0; m < BLAST_METRICS.values().length; m++) {
									BLAST_METRICS metric = BLAST_METRICS.values()[m];
									switch (metric) {
									case CROSS_HYBE_EVAL:
										dHistograms[m][l].addDataPair(Array.min(Array.toDoubleArray(summaries.get(k).getCrossHybEvalue())), parser.getNumericData()[l][markerIndex]);
										break;
									case CROSS_HYBE_LENGTH:
										dHistograms[m][l].addDataPair(Array.max(Array.toDoubleArray(summaries.get(k).getCrossHybLength())), parser.getNumericData()[l][markerIndex]);

										break;
									case CROSS_HYBE_PERCENT_MATCH:
										dHistograms[m][l].addDataPair((double) Array.max(Array.toDoubleArray(summaries.get(k).getCrossHybLength())) / proj.getArrayType().getProbeLength(), parser.getNumericData()[l][markerIndex]);
										break;
									default:
										proj.getLog().reportTimeError("Invalid Metric " + metric);

										break;
									}
								}
							}
						}
					}
				} else {
					noAppropriateMatch.add(summaries.get(k).getMarkerName() + "\t" + markerSet.getChrs()[markerIndex]);
					if (summaries.get(k).getCrossHybLength().size() > 0) {
						noAppropriateMatchMaxCrossHybe.add(Array.max(Array.toDoubleArray(summaries.get(k).getCrossHybLength())));
					} else {
						noAppropriateMatchMaxCrossHybe.add((double) 0);
					}
				}

			}
			String notOneHittersOut = ext.addToRoot(results[i].getOutput(), ".NOToneHitWonders_" + oneHitWonderDef);
			String notAppropriateMatchOut = ext.addToRoot(results[i].getOutput(), ".NOTAppropriateMatch_" + oneHitWonderDef);

			dumpMiss(log, oneHitWonderDef, notOneHittersOut, notOneHitters, notOneHittersMaxCrossHybe);
			dumpMiss(log, oneHitWonderDef, notAppropriateMatchOut, noAppropriateMatch, noAppropriateMatchMaxCrossHybe);

			String oneHitWonders = ext.addToRoot(results[i].getOutput(), ".oneHitWonders_" + oneHitWonderDef);
			String oneHitTargets = ext.addToRoot(proj.TARGET_MARKERS_FILENAME.getValue(), ".oneHitWonders_" + oneHitWonderDef);
			Files.writeArrayList(oneHitters, oneHitWonders);
			extractOneHitWondersFrom(proj.TARGET_MARKERS_FILENAME.getValue(), Array.toStringArray(oneHitters), oneHitTargets, log);
			if (otherMarkersToExtract != null) {
				for (int j = 0; j < otherMarkersToExtract.length; j++) {
					extractOneHitWondersFrom(otherMarkersToExtract[j], Array.toStringArray(oneHitters), ext.addToRoot(otherMarkersToExtract[j], ".oneHitWonders_" + oneHitWonderDef), log);
				}
			}
			String[] allPlotFiles = new String[BLAST_METRICS.values().length];
			String[] extras = new String[] { "_bin", "_count", "_average" };
			String resultsPlot = ext.addToRoot(results[i].getOutput(), ".hist.blasts");

			try {

				for (int j = 0; j < dHistograms.length; j++) {// for each blast metric
					String base = BLAST_METRICS.values()[j].toString();
					allPlotFiles[j] = resultsPlot + "." + base + ".txt";
					PrintWriter writer = new PrintWriter(new FileWriter(allPlotFiles[j]));
					String[] titles = new String[3 * parser.getNumericDataTitles().length];
					int titleIndex = 0;
					for (int j2 = 0; j2 < dHistograms[j].length; j2++) {// for each marker metric
						dHistograms[j][j2].average();
						String baseQC = base + "_V_" + parser.getNumericDataTitles()[j2];
						for (int k = 0; k < extras.length; k++) {
							titles[titleIndex] = baseQC + extras[k];
							writer.print((j2 == 0 && k == 0 ? "" : "\t") + titles[titleIndex]);
							titleIndex++;
						}

					}
					writer.println();

					for (int j2 = 0; j2 < dHistograms[j][0].getBins().length; j2++) {
						for (int k = 0; k < dHistograms[j].length; k++) {
							double bin = dHistograms[j][k].getBins()[j2];
							int counts = dHistograms[j][k].getCounts()[j2];
							double average = dHistograms[j][k].getAverages()[j2];
							writer.print((k == 0 ? "" : "\t") + bin + "\t" + counts + "\t" + average);
						}
						writer.println();
					}

					writer.close();
					ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

					for (int l = 0; l < QC_GROUPINGS.length; l++) {
						String groupPlot = ext.rootOf(allPlotFiles[j], false) + "_" + QC_TITLES[l];
						String title = "n=" + Array.sum(dHistograms[j][0].getCounts());
						ArrayList<String> ys = new ArrayList<String>();
						for (int k = 0; k < QC_GROUPINGS[l].length; k++) {
							ys.add(titles[3 * QC_GROUPINGS[l][k] + 2]);
						}
						String[] yColumns = ys.toArray(new String[ys.size()]);
						RScatter rScatterGroupAvg = new RScatter(allPlotFiles[j], groupPlot + ".rscript", ext.removeDirectoryInfo(groupPlot), groupPlot + ".pdf", titles[0], yColumns, SCATTER_TYPE.POINT, proj.getLog());
						rScatterGroupAvg.setxLabel(BLAST_METRICS.values()[j].toString());
						rScatterGroupAvg.setyLabel(QC_TITLES[l]);
						rScatterGroupAvg.setOverWriteExisting(false);
						rScatterGroupAvg.setFontsize(plotFontSize);
						rScatterGroupAvg.setTitle(title);
						rScatterGroupAvg.setyMin(yMin);
						rScatterGroupAvg.execute();
						rScatters.add(rScatterGroupAvg);

						groupPlot = ext.rootOf(allPlotFiles[j], false) + "_" + QC_TITLES[l] + "_counts";
						ys = new ArrayList<String>();
						for (int k = 0; k < QC_GROUPINGS[l].length; k++) {
							ys.add(titles[3 * QC_GROUPINGS[l][k] + 1]);
						}
						yColumns = ys.toArray(new String[ys.size()]);
						RScatter rScatterGroupCounts = new RScatter(allPlotFiles[j], groupPlot + ".rscript", ext.removeDirectoryInfo(groupPlot), groupPlot + ".pdf", titles[0], yColumns, SCATTER_TYPE.POINT, proj.getLog());
						rScatterGroupCounts.setxLabel(BLAST_METRICS.values()[j].toString());
						rScatterGroupCounts.setyMin(yMin);
						rScatterGroupCounts.setyLabel("COUNTS");
						rScatterGroupCounts.setLogTransformY(true);
						rScatterGroupCounts.setOverWriteExisting(false);
						rScatterGroupCounts.setFontsize(plotFontSize);
						rScatterGroupCounts.setTitle(title);
						rScatterGroupCounts.execute();
						rScatters.add(rScatterGroupCounts);

					}
					for (int k = 0; k < dHistograms[j].length; k++) {
						String plotRoot = ext.rootOf(allPlotFiles[j], false) + "_" + parser.getNumericDataTitles()[k];
						String xLabel = titles[3 * k];
						String yLabel = titles[3 * k + 2];
						// System.out.println("X\t"+xLabel);
						// System.out.println("Y\t"+yLabel);
						// try {
						// Thread.sleep(10000);
						// } catch (InterruptedException ie) {
						// }
						RScatter rScatter = new RScatter(allPlotFiles[j], plotRoot + ".rscript", ext.removeDirectoryInfo(plotRoot), plotRoot + ".pdf", xLabel, new String[] { yLabel }, SCATTER_TYPE.POINT, proj.getLog());
						rScatter.setxLabel(BLAST_METRICS.values()[j].toString());
						rScatter.setyLabel(parser.getNumericDataTitles()[k]);
						rScatter.setTitle(titles[3 * k]);
						rScatter.setOverWriteExisting(false);
						rScatter.execute();
						rScatters.add(rScatter);

					}
					RScatters rScattersAll = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), allPlotFiles[j] + ".rscript", allPlotFiles[j] + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, proj.getLog());

					rScattersAll.execute();

				}

			} catch (Exception e) {
				log.reportError("Error writing to " + resultsPlot);
				log.reportException(e);
			}

		}

	}

	private static void dumpMiss(Logger log, int oneHitWonderDef, String output, ArrayList<String> notOneHitters, ArrayList<Double> notOneHittersMaxCrossHybe) {
		int[] sorted = Sort.quicksort(Array.toDoubleArray(notOneHittersMaxCrossHybe));
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			for (int j = 0; j < sorted.length; j++) {
				writer.println(notOneHitters.get(sorted[j]) + "\t" + notOneHittersMaxCrossHybe.get(sorted[j]));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	public static void blastIter(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, int blastWordSize, int numThreads, boolean reportToTmp) {

		// int[] blastWordSizes = new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
		// int[] reportWordSizes = new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
		//
		//

		Logger log = proj.getLog();
		final int[] blastWordSizes = new int[] { blastWordSize };
		final int[] reportWordSizes = new int[] { 0 };
		final double appropriateMatchPercent = 100;
		final double mafFilter = 0.05;
		final double minCrossHybLength = 15;
		final double minExpect = 0.05;

		if (!Files.exists(proj.MARKER_METRICS_FILENAME.getValue())) {
			proj.getLog().reportTimeInfo("Creating file " + proj.MARKER_METRICS_FILENAME.getValue());
			MarkerMetrics.fullQC(proj, null, null, numThreads);
		}

		MarkerBlastResult[] results = new MarkerBlastResult[blastWordSizes.length * reportWordSizes.length];
		int index = 0;
		for (int i = 0; i < blastWordSizes.length; i++) {
			for (int j = 0; j < reportWordSizes.length; j++) {
				results[index] = MarkerBlast.blastEm(proj, fileSeq, type, blastWordSizes[i], reportWordSizes[j], numThreads, reportToTmp);
				index++;
			}
		}
		byte[] chrs = proj.getMarkerSet().getChrs();
		int[] pos = proj.getMarkerSet().getPositions();

		Hashtable<String, Integer> indices = proj.getMarkerIndices();

		for (int i = 0; i < results.length; i++) {
			String catResults = ext.addToRoot(results[i].getOutput(), ".all.blasts");
			if (!Files.exists(catResults)) {
				Files.cat(results[i].getTmpFiles(), catResults, new int[0], log);
			}
			Hashtable<String, Integer> trackResults = new Hashtable<String, Integer>();
			ArrayList<MarkerIterationSummary> summaries = new ArrayList<MarkerIterationSummary>();

			try {
				BufferedReader reader = Files.getAppropriateReader(catResults);
				int numEntries = 0;
				while (reader.ready()) {

					String[] line = reader.readLine().trim().split("\t");
					if (line.length == Blast.BLAST_HEADER.length && !line[0].startsWith(Blast.BLAST_HEADER[0])) {
						numEntries++;
						if (numEntries % 100000 == 0) {
							proj.getLog().reportTimeInfo("Processed " + numEntries + " blast results");
						}
						BlastResults blastResults = new BlastResults(line, log);
						String marker = blastResults.getQueryID();
						if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
							if (marker.endsWith("_A") || marker.endsWith("_B")) {
								marker = marker.substring(0, marker.length() - 2);
							} else {
								proj.getLog().reportTimeError("Query id did not end with _A or _B which is required for an AFFY array");
							}

						}
						boolean hasMarker = trackResults.containsKey(marker);
						int summaryindex = -1;
						if (hasMarker) {
							summaryindex = trackResults.get(marker);
						} else {
							summaryindex = summaries.size();
							trackResults.put(marker, summaryindex);
							summaries.add(new MarkerIterationSummary(marker));
						}

						Segment blastSeg = blastResults.getSegment();
						int markerIndex = indices.get(marker);
						Segment markerSeg = new Segment(chrs[markerIndex], pos[markerIndex] - 1, pos[markerIndex] + 1);
						boolean straight = blastResults.getMismatches() == 0 && blastResults.getGapOpens() == 0;
						boolean match = blastResults.getPercentIdentity() == appropriateMatchPercent && results[i].getSequenceSize() == blastResults.getAlignmentLength();
						if (match) {
							match = blastResults.getEvalue() < minExpect;
						}
						if (blastSeg.overlaps(markerSeg) && straight && match) {
							summaries.get(summaryindex).setHasAppropriateMatch(true);
							summaries.get(summaryindex).setEvalueAppropriateMatch(blastResults.getEvalue());
							summaries.get(summaryindex).setPercentAppropriateMatch(blastResults.getPercentIdentity());
						} else if (straight && blastResults.getEvalue() < minExpect && blastResults.getAlignmentLength() > minCrossHybLength && !blastSeg.overlaps(markerSeg)) {// we do not want to count probe A/B if probe A is perfect and probe B is not but is in the right spot
							summaries.get(summaryindex).getCrossHybEvalue().add(blastResults.getEvalue());
							summaries.get(summaryindex).getCrossHybPercentMatch().add(blastResults.getPercentIdentity());
							summaries.get(summaryindex).getCrossHybLength().add((double) blastResults.getAlignmentLength());
						}
					}
				}
				reader.close();
				String plotFileMatch = catResults + ".plot.match.txt";
				String plotFileALL = catResults + ".plot.ALL.txt";

				PrintWriter writerMatch = Files.getAppropriateWriter(plotFileMatch);
				PrintWriter writerALL = Files.getAppropriateWriter(plotFileALL);

				writerMatch.println(Array.toStr(PLOT_FILE_HEADER) + "\t" + Array.toStr(Array.subArray(MarkerMetrics.FULL_QC_HEADER, 1)));
				writerALL.println(Array.toStr(PLOT_FILE_HEADER) + "\t" + Array.toStr(Array.subArray(MarkerMetrics.FULL_QC_HEADER, 1)));
				Builder builder = new ExtProjectDataParser.Builder();
				builder.separator("\t");
				builder.sampleBased(false);
				builder.requireAll(true);
				builder.dataKeyColumnName(MarkerMetrics.FULL_QC_HEADER[0]);
				ExtProjectDataParser parser = builder.build(proj, proj.MARKER_METRICS_FILENAME.getValue());
				parser.loadData();
				int numMarkers = 0;
				int numMatch = 0;
				for (int j = 0; j < summaries.size(); j++) {
					numMarkers++;
					int markerIndex = indices.get(summaries.get(j).getMarkerName());
					if (parser.getNumericDataForTitle("MAF")[markerIndex] > mafFilter) {
						// writerALL.println(summaries.get(j).getMarkerName()+"\t"+summaries.get(j).getPercentAppropriateMatch()+"\t"+summaries.get(j).getEvalueAppropriateMatch()+"\t"+ summaries.get(j).getCrossHybEvalue().size());
						//
						// for (int j2 = 0; j2 < parser.getNumericData().length; j2++) {
						// writerALL.print("\t" + parser.getNumericData()[j2][markerIndex]);
						// }
						for (int k = 0; k < summaries.get(j).getCrossHybEvalue().size(); k++) {
							double percentMatch = summaries.get(j).getCrossHybPercentMatch().get(k);
							double eval = summaries.get(j).getCrossHybEvalue().get(k);
							double length = summaries.get(j).getCrossHybLength().get(k);

							writerALL.print(summaries.get(j).getMarkerName() + "\t" + percentMatch + "\t" + eval + "\t" + percentMatch + "\t" + eval + "\t" + length + "\t" + summaries.get(j).getCrossHybLength().size());

							for (int j2 = 0; j2 < parser.getNumericData().length; j2++) {
								writerALL.print("\t" + parser.getNumericData()[j2][markerIndex]);
							}
							writerALL.println();
						}
						// }

						// String summary = markerName;
						// summary += "\t" + percentAppropriateMatch;
						// summary += "\t" + evalueAppropriateMatch;
						// summary += "\t" + (crossHybPercentMatch.size() > 0 ? Array.mean(Array.toDoubleArray(crossHybPercentMatch), true) : "0");
						// summary += "\t" + (crossHybEvalue.size() > 0 ? Array.mean(Array.toDoubleArray(crossHybEvalue), true) : "1");
						// summary += "\t" + (crossHybLength.size() > 0 ? Array.mean(Array.toDoubleArray(crossHybLength), true) : "0");
						// summary += "\t" + crossHybEvalue.size();
						//
						// PrintWriter currentWriter = summaries.get(j).isHasAppropriateMatch() ? writerMatch : writerALL;

						// numMatch = summaries.get(j).isHasAppropriateMatch() ? numMatch + 1 : numMatch;
						// currentWriter.print(summaries.get(j).getSummary());
						//
						// for (int j2 = 0; j2 < parser.getNumericData().length; j2++) {
						// currentWriter.print("\t" + parser.getNumericData()[j2][markerIndex]);
						// }
						// currentWriter.println();
					}

				}
				log.reportTimeInfo("Found a total of " + numMarkers + " markers that were blasted and " + numMatch + " matched the correct location at percent " + appropriateMatchPercent + " and maf > " + mafFilter);
				writerMatch.close();
				ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

				for (int j = 0; j < X_COLUMNS_PLOT.length; j++) {
					for (int j2 = 0; j2 < QC_GROUPINGS.length; j2++) {
						String[] yColumns = Array.subArray(MarkerMetrics.FULL_QC_HEADER, QC_GROUPINGS[j2]);
						// rScatters.add(plot(proj, plotFileMatch, X_COLUMNS_PLOT[j], QC_TITLES[j2], yColumns));
						rScatters.add(plot(proj, plotFileALL, X_COLUMNS_PLOT[j], QC_TITLES[j2], yColumns));
					}
				}
				RScatters rScattersAll = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), plotFileMatch + ".rscript", plotFileMatch + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, proj.getLog());

				rScattersAll.execute();

			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + catResults + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + catResults + "\"");
				return;
			}

		}

		// has match and has cross hyb with % id gt than X;

	}

	private static void extractOneHitWondersFrom(String toExtractFile, String[] oneHitters, String outputFile, Logger log) {
		String[] toExtracts = HashVec.loadFileToStringArray(toExtractFile, false, new int[] { 0 }, false);
		int[] indices = ext.indexLargeFactors(toExtracts, oneHitters, true, log, false, false);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			PrintWriter writerNot = new PrintWriter(new FileWriter(ext.addToRoot(outputFile, "Not")));

			for (int i = 0; i < indices.length; i++) {
				if (indices[i] >= 0) {
					writer.println(toExtracts[i]);
				} else {
					writerNot.println(toExtracts[i]);
				}
			}
			writerNot.close();
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}
	}

	private static RScatter plot(Project proj, String plotFile, String xColumn, String title, String[] yColumns) {
		String rootOut = ext.addToRoot(plotFile, xColumn + "_" + title);
		RScatter rScatterMatch = new RScatter(plotFile, rootOut + ".rscript", ext.rootOf(rootOut), rootOut + ".jpeg", xColumn, yColumns, SCATTER_TYPE.POINT, proj.getLog());
		rScatterMatch.setTitle(ext.removeDirectoryInfo(plotFile) + "_" + xColumn);
		rScatterMatch.setxLabel(xColumn);
		rScatterMatch.setyLabel(title);
		rScatterMatch.setOverWriteExisting(true);
		rScatterMatch.execute();
		return rScatterMatch;
	}

	private static class MarkerIterationSummary {
		private String markerName;
		private boolean hasAppropriateMatch;
		private double percentAppropriateMatch;
		private double evalueAppropriateMatch;
		private ArrayList<Double> crossHybPercentMatch, crossHybEvalue, crossHybLength;

		private MarkerIterationSummary(String markerName) {
			super();
			this.markerName = markerName;
			this.hasAppropriateMatch = false;
			this.percentAppropriateMatch = Double.NaN;
			this.evalueAppropriateMatch = Double.NaN;
			this.crossHybPercentMatch = new ArrayList<Double>();
			this.crossHybEvalue = new ArrayList<Double>();
			this.crossHybLength = new ArrayList<Double>();
		}

		public String getSummary() {
			String summary = markerName;
			summary += "\t" + percentAppropriateMatch;
			summary += "\t" + evalueAppropriateMatch;
			summary += "\t" + (crossHybPercentMatch.size() > 0 ? Array.mean(Array.toDoubleArray(crossHybPercentMatch), true) : "0");
			summary += "\t" + (crossHybEvalue.size() > 0 ? Array.mean(Array.toDoubleArray(crossHybEvalue), true) : "1");
			summary += "\t" + (crossHybLength.size() > 0 ? Array.mean(Array.toDoubleArray(crossHybLength), true) : "0");
			summary += "\t" + crossHybEvalue.size();

			return summary;
		}

		public ArrayList<Double> getCrossHybPercentMatch() {
			return crossHybPercentMatch;
		}

		public ArrayList<Double> getCrossHybLength() {
			return crossHybLength;
		}

		public ArrayList<Double> getCrossHybEvalue() {
			return crossHybEvalue;
		}

		public String getMarkerName() {
			return markerName;
		}

		public boolean isHasAppropriateMatch() {
			return hasAppropriateMatch;
		}

		public void setHasAppropriateMatch(boolean hasAppropriateMatch) {
			this.hasAppropriateMatch = hasAppropriateMatch;
		}

		public void setPercentAppropriateMatch(double percentAppropriateMatch) {
			this.percentAppropriateMatch = percentAppropriateMatch;
		}

		public void setEvalueAppropriateMatch(double evalueAppropriateMatch) {
			this.evalueAppropriateMatch = evalueAppropriateMatch;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		// String fastaDb = "hg19_canonical.fa";
		String fileSeq = "HumanExome-12-v1-0-B.csv";
		int numThreads = 4;
		int blastWordSize = 30;
		int reportWordSize = 40;
		String[] otherMarkersToExtract = null;
		FILE_SEQUENCE_TYPE fSequence_TYPE = FILE_SEQUENCE_TYPE.AFFY_ANNOT;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to an Illumina manifest file  (i.e. fileSeq=" + fileSeq + " (default))\n" + "";
		usage += "   (3) number of threads to use  (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (4) word size for initial match in blast db  (i.e. blastWordSize=" + blastWordSize + " (default))\n" + "";
		usage += "   (5) number of base pairs with an exact match to report (i.e. reportWordSize=" + reportWordSize + " (default))\n" + "";
		usage += "   (6) sequence file type  (i.e. seqType=" + fSequence_TYPE + " (default))\n" + ", Options are:\n ";
		usage += "   (7) list of marker filenames to subset to one-hit wonders (i.e. marks= (no default))\n";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("fileSeq=")) {
				fileSeq = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("blastWordSize=")) {
				blastWordSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("reportWordSize=")) {
				reportWordSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("seqType=")) {
				fSequence_TYPE = FILE_SEQUENCE_TYPE.valueOf(ext.parseStringArg(args[i], FILE_SEQUENCE_TYPE.MANIFEST_FILE.toString()));
				numArgs--;
			} else if (args[i].startsWith("marks=")) {
				otherMarkersToExtract = ext.parseStringArg(args[i], "").split(",");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, false);

		try {
			blastIter(proj, fileSeq, fSequence_TYPE, blastWordSize, numThreads, true, otherMarkersToExtract);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
