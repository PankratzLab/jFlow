package one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import seq.analysis.Blast;
import seq.analysis.Blast.BlastResults;
import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
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
import common.Logger;
import common.PSF;
import common.ext;
import filesys.Segment;

public class MarkerBlastIterator {

	private static final String[] PLOT_FILE_HEADER = new String[] { "MarkerName", "PercentAppropriateMatch", "EvalAppropriateMatch", "AvgCrossHybPercentMatch", "AvgCrossHybEval", "AvgCrossHybLength", "NumCrossHyb" };
	private static final String[] X_COLUMNS_PLOT = new String[] { "AvgCrossHybPercentMatch", "AvgCrossHybEval", "AvgCrossHybLength", "NumCrossHyb" };
	private static final int[][] QC_GROUPINGS = new int[][] { { 2 }, { 3, 4, 5 }, { 6, 7 }, { 8, 9, 10 }, { 24 } };
	private static final String[] QC_TITLES = new String[] { "CallRate", "MeanClusterTheta", "DiffTheta", "SDClusterTheta", "LRR_SD" };

	public static void blastIter(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, int numThreads, boolean reportToTmp) {

		// int[] blastWordSizes = new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
		// int[] reportWordSizes = new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
		//
		//

		Logger log = proj.getLog();
		int[] blastWordSizes = new int[] { 10 };
		int[] reportWordSizes = new int[] { 0 };
		double appropriateMatchPercent = 100;
		double mafFilter = 0.05;
		double minCrossHybLength = 15;
		double minExpect = 0.05;

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
			if (Files.exists(catResults)) {
				new File(catResults).delete();
			}
			Files.cat(results[i].getTmpFiles(), catResults, new int[0], log);
			Hashtable<String, Integer> trackResults = new Hashtable<String, Integer>();
			ArrayList<MarkerIterationSummary> summaries = new ArrayList<MarkerIterationSummary>();

			try {
				BufferedReader reader = Files.getAppropriateReader(catResults);
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					if (line.length == Blast.BLAST_HEADER.length && !line[0].startsWith(Blast.BLAST_HEADER[0])) {

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
						} else if (blastResults.getEvalue() < minExpect && blastResults.getAlignmentLength() > minCrossHybLength && !blastSeg.overlaps(markerSeg)) {// we do not want to count probe A/B if probe A is perfect and probe B is not but is in the right spot
							summaries.get(summaryindex).getCrossHybEvalue().add(blastResults.getEvalue());
							summaries.get(summaryindex).getCrossHybPercentMatch().add(blastResults.getPercentIdentity());
							summaries.get(summaryindex).getCrossHybLength().add((double) blastResults.getAlignmentLength());
						}

					}
				}
				reader.close();
				String plotFileMatch = catResults + ".plot.match.txt";
				String plotFileNoMatch = catResults + ".plot.Nomatch.txt";

				PrintWriter writerMatch = Files.getAppropriateWriter(plotFileMatch);
				PrintWriter writerNoMatch = Files.getAppropriateWriter(plotFileNoMatch);

				writerMatch.println(Array.toStr(PLOT_FILE_HEADER) + "\t" + Array.toStr(Array.subArray(MarkerMetrics.FULL_QC_HEADER, 1)));
				writerNoMatch.println(Array.toStr(PLOT_FILE_HEADER) + "\t" + Array.toStr(Array.subArray(MarkerMetrics.FULL_QC_HEADER, 1)));
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

						PrintWriter currentWriter = summaries.get(j).isHasAppropriateMatch() ? writerMatch : writerNoMatch;
						numMatch = summaries.get(j).isHasAppropriateMatch() ? numMatch + 1 : numMatch;
						currentWriter.print(summaries.get(j).getSummary());

						for (int j2 = 0; j2 < parser.getNumericData().length; j2++) {
							currentWriter.print("\t" + parser.getNumericData()[j2][markerIndex]);
						}
						currentWriter.println();
					}

				}
				log.reportTimeInfo("Found a total of " + numMarkers + " markers that were blasted and " + numMatch + " matched the correct location at percent " + appropriateMatchPercent + " and maf > " + mafFilter);
				writerMatch.close();
				ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

				for (int j = 0; j < X_COLUMNS_PLOT.length; j++) {
					for (int j2 = 0; j2 < QC_GROUPINGS.length; j2++) {
						String[] yColumns = Array.subArray(MarkerMetrics.FULL_QC_HEADER, QC_GROUPINGS[j2]);
						rScatters.add(plot(proj, plotFileMatch, X_COLUMNS_PLOT[j], QC_TITLES[j2], yColumns));
						rScatters.add(plot(proj, plotFileNoMatch, X_COLUMNS_PLOT[j], QC_TITLES[j2], yColumns));

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

		public void setMarkerName(String markerName) {
			this.markerName = markerName;
		}

		public boolean isHasAppropriateMatch() {
			return hasAppropriateMatch;
		}

		public void setHasAppropriateMatch(boolean hasAppropriateMatch) {
			this.hasAppropriateMatch = hasAppropriateMatch;
		}

		public double getPercentAppropriateMatch() {
			return percentAppropriateMatch;
		}

		public void setPercentAppropriateMatch(double percentAppropriateMatch) {
			this.percentAppropriateMatch = percentAppropriateMatch;
		}

		public double getEvalueAppropriateMatch() {
			return evalueAppropriateMatch;
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
		boolean report = false;
		FILE_SEQUENCE_TYPE fSequence_TYPE = FILE_SEQUENCE_TYPE.AFFY_ANNOT;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to an Illumina manifest file  (i.e. fileSeq=" + fileSeq + " (default))\n" + "";
		usage += "   (3) number of threads to use  (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (4) word size for initial match in blast db  (i.e. blastWordSize=" + blastWordSize + " (default))\n" + "";
		usage += "   (5) number of base pairs with an exact match to report (i.e. reportWordSize=" + reportWordSize + " (default))\n" + "";
		usage += "   (6) sequence file type  (i.e. seqType=" + fSequence_TYPE + " (default))\n" + ", Options are:\n ";

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
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, false);

		blastIter(proj, fileSeq, fSequence_TYPE, numThreads, true);

	}
}
