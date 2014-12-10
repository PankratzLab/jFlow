package seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

public class Picard {
	public static final String SORT_SAM_JAR = "SortSam.jar";
	public static final String MARK_DUPLICATES_JAR = "MarkDuplicates.jar";
	public static final String BUILD_BAM_INDEX = "BuildBamIndex.jar";
	public static final String BAM_INDEX_EXT = ".bai";
	public static final String JAR = "-jar";
	public static final String DEFAULT_JAVA = "java";
	public static final String INPUT = "INPUT=";
	public static final String OUTPUT = "OUTPUT=";
	public static final String SORT_ORDER = "SORT_ORDER=";
	public static final String DEFAULT_SORT_ORDER = "coordinate";
	public static final String PICARD_LOCATION_COMMAND = "picard=";
	public static final String METRICS_FILE = "METRICS_FILE=";
	public static final String TMP_DIR = "TMP_DIR=";// necessary since it puts tmp files in strange places that run out of space

	private String picardLocation;
	private String javaLocation;
	boolean fail, verbose, overwriteExisting;
	private Logger log;

	public Picard(String picardLocation, String javaLocation, boolean overwriteExisting, boolean verbose, Logger log) {
		super();
		this.picardLocation = picardLocation;
		this.javaLocation = (javaLocation == null ? DEFAULT_JAVA : javaLocation);
		this.fail = verifyPicardLocationg(picardLocation);
		this.overwriteExisting = overwriteExisting;
		this.verbose = verbose;
		this.log = log;
	}

	public String getPicardLocation() {
		return picardLocation;
	}

	public boolean sortSam(String inputFile, String outputFile, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, picardLocation + SORT_SAM_JAR, INPUT + inputFile, OUTPUT + outputFile, SORT_ORDER + DEFAULT_SORT_ORDER, TMP_DIR + ext.parseDirectoryOfFile(outputFile) };
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { inputFile }, new String[] { outputFile }, verbose, overwriteExisting, true, (altLog == null ? log : altLog));
	}

	public boolean markDuplicates(String inputFile, String outputFile, String metricsFile, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, picardLocation + MARK_DUPLICATES_JAR, INPUT + inputFile, OUTPUT + outputFile, METRICS_FILE + metricsFile, TMP_DIR + ext.parseDirectoryOfFile(outputFile) };
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { inputFile }, new String[] { outputFile, metricsFile }, verbose, overwriteExisting, true, (altLog == null ? log : altLog));
	}

	public boolean indexBAM(String inputFile, String expectedOutput, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, picardLocation + BUILD_BAM_INDEX, INPUT + inputFile, TMP_DIR + ext.parseDirectoryOfFile(expectedOutput) };
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { inputFile }, new String[] { expectedOutput }, verbose, overwriteExisting, true, (altLog == null ? log : altLog));
	}

	private static boolean verifyPicardLocationg(String picardLocation) {
		return new File(picardLocation).exists();
	}

	public Picard_Analysis picardASam(String baseId, String fullPathToSamFile, Logger altLog) {
		boolean progress = false;
		Picard_Analysis picard_Analysis = new Picard_Analysis(baseId, fullPathToSamFile, (altLog == null ? log : altLog));
		picard_Analysis.parseInput();
		progress = sortSam(picard_Analysis.getFullPathToSamFile(), picard_Analysis.getFullPathToSortedBamFile(), picard_Analysis.getLog());
		if (progress) {
			progress = markDuplicates(picard_Analysis.getFullPathToSortedBamFile(), picard_Analysis.getFullPathToSortedDeDuppedBamFile(), picard_Analysis.getFullPathToMetricsTxt(), picard_Analysis.getLog());
			if (progress) {
				progress = indexBAM(picard_Analysis.getFullPathToSortedDeDuppedBamFile(), picard_Analysis.getFullPathToSortedDeDuppedBamFileIndex(), picard_Analysis.getLog());
				if (progress) {
					picard_Analysis.setAllThere(progress);
				}
			}
		}
		picard_Analysis.setFail(!progress);
		return picard_Analysis;
	}

	public static class Picard_Analysis {
		public static final String BAM = ".bam";
		public static final String SORTED = ".sorted";
		public static final String DE_DUPPED = ".dedup";
		public static final String INDEXED = ".bai";
		public static final String METRICS = ".metrics.txt";

		private String fullPathToSamFile, fullPathToSortedBamFile, fullPathToSortedDeDuppedBamFile, fullPathToSortedDeDuppedBamFileIndex, fullPathToMetricsTxt, baseID;
		private boolean allThere, fail;
		private Logger log;

		public Picard_Analysis(String baseId, String fullPathToSamFile, Logger log) {
			super();
			this.baseID = baseId;
			this.fullPathToSamFile = fullPathToSamFile;
			this.fail = false;
			this.log = log;
		}

		public Logger getLog() {
			return log;
		}

		public String getBaseID() {
			return baseID;
		}

		public void parseInput() {
			this.fullPathToSortedBamFile = ext.rootOf(fullPathToSamFile, false) + SORTED + BAM;
			this.fullPathToSortedDeDuppedBamFile = ext.addToRoot(fullPathToSortedBamFile, DE_DUPPED);
			this.fullPathToMetricsTxt = ext.rootOf(fullPathToSortedDeDuppedBamFile, false) + METRICS;
			this.fullPathToSortedDeDuppedBamFileIndex = ext.rootOf(fullPathToSortedDeDuppedBamFile, false) + INDEXED;
		}

		public String getFullPathToMetricsTxt() {
			return fullPathToMetricsTxt;
		}

		public String getFullPathToSamFile() {
			return fullPathToSamFile;
		}

		public String getFullPathToSortedBamFile() {
			return fullPathToSortedBamFile;
		}

		public String getFullPathToSortedDeDuppedBamFile() {
			return fullPathToSortedDeDuppedBamFile;
		}

		public String getFullPathToSortedDeDuppedBamFileIndex() {
			return fullPathToSortedDeDuppedBamFileIndex;
		}

		public boolean isAllThere() {
			return allThere;
		}

		public void setAllThere(boolean allThere) {
			this.allThere = allThere;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

	}

	public static class PicardMetricsParser {
		public static final String[] PICARD_METRICS = { "LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE" };
		private String[] picardMetricsFiles;
		private double[][] consolidatedMetrics;
		private boolean fail;

		private Logger log;

		public PicardMetricsParser(String[] picardMetricsFiles, Logger log) {
			super();
			this.picardMetricsFiles = picardMetricsFiles;
			this.consolidatedMetrics = new double[picardMetricsFiles.length][PICARD_METRICS.length];
			this.fail = false;
			this.log = log;
		}

		public void parse(String outputSummaryFile) {
			parseMetricsFiles();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(outputSummaryFile));
				writer.println("SAMPLE\t" + Array.toStr(PICARD_METRICS));
				for (int i = 0; i < consolidatedMetrics.length; i++) {
					writer.println(ext.rootOf(picardMetricsFiles[i]) + "\t" + Array.toStr(consolidatedMetrics[i]));
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + outputSummaryFile);
				log.reportException(e);
			}
		}

		private void parseMetricsFiles() {
			for (int i = 0; i < picardMetricsFiles.length; i++) {
				consolidatedMetrics[i] = getMetricsFromFile(picardMetricsFiles[i], log);
			}
		}

		public boolean isFail() {
			return fail;
		}

		private static double[] getMetricsFromFile(String picardMetricsFile, Logger log) {
			double[] curMetrics = new double[PICARD_METRICS.length];
			Arrays.fill(curMetrics, Double.NaN);
			try {
				BufferedReader reader = Files.getAppropriateReader(picardMetricsFile);
				String[] line;
				do {
					line = reader.readLine().trim().split("[\\s]+");
				} while (reader.ready() && (ext.indexFactors(new String[][] { PICARD_METRICS }, line, false, true, false, false)[0] == -1));
				if (!reader.ready()) {
					log.reportError("Error - could not find neccesary header " + Array.toStr(PICARD_METRICS) + " in file " + picardMetricsFile);
				} else {
					int[] indices = ext.indexFactors(line, PICARD_METRICS, true, false);
					line = reader.readLine().trim().split("[\\s]+");
					for (int i = 0; i < line.length; i++) {
						try {
							curMetrics[i] = Double.parseDouble(line[indices[i]]);
						} catch (NumberFormatException numberFormatException) {
							log.reportError("Error - invalid number on line " + Array.toStr(line) + " in file " + picardMetricsFile);
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + picardMetricsFile + "\" not found in current directory");

			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + picardMetricsFile + "\"");

			}
			return curMetrics;

		}
	}

}
