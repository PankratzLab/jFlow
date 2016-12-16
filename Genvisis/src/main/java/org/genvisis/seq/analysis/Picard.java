package org.genvisis.seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import com.google.common.collect.Lists;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Picard {
	public static final String PICARD_JAR = "picard.jar";
	public static final String SORT_SAM = "SortSam";
	public static final String MARK_DUPLICATES = "MarkDuplicates";
	public static final String BUILD_BAM_INDEX = "BuildBamIndex";
	public static final String BAM_INDEX_EXT = ".bai";
	public static final String JAR = "-jar";
	public static final String DEFAULT_JAVA = "java";
	public static final String INPUT = "INPUT=";
	public static final String OUTPUT = "OUTPUT=";
	public static final String SORT_ORDER = "SORT_ORDER=";
	public static final String SORTING_COLLECTION_SIZE_RATIO = "SORTING_COLLECTION_SIZE_RATIO=";
	public static final double DEFAULT_SORTING_COLLECTION_SIZE_RATIO = 0.25;

	public static final String DEFAULT_SORT_ORDER = "coordinate";
	public static final String PICARD_LOCATION_COMMAND = "picard=";
	public static final String METRICS_FILE = "METRICS_FILE=";
	public static final String TMP_DIR = "TMP_DIR=";// necessary since it puts tmp files in strange
																									// places that run out of space

	public static final String BAM = ".bam";
	public static final String DE_DUPPED = ".dedup";
	public static final String INDEXED = ".bai";
	public static final String METRICS = ".metrics.txt";

	private final String picardLocation;
	private final String javaLocation;
	boolean fail, verbose, overwriteExisting;
	private final Logger log;

	public Picard(String picardLocation, String javaLocation, boolean overwriteExisting,
								boolean verbose, Logger log) {
		super();
		this.picardLocation = picardLocation;
		this.javaLocation = (javaLocation == null ? DEFAULT_JAVA : javaLocation);
		fail = verifyPicardLocation(picardLocation);
		this.overwriteExisting = overwriteExisting;
		this.verbose = verbose;
		this.log = log;
	}

	public String getPicardLocation() {
		return picardLocation;
	}

	public boolean sortSam(String inputFile, String outputFile, Logger altLog) {
		String[] command = new String[] {	javaLocation, JAR, picardLocation + PICARD_JAR, SORT_SAM,
																			INPUT + inputFile, OUTPUT + outputFile,
																			SORT_ORDER + DEFAULT_SORT_ORDER,
																			getTMPdirectory(inputFile, outputFile)};
		return CmdLine.runCommandWithFileChecks(command, "", new String[] {inputFile},
																						new String[] {outputFile}, verbose, overwriteExisting,
																						true, (altLog == null ? log : altLog));
	}

	public boolean markDuplicates(String inputFile, String outputFile, String metricsFile,
																double memoryRatio, Logger altLog) {
		return markDuplicates(new String[] {inputFile}, outputFile, metricsFile, memoryRatio, altLog);
	}

	public boolean markDuplicates(String[] inputFiles, String outputFile, String metricsFile,
																double memoryRatio, Logger altLog) {
		ArrayList<String> commands = Lists.newArrayList(javaLocation, JAR, picardLocation + PICARD_JAR,
																										MARK_DUPLICATES, OUTPUT + outputFile,
																										METRICS_FILE + metricsFile,
																										getTMPdirectory(inputFiles[0], inputFiles[0]),
																										SORTING_COLLECTION_SIZE_RATIO + memoryRatio);
		for (int i = 0; i < inputFiles.length; i++) {
			commands.add(INPUT + inputFiles[i]);
		}
		return CmdLine.runCommandWithFileChecks(commands.toArray(new String[commands.size()]), "",
																						inputFiles, new String[] {outputFile, metricsFile},
																						verbose, overwriteExisting, true,
																						altLog == null ? log : altLog);
	}

	public boolean indexBAM(String inputFile, String expectedOutput, Logger altLog) {
		String[] command = new String[] {	javaLocation, JAR, picardLocation + PICARD_JAR,
																			BUILD_BAM_INDEX, INPUT + inputFile,
																			getTMPdirectory(inputFile, expectedOutput)};
		return CmdLine.runCommandWithFileChecks(command, "", new String[] {inputFile},
																						new String[] {expectedOutput}, verbose,
																						overwriteExisting, true,
																						(altLog == null ? log : altLog));
	}

	private static boolean verifyPicardLocation(String picardLocation) {
		return new File(picardLocation).exists();
	}

	private static String getTMPdirectory(String inputFile, String expectedOutput) {

		return TMP_DIR	+ ext.parseDirectoryOfFile(expectedOutput) + ext.rootOf(inputFile)
						+ "picard_tmp/";
	}

	public PicardAnalysis picardASam(	String baseId, String fullPathToSamFile, double memoryRatio,
																		Logger altLog) {
		boolean progress = false;
		PicardAnalysis picard_Analysis = new PicardAnalysis(baseId, fullPathToSamFile,
																												(altLog == null ? log : altLog));
		picard_Analysis.parseInput();
		progress = sortSam(	picard_Analysis.getFullPathToSamFile(),
												picard_Analysis.getFullPathToSortedBamFile(), picard_Analysis.getLog());
		if (progress) {
			progress = markDuplicates(picard_Analysis.getFullPathToSortedBamFile(),
																picard_Analysis.getFullPathToSortedDeDuppedBamFile(),
																picard_Analysis.getFullPathToMetricsTxt(), memoryRatio,
																picard_Analysis.getLog());
			if (progress) {
				progress = indexBAM(picard_Analysis.getFullPathToSortedDeDuppedBamFile(),
														picard_Analysis.getFullPathToSortedDeDuppedBamFileIndex(),
														picard_Analysis.getLog());
				if (progress) {
					picard_Analysis.setAllThere(progress);
				}
			}
		}
		picard_Analysis.setFail(!progress);
		return picard_Analysis;
	}

	public PicardMergeDedupe mergeDedupeBams(	String baseID, String[] fullPathsToBamFiles,
																						String[] fullPathsToBamFileIndices, String outputDir,
																						double memoryRatio, Logger altLog) {
		boolean progress = true;
		PicardMergeDedupe picardMergeDedupe = new PicardMergeDedupe(fullPathsToBamFiles,
																																fullPathsToBamFileIndices,
																																outputDir, baseID,
																																altLog == null ? log : altLog);
		picardMergeDedupe.parseInput();
		if (picardMergeDedupe.shouldMerge()) {
			progress = markDuplicates(fullPathsToBamFiles,
																picardMergeDedupe.getFullPathToMergedDedupedBam(),
																picardMergeDedupe.getFullPathToMetricsTxt(), memoryRatio, altLog);
			if (progress) {
				progress = indexBAM(picardMergeDedupe.getFullPathToMergedDedupedBam(),
														picardMergeDedupe.getFullPathToMergedDedupedBamFileIndex(),
														picardMergeDedupe.getLog());
			}
		}
		picardMergeDedupe.setFail(!progress);
		return picardMergeDedupe;
	}

	public static class PicardMergeDedupe {
		public static final String MERGED = ".merged";

		private final String[] fullPathsToInputBams;
		private final String[] fullPathsToInputBamIndices;
		private final String outputDir;
		private String fullPathToMergedDedupedBam;
		private String fullPathToMergedDedupedBamFileIndex;
		private String fullPathToMetricsTxt;
		private final String baseID;
		private final boolean shouldMerge;
		private boolean fail;
		private final Logger log;

		/**
		 * @param fullPathsToInputBams Paths to sorted BAM files to merge and dedupe
		 * @param fullPathsToInputBamIndices Paths to indices associated with BAMs in
		 *        fullPathsToInputBams
		 * @param outputDir Directory to output merged, deduped BAM to
		 * @param baseID
		 * @param log
		 */
		public PicardMergeDedupe(	String[] fullPathsToInputBams, String[] fullPathsToInputBamIndices,
															String outputDir, String baseID, Logger log) {
			super();
			if (fullPathsToInputBams == null	|| fullPathsToInputBamIndices == null
					|| fullPathsToInputBamIndices.length != fullPathsToInputBams.length) {
				throw new IllegalArgumentException(this.getClass().getName()
																						+ " requires fullPathsToInputBams and fullPathsToInputBamIndices of equal length");
			}
			this.fullPathsToInputBams = fullPathsToInputBams;
			this.fullPathsToInputBamIndices = fullPathsToInputBamIndices;
			this.shouldMerge = fullPathsToInputBams != null && fullPathsToInputBams.length > 1;
			this.outputDir = outputDir;
			this.baseID = baseID;
			this.log = log;
		}

		public void parseInput() {
			if (shouldMerge) {
				fullPathToMergedDedupedBam = outputDir + baseID + MERGED + BAM;
				fullPathToMergedDedupedBamFileIndex =
																						ext.rootOf(fullPathToMergedDedupedBam, false) + INDEXED;
				fullPathToMetricsTxt = ext.rootOf(fullPathToMergedDedupedBam, false) + METRICS;
			} else {
				fullPathToMergedDedupedBam = fullPathsToInputBams[0];
				fullPathToMergedDedupedBamFileIndex = fullPathsToInputBamIndices[0];
				fullPathToMetricsTxt = null;
			}
		}

		public boolean shouldMerge() {
			return shouldMerge;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String[] getFullPathsToInputBams() {
			return fullPathsToInputBams;
		}

		public String getOutputDir() {
			return outputDir;
		}

		public String getFullPathToMergedDedupedBam() {
			return fullPathToMergedDedupedBam;
		}

		public String getFullPathToMergedDedupedBamFileIndex() {
			return fullPathToMergedDedupedBamFileIndex;
		}

		public String getFullPathToMetricsTxt() {
			return fullPathToMetricsTxt;
		}

		public String getBaseID() {
			return baseID;
		}

		public Logger getLog() {
			return log;
		}


	}


	public static class PicardAnalysis {
		public static final String SORTED = ".sorted";

		private final String fullPathToSamFile;
		private String fullPathToSortedBamFile;
		private String fullPathToSortedDeDuppedBamFile;
		private String fullPathToSortedDeDuppedBamFileIndex;
		private String fullPathToMetricsTxt;
		private final String baseID;
		private String newBaseID;
		private boolean allThere, fail;
		private final Logger log;
		private String barcode;

		public PicardAnalysis(String baseId, String fullPathToSamFile, Logger log) {
			super();
			baseID = baseId;
			this.fullPathToSamFile = fullPathToSamFile;
			fail = false;
			this.log = log;
		}

		public Logger getLog() {
			return log;
		}

		public String getBaseID() {
			return baseID;
		}

		public String getNewBaseID() {
			return newBaseID;
		}

		public void setNewBaseID(String newBaseId) {
			this.newBaseID = newBaseId;
		}

		public void parseInput() {
			fullPathToSortedBamFile = ext.rootOf(fullPathToSamFile, false) + SORTED + BAM;
			fullPathToSortedDeDuppedBamFile = ext.addToRoot(fullPathToSortedBamFile, DE_DUPPED);
			fullPathToMetricsTxt = ext.rootOf(fullPathToSortedDeDuppedBamFile, false) + METRICS;
			fullPathToSortedDeDuppedBamFileIndex = ext.rootOf(fullPathToSortedDeDuppedBamFile, false)
																							+ INDEXED;
			if (baseID.split(BWA_Analysis.FileNameParser.SPLIT).length != 3) {
				barcode = "";
				log.reportTimeWarning("The current baseId "	+ baseID + " did not have 3 "
															+ BWA_Analysis.FileNameParser.SPLIT
															+ " - delimited fields, assuming no barcodes are present in the ids");
			} else {
				barcode = baseID.split(BWA_Analysis.FileNameParser.SPLIT)[2];
			}
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

		public String getBarcode() {
			return barcode;
		}

		public void setBarcode(String barcode) {
			this.barcode = barcode;
		}

	}

	public static class PicardMetricsParser {
		public static final String[] PICARD_METRICS =
																								{	"LIBRARY", "UNPAIRED_READS_EXAMINED",
																									"READ_PAIRS_EXAMINED", "UNMAPPED_READS",
																									"UNPAIRED_READ_DUPLICATES",
																									"READ_PAIR_DUPLICATES",
																									"READ_PAIR_OPTICAL_DUPLICATES",
																									"PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE"};
		private final String[] picardMetricsFiles;
		private final double[][] consolidatedMetrics;
		private final boolean fail;

		private final Logger log;

		public PicardMetricsParser(String[] picardMetricsFiles, Logger log) {
			super();
			this.picardMetricsFiles = picardMetricsFiles;
			consolidatedMetrics = new double[picardMetricsFiles.length][PICARD_METRICS.length];
			fail = false;
			this.log = log;
		}

		public void parse(String outputSummaryFile) {
			parseMetricsFiles();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(outputSummaryFile));
				writer.println("SAMPLE\t" + Array.toStr(PICARD_METRICS));
				for (int i = 0; i < consolidatedMetrics.length; i++) {
					writer.println(ext.rootOf(picardMetricsFiles[i])	+ "\t"
													+ Array.toStr(consolidatedMetrics[i]));
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
				} while (reader.ready() && (ext.indexFactors(	new String[][] {PICARD_METRICS}, line, false,
																											true, false, false)[0] == -1));
				if (!reader.ready()) {
					log.reportError("Error - could not find neccesary header "	+ Array.toStr(PICARD_METRICS)
													+ " in file " + picardMetricsFile);
				} else {
					int[] indices = ext.indexFactors(PICARD_METRICS, line, true, false);
					line = reader.readLine().trim().split("[\\s]+");
					for (int i = 0; i < indices.length; i++) {
						try {
							curMetrics[i] = Double.parseDouble(line[indices[i]]);
						} catch (NumberFormatException numberFormatException) {
							log.reportError("Error - invalid number on line "	+ Array.toStr(line) + " in file "
															+ picardMetricsFile);
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
