package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ChrPositionMap;
import org.genvisis.common.Files;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SegmentLists;
import org.genvisis.filesys.SnpMarkerSet;

/**
 * Program to summarize p-values from a Genome-Wide Association Study (GWAS) and bins marker
 * depending on their p-values and depending upon whether they are inside or outside of a DNAse
 * hypersensitive (DNAseHS) region. These segments are loaded from BED files and the name of the
 * celltype is taken from the first part of the filename. The p-value column is not required to be
 * in a consistent column number. As a result the program creates a file with the following
 * tab-delimited columns: celltype, bin1(0.1<pvalue<1.0), bin2(0.01<pvalue<0.1),
 * bin3(0.001<pvalue<0.01)... Then this output is plotted using LinePlot
 *
 * Author: Rohit Sinha Date: 10/24/13 Time: 2:48 AM
 */
public class DnaseEnrichment {

	private final static String[][] P_VALUE_FILE_HEADERS = {{"MarkerName", "Marker", "SNP"}, {"Chr"},
																													{"Position", "Pos"},
																													{"P-value", "p-val"}};
	private final static String INSIDE_REGION = "insideRegion";
	private final static String TOTAL_MARKERS = "totalMarkers";
	private final static int FAILURE = 1;
	private final static String BED_FILE_EXTENTION = ".bed";
	private final static String LD_FILES_EXTENTION = ".ld";
	private final static String OUTPUT_DELIMITER = "\t";
	private final static String DATA_SEPARATOR = ":";
	private final static Logger LOGGER = Logger.getLogger(DnaseEnrichment.class.getName()); // logger
																																													// for
																																													// this
																																													// class
	private static String OUTPUT_FILENAME = "DnaseEnrichment.xln";
	private final static String BED_FILE_CHR_MAP_PART_COUNT_FILENAME = "BedFileChrMapCount.xln";
	private final static String BED_FILE_CHR_MAP_FOLDER = "BedChrPositionMap";
	private final static String OUTPUT_FOLDER = "DnaseEnrichmentOutput";
	private final static Pattern WHILTE_SPACE_PATTERN = Pattern.compile("[\\s]+");
	private static int maxBinSize = 0;
	private static String CHR_MAP_FILE_ID = "_chr_map";
	private static String REFERENCE_MAP_FILENAME = "plink.bim";
	private static String FILTERED_FILE_EXTENSION = ".flt";
	private static String plinkFile = null;
	private static String ldDir = null;
	private static String bedDir = null;
	private static String[] bedFileList;
	private static int ldLines = 0;
	private static int numThreads = 1; // number of threads to run. default is 1
	private static boolean performLD = true;

	private static ArrayList<HashSet<String>> dhsregionsHashSetList; // global static variable to hold
																																	 // DHS regions
	private static Hashtable<String, Integer> bedFileChrMapPartCount = new Hashtable<String, Integer>(); // global
																																																			 // static
																																																			 // variable
																																																			 // for
																																																			 // ChrPositionMap

	/**
	 * WorkerThreads which process different LD files and create ChrPositionMap concurrently
	 */
	public class WorkerThread implements Runnable {
		private final String ldFilePath; // path of the ld file
		private final int chrNum; // chr number

		public WorkerThread(String s, int chrNum) {
			ldFilePath = s;
			this.chrNum = chrNum;
		}

		@Override
		public void run() {
			LOGGER.info(Thread.currentThread().getName() + ": Building ChrPositionMap for chr file: "
									+ ldFilePath);
			generateChrPositionMap(ldFilePath, chrNum); // generate the chrPositionMap for all the bed
																									// files using this chr file
			LOGGER.info(Thread.currentThread().getName()
									+ ": Completed building ChrPositionMap for file: " + ldFilePath);

		}
	}

	/**
	 * Function to print the command line usage.
	 *
	 * @return a {@link String}: the command line usage help
	 */
	private static String showCommandLineHelp() {

		return ("\n" + "gwas.DnaseEnrichment requires six arguments\n"
						+ "\t(1) directory with bed files (bedDir)\n"
						+ "\t(2) filename with chr/pos/pvals (pValueFile)\n"
						+ "\t(3) The bim file (pLinkFile)\n" + "\t(4) " + "ld files directory (ldDir)\n"
						+ "\t(5)number of lines on which ld files will be split (ldLines)\n"
						+ "\t(6)number of threads (numThreads: defaults=1)\n"
						+ "\t(7)Perform LD calculation ? (ld= default: true)\n" + "");

	}

	public static void main(String[] args) {

		ArrayList<OutputFileFormat> overlapStats = null;
		int numArgs = args.length;
		String filename = null;
		boolean overWrite = false;

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(showCommandLineHelp());
				System.exit(1);
			} else if (arg.toLowerCase().startsWith("beddir=")) {
				bedDir = ext.verifyDirFormat(arg.split("=")[1]);
				numArgs--;
			} else if (arg.toLowerCase().startsWith("pvaluefile=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.toLowerCase().startsWith("plinkfile=")) {
				plinkFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.toLowerCase().startsWith("lddir=")) {
				ldDir = ext.verifyDirFormat(arg.split("=")[1]);
				numArgs--;
			} else if (arg.toLowerCase().startsWith("ldlines=")) {
				ldLines = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.toLowerCase().startsWith("overwrite=")) {
				overWrite = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.toLowerCase().startsWith("numthreads=")) {
				numThreads = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.toLowerCase().startsWith("ld=")) {
				performLD = ext.parseBooleanArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0 || bedDir == null || plinkFile == null || ldDir == null || ldLines == 0
				|| filename == null) {
			System.err.println(showCommandLineHelp());
			System.exit(1);
		}

		// get names of all the bed files
		bedFileList = Files.list(bedDir, BED_FILE_EXTENTION, false);
		Arrays.sort(bedFileList);

		// Print out parameters to user
		LOGGER.info("Starting DnaseEnrichment Processing");
		LOGGER.info("LD Directory is: " + ldDir);
		LOGGER.info("Bed Directory is: " + bedDir);
		LOGGER.info("Plink File location is: " + plinkFile);
		LOGGER.info("LD Line Count is : " + ldLines);
		LOGGER.info("Num of threads to use while making ChrPositionMap: " + numThreads);
		LOGGER.info("Over write ChrPositionMap is set to: " + overWrite);
		LOGGER.info("Over write ChrPositionMap is set to: " + performLD);

		if (performLD) {
			// Filter ld files by storing only the new required columns
			LOGGER.info("Filtering LD files. Please wait...");
			filterAllLDFiles();
			if (Files.exists(bedDir + BED_FILE_CHR_MAP_FOLDER + File.separator
											 + BED_FILE_CHR_MAP_PART_COUNT_FILENAME)) {
				LOGGER.info("We found ChrPositionMap for bed files.");
				LOGGER.info("Note: If you feel like your ld files or bed files have changed from the last run then "
										+ ""
										+ "we strongly suggest you to stop here and delete the ChrMapPosition directory and run the program again");
				LOGGER.info("Chr Position Map File Directory: " + bedDir + BED_FILE_CHR_MAP_FOLDER
										+ File.separator);
				LOGGER.info("overWrite parameter is set to: " + overWrite);

				if (overWrite || (!Files.exists(bedDir + BED_FILE_CHR_MAP_FOLDER + File.separator
																				+ BED_FILE_CHR_MAP_PART_COUNT_FILENAME))) {
					LOGGER.info("Finding DHS Region markers...");
					dhsregionsHashSetList = findDHSRegionMarkers();
					LOGGER.info("Starting to build ChrPositionMap...");
					DnaseEnrichment dnaseEnrichmentObject = new DnaseEnrichment();
					dnaseEnrichmentObject.runWorkers();
					writeBedFileChrMapPartCount(bedFileChrMapPartCount);
					dhsregionsHashSetList.clear();
				} else {
					bedFileChrMapPartCount = readBedFileChrMapPartCount();
				}

				// find the overlapping position counts
				LOGGER.info("Building statistics. Please wait...");
				overlapStats = findOverlapRegions(bedDir, filename, bedFileChrMapPartCount);
			}

		} else {
			System.out.println("Skipping filtering LD files as perform LD is specified as: " + performLD);
			// find the overlapping position counts
			LOGGER.info("Building statistics. Please wait...");
			overlapStats = findOverlapRegions(bedDir, filename, null);
		}

		// Overlap statistics
		LOGGER.info(overlapStats.toString());

		// write to output file
		writeOutputFile(overlapStats);

	}

	/**
	 * Function to run worker threads concurrently to build chrPositionMap
	 */
	public void runWorkers() {
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		String[] ldFilesList = Files.list(ldDir, LD_FILES_EXTENTION, false);
		for (int i = 0; i < ldFilesList.length; i++) {
			Runnable worker = new WorkerThread(ldFilesList[i], i);
			executor.execute(worker);
		}
		executor.shutdown();
		try {
			executor.awaitTermination(7, TimeUnit.DAYS);
		} catch (InterruptedException e) {
		}
		LOGGER.info("All threads terminated successfully: Processed all chr files to build ChrPositionMap");
	}

	/**
	 * Function to write overlap output to output file.
	 *
	 * @param overlapStats : an {@link ArrayList} of {@link OutputFileFormat} which is the overlap
	 *        statistics
	 */
	private static void writeOutputFile(List<OutputFileFormat> overlapStats) {
		FileWriter fstream;
		BufferedWriter out;

		List<List<Object>> resultArrayList = new ArrayList<List<Object>>();

		for (OutputFileFormat curRecord : overlapStats) {
			ArrayList<Object> resultArray = new ArrayList<Object>();
			resultArray.add(curRecord.file);

			for (double[] element : curRecord.ratio) {
				// to avoid divide by 0 error
				double ratio = element[0] / (element[1] != 0 ? element[1] : 1);
				// write all the ratio
				resultArray.add(String.valueOf(ratio) + DATA_SEPARATOR + String.valueOf(element[0]) + "/"
												+ String.valueOf(element[1]));
			}
			resultArrayList.add(resultArray);
		}
		resultArrayList = ArrayUtils.transpose(resultArrayList);
		try {
			String outputFilePath = bedDir + OUTPUT_FOLDER + File.separator + OUTPUT_FILENAME;

			File theDir = new File(bedDir + OUTPUT_FOLDER + File.separator);
			// if the directory does not exist, create it
			if (!theDir.exists()) {
				System.out.println("creating directory: " + bedDir + OUTPUT_FOLDER + File.separator);
				boolean result = theDir.mkdir();

				if (result) {
					System.out.println("DIR created");
				}
			}
			fstream = new FileWriter(outputFilePath, false);
			out = new BufferedWriter(fstream);
			for (List<Object> curRecord : resultArrayList) {
				for (Object curString : curRecord) {
					out.write(curString.toString());
					out.write(OUTPUT_DELIMITER);
				}
				out.newLine();
			}
			out.close();
			LOGGER.info("Output File: " + outputFilePath);
		} catch (IOException e) {
			LOGGER.info("Unable to write to temp output file");
		}
	}

	/**
	 * Function to find overlapping positions in segments
	 *
	 * @param dir : the directory path containing the segments file
	 * @param filename : the filename of the p-value file
	 * @return a Map<String, Map<Integer, Map<String, Integer>>> which shows the inside and outside
	 *         the region stats for all the segments file
	 */
	private static ArrayList<OutputFileFormat> findOverlapRegions(String dir, String filename,
																																Hashtable<String, Integer> bedFileChrMapPartCount) {

		ArrayList<OutputFileFormat> thisFileOutput = new ArrayList<OutputFileFormat>();

		// read the pvalues records from the pvalue file in memory
		ArrayList<PValueFileFormat> pValueRecords = readPValueFile(filename);
		OUTPUT_FILENAME = ext.rootOf(filename) + "_" + OUTPUT_FILENAME;

		LOGGER.info("Starting to process segment files...");

		for (String element : bedFileList) {
			LOGGER.info("Processing: " + element);
			Segment[][] segs = getSegments(dir + element);
			TreeMap<Integer, Map<String, Long>> overlapStats = countOverlaps(segs, pValueRecords, element,
																																			 bedFileChrMapPartCount);
			double[][] ratioList = new double[maxBinSize + 1][2]; // array for holding numerator and
																														// denominator
			double cumBinMarkers = 0;
			double cumBinMarkersInside = 0;
			int key;
			NavigableSet<Integer> keySet = overlapStats.descendingKeySet();
			// get the numerator for all the bins cumulatively
			for (Integer integer : keySet) {
				key = integer;
				Map<String, Long> value = overlapStats.get(key);

				double totalMarkersInside = value.containsKey(INSIDE_REGION) ? value.get(INSIDE_REGION) : 0;
				double totalMarkers = value.containsKey(TOTAL_MARKERS) ? value.get(TOTAL_MARKERS) : 0;

				cumBinMarkers += totalMarkers;
				cumBinMarkersInside += totalMarkersInside;

				ratioList[key][0] = cumBinMarkersInside / cumBinMarkers;
			}
			// add the same denominator which we calculated just now for all the bins
			for (Map.Entry<Integer, Map<String, Long>> entry : overlapStats.entrySet()) {
				ratioList[entry.getKey()][1] = cumBinMarkersInside / cumBinMarkers;
			}
			thisFileOutput.add(new OutputFileFormat(element, ratioList));
		}
		return thisFileOutput;
	}

	/**
	 * Function to find the makers in DHS region for all the bedfiles
	 *
	 * @return {@link ArrayList} of {@link HashSet} of String which are the marker names in DHS region
	 */
	private static ArrayList<HashSet<String>> findDHSRegionMarkers() {

		ArrayList<HashSet<String>> dhsregionsHashSetList = new ArrayList<HashSet<String>>();
		PlinkFile plinkFileContents = null;

		String dir = ext.parseDirectoryOfFile(plinkFile);
		if (new File(dir + REFERENCE_MAP_FILENAME).exists()) {
			SnpMarkerSet markerSet = new SnpMarkerSet(plinkFile, false,
																								new org.genvisis.common.Logger(null));
			plinkFileContents = new PlinkFile(markerSet.getMarkerNames(), markerSet.getChrs(),
																				markerSet.getPositions());

			ArrayList<Segment> plinkSegments = createSegmentsFromPlinkFile(plinkFileContents);
			LOGGER.info("Starting to find DHS regions for bed files");
			for (String element : bedFileList) {
				LOGGER.info("Processing for DHS region: " + element);
				Segment[][] segs = getSegments(bedDir + element);

				HashSet<String> ldMarkerHashSet = new HashSet<String>();

				for (int j = 0; j < plinkSegments.size(); j++) {
					// if there is a overlap
					if (Segment.binarySearchForOverlap(plinkSegments.get(j),
																						 segs[plinkSegments.get(j).getChr()]) != -1) {
						ldMarkerHashSet.add(plinkFileContents.markerNames[j]);
					}
				}
				dhsregionsHashSetList.add(ldMarkerHashSet);
			}
		} else {
			LOGGER.warning("Error - could not find plink.bim; required to define segments and perform a crucial datacheck");
		}
		return dhsregionsHashSetList;
	}

	/**
	 * Function to create segments from {@link PlinkFile}
	 *
	 * @param plinkFileContents {@link PlinkFile} containing plink file contents
	 * @return {@link ArrayList} of {@link Segment}
	 */
	private static ArrayList<Segment> createSegmentsFromPlinkFile(PlinkFile plinkFileContents) {
		// get chr, position and marker names from the plink file
		LOGGER.info("Reading plink file. This might take couple of seconds. Please wait ...");
		ArrayList<Segment> plinkSegments = new ArrayList<Segment>();

		// if something was wrong in the plink file data fields inform the user
		if (plinkFileContents.markerNames.length != plinkFileContents.chrs.length
				&& plinkFileContents.chrs.length == plinkFileContents.positions.length) {
			LOGGER.log(Level.WARNING, "size of markers, chr, positions are not equal in the plink file");
		}

		// if reading plink file was successful then create segments out of it
		if (plinkFileContents != null) {
			for (int i = 0; i < plinkFileContents.markerNames.length; i++) {
				// create a segment from this line
				Segment curSeg = new Segment(plinkFileContents.chrs[i], plinkFileContents.positions[i],
																		 plinkFileContents.positions[i]);
				plinkSegments.add(curSeg);
			}
		}
		return plinkSegments;
	}

	/**
	 * Function to {@link ChrPositionMap} for all the bed files
	 */
	private static void generateChrPositionMap(String curLdFilePath, int chrNum) {
		ChrPositionMap chrPositionMap;
		HashSet<Integer> value;
		int partCount = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(ldDir
																																+ Files.removeExtention(curLdFilePath)
																																+ FILTERED_FILE_EXTENSION));
			reader.readLine(); // skill the first header line

			while (reader.ready()) {
				Hashtable<String, ChrPositionMap> chrPositionMapList = new Hashtable<String, ChrPositionMap>();
				LOGGER.info("Processing: " + curLdFilePath + Files.SERIALIZED_FILE_EXTENSION);

				for (int ldLineCounter = 0; ldLineCounter < ldLines && reader.ready(); ldLineCounter++) {
					String curline = reader.readLine();
					String[] curlineParams = WHILTE_SPACE_PATTERN.split(curline.trim());
					for (int j = 0; j < dhsregionsHashSetList.size(); j++) {

						if (chrPositionMapList.containsKey(bedFileList[j])) {
							chrPositionMap = chrPositionMapList.get(bedFileList[j]);
						} else {
							chrPositionMapList.put(bedFileList[j], chrPositionMap = new ChrPositionMap());

						}
						// if any of the maker is present then store both the positions
						if (dhsregionsHashSetList.get(j).contains(curlineParams[2])
								|| dhsregionsHashSetList.get(j).contains(curlineParams[4])) {
							if (chrPositionMap.getChrPositionMap().containsKey(Byte.valueOf(curlineParams[0]))) {
								value = chrPositionMap.getChrPositionMap().get(Byte.valueOf(curlineParams[0]));
							} else {
								chrPositionMap.getChrPositionMap().put(Byte.valueOf(curlineParams[0]),
																											 value = new HashSet<Integer>());
							}
							value.add(Integer.parseInt(curlineParams[1]));
							value.add(Integer.parseInt(curlineParams[3]));
						}
					}
				}
				LOGGER.info("Dumping ChrPositonMap with: " + curLdFilePath + Files.SERIALIZED_FILE_EXTENSION
										+ "\t " + "Part Count is: " + partCount + "Please wait ...");
				writeChrPositionMap(chrPositionMapList, chrNum, partCount);
				partCount++;

			}
			reader.close();
			bedFileChrMapPartCount.put(curLdFilePath, partCount);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Function to get filepath for {@link ChrPositionMap}
	 *
	 * @param bedFilename the name of the {@link ChrPositionMap} bed file
	 * @param chrNum the chrNum obtained from the ld file
	 * @return the file path to where this {@link ChrPositionMap} should be written
	 */
	private static String getChrPosMapSerFilePath(String bedFilename, int chrNum, int partCount) {
		File theDir = new File(bedDir + BED_FILE_CHR_MAP_FOLDER + File.separator);
		// if the directory does not exist, create it
		if (!theDir.exists()) {
			System.out.println("creating directory: " + bedDir + BED_FILE_CHR_MAP_FOLDER);
			boolean result = theDir.mkdir();

			if (result) {
				System.out.println("DIR created");
			}
		}
		String newFilepath = Files.removeExtention(bedDir + BED_FILE_CHR_MAP_FOLDER + File.separator
																							 + bedFilename)
												 + CHR_MAP_FILE_ID + chrNum + "_part" + partCount + BED_FILE_EXTENTION;
		return Files.getSerializedFilepath(newFilepath);
	}

	/**
	 * Function to write ChrPositionMap for different bedfiles
	 *
	 * @param chrPositionMapList {@link ArrayList} containing
	 *        {@link org.genvisis.common.ChrPositionMap} for all the bedfiles
	 * @param chrNum the chrNum obtained from the ld file
	 */
	private static void writeChrPositionMap(Hashtable<String, ChrPositionMap> chrPositionMapList,
																					int chrNum, int partCount) {
		for (String thisBedFile : chrPositionMapList.keySet()) {
			try {
				chrPositionMapList.get(thisBedFile)
													.writeToFile(getChrPosMapSerFilePath(thisBedFile, chrNum, partCount));
			} catch (IOException e) {
				LOGGER.log(Level.WARNING, "Unable to write ChrPositionMap to file" + e.getMessage(), e);
			}
		}
	}

	/**
	 * Function to write bedFileChrMapPartCount to a file for later use during re-run on the same ld
	 * files
	 *
	 * @param bedFileChrMapPartCount : the bedFileChrMapPartCount to be written
	 */
	private static void writeBedFileChrMapPartCount(Hashtable<String, Integer> bedFileChrMapPartCount) {
		if (bedFileChrMapPartCount != null) {
			try {
				String bedFileChrMapPartCountFilePath = bedDir + BED_FILE_CHR_MAP_FOLDER + File.separator
																								+ BED_FILE_CHR_MAP_PART_COUNT_FILENAME;
				FileWriter fstream = new FileWriter(bedFileChrMapPartCountFilePath, false);
				BufferedWriter out = new BufferedWriter(fstream);
				for (String thisFile : bedFileChrMapPartCount.keySet()) {
					out.write(thisFile);
					out.write(OUTPUT_DELIMITER);
					out.write(String.valueOf(bedFileChrMapPartCount.get(thisFile)));
					out.newLine();
				}
				out.close();
				LOGGER.info("Bed Chr Map File Parts were writtent to file: "
										+ bedFileChrMapPartCountFilePath);
			} catch (IOException e) {
				LOGGER.info("Unable to write to Bed Chr Map File Parts file");
			}
		}
	}

	/**
	 * Function to read the bedFileChrMaps if it exists
	 *
	 * @return a Hashtable of ChrPositionMapPartCount for different files
	 */
	private static Hashtable<String, Integer> readBedFileChrMapPartCount() {
		String curLine;
		String[] curLineParams;
		Hashtable<String, Integer> bedFileChrMapPartCount = new Hashtable<String, Integer>();
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(bedDir + BED_FILE_CHR_MAP_FOLDER + File.separator
																								 + BED_FILE_CHR_MAP_PART_COUNT_FILENAME));
			while (reader.ready()) {
				curLine = reader.readLine();
				curLineParams = curLine.split("\\t");
				bedFileChrMapPartCount.put(curLineParams[0], Integer.parseInt(curLineParams[1]));
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return bedFileChrMapPartCount;
	}

	/**
	 * Function to check if filtered ldFile exists or not and if not then filter and write
	 */
	private static void filterAllLDFiles() {

		String[] ldFilesList = Files.list(ldDir, LD_FILES_EXTENTION, false);

		for (String element : ldFilesList) {
			String lDFilepath = ldDir + element;
			String filteredLDFilepath = Files.removeExtention(lDFilepath) + FILTERED_FILE_EXTENSION;
			if (Files.exists(filteredLDFilepath)) {
				LOGGER.info("Filtered serialized file found for: " + filteredLDFilepath);
			} else {
				LOGGER.info("Filtered serialized file not found. Reading: " + filteredLDFilepath);
				// read the segments list
				LOGGER.info("Filtering and writing serialized file: " + filteredLDFilepath);
				writeFilteredLDFile(lDFilepath);
				LOGGER.info("Completed writing: " + filteredLDFilepath);
			}
		}
	}

	/**
	 * Function to filter and write a ldFile data
	 *
	 * @param filename the name of the file which has to be filtered
	 */
	private static void writeFilteredLDFile(String filename) {
		String curLine;
		String[] curLineParams;
		try {
			BufferedReader inFile = new BufferedReader(new FileReader(filename));
			BufferedWriter outFile = new BufferedWriter(new FileWriter(Files.removeExtention(filename)
																																 + FILTERED_FILE_EXTENSION, false));
			inFile.readLine(); // skill the first header line

			while (inFile.ready()) {
				curLine = inFile.readLine();
				curLineParams = WHILTE_SPACE_PATTERN.split(curLine.trim());

				if (ext.isValidDouble(curLineParams[6]) && Double.parseDouble(curLineParams[6]) > 0.5) {
					outFile.write(curLineParams[0]);
					outFile.write(OUTPUT_DELIMITER);
					outFile.write(curLineParams[1]);
					outFile.write(OUTPUT_DELIMITER);
					outFile.write(curLineParams[2]);
					outFile.write(OUTPUT_DELIMITER);
					outFile.write(curLineParams[4]);
					outFile.write(OUTPUT_DELIMITER);
					outFile.write(curLineParams[5]);
					outFile.newLine();
				}
			}
			outFile.close();
			inFile.close();
		} catch (FileNotFoundException e) {
			LOGGER.warning("Not able to find file: " + filename + e.getMessage() + e.getStackTrace());
		} catch (IOException e) {
			LOGGER.warning("IOException while serialzing ld file: " + filename + e.getMessage()
										 + e.getStackTrace());
		}
	}

	/**
	 * Function to find overlaps from a given pvalues list and segment list
	 *
	 * @param segs : a {@link Segment}[][] containing all the segments read from the file arranged
	 *        with chr
	 * @param pValueRecords : {@link ArrayList} of {@link PValueFileFormat} containing all the
	 *        p-values records
	 * @return a Map<Integer, Map<String, Integer>> containing stats for different bin for this
	 *         segment file
	 */
	private static TreeMap<Integer, Map<String, Long>> countOverlaps(Segment[][] segs,
																																	 List<PValueFileFormat> pValueRecords,
																																	 String bedFilename,
																																	 Hashtable<String, Integer> bedFileChrMapPartCount) {
		boolean insideRegion;
		int pValueBin;
		Map<String, Long> value;
		TreeMap<Integer, Map<String, Long>> overlapStats = new TreeMap<Integer, Map<String, Long>>();
		String[] ldFilesList = Files.list(ldDir, LD_FILES_EXTENTION, false);
		ChrPositionMap chrPositionMap = null;

		if (performLD) {
			// get this bed file chr position map
			ChrPositionMap chrPositionMapTemp = new ChrPositionMap();
			chrPositionMap = new ChrPositionMap();
			try {
				for (int i = 0; i < ldFilesList.length; i++) {
					int partCount = bedFileChrMapPartCount.get(ldFilesList[i]);
					for (int j = 0; j < partCount; j++) {
						chrPositionMapTemp.readFromFile(getChrPosMapSerFilePath(bedFilename, i, j));
						chrPositionMap.getChrPositionMap().putAll(chrPositionMapTemp.getChrPositionMap());
					}
				}

			} catch (IOException e) {
				LOGGER.log(Level.WARNING,
									 "Unable to read the serialized chr position map file: " + e.getMessage(), e);
			} catch (ClassNotFoundException e) {
				LOGGER.log(Level.WARNING, "Unable to find chr position map class" + e.getMessage(), e);
			}
		}

		for (PValueFileFormat curRecord : pValueRecords) {
			if (segs[curRecord.chr] != null) {
				// the p-value bin (-LOG10(pvalue) rounded down to nearest integer
				pValueBin = (-(int) Math.floor(Math.log10(curRecord.pValue)));
				// keep a count of maximum bin size encountered so far
				if (maxBinSize < pValueBin) {
					maxBinSize = pValueBin;
				}
				// if we already have values for this bin
				if (overlapStats.containsKey(pValueBin)) {
					// get the existing value
					value = overlapStats.get(pValueBin);
				} else {
					// else create new
					overlapStats.put(pValueBin, value = new HashMap<String, Long>());
				}
				// if there are segments for this chr

				// find the overlap using binary search
				insideRegion = Segment.binarySearchForOverlap(curRecord.seg, segs[curRecord.chr]) != -1;

				if (performLD) {
					// if not inside region then check in chr position map
					if (!insideRegion) {
						if (chrPositionMap.getChrPositionMap().containsKey(curRecord.seg.getChr())) {
							HashSet<Integer> thisBedChrPosHashSet = chrPositionMap.getChrPositionMap()
																																		.get(curRecord.seg.getChr());
							if (thisBedChrPosHashSet.contains(curRecord.seg.getStart())) {
								insideRegion = true;
							}
						}
					}
				}

				if (value.containsKey(TOTAL_MARKERS)) {
					value.put(TOTAL_MARKERS, value.get(TOTAL_MARKERS) + 1);
				} else {
					value.put(TOTAL_MARKERS, 1L);
				}

				// if it is inside the region
				if (insideRegion) {
					// if there is already a count then increase it else initialize with 1
					if (value.containsKey(INSIDE_REGION)) {
						value.put(INSIDE_REGION, (value.get(INSIDE_REGION)) + 1);
					} else {
						value.put(INSIDE_REGION, 1L);
					}
				}
			}
		} // end of while
		return overlapStats;
	}

	/**
	 * Function to read the p-value into memory
	 *
	 * @param pValueFilepath : the path to the p-value file
	 * @return a {@link ArrayList} of {@link PValueFileFormat} containing all the pvalue records
	 */
	private static ArrayList<PValueFileFormat> readPValueFile(String pValueFilepath) {
		int[] headerIndices = findHeaderIndices(pValueFilepath);
		String curLine;
		String[] curLineParams;
		byte curChr;
		Segment curSeg;
		double curPValue;

		ArrayList<PValueFileFormat> result = new ArrayList<PValueFileFormat>();

		try {
			BufferedReader reader = new BufferedReader(new FileReader(pValueFilepath));
			// skip the headers in file
			curLine = reader.readLine();
			while (reader.ready()) {
				curLine = reader.readLine();
				curLineParams = WHILTE_SPACE_PATTERN.split(curLine.trim());
				// if p-value is a valid number. This also filters records containing pvalue "NA"
				if (ext.isValidDouble(curLineParams[headerIndices[3]])) {
					// chr value
					curChr = Positions.chromosomeNumber(curLineParams[headerIndices[1]]);
					// current segment
					curSeg = new Segment(curChr, Integer.parseInt(curLineParams[headerIndices[2]]),
															 Integer.parseInt(curLineParams[headerIndices[2]]));
					curPValue = Double.parseDouble(curLineParams[headerIndices[3]]);
					// add this record
					result.add(new PValueFileFormat(curChr, curSeg, curPValue));
				}
			}
			reader.close();
		} catch (IOException e) {
			LOGGER.info("Unable to read p-value file" + e.getMessage());
			System.exit(FAILURE);
		}
		Collections.sort(result); // sort on the basis of chr
		return result;
	}

	/**
	 * Function to get all the segments from the file as a {@link Segment}[][] Gets the segments from
	 * the serialized file if exists else reads the main segmentFile gets the segments and writes a
	 * serialized output for future use.
	 *
	 * @param segmentFile : the filename
	 * @return a {@link Segment}[][] containing all the segments read from the file arranged with chr
	 */
	private static Segment[][] getSegments(String segmentFile) {
		Segment[][] segs;
		String serializedSegFile = Files.getSerializedFilepath(segmentFile);
		if (Files.exists(serializedSegFile)) {
			LOGGER.info("Serialized file found. Reading: " + serializedSegFile);
			// read the already existing serialized segment file
			segs = SegmentLists.load(serializedSegFile, false).getLists();

		} else {
			LOGGER.info("Serialized file not found. Reading: " + segmentFile);
			// read the segments list
			segs = SegmentLists.parseSegmentList(segmentFile, 0, 1, 2, false).getLists();
			// serialize and write for future use
			LOGGER.info("Creating serialized file. Writing: " + serializedSegFile);
			new SegmentLists(segs).serialize(serializedSegFile);
		}
		return segs;
	}

	/**
	 * Function to get the indices of required headers from the file
	 *
	 * @param filename : the name of the file
	 * @return a int[] which contains the indices of required headers
	 */
	private static int[] findHeaderIndices(String filename) {
		String[] headers;
		// get the headers for pvalues file
		headers = Files.getHeaderOfFile(filename, new org.genvisis.common.Logger(null));
		// get header indices
		return (ext.indexFactors(new String[][] {P_VALUE_FILE_HEADERS[0], P_VALUE_FILE_HEADERS[1],
																						 P_VALUE_FILE_HEADERS[2], P_VALUE_FILE_HEADERS[3]},
														 headers, false, true, true, false));
	}

	private static class OutputFileFormat implements Serializable {
		private static final long serialVersionUID = 1L;

		String file;
		double[][] ratio;

		public OutputFileFormat(String file, double[][] ratio) {
			this.file = file;
			this.ratio = ratio;
		}

		@Override
		public String toString() {
			StringBuilder result = new StringBuilder();
			result.append(file);
			result.append(Arrays.toString(ratio));
			return result.toString();

		}
	}

	private static class PValueFileFormat implements Serializable, Comparable<PValueFileFormat> {
		private static final long serialVersionUID = 1L;

		byte chr;
		Segment seg;
		double pValue;

		private PValueFileFormat(byte chr, Segment seg, double pValue) {
			this.chr = chr;
			this.seg = seg;
			this.pValue = pValue;
		}

		@Override
		public int compareTo(PValueFileFormat o) {
			byte compareQuantity = o.chr;

			// ascending order
			return chr - compareQuantity;
		}
	}

	private static class PlinkFile implements Serializable {
		private static final long serialVersionUID = 1L;

		String[] markerNames;
		byte[] chrs;
		int[] positions;

		private PlinkFile(String[] markerNames, byte[] chrs, int[] positions) {
			this.markerNames = markerNames;
			this.chrs = chrs;
			this.positions = positions;
		}
	}
}
