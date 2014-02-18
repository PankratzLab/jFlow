package gwas;

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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import common.Array;
import common.ChrPositionMap;
import common.Files;
import common.Positions;
import common.ext;
import filesys.Segment;
import filesys.SegmentLists;
import filesys.SnpMarkerSet;

/**
 * Program to summarize p-values from a Genome-Wide Association Study (GWAS) and bins marker depending on their p-values and depending upon whether they are inside or outside of a DNAse hypersensitive (DNAseHS) region. These segments are loaded from BED files and the name of the celltype is taken from the first part of the filename. The p-value column is not required to be in a consistent column number. As a result the program creates a file with the following tab-delimited columns: celltype, bin1(0.1<pvalue<1.0), bin2(0.01<pvalue<0.1), bin3(0.001<pvalue<0.01)... Then this output is plotted using LinePlot
 * 
 * Author: Rohit Sinha Date: 10/24/13 Time: 2:48 AM
 */
public class DnaseEnrichment {

	private final static String[][] P_VALUE_FILE_HEADERS = { { "MarkerName", "Marker", "SNP" }, { "Chr" }, { "Position", "Pos" }, { "P-value", "p-val" } };
	private final static String INSIDE_REGION = "insideRegion";
	private final static String TOTAL_MARKERS = "totalMarkers";
	private final static int FAILURE = 1;
	private final static String BED_FILE_EXTENTION = ".bed";
	private final static String LD_FILES_EXTENTION = ".ld";
	private final static String OUTPUT_DELIMITER = "\t";
	private final static String DATA_SEPARATOR = ":";
	private final static Logger LOGGER = Logger.getLogger(DnaseEnrichment.class.getName()); // logger for this class
	private final static String OUTPUT_FILENAME = "DnaseEnrichment.xln";
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
	private static String bed_dir = null;
	private static String[] bedFileList;
	private static int ldLines;
	private static int numThreads = 1; // number of threads to run. default is 1

	private static ArrayList<HashSet<String>> dhsregionsHashSetList; // global static variable to hold DHS regions
	private static Hashtable<String, Integer> bedFileChrMapPartCount = new Hashtable<String, Integer>(); // global static variable for ChrPositionMap

	/**
	 * WorkerThreads which process different LD files and create ChrPositionMap concurrently
	 */
	public class WorkerThread implements Runnable {
		private String ldFilePath; // path of the ld file
		private int chrNum; // chr number

		public WorkerThread(String s, int chrNum) {
			this.ldFilePath = s;
			this.chrNum = chrNum;
		}

		@Override
		public void run() {
			LOGGER.info(Thread.currentThread().getName() + ": Building ChrPositionMap for chr file: " + ldFilePath);
			generateChrPositionMap(ldFilePath, chrNum); // generate the chrPositionMap for all the bed files using this chr file
			LOGGER.info(Thread.currentThread().getName() + ": Completed building ChrPositionMap for file: " + ldFilePath);

		}
	}

	/**
	 * Function to print the command line usage.
	 * 
	 * @return a {@link String}: the command line usage help
	 */
	private static String showCommandLineHelp() {

		return ("\n" + "gwas.DnaseEnrichment requires six arguments\n" + "\t(1) directory with bed files (bedDir)\n" + "\t(2) filename with chr/pos/pvals (pValueFile)\n" + "\t(3) The bim file (pLinkFile)\n" + "\t(4) " + "ld files directory (ldDir)\n" + "\t(5)number of lines on which ld files will be split (ldLines)\n" + "\t(6)number of threads (numThreads: defaults=1)\n" + "");

	}

	public static void main(String[] args) {

		ArrayList<OutputFileFormat> overlapStats;
		int numArgs = args.length;
		String filename = null;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(showCommandLineHelp());
				System.exit(1);
			} else if (args[i].toLowerCase().startsWith("beddir=")) {
				bed_dir = ext.verifyDirFormat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("pvaluefile=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("plinkfile=")) {
				plinkFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("lddir=")) {
				ldDir = ext.verifyDirFormat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("ldlines=")) {
				ldLines = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("numthreads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0 && numArgs != 1) {
			System.err.println(showCommandLineHelp());
			System.exit(1);
		}

		// get names of all the bed files
		bedFileList = Files.list(bed_dir, BED_FILE_EXTENTION, false);
		Arrays.sort(bedFileList);

		// Print out parameters to user
		LOGGER.info("Starting DnaseEnrichment Processing");
		LOGGER.info("LD Directory is: " + ldDir);
		LOGGER.info("Bed Directory is: " + bed_dir);
		LOGGER.info("Plink File location is: " + plinkFile);
		LOGGER.info("LD Line Count is : " + ldLines);
		LOGGER.info("Num of threads to use while making ChrPositionMap: " + numThreads);

		// Filter ld files by storing only the new required columns
		LOGGER.info("Filtering LD files. Please wait...");
		filterAllLDFiles();

		// build ChrPositionMap using the ld files for each bed file
		LOGGER.info("Building ChrPositionMap. Please wait...");

		if (Files.exists(bed_dir + BED_FILE_CHR_MAP_FOLDER + File.separator + BED_FILE_CHR_MAP_PART_COUNT_FILENAME)) {
			System.out.println("We found ChrPositionMap for bed files. Do you want to continue with these files (y/n)" + "" + "" + "?");
			System.out.println("Note: If you feel like you ld files or bed files have changed from the last run then " + "" + "we strongly suggest you to stop here and delete the ChrMapPosition directory " + "and " + "run " + "the program" + " " + "again" + "" + "" + "" + "" + "");
			System.out.println("Chr Position Map File Directory: " + bed_dir + BED_FILE_CHR_MAP_FOLDER + File.separator);
			System.out.print("DO you to continue Y/N: ");
			Scanner scan = new Scanner(System.in); // scanner for input
			String userChoice = scan.nextLine();
			scan.close();
			if (userChoice.toLowerCase().equals("y")) {
				bedFileChrMapPartCount = readBedFileChrMapPartCount();
			} else if (userChoice.toLowerCase().equals("n")) {
				System.out.println("Exiting on user selection");
				System.exit(0);
			} else {
				System.out.println("Invalid choice. Exiting.");
				System.exit(FAILURE);
			}
		} else {
			dhsregionsHashSetList = findDHSRegionMarkers();
			DnaseEnrichment dnaseEnrichmentObject = new DnaseEnrichment();
			dnaseEnrichmentObject.runWorkers();
			writeBedFileChrMapPartCount(bedFileChrMapPartCount);
		}

		// find the overlapping position counts
		LOGGER.info("Building statistics. Please wait...");
		overlapStats = findOverlapRegions(bed_dir, filename, bedFileChrMapPartCount);

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
		while (!executor.isTerminated())
			;
		LOGGER.info("All threads terminated successfully: Processed all chr files to build ChrPositionMap");
	}

	/**
	 * Function to write overlap output to output file.
	 * 
	 * @param overlapStats
	 *            : an {@link ArrayList} of {@link OutputFileFormat} which is the overlap statistics
	 */
	private static void writeOutputFile(ArrayList<OutputFileFormat> overlapStats) {
		FileWriter fstream;
		BufferedWriter out;

		List<List<Object>> resultArrayList = new ArrayList<List<Object>>();

		for (OutputFileFormat curRecord : overlapStats) {
			ArrayList<Object> resultArray = new ArrayList<Object>();
			resultArray.add(curRecord.file);

			for (int i = 0; i < curRecord.ratio.length; i++) {
				// to avoid divide by 0 error
				double ratio = curRecord.ratio[i][0] / (curRecord.ratio[i][1] != 0 ? curRecord.ratio[i][1] : 1);
				// write all the ratio
				resultArray.add(String.valueOf(ratio) + DATA_SEPARATOR + String.valueOf(curRecord.ratio[i][0]) + "/" + String.valueOf(curRecord.ratio[i][1]));
			}
			resultArrayList.add(resultArray);
		}
		resultArrayList = Array.transpose(resultArrayList);
		try {
			String outputFilePath = bed_dir + OUTPUT_FOLDER + File.separator + OUTPUT_FILENAME;

			File theDir = new File(bed_dir + OUTPUT_FOLDER + File.separator);
			// if the directory does not exist, create it
			if (!theDir.exists()) {
				System.out.println("creating directory: " + bed_dir + OUTPUT_FOLDER + File.separator);
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
	 * @param dir
	 *            : the directory path containing the segments file
	 * @param filename
	 *            : the filename of the p-value file
	 * @return a Map<String, Map<Integer, Map<String, Integer>>> which shows the inside and outside the region stats for all the segments file
	 */
	private static ArrayList<OutputFileFormat> findOverlapRegions(String dir, String filename, Hashtable<String, Integer> bedFileChrMapPartCount) {

		ArrayList<OutputFileFormat> thisFileOutput = new ArrayList<OutputFileFormat>();

		// read the pvalues records from the pvalue file in memory
		ArrayList<PValueFileFormat> pValueRecords = readPValueFile(filename);

		LOGGER.info("Starting to process segment files...");

		for (int i = 0; i < bedFileList.length; i++) {
			LOGGER.info("Processing: " + bedFileList[i]);
			Segment[][] segs = getSegments(dir + bedFileList[i]);
			TreeMap<Integer, Map<String, Long>> overlapStats = countOverlaps(segs, pValueRecords, bedFileList[i], bedFileChrMapPartCount);
			double[][] ratioList = new double[maxBinSize + 1][2]; // array for holding numerator and denominator
			double cumBinMarkers = 0;
			double cumBinMarkersInside = 0;
			int key;
			NavigableSet<Integer> keySet = overlapStats.descendingKeySet();
			// get the numerator for all the bins cumulatively
			for (Iterator<Integer> iter = keySet.iterator(); iter.hasNext();) {
				key = iter.next();
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
			thisFileOutput.add(new OutputFileFormat(bedFileList[i], ratioList));
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
			SnpMarkerSet markerSet = new SnpMarkerSet(plinkFile, false, new common.Logger(null));
			plinkFileContents = new PlinkFile(markerSet.getMarkerNames(), markerSet.getChrs(), markerSet.getPositions());
		} else {
			LOGGER.warning("Error - could not find plink.bim; required to define segments" + " " + "and perform a crucial " + "datacheck");
		}

		ArrayList<Segment> plinkSegments = createSegmentsFromPlinkFile(plinkFileContents);
		LOGGER.info("Starting to find DHS regions for bed files");
		for (int i = 0; i < bedFileList.length; i++) {
			LOGGER.info("Processing for DHS region: " + bedFileList[i]);
			Segment[][] segs = getSegments(bed_dir + bedFileList[i]);

			HashSet<String> ldMarkerHashSet = new HashSet<String>();

			for (int j = 0; j < plinkSegments.size(); j++) {
				// if there is a overlap
				if (Segment.binarySearchForOverlap(plinkSegments.get(j), segs[plinkSegments.get(j).getChr()]) != -1) {
					ldMarkerHashSet.add(plinkFileContents.markerNames[j]);
				}
			}
			dhsregionsHashSetList.add(ldMarkerHashSet);
		}
		return dhsregionsHashSetList;
	}

	/**
	 * Function to create segments from {@link PlinkFile}
	 * 
	 * @param plinkFileContents
	 *            {@link PlinkFile} containing plink file contents
	 * @return {@link ArrayList} of {@link Segment}
	 */
	private static ArrayList<Segment> createSegmentsFromPlinkFile(PlinkFile plinkFileContents) {
		// get chr, position and marker names from the plink file
		LOGGER.info("Reading plink file. This might take couple of seconds. Please wait ...");
		ArrayList<Segment> plinkSegments = new ArrayList<Segment>();

		// if something was wrong in the plink file data fields inform the user
		if (plinkFileContents.markerNames.length != plinkFileContents.chrs.length && plinkFileContents.chrs.length == plinkFileContents.positions.length) {
			LOGGER.log(Level.WARNING, "size of markers, chr, positions are not equal in the plink file");
		}

		// if reading plink file was successful then create segments out of it
		if (plinkFileContents != null) {
			for (int i = 0; i < plinkFileContents.markerNames.length; i++) {
				// create a segment from this line
				Segment curSeg = new Segment(plinkFileContents.chrs[i], plinkFileContents.positions[i], plinkFileContents.positions[i]);
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
			BufferedReader reader = new BufferedReader(new FileReader(ldDir + Files.removeExtention(curLdFilePath) + FILTERED_FILE_EXTENSION));
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
						if (dhsregionsHashSetList.get(j).contains(curlineParams[2]) || dhsregionsHashSetList.get(j).contains(curlineParams[4])) {
							if (chrPositionMap.getChrPositionMap().containsKey(Byte.valueOf(curlineParams[0]))) {
								value = chrPositionMap.getChrPositionMap().get(Byte.valueOf(curlineParams[0]));
							} else {
								chrPositionMap.getChrPositionMap().put(Byte.valueOf(curlineParams[0]), value = new HashSet<Integer>());
							}
							value.add(Integer.parseInt(curlineParams[1]));
							value.add(Integer.parseInt(curlineParams[3]));
						}
					}
				}
				LOGGER.info("Dumping ChrPositonMap with: " + curLdFilePath + Files.SERIALIZED_FILE_EXTENSION + "\t " + "Part Count is: " + partCount + "Please wait ...");
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
	 * @param bedFilename
	 *            the name of the {@link ChrPositionMap} bed file
	 * @param chrNum
	 *            the chrNum obtained from the ld file
	 * @return the file path to where this {@link ChrPositionMap} should be written
	 */
	private static String getChrPosMapSerFilePath(String bedFilename, int chrNum, int partCount) {
		File theDir = new File(bed_dir + BED_FILE_CHR_MAP_FOLDER + File.separator);
		// if the directory does not exist, create it
		if (!theDir.exists()) {
			System.out.println("creating directory: " + bed_dir + BED_FILE_CHR_MAP_FOLDER);
			boolean result = theDir.mkdir();

			if (result) {
				System.out.println("DIR created");
			}
		}
		String newFilepath = Files.removeExtention(bed_dir + BED_FILE_CHR_MAP_FOLDER + File.separator + bedFilename) + CHR_MAP_FILE_ID + chrNum + "_part" + partCount + BED_FILE_EXTENTION;
		return Files.getSerializedFilepath(newFilepath);
	}

	/**
	 * Function to write ChrPositionMap for different bedfiles
	 * 
	 * @param chrPositionMapList
	 *            {@link ArrayList} containing {@link common.ChrPositionMap} for all the bedfiles
	 * @param chrNum
	 *            the chrNum obtained from the ld file
	 */
	private static void writeChrPositionMap(Hashtable<String, ChrPositionMap> chrPositionMapList, int chrNum, int partCount) {
		for (String thisBedFile : chrPositionMapList.keySet()) {
			try {
				chrPositionMapList.get(thisBedFile).writeToFile(getChrPosMapSerFilePath(thisBedFile, chrNum, partCount));
			} catch (IOException e) {
				LOGGER.log(Level.WARNING, "Unable to write ChrPositionMap to file" + e.getMessage(), e);
			}
		}
	}

	/**
	 * Function to write bedFileChrMapPartCount to a file for later use during re-run on the same ld files
	 * 
	 * @param bedFileChrMapPartCount
	 *            : the bedFileChrMapPartCount to be written
	 */
	private static void writeBedFileChrMapPartCount(Hashtable<String, Integer> bedFileChrMapPartCount) {
		if (bedFileChrMapPartCount != null) {
			try {
				String bedFileChrMapPartCountFilePath = bed_dir + BED_FILE_CHR_MAP_FOLDER + File.separator + BED_FILE_CHR_MAP_PART_COUNT_FILENAME;
				FileWriter fstream = new FileWriter(bedFileChrMapPartCountFilePath, false);
				BufferedWriter out = new BufferedWriter(fstream);
				for (String thisFile : bedFileChrMapPartCount.keySet()) {
					out.write(thisFile);
					out.write(OUTPUT_DELIMITER);
					out.write(String.valueOf(bedFileChrMapPartCount.get(thisFile)));
					out.newLine();
				}
				out.close();
				LOGGER.info("Bed Chr Map File Parts were writtent to file: " + bedFileChrMapPartCountFilePath);
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
			reader = new BufferedReader(new FileReader(bed_dir + BED_FILE_CHR_MAP_FOLDER + File.separator + BED_FILE_CHR_MAP_PART_COUNT_FILENAME));
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

		for (int i = 0; i < ldFilesList.length; i++) {
			String lDFilepath = ldDir + ldFilesList[i];
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
	 * @param filename
	 *            the name of the file which has to be filtered
	 */
	private static void writeFilteredLDFile(String filename) {
		String curLine;
		String[] curLineParams;
		try {
			BufferedReader inFile = new BufferedReader(new FileReader(filename));
			BufferedWriter outFile = new BufferedWriter(new FileWriter(Files.removeExtention(filename) + FILTERED_FILE_EXTENSION, false));
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
			LOGGER.warning("IOException while serialzing ld file: " + filename + e.getMessage() + e.getStackTrace());
		}
	}

	/**
	 * Function to find overlaps from a given pvalues list and segment list
	 * 
	 * @param segs
	 *            : a {@link Segment}[][] containing all the segments read from the file arranged with chr
	 * @param pValueRecords
	 *            : {@link ArrayList} of {@link PValueFileFormat} containing all the p-values records
	 * @return a Map<Integer, Map<String, Integer>> containing stats for different bin for this segment file
	 */
	private static TreeMap<Integer, Map<String, Long>> countOverlaps(Segment[][] segs, ArrayList<PValueFileFormat> pValueRecords, String bedFilename, Hashtable<String, Integer> bedFileChrMapPartCount) {
		boolean insideRegion;
		int pValueBin;
		Map<String, Long> value;
		TreeMap<Integer, Map<String, Long>> overlapStats = new TreeMap<Integer, Map<String, Long>>();
		String[] ldFilesList = Files.list(ldDir, LD_FILES_EXTENTION, false);

		// get this bed file chr position map
		ChrPositionMap chrPositionMapTemp = new ChrPositionMap();
		ChrPositionMap chrPositionMap = new ChrPositionMap();
		try {
			for (int i = 0; i < ldFilesList.length; i++) {
				int partCount = bedFileChrMapPartCount.get(ldFilesList[i]);
				for (int j = 0; j < partCount; j++) {
					chrPositionMapTemp.readFromFile(getChrPosMapSerFilePath(bedFilename, i, j));
					chrPositionMap.getChrPositionMap().putAll(chrPositionMapTemp.getChrPositionMap());
				}
			}

		} catch (IOException e) {
			LOGGER.log(Level.WARNING, "Unable to read the serialized chr position map file: " + e.getMessage(), e);
		} catch (ClassNotFoundException e) {
			LOGGER.log(Level.WARNING, "Unable to find chr position map class" + e.getMessage(), e);
		}

		for (PValueFileFormat curRecord : pValueRecords) {
			if (segs[curRecord.chr] != null) {
				// the p-value bin (-LOG10(pvalue) rounded down to nearest integer
				pValueBin = (-(int) Math.floor(Math.log10(curRecord.pValue)));
				// keep a count of maximum bin size encountered so far
				if (maxBinSize < pValueBin)
					maxBinSize = pValueBin;
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

				// if not inside region then check in chr position map
				if (!insideRegion) {
					if (chrPositionMap.getChrPositionMap().containsKey(curRecord.seg.getChr())) {
						HashSet<Integer> thisBedChrPosHashSet = chrPositionMap.getChrPositionMap().get(curRecord.seg.getChr());
						if (thisBedChrPosHashSet.contains(curRecord.seg.getStart())) {
							insideRegion = true;
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
	 * @param pValueFilepath
	 *            : the path to the p-value file
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
					curSeg = new Segment(curChr, Integer.parseInt(curLineParams[headerIndices[2]]), Integer.parseInt(curLineParams[headerIndices[2]]));
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
	 * Function to get all the segments from the file as a {@link Segment}[][] Gets the segments from the serialized file if exists else reads the main segmentFile gets the segments and writes a serialized output for future use.
	 * 
	 * @param segmentFile
	 *            : the filename
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
	 * @param filename
	 *            : the name of the file
	 * @return a int[] which contains the indices of required headers
	 */
	private static int[] findHeaderIndices(String filename) {
		String[] headers;
		// get the headers for pvalues file
		headers = Files.getHeaderOfFile(filename, new common.Logger(null));
		// get header indices
		return (ext.indexFactors(new String[][] { P_VALUE_FILE_HEADERS[0], P_VALUE_FILE_HEADERS[1], P_VALUE_FILE_HEADERS[2], P_VALUE_FILE_HEADERS[3] }, headers, false, true, true, false));
	}

	private static class OutputFileFormat implements Serializable {
		private static final long serialVersionUID = 1L;

		String file;
		double[][] ratio;

		public OutputFileFormat(String file, double[][] ratio) {
			this.file = file;
			this.ratio = ratio;
		}

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
			return this.chr - compareQuantity;
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

	// This part of code was developed earlier but due to change in requirements we are not using it anymore.
	// /**
	// * generateChrPositionMap function which uses HeapStats to see when the free memory is low and then dumps the
	// * ChrPositionMap to disk to save memory.
	// */
	// private static void generateChrPositionMap() {
	// String[] ldFilesList = Files.list(ldDir, LD_FILES_EXTENTION, false);
	// ChrPositionMap chrPositionMap;
	// HashSet<Integer> value;
	// Hashtable<String, ChrPositionMap> chrPositionMapList = new Hashtable<String, ChrPositionMap>();
	// ArrayList<HashSet<String>> dhsregionsHashSetList = findDHSRegionMarkers();
	// int count = 0;
	//
	// for (int i = 0; i < ldFilesList.length; i++) {
	// LOGGER.info("Processing: " + ldFilesList[i] + Files.SERIALIZED_FILE_EXTENSION);
	// try {
	// BufferedReader reader = new BufferedReader(new FileReader(ldDir + ldFilesList[i] + Files.SERIALIZED_FILE_EXTENSION));
	// reader.readLine(); // skill the first header line
	//
	// while (reader.ready()) {
	// count++;
	// if (count > 100000) {
	// LOGGER.info("Finding ChrPositonMap with: " + ldFilesList[i] + Files.SERIALIZED_FILE_EXTENSION + "\t Please wait ...");
	// count = 0;
	// }
	// String curline = reader.readLine();
	// String[] curlineParams = WHILTE_SPACE_PATTERN.split(curline.trim());
	// if (heapStats.getFreeMemory() < 100) {
	// LOGGER.info("Memory getting low. Dumping ChrPositionMap to disk." + "Free: " + heapStats.getFreeMemory());
	//
	// writeChrPositionMap(chrPositionMapList, CHR_POSITION_FILE_PART_COUNT++);
	// for (String thisKey : chrPositionMapList.keySet()) {
	// ChrPositionMap thisObject = chrPositionMapList.get(thisKey);
	// thisObject.clearAll();
	// }
	// chrPositionMapList.clear();
	// LOGGER.info("free before sleep. Dumping ChrPositionMap to disk." + "Free: " + heapStats.getFreeMemory());
	// try {
	// Thread.sleep(5000);
	// } catch (InterruptedException e) {
	// e.printStackTrace();
	// }
	// chrPositionMapList = new Hashtable<String, ChrPositionMap>();
	// LOGGER.info("Freed mem. Dumping ChrPositionMap to disk." + "Free: " + heapStats.getFreeMemory());
	// }
	// for (int j = 0; j < dhsregionsHashSetList.size(); j++) {
	//
	// if (chrPositionMapList.containsKey(bedFileList[j])) {
	// chrPositionMap = chrPositionMapList.get(bedFileList[j]);
	// } else {
	// chrPositionMapList.put(bedFileList[j], chrPositionMap = new ChrPositionMap());
	//
	// }
	// if (dhsregionsHashSetList.get(j).contains(curlineParams[2]) || dhsregionsHashSetList.get(j).contains(curlineParams[4])) {
	// if (chrPositionMap.getChrPositionMap().containsKey(Byte.valueOf(curlineParams[0]))) {
	// value = chrPositionMap.getChrPositionMap().get(Byte.valueOf(curlineParams[0]));
	// } else {
	// chrPositionMap.getChrPositionMap().put(Byte.valueOf(curlineParams[0]), value = new HashSet<Integer>());
	// }
	// value.add(Integer.parseInt(curlineParams[1]));
	// value.add(Integer.parseInt(curlineParams[3]));
	// }
	// }
	// }
	// reader.close();
	// } catch (IOException e) {
	// e.printStackTrace();
	// }
	// }
	// writeChrPositionMap(chrPositionMapList, CHR_POSITION_FILE_PART_COUNT++);
	// LOGGER.info("Total Chr parts: " + CHR_POSITION_FILE_PART_COUNT);
	// }

} // end of DnaseEnrichment class