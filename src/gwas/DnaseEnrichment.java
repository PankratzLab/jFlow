package gwas;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import common.Array;
import common.Files;
import common.Positions;
import common.ext;

import filesys.Segment;
import filesys.SegmentLists;

/**
 * Program to summarize p-values from a Genome-Wide Association Study (GWAS) and bins marker depending on their p-values and depending upon whether they are inside or outside of a DNAse hypersensitive (DNAseHS) region. These segments are loaded from BED files and the name of the celltype is taken from the first part of the filename. The p-value column is not required to be in a consistent column number. As a result the program creates a file with the following tab-delimited columns: celltype, bin1(0.1<pvalue<1.0), bin2(0.01<pvalue<0.1), bin3(0.001<pvalue<0.01)... Then this output is plotted using LinePlot
 * <p/>
 * Author: Rohit Sinha Date: 10/24/13 Time: 2:48 AM
 */
public class DnaseEnrichment {

	private final static String[][] P_VALUE_FILE_HEADERS = { { "MarkerName", "Marker", "SNP" }, { "Chr" }, { "Position", "Pos" }, { "P-value", "p-val" } };
	private final static String INSIDE_REGION = "insideRegion";
	private final static String OUTSIDE_REGION = "outsideRegion";
	private final static int FAILURE = 1;
	private final static String DNAseHS_FILE_EXTENSION = ".bed";
	private final static String OUTPUT_DELIMITER = "\t";
	private final static String DATA_SEPARATOR = ":";
	private final static Logger LOGGER = Logger.getLogger(DnaseEnrichment.class.getName()); // logger for this class
	private final static String OUTPUT_FILENAME = "DnaseEnrichment.xln";
	private static int maxBinSize = 0;

	/**
	 * Function to print the command line usage.
	 * 
	 * @return a {@link String}: the command line usage help
	 */
	private static String showCommandLineHelp() {

		String usage = "\n" + "gwas.DnaseEnrichment requires two arguments\n" + "   (1) directory with bed files\n" + "   (2) filename with chr/pos/pvals\n" + "";
		return usage;
	}

	public static void main(String[] args) {

		ArrayList<OutputFileFormat> overlapStats;
		int numArgs = args.length;
		String bed_dir = null;
		String filename = null;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(showCommandLineHelp());
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				bed_dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(showCommandLineHelp());
			System.exit(1);
		}

		// find the overlapping position counts
		overlapStats = findOverlapRegions(bed_dir, filename);

		// Overlap statistics
		LOGGER.info(overlapStats.toString());

		// write to output file
		writeOutputFile(overlapStats);

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
			String outputFilePath = System.getProperty("user.dir") + File.separator + OUTPUT_FILENAME;
			fstream = new FileWriter(outputFilePath, false);
			out = new BufferedWriter(fstream);
			for (List<Object> curRecord : resultArrayList){
				for(Object curString : curRecord){
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
	private static ArrayList<OutputFileFormat> findOverlapRegions(String dir, String filename) {

		ArrayList<OutputFileFormat> thisFileOutput = new ArrayList<OutputFileFormat>();
		String[] filesList = Files.list(dir, DNAseHS_FILE_EXTENSION, false);

		// read the pvalues records from the pvalue file in memory
		ArrayList<PValueFileFormat> pValueRecords = readPValueFile(filename);

		LOGGER.info("Starting to process segment files...");
		for (int i = 0; i < filesList.length; i++) {
			LOGGER.info("Processing: " + filesList[i]);
			Segment[][] segs = getSegments(dir + filesList[i]);
			Map<Integer, Map<String, Integer>> overlapStats = countOverlaps(segs, pValueRecords);
			double[][] ratioList = new double[maxBinSize + 1][2]; // array for holding numerator and denominator
			for (Map.Entry<Integer, Map<String, Integer>> entry : overlapStats.entrySet()) {
				int key = entry.getKey();
				Map<String, Integer> value = entry.getValue();

				double numerator = value.containsKey(INSIDE_REGION) ? value.get(INSIDE_REGION) : 0;
				double denominator = value.containsKey(OUTSIDE_REGION) ? value.get(OUTSIDE_REGION) : 0;
				ratioList[key][0] = numerator;
				ratioList[key][1] = denominator;
			}
			thisFileOutput.add(new OutputFileFormat(filesList[i], ratioList));
		}
		return thisFileOutput;
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
	private static Map<Integer, Map<String, Integer>> countOverlaps(Segment[][] segs, ArrayList<PValueFileFormat> pValueRecords) {
		boolean insideRegion;
		int pValueBin;
		Map<String, Integer> value;
		Map<Integer, Map<String, Integer>> overlapStats = new HashMap<Integer, Map<String, Integer>>();

		for (PValueFileFormat curRecord : pValueRecords) {
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
				overlapStats.put(pValueBin, value = new HashMap<String, Integer>());
			}
			// if there are segments for this chr
			if (segs[curRecord.chr] != null) {
				// find the overlap using binary search
				insideRegion = Segment.binarySearchForOverlap(curRecord.seg, segs[curRecord.chr]) != -1;
				// if it is inside the region
				if (insideRegion) {
					// if there is already a count then increase it else initialize with 1
					if (value.containsKey(INSIDE_REGION)) {
						value.put(INSIDE_REGION, (value.get(INSIDE_REGION)) + 1);
					} else {
						value.put(INSIDE_REGION, 1);
					}
				} else {
					if (value.containsKey(OUTSIDE_REGION)) {
						value.put(OUTSIDE_REGION, (value.get(OUTSIDE_REGION)) + 1);
					} else {
						value.put(OUTSIDE_REGION, 1);
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
				curLineParams = curLine.trim().split("[\\s]+");
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
		Segment[][] segs = null;
		String serializedSegFile = Files.getSerializedFilepath(segmentFile);
		LOGGER.info("The serialize  filepath: " + serializedSegFile);
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

	private static class PValueFileFormat implements Serializable {
		private static final long serialVersionUID = 1L;

		byte chr;
		Segment seg;
		double pValue;

		private PValueFileFormat(byte chr, Segment seg, double pValue) {
			this.chr = chr;
			this.seg = seg;
			this.pValue = pValue;
		}
	}
} // end of DnaseEnrichment class
