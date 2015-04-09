package cyto;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import cnv.analysis.BeastScore;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;

/**
 * Class to compare CytoCNVariants to other segments (common CNPS, un-reported, reported), currently only supports three files
 * <p>
 * It is not necessary to have a valid project file, however, a project file is necessary to compute beast scores
 * 
 */
public class CytoCompare {
	private static final String[] OUTPUT_HEADER = { "REGION", "UCSC_LINK", "ISCN", "INTERPRETATION", "GENE(s)", "AVERAGE_LOG_RATIO", "NUMBER_OF_PROBES", "CNP_OVERLAP", "REPORTED_OVERLAP", "UNREPORTED_OVERLAP", "NUM_UNREPORTED", "BEAST_SCORE" };
	private static final String[] TRACK_HEADERS = { "track name=" };
	private static final String[] CHR = { "chr" };
	private static final int[] INDICES_TO_LOAD = { 0, 1, 2 };
	public static final String[] SPLITS = { "\t" };
	public static final String EXT = ".xls";
	public static final String DEFAULT_LOG = "workbench.log";
	public static final double DEFAULT_SCORE_THRESHOLD = 0.5;

	/**
	 * @param proj
	 *            not necessary, except for Beast scores
	 * @param dir
	 *            where the cytoCNVariantFiles are located
	 * @param cytoCNVariantFiles
	 *            String[] file names (not full paths)
	 * @param cnpFile
	 *            full path
	 * @param reportedOverlapFile
	 *            full path
	 * @param unreportedOverlapFile
	 *            full path
	 * @param scoreThreshold
	 *            score threshold for the overlap scoring {@link Segment#overlapScore}
	 * @param outputDir
	 *            output directory
	 * @param computeBeast
	 *            compute the Beast score (requires project and sample data)
	 * @param log
	 * 
	 * @return String[][] containing the full paths to the output files after parsing cytoCNVariantFiles
	 */
	public static String[][] compare(Project proj, String dir, String[] cytoCNVariantFiles, String cnpFile, String reportedOverlapFile, String unreportedOverlapFile, double scoreThreshold, String outputDir, boolean computeBeast, Logger log) {
		return compare(proj, toFullPaths(dir, cytoCNVariantFiles), loadsegs(cnpFile, log), loadsegs(reportedOverlapFile, log), loadsegs(unreportedOverlapFile, log), scoreThreshold, outputDir, computeBeast, log);
	}

	/**
	 * Same as above, but full paths
	 */
	public static String[][] compare(Project proj, String[] cytoCNVariantFiles, String cnpFile, String reportedOverlapFile, String unreportedOverlapFile, double scoreThreshold, String outputDir, boolean computeBeast, Logger log) {
		return compare(proj, cytoCNVariantFiles, loadsegs(cnpFile, log), loadsegs(reportedOverlapFile, log), loadsegs(unreportedOverlapFile, log), scoreThreshold, outputDir, computeBeast, log);
	}

	/**
	 * Same as above, but full paths, and files loaded as segments
	 */
	public static String[][] compare(Project proj, String[] cytoCNVariantFiles, Segment[] cnpSegs, Segment[] reportSegs, Segment[] unReportSegs, double scoreThreshold, String outputDir, boolean computeBeast, Logger log) {
		String[][] alloutputs = new String[cytoCNVariantFiles.length][];
		for (int i = 0; i < cytoCNVariantFiles.length; i++) {
			if (!Files.exists(cytoCNVariantFiles[i])) {
				log.reportError("Error - could not find the file " + cytoCNVariantFiles[i]);
				return alloutputs;
			} else {
				CytoCNVariant[] cytoCNVariant = CytoCNVariant.loadCytoCNVariant(cytoCNVariantFiles[i], log);
				alloutputs[i] = compare(proj, cytoCNVariant, cnpSegs, reportSegs, unReportSegs, scoreThreshold, outputDir, computeBeast, log);
			}
		}
		return alloutputs;
	}

	/**
	 * For each individual, we summarize the overlaps (if available) to separate files. Currently file names are determined from the FID/IID of the CytoCNVariants
	 * <p>
	 * If a project and sample data is available (with parsed samples) and computeBeast==true, we compute the beast score.
	 */
	public static String[] compare(Project proj, CytoCNVariant[] cytoCNVariants, Segment[] cnpSegs, Segment[] reportSegs, Segment[] unReportSegs, double scoreThreshold, String outputDir, boolean computeBeast, Logger log) {
		CytoCNVariant[][] cytoCNVariantInds = CytoCNVariant.toIndividuals(cytoCNVariants, log);
		String[] outputs = new String[cytoCNVariantInds.length];
		BeastScore[] beastScores = beast(proj, cytoCNVariantInds, computeBeast, log);
		for (int i = 0; i < cytoCNVariantInds.length; i++) {
			if (cytoCNVariantInds[i].length > 0) {
				String ind = cytoCNVariantInds[i][0].getFamilyID() + "_" + cytoCNVariantInds[i][0].getIndividualID();
				String output = ind + EXT;
				outputs[i] = outputDir + output;
				checkFile(outputDir, output, log);

				PrintWriter writer = Files.getAppropriateWriter(outputs[i]);
				writer.println(Array.toStr(OUTPUT_HEADER));

				for (int j = 0; j < cytoCNVariantInds[i].length; j++) {
					if (ind.equals(cytoCNVariantInds[i][j].getFamilyID() + "_" + cytoCNVariantInds[i][j].getIndividualID())) {

						String report = cytoCNVariantInds[i][j].getMyReport();
						report += getOverLapSummary(cytoCNVariantInds[i][j], cnpSegs, scoreThreshold, false, log);
						report += getOverLapSummary(cytoCNVariantInds[i][j], reportSegs, scoreThreshold, false, log);
						report += getOverLapSummary(cytoCNVariantInds[i][j], unReportSegs, scoreThreshold, true, log);

						if (beastScores == null || beastScores[i] == null) {
							report += SPLITS[0] + ext.MISSING_VALUES[2];
						} else {
							report += SPLITS[0] + beastScores[i].getBeastScores()[j];
						}

						writer.println(report);
					} else {
						log.reportError("Error - the array of cyto aberrations contained mismatched IDs, this should not happen");
					}
				}
				writer.close();
			}
		}
		return outputs;
	}

	/**
	 * Computes the beast score for an aberration (if computeBeast is flagged, and a valid project/sample data file exists
	 * 
	 * @param proj
	 * @param cNVariantInds
	 *            CNVariant[][] organized as cNVariantInds[sample0][variantsForSample0]
	 * @param computeBeast
	 *            flag to even try
	 * @param log
	 * @return null if failed or don't want to do it, else the beast scores
	 */
	private static BeastScore[] beast(Project proj, CNVariant[][] cNVariantInds, boolean computeBeast, Logger log) {
		BeastScore[] beastScores = null;
		if (computeBeast) {
			if (proj == null) {
				log.reportError("Warning - a valid project file is required to compute beast scores, skipping beast scores :(");
			} else if (!Files.exists(proj.getFilename(proj.MARKERLOOKUP_FILENAME))) {
				log.reportError("Warning - a marker lookup file is needed to compute beast scores, skipping beast scores :(");
			} else if (!Files.exists(proj.getFilename(proj.SAMPLE_DATA_FILENAME))) {
				log.reportError("Warning - a sample data file is needed to lookup samples by FID/IID, skipping beast scores :(");
			} else {
				log.report(ext.getTime() + " Info - computing beast scores for " + cNVariantInds.length + " individuals");
				beastScores = BeastScore.beastInds(proj, cNVariantInds);
			}
		}
		return beastScores;
	}

	/**
	 * If we have valid segments to compare to, we compute the overlap
	 * 
	 * @param cytoCNVariantIndsSegment
	 * @param segs
	 * @param scoreThreshold
	 * @param reportNumOverlapThreshold
	 * @param log
	 * @return
	 */
	private static String getOverLapSummary(CytoCNVariant cytoCNVariantIndsSegment, Segment[] segs, double scoreThreshold, boolean reportNumOverlapThreshold, Logger log) {
		String summary = "";
		Segment.SegmentCompare segmentCompare;
		if (segs != null) {
			segmentCompare = compareIndsToTarget(cytoCNVariantIndsSegment, segs, scoreThreshold, log);
			summary = SPLITS[0] + segmentCompare.getMaximumOverlapScore() + (reportNumOverlapThreshold ? SPLITS[0] + segmentCompare.getNumOverlapingPastThreshold() : "");
		} else {
			summary = SPLITS[0] + ext.MISSING_VALUES[2] + (reportNumOverlapThreshold ? SPLITS[0] + ext.MISSING_VALUES[2] : "");
		}
		return summary;
	}

	/**
	 * @param individualCytoCNVariant
	 *            the variant to compare
	 * @param compareSegs
	 *            the segments to compare to (find the max score, maxscore segment, and num passing scoreThreshold
	 * @param scoreThreshold
	 * @param log
	 * @return
	 */
	private static Segment.SegmentCompare compareIndsToTarget(CytoCNVariant cytoCNVariantIndsSegment, Segment[] compareSegs, double scoreThreshold, Logger log) {
		Segment.SegmentCompare cytoAgilentCompare = cytoCNVariantIndsSegment.new SegmentCompare(compareSegs, scoreThreshold, log);
		cytoAgilentCompare.compare();
		return cytoAgilentCompare;
	}

	/**
	 * We allow null input and output
	 * 
	 * @param segFile
	 * @param log
	 * @return
	 */
	public static Segment[] loadsegs(String segFile, Logger log) {
		Segment[] segs = null;
		if (segFile != null && !Files.exists(segFile)) {
			log.reportError("Warning - could not find the file " + segFile + " skipping this comparison");
		} else if (segFile != null) {
			segs = loadCytoSegments(segFile, log);
		}
		return segs;
	}

	/**
	 * The input segment files have some goofy spacing going on (blank lines, "chr9 \t", etc), to this end, we use this method instead of Segment.loadRegions... Need to check line length and parse a bit more carefully
	 */
	private static Segment[] loadCytoSegments(String filename, Logger log) {
		ArrayList<Segment> segs = new ArrayList<Segment>();
		if (!Files.exists(filename)) {
			log.reportError("Error - could not load file " + filename + " it does not exist at the current path");
		} else {
			boolean skipFirstLine = true;
			String[] header = Files.getHeaderOfFile(filename, log);
			if (ext.indexOfStr(TRACK_HEADERS[0], header, false, true) == -1 && ext.indexOfStr(CHR[0], header) >= 0) {
				skipFirstLine = false;
			}
			try {
				BufferedReader reader = Files.getAppropriateReader(filename);
				if (skipFirstLine) {
					reader.readLine();
				}
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split(SPLITS[0]);
					if (line.length >= 3) {
						try {
							byte chr = Positions.chromosomeNumber(line[INDICES_TO_LOAD[0]].replaceAll(" ", ""), log);
							int start = Integer.parseInt(line[INDICES_TO_LOAD[1]]);
							int stop = Integer.parseInt(line[INDICES_TO_LOAD[2]]);
							segs.add(new Segment(chr, start, stop));
						} catch (NumberFormatException nfe) {
							log.reportError("Error - could not parse the segment on line " + Array.toStr(line));
						}
					}
				}
			} catch (FileNotFoundException e) {
				log.reportError("Error - could not find file " + filename);
				log.reportException(e);
			} catch (IOException e) {
				log.reportError("Error - reading file" + filename);
				log.reportException(e);
			}
		}
		log.report(ext.getTime() + " Info - found " + segs.size() + (segs.size() > 1 ? " segments " : " segment ") + "in " + filename);
		return segs.toArray(new Segment[segs.size()]);
	}

	private static void checkFile(String dir, String file, Logger log) {
		if (Files.exists(dir + file)) {
			log.reportError("Warning - the file " + dir + file + " already exists, backing it up");
			Files.backup(file, dir, dir);
		}
	}

	/**
	 * @param dir
	 *            add to each fo the files
	 * @param files
	 *            the files to add
	 * @return new array of full paths
	 */

	private static String[] toFullPaths(String dir, String[] files) {
		String[] fullPaths = new String[files.length];
		for (int i = 0; i < files.length; i++) {
			fullPaths[i] = dir + files[i];
		}
		return fullPaths;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "D:/data/Hirch_CYTO/Hirch_Cytolab/";
		String cytoCNVariantFile = "5 6 7 8 workbench.xls";
		String cnpFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/HG19 CNV edit for AGW.txt";
		String reportedOverlapFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/HG19 Reported 2012.05.22.txt";
		String unreportedOverlapFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/HG19 Unreported 2012.05.22-2.txt";
		String outputDir = "D:/data/Hirch_CYTO/Hirch_Cytolab/";
		String logFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/Cytolog.log";
		double scoreThreshold = DEFAULT_SCORE_THRESHOLD;
		String filename = null;
		filename = "C:/workspace/Genvisis/projects/HirchCyto.properties";
		boolean computeBeast = true;
		Project proj = null;

		String usage = "\n" + "jlDev.ParseAgilent requires 0-7 arguments\n";
		usage += "   (1) cyto call file (i.e. filename=" + cytoCNVariantFile + " (default))\n";
		usage += "   (2) a file of CNP locations to compute overlap (i.e. cnpFile=" + cnpFile + " (default))\n";
		usage += "   (3) a file of reported aberrations to compute overlap (i.e. reportedOverlapFile=" + reportedOverlapFile + " (default))\n";
		usage += "   (4) a file of un-reported aberrations to compute overlap (i.e. unreportedOverlapFile=" + unreportedOverlapFile + " (default))\n";
		usage += "   (5) the threshold to start counting aberration overlap (i.e. scoreThreshold=" + scoreThreshold + " (default))\n";
		usage += "   (6) the output directory, output will be named according to the sample names in the cyto call file (i.e. outputDir=" + outputDir + " (default))\n";
		usage += "   OPTIONAL: ";
		usage += "   (8) compute beast scores (a valid sample data file must exist, samples must be parsed, and a project file must be provided) in addition to comparing segments (i.e.-beast (not the default))\n";
		usage += "   (7) filename of the project (Only neccesary if computing Beast Scores (i.e. filename=" + filename + " (no default))\n";
		usage += "   (8) name of the log file (i.e. logFile=" + logFile + " (default))\n";
		usage += "   (9) directory where cyto call file is located (i.e. dir=" + dir + " (default))\n";

		usage += "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filename=")) {
				cytoCNVariantFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnpFile=")) {
				cnpFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("reportedOverlapFile=")) {
				reportedOverlapFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("unreportedOverlapFile=")) {
				unreportedOverlapFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("scoreThreshold=")) {
				scoreThreshold = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("outputDir=")) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filename=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-beast")) {
				computeBeast = true;
				numArgs--;
			} else if (args[i].startsWith("logFile=")) {
				logFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (filename != null) {
			proj = new Project(filename, false);
		}
		Logger log = new Logger(logFile);
		compare(proj, dir, new String[] { cytoCNVariantFile }, cnpFile, reportedOverlapFile, unreportedOverlapFile, scoreThreshold, outputDir, computeBeast, log);
	}
}
