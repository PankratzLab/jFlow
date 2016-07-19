package cyto;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Hashtable;

import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.manage.SourceFileParser;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;

/**
 * Parser for cytogenic files, can be called on a directory, or can be passed a String[] of files to parse. Output files have new extension ".genvisis", and the project file is updated.
 * 
 */
public class CytoAgilentParse {
	public static final String[] SYSTEMATIC_NAME = { "SystematicName" };
	public static final String[] MARKER_POSITION_HEADER = { "Marker", "Chr", "Position" };
	public static final String[] DATA_TO_GRAB = { "ProbeName", "gProcessedSignal", "rProcessedSignal", "LogRatio", "gMedianSignal", "rMedianSignal" };
	public static final String[] CONVERT_TO = { "SNP Name", Sample.DATA_FIELDS[7][0], Sample.DATA_FIELDS[8][0], Sample.DATA_FIELDS[3][0], Sample.DATA_FIELDS[4][0], Sample.GENOTYPE_FIELDS[0][0], Sample.GENOTYPE_FIELDS[1][0], Sample.GENOTYPE_FIELDS[2][0], Sample.GENOTYPE_FIELDS[3][0] };

	private static final String[] SPLITS = { "\t" };
	private static final String CHR = "chr";
	private static final String GENVISIS_EXT = ".genvisis";
	private static final float SCALE_FACTOR = 2000f;

	/**
	 * Will parse all appropriate files in the projects source directory
	 * 
	 * @param proj
	 * @param log
	 * @return the list of parsed files
	 */
	public static String[] parseCytoToGenvisis(Project proj, Logger log) {
		String[] files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		if (files.length < 1) {
			log.reportError("Error - did not find any files to parse, please make sure the filename extension is set in the project properties file " + proj.getPropertyFilename());
			return new String[0];
		} else {
			for (int i = 0; i < files.length; i++) {
				files[i] = proj.SOURCE_DIRECTORY.getValue(false, true) + files[i];
			}
		}
		return parseCytoToGenvisis(proj, files, log);
	}

	/**
	 * Main function to generate .genvisis files that will pass through ParseIllumina.
	 * <p>
	 * The first sample is used to generate a markerPositions file. Returned are the output files, which are just the input files with a ".genvisis" extension.
	 * <p>
	 * Project properties are modified (SOURCE_FILENAME_EXTENSION -> GENVISIS_EXT, ID_HEADER -> ParseIllumina.FILENAME_AS_ID_OPTION, PARSE_AT_AT_SYMBOL ->FALSE ,SOURCE_FILE_DELIMITER -> TAB)
	 * <p>
	 * Went with intermediate files to avoid a new parser, there are also some quirks
	 * <p>
	 * Warning - files will be overwritten so appropriate checks should be made
	 * 
	 * @param proj
	 * @param filesToParse
	 *            (String[] of full paths)
	 * @param log
	 * @return
	 */
	public static String[] parseCytoToGenvisis(Project proj, String[] filesToParse, Logger log) {
		boolean createdMarkerPostions = false;
		String[] parsedFiles = new String[filesToParse.length];
		for (int i = 0; i < filesToParse.length; i++) {
			String genOutput = proj.SOURCE_DIRECTORY.getValue(true, true) + ext.rootOf(filesToParse[i]) + GENVISIS_EXT;
			parsedFiles[i] = genOutput;
			// TODO change to not , currently this overwrites existing files?
			// if (!Files.exists(genOutput)) {
			if (parseCytoToGenvisis(proj, filesToParse[i], genOutput, SCALE_FACTOR, log)) {
//				if ((!Files.exists(proj.getFilename(proj.MARKERSET_FILENAME, null, false, false)) && !Files.exists(proj.getFilename(proj.MARKER_POSITION_FILENAME, null, false, false))) || !createdMarkerPostions) {
				if ((!Files.exists(proj.MARKERSET_FILENAME.getValue(false, false)) && !Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false))) || !createdMarkerPostions) {
//					generateMarkerPositions(filesToParse[i], proj.getFilename(proj.MARKER_POSITION_FILENAME), log);
					generateMarkerPositions(filesToParse[i], proj.MARKER_POSITION_FILENAME.getValue(true, false), log);
					createdMarkerPostions = true;
				}
			} else {
				log.report("Error - could not parse file " + filesToParse[i] + ", if this is not an actual sample file, then this is not actually an error and will simply be skipped");
			}
			// } else {
			// log.report("Warning - skipping parsing of " + filesToParse[i] + ", the output file " + genOutput + " already exists");
			// }
		}

		proj.setProperty(proj.SOURCE_FILENAME_EXTENSION, GENVISIS_EXT);
		proj.setProperty(proj.ID_HEADER, SourceFileParser.FILENAME_AS_ID_OPTION);
//		proj.setProperty(proj.PARSE_AT_AT_SYMBOL, "FALSE");
		proj.setProperty(proj.PARSE_AT_AT_SYMBOL, Boolean.FALSE);
		proj.setProperty(proj.SOURCE_FILE_DELIMITER, Project.SOURCE_FILE_DELIMITERS.TAB);//"TAB");
		proj.saveProperties();
		return parsedFiles;
	}

	/**
	 * Generates an intermediate file with a .genvisis extension to be parsed with parseIllumina
	 * 
	 * @param proj
	 *            current project
	 * @param fileToParse
	 *            the file (full path)
	 * @param output
	 *            yup
	 * @param scale
	 *            the scale factor to scale gProcessedSignal and rProcessedSignal ( which are used as X,Y)
	 * @param log
	 * @return
	 */
	public static boolean parseCytoToGenvisis(Project proj, String fileToParse, String output, float scale, Logger log) {
		boolean valid = true;
		try {
			BufferedReader reader = Files.getAppropriateReader(fileToParse);
			String[] line;
			int[] indices;
			int count = 0;
			Hashtable<String, String> track = new Hashtable<String, String>();
			do {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				indices = ext.indexFactors(DATA_TO_GRAB, line, true, log, false, false);
			} while (reader.ready() && (indices[0] == -1));

			for (int i = 0; i < indices.length; i++) {
				if (indices[i] == -1) {
					log.reportError("Error - could not find neccesary column " + DATA_TO_GRAB[i] + " in file " + fileToParse);
					valid = false;
				}
			}
			if (!reader.ready()) {
				log.reportError("Error - could not neccesary columns with headers " + Array.toStr(DATA_TO_GRAB));
				valid = false;

			} else if (valid) {
				PrintWriter writer = Files.getAppropriateWriter(output);
				writer.println(Array.toStr(CONVERT_TO));
				while (reader.ready()) {
					line = reader.readLine().trim().split(SPLITS[0], -1);
					if (!track.containsKey(line[indices[0]])) {
						String[] converted = convert(line, indices, scale, log);
						writer.println(Array.toStr(converted));
						track.put(line[indices[0]], line[indices[0]]);
						count++;
					}
				}
				writer.close();
			}
			reader.close();
			log.report(ext.getTime() + " Info - found " + count + " unique markers");
		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find file " + fileToParse);
			log.reportException(e);
		} catch (IOException e) {
			log.reportError("Error - could not read file " + fileToParse);
			log.reportException(e);
		}
		return valid;
	}

	/**
	 * the genotypes being set to "A" are place holders for the 180K + SNP chips to come
	 */
	private static String[] convert(String[] input, int[] indices, float scale, Logger log) {
		String[] conversion = new String[CONVERT_TO.length];
		// Snp name
		conversion[0] = input[indices[0]];
		// BAF custom =
		conversion[1] = getBAF(input[indices[4]], input[indices[5]], input[indices[3]]) + "";
		// LRR (log10(r/g)
		conversion[2] = input[indices[3]];
		// X is the gProcessedSignal
		conversion[3] = scaleXY(input[indices[1]], scale) + "";
		// Y is the rProcessedSignal
		conversion[4] = scaleXY(input[indices[2]], scale) + "";
		conversion[5] = "A";
		conversion[6] = "A";
		conversion[7] = "A";
		conversion[8] = "A";
		return conversion;
	}

	private static float getBAF(String gMedianSignal, String rMedianSignal, String logRRatio) {
		return getBAF(Float.parseFloat(gMedianSignal), Float.parseFloat(rMedianSignal), Float.parseFloat(logRRatio));
	}

	/**
	 * A cute litte BAF fake thing. First we compute the log Ratio using gMedianSignal and rMedianSignal (as opposed to gProcessedSignal and rProcessedSignal). If this value is less than 0, we set it to 0. Else, we limit it to a max of 1(as BAF must be between 0 and 1)
	 * <p>
	 * This definitely truncates the values, but you can still get a sense of the aberration i.e, this value should follow the actual Log R Ratio
	 * 
	 */
	private static float getBAF(float gMedianSignal, float rMedianSignal, float logRRatio) {
		float medianLogRRatio = (float) Math.log10(rMedianSignal / gMedianSignal);
		if (medianLogRRatio < 0) {
			return 0;

		}
		return Math.min(Math.abs(medianLogRRatio), 1);
	}

	/**
	 * Since the gProcessedSignal and rProcessedSignal are often too large
	 */
	private static float scaleXY(String val, float scale) {
		return scaleXY(Float.parseFloat(val), scale);
	}

	/**
	 * gProcessedSignal and rProcessedSignal are set to X and Y. Since they are often quite large, we scale them.
	 */
	private static float scaleXY(float val, float scale) {
		return val / scale;
	}

	/**
	 * Generate a marker positions file from a cyto input file. The Systematic Name column is used to parse positions from the UCSC format
	 * <p>
	 * If the entry for a probeName in the Systematic Name file does not contain a valid UCSC location, the position is set to chr =0, position=0
	 * <p>
	 * For valid UCSC entries, the start of the the region is taken to be the position
	 */
	public static void generateMarkerPositions(String input, String output, Logger log) {
		log.report(ext.getTime() + " Info - generating a marker position file from " + input);
		log.report(ext.getTime() + " Info - positions will be parsed from UCSC locations only, using only the start position of the probe");
		log.report(ext.getTime() + " Info - probes without a valid UCSC location will be set to chr=0, position=0");

		try {
			BufferedReader reader = Files.getAppropriateReader(input);
			String[] line;
			int systematicIndex = -1;
			int probeNameIndex = -1;
			int count = 0;
			do {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				systematicIndex = ext.indexOfStr(SYSTEMATIC_NAME[0], line, true, true);
				probeNameIndex = ext.indexOfStr(DATA_TO_GRAB[0], line, true, true);
			} while (reader.ready() && (systematicIndex == -1 || probeNameIndex == -1));

			if (!reader.ready()) {
				log.reportError("Error - could not neccesary columns with headers " + Array.toStr(SYSTEMATIC_NAME) + " and " + DATA_TO_GRAB[0]);
				return;
			}
			PrintWriter writer = Files.getAppropriateWriter(output);
			writer.println(Array.toStr(MARKER_POSITION_HEADER));
			while (reader.ready()) {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				String potentialUCSC = line[systematicIndex];
				String probeName = line[probeNameIndex];
				int[] loc = getLoc(potentialUCSC, log);
				writer.println(probeName + SPLITS[0] + loc[0] + SPLITS[0] + loc[1]);
				count++;
			}
			writer.close();
			reader.close();
			log.report(ext.getTime() + " Info - generated positions file for " + count + " markers");
		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find file " + input);
			log.reportException(e);
		} catch (IOException e) {
			log.reportError("Error - could not read file " + input);
			log.reportException(e);
		}
	}

	/**
	 * @param potentialUCSC
	 *            if the potential UCSC position starts with "chr", we take it. Other wise we set the positon to chr=0,position =0;
	 * @param log
	 * @return
	 */
	private static int[] getLoc(String potentialUCSC, Logger log) {
		int[] loc = new int[3];
		if (potentialUCSC.startsWith(CHR)) {
			loc = Positions.parseUCSClocation(potentialUCSC);
		} else {
			Arrays.fill(loc, 0);
		}
		return loc;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/workspace/Genvisis/projects/HirchCyto.properties";
		String logFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/testlog.log";

		String usage = "\n" + "jlDev.CytoCNVariant requires 0-2 arguments\n";
		usage += "   (1) name of the project filename (i.e. filename=" + filename + " (default))\n";
		usage += "   (2) name of the log file (i.e. logFile=" + logFile + " (default))\n";
		usage += "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("filename=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("logFile=")) {
				logFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Invalid argument " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Logger log = new Logger(logFile);
		Project proj = new Project(filename, false);
		parseCytoToGenvisis(proj, log);
	}
}
