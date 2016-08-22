package org.genvisis.seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class BWA_Analysis {
	public static class BWA_AnalysisIndividual {
		private static final String OUTPUT_SEP = "_";

		public static String[] getBatchesByLane(BWA_AnalysisIndividual[] bwIndividuals) {
			Hashtable<String, ArrayList<BWA_AnalysisIndividual>> track =
			                                                           new Hashtable<String, ArrayList<BWA_AnalysisIndividual>>();
			ArrayList<String> unique = new ArrayList<String>();
			for (BWA_AnalysisIndividual bwIndividual : bwIndividuals) {
				String baseId = bwIndividual.getID();
				if (!track.containsKey(baseId)) {
					track.put(baseId, new ArrayList<BWA_AnalysisIndividual>());
					unique.add(baseId);
				}
				track.get(baseId).add(bwIndividual);
			}
			String[] batchedByLane = new String[unique.size()];
			for (int i = 0; i < unique.size(); i++) {
				ArrayList<BWA_AnalysisIndividual> current = track.get(unique.get(i));
				batchedByLane[i] = "";
				for (int j = 0; j < current.size(); j++) {
					if (j == 0) {
						batchedByLane[i] = current.get(j).getAvailableFiles("\t");
					} else {
						batchedByLane[i] += "\n" + current.get(j).getAvailableFiles("\t");
					}
				}
			}
			return batchedByLane;
		}

		private final BWA bwa;
		private final String ID;
		private final String Lane;
		private String library;
		private final String barcode;
		private final String referenceGenomeFasta;
		private String readFQ1;
		private String readFQ2;
		private String output;
		private final String outputDir;
		private boolean fail;
		private boolean success;

		private final Logger log;

		public BWA_AnalysisIndividual(BWA bwa, String outputDir, String iD, String lane, String library,
		                              String barcode, String referenceGenomeFasta, Logger log) {
			super();
			this.bwa = bwa;
			this.outputDir = outputDir;
			ID = iD;
			Lane = lane;
			this.library = library;
			this.barcode = barcode;
			this.referenceGenomeFasta = referenceGenomeFasta;
			fail = false;
			this.log = log;
			output = outputDir + ID + OUTPUT_SEP + barcode + OUTPUT_SEP + Lane + OUTPUT_SEP + library
			         + SAM_EXT;
			success = false;
		}

		public boolean analyze(int numMemThreads, Logger altLog) {
			if (!fail) {
				fail = !hasBothFiles();
				if (!fail) {
					new File(outputDir).mkdirs();
					success = bwa.bwaMEM(referenceGenomeFasta, readFQ1, readFQ2, output, getReadGroup(),
					                     numMemThreads, (altLog == null ? log : altLog));
				} else {
					log.reportError("could not find both files for ID" + ID);
				}
			}
			return success;
		}

		public void assignFile(String readFQ) {
			if (readFQ1 == null) {
				readFQ1 = readFQ;
			} else if (readFQ2 == null) {
				readFQ2 = readFQ;
			} else {
				log.reportError("internal error, adding too many samples to the bwa analysis");
				fail = true;
			}
		}

		public String getAvailableFiles(String sep) {
			String available = "";
			if (hasBothFiles()) {
				available += readFQ1;
				available += sep + readFQ2;
			} else {
				if (readFQ1 != null) {
					available += sep + readFQ1;
				}
				if (readFQ2 != null) {
					available += sep + readFQ2;
				}
			}

			return available;
		}

		public String getID() {
			return ID;
		}

		public String getLane() {
			return Lane;
		}

		public String getLibrary() {
			return library;
		}

		public String getOutput() {
			return output;
		}

		// @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1
		public String getReadGroup() {
			String RG = "\"@RG";
			RG += "\\tID:" + ID + "_" + barcode + "_" + Lane + "_" + library;// unique for this run of the
			                                                                 // sample
			RG += "\\tSM:" + ID;// unique for sample
			RG += "\\tPL:ILLUMINA";// not used, currently TODO
			RG += "\\tLB:" + library;// library prep, I think this is correct TODO
			RG += "\"";
			return RG;
		}

		public boolean hasBothFiles() {
			return readFQ1 != null && readFQ2 != null;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public void setLibrary(String library) {
			this.library = library;
		}

		public void setOutput(String output) {
			this.output = output;
		}

	}
	public static class FileNameParser {
		public static final String SPLIT = "_";
		private final String fileName;
		private final String[] split;
		private String ID;
		private String lane;
		private String batch;
		private String barcode;
		private boolean valid;
		private final Logger log;

		public FileNameParser(String fileName, Logger log) {
			super();
			this.fileName = fileName;
			split = ext.rootOf(fileName.trim()).split(SPLIT);
			valid = true;
			this.log = log;
		}

		public String getBarcode() {
			return barcode;
		}

		public String getBatch() {
			return batch;
		}

		public String getID() {
			return ID;
		}

		public String getIDLane() {
			return ID + SPLIT + lane;
		}

		public String getLane() {
			return lane;
		}

		public boolean isValid() {
			return valid;
		}

		public void parse() {
			if (split.length < 4) {
				log.reportError("could not parse filename " + fileName
				                + " to ID, lane, barcode, and batch");
				valid = false;
			} else {
				lane = split[split.length - 3];
				batch = split[split.length - 1];
				barcode = split[split.length - 4];
				ID = Array.toStr(Array.subArray(split, 0, split.length - 4), SPLIT);// we do not include
				                                                                    // barcode in the id,
				                                                                    // instead adding it to
				                                                                    // the RG
			}
		}

	}
	private static class WorkerBWA_Analysis implements Callable<Boolean> {
		private final BWA_AnalysisIndividual bwAnalysisIndividual;
		private final boolean verbose;
		private final int numMemThreads;
		private final Logger log;

		public WorkerBWA_Analysis(BWA_AnalysisIndividual bwAnalysisIndividual, int numMemThreads,
		                          boolean verbose, Logger log) {
			super();
			this.bwAnalysisIndividual = bwAnalysisIndividual;
			this.verbose = verbose;
			this.numMemThreads = numMemThreads;
			this.log = log;
		}

		@Override
		public Boolean call() {// acts like run
			boolean success = false;
			if (!bwAnalysisIndividual.isFail()) {
				if (verbose) {
					log.report(ext.getTime() + "Info - running bwa mem on thread"
					           + Thread.currentThread().getName() + " for "
					           + bwAnalysisIndividual.getAvailableFiles("\n"));
				}
				success = bwAnalysisIndividual.analyze(numMemThreads, log);
				if (verbose) {
					log.report(ext.getTime() + "Info - finished running bwa mem on thread"
					           + Thread.currentThread().getName() + " for "
					           + bwAnalysisIndividual.getAvailableFiles("\n"));
				}
			} else {
				log.reportError("initializing has failed for analysis "
				                + bwAnalysisIndividual.getAvailableFiles("\n"));
			}
			return success;
		}
	}

	public static final String INDEXED_REFERENCE_GENOME_EXT = ".sa";
	public static final String SAM_EXT = ".sam";

	public static final String[] FQ_EXTS = {".fq", ".fastq"};
	public static final String PAIRED_SUMMARY = "fqLaneReadPaired.txt";
	public static final String SPACE = " ";
	public static final String ROOT_INPUT_COMMAND = "rootInputDir=";
	public static final String ROOT_OUTPUT_COMMAND = "rootOutputDir=";
	public static final String FILE_OF_SAMPLE_PAIRS_COMMAND = "fileOfSamplePairs=";
	public static final String REFERENCE_GENOME_COMMAND = "referenceGenomeFasta=";
	public static final String BWA_LOCATION_COMMAND = "bwaLocation=";
	public static final String QUIET_COMMAND = "-quiet";
	public static final String NUM_BETWEEN_THREADS_COMMAND = "numBetweenThreads=";
	public static final String NUM_WITHIN_THREADS_COMMAND = "numWithinThreads=";
	public static final String LOG_FILE_COMMAND = "logFile=";
	public static final String BATCH_COMMAND = "-batch";

	public static final String NUMBATCHES_COMMAND = "numBatches=";
	public static final String OVERWRITE_EXISTING_COMMAND = "-overwrite";
	public static final String BASE_NAME_COMMAND = "baseName=";

	private final String rootInputDir;
	private final String rootOutputDir;
	private final String referenceGenomeFasta;
	private String[] inputFiles;

	private BWA_AnalysisIndividual[] bwAnalysisIndividuals;

	private final BWA bwa;

	private boolean fail;

	private final boolean verbose;

	private final int numWithinSampleThreads;

	private final int numBetweenSampleThreads;

	private final Logger log;

	public BWA_Analysis(String rootInputDir, String rootOutputDir, String referenceGenomeFasta,
	                    boolean verbose, int numWithinSampleThreads, int numBetweenSampleThreads,
	                    BWA bwa, Logger log) {
		super();
		this.rootInputDir = rootInputDir;
		this.rootOutputDir = rootOutputDir;
		new File(rootOutputDir).mkdirs();
		this.referenceGenomeFasta = referenceGenomeFasta;
		this.verbose = verbose;
		this.numWithinSampleThreads = numWithinSampleThreads;
		this.numBetweenSampleThreads = numBetweenSampleThreads;
		this.bwa = bwa;
		fail = bwa.isFail();
		this.log = log;
	}

	public boolean analyzeBWA_MEM() {
		boolean success = false;
		if (!fail) {
			fail = !verifyReferenceGenome();
			verifyAnalsyisInds();
			if (!fail) {
				ExecutorService executor = Executors.newFixedThreadPool(numBetweenSampleThreads);
				Hashtable<String, Future<Boolean>> tmpResults = new Hashtable<String, Future<Boolean>>();
				for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
					Logger tmpLog = new Logger(ext.rootOf(log.getFilename(), false) + "_BWA_ID_"
					                           + bwAnalysisIndividuals[i].getID() + "_batch"
					                           + bwAnalysisIndividuals[i].getLibrary() + ".log");
					tmpResults.put(i + "",
					               executor.submit(new WorkerBWA_Analysis(bwAnalysisIndividuals[i],
					                                                      numWithinSampleThreads, verbose,
					                                                      tmpLog)));
				}
				for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
					try {
						if (!tmpResults.get(i + "").get()) {
							log.reportError("failed bwa for "
							                + bwAnalysisIndividuals[i].getAvailableFiles("\n"));
							fail = true;
						}
					} catch (InterruptedException e) {
						log.reportError("Could not complete running bwa mem on internal index " + i);
						log.reportException(e);
						fail = true;
					} catch (ExecutionException e) {
						log.reportError("Could not complete running bwa mem on internal index " + i);
						log.reportException(e);
						fail = true;
					}
				}
				executor.shutdown();
				try {
					executor.awaitTermination(10, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					log.reportException(e);
				}
			}
		}
		return success;
	}

	public void batch(int numBatches, int memoryInMB, int wallTimeInHours, String baseName) {
		String[] allFilesMatched = new String[bwAnalysisIndividuals.length];
		for (int i = 0; i < allFilesMatched.length; i++) {
			allFilesMatched[i] = bwAnalysisIndividuals[i].getAvailableFiles("\t");
		}
		String[][] batchedMatchedFiles = Array.splitUpStringArray(allFilesMatched, numBatches, log);

		String[][] batches = new String[batchedMatchedFiles.length][1];
		for (int i = 0; i < batches.length; i++) {
			batches[i][0] = "batch_" + i + "_" + baseName;
			Files.writeArray(batchedMatchedFiles[i], rootOutputDir + batches[i][0] + ".txt");
		}
		// String rootInputDir, String rootOutputDir, String referenceGenomeFasta, String bwaLocation,
		// String fileOfSamplePairs, boolean overwriteExisting, boolean verbose, int numMemThreads, int
		// numSampleThreads, boolean batch, Logger log) {

		String command = "load module java\njava -Xmx" + memoryInMB + "m -jar "
		                 + org.genvisis.common.PSF.Java.GENVISIS + " seq.BWA_Analysis "
		                 + ROOT_INPUT_COMMAND + rootInputDir + SPACE + ROOT_OUTPUT_COMMAND
		                 + rootOutputDir + SPACE;
		command += REFERENCE_GENOME_COMMAND + referenceGenomeFasta + SPACE + BWA_LOCATION_COMMAND
		           + bwa.getBwaLocation() + SPACE;
		command += NUM_BETWEEN_THREADS_COMMAND + numWithinSampleThreads + SPACE
		           + FILE_OF_SAMPLE_PAIRS_COMMAND + rootOutputDir + "[%0].txt" + SPACE
		           + NUM_WITHIN_THREADS_COMMAND + numBetweenSampleThreads;
		Files.qsub("BWA_MEM" + baseName, command, batches, memoryInMB, wallTimeInHours,
		           numWithinSampleThreads * numBetweenSampleThreads);
	}

	public BWA getBwa() {
		return bwa;
	}

	public BWA_AnalysisIndividual[] getBwAnalysisIndividuals() {
		return bwAnalysisIndividuals;
	}

	public String[] getInputFiles() {
		return inputFiles;
	}

	public Logger getLog() {
		return log;
	}

	public int getNumBetweenSampleThreads() {
		return numBetweenSampleThreads;
	}

	public int getNumWithinSampleThreads() {
		return numWithinSampleThreads;
	}

	public String getReferenceGenomeFasta() {
		return referenceGenomeFasta;
	}

	public String getRootInputDir() {
		return rootInputDir;
	}

	public String getRootOutputDir() {
		return rootOutputDir;
	}

	public void init(String fileOfSamplePairs) {
		if (fileOfSamplePairs == null) {
			gatherInputFilesFromRootDirectory();
			parseAnalysisIndividuals();
		} else {
			gatherInputFilesFromFile(fileOfSamplePairs);
		}
	}

	public boolean isFail() {
		return fail;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setBwAnalysisIndividuals(BWA_AnalysisIndividual[] bwAnalysisIndividuals) {
		this.bwAnalysisIndividuals = bwAnalysisIndividuals;
	}

	public void setFail(boolean fail) {
		this.fail = fail;
	}

	private void gatherInputFilesFromFile(String fileOfSamplePairs) {

		try {
			int numSamples = Files.countLines(fileOfSamplePairs, 0);
			bwAnalysisIndividuals = new BWA_AnalysisIndividual[numSamples];// this a per lane way of doing
			                                                               // it

			BufferedReader reader = Files.getAppropriateReader(fileOfSamplePairs);
			int index = 0;
			Hashtable<String, Integer> track = new Hashtable<String, Integer>();
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("[\\s]+");
				if (line.length != 2) {
					log.reportError(fileOfSamplePairs
					                + " must be two columns only, with lane matched samples in each");
					fail = true;
				} else {
					FileNameParser fileNameParser1 = new FileNameParser(line[0], log);
					fileNameParser1.parse();
					FileNameParser fileNameParser2 = new FileNameParser(line[1], log);
					fileNameParser2.parse();
					if (fileNameParser1.isValid() && fileNameParser2.isValid()) {
						if (!fileNameParser1.getLane().equals(fileNameParser2.getLane())) {
							log.reportError("Warning - the determined lane for the two samples "
							                + Array.toStr(line)
							                + " did not match up, please make sure this is what you want to do");
						}
						if (!fileNameParser1.getID().equals(fileNameParser2.getID())) {
							log.reportError("the determined root ID for the two samples "
							                + Array.toStr(line) + " did not match up");
							fail = true;
						} else {
							bwAnalysisIndividuals[index] = getAnalysisIndFromFileParser(fileNameParser1);
							bwAnalysisIndividuals[index].assignFile(line[0]);
							bwAnalysisIndividuals[index].assignFile(line[1]);
							if (track.containsKey(bwAnalysisIndividuals[index].getOutput())) {
								int num = track.get(bwAnalysisIndividuals[index].getOutput());
								bwAnalysisIndividuals[index].setOutput(ext.addToRoot(bwAnalysisIndividuals[index].getOutput(),
								                                                     "rep" + num));
								bwAnalysisIndividuals[index].setLibrary(bwAnalysisIndividuals[index].getLibrary()
								                                        + "rep" + num);
								track.put(bwAnalysisIndividuals[index].getOutput(), (num + 1));
							} else {
								track.put(bwAnalysisIndividuals[index].getOutput(), 1);
							}
							System.out.println(bwAnalysisIndividuals[index].getReadGroup());

							index++;

						}
					} else {
						fail = true;
					}
				}
			}

		} catch (FileNotFoundException e) {
			log.reportError("could not find file " + fileOfSamplePairs);
			log.reportException(e);
			e.printStackTrace();
		} catch (IOException e) {
			log.reportError("could not read file " + fileOfSamplePairs);
			log.reportException(e);
		}
	}

	private void gatherInputFilesFromRootDirectory() {
		if (verbose) {
			log.report(ext.getTime() + " Info - gathering samples by lane from " + rootInputDir);
		}
		for (String element : FQ_EXTS) {
			String[] tmpFiles = Files.list(rootInputDir, element, false);
			if (tmpFiles != null && tmpFiles.length > 0) {
				if (verbose) {
					log.report("Info - found " + tmpFiles.length + " of type " + element);
				}
				inputFiles = Files.toFullPaths(tmpFiles, rootInputDir);
				break;
			}
		}
	}

	private BWA_AnalysisIndividual getAnalysisIndFromFileParser(FileNameParser fileNameParser) {
		return new BWA_AnalysisIndividual(bwa, rootOutputDir, fileNameParser.getID(),
		                                  fileNameParser.getLane(), fileNameParser.getBatch(),
		                                  fileNameParser.getBarcode(), referenceGenomeFasta, log);
	}

	private void parseAnalysisIndividuals() {
		if (!fail) {
			fail = !verifyInputFiles();
			if (!fail) {
				Hashtable<String, Integer> track = new Hashtable<String, Integer>();
				bwAnalysisIndividuals = new BWA_AnalysisIndividual[inputFiles.length / 2];
				int currIndex = 0;
				for (String inputFile : inputFiles) {
					FileNameParser fileNameParser = new FileNameParser(inputFile, log);
					fileNameParser.parse();
					if (fileNameParser.isValid()) {
						if (!track.containsKey(fileNameParser.getIDLane())) {
							track.put(fileNameParser.getIDLane(), currIndex);
							bwAnalysisIndividuals[currIndex] = getAnalysisIndFromFileParser(fileNameParser);
							bwAnalysisIndividuals[currIndex].assignFile(inputFile);
							currIndex++;
						} else {
							bwAnalysisIndividuals[track.get(fileNameParser.getIDLane())].assignFile(inputFile);
						}
					} else {
						log.reportError("file name " + inputFile
						                + " could not be parsed according to our assumptions");
					}
				}
				if (currIndex != bwAnalysisIndividuals.length) {
					log.reportError("could not match all files for input set\n"
					                + Array.toStr(inputFiles, "\n"));
					fail = true;
				}
			}
		}
	}

	private void verifyAnalsyisInds() {
		for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
			if (!bwAnalysisIndividuals[i].hasBothFiles()) {
				log.reportError("internal index " + i + " only had the following files:"
				                + bwAnalysisIndividuals[i].getAvailableFiles("\n"));
				fail = true;
			}
		}
	}

	private boolean verifyInputFiles() {
		if (inputFiles != null && inputFiles.length % 2 == 0) {
			return true;
		} else {
			log.reportError("found an odd number of input files, currently they must be paired reads");
			return false;
		}
	}

	private boolean verifyReferenceGenome() {
		boolean verified = true;
		if (!Files.exists(referenceGenomeFasta)) {
			log.reportError("could not find reference genome file" + referenceGenomeFasta);
			verified = false;
		}
		return verified;
	}

	public static void run(String rootInputDir, String rootOutputDir, String referenceGenomeFasta,
	                       String bwaLocation, String fileOfSamplePairs, boolean overwriteExisting,
	                       boolean verbose, int numMemThreads, int numSampleThreads, boolean batch,
	                       int numBatches, int memoryInMB, int wallTimeInHours, String baseName,
	                       Logger log) {
		BWA bwa = new BWA(bwaLocation, overwriteExisting, verbose, log);
		BWA_Analysis bwa_Analysis =
		                          new BWA_Analysis(rootInputDir, rootOutputDir, referenceGenomeFasta,
		                                           verbose, numMemThreads, numSampleThreads, bwa, log);
		bwa_Analysis.init(fileOfSamplePairs);
		if (batch) {
			bwa_Analysis.batch(numBatches, memoryInMB, wallTimeInHours, baseName);
		} else {
			bwa_Analysis.analyzeBWA_MEM();
		}
	}

	public static void test(String rootDir, String rootOut, String ref) {

		Logger log = new Logger("/home/pankrat2/shared/testGATK/testing.log");
		BWA bwa = new BWA("/home/pankrat2/lanej/bin/bwa", true, true, log);
		// bwa.indexReferenceGenome(ref);
		BWA_Analysis bwa_Analysis = new BWA_Analysis(rootDir, rootOut, ref, true, 8, 1, bwa, log);
		System.out.println("HI");

		bwa_Analysis.init(null);
		System.out.println("HI");

		bwa_Analysis.analyzeBWA_MEM();
		System.out.println("END");

	}

	public static void main(String[] args) {
		// public BWA_Analysis(String rootInputDir, String rootOutputDir, String referenceGenomeFasta,
		// BWA bwa, boolean verbose, int numMemThreads, int numSampleThreads, Logger log) {
		int numArgs = args.length;
		String rootInputDir = null;
		String rootOutputDir = null;
		String referenceGenomeFasta = null;
		String bwaLocation = "";
		String fileOfSamplePairs = null;
		boolean verbose = true;
		boolean batch = false;
		int numBatches = 5;
		int memoryInMB = 23000;
		int wallTimeInHours = 48;
		int numMemThreads = 8;
		int numSampleThreads = 1;
		boolean overwriteExisting = false;
		String baseName = "";

		String logFile = "bwaMem.log";

		String usage = "\n" + "seq.BWA_Analysis requires 2 argument\n";
		usage += "   (1) root input directory (i.e. " + ROOT_INPUT_COMMAND + rootInputDir
		         + " (no default))\n" + "";
		usage += "   (2) root output directory (i.e. " + ROOT_OUTPUT_COMMAND + rootOutputDir
		         + " (no default))\n" + "";
		usage += "   (3) tab-delimited file with no header of paired .fastq (i.e. "
		         + FILE_OF_SAMPLE_PAIRS_COMMAND + fileOfSamplePairs + " (optional, no default))\n" + "";
		usage += "   (4) the full path to a  reference genome in fasta format (i.e."
		         + REFERENCE_GENOME_COMMAND + referenceGenomeFasta + " (no default))\n" + "";
		usage += "   (5) the full path to the bwa executable (i.e. " + BWA_LOCATION_COMMAND
		         + bwaLocation + " (no default, defualts to systems path))\n" + "";
		usage += "   (6) run in quiet mode (i.e. " + QUIET_COMMAND + " (not tbe default))\n" + "";
		usage += "   (7) number of threads for bwa mem (i.e." + NUM_BETWEEN_THREADS_COMMAND
		         + numMemThreads + " (default))\n" + "";
		usage += "   (8) filename for a log (i.e. " + LOG_FILE_COMMAND + logFile + " (default))\n" + "";
		usage += "   (9) set up a batch analysis for the root input directory for a log (i.e. "
		         + BATCH_COMMAND + " (not the default))\n" + "";
		usage += "   (10) number of batches for a batched analysis (i.e. " + NUMBATCHES_COMMAND
		         + numBatches + " (the default))\n" + "";
		usage += "   (11) over-write exsiting files (i.e. " + OVERWRITE_EXISTING_COMMAND
		         + " (not the default))\n" + "";
		usage +=
		      "   (11) base-name for batch analysis (i.e. " + BASE_NAME_COMMAND + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith(ROOT_INPUT_COMMAND)) {
				rootInputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(ROOT_OUTPUT_COMMAND)) {
				rootOutputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(FILE_OF_SAMPLE_PAIRS_COMMAND)) {
				fileOfSamplePairs = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(REFERENCE_GENOME_COMMAND)) {
				referenceGenomeFasta = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(BWA_LOCATION_COMMAND)) {
				bwaLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(LOG_FILE_COMMAND)) {
				logFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(BASE_NAME_COMMAND)) {
				baseName = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(NUM_BETWEEN_THREADS_COMMAND)) {
				numMemThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUM_WITHIN_THREADS_COMMAND)) {
				numSampleThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUMBATCHES_COMMAND)) {
				numBatches = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("memoryInMB=")) {
				memoryInMB = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("wallTimeInHours=")) {
				wallTimeInHours = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(QUIET_COMMAND)) {
				verbose = false;
				numArgs--;
			} else if (arg.startsWith(BATCH_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (arg.startsWith(OVERWRITE_EXISTING_COMMAND)) {
				batch = true;
				numArgs--;
			} else {
				System.err.println("invalid argument: " + arg);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (rootOutputDir == null) {
			rootOutputDir = rootInputDir;
		}
		Logger log = new Logger(rootOutputDir + "bwa.log");
		run(rootInputDir, rootOutputDir, referenceGenomeFasta, bwaLocation, fileOfSamplePairs,
		    overwriteExisting, verbose, numMemThreads, numSampleThreads, batch, numBatches, memoryInMB,
		    wallTimeInHours, baseName, log);
	}
}

// private boolean verifyReferenceGenome() {
// boolean verfied = true;
// if (referenceGenomeFasta == null && referenceGenomeIndexed == null) {
// log.reportError("a reference (indexed or un-indexed) genome must be supplied");
// verfied = false;
// } else {
// if (referenceGenomeIndexed != null) {
// if (Files.exists(referenceGenomeIndexed)) {
// if (referenceGenomeIndexed.endsWith(INDEXED_REFERENCE_GENOME_EXT)) {
// if (verbose) {
// log.report(ext.getTime() + " Info - using pre-indexed reference genome " +
// referenceGenomeIndexed);
// }
// } else {
// log.reportError("Warning - pre-indexed reference genome files for bwa typically end with " +
// INDEXED_REFERENCE_GENOME_EXT + ", " + referenceGenomeIndexed + " did not");
// }
// } else {
// verfied = false;
// log.reportError("indexed reference genome was specified, but did not exist at" +
// referenceGenomeIndexed);
// }
// } else if (referenceGenomeFasta != null) {
// referenceGenomeIndexed = ext.rootOf(referenceGenomeFasta, false) + INDEXED_REFERENCE_GENOME_EXT;
// if (Files.exists(referenceGenomeIndexed)) {
// if (verbose) {
// log.report(ext.getTime() + " Info - detected pre-indexed reference genome " +
// referenceGenomeIndexed);
// }
// } else {
// if (!Files.exists(referenceGenomeFasta)) {
// log.reportError("reference genome was specified, but did not exist at" +
// referenceGenomeFasta);
// verfied = false;
// } else {
// new BWA(bwaLocation, false, verbose, log).indexReferenceGenome(referenceGenomeFasta);
// if (!Files.exists(referenceGenomeIndexed)) {
// log.reportError("could not create indexed reference genome" + referenceGenomeIndexed);
// verfied = false;
// }
// }
// }
// }
// }
// return verfied;
// }
