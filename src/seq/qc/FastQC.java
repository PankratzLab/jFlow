package seq.qc;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import seq.analysis.BWA_Analysis;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

public class FastQC {
	public static final String FASTA_QC_END_MODULE = ">>END_MODULE";
	public static final String FASTA_QC_BEGIN_MODULE = ">>";

	public static final String FAST_QC_RESULTS_FILE = "fastqc_data.txt";
	public static final String FAST_QC_RESULTS_FILENAME_TAG = "Filename";
	public static final String FAST_SUMMARY_FILE = ".fastqc.allSamples.summary.txt";

	public static final String ZIP = ".zip";
	public static final String[] MODULES_TO_EXTRACT = { ">>Sequence Duplication Levels", ">>Adapter Content" };
	public static final String FASTA_QC_LOCATION_COMMAND = "fastaQCLocation=";
	public static final String NUM_THREADS_COMMAND = "numThreads=";
	private String fastQCLocation;
	private String[] rootInputDirs;
	private String rootOutputDir;
	private String[] fastaQFiles;
	private int numThreads;
	private boolean verbose, fail, overWriteExistingOutput;
	private Logger log;

	public FastQC(String fastQCLocation, String[] rootInputDirs, String rootOutputDir, int numThreads, boolean verbose, boolean overWriteExistingOutput, Logger log) {
		super();
		this.fastQCLocation = fastQCLocation;
		this.rootInputDirs = rootInputDirs;
		this.rootOutputDir = (rootOutputDir == null ? rootInputDirs[0] : rootOutputDir);
		this.numThreads = numThreads;
		this.verbose = verbose;
		this.overWriteExistingOutput = overWriteExistingOutput;
		this.fail = !verify(new String[] { fastQCLocation });
		this.log = log;
		if (fail) {
			log.reportError("Error - could not find fastQC location " + fastQCLocation);
		}
	}

	public void fastaQC() {
		gatherFastaQFiles();
		if (!fail) {
			ExecutorService executor = Executors.newFixedThreadPool(numThreads);
			Hashtable<String, Future<Boolean>> tmpResults = new Hashtable<String, Future<Boolean>>();
			for (int i = 0; i < fastaQFiles.length; i++) {
				tmpResults.put(i + "", executor.submit(new FastQCWorkerThread(fastQCLocation, fastaQFiles[i], rootOutputDir, verbose, overWriteExistingOutput, log)));
			}
			for (int i = 0; i < fastaQFiles.length; i++) {
				try {
					boolean success = tmpResults.get(i + "").get();
					if (!success) {
						fail = success;// all results did not come through
					}
				} catch (InterruptedException e) {
					log.reportError("Error - could running fastaQC on internal index " + i);
					log.reportException(e);
					fail = true;
				} catch (ExecutionException e) {
					log.reportError("Error - could running fastaQC on internal index " + i);
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

	public void parseResultsFiles(boolean allModules) {
		if (!fail) {
			String[] zipFiles = Files.toFullPaths(Files.list(rootOutputDir, ZIP, false), rootOutputDir);
			if (zipFiles == null || zipFiles.length < 1) {
				log.reportError("Error - did not find any files with extension" + ZIP + " in output directory " + rootOutputDir);
				fail = true;
			} else {
				FastaQCModuleResults[][] fastaQCModuleResults = new FastaQCModuleResults[zipFiles.length][];
				for (int i = 0; i < zipFiles.length; i++) {
					fastaQCModuleResults[i] = parseFastaQCResult(zipFiles[i], allModules, log);
					if (!allModules && fastaQCModuleResults[i].length != MODULES_TO_EXTRACT.length) {
						log.reportError("Error - did not detect all neccesary modules in zip directory " + zipFiles[i]);
						fail = true;
					}
				}
				if (!fail) {
					for (int i = 0; i < fastaQCModuleResults[0].length; i++) {
						String currentModule = fastaQCModuleResults[0][i].getModuleTitleFormatted();
						String currentOutput = rootOutputDir + currentModule + FAST_SUMMARY_FILE;
						try {
							PrintWriter writer = new PrintWriter(new FileWriter(currentOutput));
							String[] header = fastaQCModuleResults[0][i].getModuleHeader();
							writer.println("Sample\tInternalKey\t" + Array.toStr(header));
							for (int j = 0; j < fastaQCModuleResults.length; j++) {
								String[] allKeys = fastaQCModuleResults[j][i].getAllArrayKeys();
								for (int k = 0; k < allKeys.length; k++) {
									writer.println(fastaQCModuleResults[j][i].getSourceFile() + "\t" + allKeys[k] + "\t" + Array.toStr(fastaQCModuleResults[j][i].getValueForKey(allKeys[k])));
								}
							}
							writer.println();
							writer.close();
						} catch (Exception e) {
							log.reportError("Error writing to " + currentOutput);
							log.reportException(e);
						}
					}

				}
			}
		}
	}

	private static FastaQCModuleResults[] parseFastaQCResult(String fullPathFastQCZipFile, boolean allModules, Logger log) {
		ArrayList<FastaQCModuleResults> fastaQCModuleResults = new ArrayList<FastQC.FastaQCModuleResults>(10);
		try {
			ZipFile zipFile = new ZipFile(fullPathFastQCZipFile);
			Enumeration<? extends ZipEntry> entries = zipFile.entries();
			while (entries.hasMoreElements()) {
				ZipEntry entry = entries.nextElement();
				if (ext.removeDirectoryInfo(entry.getName()).equals(FAST_QC_RESULTS_FILE)) {
					InputStream stream = (InputStream) zipFile.getInputStream(entry);
					BufferedReader reader = new BufferedReader(new InputStreamReader(stream, "UTF-8"));
					boolean scanningToModule = true;
					boolean scanningToHeader = true;
					String inputFileName = "NA";
					String currentModule = "NA";
					String[] currentHeader = null;
					int currentModuleIndex = 0;
					while (reader.ready()) {
						String[] line = reader.readLine().trim().split("\t");
						if (!scanningToModule) {
							scanningToHeader = line[0].startsWith("#");
						}
						if (line[0].equals(FAST_QC_RESULTS_FILENAME_TAG)) {
							inputFileName = line[1];
						}
						if (scanningToHeader || scanningToModule) {
							int index = ext.indexOfStr(line[0], MODULES_TO_EXTRACT, true, true);
							if (scanningToModule && (index >= 0 || (allModules && line[0].startsWith(FASTA_QC_BEGIN_MODULE) && !line[0].equals(FASTA_QC_END_MODULE)))) {
								currentModule = line[0];
								scanningToModule = false;
							}
							if (!scanningToModule && scanningToHeader && line[0].startsWith("#")) {
								currentHeader = line;
								currentHeader[0] = currentHeader[0].substring(1);
								scanningToHeader = false;
							}
						} else if (line[0].equals(FASTA_QC_END_MODULE)) {
							scanningToHeader = true;
							scanningToModule = true;
							if (currentModuleIndex != fastaQCModuleResults.size()) {// nothing to gather
								currentModuleIndex++;
							}
						} else {
							if (currentModuleIndex == fastaQCModuleResults.size()) {
								fastaQCModuleResults.add(new FastaQCModuleResults(inputFileName, currentModule, currentHeader, log));
							}

							fastaQCModuleResults.get(currentModuleIndex).add(line[0], line);
						}

					}
					reader.close();
				}
			}
			zipFile.close();
		} catch (IOException e) {
			log.reportError("Error - could not access zip file " + fullPathFastQCZipFile);
			log.reportException(e);
			e.printStackTrace();
		}
		return fastaQCModuleResults.toArray(new FastaQCModuleResults[fastaQCModuleResults.size()]);
	}

	private void gatherFastaQFiles() {
		if (!fail) {
			fail = !verify(rootInputDirs);
			if (!fail) {
				if (verbose) {
					log.report(ext.getTime() + " Info - gathering files from " + Array.toStr(rootInputDirs, "\n"));
				}
				ArrayList<String> tmpAllFiles = new ArrayList<String>();

				for (int i = 0; i < BWA_Analysis.FQ_EXTS.length; i++) {
					for (int j = 0; j < rootInputDirs.length; j++) {

						String[] tmpFiles = Files.list(rootInputDirs[j], BWA_Analysis.FQ_EXTS[i], false);
						if (tmpFiles != null && tmpFiles.length > 0) {
							if (verbose) {
								log.report("Info - found " + tmpFiles.length + " of type " + BWA_Analysis.FQ_EXTS[i] + " in " + rootInputDirs[j]);
							}
							tmpFiles = Files.toFullPaths(tmpFiles, rootInputDirs[j]);
							for (int k = 0; k < tmpFiles.length; k++) {
								tmpAllFiles.add(tmpFiles[k]);
							}
						}
					}
				}
				this.fastaQFiles = tmpAllFiles.toArray(new String[tmpAllFiles.size()]);
			} else {
				log.reportError("Error - could not find input location(s) " + Array.toStr(rootInputDirs, "\n"));
			}
		}
		if (fastaQFiles == null || fastaQFiles.length < 1) {
			fail = true;
			log.reportError("Error - could not find any input files in " + Array.toStr(rootInputDirs, "\n") + "\nwith any of the following extensions:" + Array.toStr(BWA_Analysis.FQ_EXTS, "\n"));
		} else {
			log.report(ext.getTime() + " Info - found " + fastaQFiles.length + " file(s) to QC");
		}
	}

	public boolean isFail() {
		return fail;
	}

	private static boolean verify(String[] dirs) {
		for (int i = 0; i < dirs.length; i++) {
			if (!Files.exists(dirs[i])) {
				return false;
			}
		}
		return true;
	}

	private static class FastQCWorkerThread implements Callable<Boolean> {
		// private static final String F = "-f";
		private static final String O = "-o";
		private static final String FASTQC_ZIP = "_fastqc.zip";
		private static final String FASTQC_HTML = "_fastqc.html";

		private String fastQCLocation;
		private String fastaQFile;
		private String rootOutputDir;
		private String[] command;
		private boolean verbose, overWriteExistingOutput;
		private Logger log;

		public FastQCWorkerThread(String fastQCLocation, String fastaQFile, String rootOutputDir, boolean verbose, boolean overWriteExistingOutput, Logger log) {
			super();
			this.fastQCLocation = fastQCLocation;
			this.fastaQFile = fastaQFile;
			this.rootOutputDir = rootOutputDir;
			this.overWriteExistingOutput = overWriteExistingOutput;
			this.verbose = verbose;
			this.log = log;
		}

		@Override
		public Boolean call() throws Exception {
			this.command = new String[] { fastQCLocation, fastaQFile, O, rootOutputDir };
			String baseOutput = rootOutputDir + ext.rootOf(fastaQFile);
			if (verbose) {
				log.report(ext.getTime() + " Info - beginning fastaQC for " + fastaQFile);
			}
			boolean success = CmdLine.runCommandWithFileChecks(command, "", new String[] { fastaQFile }, new String[] { baseOutput + FASTQC_ZIP, baseOutput + FASTQC_HTML }, verbose, overWriteExistingOutput, true, log);
			if (verbose && success) {
				log.report(ext.getTime() + " Info - finished fastaQC for " + fastaQFile);
			} else if (!success) {
				log.reportError(ext.getTime() + " Error - could not complete fastaQC for " + fastaQFile);
			}
			return success;
		}
	}

	public static void run(String fastQCLocation, String[] rootInputDirs, String rootOutputDir, int numThreads, boolean verbose, boolean overWriteExistingOutput, boolean batch, Logger log) {
		FastQC fastQC = new FastQC(fastQCLocation, rootInputDirs, rootOutputDir, numThreads, verbose, overWriteExistingOutput, log);
		if (!fastQC.isFail()) {
			fastQC.fastaQC();
		}
		if (!fastQC.isFail()) {

			fastQC.parseResultsFiles(true);
		}
	}

	public static class FastaQCModuleResults {
		private String moduleTitle, sourceFile;
		private String[] moduleHeader;
		private Hashtable<String, String[]> moduleValues;
		private int internalIndex;
		private ArrayList<String> allKeys;

		private Logger log;

		public FastaQCModuleResults(String sourceFile, String moduleTitle, String[] moduleHeader, Logger log) {
			super();
			this.sourceFile = sourceFile;
			this.moduleTitle = moduleTitle;
			this.moduleHeader = moduleHeader;
			this.moduleValues = new Hashtable<String, String[]>();
			this.allKeys = new ArrayList<String>();
			this.internalIndex = 0;
			this.log = log;
		}

		public void add(String key, String[] value) {

			if (key.contains("-")) {
				key = key.split("-")[0];
			}
			if (moduleValues.containsKey(key)) {
				key = internalIndex + "_" + key;
				internalIndex++;
			}
			moduleValues.put(key, value);

			allKeys.add(key);

		}

		public String[] getValueForKey(String key) {
			return getModuleValues().get(key);
		}

		public Hashtable<String, String[]> getModuleValues() {
			return moduleValues;
		}

		public Logger getLog() {
			return log;
		}

		public String getSourceFile() {
			return sourceFile;
		}

		public String[] getAllArrayKeys() {
			return getAllKeys().toArray(new String[getAllKeys().size()]);
		}

		public ArrayList<String> getAllKeys() {
			return allKeys;
		}

		public String getModuleTitleFormatted() {
			return ext.replaceWithLinuxSafeCharacters(getModuleTitle().replaceAll(">", ""), true);
		}

		public String getModuleTitle() {
			return moduleTitle;
		}

		public String[] getModuleHeader() {
			return moduleHeader;
		}

		// private

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		ArrayList<String> rootInputDirs = new ArrayList<String>();
		String rootOutputDir = null;
		boolean overwriteExisting = false;

		String fastQCLocation = "";
		boolean verbose = true;
		boolean batch = false;
		int numThreads = 8;
		// int memoryInMB = 23000;
		// int wallTimeInHours = 48;

		String logFile = "FastaQC.log";

		String usage = "\n" + "seq.FastaQC requires 2 argument\n";
		usage += "   (1) root input directory (if multiple, use the argument for each directory) (i.e. " + BWA_Analysis.ROOT_INPUT_COMMAND + " (no default))\n" + "";
		usage += "   (2) root output directory (i.e. " + BWA_Analysis.ROOT_OUTPUT_COMMAND + rootOutputDir + " (default))\n" + "";
		usage += "   (3) the full path to the fastaQC executable (i.e. " + FASTA_QC_LOCATION_COMMAND + fastQCLocation + " (defualts to systems path))\n" + "";
		usage += "   (4) run in quiet mode (i.e. " + BWA_Analysis.QUIET_COMMAND + " (not tbe default))\n" + "";
		usage += "   (5) number of threads for fastaQC (i.e." + NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (6) filename for a log (i.e. " + BWA_Analysis.LOG_FILE_COMMAND + logFile + " (default))\n" + "";
		usage += "   (7) over-write exsiting files (i.e. " + BWA_Analysis.OVERWRITE_EXISTING_COMMAND + " (not the default))\n" + "";

		// usage += "   (7) set up an analysis for a compute node (i.e. " + BWA_Analysis.BATCH_COMMAND + " (not the default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(BWA_Analysis.ROOT_INPUT_COMMAND)) {
				rootInputDirs.add(ext.parseStringArg(args[i], ""));
				numArgs--;
			} else if (args[i].startsWith(BWA_Analysis.ROOT_OUTPUT_COMMAND)) {
				rootOutputDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(FASTA_QC_LOCATION_COMMAND)) {
				fastQCLocation = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(BWA_Analysis.LOG_FILE_COMMAND)) {
				logFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("memoryInMB=")) {
				// memoryInMB = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("wallTimeInHours=")) {
				// wallTimeInHours = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(BWA_Analysis.QUIET_COMMAND)) {
				verbose = false;
				numArgs--;
			} else if (args[i].startsWith(BWA_Analysis.BATCH_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (args[i].startsWith(BWA_Analysis.OVERWRITE_EXISTING_COMMAND)) {
				overwriteExisting = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Logger log = new Logger(rootOutputDir + logFile);
		run(fastQCLocation, rootInputDirs.toArray(new String[rootInputDirs.size()]), rootOutputDir, numThreads, verbose, overwriteExisting, batch, log);
	}

}
