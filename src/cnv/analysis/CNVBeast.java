package cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * Wrapper for the Beast software package as available here http://people.duke.edu/~asallen/Software.html
 * <p>
 * Currently will parse data (sample-based) to beast formatted data and configuration files, execute them, and consolidate all results in a .cnv file
 * <p>
 * Warning - cnv.beast creates tons of intermediate files and can easily reach .5 million. To estimate the number, take 26*3*number of samples
 * 
 *
 */
public class CNVBeast {
	public static final String OVERWRITE_OPTION = "-overwrite";
	public static final String BEAST_DELIM = "\\s+";

	public static final String[] BEAST_DATA_HEADER = { "probeset_id", "chr", "position", "LRR" };
	public static final String[] BEAST_RESULTS_HEADER = { "start", "end", "height", "length", "type", "score", "start_probe", "end_probe" };
	public static final String[] STRANGE_BEAST_SOURCE = { "***" };

	private Project proj;
	private String fullPathToBeastExe;
	private boolean[] samplesToAnalyze;
	private boolean overWriteExistingFiles;
	private int numThreads;
	private String analysisDirectory;
	private String outputCNVFile;

	/**
	 * @param proj
	 *            Current project
	 * @param fullPathToBeastExe
	 *            will only be valid on Windows
	 * @param samplesToAnalyze
	 *            only these samples will be exported
	 * @param overWriteExistingFiles
	 *            create new data files
	 * @param analysisDirectory
	 *            sub-directory under the project directory
	 * @param outputCNVFile
	 *            relative to the project and analysis directory
	 * @param numThreads
	 *            number of threads for parsing and analysis
	 */
	public CNVBeast(Project proj, String fullPathToBeastExe, boolean[] samplesToAnalyze, boolean overWriteExistingFiles, String analysisDirectory, String outputCNVFile, int numThreads) {
		super();
		this.proj = proj;
		this.fullPathToBeastExe = fullPathToBeastExe;
		this.samplesToAnalyze = samplesToAnalyze;
		this.overWriteExistingFiles = overWriteExistingFiles;
		this.analysisDirectory = analysisDirectory;
		this.numThreads = numThreads;
		this.outputCNVFile = outputCNVFile;
	}

	// TODO, implement a marker based method to generate corrected intensity data
	public void analyzeSampleBased() {
		String[] samples = proj.getSamples();
		if (samplesToAnalyze != null && samples.length != samplesToAnalyze.length) {
			proj.getLog().reportError("Error - mismatched array sizes for project's samples and sample mask");
		} else {
			MarkerSet markerSet = proj.getMarkerSet();
			int numSamples = (samplesToAnalyze != null ? Array.booleanArraySum(samplesToAnalyze) : samples.length);
			BeastConfig[][] sampleConfigs = new BeastConfig[numSamples][];
			int configIndex = 0;
			proj.getLog().report(ext.getTime() + " Initializing configurations for " + numSamples + (numSamples == 1 ? " sample" : " samples"));
			for (int i = 0; i < samples.length; i++) {
				if (i % 100 == 0) {
					proj.getLog().report(ext.getTime() + " Info - initializing sample " + (i + 1) + "/" + samples.length);
				}
				if (samplesToAnalyze == null || samplesToAnalyze[i]) {
					sampleConfigs[configIndex] = generateConfigs(proj, samples[i], markerSet, proj.getProjectDir() + analysisDirectory, fullPathToBeastExe, overWriteExistingFiles);
					configIndex++;
				}
			}
			BeastConfig[] allConfigs = parseAllSamples(proj, sampleConfigs, markerSet, numThreads);
			allConfigs = analyzeConfigs(proj, allConfigs, numThreads);
			CNVariant[][] allCNVs = parseBeastResults(proj, allConfigs);
			reportAllCNVs(proj, allCNVs, proj.getProjectDir() + analysisDirectory + outputCNVFile);
		}
	}

	private static void reportAllCNVs(Project proj, CNVariant[][] allCNVs, String fullPathToOutput) {
		proj.getLog().report(ext.getTime() + " Info - reporting all cnvs to " + fullPathToOutput);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i < allCNVs.length; i++) {
				if (allCNVs[i] != null) {
					for (int j = 0; j < allCNVs[i].length; j++) {
						writer.println(allCNVs[i][j].toPlinkFormat());
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + fullPathToOutput);
			proj.getLog().reportException(e);
		}
	}

	private static CNVariant[][] parseBeastResults(Project proj, BeastConfig[] allConfigs) {
		SampleData sampleData = proj.getSampleData(0, false);
		CNVariant[][] allCNVs = new CNVariant[allConfigs.length][];
		for (int i = 0; i < allConfigs.length; i++) {
			if (allConfigs[i].hasSummaryFile()) {
				allCNVs[i] = loadBeastResultsToCNVariant(proj, sampleData, allConfigs[i]);
			} else {
				allCNVs[i] = new CNVariant[0];
				proj.getLog().reportError("Error - the summary file for sample " + allConfigs[i].getSummaryFile() + " does not exist and will not be included in the final output");
			}
		}
		return allCNVs;
	}

	/**
	 * Finds each beast result file (for each sample and each chromosome) and loads it to a {@link CNVariant}[][]
	 */
	private static CNVariant[] loadBeastResultsToCNVariant(Project proj, SampleData sampleData, BeastConfig config) {
		ArrayList<CNVariant> indCNVariants = new ArrayList<CNVariant>();
		String[] ind = sampleData.lookup(config.getSample());
		String summaryFile = config.getSummaryFile();
		byte chr = (byte) config.getAnalysisChr();
		if (ind != null) {
			String fid = ind[1].split("\t")[0];
			String iid = ind[1].split("\t")[1];
			try {
				BufferedReader reader = Files.getAppropriateReader(summaryFile);
				String[] beastHeader = reader.readLine().trim().split(BEAST_DELIM);
				int[] indices = ext.indexFactors(beastHeader, BEAST_RESULTS_HEADER, true, false);
				if (Array.countIf(indices, -1) == 0) {
					while (reader.ready()) {
						String[] beastLine = reader.readLine().trim().split(BEAST_DELIM);
						try {
							int start = Integer.parseInt(beastLine[indices[0]]);
							int stop = Integer.parseInt(beastLine[indices[1]]);
							int cn = determineCN(Double.parseDouble(beastLine[indices[2]]));
							int numMarkers = Integer.parseInt(beastLine[indices[3]]);
							int source = -99;
							// TODO, figure out the weird types, and the case where numMarkers and types seem to be switched
							if (ext.indexOfStr(beastLine[indices[4]], STRANGE_BEAST_SOURCE) < 0) {// sometimes this entry is goofy
								source = Integer.parseInt(beastLine[indices[4]]);
							}
							double score = Double.parseDouble(beastLine[indices[5]]);
							indCNVariants.add(new CNVariant(fid, iid, chr, start, stop, cn, score, numMarkers, source));
						} catch (NumberFormatException nfe) {
							proj.getLog().reportError("Error - results file " + summaryFile + " had an improper number on line " + Array.toStr(beastLine));
						}
					}
					reader.close();
				} else {
					proj.getLog().reportError("Error - could not detect proper header in results file " + summaryFile);
				}
			} catch (FileNotFoundException fnfe) {
				proj.getLog().reportError("Error: file \"" + summaryFile + "\" not found in current directory");

			} catch (IOException ioe) {
				proj.getLog().reportError("Error reading file \"" + summaryFile + "\"");
			}
		} else {
			proj.getLog().reportError("Error - could look up sample " + config.getSample() + " in sample data file " + proj.getFilename(Project.SAMPLE_DATA_FILENAME));
		}
		return indCNVariants.toArray(new CNVariant[indCNVariants.size()]);
	}

	private static int determineCN(double height) {
		return height > 0 ? 3 : 1;
	}

	private static BeastConfig[] parseAllSamples(Project proj, BeastConfig[][] configs, MarkerSet markerSet, int numThreads) {
		ArrayList<BeastConfig> allConfigs = new ArrayList<BeastConfig>();
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		ArrayList<Future<BeastConfig[]>> tmpResults = new ArrayList<Future<BeastConfig[]>>();
		for (int i = 0; i < configs.length; i++) {// need to submit the jobs first
			if (configs[i] != null && configs.length > 0) {
				BeastSampleParserThread worker = new BeastSampleParserThread(proj, configs[i], configs[i][0].getSample(), markerSet.getIndicesByChr(), markerSet.getMarkerNames(), markerSet.getPositions());
				tmpResults.add(executor.submit(worker));
			}
		}
		for (int i = 0; i < configs.length; i++) {
			try {
				BeastConfig[] tmpConfigs = tmpResults.get(i).get();// get is only applied after the job has finished
				for (int j = 0; j < tmpConfigs.length; j++) {
					allConfigs.add(tmpConfigs[j]);
				}
			} catch (InterruptedException e) {
				proj.getLog().reportException(e);
			} catch (ExecutionException e) {
				proj.getLog().reportException(e);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(7, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			return allConfigs.toArray(new BeastConfig[allConfigs.size()]);
		}
		return allConfigs.toArray(new BeastConfig[allConfigs.size()]);
	}

	private static BeastConfig[] analyzeConfigs(Project proj, BeastConfig[] configs, int numThreads) {
		BeastConfig[] analyzedConfigs = new BeastConfig[configs.length];
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		ArrayList<Future<BeastConfig>> tmpResults = new ArrayList<Future<BeastConfig>>();
		for (int i = 0; i < configs.length; i++) {
			BeastAnalysisThread worker = new BeastAnalysisThread(proj, configs[i]);
			tmpResults.add(executor.submit(worker));
		}
		for (int i = 0; i < configs.length; i++) {
			try {
				analyzedConfigs[i] = tmpResults.get(i).get();
			} catch (InterruptedException e) {
				proj.getLog().reportException(e);
			} catch (ExecutionException e) {
				proj.getLog().reportException(e);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(7, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			return analyzedConfigs;
		}
		return analyzedConfigs;
	}

	private static BeastConfig[] generateConfigs(Project proj, String sample, MarkerSet markerSet, String analysisDirectoryFullPath, String fullPathToBeastExe, boolean overWriteExistingFiles) {
		int[][] indicesByChr = markerSet.getIndicesByChr();
		ArrayList<BeastConfig> configs = new ArrayList<BeastConfig>(indicesByChr.length);
		for (int i = 1; i < indicesByChr.length; i++) {// skip 0
			if (hasDataAt(indicesByChr, i)) {
				String baseName = sample + "_chr" + i;
				String subDir = "chr" + i + "/";
				BeastConfig sampChrConfig = new BeastConfig(sample, baseName, fullPathToBeastExe, analysisDirectoryFullPath + subDir, overWriteExistingFiles, i, markerSet.getMarkerNames().length, proj.getLog());
				configs.add(sampChrConfig);
			}
		}
		return configs.toArray(new BeastConfig[configs.size()]);
	}

	private static boolean hasDataAt(int[][] indicesByChr, int index) {
		return indicesByChr[index] != null && indicesByChr[index].length > 0;
	}

	private static class BeastSampleParserThread implements Callable<BeastConfig[]> {
		private Project proj;
		private BeastConfig[] configs;
		private String sampleToParse;
		private int[][] indicesByChr;
		private String[] markerNames;
		private int[] positions;

		public BeastSampleParserThread(Project proj, BeastConfig[] configs, String sampleToParse, int[][] indicesByChr, String[] markerNames, int[] positions) {
			super();
			this.proj = proj;
			this.configs = configs;
			this.sampleToParse = sampleToParse;
			this.indicesByChr = indicesByChr;
			this.markerNames = markerNames;
			this.positions = positions;
		}

		@Override
		public BeastConfig[] call() throws Exception {
			if (verifySameSample(sampleToParse, configs, proj.getLog())) {
				Sample samp = null;
				float[] lrrs = null;
				for (int i = 0; i < configs.length; i++) {
					if (!Files.exists(configs[i].getDataFile())) {
						if (samp == null || lrrs == null) {// only load once, and only if needed
							samp = proj.getFullSampleFromRandomAccessFile(sampleToParse);
							lrrs = samp.getLRRs();
						}
						try {
							PrintWriter writer = new PrintWriter(new FileWriter(configs[i].getDataFile()));
							writer.println(Array.toStr(BEAST_DATA_HEADER));
							int analysisChr = configs[i].getAnalysisChr();
							for (int j = 0; j < indicesByChr[analysisChr].length; j++) {
								writer.println(markerNames[indicesByChr[analysisChr][j]] + "\t" + i + "\t" + positions[indicesByChr[analysisChr][j]] + "\t" + lrrs[indicesByChr[analysisChr][j]]);
							}
							writer.close();
						} catch (Exception e) {
							proj.getLog().reportError("Error writing to " + configs[i].getDataFile());
							proj.getLog().reportException(e);
						}
					} else {
						proj.getLog().report(ext.getTime() + " Info - skipping parsing for sample " + configs[0].getSample() + " chromosome:," + configs[i].getAnalysisChr() + " data file " + configs[0].getDataFile() + " already exists");
					}
				}
				return configs;
			} else {
				return null;
			}
		}

		private static boolean verifySameSample(String sample, BeastConfig[] configs, Logger log) {
			for (int i = 0; i < configs.length; i++) {
				if (configs[i].getSample() == null || sample == null || !sample.equals(configs[i].getSample())) {
					log.reportError("Error - mismatched samples in config array, this should not happen");
					return false;
				}
			}
			return true;
		}

	}

	private static class BeastAnalysisThread implements Callable<BeastConfig> {
		private Project proj;
		private BeastConfig config;

		public BeastAnalysisThread(Project proj, BeastConfig config) {
			super();
			this.proj = proj;
			this.config = config;
		}

		@Override
		public BeastConfig call() throws Exception {
			long time = System.currentTimeMillis();
			proj.getLog().report(ext.getTime() + " Info - beginning analysis for sample " + config.getSample() + " on chromosome " + config.getAnalysisChr() + " with thread" + Thread.currentThread().getName());
			config.executeConfig();
			proj.getLog().report(ext.getTime() + " Info - finished analysis in " + ext.getTimeElapsed(time) + " for sample " + config.getSample() + " on chromosome " + config.getAnalysisChr());
			return config;
		}
	}

	private static class BeastConfig {

		private static final int DEFAULT_MIN_GAP = 6;
		private static final int DEFAULT_MAX_GAP = 30;
		// private static final int DEFAULT_BLOCK_SIZE = 50000;
		private static final int DEFAULT_NUM_DELIM_BEFORE_LOC = 2;
		private static final int DEFAULT_NUM_DELIM_BEFORE_LRR = 3;
		private static final double DEFAULT_MIN_RATIO = 0.25;
		private static final double DEFAULT_SCORE_EXPONENT = 0.5;
		private static final int DEFAULT_ASCII_DELIM = 9;// tab
		private static final boolean DEFAULT_QUANTILE = true;
		private static final boolean DEFAULT_OUTPUT_PREDRES = false;
		private static final boolean DEFAULT_SORTED_FILE = true;
		private static final String[] OUTPUT_EXT = { ".config", ".dat", ".summary.res", ".pred.res" };

		private int minGap, maxGap, numDelimBeforeLoc, numDelimBeforeLRR, blockSize, asciiDelim, analysisChr;
		private double minRatio, scoreExponent;
		private boolean quantile, outputPredicted, overWriteExistingFiles, sortedFile;

		private String sample, baseName, summaryFile, predictedResFile, dataFile, configFile, fullPathToBeastExe, analysisDirectoryFullPath;
		private Logger log;

		/**
		 * @param minGap
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param maxGap
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param blockSize
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param numDelimBeforeLoc
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param numDelimBeforeLRR
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param minRatio
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param scoreExponent
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param asciiDelim
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param quantileNormalize
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param outputPredicted
		 *            refer to http://people.duke.edu/~asallen/Software.html
		 * @param baseName
		 *            the base name of the analysis, the {@link #OUTPUT_EXT} extensions will be tagged on
		 * @param analysisDirectoryFullPath
		 *            directory where results will be stored
		 * @param fullPathToBeastExe
		 * @param overWriteExistingFiles
		 *            overwrite the results of a previous run
		 * @param analysisChr
		 *            used for tracking since beast does not require it
		 * @param log
		 */
		public BeastConfig(int minGap, int maxGap, int blockSize, int numDelimBeforeLoc, int numDelimBeforeLRR, double minRatio, boolean sortedFile, double scoreExponent, int asciiDelim, boolean quantileNormalize, boolean outputPredicted, String baseName, String analysisDirectoryFullPath, String fullPathToBeastExe, boolean overWriteExistingFiles, int analysisChr, String sample, Logger log) {
			super();
			this.minGap = minGap;
			this.maxGap = maxGap;
			this.blockSize = blockSize;
			this.numDelimBeforeLoc = numDelimBeforeLoc;
			this.numDelimBeforeLRR = numDelimBeforeLRR;
			this.minRatio = minRatio;
			this.sortedFile = sortedFile;
			this.scoreExponent = scoreExponent;
			this.asciiDelim = asciiDelim;
			this.quantile = quantileNormalize;
			this.outputPredicted = outputPredicted;
			this.baseName = baseName;
			this.fullPathToBeastExe = fullPathToBeastExe;
			this.analysisDirectoryFullPath = analysisDirectoryFullPath;
			this.overWriteExistingFiles = overWriteExistingFiles;
			this.analysisChr = analysisChr;
			this.sample = sample;
			this.log = log;
			initOutputs();

		}

		/**
		 * Constructor to populate with default analysis params for a tab-delimited input files produced by {@link CNVBeast}
		 */
		public BeastConfig(String sample, String baseName, String fullPathToBeastExe, String analysisDirectoryFullPath, boolean overWriteExistingFiles, int analysisChr, int blockSize, Logger log) {
			this(DEFAULT_MIN_GAP, DEFAULT_MAX_GAP, blockSize, DEFAULT_NUM_DELIM_BEFORE_LOC, DEFAULT_NUM_DELIM_BEFORE_LRR, DEFAULT_MIN_RATIO, DEFAULT_SORTED_FILE, DEFAULT_SCORE_EXPONENT, DEFAULT_ASCII_DELIM, DEFAULT_QUANTILE, DEFAULT_OUTPUT_PREDRES, baseName, analysisDirectoryFullPath, fullPathToBeastExe, overWriteExistingFiles, analysisChr, sample, log);
		}

		/**
		 * This executes cnv.beast.exe in a given directory and configuration
		 * 
		 * @return true if successfully produced output
		 */
		public boolean executeConfig() {
			writeConfigFile();
			boolean complete = false;
			if (readyToGo(overWriteExistingFiles)) {
				CmdLine.run(fullPathToBeastExe + " " + configFile, analysisDirectoryFullPath);
			}
			if (!Files.exists(summaryFile)) {
				log.reportError("Error - the result file " + summaryFile + " could not be found (the analysis failed)");
			} else {
				complete = true;
			}
			return complete;
		}

		/**
		 * Make sure the neccesary pieces are already available
		 */
		private boolean readyToGo(boolean overWriteExistingFiles) {
			boolean ready = true;
			if (!Files.exists(fullPathToBeastExe)) {
				log.reportError("Error - could not find the beast executable " + fullPathToBeastExe);
				ready = false;
			}
			if (!Files.exists(configFile)) {
				log.reportError("Error - could not find the configuration file " + configFile);
				ready = false;
			}
			if (!Files.exists(dataFile)) {
				log.reportError("Error - could not find the data file " + dataFile);
				ready = false;
			}
			if (hasSummaryFile() && !overWriteExistingFiles) {
				log.reportError("Warning - the summary file " + summaryFile + " exists and " + OVERWRITE_OPTION + " was not flagged, skipping this beast run");
				ready = false;
			}
			return ready;
		}

		public boolean hasSummaryFile() {
			return Files.exists(summaryFile);
		}

		public int getAnalysisChr() {
			return analysisChr;
		}

		public String getSample() {
			return sample;
		}

		/**
		 * @return a String with the parameters for a beast configuration file
		 */
		public String getConfigs() {
			String config = "";
			config += minGap + "\n";
			config += maxGap + "\n";
			config += minRatio + "\n";
			config += blockSize + "\n";
			config += numDelimBeforeLoc + "\n";
			config += numDelimBeforeLRR + "\n";
			config += dataFile + "\n";
			config += (sortedFile ? 1 : 0) + "\n";
			config += scoreExponent + "\n";
			config += asciiDelim + "\n";
			config += (quantile ? 1 : 0) + "\n";
			config += (outputPredicted ? 1 : 0) + "\n";
			config += summaryFile + "\n";
			config += predictedResFile;
			return config;
		}

		public String getConfigFile() {
			return configFile;
		}

		public void setConfigFile(String configFile) {
			this.configFile = configFile;
		}

		public String getSummaryFile() {
			return summaryFile;
		}

		public void setSummaryFile(String summaryFile) {
			this.summaryFile = summaryFile;
		}

		// public String getPredictedResFile() {
		// return predictedResFile;
		// }

		public void setPredictedResFile(String predictedResFile) {
			this.predictedResFile = predictedResFile;
		}

		public String getDataFile() {
			return dataFile;
		}

		public void setDataFile(String dataFile) {
			this.dataFile = dataFile;
		}

		private void writeConfigFile() {
			initOutputs();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(getConfigFile()));
				writer.print(getConfigs());
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + dataFile);
				log.reportException(e);
			}
		}

		private void initOutputs() {
			new File(analysisDirectoryFullPath).mkdirs();
			setSummaryFile(analysisDirectoryFullPath + baseName + OUTPUT_EXT[2]);
			if (outputPredicted) {
				setPredictedResFile(analysisDirectoryFullPath + baseName + OUTPUT_EXT[3]);
			} else {
				setPredictedResFile("NA");
			}
			setDataFile(analysisDirectoryFullPath + baseName + OUTPUT_EXT[1]);
			setConfigFile(analysisDirectoryFullPath + baseName + OUTPUT_EXT[0]);
		}

	}

	private static boolean[] getSubset(Project proj, String samplesToAnalyzeFile) {
		boolean[] samplesToUse = null;
		if (samplesToAnalyzeFile == null) {
			return samplesToUse;
		} else {
			samplesToUse = new boolean[proj.getSamples().length];
			Arrays.fill(samplesToUse, false);
			String[] samps = HashVec.loadFileToStringArray(proj.getProjectDir() + samplesToAnalyzeFile, false, new int[] { 0 }, true);
			int[] indices = ext.indexLargeFactors(samps, proj.getSamples(), true, proj.getLog(), true, false);
			for (int i = 0; i < indices.length; i++) {
				if (indices[i] >= 0) {
					samplesToUse[indices[i]] = true;
				} else {
					proj.getLog().report("Warning - " + samps[i] + " was not found in the current project, skipping");
				}
			}
			return samplesToUse;
		}
	}

	public static void analyze(Project proj, String fullPathToBeastExe, boolean overWriteExistingFiles, String analysisDirectory, String outputCNVFile, String samplesToAnalyzeFile, int numThreads) {
		CNVBeast cnvBeast = new CNVBeast(proj, fullPathToBeastExe, getSubset(proj, samplesToAnalyzeFile), overWriteExistingFiles, analysisDirectory, outputCNVFile, numThreads);
		cnvBeast.analyzeSampleBased();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String fullPathToBeastExe = "C:/bin/BEAST/cnv.beast.exe";
		boolean overWriteExistingFiles = false;
		String analysisDirectory = "BEAST/";
		String samplesToAnalyzeFile = null;
		String outputCNVFile = "beast.cnv";
		String logfile = null;
		int numThreads = 1;
		// TODO, add customizations to the config file, currently running with beast defaults
		String usage = "\n" + "cnv.analysis.CNVBeast requires 1 argument\n";
		usage += "   (1) project filename (i.e. file=" + filename + " (no default))\n" + "";
		usage += "   (2) full path to beast.exe (i.e. beast=" + fullPathToBeastExe + " (default))\n" + "";// if you can get this to work on linux, this same argument could be used in that environment
		usage += "   (3) analysis subdirectory under the project directory (i.e. dir=" + analysisDirectory + " (default))\n" + "";
		usage += "   (4) overwrite existing result files (i.e. " + OVERWRITE_OPTION + " (default))\n" + "";
		usage += "   (5) int number of threads to use (i.e. numThreads=" + numThreads + " (default))\n" + "";
		usage += "   (6) output filename (relative to the project directory, and analysis directory if supplied) (i.e. output=" + outputCNVFile + " (default))\n" + "";
		usage += "   (7) name of file containing DNA id of individuals to analyze (relative to the project directory) (i.e. samples=" + samplesToAnalyzeFile + " (no default))\n" + "";
		usage += "   (8) name of a log file (i.e. log=" + logfile + " (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("beast=")) {
				fullPathToBeastExe = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("samples=")) {
				samplesToAnalyzeFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				analysisDirectory = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				outputCNVFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(OVERWRITE_OPTION)) {
				overWriteExistingFiles = true;
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, logfile, false);
			analyze(proj, fullPathToBeastExe, overWriteExistingFiles, analysisDirectory, outputCNVFile, samplesToAnalyzeFile, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
