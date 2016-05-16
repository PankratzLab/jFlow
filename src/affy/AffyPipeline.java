package affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.CentroidBuilder;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import cnv.manage.Markers;
import cnv.manage.SourceFileParser;
import cnv.manage.TransposeData;
import cnv.var.SampleData;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.ext;

/**
 * @author lane0212
 *
 *         Maybe the last affy raw data processor.
 *         <p>
 *         Takes us from .CEL -> to Genvisis. Requires affy power tools and related library files
 *         <p>
 * 
 *         Basically follows http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/
 *
 */
public class AffyPipeline {

	private static final String AFFY_CEL_LIST_HEADER = "cel_files";
	private static final String AFFY_PROBELIST_HEADER = "probeset_id";

	private String aptExeDir;// holds "apt-geno-qc", "apt-probeset-genotype", "apt-probeset-summarize" etc
	private String aptLibDir;
	private boolean full;
	private Logger log;

	public AffyPipeline(String aptExeDir, String aptLibDir, boolean full,Logger log) {
		super();
		this.aptExeDir = aptExeDir;
		this.aptLibDir = aptLibDir;
		this.full=full;
		this.log = log;
		if (!new File(aptExeDir).exists()) {
			log.reportTimeError(aptExeDir + " did not exist (Affy exe directory)");
			throw new IllegalArgumentException();
		}
		if (!new File(aptLibDir).exists()) {
			log.reportTimeError(aptLibDir + " did not exist (Affy lib directory)");
			throw new IllegalArgumentException();
		}
		validatePreReq();
		log.reportTimeInfo("Validated Affy Power Tools initialization");
	}

	private void validatePreReq() {

		for (int i = 0; i < AFFY_LIB_FILES.values().length; i++) {
			if (!Files.exists(aptLibDir + AFFY_LIB_FILES.values()[i].getLibFile(full))) {
				log.reportTimeError(aptLibDir + AFFY_LIB_FILES.values()[i].getLibFile(full) + " did not exist");
				throw new IllegalArgumentException();
			}
		}

		for (int i = 0; i < AFFY_ANALYSIS_TYPES.values().length; i++) {
			if (!Files.exists(aptExeDir + AFFY_ANALYSIS_TYPES.values()[i].getExe())) {
				log.reportTimeError(aptExeDir + AFFY_ANALYSIS_TYPES.values()[i].getExe() + " did not exist");
				throw new IllegalArgumentException();
			}
		}

	}

	/**
	 * Required library files for our pipeline
	 *
	 */
	private static enum AFFY_LIB_FILES {
		GW6_CDF("GenomeWideSNP_6.cdf"),
		GW6_BIRDSEED_MODELS("GenomeWideSNP_6.birdseed-v2.models"),
		GW6_SPECIAL_SNPS("GenomeWideSNP_6.specialSNPs"),
		GW6_CHRX("GenomeWideSNP_6.chrXprobes"),
		GW6_CHRY("GenomeWideSNP_6.chrYprobes");

		private String libFile;

		private AFFY_LIB_FILES(String libFile) {
			this.libFile = libFile;
		}

		public String getLibFile(boolean full) {
			if (full) {
				switch (this) {
				case GW6_CDF:
				case GW6_SPECIAL_SNPS:
					return ext.addToRoot(libFile, ".Full");
				case GW6_CHRY:
				case GW6_BIRDSEED_MODELS:
				case GW6_CHRX:
				default:
					break;

				}
			}
			return libFile;

		}
	}

	/**
	 * The basic analyses we perform
	 *
	 */
	private static enum AFFY_ANALYSIS_TYPES {

		GENERATE_PROBE_LIST("apt-probeset-summarize"),
		GENOTYPE("apt-probeset-genotype"),
		NORMALIZE("apt-probeset-summarize");

		private String exe;

		private AFFY_ANALYSIS_TYPES(String exe) {
			this.exe = exe;
		}

		public String getExe() {
			if (Files.isWindows()) {
				return exe + ".exe";
			} else {
				return exe; // *nix systems are much superior
			}
		}

	}

	private String generateCelList(String[] celFiles, String outDir, String analysisName) {
		String celListFile = outDir + analysisName + ".celList.txt";
		String[] toWrite = new String[] { AFFY_CEL_LIST_HEADER };
		toWrite = Array.concatAll(toWrite, celFiles);
		Files.writeList(toWrite, celListFile);
		return celListFile;
	}

	private static class Probesets {
		private String snpOnlyFile;
		private String allFile;
		private boolean fail;

		public Probesets(String snpOnlyFile, String allFile) {
			super();
			this.snpOnlyFile = snpOnlyFile;
			this.allFile = allFile;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getSnpOnlyFile() {
			return snpOnlyFile;
		}

		public String getAllFile() {
			return allFile;
		}

	}

	/**
	 * run a simple command to extract all probesets, return a formatted file to use for downstream analysis
	 */
	private Probesets getAnalysisProbesetList(String celFile, String outDir, String analysisName, String markerPositionFile) {
		String smallCelList = outDir + "tmp" + ext.getTimestampForFilename();
		String currentAnalysis = analysisName + "_probelistGenerator";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(smallCelList));
			writer.println(AFFY_CEL_LIST_HEADER);
			writer.println(celFile);
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to .cel list " + smallCelList);
			e.printStackTrace();
		}

		ArrayList<String> psetCommand = new ArrayList<String>();
		psetCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.GENERATE_PROBE_LIST.getExe());
		psetCommand.add("--cdf-file");
		psetCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile(full));
		psetCommand.add("--analysis");
		psetCommand.add("quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true");
		psetCommand.add("--out-dir");
		psetCommand.add(outDir);
		psetCommand.add("--cel-files");
		psetCommand.add(smallCelList);
		psetCommand.add("--set-analysis-name");
		psetCommand.add(currentAnalysis);

		log.report(ext.getTime() + " Info - running a command to extract probeset ids: " + psetCommand);

		String probeResults = outDir + currentAnalysis + ".summary.txt";
		boolean progress = CmdLine.runCommandWithFileChecks(Array.toStringArray(psetCommand), "", new String[] { celFile }, new String[] { probeResults }, true, false, false, log);

		log.reportTimeInfo("Parsing " + probeResults + " to obtain all probesIds");
		ArrayList<String> probesetIdsAll = new ArrayList<String>(1800000);
		ArrayList<String> probesetIdsSNP = new ArrayList<String>(925000);

		probesetIdsAll.add(AFFY_PROBELIST_HEADER);
		probesetIdsSNP.add(AFFY_PROBELIST_HEADER);
		String tmpMarkerSet = outDir + analysisName + "tmpMarkerSet.ser";
		Markers.orderMarkers(null, markerPositionFile, tmpMarkerSet, log);
		MarkerSet markerSet = MarkerSet.load(tmpMarkerSet, false);
		String[] names = markerSet.getMarkerNames();

		HashMap<String, String> track = new HashMap<String, String>();
		ArrayList<String> markersNotUsed = new ArrayList<String>();

		for (int i = 0; i < names.length; i++) {
			track.put(names[i], names[i]);
		}
		try {
			BufferedReader reader = Files.getAppropriateReader(probeResults);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				if (line[0].startsWith("CN_") || line[0].startsWith("SNP_") || line[0].startsWith("AFFX-SNP")) {
					String pId = line[0];
					if (!pId.endsWith("-B")) {
						String probeParsed = pId.replaceAll("-A", "");
						if (track.containsKey(probeParsed)) {
							probesetIdsAll.add(probeParsed);
							if (!pId.startsWith("CN_")) {
								probesetIdsSNP.add(probeParsed);
							}
						} else {
							markersNotUsed.add(probeParsed);
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + probeResults + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + probeResults + "\"");
			return null;
		}

		log.reportTimeInfo("Detected " + (probesetIdsAll.size() - 1) + " total probesets with " + (probesetIdsSNP.size() - 1) + " SNP_ probesets");
		if (markersNotUsed.size() > 0) {
			log.reportTimeInfo(markersNotUsed.size() + " markers where skipped");
			String pIdSkipFile = outDir + analysisName + ".probesetIdsSkipped.txt";
			Files.writeList(Array.toStringArray(markersNotUsed), pIdSkipFile);

		}
		String pIdAllFile = outDir + analysisName + ".probesetIdsAll.txt";
		String pIdSnpFile = outDir + analysisName + ".probesetIdsSNPS.txt";
		Files.writeList(Array.toStringArray(probesetIdsAll), pIdAllFile);
		Files.writeList(Array.toStringArray(probesetIdsSNP), pIdSnpFile);

		new File(smallCelList).delete();
		Probesets probesets = new Probesets(pIdSnpFile, pIdAllFile);
		probesets.setFail(!progress);
		return probesets;
	}

	private static class NormalizationResult {
		private String quantNormFile;
		private boolean fail;

		public NormalizationResult(String quantNormFile) {
			super();
			this.quantNormFile = quantNormFile;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getQuantNormFile() {
			return quantNormFile;
		}

	}

	private NormalizationResult normalize(String celListFile, String pIDFile, String analysisName, String outDirRoot, String targetSketch) {
		String outCurrent = outDirRoot + analysisName + "_Normalization/";
		new File(outCurrent).mkdirs();
		ArrayList<String> normalizeCommand = new ArrayList<String>();
		normalizeCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.NORMALIZE.getExe());
		normalizeCommand.add("-cdf-file");
		normalizeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile(full));
		normalizeCommand.add("--analysis");
		normalizeCommand.add("quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true");
		normalizeCommand.add("--target-sketch");
		normalizeCommand.add(targetSketch);
		normalizeCommand.add("--out-dir");
		normalizeCommand.add(outCurrent);
		normalizeCommand.add("--cel-files");
		normalizeCommand.add(celListFile);
		normalizeCommand.add("--probeset-ids");
		normalizeCommand.add(pIDFile);
		normalizeCommand.add("--set-analysis-name");
		normalizeCommand.add(analysisName);

		String quantNormFile = outCurrent + analysisName + ".summary.txt";
		String report = outCurrent + analysisName + ".report.txt";

		boolean progress = CmdLine.runCommandWithFileChecks(Array.toStringArray(normalizeCommand), "", null, new String[] { quantNormFile, report }, true, false, false, log);
		NormalizationResult normalizationResult = new NormalizationResult(quantNormFile);
		normalizationResult.setFail(!progress);
		return normalizationResult;
	}

	private static class GenotypeResult {
		private String callFile;
		private String confFile;
		private boolean failed;

		public boolean isFailed() {
			return failed;
		}

		public void setFailed(boolean failed) {
			this.failed = failed;
		}

		public GenotypeResult(String callFile, String confFile) {
			super();
			this.callFile = callFile;
			this.confFile = confFile;
		}

		public String getCallFile() {
			return callFile;
		}

		public String getConfFile() {
			return confFile;
		}

	}

	private GenotypeResult genotype(String celListFile, String pIDFile, String analysisName, String outDirRoot) {

		String outCurrent = outDirRoot + analysisName + "_Genotypes/";
		new File(outCurrent).mkdirs();
		ArrayList<String> genotypeCommand = new ArrayList<String>();
		genotypeCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.GENOTYPE.getExe());
		genotypeCommand.add("-c");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile(full));
		genotypeCommand.add("--table-output");
		genotypeCommand.add("true");
		genotypeCommand.add("-a");
		genotypeCommand.add("birdseed-v2");
		genotypeCommand.add("--set-gender-method");
		genotypeCommand.add("cn-probe-chrXY-ratio");
		genotypeCommand.add("--read-models-birdseed");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_BIRDSEED_MODELS.getLibFile(full));
		genotypeCommand.add("--special-snps");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_SPECIAL_SNPS.getLibFile(full));
		genotypeCommand.add("--chrX-probes");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CHRX.getLibFile(full));
		genotypeCommand.add("--chrY-probes");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CHRY.getLibFile(full));
		genotypeCommand.add("--probeset-ids");
		genotypeCommand.add(pIDFile);
		genotypeCommand.add("--set-analysis-name");
		genotypeCommand.add(analysisName);
		genotypeCommand.add("-out-dir");
		genotypeCommand.add(outCurrent);
		genotypeCommand.add("--cel-files");
		genotypeCommand.add(celListFile);

		String callFile = outCurrent + analysisName + ".calls.txt";
		String confFile = outCurrent + analysisName + ".confidences.txt";
		String reportFile = outCurrent + analysisName + ".report.txt";

		GenotypeResult genotypeResult = new GenotypeResult(callFile, confFile);
		String[] output = new String[] { genotypeResult.getCallFile(), genotypeResult.getConfFile(), reportFile };
		boolean progress = CmdLine.runCommandWithFileChecks(Array.toStringArray(genotypeCommand), "", null, output, true, false, false, log);
		genotypeResult.setFailed(!progress);
		return genotypeResult;
	}

	private static void validateCelSelection(String[] celFiles, Logger log) {
		HashMap<String, String> uniq = new HashMap<String, String>();
		boolean error = false;
		for (int i = 0; i < celFiles.length; i++) {
			if (uniq.containsKey(ext.removeDirectoryInfo(celFiles[i]))) {
				log.reportTimeError(ext.removeDirectoryInfo(celFiles[i]) + " was seen multiple times, perhaps across directories");
				error = true;
			} else {
				uniq.put(ext.removeDirectoryInfo(celFiles[i]), ext.removeDirectoryInfo(celFiles[i]));
			}
			if (celFiles[i].length() != celFiles[i].replaceAll("[\\s]+", "").length()) {
				log.reportTimeError(celFiles[i] + " contained whitespace");
			}
		}
		if (error) {
			log.reportTimeError("Invalid cel file list, apt will not run");
			throw new IllegalArgumentException();
		}
	}

	public static void run(String aptExeDir, String aptLibDir, String cels, String outDir, String quantNormTarget, String analysisName, String markerPositions, int markerBuffer, int maxWritersOpen,boolean full, int numThreads) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "affyPipeline.log");
		String[] celFiles;
		if (Files.isDirectory(cels)) {
			celFiles = Files.list(cels, null, ".cel", false, false, true);
		} else {
			celFiles = HashVec.loadFileToStringArray(cels, false, new int[] { 0 }, true);
		}
		validateCelSelection(celFiles, log);
		log.reportTimeInfo("Found " + celFiles.length + " .cel files to process");
		if (markerPositions == null || !Files.exists(markerPositions)) {
			log.reportTimeError("Could not find marker position file " + markerPositions);
			throw new IllegalArgumentException();
		}
		if (quantNormTarget == null || !Files.exists(quantNormTarget)) {
			log.reportTimeError("A valid target sketch file is required, and available from http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz");
			throw new IllegalArgumentException();
		}
		if (full) {
			log.reportTimeInfo("Running with full affymetrix cdf");
		} else {
			log.reportTimeInfo("Running with default affymetrix cdf, use the \"-full\" command to use the full version");
		}

		AffyPipeline pipeline = new AffyPipeline(aptExeDir, aptLibDir, full, log);
		Probesets probeSets = pipeline.getAnalysisProbesetList(celFiles[0], outDir, analysisName, markerPositions);
		if (!probeSets.isFail()) {
			String celListFile = pipeline.generateCelList(celFiles, outDir, analysisName);
			GenotypeResult genotypeResult = pipeline.genotype(celListFile, probeSets.getSnpOnlyFile(), analysisName, outDir);
			if (!genotypeResult.isFailed()) {

				NormalizationResult normalizationResult = pipeline.normalize(celListFile, probeSets.getAllFile(), analysisName, outDir, quantNormTarget);

				if (!normalizationResult.isFail()) {
					String tmpDir = outDir + analysisName + "_TMP/";

					String outSnpSrc = tmpDir + "SNP_Src/";
					new File(outSnpSrc).mkdirs();

					AffySNP6Tables AS6T = new AffySNP6Tables(outSnpSrc, genotypeResult.getCallFile(), genotypeResult.getConfFile(), normalizationResult.getQuantNormFile(), maxWritersOpen, log);
					AS6T.parseSNPTables(markerBuffer);

					String outCNSrc = tmpDir + "CN_Src/";
					new File(outCNSrc).mkdirs();

					AffySNP6Tables AS6TCN = new AffySNP6Tables(outCNSrc, normalizationResult.getQuantNormFile(), maxWritersOpen, log);
					AS6TCN.parseCNTable(markerBuffer);
					log.reportTimeInfo("Generating Genvisis project in " + outDir);

					String projectFile = outDir + analysisName + ".properties";
					Files.write("PROJECT_DIRECTORY=" + outDir, projectFile);
					Project proj = new Project(projectFile, false);

					proj.PROJECT_DIRECTORY.setValue(outDir);
					proj.SOURCE_DIRECTORY.setValue(analysisName + "_00src");
					proj.PROJECT_NAME.setValue(analysisName);
					proj.XY_SCALE_FACTOR.setValue((double) 100);

					proj.ARRAY_TYPE.setValue(ARRAY.AFFY_GW6);

					proj.SOURCE_FILENAME_EXTENSION.setValue(".txt.gz");
					proj.SOURCE_FILE_DELIMITER.setValue(SOURCE_FILE_DELIMITERS.TAB);
					proj.ID_HEADER.setValue("[FILENAME_ROOT]");
					proj.LONG_FORMAT.setValue(false);
					MergeChp.combineChpFiles(tmpDir, numThreads, "", ".txt", proj.SOURCE_DIRECTORY.getValue(true, true), log);
					if (Files.exists(markerPositions)) {
						if (!proj.MARKER_POSITION_FILENAME.getValue().equals(markerPositions)) {
							Files.copyFileUsingFileChannels(markerPositions, proj.MARKER_POSITION_FILENAME.getValue(), log);
						}
						// proj.MARKER_POSITION_FILENAME.setValue(markerPositions);
						if (proj.getSourceFileHeaders(true) == null || proj.getSourceFileHeaders(true).size() == 0 || Files.exists(proj.PROJECT_DIRECTORY.getValue() + Project.HEADERS_FILENAME)) {
							new File(proj.PROJECT_DIRECTORY.getValue() + Project.HEADERS_FILENAME).delete();
						}
						// proj.setSourceFileHeaders(proj.getSourceFileHeaders(true));
						proj.saveProperties();
						proj = new Project(projectFile, false);
						SourceFileParser.createFiles(proj, numThreads);
						TransposeData.transposeData(proj, 2000000000, false);
						CentroidCompute.computeAndDumpCentroids(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), new CentroidBuilder(), numThreads, 2);
						Centroids.recompute(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), false, numThreads);
						TransposeData.transposeData(proj, 2000000000, false);
						SampleData.createMinimalSampleData(proj);
					} else {
						log.reportTimeWarning("Missing file " + markerPositions);
						log.reportTimeWarning("Please provide the marker position at the command line, or use the Genvisis gui to finish parsing your affy project");
					}
					proj.saveProperties();

				}
			}
		}

	}

	public static void main(String[] args) {
		String analysisName = "Genvisis_affy_pipeline";
		String cels = "~/Affy6/cels/";
		String targetSketch = "~/resources/hapmap.quant-norm.normalization-target.txt";
		String aptExeDir = "~/apt-1.18.0-x86_64-intel-linux/bin/";
		String aptLibDir = "~/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/";
		String outDir = "~/Affy6_results/";
		String markerPositions = "~/resources/hg18.affy6.markerPostions";
		int numThreads = 1;
		int markerBuffer = 100;
		int maxWritersOpen = 1000000;
		int numArgs = args.length;
		boolean full = false;

		String usage = "\n" +
				"affy.AffyPipeline requires 0-1 arguments\n" +
				"   (1) analysis name (i.e. analysisName=" + analysisName + " (default))\n" +
				"   (2) a directory or full path to a file containing .cel files for analyiss (i.e. cels=" + cels + " (default))\n" +
				"   (3) a target sketch file (such as hapmap.quant-norm.normalization-target.txt Available at ) (i.e. sketch=" + targetSketch + " (default))\n" +
				"   (4) directory with Affy Power Tools executables (should contain apt-probeset-genotype, etc. Available at http://www.affymetrix.com/) (i.e. aptExeDir=" + aptExeDir + " (default))\n" +
				"   (5) directory with Affy Power Tools library files (should contain GenomeWideSNP_6.cdf, etc. Available at http://www.affymetrix.com/) (i.e. aptLibDir=" + aptLibDir + " (default))\n" +
				"   (6) output directory  (i.e. outDir=" + outDir + " (default))\n" +
				"   (7) full path to a file of marker positions (i.e. markerPositions=" + markerPositions + " (default))\n" +
				"   (8) optional: number of threads (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + "=" + numThreads + " (default))\n" +
				"   (9) optional: number of markers to buffer when splitting files (i.e. markerBuffer=" + markerBuffer + " (default))\n" +
				"   (10) optional: maximum number of writers to open, if this is less than the sample size parsing will slow drastically (i.e. maxWritersOpen=" + maxWritersOpen + " (default))\n" +
				"   (11) optional: use the full affymetrix cdf, which contains more mitochondrial probesets (i.e. -full (not the default))\n" +

				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("analysisName=")) {
				analysisName = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("markerPositions=")) {
				markerPositions = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("cels=")) {
				cels = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("sketch=")) {
				targetSketch = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("aptExeDir=")) {
				aptExeDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("aptLibDir=")) {
				aptLibDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("maxWritersOpen=")) {
				maxWritersOpen = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("markerBuffer")) {
				markerBuffer = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-full")) {
				full = true;
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
			run(aptExeDir, aptLibDir, cels, outDir, targetSketch, analysisName, markerPositions, markerBuffer, maxWritersOpen, full, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
