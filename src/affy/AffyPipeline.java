package affy;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

/**
 * @author lane0212
 *
 *         Maybe the last affy raw data processor. Takes us from .CEL -> to Genvisis. Requires affy power tools and related library files
 *
 */
public class AffyPipeline {

	private static final String AFFY_CEL_LIST_HEADER = "cel_files";
	private static final String AFFY_PROBELIST_HEADER = "probeset_id";

	private String aptExeDir;// holds "apt-geno-qc", "apt-probeset-genotype", "apt-probeset-summarize" etc
	private String aptLibDir;
	private Logger log;

	public AffyPipeline(String aptExeDir, String aptLibDir, Logger log) {
		super();
		this.aptExeDir = aptExeDir;
		this.aptLibDir = aptLibDir;
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
			if (!Files.exists(aptLibDir + AFFY_LIB_FILES.values()[i].getLibFile())) {
				log.reportTimeError(aptLibDir + AFFY_LIB_FILES.values()[i].getLibFile() + " did not exist");
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

		public String getLibFile() {
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
			return exe;
		}

	}

	private String getFullProbesetList(String celFile, String outDir, String analysisName, Logger log) {
		String smallCelList = outDir + "tmp" + ext.getTimestampForFilename();
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(smallCelList));
			writer.println(AFFY_CEL_LIST_HEADER);
			writer.println(celFile);
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to .cel list " + smallCelList);
			e.printStackTrace();
			System.exit(1);
		}

		ArrayList<String> psetCommand = new ArrayList<String>();
		psetCommand.add(aptExeDir);
		psetCommand.add("--cdf-file");
		psetCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile());
		psetCommand.add("--analysis");
		psetCommand.add("--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true");
		psetCommand.add("--out-dir");
		psetCommand.add("outDir");
		psetCommand.add("--cel-files");
		psetCommand.add("smallCelList");
		psetCommand.add("--set-analysis-name");
		psetCommand.add("analysisName");

		CmdLine.runCommandWithFileChecks(Array.toStringArray(psetCommand), "", new String[] { celFile }, null, true, false, false, log);

		// --out-dir " + outDir + " --cel-files " + smallCelList + " " + analysisName;
		// String psetCommand = aptExeDir + AFFY_ANALYSIS_TYPES[2] + " -a rma --cdf-file " + affyCDF + " -o " + affyResultsDir + " --cel-files " + smallCelList;
		log.report(ext.getTime() + " Info - running a command to extract probeset ids: " + psetCommand);
		// CmdLine.run(psetCommand, aptExeDir);
		// String toExtractFile = getMatchedFiles(affyResultsDir, log, AFFY_ANALYSIS_OUTPUTS_PREFIXES[2] + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0])[0];
		// log.report(ext.getTime() + " Info - extracting probesets from " + toExtractFile);
		// extractProbesets(affyResultsDir, markerFile, toExtractFile, log);
		// log.report(ext.getTime() + " Info - cleaning up files...");
		// deleteFile(smallCelList, log);
		// // delete summarize .summary
		// deleteFile(toExtractFile, log);
		// // delete summarize log
		// deleteFile(affyResultsDir + AFFY_ANALYSIS_TYPES[2] + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[4], log);
		// // delete summarize .report
		// deleteFile(affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[2] + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[1], log);
		return null;
	}

	public static void run(String aptExeDir, String aptLibDir, String celDir, String outDir, String quantNormTarget) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "affyPipeline.log");
		String[] celFiles = Files.list(celDir, null, ".cel", false, false, true);

		AffyPipeline pipeline = new AffyPipeline(aptExeDir, aptLibDir, log);
	}

	public static void main(String[] args) {
		String celDir = "/scratch.global/lanej/Affy6_1000g/cels/";
		String quantNormTarget = null;
		String aptExeDir = "/home/pankrat2/public/bin/affyPowerTools/apt-1.18.0-x86_64-intel-linux/bin/";
		String aptLibDir = "/home/pankrat2/public/bin/affyPowerTools/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/";
		String outDir = "/scratch.global/lanej/Affy6_1000g/";

		run(aptExeDir, aptLibDir, celDir, outDir, quantNormTarget);
	}

}
