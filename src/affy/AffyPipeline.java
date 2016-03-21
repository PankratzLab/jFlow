package affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
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

	private String generateCelList(String[] celFiles, String outDir, String analysisName) {
		String celListFile = outDir + analysisName + ".celList.txt";
		String[] toWrite = new String[] { AFFY_CEL_LIST_HEADER };
		toWrite = Array.concatAll(toWrite, celFiles);
		Files.writeList(toWrite, celListFile);
		return celListFile;
	}

	/**
	 * run a simple command to extract all probesets, return a formatted file to use for downstream analysis
	 */
	private String getFullProbesetList(String celFile, String outDir, String analysisName) {
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
			System.exit(1);
		}

		ArrayList<String> psetCommand = new ArrayList<String>();
		psetCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.GENERATE_PROBE_LIST.getExe());
		psetCommand.add("--cdf-file");
		psetCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile());
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
		CmdLine.runCommandWithFileChecks(Array.toStringArray(psetCommand), "", new String[] { celFile }, new String[] { probeResults }, true, false, false, log);

		log.reportTimeInfo("Parsing " + probeResults + " to obtain all probesIds");
		ArrayList<String> probesetIds = new ArrayList<String>(1800000);
		try {
			BufferedReader reader = Files.getAppropriateReader(probeResults);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				if (line[0].startsWith("CN_") || line[0].startsWith("SNP_") || line[0].startsWith("AFFX-SNP")) {
					String pId = line[0];
					if (!pId.endsWith("-B")) {
						probesetIds.add(pId);
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

		log.reportTimeInfo("Detected " + probesetIds.size() + " total probesets");
		String pIdFile = outDir + analysisName + ".probesetIds.txt";
		Files.writeList(Array.toStringArray(probesetIds), pIdFile);
		new File(smallCelList).delete();
		return pIdFile;
	}

	private void genotype(String celListFile, String pIDFile, String analysisName, String outDir) {
		
		String outCurrent = outDir+analysisName+"_Genotyping/";
		new File(outCurrent).mkdirs();
		ArrayList<String > genotypeCommand =new ArrayList<String>();
		genotypeCommand.add(aptExeDir+AFFY_ANALYSIS_TYPES.GENOTYPE.getExe());
		genotypeCommand.add("-c");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile());
		genotypeCommand.add("--table-output");
		genotypeCommand.add("true");
		genotypeCommand.add("-a");
		genotypeCommand.add("birdseed-v2");
		genotypeCommand.add("--set-gender-method");
		genotypeCommand.add("cn-probe-chrXY-ratio");
		genotypeCommand.add("--read-models-birdseed");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_BIRDSEED_MODELS.getLibFile());
		genotypeCommand.add("--special-snps");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_SPECIAL_SNPS.getLibFile());
		genotypeCommand.add("--chrX-probes");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CHRX.getLibFile());
		genotypeCommand.add("--chrY-probes");
		genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CHRY.getLibFile());
		genotypeCommand.add("--probeset-ids");
		genotypeCommand.add(pIDFile);
		genotypeCommand.add("--set-analysis-name");
		genotypeCommand.add(analysisName);
		genotypeCommand.add("-out-dir");
		genotypeCommand.add(outCurrent);
		genotypeCommand.add("--cel-files");
		genotypeCommand.add(celListFile);
		CmdLine.runCommandWithFileChecks(Array.toStringArray(genotypeCommand), "", null, null, true, false, false, log);

		

	}

	public static void run(String aptExeDir, String aptLibDir, String celDir, String outDir, String quantNormTarget, String analysisName) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "affyPipeline.log");
		String[] celFiles = Files.list(celDir, null, ".cel", false, false, true);
		log.reportTimeInfo("Found " + celFiles.length + " .cel files to process");
		AffyPipeline pipeline = new AffyPipeline(aptExeDir, aptLibDir, log);
		String piDFile = pipeline.getFullProbesetList(celFiles[0], outDir, analysisName);
		String celListFile = pipeline.generateCelList(celFiles, outDir, analysisName);

	}

	public static void main(String[] args) {
		String analysisName = "Genvisis_affy_pipeline";
		String celDir = "/scratch.global/lanej/Affy6_1000g/cels/";
		String quantNormTarget = null;
		String aptExeDir = "/home/pankrat2/public/bin/affyPowerTools/apt-1.18.0-x86_64-intel-linux/bin/";
		String aptLibDir = "/home/pankrat2/public/bin/affyPowerTools/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/";
		String outDir = "/scratch.global/lanej/Affy6_1000g/";

		run(aptExeDir, aptLibDir, celDir, outDir, quantNormTarget, analysisName);
	}

}
