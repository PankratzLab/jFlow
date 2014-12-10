package seq.analysis;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

public class GATK {
	public static final String GATK_LOCATION_COMMAND = "gatk=";
	public static final String KNOWN_SITES_SNP_LOCATION_COMMAND = "knownSnps=";
	public static final String KNOWN_SITES_INDEL_LOCATION_COMMAND = "knownIndels=";
	public static final String KNOWN_SITES_SPLITTER = ",";

	public static final String DEFAULT_JAVA = "java";
	public static final String JAR = "-jar";
	public static final String GENOME_ANALYSIS_TK = "GenomeAnalysisTK.jar";

	public static final String T = "-T";
	public static final String REALIGNER_TARGET_CREATOR = "RealignerTargetCreator";
	public static final String INDEL_REALIGNER = "IndelRealigner";
	public static final String BASE_RECALIBRATOR = "BaseRecalibrator";
	public static final String ANALYZE_COVARIATES = "AnalyzeCovariates";
	public static final String HAPLOTYPE_CALLER = "HaplotypeCaller";

	public static final String PRINT_READS = "PrintReads";
	public static final String[] LOAD_R = { "module load R", "R" };
	public static final String RSCRIPT = "Rscript";
	public static final String NT = "-nt";
	public static final String NCT = "-nct";

	public static final String ERC_MODE = "--emitRefConfidence GVCF";
	public static final String VARIANT_INDEX_LINEAR = "--variant_index_type LINEAR";
	public static final String VARIANT_INDEX_PARAMETER = "--variant_index_parameter 128000";
	public static final String DB_SNP = "--dbsnp";
	public static final String DB_SNP_FILE = "dbsnp";

	public static final String BEFORE = "-before";
	public static final String AFTER = "-after";
	public static final String PLOTS = "-plots";
	public static final String BQSR = "-BQSR ";
	public static final String KNOWN = "-known";
	public static final String KNOWN_SITES = "-knownSites";

	public static final String R = "-R";
	public static final String I = "-I";
	// public static final String L = "-L";
	// public static final String DEFAULT_L = "20";

	public static final String TARGET_INTERVALS = "-targetIntervals";
	public static final String O = "-o";

	public static final String VCF = ".vcf";
	public static final String GVCF = ".gvcf";

	private String GATKLocation, referenceGenomeFasta;
	private String[] knownSitesSnpFile, knownSitesIndelFile;
	private String javaLocation;
	private boolean fail, verbose, overWriteExisting;
	private Logger log;

	public GATK(String gATKLocation, String referenceGenomeFasta, String javaLocation, String[] knownSitesSnpFile, String[] knownSitesIndelFile, boolean verbose, boolean overWriteExisting, Logger log) {
		super();
		this.GATKLocation = gATKLocation;
		this.referenceGenomeFasta = referenceGenomeFasta;
		this.javaLocation = (javaLocation == null ? DEFAULT_JAVA : javaLocation);
		this.knownSitesSnpFile = knownSitesSnpFile;
		this.knownSitesIndelFile = knownSitesIndelFile;
		this.fail = verifyGATKLocation();
		this.verbose = verbose;
		this.overWriteExisting = overWriteExisting;
		this.log = log;
	}

	public boolean isFail() {
		return fail;
	}

	public String getGATKLocation() {
		return GATKLocation;
	}

	public String[] getKnownSitesSnpFile() {
		return knownSitesSnpFile;
	}

	public String[] getKnownSitesIndelFile() {
		return knownSitesIndelFile;
	}

	public BaseRecalibration recalibrateABam(String baseId, String realigned_dedup_reads_bam, Logger altLog) {
		BaseRecalibration baseRecalibration = new BaseRecalibration(baseId, realigned_dedup_reads_bam, (altLog == null ? log : altLog));
		baseRecalibration.parseInput();
		boolean progress = determineBaseCovariation(realigned_dedup_reads_bam, baseRecalibration.getBqsr_before(), baseRecalibration.getLog());
		if (progress) {
			progress = secondPassBaseCovariation(realigned_dedup_reads_bam, baseRecalibration.getBqsr_before(), baseRecalibration.getBqsr_post(), altLog);
			if (progress) {
				progress = analyzeBaseCovariation(baseRecalibration.getBqsr_before(), baseRecalibration.getBqsr_post(), baseRecalibration.getRecalibration_plots(), altLog);
				if (progress) {
					progress = applyBaseRecalibration(realigned_dedup_reads_bam, baseRecalibration.getBqsr_before(), baseRecalibration.getRrd_bam(), altLog);
				}
			}
		}
		baseRecalibration.setFail(!progress);
		return baseRecalibration;
	}

	public IndelPrep realignABam(String baseId, String dedup_reads_bam, Logger altLog) {
		boolean progress = false;
		IndelPrep indelPrep = new IndelPrep(baseId, dedup_reads_bam, (altLog == null ? log : altLog));
		indelPrep.parseInput();
		progress = determineTargetIndels(indelPrep.getDedup_reads_bam(), indelPrep.getTargetIntervalsList(), indelPrep.getLog());
		if (progress) {
			progress = realignTargetIndels(indelPrep.getDedup_reads_bam(), indelPrep.getTargetIntervalsList(), indelPrep.getRealigned_dedup_reads_bam(), indelPrep.getLog());
			if (progress) {
				indelPrep.setAllThere(progress);
			}
		}
		indelPrep.setFail(!progress);
		return indelPrep;
	}

	public SingleSampleHaplotypeCaller haplotypeCallABam(String baseId, String inputBam, Logger altLog) {
		boolean progress = false;
		SingleSampleHaplotypeCaller haplotypeCaller = new SingleSampleHaplotypeCaller(inputBam, baseId, altLog);
		haplotypeCaller.parseInput();
		progress = singleSampleAllSitesCall(haplotypeCaller.getInputBam(), haplotypeCaller.getOutputGVCF(), haplotypeCaller.getLog());
		haplotypeCaller.setFail(progress);
		return haplotypeCaller;
	}

	public String getReferenceGenomeFasta() {
		return referenceGenomeFasta;
	}

	private boolean verifyGATKLocation() {
		boolean verified = false;
		if (Files.exists(GATKLocation)) {
			verified = true;
		} else {
			fail = true;
		}
		return verified;
	}

	private boolean determineTargetIndels(String dedup_reads_bam, String output, Logger altLog) {
		boolean useKnownIndels = knownSitesIndelFile == null ? false : Files.exists("", knownSitesIndelFile);
		if (!useKnownIndels && verbose) {
			if (knownSitesIndelFile == null) {
				log.report("Warning - known indel file(s) were not provided, skipping known indel realignment");
			} else {
				log.report("Warning - could not find all of the following known indel files:\n" + Array.toStr(knownSitesIndelFile, "\n"));
			}
		}
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, REALIGNER_TARGET_CREATOR, R, referenceGenomeFasta, I, dedup_reads_bam, O, output };

		if (useKnownIndels) {
			command = parseAndAddToCommand(command, KNOWN, knownSitesIndelFile);
		}
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, dedup_reads_bam }, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private boolean realignTargetIndels(String dedup_reads_bam, String targetIntervalFile, String output, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, INDEL_REALIGNER, R, referenceGenomeFasta, I, dedup_reads_bam, TARGET_INTERVALS, targetIntervalFile, O, output };
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, dedup_reads_bam, targetIntervalFile }, new String[] { output }, verbose, overWriteExisting, false, (altLog == null ? log : altLog));
	}

	private boolean determineBaseCovariation(String realigned_dedup_reads_bam, String output, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, BASE_RECALIBRATOR, R, referenceGenomeFasta, I, realigned_dedup_reads_bam, O, output };
		if (checkKnowns()) {
			String[] neccesaryInputFiles = new String[] { referenceGenomeFasta, realigned_dedup_reads_bam };
			neccesaryInputFiles = handleKnownSites(neccesaryInputFiles, command);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
			return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles, new String[] { output }, verbose, overWriteExisting, false, (altLog == null ? log : altLog));
		} else {
			return false;
		}
	}

	private boolean secondPassBaseCovariation(String realigned_dedup_reads_bam, String bqsrFile, String output, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, BASE_RECALIBRATOR, R, referenceGenomeFasta, I, realigned_dedup_reads_bam, O, output, BQSR, bqsrFile };
		if (checkKnowns()) {
			String[] neccesaryInputFiles = new String[] { referenceGenomeFasta, realigned_dedup_reads_bam, bqsrFile };
			neccesaryInputFiles = handleKnownSites(neccesaryInputFiles, command);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
			return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles, new String[] { output }, verbose, overWriteExisting, false, (altLog == null ? log : altLog));
		} else {
			return false;
		}
	}

	private boolean analyzeBaseCovariation(String before_recal_data, String after_recal_data, String output, Logger altLog) {

		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, ANALYZE_COVARIATES, R, referenceGenomeFasta, BEFORE, before_recal_data, AFTER, after_recal_data, PLOTS, output };
		if (!CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, before_recal_data, after_recal_data }, new String[] { output }, verbose, overWriteExisting, false, (altLog == null ? log : altLog))) {
			altLog.reportError("Often this command fails due to not finding an R installation, R is needed to analyze the base covariation and is required for this pipeline");
			altLog.reportError("	 Please add the R and Rscript directory to your environment ${PATH}, or module load R if using a compute cluster");
			altLog.reportError("     Often the R library ggplot2 is unavailable, please install to generate plots using \"install.packages('ggplot2', dependencies = TRUE)\" on the R command line (you may need gplots, gsalib, and reshape as well)");
			return false;
		}
		return true;

	}

	private boolean applyBaseRecalibration(String realigned_dedup_reads_bam, String bqsrFile, String output, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, PRINT_READS, R, referenceGenomeFasta, I, realigned_dedup_reads_bam, BQSR, bqsrFile, O, output };
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, realigned_dedup_reads_bam, bqsrFile }, new String[] { output }, verbose, overWriteExisting, false, (altLog == null ? log : altLog));
	}

	private boolean singleSampleAllSitesCall(String bamFile, String output, Logger altLog) {
		String dbSnpFile = null;
		String[] input = new String[] { referenceGenomeFasta, bamFile };
		if (knownSitesSnpFile != null && knownSitesSnpFile.length > 0) {
			for (int i = 0; i < knownSitesSnpFile.length; i++) {
				if (knownSitesSnpFile[i].contains("dbsnp")) {
					dbSnpFile = knownSitesSnpFile[i];
				}
			}
		}
		if (dbSnpFile == null) {
			log.reportError("Warning - did not detect a file containing \"" + DB_SNP_FILE + "\" in the file name, will not annotate variants");
		} else {
			if (verbose) {
				log.report(ext.getTime() + " Info - will annotate variants from " + bamFile + " with db snp file " + dbSnpFile);
			}
			input = Array.concatAll(input, new String[] { dbSnpFile });
		}

		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, HAPLOTYPE_CALLER, R, referenceGenomeFasta, I, bamFile, ERC_MODE, VARIANT_INDEX_LINEAR, VARIANT_INDEX_PARAMETER, dbSnpFile == null ? "" : DB_SNP + " " + dbSnpFile };
		return CmdLine.runCommandWithFileChecks(command, "", input, new String[] { output }, verbose, overWriteExisting, false, (altLog == null ? log : altLog));
	}

	public static class BaseRecalibration {
		private static final String RECAL_DATA = ".recal_data.table";
		private static final String RECAL = ".recal";

		private static final String POST = ".post";
		private static final String RECALIBRATION_PLOTS = ".recalibration_plots.pdf";

		private String realigned_dedup_reads_bam, rrd_bam, bqsr_before, bqsr_post, recalibration_plots, baseId;
		private boolean allThere, fail;
		private Logger log;

		public BaseRecalibration(String baseId, String realigned_dedup_reads_bam, Logger log) {
			super();
			this.baseId = baseId;
			this.realigned_dedup_reads_bam = realigned_dedup_reads_bam;
			this.allThere = false;
			this.fail = false;
			this.log = log;
		}

		public void parseInput() {
			this.rrd_bam = ext.addToRoot(realigned_dedup_reads_bam, RECAL);
			this.bqsr_before = ext.rootOf(realigned_dedup_reads_bam, false) + RECAL_DATA;
			this.bqsr_post = ext.addToRoot(bqsr_before, POST);
			this.recalibration_plots = ext.rootOf(realigned_dedup_reads_bam, false) + RECALIBRATION_PLOTS;
		}

		public String getRealigned_dedup_reads_bam() {
			return realigned_dedup_reads_bam;
		}

		public String getRrd_bam() {
			return rrd_bam;
		}

		public String getBqsr_before() {
			return bqsr_before;
		}

		public String getBqsr_post() {
			return bqsr_post;
		}

		public String getRecalibration_plots() {
			return recalibration_plots;
		}

		public boolean isAllThere() {
			return allThere;
		}

		public void setAllThere(boolean allThere) {
			this.allThere = allThere;
		}

		public Logger getLog() {
			return log;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getBaseId() {
			return baseId;
		}

	}

	public static class SingleSampleHaplotypeCaller {
		private String inputBam, outputGVCF, baseId;
		private boolean fail, allThere;
		private Logger log;

		public SingleSampleHaplotypeCaller(String inputBam, String baseId, Logger log) {
			super();
			this.inputBam = inputBam;
			this.baseId = baseId;
			this.log = log;
		}

		public void parseInput() {
			this.outputGVCF = ext.rootOf(inputBam, false) + GVCF;
		}

		public String getBaseId() {
			return baseId;
		}

		public void setBaseId(String baseId) {
			this.baseId = baseId;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public Logger getLog() {
			return log;
		}

		public boolean isAllThere() {
			return allThere;
		}

		public void setAllThere(boolean allThere) {
			this.allThere = allThere;
		}

		public String getInputBam() {
			return inputBam;
		}

		public String getOutputGVCF() {
			return outputGVCF;
		}

	}

	public static class IndelPrep {
		private static final String TARGET_INTERVALS = ".targetIntervals.list";
		private static final String REALIGNED = ".realigned";
		private String dedup_reads_bam, targetIntervalsList, realigned_dedup_reads_bam, baseId;
		private boolean allThere, fail;
		private Logger log;

		public IndelPrep(String baseId, String dedup_reads_bam, Logger log) {
			super();
			this.baseId = baseId;
			this.dedup_reads_bam = dedup_reads_bam;
			this.allThere = false;
			this.fail = false;
			this.log = log;
		}

		public void parseInput() {
			this.targetIntervalsList = ext.rootOf(dedup_reads_bam, false) + TARGET_INTERVALS;
			this.realigned_dedup_reads_bam = ext.addToRoot(dedup_reads_bam, REALIGNED);
		}

		public String getDedup_reads_bam() {
			return dedup_reads_bam;
		}

		public String getTargetIntervalsList() {
			return targetIntervalsList;
		}

		public String getBaseId() {
			return baseId;
		}

		public String getRealigned_dedup_reads_bam() {
			return realigned_dedup_reads_bam;
		}

		public boolean isAllThere() {
			return allThere;
		}

		public void setAllThere(boolean allThere) {
			this.allThere = allThere;
		}

		public Logger getLog() {
			return log;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

	}

	private static String[] parseAndAddToCommand(String[] command, String commandToAdd, String[] values) {
		String[] knowns;
		knowns = new String[values.length * 2];
		for (int i = 0; i < values.length; i++) {
			knowns[2 * i] = commandToAdd;
			knowns[2 * i + 1] = values[i];
		}
		command = Array.concatAll(command, knowns);
		return command;
	}

	private String[] handleKnownSites(String[] neccesaryInputFiles, String[] command) {
		neccesaryInputFiles = Array.concatAll(neccesaryInputFiles, knownSitesIndelFile);
		neccesaryInputFiles = Array.concatAll(neccesaryInputFiles, knownSitesSnpFile);
		return neccesaryInputFiles;
	}

	private boolean checkKnowns() {
		if (knownSitesSnpFile == null || knownSitesIndelFile == null) {
			log.reportError("Error - known indels and snps must be provided for " + BASE_RECALIBRATOR);
			return false;
		} else {
			return true;
		}
	}

}
