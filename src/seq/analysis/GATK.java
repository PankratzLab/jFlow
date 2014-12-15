package seq.analysis;

import seq.analysis.GATK_Genotyper.JointGATKGenotyper;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

public class GATK {
	public static final String GATK_LOCATION_COMMAND = "gatk=";
	public static final String KNOWN_SITES_SNP_LOCATION_COMMAND = "knownSnps=";
	public static final String KNOWN_SITES_INDEL_LOCATION_COMMAND = "knownIndels=";
	public static final String SPLIT = ",";

	public static final String DEFAULT_JAVA = "java";
	public static final String JAR = "-jar";
	public static final String GENOME_ANALYSIS_TK = "GenomeAnalysisTK.jar";

	public static final String T = "-T";
	public static final String REALIGNER_TARGET_CREATOR = "RealignerTargetCreator";
	public static final String INDEL_REALIGNER = "IndelRealigner";
	public static final String BASE_RECALIBRATOR = "BaseRecalibrator";
	public static final String ANALYZE_COVARIATES = "AnalyzeCovariates";
	public static final String HAPLOTYPE_CALLER = "HaplotypeCaller";
	public static final String GENOTYPEGVCFS = "GenotypeGVCFs";
	public static final String VARIANT_RECALIBRATOR = "VariantRecalibrator";
	public static final String APPLY_RECALIBRATION = "ApplyRecalibration";

	public static final String PRINT_READS = "PrintReads";
	public static final String[] LOAD_R = { "module load R", "R" };
	public static final String RSCRIPT = "Rscript";
	public static final String NT = "-nt";
	public static final String NCT = "-nct";

	public static final String ERC_MODE = "-ERC";
	public static final String GVCF_MODE = "GVCF";
	public static final String VARIANT_INDEX_TYPE = "-variant_index_type";
	public static final String LINEAR = "LINEAR";
	public static final String VARIANT_INDEX_PARAMETER = "-variant_index_parameter";
	public static final String VARIANT_INDEX_DEFAULT = "128000";

	public static final String DB_SNP = "--dbsnp";
	public static final String DB_SNP_FILE = "dbsnp";

	public static final String BEFORE = "-before";
	public static final String AFTER = "-after";
	public static final String PLOTS = "-plots";
	public static final String BQSR = "-BQSR ";
	public static final String KNOWN = "-known";
	public static final String KNOWN_SITES = "-knownSites";
	public static final String VARIANT = "--variant";
	public static final String R = "-R";
	public static final String I = "-I";
	public static final String INPUT = "-input";

	public static final String MAX_GAUSSIANS = "--maxGaussians";
	public static final String DEFAULT_MAX_GAUSSIANS = "4";
	public static final String TS_FILTER_LEVEL = "--ts_filter_level";
	public static final String DEFUALT_TS_FILTER_LEVEL = "99.0";
	public static final String AN = "-an";
	public static final String[] ANS = { "DP", "FS", "MQRankSum", "ReadPosRankSum" };
	public static final String AN_QD = "QD";
	public static final String MODE = "-mode";
	public static final String SNP = "SNP";
	public static final String INDEL = "INDEL";

	public static final String RECAL_FILE = "-recalFile";
	public static final String TRANCHES_FILE = "-tranchesFile";
	public static final String R_SCRIPT_FILE = "-rscriptFile";

	// public static final String L = "-L";
	// public static final String DEFAULT_L = "20";

	public static final String TARGET_INTERVALS = "-targetIntervals";
	public static final String O = "-o";

	public static final String VCF = ".vcf";
	public static final String GVCF = ".gvcf";
	public static final String RESOURCE = "-resource:";
	public static final String[] RESOURCES = { "hapmap", "omni", "1000G", "dbsnp" };
	public static final String KNOWN_RESOURCE = "known=";
	public static final String[] KNOWN_RESOURCES = { "false", "false", "false", "true" };
	public static final String TRAINING = "training=";
	public static final String[] TRAININGS = { "true", "true", "true", "false" };
	public static final String TRUTH = "truth=";
	public static final String[] TRUTHS = { "true", "true", "false", "false" };
	public static final String PRIOR = "prior=";
	public static final String[] PRIORS = { "15.0", "12.0", "10.0", "2.0" };

	public static final String TRANCHE = "-tranche";
	public static final String[] TRANCHES = { "100.0", "99.9", "99.0", "90.0" };

	public static final String INDEL_RESOURCE_FULL = "-resource:mills,known=true,training=true,truth=true,prior=12.0";

	private String GATKLocation, referenceGenomeFasta;
	private String[] knownSitesSnpFile, knownSitesIndelFile;
	private String javaLocation;
	private boolean fail, verbose, overWriteExisting;
	private Logger log;
	private String hapMapTraining;
	private String omniTraining;
	private String thousandGTraining;
	private String dbSnpTraining;
	private String millsIndelTraining;

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

	public String getMillsIndelTraining() {
		return millsIndelTraining;
	}

	public void setMillsIndelTraining(String millsIndelTraining) {
		this.millsIndelTraining = millsIndelTraining;
	}

	public String[] getKnownSitesIndelFile() {
		return knownSitesIndelFile;
	}

	public String getHapMapTraining() {
		return hapMapTraining;
	}

	public void setHapMapTraining(String hapMapTraining) {
		this.hapMapTraining = hapMapTraining;
	}

	public String getOmniTraining() {
		return omniTraining;
	}

	public void setOmniTraining(String omniTraining) {
		this.omniTraining = omniTraining;
	}

	public String getThousandGTraining() {
		return thousandGTraining;
	}

	public void setThousandGTraining(String thousandGTraining) {
		this.thousandGTraining = thousandGTraining;
	}

	public String getDbSnpTraining() {
		return dbSnpTraining;
	}

	public void setDbSnpTraining(String dbSnpTraining) {
		this.dbSnpTraining = dbSnpTraining;
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

	public SingleSampleHaplotypeCaller haplotypeCallABam(String baseId, String inputBam, int numWithinSampleThreads, Logger altLog) {
		boolean progress = false;
		SingleSampleHaplotypeCaller haplotypeCaller = new SingleSampleHaplotypeCaller(inputBam, baseId, altLog);
		haplotypeCaller.parseInput();
		progress = singleSampleAllSitesCall(haplotypeCaller.getInputBam(), haplotypeCaller.getOutputGVCF(), numWithinSampleThreads, haplotypeCaller.getLog());
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
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, dedup_reads_bam, targetIntervalFile }, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private boolean determineBaseCovariation(String realigned_dedup_reads_bam, String output, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, BASE_RECALIBRATOR, R, referenceGenomeFasta, I, realigned_dedup_reads_bam, O, output };
		if (checkKnowns()) {
			String[] neccesaryInputFiles = new String[] { referenceGenomeFasta, realigned_dedup_reads_bam };
			neccesaryInputFiles = handleKnownSites(neccesaryInputFiles, command);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
			return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
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
			return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
		} else {
			return false;
		}
	}

	private boolean analyzeBaseCovariation(String before_recal_data, String after_recal_data, String output, Logger altLog) {

		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, ANALYZE_COVARIATES, R, referenceGenomeFasta, BEFORE, before_recal_data, AFTER, after_recal_data, PLOTS, output };
		if (!CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, before_recal_data, after_recal_data }, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog))) {
			altLog.reportError("Often this command fails due to not finding an R installation, R is needed to analyze the base covariation and is required for this pipeline");
			altLog.reportError("	 Please add the R and Rscript directory to your environment ${PATH}, or module load R if using a compute cluster");
			altLog.reportError("     Often the R library ggplot2 is unavailable, please install to generate plots using \"install.packages('ggplot2', dependencies = TRUE)\" on the R command line (you may need gplots, gsalib, and reshape as well)");
			return false;
		}
		return true;

	}

	private boolean applyBaseRecalibration(String realigned_dedup_reads_bam, String bqsrFile, String output, Logger altLog) {
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, PRINT_READS, R, referenceGenomeFasta, I, realigned_dedup_reads_bam, BQSR, bqsrFile, O, output };
		return CmdLine.runCommandWithFileChecks(command, "", new String[] { referenceGenomeFasta, realigned_dedup_reads_bam, bqsrFile }, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private boolean singleSampleAllSitesCall(String bamFile, String output, int numWithinSampleThreads, Logger altLog) {
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

		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, HAPLOTYPE_CALLER, R, referenceGenomeFasta, I, bamFile, ERC_MODE, GVCF_MODE, VARIANT_INDEX_TYPE, LINEAR, VARIANT_INDEX_PARAMETER, VARIANT_INDEX_DEFAULT, dbSnpFile == null ? "" : DB_SNP, dbSnpFile == null ? "" : dbSnpFile, O, output, NCT, numWithinSampleThreads + "" };
		return CmdLine.runCommandWithFileChecks(command, "", input, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private String[] getCurrentResourceBundle() {
		String[] resourceArray = new String[RESOURCES.length * 2];
		String[] currentResourceBundle = new String[] { getHapMapTraining(), getOmniTraining(), getThousandGTraining(), getDbSnpTraining() };
		int index = 0;
		for (int i = 0; i < currentResourceBundle.length; i++) {
			resourceArray[index] = "";
			resourceArray[index] += RESOURCE + RESOURCES[i] + SPLIT;
			resourceArray[index] += KNOWN_RESOURCE + KNOWN_RESOURCES[i] + SPLIT;
			resourceArray[index] += TRAINING + TRAININGS[i] + SPLIT;
			resourceArray[index] += TRUTH + TRUTHS[i] + SPLIT;
			resourceArray[index] += PRIOR + PRIORS[i];
			index++;
			resourceArray[index] = currentResourceBundle[i];
			index++;
		}
		return resourceArray;
	}

	public JointGATKGenotyper recalibrateAVCF(JointGATKGenotyper jGatkGenotyper, int numThreads, Logger log) {
		boolean progress = !jGatkGenotyper.isFail();
		if (progress) {
			progress = buildSNPRecalibrationModel(jGatkGenotyper.getRawVCF(), jGatkGenotyper.getRecalSNPFile(), jGatkGenotyper.getTranchesSNPFile(), jGatkGenotyper.getRscriptSNPFile(), numThreads, log);
			if (progress) {
				progress = applySNPRecalibrationModel(jGatkGenotyper.getRawVCF(), jGatkGenotyper.getRecalSNPFile(), jGatkGenotyper.getTranchesSNPFile(), jGatkGenotyper.getRecalSNP_VCF_File(), numThreads, log);
				if (progress) {
					buildINDELRecalibrationModel(jGatkGenotyper.getRecalSNP_VCF_File(), jGatkGenotyper.getRecalINDELFile(), jGatkGenotyper.getTranchesINDELFile(), jGatkGenotyper.getRscriptINDELFile(), numThreads, log);
					if (progress) {
						applyINDELRecalibrationModel(jGatkGenotyper.getRecalSNP_VCF_File(), jGatkGenotyper.getRecalINDELFile(), jGatkGenotyper.getTranchesINDELFile(), jGatkGenotyper.getRecalSNP_Indel_VCF_File(), numThreads, log);
					}
				}
			}
		}
		jGatkGenotyper.setFail(!progress);
		return jGatkGenotyper;
	}

	private boolean buildSNPRecalibrationModel(String inputVCF, String recalFile, String tranchesFile, String rscriptFile, int numThreads, Logger altLog) {
		String[] inputs = new String[] { inputVCF, getHapMapTraining(), getOmniTraining(), getThousandGTraining(), getDbSnpTraining() };
		String[] ouputs = new String[] { recalFile, tranchesFile, rscriptFile };
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, VARIANT_RECALIBRATOR, R, referenceGenomeFasta, INPUT, inputVCF, MODE, SNP, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile, R_SCRIPT_FILE, rscriptFile, NCT, numThreads + "" };
		command = Array.concatAll(command, buildAns(true), getCurrentResourceBundle(), buildTranches());
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose, overWriteExisting, false, (altLog == null ? log : altLog));
	}

	private boolean applySNPRecalibrationModel(String inputVCF, String recalFile, String tranchesFile, String output, int numThreads, Logger altLog) {
		String[] inputs = new String[] { inputVCF, recalFile, tranchesFile };
		String[] ouputs = new String[] { output };
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, APPLY_RECALIBRATION, R, referenceGenomeFasta, INPUT, inputVCF, MODE, SNP, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile, TS_FILTER_LEVEL, DEFUALT_TS_FILTER_LEVEL, O, output, NCT, numThreads + "" };
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private boolean buildINDELRecalibrationModel(String inputVCF, String recalFile, String tranchesFile, String rscriptFile, int numThreads, Logger altLog) {
		String[] inputs = new String[] { inputVCF, getMillsIndelTraining() };
		String[] ouputs = new String[] { recalFile, tranchesFile, rscriptFile };
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, VARIANT_RECALIBRATOR, R, referenceGenomeFasta, INPUT, inputVCF, MODE, INDEL, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile, R_SCRIPT_FILE, rscriptFile, MAX_GAUSSIANS, DEFAULT_MAX_GAUSSIANS, INDEL_RESOURCE_FULL, getMillsIndelTraining(), NCT, numThreads + "" };
		command = Array.concatAll(command, buildAns(false), buildTranches());
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private boolean applyINDELRecalibrationModel(String inputVCF, String recalFile, String tranchesFile, String output, int numThreads, Logger altLog) {
		String[] inputs = new String[] { inputVCF, recalFile, tranchesFile };
		String[] ouputs = new String[] { output };
		// NO DQ!
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, APPLY_RECALIBRATION, R, referenceGenomeFasta, INPUT, inputVCF, MODE, INDEL, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile, TS_FILTER_LEVEL, DEFUALT_TS_FILTER_LEVEL, O, output, NCT, numThreads + "" };
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	public boolean jointGenotypeGVCFs(String[] inputGVCFs, String output, int numWithinSampleThreads, Logger altLog) {
		String[] inputs = new String[] { referenceGenomeFasta };
		inputs = Array.concatAll(inputs, inputGVCFs);
		String[] inputGVCFArgs = new String[inputGVCFs.length * 2];
		int index = 0;
		for (int i = 0; i < inputGVCFs.length; i++) {
			inputGVCFArgs[index] = VARIANT;
			index++;
			inputGVCFArgs[index] = inputGVCFs[i];
			index++;
		}
		String[] command = new String[] { javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T, GENOTYPEGVCFS, R, referenceGenomeFasta, O, output, NT, numWithinSampleThreads + "" };
		command = Array.concatAll(command, inputGVCFArgs);

		// System.out.println(Array.toStr(command));
		// System.exit(1);
		return CmdLine.runCommandWithFileChecks(command, "", inputs, new String[] { output }, verbose, overWriteExisting, true, (altLog == null ? log : altLog));
	}

	private static String[] buildAns(boolean SNP) {
		String[] ans = new String[ANS.length * 2];
		int index = 0;
		for (int i = 0; i < ANS.length; i++) {
			System.out.println(ans.length + "\t" + index);
			ans[index] = AN;
			index++;
			ans[index] = ANS[i];
			index++;
		}
		if (SNP) {
			ans = Array.concatAll(ans, new String[] { AN, AN_QD });
		}
		return ans;
	}

	private static String[] buildTranches() {
		String[] tranches = new String[TRANCHES.length * 2];
		int index = 0;
		for (int i = 0; i < TRANCHES.length; i++) {
			tranches[index] = TRANCHE;
			index++;
			tranches[index] = TRANCHES[i];
			index++;
		}
		return tranches;
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
