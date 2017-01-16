package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.GATK_Genotyper.JointGATKGenotyper;
import org.genvisis.seq.analysis.SNPEFF.SnpEffResult;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFTumorNormalOps;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.googlecode.charts4j.collect.Maps;

public class GATK {
	public static final String GATK_LOCATION_COMMAND = "gatk=";
	public static final String TARGETED_REGION_COMMAND = "seqTarget=";
	public static final String DEFAULT_GATK = "/home/pankrat2/public/bin/GATK/";
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
	public static final String COMBINE_VARIANTS = "CombineVariants";
	public static final String VARIANT_RECALIBRATOR = "VariantRecalibrator";
	public static final String APPLY_RECALIBRATION = "ApplyRecalibration";
	public static final String VARIANT_ANNOTATOR = "VariantAnnotator";
	public static final String MUTECT2 = "MuTect2";

	public static final String PRINT_READS = "PrintReads";
	public static final String[] LOAD_R = {"module load R", "R"};
	public static final String RSCRIPT = "Rscript";
	public static final String NT = "-nt";
	public static final String NCT = "-nct";

	public static final String KNOWN_SITES_SNP_LOCATION_COMMAND = "knownSnps=";
	public static final String KNOWN_SITES_INDEL_LOCATION_COMMAND = "knownIndels=";
	public static final String ERC_MODE = "-ERC";
	public static final String GVCF_MODE = "GVCF";
	public static final String VARIANT_INDEX_TYPE = "-variant_index_type";
	public static final String LINEAR = "LINEAR";
	public static final String VARIANT_INDEX_PARAMETER = "-variant_index_parameter";
	public static final String VARIANT_INDEX_DEFAULT = "128000";
	public static final String MIN_N = "-minN";
	public static final String SET_KEY = "--setKey";
	public static final String FILTERED_ARE_UNCALLED = "--filteredAreUncalled";
	public static final String FILTERED_RECORDS_MERGE_TYPE = "--filteredrecordsmergetype";
	public static final String KEEP_IF_ANY_UNFILTERED = "KEEP_IF_ANY_UNFILTERED";
	public static final String DB_SNP = "--dbsnp";
	public static final String DB_SNP_FILE = "dbsnp";
	public static final String COSMIC = "--cosmic";
	public static final String PON = "-PON";

	public static final String GENOTYPING_MODE = "--genotyping_mode";
	public static final String DISCOVERY = "DISCOVERY";
	public static final String STAND_EMIT_CONF = "-stand_emit_conf";
	public static final String DEFAULT_STAND_EMIT_CONF = "10";
	public static final String STAND_CALL_CONF = "-stand_call_conf";
	public static final String DEFAULT_STAND_CALL_CONF = "30";

	public static final String BEFORE = "-before";
	public static final String AFTER = "-after";
	public static final String PLOTS = "-plots";
	public static final String BQSR = "-BQSR ";
	public static final String KNOWN = "-known";
	public static final String KNOWN_SITES = "-knownSites";
	public static final String VARIANT = "--variant";
	public static final String R = "-R";
	public static final String A = "-A";
	public static final String E = "-E";
	public static final String SNP_EFF = "SnpEff";
	public static final String SNP_EFF_FILE = "--snpEffFile";
	public static final String L = "-L";
	public static final String INTERVAL_PADDING = "-ip";
	public static final int DEFAULT_INTERVAL_PADDING = 100;
	public static final String V = "-V";

	public static final String I = "-I";
	public static final String I_TUMOR = "-I:tumor";
	public static final String I_NORMAL = "-I:normal";

	public static final String INPUT = "-input";

	public static final String MAX_GAUSSIANS = "--maxGaussians";
	public static final String DEFAULT_INDEL_MAX_GAUSSIANS = "4";
	// TODO: Determine if macGuassians 4 should apply to Indels (Question @ GATK)
	public static final String DEFAULT_TARGETED_MAX_GAUSSIANS = "4";
	public static final String TS_FILTER_LEVEL = "--ts_filter_level";
	public static final String DEFUALT_TS_FILTER_LEVEL_SNP = "99.0";
	public static final String DEFUALT_TS_FILTER_LEVEL_INDEL = "99.0";
	public static final String ARTIFACT_DETECTION_MODE = "--artifact_detection_mode";

	public static final String AN = "-an";
	// from https://www.broadinstitute.org/gatk/guide/article?id=1259
	// accessed 2016-08-29 (last updated 2016-06-28)
	public static final Set<String> ANS_BASE = ImmutableSet.of(	"QD", "FS", "SOR", "MQRankSum",
																															"ReadPosRankSum");
	public static final Set<String> ANS_SNP_ADDITIONS = ImmutableSet.of("MQ");
	public static final Set<String> ANS_INDEL_ADDITIONS = ImmutableSet.of();
	public static final Set<String> ANS_GENOME_ADDITIONS = ImmutableSet.of("DP");
	public static final Set<String> ANS_EXOME_ADDITIONS = ImmutableSet.of();
	public static final Set<String> ANS_INBREEDING_ADDITIONS = ImmutableSet.of("InbreedingCoeff");

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

	public static final String GZ = ".gz";
	public static final String GZ_INDEX = ".tbi";

	public static final String VCF = ".vcf";
	public static final String VCF_INDEX = ".idx";

	public static final String G = ".g";
	public static final String GVCF = G + VCF;

	public static final Set<String> VCF_EXTENSIONS = ImmutableSet.of(	GZ, GZ_INDEX, VCF, VCF_INDEX,
																																		GVCF, G);
	private static final String RESOURCE_ARG = "-resource:";
	
	public static final String TRANCHE = "-tranche";
	public static final String[] TRANCHES = {"100.0", "99.9", "99.0", "90.0"};

	private final String gatkLocation;
	private final String referenceGenomeFasta;
	private String[] knownSitesSnpFile;
	private String[] knownSitesIndelFile;
	private String dbSnpKnownSites;
	private String cosmicKnownSites;
	private final String javaLocation;
	private int memoryInMB;
	private boolean fail;
	private final boolean verbose;
	private final boolean overWriteExistingOutput;
	private final Logger log;
	private Map<RESOURCE, String> trainingResources;
	private String regionsFile;
	private String supportingSnps;
	private SEQ_TARGET seqTarget;

	private GATK(	String gATKLocation, String referenceGenomeFasta, String javaLocation,
								int memoryInMB, String dbSnpKnownSites, String regionsFile, SEQ_TARGET seqTarget,
								String cosmicKnownSites, boolean verbose, boolean overWriteExisting, Logger log) {
		this(	gATKLocation, referenceGenomeFasta, regionsFile, seqTarget, javaLocation, memoryInMB,
					verbose, overWriteExisting, log);
		this.dbSnpKnownSites = dbSnpKnownSites;
		this.cosmicKnownSites = cosmicKnownSites;
	}

	public GATK(String gATKLocation, String referenceGenomeFasta, String regionsFile, SEQ_TARGET seqTarget,
							int memoryInMB, boolean verbose, boolean overWriteExisting, Logger log) {
		this(gATKLocation, referenceGenomeFasta, regionsFile, seqTarget, null, memoryInMB, verbose, overWriteExisting, log);
	}

	public GATK(String GATKLocation, String referenceGenomeFasta, String regionsFile, SEQ_TARGET seqTarget,
							String javaLocation, int memoryInMB, boolean verbose, boolean overWriteExisting, Logger log) {
		this.gatkLocation = GATKLocation;
		this.referenceGenomeFasta = referenceGenomeFasta;
		this.regionsFile = regionsFile;
		this.seqTarget = seqTarget;
		this.javaLocation = javaLocation == null ? DEFAULT_JAVA : javaLocation;
		this.memoryInMB = memoryInMB;
		this.verbose = verbose;
		overWriteExistingOutput = overWriteExisting;
		this.log = log;
		fail = verifyGATKLocation();
	}

	public GATK(String gATKLocation, String referenceGenomeFasta, String regionsFile, SEQ_TARGET seqTarget,
							String javaLocation, int memoryInMB, String[] knownSitesSnpFile,
							String[] knownSitesIndelFile, boolean verbose, boolean overWriteExisting, Logger log) {
		this(	gATKLocation, referenceGenomeFasta, regionsFile, seqTarget, javaLocation, memoryInMB,
					verbose, overWriteExisting, log);
		this.knownSitesSnpFile = knownSitesSnpFile;
		this.knownSitesIndelFile = knownSitesIndelFile;
	}

	public static class Mutect extends GATK {
		public Mutect(String gATKLocation, String referenceGenomeFasta, int memoryInMB,
									String dbSnpKnownSites, String regionsFile, SEQ_TARGET seqTarget,
									String cosmicKnownSites, boolean verbose, boolean overWriteExisting, Logger log) {
			this(	gATKLocation, referenceGenomeFasta, null, memoryInMB, dbSnpKnownSites, regionsFile,
						seqTarget, cosmicKnownSites, verbose, overWriteExisting, log);
		}

		public Mutect(String gATKLocation, String referenceGenomeFasta, String javaLocation,
									int memoryInMB, String dbSnpKnownSites, String regionsFile,
									SEQ_TARGET seqTarget, String cosmicKnownSites, boolean verbose, boolean overWriteExisting, Logger log) {
			super(gATKLocation, referenceGenomeFasta, javaLocation, memoryInMB, dbSnpKnownSites,
						regionsFile, seqTarget, cosmicKnownSites, verbose, overWriteExisting, log);
		}
	}

	/**
	 * 
	 * @param filename
	 * @return
	 */
	public static final String getVcfIndex(String filename) {
		if (filename.endsWith(GZ)) {
			return filename + GZ_INDEX;
		} else {
			return filename + VCF_INDEX;
		}
	}

	/**
	 * Gets the root of a VCF filename, removing the directory info and any extensions found in
	 * {@link org.genvisis.seq.manage.VCFOps.VCF_EXTENSIONS}
	 * 
	 * @param filename a VCF filename
	 * @param trimDirectoryInfo true to also remove the directory info
	 * @return the root of the filename
	 * 
	 */
	public static final String getVcfRoot(String filename) {
		return getVcfRoot(filename, true);
	}

	/**
	 * Gets the root of a VCF filename, removing any extensions found in
	 * {@link org.genvisis.seq.manage.VCFOps.VCF_EXTENSIONS}
	 * 
	 * @param filename a VCF filename
	 * @param trimDirectoryInfo true to also remove the directory info
	 * @return the root of the filename
	 * 
	 */
	public static final String getVcfRoot(String filename, boolean trimDirectoryInfo) {
		if (trimDirectoryInfo) {
			filename = ext.removeDirectoryInfo(filename);
		}
		while (filename.lastIndexOf('.') > 0
						&& VCF_EXTENSIONS.contains(filename.substring(filename.lastIndexOf('.')))) {
			filename = filename.substring(0, filename.lastIndexOf('.'));
		}
		return filename;
	}

	public String getRegionsFile() {
		return regionsFile;
	}

	public SEQ_TARGET getSeqTarget() {
		return seqTarget;
	}

	public void setSupportingSnps(String supportingSnps) {
		this.supportingSnps = supportingSnps;
	}

	public boolean isFail() {
		return fail;
	}

	public String getGATKLocation() {
		return gatkLocation;
	}

	public String[] getKnownSitesSnpFile() {
		return knownSitesSnpFile;
	}

	public String getDbSnpKnownSites() {
		return dbSnpKnownSites;
	}

	public String getCosmicKnownSites() {
		return cosmicKnownSites;
	}

	public void setCosmicKnownSites(String cosmicKnownSites) {
		this.cosmicKnownSites = cosmicKnownSites;
	}

	public String[] getKnownSitesIndelFile() {
		return knownSitesIndelFile;
	}

	public Map<RESOURCE, String> getTrainingResources() {
		return trainingResources;
	}

	public void setTrainingResources(Map<RESOURCE, String> trainingResources) {
		this.trainingResources = ImmutableMap.copyOf(trainingResources);
	}

	public BaseRecalibration recalibrateABam(String baseId, String dedup_reads_bam, Logger altLog) {
		BaseRecalibration baseRecalibration = new BaseRecalibration(baseId, dedup_reads_bam,
																																altLog == null ? log : altLog);
		baseRecalibration.parseInput();
		boolean progress = determineBaseCovariation(dedup_reads_bam, baseRecalibration.getBqsr_before(),
																								baseRecalibration.getLog());
		if (progress) {
			progress = secondPassBaseCovariation(	dedup_reads_bam, baseRecalibration.getBqsr_before(),
																						baseRecalibration.getBqsr_post(), altLog);
			if (progress) {
				progress = analyzeBaseCovariation(baseRecalibration.getBqsr_before(),
																					baseRecalibration.getBqsr_post(),
																					baseRecalibration.getRecalibration_plots(), altLog);
				if (progress) {
					progress = applyBaseRecalibration(dedup_reads_bam, baseRecalibration.getBqsr_before(),
																						baseRecalibration.getRrd_bam(), altLog);
				}
			}
		}
		baseRecalibration.setFail(!progress);
		return baseRecalibration;
	}

	public IndelPrep realignABam(String baseId, String dedup_reads_bam, Logger altLog) {
		boolean progress = false;
		IndelPrep indelPrep = new IndelPrep(baseId, dedup_reads_bam, altLog == null ? log : altLog);
		indelPrep.parseInput();
		progress = determineTargetIndels(	indelPrep.getDedup_reads_bam(),
																			indelPrep.getTargetIntervalsList(), indelPrep.getLog());
		if (progress) {
			progress = realignTargetIndels(	indelPrep.getDedup_reads_bam(),
																			indelPrep.getTargetIntervalsList(),
																			indelPrep.getRealigned_dedup_reads_bam(), indelPrep.getLog());
			if (progress) {
				indelPrep.setAllThere(progress);
			}
		}
		indelPrep.setFail(!progress);
		return indelPrep;
	}

	public SingleSampleHaplotypeCaller haplotypeCallABam(	String baseId, String inputBam,
																												int numWithinSampleThreads, Logger altLog) {
		boolean progress = false;
		SingleSampleHaplotypeCaller haplotypeCaller = new SingleSampleHaplotypeCaller(inputBam, baseId,
																																									altLog);
		haplotypeCaller.parseInput();
		progress = singleSampleAllSitesCall(haplotypeCaller.getInputBam(),
																				haplotypeCaller.getOutputGVCF(), numWithinSampleThreads,
																				haplotypeCaller.getLog());
		haplotypeCaller.setFail(!progress);
		return haplotypeCaller;
	}

	public SnpEffResult annotateAVcfWithSnpEFF(SnpEffResult snpEffResult, boolean addDBSNP) {
		if (!snpEffResult.isFail()) {
			boolean progress = addSnpEffAnnotation(	snpEffResult.getInputVCF(),
																							snpEffResult.getOutputSnpEffVCF(),
																							snpEffResult.getOutputGatkSnpEffVCF(), addDBSNP,
																							snpEffResult.getLog());
			snpEffResult.setFail(!progress);
		} else {
			log.reportError("Error - could not annotate input vcf "	+ snpEffResult.getInputVCF()
											+ " with SNPEFF results");
		}
		return snpEffResult;
	}

	public String getReferenceGenomeFasta() {
		return referenceGenomeFasta;
	}

	public int getMemoryInMB() {
		return memoryInMB;
	}

	public void setMemoryInMB(int memoryInMB) {
		this.memoryInMB = memoryInMB;
	}

	private boolean verifyGATKLocation() {
		boolean verified = false;
		if (Files.exists(gatkLocation)) {
			verified = true;
		} else {
			fail = true;
		}
		return verified;
	}


	public boolean annotateWithAnotherVCF(String inputVcf, String annoVcf, String outVCF,
																				String[] annotations, String resourceName,
																				int numThreads) {
		String[] inputs = new String[] {inputVcf, annoVcf};
		String[] outputs = new String[] {outVCF};
		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add(VARIANT_ANNOTATOR);
		command.add(R);
		command.add(referenceGenomeFasta);
		command.add(V);
		command.add(inputVcf);
		command.addAll(intervalCommands());
		command.add(O);
		command.add(outVCF);
		command.add(RESOURCE_ARG + resourceName);
		command.add(annoVcf);
		for (String annotation : annotations) {
			command.add(E);
			command.add(resourceName + "." + annotation);
		}
		if (numThreads > 1) {
			command.add(NT);
			command.add(numThreads + "");
		}
		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs,
																						verbose, overWriteExistingOutput, false, log);
	}

	private boolean determineTargetIndels(String dedup_reads_bam, String output, Logger altLog) {
		boolean useKnownIndels = knownSitesIndelFile == null	? false
																													: Files.exists("", knownSitesIndelFile);
		if (!useKnownIndels && verbose) {
			if (knownSitesIndelFile == null) {
				log.report("Warning - known indel file(s) were not provided, skipping known indel realignment");
			} else {
				log.report("Warning - could not find all of the following known indel files:\n"
										+ Array.toStr(knownSitesIndelFile, "\n"));
			}
		}
		List<String> command = Lists.newArrayList(javaLocation,
																							PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																							gatkLocation + GENOME_ANALYSIS_TK, T,
																							REALIGNER_TARGET_CREATOR, R, referenceGenomeFasta, I,
																							dedup_reads_bam, O, output);

		if (useKnownIndels) {
			command = parseAndAddToCommand(command, KNOWN, knownSitesIndelFile);
		}
		command.addAll(intervalCommands());
		return CmdLine.runCommandWithFileChecks(command, "",
																						ImmutableList.of(referenceGenomeFasta, dedup_reads_bam),
																						ImmutableList.of(output), verbose,
																						overWriteExistingOutput, true,
																						altLog == null ? log : altLog);
	}

	private boolean realignTargetIndels(String dedup_reads_bam, String targetIntervalFile,
																			String output, Logger altLog) {
		String[] command = new String[] {	javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, INDEL_REALIGNER, R,
																			referenceGenomeFasta, I, dedup_reads_bam, TARGET_INTERVALS,
																			targetIntervalFile, O, output};
		return CmdLine.runCommandWithFileChecks(command, "",
																						new String[] {referenceGenomeFasta, dedup_reads_bam,
																													targetIntervalFile},
																						new String[] {output}, verbose, overWriteExistingOutput,
																						true, altLog == null ? log : altLog);
	}

	private boolean determineBaseCovariation(String dedup_reads_bam, String output, Logger altLog) {
		List<String> command = Lists.newArrayList(javaLocation,
																							PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																							gatkLocation + GENOME_ANALYSIS_TK, T,
																							BASE_RECALIBRATOR, R, referenceGenomeFasta, I,
																							dedup_reads_bam, O, output);
		if (checkKnowns()) {
			List<String> neccesaryInputFiles = ImmutableList.of(referenceGenomeFasta, dedup_reads_bam);
			neccesaryInputFiles = handleKnownSites(neccesaryInputFiles);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
			command.addAll(intervalCommands());
			return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles,
																							ImmutableList.of(output), verbose,
																							overWriteExistingOutput, true,
																							altLog == null ? log : altLog);
		} else {
			return false;
		}
	}

	private boolean secondPassBaseCovariation(String dedup_reads_bam, String bqsrFile, String output,
																						Logger altLog) {

		List<String> command = Lists.newArrayList(javaLocation,
																							PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																							gatkLocation + GENOME_ANALYSIS_TK, T,
																							BASE_RECALIBRATOR, R, referenceGenomeFasta, I,
																							dedup_reads_bam, O, output, BQSR, bqsrFile);
		if (checkKnowns()) {
			List<String> neccesaryInputFiles = ImmutableList.of(referenceGenomeFasta, dedup_reads_bam,
																													bqsrFile);
			neccesaryInputFiles = handleKnownSites(neccesaryInputFiles);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
			command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
			command.addAll(intervalCommands());
			return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles,
																							ImmutableList.of(output), verbose,
																							overWriteExistingOutput, true,
																							altLog == null ? log : altLog);
		} else {
			return false;
		}
	}

	private boolean analyzeBaseCovariation(	String before_recal_data, String after_recal_data,
																					String output, Logger altLog) {

		String[] command = new String[] {	javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, ANALYZE_COVARIATES, R,
																			referenceGenomeFasta, BEFORE, before_recal_data, AFTER,
																			after_recal_data, PLOTS, output};
		if (!CmdLine.runCommandWithFileChecks(command, "",
																					new String[] {referenceGenomeFasta, before_recal_data,
																												after_recal_data},
																					new String[] {output}, verbose, overWriteExistingOutput,
																					true, altLog == null ? log : altLog)) {
			altLog.reportError("Often this command fails due to not finding an R installation, R is needed to analyze the base covariation and is required for this pipeline");
			altLog.reportError("	 Please add the R and Rscript directory to your environment ${PATH}, or module load R if using a compute cluster");
			altLog.reportError("     Often the R library ggplot2 is unavailable, please install to generate plots using \"install.packages('ggplot2', dependencies = TRUE)\" on the R command line (you may need gplots, gsalib, and reshape as well)");
			return false;
		}
		return true;

	}

	private boolean applyBaseRecalibration(	String dedup_reads_bam, String bqsrFile, String output,
																					Logger altLog) {
		String[] command = new String[] {	javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, PRINT_READS, R,
																			referenceGenomeFasta, I, dedup_reads_bam, BQSR, bqsrFile, O,
																			output};
		return CmdLine.runCommandWithFileChecks(command, "",
																						new String[] {referenceGenomeFasta, dedup_reads_bam,
																													bqsrFile},
																						new String[] {output}, verbose, overWriteExistingOutput,
																						true, altLog == null ? log : altLog);
	}

	private boolean singleSampleAllSitesCall(	String bamFile, String output,
																						int numWithinSampleThreads, Logger altLog) {
		String dbSnpFile = null;
		List<String> input = Lists.newArrayList(referenceGenomeFasta, bamFile);
		if (knownSitesSnpFile != null && knownSitesSnpFile.length > 0) {
			for (String element : knownSitesSnpFile) {
				if (element.contains("dbsnp")) {
					dbSnpFile = element;
				}
			}
		}
		if (dbSnpFile == null) {
			log.reportError("Warning - did not detect a file containing \""	+ DB_SNP_FILE
											+ "\" in the file name, will not annotate variants");
		} else {
			if (verbose) {
				log.report(ext.getTime()	+ " Info - will annotate variants from " + bamFile
										+ " with db snp file " + dbSnpFile);
			}
			input.add(dbSnpFile);
		}

		List<String> command = Lists.newArrayList(javaLocation,
																							PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																							gatkLocation + GENOME_ANALYSIS_TK, T,
																							HAPLOTYPE_CALLER, R, referenceGenomeFasta, I, bamFile,
																							O, output, ERC_MODE, GVCF_MODE, GENOTYPING_MODE,
																							DISCOVERY, STAND_EMIT_CONF, DEFAULT_STAND_EMIT_CONF,
																							STAND_CALL_CONF, DEFAULT_STAND_CALL_CONF,
																							dbSnpFile == null ? "" : DB_SNP,
																							dbSnpFile == null ? "" : dbSnpFile, NCT,
																							Integer.toString(numWithinSampleThreads));
		command.addAll(intervalCommands());
		return CmdLine.runCommandWithFileChecks(command, "", input,
																						ImmutableList.of(output, getVcfIndex(output)), verbose,
																						overWriteExistingOutput, false,
																						altLog == null ? log : altLog);
	}

	/**
	 * @param vcfs
	 * @param outputVcf
	 * @param minN records with lt this number will not be combined
	 * @param log
	 * @return
	 */
	public boolean combinePonVcfs(String[] vcfs, String outputVcf, int minN, Logger log) {
		String[] input = new String[] {referenceGenomeFasta, regionsFile};
		input = Array.concatAll(input, vcfs);
		String[] outputs = new String[] {outputVcf, getVcfIndex(outputVcf)};

		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add(COMBINE_VARIANTS);
		command.add(R);
		command.add(referenceGenomeFasta);
		for (String vcf2 : vcfs) {
			command.add(V);
			command.add(vcf2);
		}
		command.add(MIN_N);
		command.add(minN + "");
		command.add(SET_KEY);
		command.add("\"null\"");
		command.add(FILTERED_ARE_UNCALLED);
		command.add(FILTERED_RECORDS_MERGE_TYPE);
		command.add(KEEP_IF_ANY_UNFILTERED);
		command.addAll(intervalCommands());
		command.add(O);
		command.add(outputVcf);
		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
																						verbose, overWriteExistingOutput, false, log);
	}

	public MutectTumorNormal callTumor(	String normalBam, String tumorBam, String outputVCF,
																			String pon, boolean rename, Logger log) {
		String[] input = new String[] {	referenceGenomeFasta, normalBam, tumorBam, dbSnpKnownSites,
																		regionsFile, cosmicKnownSites};
		if (pon != null) {
			input = Array.concatAll(input, new String[] {pon});
		} else {
			log.reportTimeWarning("Running tumor normal calling without PON");
		}
		String[] outputs = new String[] {outputVCF, getVcfIndex(outputVCF)};

		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add(MUTECT2);
		command.add(DB_SNP);
		command.add(dbSnpKnownSites);
		command.add(COSMIC);
		command.add(cosmicKnownSites);
		command.add(R);
		command.add(referenceGenomeFasta);
		command.add(I_NORMAL);
		command.add(normalBam);
		command.add(I_TUMOR);
		command.add(tumorBam);
		if (pon != null) {
			command.add(PON);
			command.add(pon);
		}
		command.addAll(intervalCommands());
		command.add(O);
		command.add(outputVCF);
		boolean progress = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input,
																												outputs, verbose, overWriteExistingOutput,
																												false, log);

		MutectTumorNormal mutectTumorNormal = new MutectTumorNormal(normalBam, tumorBam, outputVCF,
																																!progress);
		if (progress && rename && (!Files.exists(mutectTumorNormal.getReNamedFilteredVCF())
																|| !Files.exists(mutectTumorNormal.getReNamedOutputVCF()))) {
			log.reportTimeInfo("Re-naming samples in file " + outputVCF);
			String normalSamp = BamOps.getSampleName(normalBam, log);
			String tumorSamp = BamOps.getSampleName(tumorBam, log);
			VCFTumorNormalOps.renameTumorNormalVCF(	outputVCF, tumorSamp, normalSamp,
																							mutectTumorNormal.getReNamedOutputVCF(),
																							mutectTumorNormal.getReNamedFilteredVCF(), log);
		}
		return mutectTumorNormal;
	}

	public Mutect2Normal generateMutect2Normal(	String bamFile, String outputVcf,
																							int numWithinSampleThreads, Logger log) {
		boolean progress = mutect2NormalABam(bamFile, outputVcf, numWithinSampleThreads, log);
		return new Mutect2Normal(bamFile, outputVcf, !progress);
	}

	public static class GenotypeRefiner {
		private final String ped;
		private final String baseVCF;
		private final String cgpVCF;
		private final String filtVCF;
		private final String denovoVCF;
		private final String outputDir;
		private final Logger log;
		private boolean fail;

		public GenotypeRefiner(String ped, String baseVCF, String outputDir, Logger log) {
			super();
			this.ped = ped;
			this.baseVCF = baseVCF;
			this.outputDir = outputDir;
			cgpVCF = outputDir + VCFOps.getAppropriateRoot(baseVCF, true) + "postCGP.vcf.gz";
			filtVCF = outputDir + VCFOps.getAppropriateRoot(baseVCF, true) + "postCGP.Gfiltered.vcf.gz";
			denovoVCF = outputDir	+ VCFOps.getAppropriateRoot(baseVCF, true)
									+ "postCGP.Gfiltered.deNovos.vcf.gz";
			this.log = log;
			fail = false;

		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getPed() {
			return ped;
		}

		public String getBaseVCF() {
			return baseVCF;
		}

		public String getCgpVCF() {
			return cgpVCF;
		}

		public String getFiltVCF() {
			return filtVCF;
		}

		public String getDenovoVCF() {
			return denovoVCF;
		}

		public String getOutputDir() {
			return outputDir;
		}

		public Logger getLog() {
			return log;
		}

		public boolean isFail() {
			return fail;
		}

	}

	/**
	 * Running https://www.broadinstitute.org/gatk/guide/tagged?tag=denovo
	 */
	public GenotypeRefiner refineGenotypes(String vcf, String ped, String outputDir, Logger log) {
		new File(outputDir).mkdirs();
		GenotypeRefiner refiner = new GenotypeRefiner(ped, vcf, outputDir, log);
		boolean progress = derivePosteriorProbabilities(refiner.getBaseVCF(), refiner.getCgpVCF(), ped,
																										outputDir, log);

		if (progress) {
			progress = filterLowQualGenotpyes(refiner.getCgpVCF(), refiner.getFiltVCF(),
																				refiner.getOutputDir(), log);
			if (progress) {
				annotateDenovo(	refiner.getFiltVCF(), refiner.getDenovoVCF(), refiner.getPed(),
												refiner.getOutputDir(), log);
			}
		}
		refiner.setFail(!progress);
		return refiner;
	}

	private boolean annotateDenovo(	String inputVCF, String outputVCF, String ped, String outputDir,
																	Logger log) {

		new File(outputDir).mkdirs();
		String[] input = new String[] {	referenceGenomeFasta, inputVCF, ped, supportingSnps,
																		cosmicKnownSites};
		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add("VariantAnnotator");
		command.add(R);
		command.add(referenceGenomeFasta);
		command.add("-A");
		command.add("PossibleDeNovo");
		command.add("-ped");
		command.add(ped);
		command.add(V);
		command.add(inputVCF);
		command.add(O);
		command.add(outputVCF);
		String[] outputs = new String[] {outputVCF, outputVCF + GZ_INDEX};
		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
																						verbose, overWriteExistingOutput, false, log);
	}

	private boolean filterLowQualGenotpyes(	String inputVCF, String outputVCF, String outputDir,
																					Logger log) {

		new File(outputDir).mkdirs();
		String[] input = new String[] {	referenceGenomeFasta, inputVCF, supportingSnps,
																		cosmicKnownSites};
		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add("VariantFiltration");
		command.add(R);
		command.add(referenceGenomeFasta);
		command.add("-G_filter");
		command.add("\"GQ < 20.0\"");
		command.add("-G_filterName");
		command.add("lowGQ");
		command.add(V);
		command.add(inputVCF);
		command.add(O);
		command.add(outputVCF);
		String[] outputs = new String[] {outputVCF, outputVCF + GZ_INDEX};
		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
																						verbose, overWriteExistingOutput, false, log);
	}

	private boolean derivePosteriorProbabilities(	String inputVCF, String outputVCF, String ped,
																								String outputDir, Logger log) {

		new File(outputDir).mkdirs();
		String[] input = new String[] {referenceGenomeFasta, inputVCF, ped, supportingSnps};
		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add("CalculateGenotypePosteriors");
		command.add(R);
		command.add(referenceGenomeFasta);

		command.add("--supporting");
		command.add(supportingSnps);
		command.add("-ped");
		command.add(ped);
		command.add(V);
		command.add(inputVCF);
		command.add(O);
		command.add(outputVCF);
		String[] outputs = new String[] {outputVCF, outputVCF + GZ_INDEX};
		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
																						verbose, overWriteExistingOutput, false, log);
	}

	// TODO
	// https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php
	private boolean mutect2NormalABam(String bamFile, String outputVcf, int numWithinSampleThreads,
																		Logger log) {
		String[] input = new String[] {	referenceGenomeFasta, bamFile, dbSnpKnownSites, regionsFile,
																		cosmicKnownSites};
		ArrayList<String> command = new ArrayList<String>();
		command.add(javaLocation);
		command.add(JAR);
		command.add(gatkLocation + GENOME_ANALYSIS_TK);
		command.add(T);
		command.add(MUTECT2);
		command.add(R);
		command.add(referenceGenomeFasta);
		command.add(I_TUMOR);
		command.add(bamFile);
		// command.add(DB_SNP);
		// command.add(dbSnpKnownSites);
		// command.add(COSMIC);
		// command.add(cosmicKnownSites);
		command.add(ARTIFACT_DETECTION_MODE);
		command.addAll(intervalCommands());
		command.add(O);
		command.add(outputVcf);
		if (numWithinSampleThreads > 1) {
			command.add(NCT);
			command.add(Integer.toString(numWithinSampleThreads));
		}

		String[] outputs = new String[] {outputVcf, getVcfIndex(outputVcf)};
		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
																						verbose, overWriteExistingOutput, false, log);
	}

	private boolean addSnpEffAnnotation(String inputVCF, String snpEffVcf, String outputVCF,
																			boolean addDBSNP, Logger log) {
		boolean progress = true;
		List<String> inputFiles = Lists.newArrayList(inputVCF, snpEffVcf);
		List<String> outputFiles = Lists.newArrayList(getVcfIndex(outputVCF), outputVCF);
		List<String> command = Lists.newArrayList(javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, VARIANT_ANNOTATOR, R,
																			referenceGenomeFasta, A, SNP_EFF, VARIANT, inputVCF,
																			SNP_EFF_FILE, snpEffVcf, L, inputVCF, O, outputVCF);
		if (addDBSNP) {
			String dbSNPResource = trainingResources.get(RESOURCE.DBSNP);
			command.add(DB_SNP);
			command.add(dbSNPResource);
			inputFiles.add(dbSNPResource);
		}

		progress = CmdLine.runCommandWithFileChecks(command, "", inputFiles, outputFiles, verbose,
																								overWriteExistingOutput, true, log);
		return progress;
	}

	public JointGATKGenotyper recalibrateAVCF(final JointGATKGenotyper jGatkGenotyper, int numThreads,
																						Logger log) {
		boolean progress = !jGatkGenotyper.isFail();
		if (progress) {
			progress = buildSNPRecalibrationModel(jGatkGenotyper.getRawVCF(),
																						jGatkGenotyper.getRecalSNPFile(),
																						jGatkGenotyper.getTranchesSNPFile(),
																						jGatkGenotyper.getRscriptSNPFile(),
																						jGatkGenotyper.isIgnoreInbreeding(),
																						numThreads, log);
			if (progress) {
				progress = applySNPRecalibrationModel(jGatkGenotyper.getRawVCF(),
																							jGatkGenotyper.getRecalSNPFile(),
																							jGatkGenotyper.getTranchesSNPFile(),
																							jGatkGenotyper.getRecalSNP_VCF_File(), numThreads,
																							log);
				if (progress) {
					buildINDELRecalibrationModel(	jGatkGenotyper.getRecalSNP_VCF_File(),
																				jGatkGenotyper.getRecalINDELFile(),
																				jGatkGenotyper.getTranchesINDELFile(),
																				jGatkGenotyper.getRscriptINDELFile(),
																				jGatkGenotyper.isIgnoreInbreeding(),
																				numThreads, log);
					if (progress) {
						applyINDELRecalibrationModel(	jGatkGenotyper.getRecalSNP_VCF_File(),
																					jGatkGenotyper.getRecalINDELFile(),
																					jGatkGenotyper.getTranchesINDELFile(),
																					jGatkGenotyper.getRecalSNP_Indel_VCF_File(), numThreads,
																					log);
					}
				}
			}
		}
		jGatkGenotyper.setFail(!progress);
		return jGatkGenotyper;
	}

	private boolean buildSNPRecalibrationModel(	String inputVCF, String recalFile, String tranchesFile,
																							String rscriptFile, boolean ignoreInbreeding,
																							int numThreads, Logger altLog) {
		List<String> inputs = Lists.newArrayList(inputVCF);
		for (RESOURCE training : RESOURCE.SNP_TRAINING_RESOURCES) {
			inputs.add(trainingResources.get(training));
		}
		List<String> ouputs = ImmutableList.of(recalFile, tranchesFile, rscriptFile);
		List<String> command = Lists.newArrayList(javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, VARIANT_RECALIBRATOR, R,
																			referenceGenomeFasta, INPUT, inputVCF, MODE, SNP,
																			TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
																			R_SCRIPT_FILE, rscriptFile);
		if (seqTarget == SEQ_TARGET.TARGETED) {
			command.add(MAX_GAUSSIANS);
			command.add(DEFAULT_TARGETED_MAX_GAUSSIANS);
		}
		command.addAll(buildAns(true, getSeqTarget(), ignoreInbreeding, log));
		command.addAll(buildTranches());
		for (RESOURCE training : RESOURCE.SNP_TRAINING_RESOURCES) {
			command.add(training.getArgument());
			command.add(trainingResources.get(training));
		}
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
																						overWriteExistingOutput, false,
																						altLog == null ? log : altLog);
	}

	private boolean applySNPRecalibrationModel(	String inputVCF, String recalFile, String tranchesFile,
																							String output, int numThreads, Logger altLog) {
		String[] inputs = new String[] {inputVCF, recalFile, tranchesFile};
		String[] ouputs = new String[] {output};
		String[] command = new String[] {	javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, APPLY_RECALIBRATION, R,
																			referenceGenomeFasta, INPUT, inputVCF, MODE, SNP,
																			TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
																			TS_FILTER_LEVEL, DEFUALT_TS_FILTER_LEVEL_SNP, O, output};
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
																						overWriteExistingOutput, true,
																						altLog == null ? log : altLog);
	}

	private boolean buildINDELRecalibrationModel(	String inputVCF, String recalFile,
																								String tranchesFile, String rscriptFile,
																								boolean ignoreInbreeding, int numThreads,
																								Logger altLog) {
		List<String> inputs = Lists.newArrayList(inputVCF);
		for (RESOURCE training : RESOURCE.INDEL_TRAINING_RESOURCES) {
			inputs.add(trainingResources.get(training));
		}
		List<String> ouputs = Lists.newArrayList(recalFile, tranchesFile, rscriptFile);
		List<String> command = Lists.newArrayList(	javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
		      																			gatkLocation + GENOME_ANALYSIS_TK, T, VARIANT_RECALIBRATOR, R,
		      																			referenceGenomeFasta, INPUT, inputVCF, MODE, INDEL,
		      																			TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
		      																			R_SCRIPT_FILE, rscriptFile, MAX_GAUSSIANS,
		      																			DEFAULT_INDEL_MAX_GAUSSIANS);
		command.addAll(buildAns(false, getSeqTarget(), ignoreInbreeding, log));
		command.addAll(buildTranches());
		for (RESOURCE training : RESOURCE.INDEL_TRAINING_RESOURCES) {
			command.add(training.getArgument());
			command.add(trainingResources.get(training));
		}
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
																						overWriteExistingOutput, true,
																						altLog == null ? log : altLog);
	}

	private boolean applyINDELRecalibrationModel(	String inputVCF, String recalFile,
																								String tranchesFile, String output, int numThreads,
																								Logger altLog) {
		String[] inputs = new String[] {inputVCF, recalFile, tranchesFile};
		String[] ouputs = new String[] {output};
		// NO DQ!
		String[] command = new String[] {	javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, APPLY_RECALIBRATION, R,
																			referenceGenomeFasta, INPUT, inputVCF, MODE, INDEL,
																			TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
																			TS_FILTER_LEVEL, DEFUALT_TS_FILTER_LEVEL_INDEL, O, output};
		return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
																						overWriteExistingOutput, true,
																						altLog == null ? log : altLog);
	}

	/**
	 * @param inputGVCFs
	 * @param output
	 * @param numWithinSampleThreads
	 * @param altLog
	 * @return
	 */
	public boolean jointGenotypeGVCFs(String[] inputGVCFs, String output, int numWithinSampleThreads,
																		Logger altLog) {
		List<String> inputs = Lists.newArrayList();
		inputs.add(referenceGenomeFasta);
		for (String inputGVCF : inputGVCFs) {
			inputs.add(inputGVCF);
			inputs.add(getVcfIndex(inputGVCF));
		}
		List<String> inputGVCFArgs = Lists.newArrayList();
		for (String inputGVCF : inputGVCFs) {
			inputGVCFArgs.add(VARIANT);
			inputGVCFArgs.add(inputGVCF);
		}
		List<String> command = Lists.newArrayList(javaLocation,
																							PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																							gatkLocation + GENOME_ANALYSIS_TK, T, GENOTYPEGVCFS,
																							R, referenceGenomeFasta, O, output, NT,
																							Integer.toString(numWithinSampleThreads));
		command.addAll(inputGVCFArgs);
		return CmdLine.runCommandWithFileChecks(command, "", inputs,
																						ImmutableList.of(output, getVcfIndex(output)), verbose,
																						overWriteExistingOutput, true,
																						altLog == null ? log : altLog);
	}

	/**
	 * @param vcfs these vcfs will be merged to the output file
	 * @param output
	 * @param log
	 * @return true on successful merge
	 * @deprecated Use {@link #mergeVCFs(Collection,String,int,boolean,Logger)} instead
	 */
	public boolean mergeVCFs(	String[] vcfs, String output, int numthreads, boolean skipReporting,
														Logger log) {
															return mergeVCFs(Lists.newArrayList(vcfs), output, numthreads, skipReporting, log);
														}

	/**
	 * @param vcfs these vcfs will be merged to the output file
	 * @param output
	 * @param log
	 * @return true on successful merge
	 */
	public boolean mergeVCFs(	Collection<String> vcfs, String output, int numthreads, boolean skipReporting,
														Logger log) {
		List<String> command = Lists.newArrayList(javaLocation, PSF.Java.buildXmxString(getMemoryInMB()), JAR,
																			gatkLocation + GENOME_ANALYSIS_TK, T, COMBINE_VARIANTS, R,
																			referenceGenomeFasta, O, output, "-genotypeMergeOptions",
																			"UNIQUIFY");
		if (numthreads > 1) {
			command.add(NT);
			command.add(Integer.toString(numthreads));
		}
		if (regionsFile != null) {
			command.addAll(intervalCommands());
		}
		for (String vcf2 : vcfs) {
			command.add(VARIANT);
			command.add(vcf2);
		}
		return CmdLine.runCommandWithFileChecks(command, "", vcfs,
																						Lists.newArrayList(output, getVcfIndex(output)), verbose,
																						overWriteExistingOutput, skipReporting, log);
	}


	/**
	 * 
	 * @param SNP true for SNP recal, false for indel
	 * @param exome true for exome recal, false for genome
	 * @param ignoreInbreeding true to leave out InbreedingCoeff an, false to include
	 * @return
	 */
	private static List<String> buildAns(	boolean SNP, SEQ_TARGET seqTarget, boolean ignoreInbreeding,
																		Logger log) {
		Set<String> ansToUse = Sets.newHashSet(ANS_BASE);
		if (SNP) {
			ansToUse.addAll(ANS_SNP_ADDITIONS);
		} else {
			ansToUse.addAll(ANS_INDEL_ADDITIONS);
		}
		switch (seqTarget) {
			case TARGETED:
				log.reportTimeWarning("Using Exome recalibration annotationa for custom targeted sequencing");
				ansToUse.addAll(ANS_EXOME_ADDITIONS);
				break;
			case EXOME:
				ansToUse.addAll(ANS_EXOME_ADDITIONS);
				break;
			case GENOME:
				ansToUse.addAll(ANS_GENOME_ADDITIONS);
				break;
			default:
				log.reportError("Unrecognized Sequencing target, using Exome recalibration annotations");
				ansToUse.addAll(ANS_EXOME_ADDITIONS);
				break;
		}
		if (!ignoreInbreeding) {
			ansToUse.addAll(ANS_INBREEDING_ADDITIONS);
		}
		
		List<String> ans = Lists.newArrayList();
		for (String an : ansToUse) {
			ans.add(AN);
			ans.add(an);
		}
		return ans;
	}

	private static List<String> buildTranches() {
		List<String> tranches = Lists.newArrayList();
		for (String element : TRANCHES) {
			tranches.add(TRANCHE);
			tranches.add(element);
		}
		return tranches;
	}
	
	public static enum RESOURCE {
		HAPMAP("hapmap", false, true, true, 15.0),
		OMNI("omni", false, true, true, 12.0),
		G1K("1000G", false, true, false, 10.0),
		DBSNP("dbsnp", true, false, false, 2.0),
		MILLS("mills", false, true, true, 12.0);
		
		public static final String KNOWN_RESOURCE = "known=";
		public static final String TRAINING = "training=";
		public static final String TRUTH = "truth=";
		public static final String PRIOR = "prior=";
		
		public static final Set<RESOURCE> SNP_TRAINING_RESOURCES = ImmutableSet.of(HAPMAP, OMNI, G1K, DBSNP);
		public static final Set<RESOURCE> INDEL_TRAINING_RESOURCES = ImmutableSet.of(DBSNP, MILLS);
		
		private static final Map<String, RESOURCE> NAME_MAP;
		
		static {
			Map<String, RESOURCE> buildNameMap = Maps.newHashMap();
			for (RESOURCE resource : RESOURCE.values()) {
				buildNameMap.put(resource.name, resource);
			}
			NAME_MAP = ImmutableMap.copyOf(buildNameMap);
		}
		
		private String name;
		private boolean known;
		private boolean training;
		private boolean truth;
		private double prior;
		/**
		 * @param name
		 * @param known
		 * @param training
		 * @param truth
		 * @param prior
		 */
		private RESOURCE(String name, boolean known, boolean training, boolean truth, double prior) {
			this.name = name;
			this.known = known;
			this.training = training;
			this.truth = truth;
			this.prior = prior;
		}
		
		public String getName() {
			return name;
		}

		public boolean isKnown() {
			return known;
		}

		public boolean isTraining() {
			return training;
		}

		public boolean isTruth() {
			return truth;
		}

		public double getPrior() {
			return prior;
		}

		public String getArgument() {
			List<String> pieces = Lists.newArrayList();
			pieces.add(RESOURCE_ARG + name);
			pieces.add(KNOWN_RESOURCE + Boolean.toString(known));
			pieces.add(TRAINING + Boolean.toString(training));
			pieces.add(TRUTH + Boolean.toString(truth));
			pieces.add(PRIOR + Double.toString(prior));
			return Joiner.on(SPLIT).join(pieces);
		}
		
		public static Set<String> names() {
			return NAME_MAP.keySet();
		}
		
		public static RESOURCE getResourceByName(String name) {
			return NAME_MAP.get(name);
		}
	}

	public static enum SEQ_TARGET {
																	TARGETED, EXOME, GENOME;
	}

	public static class BaseRecalibration {
		private static final String RECAL_DATA = ".recal_data.table";
		private static final String RECAL = ".recal";

		private static final String POST = ".post";
		private static final String RECALIBRATION_PLOTS = ".recalibration_plots.pdf";

		private final String dedup_reads_bam;
		private String rrd_bam;
		private String bqsr_before;
		private String bqsr_post;
		private String recalibration_plots;
		private String baseID;
		private boolean allThere, fail;
		private final Logger log;

		public BaseRecalibration(String baseId, String dedup_reads_bam, Logger log) {
			super();
			this.baseID = baseId;
			this.dedup_reads_bam = dedup_reads_bam;
			allThere = false;
			fail = false;
			this.log = log;
		}

		public void parseInput() {
			rrd_bam = ext.addToRoot(dedup_reads_bam, RECAL);
			bqsr_before = ext.rootOf(dedup_reads_bam, false) + RECAL_DATA;
			bqsr_post = ext.addToRoot(bqsr_before, POST);
			recalibration_plots = ext.rootOf(dedup_reads_bam, false) + RECALIBRATION_PLOTS;
		}

		public String getDedup_reads_bam() {
			return dedup_reads_bam;
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

		public String getBaseID() {
			return baseID;
		}

		public void setBaseID(String baseId) {
			this.baseID = baseId;
		}

		public void setRrd_bam(String rrd_bam) {
			this.rrd_bam = rrd_bam;
		}
	}

	public static class SingleSampleHaplotypeCaller {
		private final String inputBam;
		private String outputGVCF;
		private String baseId;
		private boolean fail, allThere;
		private final Logger log;

		public SingleSampleHaplotypeCaller(String inputBam, String baseId, Logger log) {
			super();
			this.inputBam = inputBam;
			this.baseId = baseId;
			this.log = log;
		}

		public void parseInput() {
			outputGVCF = ext.rootOf(inputBam, false) + GVCF + GZ;
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
		private final String dedup_reads_bam;
		private String targetIntervalsList;
		private String realigned_dedup_reads_bam;
		private final String baseId;
		private boolean allThere, fail;
		private final Logger log;

		public IndelPrep(String baseId, String dedup_reads_bam, Logger log) {
			super();
			this.baseId = baseId;
			this.dedup_reads_bam = dedup_reads_bam;
			allThere = false;
			fail = false;
			this.log = log;
		}

		public void parseInput() {
			targetIntervalsList = ext.rootOf(dedup_reads_bam, false) + TARGET_INTERVALS;
			realigned_dedup_reads_bam = ext.addToRoot(dedup_reads_bam, REALIGNED);
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

	public static class Mutect2Normal {
		private final String normalBam;
		private final String outputVCF;
		private final boolean fail;

		public Mutect2Normal(String normalBam, String outputVCF, boolean fail) {
			super();
			this.normalBam = normalBam;
			this.outputVCF = outputVCF;
			this.fail = fail;
		}

		public boolean isFail() {
			return fail;
		}

		public String getOutputVCF() {
			return outputVCF;
		}

		public String getNormalBam() {
			return normalBam;
		}

	}

	public static class MutectTumorNormal {
		private final String normalBam;
		private final String tumorBam;
		private final String outputVCF;
		private String reNamedOutputVCF;
		private final String reNamedFilteredVCF;
		private final boolean fail;

		public MutectTumorNormal(String normalBam, String tumorBam, String outputVCF, boolean fail) {
			super();
			this.normalBam = normalBam;
			this.tumorBam = tumorBam;
			this.outputVCF = outputVCF;
			reNamedOutputVCF = VCFOps.getAppropriateRoot(outputVCF, false) + ".renamed.vcf.gz";
			reNamedFilteredVCF = VCFOps.getAppropriateRoot(outputVCF, false) + ".renamed.filtered.vcf.gz";
			this.fail = fail;
		}

		public String getReNamedOutputVCF() {
			return reNamedOutputVCF;
		}

		public void setReNamedOutputVCF(String reNamedOutputVCF) {
			this.reNamedOutputVCF = reNamedOutputVCF;
		}

		public String getReNamedFilteredVCF() {
			return reNamedFilteredVCF;
		}

		public String getNormalBam() {
			return normalBam;
		}

		public String getTumorBam() {
			return tumorBam;
		}

		public String getOutputVCF() {
			return outputVCF;
		}

		public boolean isFail() {
			return fail;
		}

	}

	// public static class Mutect2 {
	// private String normalBam;
	// private String tumorBam;
	// private String outputVCF;
	//
	// }

	/**
	 * 
	 * @param command existing List of commands
	 * @param commandToAdd command to add with each value
	 * @param values values to be added, each with commandToAdd
	 * @return the new List
	 */
	private static List<String> parseAndAddToCommand(	final List<String> command, String commandToAdd,
																										String... values) {
		List<String> combinedCommands = Lists.newArrayList(command);
		for (String value : values) {
			combinedCommands.add(commandToAdd);
			combinedCommands.add(value);
		}
		return combinedCommands;
	}
	
	private List<String> intervalCommands() {
		List<String> command = Lists.newArrayList();
		if (regionsFile == null) {
			log.reportTimeWarning("No regions file specified, necessary GATK calls will not be limited to"
														+ " targeted regions. This is only recommended for Whole Genome data.");
		} else {
			command.add(L);
			command.add(regionsFile);
			command.add(INTERVAL_PADDING);
			command.add(Integer.toString(DEFAULT_INTERVAL_PADDING));
		}
		return command;
	}

	private List<String> handleKnownSites(final List<String> neccesaryInputFiles) {
		return Array.concactAllToList(neccesaryInputFiles, knownSitesIndelFile, knownSitesSnpFile);
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
