package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.GATK_Genotyper.JointGATKGenotyper;
import org.genvisis.seq.analysis.SNPEFF.SnpEffResult;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFTumorNormalOps;

public class GATK {
  public static final String GATK_LOCATION_COMMAND = "gatk=";
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
  public static final String V = "-V";

  public static final String I = "-I";
  public static final String I_TUMOR = "-I:tumor";
  public static final String I_NORMAL = "-I:normal";

  public static final String INPUT = "-input";

  public static final String MAX_GAUSSIANS = "--maxGaussians";
  public static final String DEFAULT_MAX_GAUSSIANS = "4";
  public static final String TS_FILTER_LEVEL = "--ts_filter_level";
  public static final String DEFUALT_TS_FILTER_LEVEL_SNP = "99.5";
  public static final String DEFUALT_TS_FILTER_LEVEL_INDEL = "99.0";
  public static final String ARTIFACT_DETECTION_MODE = "--artifact_detection_mode";

  public static final String AN = "-an";
  // from https://www.broadinstitute.org/gatk/guide/article?id=1259
  // date = 12-17-14
  public static final String[] ANS_SNP = {"QD", "MQ", "MQRankSum", "ReadPosRankSum", "FS", "SOR",
                                          "InbreedingCoeff"};// NO DP for
                                                             // Exomes
  public static final String[] ANS_INDEL = {"QD", "FS", "SOR", "ReadPosRankSum", "MQRankSum",
                                            "InbreedingCoeff"};// NO DP for Exomes

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
  public static final String VCF_INDEX = ".idx";
  public static final String VCF_GZ_INDEX = ".tbi";

  public static final String GVCF = ".gvcf";
  public static final String RESOURCE = "-resource:";
  public static final String[] RESOURCES = {"hapmap", "omni", "1000G", "dbsnp"};
  public static final String KNOWN_RESOURCE = "known=";
  public static final String[] KNOWN_RESOURCES = {"false", "false", "false", "true"};
  public static final String TRAINING = "training=";
  public static final String[] TRAININGS = {"true", "true", "true", "false"};
  public static final String TRUTH = "truth=";
  public static final String[] TRUTHS = {"true", "true", "false", "false"};
  public static final String PRIOR = "prior=";
  public static final String[] PRIORS = {"15.0", "12.0", "10.0", "2.0"};

  public static final String TRANCHE = "-tranche";
  public static final String[] TRANCHES = {"100.0", "99.9", "99.5", "99.0", "90.0"};

  public static final String INDEL_RESOURCE_FULL_MILLS =
                                                       "-resource:mills,known=false,training=true,truth=true,prior=12.0";
  public static final String INDEL_RESOURCE_FULL_DBSNP =
                                                       "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0";

  private final String GATKLocation, referenceGenomeFasta;
  private String[] knownSitesSnpFile, knownSitesIndelFile;
  private String dbSnpKnownSites;
  private String cosmicKnownSites;
  private final String javaLocation;
  private boolean fail;
  private final boolean verbose;
  private final boolean overWriteExistingOutput;
  private final Logger log;
  private String hapMapTraining;
  private String omniTraining;
  private String thousandGTraining;
  private String dbSnpTraining;
  private String millsIndelTraining;
  private String regionsFile;
  private String supportingSnps;

  /**
   * Mutect constructor
   */
  private GATK(String gATKLocation, String referenceGenomeFasta, String javaLocation,
               String dbSnpKnownSites, String regionsFile, String cosmicKnownSites, boolean verbose,
               boolean overWriteExisting, Logger log) {
    this(gATKLocation, referenceGenomeFasta, javaLocation, verbose, overWriteExisting, log);
    this.dbSnpKnownSites = dbSnpKnownSites;
    this.cosmicKnownSites = cosmicKnownSites;
    this.regionsFile = regionsFile;
  }

  public GATK(String gATKLocation, String referenceGenomeFasta, boolean verbose,
              boolean overWriteExisting, Logger log) {
    this(gATKLocation, referenceGenomeFasta, null, verbose, overWriteExisting, log);
  }

  public GATK(String GATKLocation, String referenceGenomeFasta, String javaLocation,
              boolean verbose, boolean overWriteExisting, Logger log) {
    this.GATKLocation = GATKLocation;
    this.referenceGenomeFasta = referenceGenomeFasta;
    this.javaLocation = javaLocation == null ? DEFAULT_JAVA : javaLocation;
    this.verbose = verbose;
    overWriteExistingOutput = overWriteExisting;
    this.log = log;
    fail = verifyGATKLocation();
  }

  public GATK(String gATKLocation, String referenceGenomeFasta, String javaLocation,
              String[] knownSitesSnpFile, String[] knownSitesIndelFile, boolean verbose,
              boolean overWriteExisting, Logger log) {
    this(gATKLocation, referenceGenomeFasta, javaLocation, verbose, overWriteExisting, log);
    this.knownSitesSnpFile = knownSitesSnpFile;
    this.knownSitesIndelFile = knownSitesIndelFile;
  }

  public static class Mutect extends GATK {
    public Mutect(String gATKLocation, String referenceGenomeFasta, String dbSnpKnownSites,
                  String regionsFile, String cosmicKnownSites, boolean verbose,
                  boolean overWriteExisting, Logger log) {
      this(gATKLocation, referenceGenomeFasta, null, dbSnpKnownSites, regionsFile, cosmicKnownSites,
           verbose, overWriteExisting, log);
    }

    public Mutect(String gATKLocation, String referenceGenomeFasta, String javaLocation,
                  String dbSnpKnownSites, String regionsFile, String cosmicKnownSites,
                  boolean verbose, boolean overWriteExisting, Logger log) {
      super(gATKLocation, referenceGenomeFasta, javaLocation, dbSnpKnownSites, regionsFile,
            cosmicKnownSites, verbose, overWriteExisting, log);
    }
  }

  public String getRegionsFile() {
    return regionsFile;
  }

  public void setSupportingSnps(String supportingSnps) {
    this.supportingSnps = supportingSnps;
  }

  public void setRegionsFile(String regionsFile) {
    this.regionsFile = regionsFile;
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

  public String getDbSnpKnownSites() {
    return dbSnpKnownSites;
  }

  public String getCosmicKnownSites() {
    return cosmicKnownSites;
  }

  public void setCosmicKnownSites(String cosmicKnownSites) {
    this.cosmicKnownSites = cosmicKnownSites;
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

  public BaseRecalibration recalibrateABam(String baseId, String realigned_dedup_reads_bam,
                                           Logger altLog) {
    BaseRecalibration baseRecalibration = new BaseRecalibration(baseId, realigned_dedup_reads_bam,
                                                                (altLog == null ? log : altLog));
    baseRecalibration.parseInput();
    boolean progress = determineBaseCovariation(realigned_dedup_reads_bam,
                                                baseRecalibration.getBqsr_before(),
                                                baseRecalibration.getLog());
    if (progress) {
      progress = secondPassBaseCovariation(realigned_dedup_reads_bam,
                                           baseRecalibration.getBqsr_before(),
                                           baseRecalibration.getBqsr_post(), altLog);
      if (progress) {
        progress = analyzeBaseCovariation(baseRecalibration.getBqsr_before(),
                                          baseRecalibration.getBqsr_post(),
                                          baseRecalibration.getRecalibration_plots(), altLog);
        if (progress) {
          progress = applyBaseRecalibration(realigned_dedup_reads_bam,
                                            baseRecalibration.getBqsr_before(),
                                            baseRecalibration.getRrd_bam(), altLog);
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
    progress = determineTargetIndels(indelPrep.getDedup_reads_bam(),
                                     indelPrep.getTargetIntervalsList(), indelPrep.getLog());
    if (progress) {
      progress = realignTargetIndels(indelPrep.getDedup_reads_bam(),
                                     indelPrep.getTargetIntervalsList(),
                                     indelPrep.getRealigned_dedup_reads_bam(), indelPrep.getLog());
      if (progress) {
        indelPrep.setAllThere(progress);
      }
    }
    indelPrep.setFail(!progress);
    return indelPrep;
  }

  public SingleSampleHaplotypeCaller haplotypeCallABam(String baseId, String inputBam,
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
      boolean progress = addSnpEffAnnotation(snpEffResult.getInputVCF(),
                                             snpEffResult.getOutputSnpEffVCF(),
                                             snpEffResult.getOutputGatkSnpEffVCF(), addDBSNP,
                                             snpEffResult.getLog());
      snpEffResult.setFail(!progress);
    } else {
      log.reportError("Error - could not annotate input vcf " + snpEffResult.getInputVCF()
                      + " with SNPEFF results");
    }
    return snpEffResult;
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

  public boolean annotateWithAnotherVCF(String inputVcf, String annoVcf, String outVCF,
                                        String[] annotations, String resourceName, int numThreads) {
    String[] inputs = new String[] {inputVcf, annoVcf};
    String[] outputs = new String[] {outVCF};
    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
    command.add(T);
    command.add(VARIANT_ANNOTATOR);
    command.add(R);
    command.add(referenceGenomeFasta);
    command.add(V);
    command.add(inputVcf);
    command.add(O);
    command.add(outVCF);
    command.add(RESOURCE + resourceName);
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
    boolean useKnownIndels = knownSitesIndelFile == null ? false
                                                         : Files.exists("", knownSitesIndelFile);
    if (!useKnownIndels && verbose) {
      if (knownSitesIndelFile == null) {
        log.report("Warning - known indel file(s) were not provided, skipping known indel realignment");
      } else {
        log.report("Warning - could not find all of the following known indel files:\n"
                   + Array.toStr(knownSitesIndelFile, "\n"));
      }
    }
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     REALIGNER_TARGET_CREATOR, R, referenceGenomeFasta, I,
                                     dedup_reads_bam, O, output};

    if (useKnownIndels) {
      command = parseAndAddToCommand(command, KNOWN, knownSitesIndelFile);
    }
    return CmdLine.runCommandWithFileChecks(command, "",
                                            new String[] {referenceGenomeFasta, dedup_reads_bam},
                                            new String[] {output}, verbose, overWriteExistingOutput,
                                            true, (altLog == null ? log : altLog));
  }

  private boolean realignTargetIndels(String dedup_reads_bam, String targetIntervalFile,
                                      String output, Logger altLog) {
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     INDEL_REALIGNER, R, referenceGenomeFasta, I, dedup_reads_bam,
                                     TARGET_INTERVALS, targetIntervalFile, O, output};
    return CmdLine.runCommandWithFileChecks(command, "",
                                            new String[] {referenceGenomeFasta, dedup_reads_bam,
                                                          targetIntervalFile},
                                            new String[] {output}, verbose, overWriteExistingOutput,
                                            true, (altLog == null ? log : altLog));
  }

  private boolean determineBaseCovariation(String realigned_dedup_reads_bam, String output,
                                           Logger altLog) {
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     BASE_RECALIBRATOR, R, referenceGenomeFasta, I,
                                     realigned_dedup_reads_bam, O, output};
    if (checkKnowns()) {
      String[] neccesaryInputFiles = new String[] {referenceGenomeFasta, realigned_dedup_reads_bam};
      neccesaryInputFiles = handleKnownSites(neccesaryInputFiles, command);
      command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
      command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
      return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles,
                                              new String[] {output}, verbose,
                                              overWriteExistingOutput, true,
                                              (altLog == null ? log : altLog));
    } else {
      return false;
    }
  }

  private boolean secondPassBaseCovariation(String realigned_dedup_reads_bam, String bqsrFile,
                                            String output, Logger altLog) {
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     BASE_RECALIBRATOR, R, referenceGenomeFasta, I,
                                     realigned_dedup_reads_bam, O, output, BQSR, bqsrFile};
    if (checkKnowns()) {
      String[] neccesaryInputFiles = new String[] {referenceGenomeFasta, realigned_dedup_reads_bam,
                                                   bqsrFile};
      neccesaryInputFiles = handleKnownSites(neccesaryInputFiles, command);
      command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesIndelFile);
      command = parseAndAddToCommand(command, KNOWN_SITES, knownSitesSnpFile);
      return CmdLine.runCommandWithFileChecks(command, "", neccesaryInputFiles,
                                              new String[] {output}, verbose,
                                              overWriteExistingOutput, true,
                                              (altLog == null ? log : altLog));
    } else {
      return false;
    }
  }

  private boolean analyzeBaseCovariation(String before_recal_data, String after_recal_data,
                                         String output, Logger altLog) {

    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     ANALYZE_COVARIATES, R, referenceGenomeFasta, BEFORE,
                                     before_recal_data, AFTER, after_recal_data, PLOTS, output};
    if (!CmdLine.runCommandWithFileChecks(command, "",
                                          new String[] {referenceGenomeFasta, before_recal_data,
                                                        after_recal_data},
                                          new String[] {output}, verbose, overWriteExistingOutput,
                                          true, (altLog == null ? log : altLog))) {
      altLog.reportError("Often this command fails due to not finding an R installation, R is needed to analyze the base covariation and is required for this pipeline");
      altLog.reportError("	 Please add the R and Rscript directory to your environment ${PATH}, or module load R if using a compute cluster");
      altLog.reportError("     Often the R library ggplot2 is unavailable, please install to generate plots using \"install.packages('ggplot2', dependencies = TRUE)\" on the R command line (you may need gplots, gsalib, and reshape as well)");
      return false;
    }
    return true;

  }

  private boolean applyBaseRecalibration(String realigned_dedup_reads_bam, String bqsrFile,
                                         String output, Logger altLog) {
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     PRINT_READS, R, referenceGenomeFasta, I,
                                     realigned_dedup_reads_bam, BQSR, bqsrFile, O, output};
    return CmdLine.runCommandWithFileChecks(command, "",
                                            new String[] {referenceGenomeFasta,
                                                          realigned_dedup_reads_bam, bqsrFile},
                                            new String[] {output}, verbose, overWriteExistingOutput,
                                            true, (altLog == null ? log : altLog));
  }

  private boolean singleSampleAllSitesCall(String bamFile, String output,
                                           int numWithinSampleThreads, Logger altLog) {
    String dbSnpFile = null;
    String[] input = new String[] {referenceGenomeFasta, bamFile};
    if (knownSitesSnpFile != null && knownSitesSnpFile.length > 0) {
      for (String element : knownSitesSnpFile) {
        if (element.contains("dbsnp")) {
          dbSnpFile = element;
        }
      }
    }
    if (dbSnpFile == null) {
      log.reportError("Warning - did not detect a file containing \"" + DB_SNP_FILE
                      + "\" in the file name, will not annotate variants");
    } else {
      if (verbose) {
        log.report(ext.getTime() + " Info - will annotate variants from " + bamFile
                   + " with db snp file " + dbSnpFile);
      }
      input = Array.concatAll(input, new String[] {dbSnpFile});
    }

    String[] command =
                     new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                   HAPLOTYPE_CALLER, R, referenceGenomeFasta, I, bamFile, ERC_MODE,
                                   GVCF_MODE, VARIANT_INDEX_TYPE, LINEAR, VARIANT_INDEX_PARAMETER,
                                   VARIANT_INDEX_DEFAULT, dbSnpFile == null ? "" : DB_SNP,
                                   dbSnpFile == null ? "" : dbSnpFile, O, output, NCT,
                                   numWithinSampleThreads + ""};
    return CmdLine.runCommandWithFileChecks(command, "", input,
                                            new String[] {output, output + VCF_INDEX}, verbose,
                                            overWriteExistingOutput, true,
                                            (altLog == null ? log : altLog));
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
    String[] outputs = new String[] {outputVcf, outputVcf + VCF_INDEX};

    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
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
    command.add(L);
    command.add(regionsFile);
    command.add(O);
    command.add(outputVcf);
    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
                                            verbose, overWriteExistingOutput, false, log);
  }

  public MutectTumorNormal callTumor(String normalBam, String tumorBam, String outputVCF,
                                     String pon, boolean rename, Logger log) {
    String[] input = new String[] {referenceGenomeFasta, normalBam, tumorBam, dbSnpKnownSites,
                                   regionsFile, cosmicKnownSites};
    if (pon != null) {
      input = Array.concatAll(input, new String[] {pon});
    } else {
      log.reportTimeWarning("Running tumor normal calling without PON");
    }
    String[] outputs =
                     new String[] {outputVCF, outputVCF + (outputVCF.endsWith(".gz") ? VCF_GZ_INDEX
                                                                                     : VCF_INDEX)};

    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
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
    command.add(L);
    command.add(regionsFile);
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
      String normalSamp = BamOps.getSampleName(normalBam);
      String tumorSamp = BamOps.getSampleName(tumorBam);
      VCFTumorNormalOps.renameTumorNormalVCF(outputVCF, tumorSamp, normalSamp,
                                             mutectTumorNormal.getReNamedOutputVCF(),
                                             mutectTumorNormal.getReNamedFilteredVCF(), log);
    }
    return mutectTumorNormal;
  }

  public Mutect2Normal generateMutect2Normal(String bamFile, String outputVcf,
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
      denovoVCF = outputDir + VCFOps.getAppropriateRoot(baseVCF, true)
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
        annotateDenovo(refiner.getFiltVCF(), refiner.getDenovoVCF(), refiner.getPed(),
                       refiner.getOutputDir(), log);
      }
    }
    refiner.setFail(!progress);
    return refiner;
  }

  private boolean annotateDenovo(String inputVCF, String outputVCF, String ped, String outputDir,
                                 Logger log) {

    new File(outputDir).mkdirs();
    String[] input = new String[] {referenceGenomeFasta, inputVCF, ped, supportingSnps,
                                   cosmicKnownSites};
    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
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
    String[] outputs = new String[] {outputVCF, outputVCF + VCF_GZ_INDEX};
    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
                                            verbose, overWriteExistingOutput, false, log);
  }

  private boolean filterLowQualGenotpyes(String inputVCF, String outputVCF, String outputDir,
                                         Logger log) {

    new File(outputDir).mkdirs();
    String[] input =
                   new String[] {referenceGenomeFasta, inputVCF, supportingSnps, cosmicKnownSites};
    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
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
    String[] outputs = new String[] {outputVCF, outputVCF + VCF_GZ_INDEX};
    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
                                            verbose, overWriteExistingOutput, false, log);
  }

  private boolean derivePosteriorProbabilities(String inputVCF, String outputVCF, String ped,
                                               String outputDir, Logger log) {

    new File(outputDir).mkdirs();
    String[] input = new String[] {referenceGenomeFasta, inputVCF, ped, supportingSnps};
    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
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
    String[] outputs = new String[] {outputVCF, outputVCF + VCF_GZ_INDEX};
    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
                                            verbose, overWriteExistingOutput, false, log);
  }

  // TODO
  // https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php
  private boolean mutect2NormalABam(String bamFile, String outputVcf, int numWithinSampleThreads,
                                    Logger log) {
    String[] input = new String[] {referenceGenomeFasta, bamFile, dbSnpKnownSites, regionsFile,
                                   cosmicKnownSites};
    ArrayList<String> command = new ArrayList<String>();
    command.add(javaLocation);
    command.add(JAR);
    command.add(GATKLocation + GENOME_ANALYSIS_TK);
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
    command.add(L);
    command.add(regionsFile);
    command.add(O);
    command.add(outputVcf);
    if (numWithinSampleThreads > 1) {
      command.add(NCT);
      command.add(numWithinSampleThreads + "");
    }

    String[] outputs = new String[] {outputVcf, outputVcf + VCF_INDEX};
    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs,
                                            verbose, overWriteExistingOutput, false, log);
  }

  private boolean addSnpEffAnnotation(String inputVCF, String snpEffVcf, String outputVCF,
                                      boolean addDBSNP, Logger log) {
    boolean progress = true;
    String[] inputFiles = new String[] {inputVCF, snpEffVcf};
    String[] outputFiles = new String[] {inputVCF + VCF_INDEX, outputVCF};
    String[] command =
                     new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                   VARIANT_ANNOTATOR, R, referenceGenomeFasta, A, SNP_EFF, VARIANT,
                                   inputVCF, SNP_EFF_FILE, snpEffVcf, L, inputVCF, O, outputVCF};
    if (addDBSNP) {
      String[] dbSnp = new String[] {DB_SNP, getDbSnpTraining()};
      command = Array.concatAll(command, dbSnp);
      inputFiles = Array.concatAll(inputFiles, new String[] {getDbSnpTraining()});
    }

    progress = CmdLine.runCommandWithFileChecks(command, "", inputFiles, outputFiles, verbose,
                                                overWriteExistingOutput, true, log);
    return progress;
  }

  private String[] getCurrentResourceBundle() {
    String[] resourceArray = new String[RESOURCES.length * 2];
    String[] currentResourceBundle = new String[] {getHapMapTraining(), getOmniTraining(),
                                                   getThousandGTraining(), getDbSnpTraining()};
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

  public JointGATKGenotyper recalibrateAVCF(final JointGATKGenotyper jGatkGenotyper, int numThreads,
                                            Logger log) {
    boolean progress = !jGatkGenotyper.isFail();
    if (progress) {
      progress = buildSNPRecalibrationModel(jGatkGenotyper.getRawVCF(),
                                            jGatkGenotyper.getRecalSNPFile(),
                                            jGatkGenotyper.getTranchesSNPFile(),
                                            jGatkGenotyper.getRscriptSNPFile(), numThreads, log);
      if (progress) {
        progress =
                 applySNPRecalibrationModel(jGatkGenotyper.getRawVCF(),
                                            jGatkGenotyper.getRecalSNPFile(),
                                            jGatkGenotyper.getTranchesSNPFile(),
                                            jGatkGenotyper.getRecalSNP_VCF_File(), numThreads, log);
        if (progress) {
          buildINDELRecalibrationModel(jGatkGenotyper.getRecalSNP_VCF_File(),
                                       jGatkGenotyper.getRecalINDELFile(),
                                       jGatkGenotyper.getTranchesINDELFile(),
                                       jGatkGenotyper.getRscriptINDELFile(), numThreads, log);
          if (progress) {
            applyINDELRecalibrationModel(jGatkGenotyper.getRecalSNP_VCF_File(),
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

  private boolean buildSNPRecalibrationModel(String inputVCF, String recalFile, String tranchesFile,
                                             String rscriptFile, int numThreads, Logger altLog) {
    String[] inputs = new String[] {inputVCF, getHapMapTraining(), getOmniTraining(),
                                    getThousandGTraining(), getDbSnpTraining()};
    String[] ouputs = new String[] {recalFile, tranchesFile, rscriptFile};
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     VARIANT_RECALIBRATOR, R, referenceGenomeFasta, INPUT, inputVCF,
                                     MODE, SNP, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
                                     R_SCRIPT_FILE, rscriptFile};
    command = Array.concatAll(command, buildAns(true), getCurrentResourceBundle(), buildTranches());
    return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
                                            overWriteExistingOutput, false,
                                            (altLog == null ? log : altLog));
  }

  private boolean applySNPRecalibrationModel(String inputVCF, String recalFile, String tranchesFile,
                                             String output, int numThreads, Logger altLog) {
    String[] inputs = new String[] {inputVCF, recalFile, tranchesFile};
    String[] ouputs = new String[] {output};
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     APPLY_RECALIBRATION, R, referenceGenomeFasta, INPUT, inputVCF,
                                     MODE, SNP, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
                                     TS_FILTER_LEVEL, DEFUALT_TS_FILTER_LEVEL_SNP, O, output};
    return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
                                            overWriteExistingOutput, true,
                                            (altLog == null ? log : altLog));
  }

  private boolean buildINDELRecalibrationModel(String inputVCF, String recalFile,
                                               String tranchesFile, String rscriptFile,
                                               int numThreads, Logger altLog) {
    String[] inputs = new String[] {inputVCF, getMillsIndelTraining()};
    String[] ouputs = new String[] {recalFile, tranchesFile, rscriptFile};
    String[] command =
                     new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                   VARIANT_RECALIBRATOR, R, referenceGenomeFasta, INPUT, inputVCF,
                                   MODE, INDEL, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
                                   R_SCRIPT_FILE, rscriptFile, MAX_GAUSSIANS, DEFAULT_MAX_GAUSSIANS,
                                   INDEL_RESOURCE_FULL_MILLS, getMillsIndelTraining(),
                                   INDEL_RESOURCE_FULL_DBSNP, getDbSnpTraining()};
    command = Array.concatAll(command, buildAns(false), buildTranches());
    return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
                                            overWriteExistingOutput, true,
                                            (altLog == null ? log : altLog));
  }

  private boolean applyINDELRecalibrationModel(String inputVCF, String recalFile,
                                               String tranchesFile, String output, int numThreads,
                                               Logger altLog) {
    String[] inputs = new String[] {inputVCF, recalFile, tranchesFile};
    String[] ouputs = new String[] {output};
    // NO DQ!
    String[] command =
                     new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                   APPLY_RECALIBRATION, R, referenceGenomeFasta, INPUT, inputVCF,
                                   MODE, INDEL, TRANCHES_FILE, tranchesFile, RECAL_FILE, recalFile,
                                   TS_FILTER_LEVEL, DEFUALT_TS_FILTER_LEVEL_INDEL, O, output};
    return CmdLine.runCommandWithFileChecks(command, "", inputs, ouputs, verbose,
                                            overWriteExistingOutput, true,
                                            (altLog == null ? log : altLog));
  }

  public boolean jointGenotypeGVCFs(String[] inputGVCFs, String output, int numWithinSampleThreads,
                                    Logger altLog) {
    String[] inputs = new String[] {referenceGenomeFasta};
    inputs = Array.concatAll(inputs, inputGVCFs);
    String[] inputGVCFArgs = new String[inputGVCFs.length * 2];
    int index = 0;
    for (String inputGVCF : inputGVCFs) {
      inputGVCFArgs[index] = VARIANT;
      index++;
      inputGVCFArgs[index] = inputGVCF;
      index++;
    }
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     GENOTYPEGVCFS, R, referenceGenomeFasta, O, output, NT,
                                     numWithinSampleThreads + ""};
    command = Array.concatAll(command, inputGVCFArgs);
    return CmdLine.runCommandWithFileChecks(command, "", inputs, new String[] {output}, verbose,
                                            overWriteExistingOutput, true,
                                            (altLog == null ? log : altLog));
  }

  /**
   * @param vcfs these vcfs will be merged to the output file
   * @param output
   * @param log
   * @return true on successful merge
   */
  public boolean mergeVCFs(String[] vcfs, String output, int numthreads, boolean skipReporting,
                           Logger log) {
    String[] command = new String[] {javaLocation, JAR, GATKLocation + GENOME_ANALYSIS_TK, T,
                                     COMBINE_VARIANTS, R, referenceGenomeFasta, O, output,
                                     "-genotypeMergeOptions", "UNIQUIFY"};
    if (numthreads > 1) {
      command = Array.concatAll(command, new String[] {NT, numthreads + ""});
    }
    String[] inputArgVCF = new String[vcfs.length * 2];

    int index = 0;
    for (String vcf2 : vcfs) {
      inputArgVCF[index] = VARIANT;
      index++;
      inputArgVCF[index] = vcf2;
      index++;
    }
    command = Array.concatAll(command, inputArgVCF);
    return CmdLine.runCommandWithFileChecks(command, "", vcfs,
                                            new String[] {output,
                                                          output.endsWith(".gz") ? output
                                                                                   + VCF_GZ_INDEX
                                                                                 : output
                                                                                   + VCF_INDEX},
                                            verbose, overWriteExistingOutput, skipReporting, log);
  }

  private static String[] buildAns(boolean SNP) {
    String[] ansToUse = SNP ? ANS_SNP : ANS_INDEL;
    String[] ans = new String[ansToUse.length * 2];
    int index = 0;
    for (String element : ansToUse) {
      ans[index] = AN;
      index++;
      ans[index] = element;
      index++;
    }
    return ans;
  }

  private static String[] buildTranches() {
    String[] tranches = new String[TRANCHES.length * 2];
    int index = 0;
    for (String element : TRANCHES) {
      tranches[index] = TRANCHE;
      index++;
      tranches[index] = element;
      index++;
    }
    return tranches;
  }

  public static class BaseRecalibration {
    private static final String RECAL_DATA = ".recal_data.table";
    private static final String RECAL = ".recal";

    private static final String POST = ".post";
    private static final String RECALIBRATION_PLOTS = ".recalibration_plots.pdf";

    private final String realigned_dedup_reads_bam;
    private String rrd_bam;
    private String bqsr_before;
    private String bqsr_post;
    private String recalibration_plots;
    private String baseId;
    private String barcode;
    private String newBaseId;
    private boolean allThere, fail;
    private final Logger log;

    public BaseRecalibration(String baseId, String realigned_dedup_reads_bam, Logger log) {
      super();
      this.baseId = baseId;
      this.realigned_dedup_reads_bam = realigned_dedup_reads_bam;
      allThere = false;
      fail = false;
      this.log = log;
    }

    public void parseInput() {
      rrd_bam = ext.addToRoot(realigned_dedup_reads_bam, RECAL);
      bqsr_before = ext.rootOf(realigned_dedup_reads_bam, false) + RECAL_DATA;
      bqsr_post = ext.addToRoot(bqsr_before, POST);
      recalibration_plots = ext.rootOf(realigned_dedup_reads_bam, false) + RECALIBRATION_PLOTS;
      if (baseId.split(BWA_Analysis.FileNameParser.SPLIT).length != 3) {
        barcode = "";
        log.reportTimeWarning("The current baseId " + baseId + " did not have 3 "
                              + BWA_Analysis.FileNameParser.SPLIT
                              + " - delimited fields, assuming no barcodes are present in the ids");
      } else {
        barcode = baseId.split(BWA_Analysis.FileNameParser.SPLIT)[2];
      }
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

    public void setBaseId(String baseId) {
      this.baseId = baseId;
    }

    public String getNewBaseId() {
      return newBaseId;
    }

    public void setNewBaseId(String newBaseId) {
      this.newBaseId = newBaseId;
    }

    public String getBarcode() {
      return barcode;
    }

    public void setBarcode(String barcode) {
      this.barcode = barcode;
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
      outputGVCF = ext.rootOf(inputBam, false) + GVCF;
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

  private static String[] parseAndAddToCommand(String[] command, String commandToAdd,
                                               String[] values) {
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
