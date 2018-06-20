package org.genvisis.seq.analysis;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;
import org.genvisis.seq.analysis.ANNOVAR.AnnovarResults;
import org.genvisis.seq.analysis.GATK.RESOURCE;
import org.genvisis.seq.analysis.GATK.SEQ_TARGET;
import org.genvisis.seq.analysis.SNPEFF.SnpEffResult;
import org.genvisis.seq.manage.VCFOps;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class GATK_Genotyper {

  public static final String SPACE = " ";
  private final GATK gatk;
  private final SNPEFF snpeff;
  private final ANNOVAR annovar;
  private final SNPSIFT snpsift;
  private GATK.SingleSampleHaplotypeCaller[] siSampleHaplotypeCallers;
  private boolean fail, verbose;
  private final int numBetweenSampleThreads;
  private final int numWithinSampleThreads;
  private final Logger log;

  public GATK_Genotyper(GATK gatk, SNPEFF snpeff, SNPSIFT snpsift, ANNOVAR annovar,
                        int numBetweenSampleThreads, int numWithinSampleThreads, boolean verbose,
                        Logger log) {
    super();
    this.gatk = gatk;
    this.snpeff = snpeff;
    this.snpsift = snpsift;
    this.annovar = annovar;
    this.numBetweenSampleThreads = numBetweenSampleThreads;
    this.numWithinSampleThreads = numWithinSampleThreads;
    this.log = log;
  }

  public boolean isFail() {
    return fail;
  }

  public void setFail(boolean fail) {
    this.fail = fail;
  }

  public boolean isVerbose() {
    return verbose;
  }

  public void setVerbose(boolean verbose) {
    this.verbose = verbose;
  }

  public SNPEFF getSnpeff() {
    return snpeff;
  }

  public ANNOVAR getAnnovar() {
    return annovar;
  }

  public boolean runJointGenotyping(JointGATKGenotyper jGatkGenotyper) {
    boolean progress = !jGatkGenotyper.isFail();
    if (progress) {
      progress = gatk.jointGenotypeGVCFs(jGatkGenotyper.getInputGVCFs(),
                                         jGatkGenotyper.getOutputVCF(), numWithinSampleThreads,
                                         jGatkGenotyper.getLog());
      jGatkGenotyper.setFail(!progress);
    }
    return progress;
  }

  public boolean runRecalibration(final JointGATKGenotyper jGatkGenotyper) {
    boolean progress = !jGatkGenotyper.isFail();
    if (progress) {
      gatk.recalibrateAVCF(jGatkGenotyper, numWithinSampleThreads, log);
    }
    return progress;
  }

  public String annotateVCF(String inputVCF, ANNOTATION_BUILD build, MergeVCF mergeVCF,
                            ANNOVCF annoVCF) {
    if (!fail) {
      String in = inputVCF;
      String out = "";
      if (mergeVCF != null) {
        log.reportTimeInfo("Applying merge " + mergeVCF.getArg() + " prior to annotation");
        in = VCFOps.getAppropriateRoot(inputVCF, false) + ".merge_" + mergeVCF.getTag() + GATK.VCF
             + GATK.GZ;
        log.reportTimeInfo("Output merge: " + in);
        List<String> vcfsToMerge = Lists.newArrayList(mergeVCF.getVcfsToMergeWith());
        vcfsToMerge.add(0, inputVCF);
        gatk.mergeVCFs(vcfsToMerge, in, numWithinSampleThreads, false, log);
        out = in;
      }

      if (!annovar.isFail()) {
        AnnovarResults annovarResults = annovar.AnnovarAVCF(in, build, numWithinSampleThreads, log);
        in = annovarResults.getOutputVCF();
        out = annovarResults.getOutputVCF();
      }
      if (!snpeff.isFail()) {
        SnpEffResult snpEffResult = snpeff.annotateAVCF(in, build.getSnpEffBuild());
        gatk.annotateAVcfWithSnpEFF(snpEffResult, false);
        out = snpEffResult.getOutputGatkSnpEffVCF();
        if (!snpEffResult.isFail()) {
          if (annoVCF != null) {
            log.reportTimeInfo("Applying annotations from " + annoVCF.getArg());
            out = VCFOps.getAppropriateRoot(snpEffResult.getOutputGatkSnpEffVCF(), false) + ".anno_"
                  + annoVCF.getTag() + ".vcf";
            log.reportTimeInfo("Output anno: " + out);
            gatk.annotateWithAnotherVCF(snpEffResult.getOutputGatkSnpEffVCF(), annoVCF.getVcf(),
                                        out, annoVCF.getAnnos(), annoVCF.getTag(),
                                        numWithinSampleThreads);
            in = out;
          }
          log.reportTimeInfo("Since Annovar uses invalid char sequence for 1000g2014oct_* and 1000g2015aug_* (starts with number), replacing with g10002014oct_* or g10002015aug_*");
          log.reportTimeInfo("Note this is a hard-coded sed, so...");
          in = sed1000g(out, log);
          out = VCFOps.gzipAndIndex(in, log);
          // SnpSiftResult ssr = snpsift.annotateDbnsfp(in, log);
          // out = ssr.getOutputVCF();
        }
      }
      return out;
    }
    return null;
  }

  private String sed1000g(String in, Logger log) {
    String out = ext.addToRoot(in, ".sed1000g");
    String command = "cat " + in
                     + "|sed 's/1000g2014oct_/g10002014oct_/g'|sed 's/1000g2015aug_/g10002015aug_/g'>"
                     + out;
    String[] bat = CmdLine.prepareBatchForCommandLine(out + ".bat", true, log,
                                                      new String[] {command});
    if (CmdLine.runCommandWithFileChecks(bat, "", new String[] {in}, new String[] {out}, true,
                                         false, false, log)) {
      return out;
    } else {
      return null;
    }
  }

  public boolean determineTsTV(String inputVCF) {
    return !snpsift.tsTv(inputVCF, log).isFail();
  }

  public void batch(JointGATKGenotyper jointGATKGenotyper, String rootOutputDir, MergeVCF mergeVCF,
                    ANNOVCF annoVCF, boolean annotate, int memoryInMB, int wallTimeInHours,
                    String baseName) {
    String command = ArrayUtils.toStr(PSF.Load.getAllModules(), "\n");
    command += "\njava -Xmx" + memoryInMB + "m -jar ~/genvisisGATK3.6.jar "
               + this.getClass().getName() + SPACE;
    command += GATK_LanePrep.ROOT_INPUT_COMMAND + jointGATKGenotyper.getRootInputDir() + SPACE;
    command += GATK_LanePrep.ROOT_OUTPUT_COMMAND + rootOutputDir + SPACE;
    command += GATK_LanePrep.REFERENCE_GENOME_COMMAND + gatk.getReferenceGenomeFasta() + SPACE;
    command += GATK.GATK_LOCATION_COMMAND + gatk.getGATKLocation() + SPACE;
    command += NUM_THREADS + numWithinSampleThreads + SPACE;
    command += OUTPUT_COMMAND + jointGATKGenotyper.getOutput() + SPACE;
    command += GATK_LanePrep.LOG_FILE_COMMAND
               + ext.removeDirectoryInfo(jointGATKGenotyper.getLog().getFilename()) + SPACE;
    if (jointGATKGenotyper.getFileOfGVCFs() != null) {
      command += FILE_OF_GVCFS + jointGATKGenotyper.getFileOfGVCFs() + SPACE;
    }
    for (Map.Entry<GATK.RESOURCE, String> trainingResource : gatk.getTrainingResources()
                                                                 .entrySet()) {
      command += trainingResource.getKey().getName() + "=" + trainingResource.getValue() + SPACE;
    }
    if (annotate) {
      command += SNPEFF.SNP_EFF_COMMAND + snpeff.getSnpEffLocation() + SPACE;
      command += SNPSIFT.SNP_SIFT_LOCATION_COMMAND + snpsift.getSnpSiftLocation() + SPACE;
      command += ANNOVAR.ANNOVAR_COMMAND + annovar.getAnnovarLocation() + SPACE;
    } else {
      command += SNPEFF.SNP_EFF_NO_ANNO_COMMAND + SPACE;
    }
    if (jointGATKGenotyper.isIgnoreInbreeding()) {
      command += IGNORE_INBREEDING;
    }
    command += GATK.TARGETED_REGION_COMMAND + gatk.getSeqTarget().toString();
    if (mergeVCF != null) {
      command += MERGE_WITH + mergeVCF.getArg() + SPACE;
    }
    if (annoVCF != null) {
      command += EXTRA_VCF_ANNOTATIONS + annoVCF.getArg() + SPACE;
    }
    if (gatk.getRegionsFile() != null) {
      command += GATK_LanePrep.REGIONS_FILE_COMMAND + gatk.getRegionsFile();
    }
    Qsub.qsub("GATK_Genotype_" + baseName, command, memoryInMB, wallTimeInHours,
              numWithinSampleThreads);
  }

  public void runSingleSampleAllSites(String[] inputBams) {

    if (!isFail() && Files.checkAllFiles("", verbose, log, inputBams)) {
      if (inputBams != null) {
        siSampleHaplotypeCallers = new GATK.SingleSampleHaplotypeCaller[inputBams.length];
        int[] actualWithinSampleThreads = optimizeThreads(inputBams.length, numBetweenSampleThreads,
                                                          numWithinSampleThreads, log);
        WorkerHive<GATK.SingleSampleHaplotypeCaller> hive = new WorkerHive<>(numBetweenSampleThreads,
                                                                             10, log);
        WorkerSingleSampleAllSites[] workers = new WorkerSingleSampleAllSites[inputBams.length];
        for (int i = 0; i < inputBams.length; i++) {
          Logger altLog = new Logger(ext.rootOf(inputBams[i], false) + ".HC_ERC.log");
          workers[i] = new WorkerSingleSampleAllSites(gatk, inputBams[i], ext.rootOf(inputBams[i]),
                                                      actualWithinSampleThreads[i], altLog);
        }
        hive.addCallables(workers);
        hive.execute(true);
        siSampleHaplotypeCallers = hive.getResults()
                                       .toArray(new GATK.SingleSampleHaplotypeCaller[hive.getResults()
                                                                                         .size()]);
        for (int i = 0; i < workers.length; i++) {
          if (siSampleHaplotypeCallers[i].isFail()) {
            log.reportError("Failed single sample haplotype calling for "
                            + siSampleHaplotypeCallers[i].getInputBam());
            fail = true;
          }
        }
      } else {
        // TODO better check
      }
    }
  }

  private static int[] optimizeThreads(int numInputs, int numBetweenSampleThreads,
                                       int numWithinSampleThreads, Logger log) {
    int[] optimizedWithin = new int[numInputs];
    Arrays.fill(optimizedWithin, numWithinSampleThreads);
    if (numInputs < numBetweenSampleThreads) {
      int numExtra = numBetweenSampleThreads - numInputs;
      int index = 0;
      while (numExtra > 0) {
        numExtra--;
        optimizedWithin[index]++;
        index++;
        if (index == numInputs) {
          index = 0;
        }
      }
      log.report(ext.getTime()
                 + " Info - since we have extra between sample threads, we will allocate some to within sample(s) threads");
      log.report(ext.getTime() + " Info - allocated " + (numBetweenSampleThreads - numInputs)
                 + " extra thread(s) within sample(s) as follows: "
                 + ArrayUtils.toStr(optimizedWithin));
    }
    return optimizedWithin;
  }

  private static class WorkerSingleSampleAllSites implements Callable<GATK.SingleSampleHaplotypeCaller> {

    private final GATK GATK;
    private final String inputBam, baseId;
    private final int numWithinSampleThreads;
    private final Logger altLog;

    public WorkerSingleSampleAllSites(org.genvisis.seq.analysis.GATK gATK, String inputBam,
                                      String baseId, int numWithinSampleThreads, Logger altLog) {
      super();
      GATK = gATK;
      this.inputBam = inputBam;
      this.baseId = baseId;
      this.numWithinSampleThreads = numWithinSampleThreads;
      this.altLog = altLog;
    }

    @Override
    public GATK.SingleSampleHaplotypeCaller call() {
      return GATK.haplotypeCallABam(baseId, inputBam, numWithinSampleThreads, altLog);
    }
  }

  public static class JointGATKGenotyper {

    public static final String RECAL_EXT = ".recal";
    public static final String TRANCHES_EXT = ".tranches";
    public static final String RScript_EXT = ".R";

    private final String rootInputDir;
    private final String rootOutputDir;
    private final boolean ignoreInbreeding;
    private String output;
    private String rawVCF;
    private String fileOfGVCFs;
    private String recalSNP_VCF_File;
    private String recalSNP_Indel_VCF_File;
    private String recalSNPFile;
    private String tranchesSNPFile;
    private String rscriptSNPFile;
    private String recalINDELFile;
    private String tranchesINDELFile;
    private String rscriptINDELFile;
    private String[] inputGVCFs;
    private boolean fail;
    private final Logger log;

    /**
     * @param rootInputDir
     * @param rootOutputDir
     * @param output
     * @param log
     */
    public JointGATKGenotyper(String rootInputDir, String rootOutputDir, String output,
                              boolean ignoreInbreeding, Logger log) {
      super();
      this.rootInputDir = rootInputDir;
      this.rootOutputDir = rootOutputDir;
      this.output = output;
      this.ignoreInbreeding = ignoreInbreeding;
      fileOfGVCFs = null;
      this.log = log;
      fail = false;
    }

    public void setOutput(String output) {
      this.output = output;
    }

    public void init(String fileOfGVCFs) {
      this.fileOfGVCFs = fileOfGVCFs;
      initOutputs();
      if (fileOfGVCFs != null) {
        log.report(ext.getTime() + " Info - using GVCF files listed in the first column of"
                   + fileOfGVCFs);
        inputGVCFs = HashVec.loadFileToStringArray(fileOfGVCFs, false, new int[] {0}, true);
      } else if (rootInputDir == null) {
        log.reportError("Error - a file listing GVCF files was not provided and the root input directory was not provided, halting...");
        fail = true;
      } else {
        log.report(ext.getTime() + " Info - finding files with extension " + GATK.GVCF + GATK.GZ
                   + " in " + rootInputDir);
        inputGVCFs = Files.toFullPaths(Files.list(rootInputDir, GATK.GVCF + GATK.GZ), rootInputDir);

      }
      if (inputGVCFs == null || inputGVCFs.length < 1) {
        log.reportError("Error - could not find any GVCF files to joint genotype");
        fail = true;
      } else {
        log.report(ext.getTime() + " Info - using " + inputGVCFs.length
                   + " file(s) for joint Genotyping");
      }

    }

    public void initOutputs() {
      String currentRoot = rootOutputDir + GATK.getVcfRoot(output);
      rawVCF = currentRoot + GATK.VCF + GATK.GZ;

      recalSNPFile = currentRoot + "." + GATK.SNP + RECAL_EXT;
      tranchesSNPFile = currentRoot + "." + GATK.SNP + TRANCHES_EXT;
      rscriptSNPFile = currentRoot + "." + GATK.SNP + RScript_EXT;

      recalINDELFile = currentRoot + "." + GATK.INDEL + RECAL_EXT;
      tranchesINDELFile = currentRoot + "." + GATK.INDEL + TRANCHES_EXT;
      rscriptINDELFile = currentRoot + "." + GATK.INDEL + RScript_EXT;

      recalSNP_VCF_File = VCFOps.addToRoot(rawVCF, "." + GATK.SNP + RECAL_EXT);
      recalSNP_Indel_VCF_File = VCFOps.addToRoot(recalSNP_VCF_File, "." + GATK.INDEL + RECAL_EXT);
    }

    public String getOutput() {
      return output;
    }

    public String getOutputVCF() {
      return rawVCF;
    }

    public String getFileOfGVCFs() {
      return fileOfGVCFs;
    }

    public boolean isFail() {
      return fail;
    }

    public void setFail(boolean fail) {
      this.fail = fail;
    }

    public String getRootInputDir() {
      return rootInputDir;
    }

    public String[] getInputGVCFs() {
      return inputGVCFs;
    }

    public Logger getLog() {
      return log;
    }

    public String getRootOutputDir() {
      return rootOutputDir;
    }

    public static String getRecalExt() {
      return RECAL_EXT;
    }

    public static String getTranchesExt() {
      return TRANCHES_EXT;
    }

    public static String getRscriptExt() {
      return RScript_EXT;
    }

    public String getRawVCF() {
      return rawVCF;
    }

    public String getRecalSNP_VCF_File() {
      return recalSNP_VCF_File;
    }

    public String getRecalSNP_Indel_VCF_File() {
      return recalSNP_Indel_VCF_File;
    }

    public String getRecalSNPFile() {
      return recalSNPFile;
    }

    public String getTranchesSNPFile() {
      return tranchesSNPFile;
    }

    public String getRscriptSNPFile() {
      return rscriptSNPFile;
    }

    public String getRecalINDELFile() {
      return recalINDELFile;
    }

    public String getTranchesINDELFile() {
      return tranchesINDELFile;
    }

    public String getRscriptINDELFile() {
      return rscriptINDELFile;
    }

    public boolean isIgnoreInbreeding() {
      return ignoreInbreeding;
    }

  }

  public static void jointGenotype(String rootInputDir, String rootOutputDir, String output,
                                   String gATKLocation, String referenceGenomeFasta,
                                   String fileOfGVCFs, Map<RESOURCE, String> trainingResources,
                                   String snpEffLocation, String snpSiftLocation,
                                   String annovarLocation, ANNOTATION_BUILD annoBuild,
                                   String regionsFile, SEQ_TARGET seqTarget, MergeVCF mergeVCF,
                                   ANNOVCF annoVCF, boolean verbose, boolean overwriteExisting,
                                   boolean batch, boolean annotate, boolean ignoreInbreeding,
                                   boolean skipRecalibration, int numThreads, int memoryInMB,
                                   int wallTimeInHours, Logger log) {
    GATK gatk = new GATK(gATKLocation, referenceGenomeFasta, regionsFile, seqTarget, null,
                         memoryInMB, null, null, verbose, overwriteExisting, log);
    gatk.setTrainingResources(trainingResources);
    SNPEFF snpeff = new SNPEFF(snpEffLocation, verbose, overwriteExisting, log);
    SNPSIFT snpsift = new SNPSIFT(snpSiftLocation, verbose, overwriteExisting, log);
    ANNOVAR annovar = new ANNOVAR(annovarLocation, verbose, overwriteExisting, log);
    GATK_Genotyper genotyper = new GATK_Genotyper(gatk, snpeff, snpsift, annovar, 0, numThreads,
                                                  verbose, log);
    JointGATKGenotyper jGatkGenotyper = new JointGATKGenotyper(rootInputDir, rootOutputDir, output,
                                                               ignoreInbreeding, log);
    jGatkGenotyper.init(fileOfGVCFs);
    new File(rootOutputDir).mkdirs();
    if (batch) {
      genotyper.batch(jGatkGenotyper, rootOutputDir, mergeVCF, annoVCF, annotate, memoryInMB,
                      wallTimeInHours, output);
    } else {
      boolean progress = genotyper.runJointGenotyping(jGatkGenotyper);
      if (progress) {
        if (gatk.getRegionsFile() != null) {
          log.reportTimeInfo("Subsetting the vcf to regions defined in " + regionsFile + " with a "
                             + GATK.DEFAULT_INTERVAL_PADDING + " bp buffer");
          String subsetVcf = VCFOps.extractSegments(jGatkGenotyper.getRawVCF(), regionsFile,
                                                    GATK.DEFAULT_INTERVAL_PADDING, null,
                                                    rootOutputDir, false, true, 1, log);
          String newOutput = ext.removeDirectoryInfo(subsetVcf);
          jGatkGenotyper.setOutput(newOutput);
          jGatkGenotyper.initOutputs();// now starting with the subset vcf
        }
        if (skipRecalibration) {
          if (annotate) {
            genotyper.annotateVCF(jGatkGenotyper.output, annoBuild, mergeVCF, annoVCF);
          }
        } else {
          progress = genotyper.runRecalibration(jGatkGenotyper);
          if (progress) {
            if (progress && annotate) {
              genotyper.annotateVCF(jGatkGenotyper.getRecalSNP_Indel_VCF_File(), annoBuild,
                                    mergeVCF, annoVCF);
              genotyper.determineTsTV(jGatkGenotyper.getRecalSNP_Indel_VCF_File());
            }
          }
        }
      }
    }
  }

  public static String annotateOnlyWithDefualtLocations(String vcf, ANNOVCF annoVCF, int memoryInMB,
                                                        boolean verbose, boolean overwriteExisting,
                                                        Logger log) {
    return annotateOnly(vcf, AnnotationDefaults.GATK_LOC, AnnotationDefaults.REF, SEQ_TARGET.GENOME,
                        memoryInMB, null, AnnotationDefaults.SNP_EFF, null,
                        AnnotationDefaults.ANNOVAR, ANNOTATION_BUILD.HG19, annoVCF, verbose,
                        overwriteExisting, log);
  }

  public static String annotateOnly(String vcf, String gATKLocation, String referenceGenomeFasta,
                                    SEQ_TARGET seqTarget, int memoryInMB, String fileOfGVCFs,
                                    String snpEffLocation, String snpSiftLocation,
                                    String annovarLocation, ANNOTATION_BUILD annoBuild,
                                    ANNOVCF annoVCF, boolean verbose, boolean overwriteExisting,
                                    Logger log) {
    String snpSiftLoc = snpSiftLocation;
    if (snpSiftLoc == null || snpSiftLoc.equals(PSF.Ext.BLANK)) {
      snpSiftLoc = snpEffLocation;
    }
    GATK gatk = new GATK(gATKLocation, referenceGenomeFasta, null, seqTarget, null, memoryInMB,
                         null, null, verbose, overwriteExisting, log);
    SNPEFF snpeff = new SNPEFF(snpEffLocation, verbose, overwriteExisting, log);
    SNPSIFT snpsift = new SNPSIFT(snpSiftLoc, verbose, overwriteExisting, log);
    ANNOVAR annovar = new ANNOVAR(annovarLocation, verbose, overwriteExisting, log);
    GATK_Genotyper genotyper = new GATK_Genotyper(gatk, snpeff, snpsift, annovar, 0, 1, verbose,
                                                  log);
    return genotyper.annotateVCF(vcf, annoBuild, null, annoVCF);
  }

  public static String annotateOnly(String vcf, String gATKLocation, String referenceGenomeFasta,
                                    int memoryInMB, String snpEffLocation, String snpSiftLocation,
                                    String annovarLocation, ANNOTATION_BUILD annoBuild,
                                    boolean verbose, boolean overwriteExisting, Logger log) {
    return annotateOnly(vcf, gATKLocation, referenceGenomeFasta, SEQ_TARGET.GENOME, memoryInMB,
                        null, snpEffLocation, snpSiftLocation, annovarLocation, annoBuild, null,
                        verbose, overwriteExisting, log);
  }

  private static class MergeVCF {

    private final String tag;
    private final String[] vcfsToMergeWith;

    private MergeVCF(String tag, String[] vcfsToMergeWith) {
      super();
      this.tag = tag;
      this.vcfsToMergeWith = vcfsToMergeWith;
    }

    private String getArg() {
      return tag + ":" + ArrayUtils.toStr(vcfsToMergeWith, ",");

    }

    public String getTag() {
      return tag;
    }

    public String[] getVcfsToMergeWith() {
      return vcfsToMergeWith;
    }

    private static MergeVCF fromArg(String arg) {
      try {
        String[] tmp = ext.parseStringArg(arg).split(":");
        String tag = tmp[0];
        String[] vcfsToMergeWith = tmp[1].split(",");

        return new MergeVCF(tag, vcfsToMergeWith);
      } catch (Exception e) {
        System.err.println("Could not process merge arg " + arg + " format is tag:vcf1,vcf2");
        return null;
      }
    }

  }

  public static class ANNOVCF {

    private final String[] annos;
    private final String vcf;
    private final String tag;

    private ANNOVCF(String tag, String[] annos, String vcf) {
      super();
      this.tag = tag;
      this.annos = annos;
      this.vcf = vcf;
    }

    public String getArg() {
      return tag + ":" + vcf + ":" + ArrayUtils.toStr(annos, ",");

    }

    public String[] getAnnos() {
      return annos;
    }

    public String getVcf() {
      return vcf;
    }

    public String getTag() {
      return tag;
    }

    public static ANNOVCF fromArg(String arg) {
      try {
        String[] tmp = ext.parseStringArg(arg).split(":");
        String tag = tmp[0];
        String vcf = tmp[1];
        String[] annos = tmp[2].split(",");

        return new ANNOVCF(tag, annos, vcf);
      } catch (Exception e) {
        System.err.println("Could not process merge arg " + arg
                           + " format is tag:vcf.vcf:anno1,anno2");
        return null;
      }
    }

  }

  public static final String NUM_THREADS = "numThreads=";
  public static final String FILE_OF_GVCFS = "gvcfs=";
  public static final String OUTPUT_COMMAND = "output=";
  public static final String HAPMAP_TRAIN_COMMAND = "hapmap=";
  public static final String OMNI_TRAIN_COMMAND = "omni=";
  public static final String G1000_COMMAND = "1000G=";
  public static final String DBSNP_COMMMAND = "dbSNP=";
  public static final String MILLS = "mills=";
  public static final String ANNOTATE_VCF = "annotateThisVCFOnly=";
  public static final String IGNORE_INBREEDING = "-ignoreInbreeding";
  public static final String MERGE_WITH = "mergeWith=";
  public static final String EXTRA_VCF_ANNOTATIONS = "extraVCFAnno=";

  public enum ANNOTATION_BUILD {
    HG19("hg19", "hg19"), HG38("hg38", "GRCh38.86");

    private String annovarBuild;
    private String snpEffBuild;

    private ANNOTATION_BUILD(String annovarBuild, String snpEffBuild) {
      this.annovarBuild = annovarBuild;
      this.snpEffBuild = snpEffBuild;
    }

    /**
     * @return the annovarBuild
     */
    public String getAnnovarBuild() {
      return annovarBuild;
    }

    /**
     * @return the snpEffBuild
     */
    public String getSnpEffBuild() {
      return snpEffBuild;
    }

  }

  public static void main(String[] args) {

    int numArgs = args.length;
    String rootInputDir = null;
    String rootOutputDir = null;
    String referenceGenomeFasta = "";
    String gATKLocation = "";
    String output = "joint_genotypes";
    String fileOfGVCFs = null;
    boolean verbose = true;
    boolean batch = false;
    int memoryInMB = 23000;
    int wallTimeInHours = 24;
    int numThreads = 1;
    boolean overwriteExisting = false;
    Map<GATK.RESOURCE, String> trainingResources = Maps.newEnumMap(GATK.RESOURCE.class);
    String snpEffLocation = PSF.Ext.BLANK;
    String snpSiftLocation = PSF.Ext.BLANK;
    String annovarLocation = PSF.Ext.BLANK;
    ANNOTATION_BUILD annoBuild = ANNOTATION_BUILD.HG19;
    String vcfToAnnotate = null;
    boolean annotate = true;
    boolean ignoreInbreeding = false;
    String regionsFile = null;
    String logFile = "GATK_GENOTYPE.log";
    boolean skipRecalibration = false;
    SEQ_TARGET seqTarget = null;
    MergeVCF mergeVCF = null;
    ANNOVCF annoVCF = null;
    int argNum = 1;
    String usage = "\n" + "seq.GATK_Genotyper requires 2 argument\n";
    usage += "   (" + argNum++ + ") root input directory (i.e. " + GATK_LanePrep.ROOT_INPUT_COMMAND
             + rootInputDir + " (no default))\n" + "";
    usage += "   (" + argNum++ + ") root output directory (i.e. "
             + GATK_LanePrep.ROOT_OUTPUT_COMMAND + rootOutputDir + " (no default))\n" + "";
    usage += "   (" + argNum++ + ") tab-delimited file with no header of (i.e. " + FILE_OF_GVCFS
             + fileOfGVCFs + " (optional, no default))\n" + "";
    usage += "   (" + argNum++ + ") the full path to a  reference genome in fasta format (i.e."
             + GATK_LanePrep.REFERENCE_GENOME_COMMAND + referenceGenomeFasta + " (no default))\n"
             + "";
    usage += "   (" + argNum++ + ") the full path to the GATK executable (i.e. "
             + GATK.GATK_LOCATION_COMMAND + gATKLocation + " (defaults to systems path))\n" + "";

    usage += "   (" + argNum++ + ") run in quiet mode (i.e. " + GATK_LanePrep.QUIET_COMMAND
             + " (not the default))\n" + "";
    usage += "   (" + argNum++ + ") number of  threads for analysis(i.e." + NUM_THREADS + numThreads
             + " (default))\n" + "";

    usage += "   (" + argNum++ + ") filename for a log (i.e. " + GATK_LanePrep.LOG_FILE_COMMAND
             + logFile + " (default))\n" + "";
    usage += "   (" + argNum++ + ") over-write exsisting files (i.e. "
             + GATK_LanePrep.OVERWRITE_EXISTING_COMMAND + " (not the default))\n" + "";
    usage += "   (" + argNum++
             + ") set up a batch analysis for the root input directory for a log (i.e. "
             + GATK_LanePrep.BATCH_COMMAND + " (not the default))\n" + "";
    usage += "   (" + argNum++ + ") root output for analysis (i.e. " + OUTPUT_COMMAND + output
             + " ( default))\n" + "";
    for (GATK.RESOURCE training : GATK.RESOURCE.values()) {
      usage += "   (" + argNum++ + ") " + training.getName() + " Training Reference (i.e. "
               + training.getName() + "=example.vcf (no default))\n";
    }
    usage += "   (" + argNum++ + ") full path to the SNP EFF directory (i.e. "
             + SNPEFF.SNP_EFF_COMMAND + " ( no default))\n" + "";
    usage += "   (" + argNum++ + ") the build version for SNP EFF annotation (options are "
             + ArrayUtils.toStr(ANNOTATION_BUILD.values(), ", ") + " (i.e. "
             + SNPEFF.SNP_EFF_BUILD_COMMAND + annoBuild.toString() + " ( default))\n" + "";
    usage += "   (" + argNum++
             + ") full path to the SNP SIFT directory (only if different from the SNP EFF directory) (i.e. "
             + SNPSIFT.SNP_SIFT_LOCATION_COMMAND + " ( no default))\n" + "";
    usage += "   (" + argNum++ + ") do not annotate with SNP EFF/SIFT/ANNOVAR (i.e. "
             + SNPEFF.SNP_EFF_NO_ANNO_COMMAND + " ( not the default))\n" + "";
    usage += "   (" + argNum++ + ") full path to the ANNOVAR directory (i.e. "
             + ANNOVAR.ANNOVAR_COMMAND + " ( no default))\n" + "";
    usage += "   (" + argNum++
             + ") full path to a file for restricting the vcf...and subsequent recalibrations (i.e. "
             + GATK_LanePrep.REGIONS_FILE_COMMAND + " ( no default))\n" + "";

    usage += "   (" + argNum++ + ") annotate this vcf only (skipping all previous steps) (i.e. "
             + ANNOTATE_VCF + " ( no default))\n" + "";
    usage += "   (" + argNum++
             + ") Don't use Inbreeding Coefficient in variant filtering (use when <10 samples or highly related) (i.e. "
             + IGNORE_INBREEDING + " (not the default))";
    usage += "   (" + argNum++ + ") Region targeted by sequencing ("
             + ArrayUtils.toStr(SEQ_TARGET.values(), ",") + ") ( (i.e. "
             + GATK.TARGETED_REGION_COMMAND + " ( no default))\n" + "";
    usage += "   (" + argNum++
             + ") merge in these vcfs prior to annotating( like ARIC) (ex. MERGE_WITH=ARIC:aric1.vcf,aric2.vcf) ( (i.e. "
             + MERGE_WITH + " ( no default))\n" + "";
    usage += "   (" + argNum++
             + ") use another vcf to add more annotations (ex. EXTRA_VCF_ANNOTATIONS=charge.vcf:charge.maf1,charge.maf2) ( (i.e. "
             + EXTRA_VCF_ANNOTATIONS + " ( no default))\n" + "";

    usage += "   (" + argNum++ + ") restrict genotyping to a specific contig ( (i.e. " + "chrM"
             + " ( no default))\n" + "";
    usage += "   (" + argNum++
             + ") since targeted sequencing, or mtDNA genotyping will generally fail variant recalibration, skip it and see caveats here https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php ( (i.e. "
             + "-skipRecal" + " ( not the default))\n" + "";

    // usage += " (" + argNum++ + ") log file name (i.e. " + MILLS + " ( no default))\n" + "";

    for (String arg : args) {

      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(GATK_LanePrep.ROOT_INPUT_COMMAND)) {
        rootInputDir = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.ROOT_OUTPUT_COMMAND)) {
        rootOutputDir = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(FILE_OF_GVCFS)) {
        fileOfGVCFs = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.REFERENCE_GENOME_COMMAND)) {
        referenceGenomeFasta = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.LOG_FILE_COMMAND)) {
        logFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(NUM_THREADS)) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.MEMORY_MB)) {
        memoryInMB = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("wallTimeInHours=")) {
        wallTimeInHours = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.QUIET_COMMAND)) {
        verbose = false;
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.BATCH_COMMAND)) {
        batch = true;
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.OVERWRITE_EXISTING_COMMAND)) {
        overwriteExisting = true;
        numArgs--;
      } else if (arg.startsWith(GATK.GATK_LOCATION_COMMAND)) {
        gATKLocation = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(SNPEFF.SNP_EFF_COMMAND)) {
        snpEffLocation = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(SNPSIFT.SNP_SIFT_LOCATION_COMMAND)) {
        snpSiftLocation = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(SNPEFF.SNP_EFF_NO_ANNO_COMMAND)) {
        annotate = false;
        numArgs--;
      } else if (arg.startsWith("-skipRecal")) {
        skipRecalibration = true;
        numArgs--;
      } else if (arg.startsWith(OUTPUT_COMMAND)) {
        output = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(ANNOVAR.ANNOVAR_COMMAND)) {
        annovarLocation = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(ANNOTATE_VCF)) {
        vcfToAnnotate = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(IGNORE_INBREEDING)) {
        ignoreInbreeding = true;
        numArgs--;
      } else if (arg.startsWith(GATK_LanePrep.REGIONS_FILE_COMMAND)) {
        regionsFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith(GATK.TARGETED_REGION_COMMAND)) {
        try {
          seqTarget = SEQ_TARGET.valueOf(ext.parseStringArg(arg));
          numArgs--;
        } catch (IllegalArgumentException iae) {
          System.err.println(GATK.TARGETED_REGION_COMMAND + " must be one of: "
                             + ArrayUtils.toStr(SEQ_TARGET.values()));
        }
      } else if (arg.startsWith(MERGE_WITH)) {
        mergeVCF = MergeVCF.fromArg(arg);
        numArgs--;
      } else if (arg.startsWith(EXTRA_VCF_ANNOTATIONS)) {
        annoVCF = ANNOVCF.fromArg(arg);
        numArgs--;
      } else if (arg.startsWith(SNPEFF.SNP_EFF_BUILD_COMMAND)) {
        annoBuild = ANNOTATION_BUILD.valueOf(ext.parseStringArg(arg));
        numArgs--;
      } else if (ext.startsWithOneOf(arg, GATK.RESOURCE.names())) {
        String name = arg.substring(0, arg.indexOf('='));
        trainingResources.put(GATK.RESOURCE.getResourceByName(name), ext.parseStringArg(arg));
        numArgs--;
      }

      else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (rootOutputDir != null) {
      new File(rootOutputDir).mkdirs();
    }
    Logger log = new Logger(rootOutputDir + logFile);
    if (snpSiftLocation.equals(PSF.Ext.BLANK)) {
      snpSiftLocation = snpEffLocation;
    }
    if (vcfToAnnotate != null) {
      log.reportTimeInfo("Attempting to annotate " + vcfToAnnotate);
      annotateOnly(vcfToAnnotate, gATKLocation, referenceGenomeFasta, seqTarget, memoryInMB,
                   fileOfGVCFs, snpEffLocation, snpSiftLocation, annovarLocation, annoBuild, null,
                   verbose, overwriteExisting, log);
    } else {
      jointGenotype(rootInputDir, rootOutputDir, output, gATKLocation, referenceGenomeFasta,
                    fileOfGVCFs, trainingResources, snpEffLocation, snpSiftLocation,
                    annovarLocation, annoBuild, regionsFile, seqTarget, mergeVCF, annoVCF, verbose,
                    overwriteExisting, batch, annotate, skipRecalibration, ignoreInbreeding,
                    numThreads, memoryInMB, wallTimeInHours, log);
    }
  }
}
