package org.genvisis.one.spencer.ewingGATK;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.genvisis.seq.analysis.GATK;
import org.genvisis.seq.analysis.GATK.RESOURCE;
import org.genvisis.seq.analysis.GATK.SEQ_TARGET;
import org.genvisis.seq.analysis.RegNovo;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class Pipeline {

  private static final String BAM = ".bam";
  private static final String BAI = ".bai";
  private static final String PBS = ".pbs";

  private static final String QSUB = "qsub -A $MYGROUP -W group_list=$MYGROUP ";

  private static final int PANKRATZ_CORES = 24;
  private static final int PANKRATZ_MEM_GB = 252;

  private static final int SMALL_CORES = 24;
  private static final int SMALL_MEM_GB = 62;

  private static final Joiner FILENAME_JOINER = Joiner.on("_");

  private static final String BASE_DIR = "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/";
  private static final String REF_GENO = BASE_DIR + "WUSTL_RefGenome/all_sequences.fa";
  private static final String REG_NOVO_DIR = BASE_DIR + "RegNovo/";
  private static final String VCF_DIR = BASE_DIR + "vcfs/";

  private static final String RECAL_SNP_FILE = VCF_DIR + "ES_recalibrate_SNP.recal";
  private static final String RECAL_SNP_TRANCHES_FILE = VCF_DIR + "ES_recalibrate_SNP.tranches";
  private static final String RECAL_SNP_RSCRIPT = VCF_DIR + "ES_recalibrate_SNP_plots.R";

  private static final String RECAL_INDEL_FILE = VCF_DIR + "ES_recalibrate_INDEL.recal";
  private static final String RECAL_INDEL_TRANCHES_FILE = VCF_DIR + "ES_recalibrate_INDEL.tranches";
  private static final String RECAL_INDEL_RSCRIPT = VCF_DIR + "ES_recalibrate_INDEL_plots.R";

  private final String bamDir;
  private final Set<String> ids;

  public Pipeline(String bamDir) {
    this.bamDir = bamDir;
    String[] bams = Files.list(bamDir, "", ".bam", true);
    ids = Arrays.stream(bams).map(ext::rootOf).filter(id -> !id.contains("."))
                .filter(id -> !id.contains("_")).collect(Collectors.toSet());
  }

  public void generatePipeline() {
    // Map<String, String> toIndex = indexBAMs();
    // Map<String, String> toRecal = runRecal();
    // Map<String, String> toCall = runHaplotypeCaller();

    // Set<String> pbsToIndex = Sets.newTreeSet();
    //
    // for (Map.Entry<String, String> indexEntry : toIndex.entrySet()) {
    // String id = indexEntry.getKey();
    // String pbs = indexEntry.getValue();
    //
    // pbsToIndex.add(QSUB + pbs);
    // }
    // Files.writeIterable(pbsToIndex, bamDir + "../pbsScripts/pbsToIndex");

    // Map<String, String> toRecal = runRecal();
    // Map<String, String> toCall = runHaplotypeCaller();

    // runHaplotypeCallerItasca();

    // runCombineGVCFs();

    // runCombinedGVCFJointGenotyping();

    // runRegNovoPreVQSR();

    runVQSR();
    // runVariantAnnotator();

  }

  private Map<String, String> indexBAMs() {
    Map<String, String> toIndex = Maps.newHashMap();
    for (String id : ids) {
      if (!fileExists(bamIndex(id))) {
        StringJoiner cmd = new StringJoiner("\n");
        cmd.add(cdBamDir());
        cmd.add("");
        cmd.add("java -jar /home/pankrat2/public/bin/picard/picard.jar BuildBamIndex INPUT="
                + bam(id));
        cmd.add("");
        // cmd.add(QSUB + recalPBS(id));
        Qsub.qsubGb(indexPBS(id), cmd.toString(), SMALL_MEM_GB / SMALL_CORES, 2, 1);
        toIndex.put(id, indexPBS(id));
      }
    }
    return toIndex;
  }

  private Map<String, String> runRecal() {
    final int threads = 8;
    Map<String, String> toRecal = Maps.newHashMap();
    for (String id : ids) {
      // if (!fileExists(recal1(id), recalPBS(id))) {
      if (true) {
        StringJoiner cmd = new StringJoiner("\n");
        cmd.add(cdBamDir());
        cmd.add("");
        cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
           .add("\t-T BaseRecalibrator \\").add("\t-R " + REF_GENO + " \\")
           .add("\t-I " + bam(id) + " \\")
           .add("\t-knownSites /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
           .add("\t-knownSites /home/pankrat2/public/bin/ref/Mills_and_1000G_gold_standard.indels.b37.vcf \\")
           .add("\t-o " + recal1(id) + " \\").add("\t-nct 4");
        cmd.add("");
        cmd.add("qsub -q batch " + callerPBS(id));
        // cmd.add(QSUB + callerPBS(id));
        // cmd.add("");
        // cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
        // .add("\t-T BaseRecalibrator \\")
        // .add("\t-R " + REF_GENO + " \\")
        // .add("\t-I " + bam(id) + " \\")
        // .add("\t-knownSites /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
        // .add("\t-knownSites
        // /home/pankrat2/public/bin/ref/Mills_and_1000G_gold_standard.indels.b37.vcf \\")
        // .add("\t-o " + recal2(id) + " \\")
        // .add("\t-BQSR " + recal1(id))
        // .add("\t-nct " + threads);
        // cmd.add("");
        // cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
        // .add("\t-T AnalyzeCovariates \\")
        // .add("\t-R " + REF_GENO + " \\")
        // .add("\t-before " + recal1(id) + " \\")
        // .add("\t-after " + recal2(id) + " \\")
        // .add("\t-plots " + recalPlots(id));
        Qsub.qsubGb(recalPBS(id), cmd.toString(), 22, 20, threads);
        toRecal.put(id, recalPBS(id));
      }
    }
    return toRecal;
  }

  private Map<String, String> runHaplotypeCaller() {
    final int threads = 4;
    Map<String, String> toCall = Maps.newHashMap();
    for (String id : ids) {
      // if (!fileExists(callerGvcf(id), callerPBS(id))) {
      if (!fileExists(callerGvcf(id), callerGvcfIndex(id))) {
        StringJoiner cmd = new StringJoiner("\n");
        cmd.add(cdBamDir());
        cmd.add("");
        cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
           .add("\t-T HaplotypeCaller \\").add("\t-R " + REF_GENO + " \\")
           .add("\t-I " + bam(id) + " \\").add("\t-BQSR " + recal1(id) + " \\")
           .add("\t--genotyping_mode DISCOVERY \\").add("\t-o " + callerGvcf(id) + " \\")
           .add("\t--dbsnp /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
           .add("\t-ERC GVCF \\").add("\t-nct " + threads);

        Qsub.qsubGb(callerPBS(id), cmd.toString(), 20, 40, threads);
        toCall.put(id, callerPBS(id));
      }
    }
    return toCall;
  }

  private void runHaplotypeCallerItasca() {
    final int threads = 4;
    Iterator<String> idIter = ids.iterator();
    Set<String> itascaIDs = HashVec.loadFileToHashSet(BASE_DIR + "bams/itascaIDs", false);
    Set<String> pbsScripts = Sets.newHashSet();
    while (idIter.hasNext()) {
      Set<String> idsToBatch = Sets.newHashSetWithExpectedSize(3);
      int i = 0;
      while (i < 3 && idIter.hasNext()) {
        String id = idIter.next();
        if (!fileExists(callerGvcf(id)) && itascaIDs.contains(id)) {
          idsToBatch.add(id);
          i++;
        }
      }

      StringJoiner cmd = new StringJoiner("\n");
      cmd.add(cdBamDir());
      for (String id : idsToBatch) {
        cmd.add("");
        cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
           .add("\t-T HaplotypeCaller \\").add("\t-R " + REF_GENO + " \\")
           .add("\t-I " + bam(id) + " \\").add("\t-BQSR " + recal1(id) + " \\")
           .add("\t--genotyping_mode DISCOVERY \\").add("\t-o " + callerGvcf(id) + " \\")
           .add("\t--dbsnp /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
           .add("\t-ERC GVCF \\").add("\t-nct " + threads + " &");
        cmd.add("");
      }

      cmd.add("wait");
      String pbs = callerPBS(Joiner.on("_").join(idsToBatch));
      Qsub.qsubGb(pbs, cmd.toString(), 62, 40, 16);
      pbsScripts.add(QSUB + " -q sb " + pbs);
    }
    Files.writeIterable(pbsScripts, BASE_DIR + "bams/runItascaBatch");

  }

  private void runSplitGVCFs() {

    final int threads = 1;
    List<String> commands = Lists.newArrayList();
    int cmdCount = 0;
    int commandScriptCount = 0;
    for (String id : ids) {
      boolean write = false;
      StringJoiner cmd = new StringJoiner("\n");
      cmd.add(cdBamDir());
      cmd.add("");
      for (int i = 2; i < 26; i++) {

        // if (!fileExists(callerGvcf(id), callerPBS(id))) {
        String chr = Integer.toString(i);
        if (i == 23) chr = "X";
        if (i == 24) chr = "Y";
        if (i == 25) chr = "MT";
        if (!fileExists(splitGvcf(id, chr), splitGvcfIndex(id, chr))) {
          write = true;
          cmd.add("");
          cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
             .add("\t-T SelectVariants \\").add("\t-R " + REF_GENO + " \\")
             .add("\t-V " + callerGvcf(id) + " \\").add("\t-o " + splitGvcf(id, chr) + " \\")
             .add("\t-L " + chr);
          cmd.add("");

        }
      }
      if (write) {
        Qsub.qsub(splitPBS(id), cmd.toString(), 2500, 24, threads);
        commands.add("qsub " + (cmdCount < 48 ? "-q pankratz " : "") + splitPBS(id));
        if (cmdCount++ > 400) {
          Files.writeIterable(commands,
                              BASE_DIR + "pbsScripts/splitGVCFJobs_" + commandScriptCount);
          commands.clear();
          commandScriptCount++;
          cmdCount = 0;
        }
      }
    }
    Files.writeIterable(commands, BASE_DIR + "pbsScripts/splitGVCFJobs_" + commandScriptCount);
  }

  private void runCombineGVCFs() {

    final int threads = 1;

    for (int i = 1; i < 26; i++) {

      String chr = intToChr(i);
      for (int batch = 1; batch <= 6; batch++) {
        String gvcfList = BASE_DIR + "pbsScripts/chr" + chr + "_gVCFs_batch_" + batch + ".list";
        // Files.writeIterable(ids.stream().map(id -> splitGvcf(id,
        // chr)).collect(Collectors.toList()),
        // gvcfList);

        String pbsScript = BASE_DIR + "pbsScripts/CG_chr" + chr + "_gVCFs_batch_" + batch
                           + "_combineGVCFs.pbs";

        if (!fileExists(combinedGVCF(chr, batch), combinedGVCFIndex(chr, batch))) {

          StringJoiner cmd = new StringJoiner("\n");
          cmd.add(cdBamDir());
          cmd.add("");
          cmd.add("java -Xmx14G -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
             .add("\t-T CombineGVCFs \\").add("\t-R " + REF_GENO + " \\")
             .add("\t-V " + gvcfList + " \\").add("\t-o " + combinedGVCF(chr, batch) + " \\")
             .add("\t-L " + chr);
          cmd.add("");

          Qsub.qsubGb(pbsScript, cmd.toString(), 16, 96, threads);
        }

      }
    }

  }

  private static String combinedGVCF(String chr, int batch) {
    return "ES_combined_chr" + chr + "_batch_" + batch + "_presplit.g.vcf.gz";
  }

  private static String combinedGVCFIndex(String chr, int batch) {
    return combinedGVCF(chr, batch) + ".tbi";
  }

  private void runJointGenotyping() {

    final int threads = 1;
    for (int i = 1; i < 26; i++) {
      StringJoiner cmd = new StringJoiner("\n");
      cmd.add(cdBamDir());
      cmd.add("");
      final String chr = intToChr(i);

      String chrGVCFs = BASE_DIR + "chr" + chr + "_gVCFs.list";

      Files.writeIterable(ids.stream()
                             .map(s -> BASE_DIR + "gVCFs_byChr/" + s + "_chr" + chr
                                       + "_recalibrated.snps.indels.g.vcf.gz")
                             .collect(Collectors.toList()),
                          chrGVCFs);

      if (!fileExists(jointVCF(chr), jointVCFIndex(chr))) {
        cmd.add("");
        cmd.add("java -Xmx15G -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
           .add("\t-T GenotypeGVCFs \\").add("\t-R " + REF_GENO + " \\")
           .add("\t-V " + chrGVCFs + " \\").add("\t-o " + jointVCF(chr) + " \\").add("\t-L " + chr);
        cmd.add("");
        Qsub.qsubGb(jointGenoPBS(chr), cmd.toString(), 19, 440, threads);
      }
    }

  }

  private void runCombinedGVCFJointGenotyping() {

    final int threads = 1;
    for (int i = 1; i < 26; i++) {

      final String chr = intToChr(i);

      if (!Files.checkAllFiles("", combinedGVCFjointVCF(chr), combinedGVCFjointVCF(chr))) {
        boolean ready = true;
        StringJoiner cmd = new StringJoiner("\n");
        cmd.add(cdBamDir());
        cmd.add("");
        cmd.add("");
        cmd.add("java -Xmx18G -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
           .add("\t-T GenotypeGVCFs \\").add("\t-R " + REF_GENO + " \\");
        for (int batch = 1; batch <= 6; batch++) {
          cmd.add("\t-V " + combinedGVCF(chr, batch) + " \\");
          if (!fileExists(combinedGVCF(chr, batch), combinedGVCFIndex(chr, batch))) {
            System.err.println(combinedGVCFIndex(chr, batch) + " does not exist");
            ready = false;
            break;
          }
        }
        cmd.add("\t-o " + combinedGVCFjointVCF(chr) + " \\").add("\t-L " + chr);
        cmd.add("");
        if (ready) Qsub.qsubGb(combinedGVCFjointGenoPBS(chr), cmd.toString(), 20, 288, threads);
      }
    }

  }

  private void runVQSR() {
    GATK gatk = new GATK("/home/pankrat2/public/bin/GATK_3.7/", REF_GENO, null, SEQ_TARGET.GENOME,
                         null, 250000, null, null, true, false, new Logger("VQSR.log"));

    Map<GATK.RESOURCE, String> trainingResources = Maps.newEnumMap(GATK.RESOURCE.class);
    String RESOURCE_DIRECTORY = "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/WUSTL_RefGenome/";
    trainingResources.put(RESOURCE.HAPMAP, RESOURCE_DIRECTORY + "hapmap_3.3.b37.vcf");
    trainingResources.put(RESOURCE.OMNI, RESOURCE_DIRECTORY + "1000G_omni2.5.b37.vcf");
    trainingResources.put(RESOURCE.DBSNP, RESOURCE_DIRECTORY + "dbsnp_138.b37.vcf");
    trainingResources.put(RESOURCE.G1K,
                          RESOURCE_DIRECTORY + "1000G_phase1.snps.high_confidence.b37.vcf");
    trainingResources.put(RESOURCE.MILLS,
                          RESOURCE_DIRECTORY + "Mills_and_1000G_gold_standard.indels.b37.vcf");
    gatk.setTrainingResources(trainingResources);
    int threads = 24;
    Joiner cmdJoiner = Joiner.on(' ');

    String[] jointGenoVCFs = HashVec.loadFileToStringArray(BASE_DIR + "jointGenotypeVCFs.list",
                                                           false, null, false);
    StringJoiner cmd = new StringJoiner("\n");
    cmd.add(cdVCFDirectory());
    if (!Files.checkAllFiles("", RECAL_SNP_FILE, RECAL_SNP_TRANCHES_FILE, RECAL_SNP_RSCRIPT)) {
      cmd.add(cmdJoiner.join(gatk.buildSNPRecalibrationModelCommand(BASE_DIR
                                                                    + "jointGenotypeVCFs.list",
                                                                    RECAL_SNP_FILE,
                                                                    RECAL_SNP_TRANCHES_FILE,
                                                                    RECAL_SNP_RSCRIPT, false,
                                                                    threads)));
    }
    if (!Files.checkAllFiles("", snpRecaledVCF(), snpRecaledVCF() + ".tbi")) {
      cmd.add("");
      cmd.add(cmdJoiner.join(gatk.applySNPRecalibrationModelCommand(RECAL_SNP_FILE,
                                                                    RECAL_SNP_TRANCHES_FILE,
                                                                    snpRecaledVCF(), threads,
                                                                    jointGenoVCFs)));
      cmd.add("");
    }
    if (!Files.checkAllFiles("", RECAL_INDEL_FILE, RECAL_INDEL_TRANCHES_FILE,
                             RECAL_INDEL_RSCRIPT)) {
      cmd.add("");
      cmd.add(cmdJoiner.join(gatk.buildINDELRecalibrationModelCommand(snpRecaledVCF(),
                                                                      RECAL_INDEL_FILE,
                                                                      RECAL_INDEL_TRANCHES_FILE,
                                                                      RECAL_INDEL_RSCRIPT, false,
                                                                      threads)));
      cmd.add("");
    }
    if (!Files.checkAllFiles("", snpIndelRecaledVCF(), snpIndelRecaledVCF() + ".tbi")) {
      cmd.add("");
      cmd.add(cmdJoiner.join(gatk.applyINDELRecalibrationModelCommand(RECAL_INDEL_FILE,
                                                                      RECAL_INDEL_TRANCHES_FILE,
                                                                      snpIndelRecaledVCF(), threads,
                                                                      snpRecaledVCF())));
    }
    if (!cmd.toString().equals(cdVCFDirectory())) {
      Qsub.qsubGb(BASE_DIR + "VQSR_ES.pbs", cmd.toString(), 252, 175, threads);
    }
  }

  private void runVariantAnnotator() {
    for (int i = 1; i < 26; i++) {
      String chr = intToChr(i);
      StringJoiner cmd = new StringJoiner("\n");
      cmd.add(cdVCFDirectory());
      cmd.add("");
      cmd.add("java -Xmx8g -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\");
      cmd.add("-T VariantAnnotator \\");
      cmd.add("-V " + combinedGVCFjointVCF(chr) + " \\");
      cmd.add("-R " + REF_GENO + " \\");
      cmd.add("-A InbreedingCoeff \\");
      cmd.add("-ped /panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/sampleInfo/ES.ped \\");
      cmd.add("--dbsnp /panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/WUSTL_RefGenome/dbsnp_138.b37.vcf \\");
      cmd.add("-o " + combinedGVCFjointVCFExtraAnno(chr));

      Qsub.qsub(BASE_DIR + "pbsScripts/ICA_chr" + chr + "_AnnotateInbreedingCoeff.pbs",
                cmd.toString(), 10500, 100, 1);
    }
  }

  private void runRegNovoPreVQSR() {
    for (int i = 1; i <= 25; i++) {
      String chr = intToChr(i);
      StringJoiner cmd = new StringJoiner("\n");
      cmd.add(cdRegNovoDir());
      cmd.add("java -jar ~/genvisis_RegNovo.jar " + RegNovo.class.getName() + " \\");
      cmd.add("\t" + "vcf=../vcfs/ES_JointGenotypes_chr" + chr + "_combinedGVCFs.vcf.gz \\");
      cmd.add("\t" + "vpop=ES_Trios.vpop \\");
      cmd.add("\t" + "regions=chr" + chr + " \\");
      cmd.add("\t" + "outputRoot=ES_RegNovo_PreVQSR_chr" + chr + " \\");
      cmd.add("\t" + "log=ES_RegNovo_PreVQSR_chr" + chr + ".log");

      Qsub.qsub(REG_NOVO_DIR + "RNPV_chr" + chr + "_RegNovoPreVQSR.pbs", cmd.toString(), 10500, 200,
                1);
    }
  }

  private String intToChr(int i) {
    String chr = Integer.toString(i);
    if (i == 23)
      chr = "X";
    else if (i == 24)
      chr = "Y";
    else if (i == 25) chr = "MT";
    return chr;
  }

  private String cdBamDir() {
    return "cd " + bamDir;
  }

  private String cdRegNovoDir() {
    return "cd " + REG_NOVO_DIR;
  }

  private String cdVCFDirectory() {
    return "cd " + VCF_DIR;
  }

  private String bam(String id) {
    return id + BAM;
  }

  private String bamIndex(String id) {
    return id + BAI;
  }

  private String recal1(String id) {
    return id + ".recal_data.table";
  }

  private String recal2(String id) {
    return id + ".post_recal_data.table";
  }

  private String recalPlots(String id) {
    return id + ".recalibration_plots.pdf";
  }

  private String[] recalOuts(String id) {
    return new String[] {recal1(id), recal2(id), recalPlots(id)};
  }

  private String callerGvcf(String id) {
    return "../gVCFs/" + id + "_recalibrated.snps.indels.g.vcf.gz";
  }

  private String splitGvcf(String id, String chr) {
    return "../gVCFs_byChr/" + id + "_chr" + chr + "_recalibrated.snps.indels.g.vcf.gz";
  }

  private String jointVCF(String chr) {
    return BASE_DIR + "vcfs/ES_JointGenotypes_chr" + chr + "_splitGVCFs.vcf.gz";
  }

  private String jointVCFIndex(String chr) {
    return jointVCF(chr) + ".tbi";
  }

  private String combinedGVCFjointVCF(String chr) {
    return BASE_DIR + "vcfs/ES_JointGenotypes_chr" + chr + "_combinedGVCFs.vcf.gz";
  }

  private String combinedGVCFjointVCFExtraAnno(String chr) {
    return BASE_DIR + "vcfs/ES_JointGenotypes_chr" + chr + "_combinedGVCFs_withPED.vcf.gz";
  }

  private String combinedGVCFjointVCFIndex(String chr) {
    return combinedGVCFjointVCF(chr) + ".tbi";
  }

  private String callerGvcfIndex(String id) {
    return "../gVCFs/" + id + "_recalibrated.snps.indels.g.vcf.gz.tbi";
  }

  private String splitGvcfIndex(String id, String chr) {
    return splitGvcf(id, chr) + ".tbi";
  }

  private String indexPBS(String id) {
    return "../pbsScripts/I_" + id + "_index" + PBS;
  }

  private String recalPBS(String id) {
    return "../pbsScripts/R_" + id + "_recal" + PBS;
  }

  private String callerPBS(String id) {
    return "../pbsScripts/C_" + id + "_haplotypeCaller" + PBS;
  }

  private String splitPBS(String id) {
    return "../pbsScripts/S_" + id + "_splitGVCFs" + PBS;
  }

  private String jointGenoPBS(String chr) {
    return "../pbsScripts/J_chr" + chr + "_jointGenotyping" + PBS;
  }

  private String combinedGVCFjointGenoPBS(String chr) {
    return "../pbsScripts/JC_chr" + chr + "_combinedGVCFs_jointGenotyping" + PBS;
  }

  private String snpRecaledVCF() {
    return VCF_DIR + "ES_recalibrated_snps_raw_indels.vcf.gz";
  }

  private String snpIndelRecaledVCF() {
    return VCF_DIR + "ES_recalibrated_snps_indels.vcf.gz";
  }

  private boolean fileExists(String... filenames) {
    return Files.exists(bamDir, filenames);
  }

  public static void main(String[] args) {
    new Pipeline(BASE_DIR + "bams/").generatePipeline();
  }

}
