package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCFTumorNormalOps;

/**
 * Wrapper and parser/assembler for Somatic Sniper (https://github.com/genome/somatic-sniper)
 */
public class SomaticSniper {

  private static void run(GATK gatk, SomaticParams params, String bams, String vcfPop,
                          int numThreads, Logger log) {
    VcfPopulation vpop = VcfPopulation.load(vcfPop, POPULATION_TYPE.TUMOR_NORMAL, log);
    String[] bamFiles = null;
    if (Files.isDirectory(bams)) {
      bamFiles = Files.list(bams, ".bam");
    } else {
      bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] {0}, true);
    }
    log.reportTimeInfo("Detected " + bamFiles.length + " total bam files");

    TNSample[] tnSamples = matchSamples(bamFiles, params.getOutputDir(), vpop, params, log);
    String tnMatch = params.getOutputDir() + "tnMatch.txt";
    ArrayList<String> summaryMatch = new ArrayList<String>();
    summaryMatch.add("NORMAL_BAM\tTUMOR_BAM");
    for (TNSample tnSample : tnSamples) {
      summaryMatch.add(tnSample.getNormalBam() + "\t" + tnSample.getTumorBam());
    }
    Files.writeArray(ArrayUtils.toStringArray(summaryMatch), tnMatch);
    TNProducer producer = new TNProducer(tnSamples);
    WorkerTrain<TNSample> train = new WorkerTrain<TNSample>(producer, numThreads, 2, log);
    ArrayList<String> finalOuts = new ArrayList<String>();
    while (train.hasNext()) {
      TNSample tmp = train.next();
      finalOuts.add(tmp.getOutputGz());
    }
    String out = params.getOutputDir() + "tn.out.vcf.gz";
    params.getOutputDir();

    if (gatk.mergeVCFs(ArrayUtils.toStringArray(finalOuts), out, numThreads, false, log)) {

    }

  }

  private static class TNProducer extends AbstractProducer<TNSample> {

    private final TNSample[] tnSamples;
    private int index;

    public TNProducer(TNSample[] tnSamples) {
      super();
      this.tnSamples = tnSamples;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < tnSamples.length;
    }

    @Override
    public Callable<TNSample> next() {
      TNSample current = tnSamples[index];
      index++;

      return current;
    }
  }

  private static TNSample[] matchSamples(String[] bamFiles, String outputDir, VcfPopulation vpop,
                                         SomaticParams somaticParams, Logger log) {
    Hashtable<String, String> all = new Hashtable<String, String>();
    for (String bamFile : bamFiles) {
      all.put(BamOps.getSampleName(bamFile), bamFile);
    }
    ArrayList<String> analysisBams = new ArrayList<String>();
    ArrayList<TNSample> tnSamples = new ArrayList<TNSample>();
    for (String tnPair : vpop.getSubPop().keySet()) {
      Set<String> samps = vpop.getSubPop().get(tnPair);
      String tumor = null;
      String normal = null;
      for (String samp : samps) {
        if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
          tumor = samp;
        } else if (vpop.getPopulationForInd(samp,
                                            RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.NORMAL)) {
          normal = samp;
        } else {
          throw new IllegalArgumentException("Unknown types");
        }
      }

      if (!all.containsKey(tumor) || !all.containsKey(normal)) {
        throw new IllegalArgumentException("Could not find bam file for Tumor " + tumor
                                           + " or for Normal " + normal);
      } else {
        TNSample tSample = new TNSample(all.get(normal), normal, all.get(tumor), tumor,
                                        somaticParams, log);
        tnSamples.add(tSample);
        analysisBams.add(all.get(normal));
        analysisBams.add(all.get(tumor));

      }
    }
    if (analysisBams.size() < bamFiles.length) {
      Files.writeArray(ArrayUtils.toStringArray(analysisBams),
                       outputDir + ext.rootOf(vpop.getFileName() + ".analysis.bams.txt"));
    }

    log.reportTimeInfo("Matched " + tnSamples.size() + " tumor normals with bams ");
    return tnSamples.toArray(new TNSample[tnSamples.size()]);
  }

  private static class SomaticParams {

    private final String somaticSnipeLoc;
    private final String refGenome;
    private final String outputDir;
    private final String samtoolsLoc;
    private final boolean overwriteExisting;

    public SomaticParams(String somaticSnipeLoc, String refGenome, String outputDir,
                         boolean overwriteExisting) {
      super();
      this.somaticSnipeLoc = somaticSnipeLoc;
      this.refGenome = refGenome;
      this.outputDir = outputDir;
      this.overwriteExisting = overwriteExisting;
      samtoolsLoc = somaticSnipeLoc + "build/vendor/samtools/";

    }

    public String getSamtoolsLoc() {
      return samtoolsLoc;
    }

    public String getOutputDir() {
      return outputDir;
    }

    public boolean isOverwriteExisting() {
      return overwriteExisting;
    }

    public String getRefGenome() {
      return refGenome;
    }

    public String getSnipe() {
      return somaticSnipeLoc + "build/bin/" + "bam-somaticsniper";
    }

    private String[] generateCommand(String normalBam, String tumorBam, String output) {
      ArrayList<String> command = new ArrayList<String>();
      command.add(getSnipe());
      command.add("-F");
      command.add("vcf");

      command.add("-Q");
      command.add("40");
      command.add("-G");
      command.add("-L");
      command.add("-f");
      command.add(refGenome);
      command.add(tumorBam);
      command.add(normalBam);
      command.add(output);
      return ArrayUtils.toStringArray(command);
    }
  }

  private static class TNSample implements Callable<TNSample> {

    private final String normalBam;
    private final String normalSample;
    private final String tumorBam;
    private final String tumorSample;
    private final SomaticParams somaticParams;
    private final String indelPile;
    private final String indelPileFilt;
    private final String output;
    private final String outputGz;
    private final Logger log;

    public TNSample(String normalBam, String normalSample, String tumorBam, String tumorSample,
                    SomaticParams somaticParams, Logger log) {
      super();
      this.normalBam = normalBam;
      this.normalSample = normalSample;
      this.tumorBam = tumorBam;
      this.tumorSample = tumorSample;
      this.somaticParams = somaticParams;
      output = somaticParams.getOutputDir() + "T_" + tumorSample + "_N_" + normalSample + ".vcf";
      outputGz = VCFOps.getAppropriateRoot(output, false) + "_rename.vcf.gz";
      indelPile = ext.rootOf(output, false) + ".indel.pileup";
      indelPileFilt = ext.rootOf(output, false) + ".indel.pileup";

      this.log = log;
    }

    public String getOutputGz() {
      return outputGz;
    }

    public String getNormalBam() {
      return normalBam;
    }

    public String getTumorBam() {
      return tumorBam;
    }

    private boolean runSamToolsIndel() {
      String[] inputs = new String[] {tumorBam, somaticParams.getRefGenome()};
      String[] output = new String[] {indelPile};
      ArrayList<String> commandCall = new ArrayList<String>();
      commandCall.add(somaticParams.getSamtoolsLoc() + "samtools");
      commandCall.add("pileup");
      commandCall.add("�cvi");
      commandCall.add("�f");
      commandCall.add(somaticParams.getRefGenome());
      commandCall.add(tumorBam);
      commandCall.add(indelPile);

      boolean progress = CmdLine.runCommandWithFileChecks(ArrayUtils.toStringArray(commandCall), "",
                                                          inputs, output, true,
                                                          somaticParams.isOverwriteExisting(),
                                                          false, log);
      if (progress) {
        ArrayList<String> commandFilt = new ArrayList<String>();
        commandFilt.add("perl");
        commandFilt.add(somaticParams.getSamtoolsLoc() + "misc/samtools.pl");
        commandFilt.add("varFilter");
        commandFilt.add(indelPile);
        commandFilt.add(">");
        commandFilt.add(indelPileFilt);
        String[] runBat = CmdLine.prepareBatchForCommandLine(ArrayUtils.toStringArray(commandFilt),
                                                             indelPileFilt + ".sh", true, log);
        progress = CmdLine.runCommandWithFileChecks(runBat, "", new String[] {indelPile},
                                                    new String[] {indelPileFilt}, true,
                                                    somaticParams.isOverwriteExisting(), false,
                                                    log);
      }
      return progress;

    }

    @Override
    public TNSample call() throws Exception {
      String[] inputs = new String[] {somaticParams.getRefGenome(), somaticParams.getSnipe(),
                                      tumorBam, normalBam};
      String[] outputs = new String[] {output};
      String[] command = somaticParams.generateCommand(normalBam, tumorBam, output);

      boolean progress = CmdLine.runCommandWithFileChecks(command, "", inputs, outputs, true,
                                                          somaticParams.isOverwriteExisting(),
                                                          false, log);

      if (progress) {

        progress = runSamToolsIndel();
        if (!Files.exists(outputGz)) {
          VCFTumorNormalOps.renameTumorNormalVCF(output, tumorSample, normalSample, outputGz,
                                                 VCFOps.getAppropriateRoot(outputGz,
                                                                           false) + ".filtered.vcf.gz",
                                                 log);
        }

      }
      if (progress) {
        return this;
      } else {
        return null;
      }
    }
  }

  public static void main(String[] args) {
    String bams = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/TN.vpop.analysis.bams";
    String vpop = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/TN.vpop";
    String outputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/";
    String refGenome = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
    String somaticSnipLoc = "/home/pankrat2/public/bin/somaticSniper/somatic-sniper/";
    String gatkLoc = "/home/pankrat2/public/bin/GATK/";

    int numThreads = 1;
    new File(outputDir).mkdirs();
    Logger log = new Logger(outputDir + "somaticSniper.log");
    GATK gatk = new GATK(gatkLoc, refGenome, null, null, null, PSF.Ext.DEFAULT_MEMORY_MB, true,
                         false, log);
    SomaticParams params = new SomaticParams(somaticSnipLoc, gatk.getReferenceGenomeFasta(),
                                             outputDir, false);
    run(gatk, params, bams, vpop, numThreads, log);

  }

}
