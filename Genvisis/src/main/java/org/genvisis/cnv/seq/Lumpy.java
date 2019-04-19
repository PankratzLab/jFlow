/**
 * 
 */
package org.genvisis.cnv.seq;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

/**
 * wrapper for lumpy https://github.com/arq5x/lumpy-sv
 */
public class Lumpy {

  /**
   * 
   */
  private static final String SVTYPER = "svtyper";
  /**
   * 
   */
  private static final String BAM = "bam";
  /**
   * 
   */
  private static final String LUMPYEXPRESS = "lumpyexpress";

  private Lumpy() {

  }

  /**
   * @param lumpyExpressLoc full path to lumpyexpress
   * @param toRun list of {@link LumpyPrep} to run
   * @param outputVCF
   * @param log
   */
  private static void runLumpy(String lumpyExpressLoc, String lumpyExcludeBed,
                               List<PairedEndSVAnalysis> toRun, String outputVCF, Logger log) {
    // -x FILE
    StringBuilder bArg = new StringBuilder("-B ");
    StringBuilder sArg = new StringBuilder("-S ");
    StringBuilder dArg = new StringBuilder("-D ");
    List<String> outputFiles = new ArrayList<>();
    outputFiles.add(outputVCF);
    List<String> inputFiles = new ArrayList<>();
    boolean first = true;
    for (PairedEndSVAnalysis lumpyPrep : toRun) {
      inputFiles.add(lumpyPrep.getBaseBam());
      inputFiles.add(lumpyPrep.getSplitterBam());
      inputFiles.add(lumpyPrep.getDiscordantBam());

      bArg.append(first ? lumpyPrep.getBaseBam() : "," + lumpyPrep.getBaseBam());
      sArg.append(first ? lumpyPrep.getSplitterBam() : "," + lumpyPrep.getSplitterBam());
      dArg.append(first ? lumpyPrep.getDiscordantBam() : "," + lumpyPrep.getDiscordantBam());

      first = false;
    }
    List<String> cmd = new ArrayList<>();
    cmd.add(lumpyExpressLoc);
    cmd.add("-P");// output probability curves for each variant
                  // https://github.com/arq5x/lumpy-sv/issues/108?
    cmd.add(bArg.toString());
    cmd.add(sArg.toString());
    cmd.add(dArg.toString());
    cmd.add("-o");
    cmd.add(outputVCF);
    if (lumpyExcludeBed != null) {
      cmd.add("-x"); // see http://www.nature.com/ng/journal/v49/n5/full/ng.3834.html
      cmd.add(lumpyExcludeBed);
      inputFiles.add(lumpyExcludeBed);
    }

    CmdLine.runCommandWithFileChecks(cmd, "", inputFiles, outputFiles, true, false, false, true,
                                     log);

  }

  private static void run(String lumpyExpressLoc, String lumpyExcludeBed, String svtyperLoc,
                          String bam, String outDir, int threads) {
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + ext.rootOf(bam) + "_lumpy.log");
    List<String> bams = new ArrayList<>();
    bams.add(bam);
    List<PairedEndSVAnalysis> preps = LumpyPrep.runPrep(bams, outDir, threads, log);
    log.reportTimeWarning("Lumpy analysis assumes that the following are on your systempath\ngawk\nvcf-sort\n\nWindows users are gonna have a bad time");

    String outputVCF = outDir
                       + (preps.size() > 1 ? "lumpy.vcf"
                                           : BamOps.getSampleName(preps.get(0).getBaseBam(), log)
                                             + ".lumpy.vcf");
    runLumpy(lumpyExpressLoc, lumpyExcludeBed, preps, outputVCF, log);

    for (PairedEndSVAnalysis pe : preps) {
      String svtyperVCF = outDir + BamOps.getSampleName(pe.getBaseBam(), log) + "lumpy.gt.vcf";
      SVTyper.run(svtyperLoc, outputVCF, svtyperVCF, pe, log);
      sortVcf(svtyperVCF, log);
    }
  }

  private static String sortVcf(String inputVCF, Logger log) {
    String sortVCF = VCFOps.getAppropriateRoot(inputVCF, false) + ".sort.vcf";

    String bat = VCFOps.getAppropriateRoot(inputVCF, false) + ".sort.bat";
    String script = "vcf-sort " + inputVCF + " > " + sortVCF;

    CmdLine.prepareBatchForCommandLine(bat, true, log, new String[] {script});
    CmdLine.runCommandWithFileChecks(new String[] {bat}, "", new String[] {bat, inputVCF},
                                     new String[] {sortVCF}, true, false, false, log);

    VCFOps.verifyIndex(sortVCF, log);
    return sortVCF;

  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(Lumpy.class);
    c.addArgWithDefault(LUMPYEXPRESS,
                        "full path to lumpyexpress, see https://github.com/arq5x/lumpy-sv",
                        LUMPYEXPRESS);
    c.addArgWithDefault(SVTYPER,
                        "full path to svtyper script, see https://github.com/hall-lab/svtyper.git",
                        SVTYPER);

    c.addArgWithDefault("bed",
                        "(Optional) full path to a bed file to exclude calls in, see https://github.com/hall-lab/speedseq/blob/30f8dc688ca37b74a3df2cf357fc7c64eaa74b7e/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed",
                        "a.bed");
    c.addArgWithDefault(BAM, "Bam file to analyze", "a.bam");
    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "na");
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "1");

    c.parseWithExit(args);
    run(c.get(LUMPYEXPRESS), c.get("bed"), c.get(SVTYPER), c.get(BAM), c.get(CLI.ARG_OUTDIR),
        c.getI(CLI.ARG_THREADS));
  }

}

// NOTES: http://shiny.wehi.edu.au/cameron.d/sv_benchmark/ comp of callers

// https://www.biostars.org/p/243145/
// Ideally, you would filter on variant quality directly, instead of RP and SR independently but not
// many SV callers actually report a variant quality score, and none that I know of report an even
// remotely calibrated quality score.

// https://www.nature.com/articles/npjgenmed201626
