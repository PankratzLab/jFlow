/**
 * 
 */
package org.genvisis.seq.cnv;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;

/**
 * Wrapper for SBTyper - Bayesian genotyper for structural variants Using to genotype lumpy calls
 */
public class SVTyper {

  private SVTyper() {

  }

  /**
   * Run svtyper on Lumpy output
   * 
   * @param svtyperLoc full path the the SVTyper script
   * @param inputVCF vcf created by lumpy
   * @param outputVCF output of results
   * @param pairedEndSVAnalysis what to analyze
   * @param log
   */
  public static void run(String svtyperLoc, String inputVCF, String outputVCF,
                         PairedEndSVAnalysis pairedEndSVAnalysis, Logger log) {
    new File(ext.parseDirectoryOfFile(outputVCF)).mkdirs();
    String bat = outputVCF + ".bat";

    List<String> inputFiles = new ArrayList<>();
    inputFiles.add(inputVCF);
    inputFiles.add(svtyperLoc);

    inputFiles.add(pairedEndSVAnalysis.getBaseBam());
    BamOps.verifyIndex(pairedEndSVAnalysis.getBaseBam(), log);

    inputFiles.add(pairedEndSVAnalysis.getSplitterBam());
    BamOps.verifyIndex(pairedEndSVAnalysis.getSplitterBam(), log);

    inputFiles.add(bat);

    List<String> outputFiles = new ArrayList<>();
    outputFiles.add(outputVCF);

    StringBuilder cmd = new StringBuilder(svtyperLoc);
    cmd.append(" -B " + pairedEndSVAnalysis.getBaseBam());
    cmd.append(" -S " + pairedEndSVAnalysis.getSplitterBam());
    cmd.append(" -i " + inputVCF);
    cmd.append(" > " + outputVCF);
    CmdLine.prepareBatchForCommandLine(new String[] {cmd.toString()}, bat, true, log);

    List<String> run = new ArrayList<>();
    run.add(bat);

    CmdLine.runCommandWithFileChecks(run, "", inputFiles, outputFiles, true, false, false, true,
                                     log);

  }

}
