package org.genvisis.one.JL.beta;

import java.io.File;
import org.genvisis.cnv.analysis.pca.BetaOptimizer;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ext;

/**
 * @author Kitty Fixing beta optimaization...hopefully
 */
public class BetaOpDev {

  private BetaOpDev() {

  }

  private static void test() {
    Project proj = new Project("/Users/Kitty/workspace.other/Genvisis/Genvisis/projects/LLFS.properties");
    String pcFile = "/Volumes/Beta/data/LLFS/Final_ohw_ws_40_ALL1000PCs_gc_corrected_OnTheFly_SampLRR_Recomp_LRR_030CR_096.PCs.extrapolated.txt";
    String pcSamps = "/Volumes/Beta/data/LLFS/Final_ohw_ws_40_ALL1000PCs_gc_corrected_OnTheFly_SampLRR_Recomp_LRR_030CR_096.samples.USED_PC.txt";
    double[] pvals = new double[] {.001};
    proj.BLAST_ANNOTATION_FILENAME.setValue("/Users/Kitty/temp/LLFS/blast.vcfno.gz");// SSD
    String outDir = ext.parseDirectoryOfFile(pcFile) + "betaOpt2/";
    new File(outDir).mkdirs();

    new File(outDir + "WBC_SE_beta_summary.txt").delete();

    String unrelatedWhites = "/Volumes/Beta/data/LLFS/unrelateWParse/unrelatedWStudySample.txt";

    BetaOptimizer.optimize(proj, pcFile, outDir,
                           "/Users/Kitty/workspace.other/Genvisis/Genvisis/target/resources/MitoCN/",
                           unrelatedWhites, pcSamps, pvals, 100, .98, 10, 0.0000000001, 3);
    // 0.000001
  }

  public static void main(String[] args) {

    test();
  }

}
