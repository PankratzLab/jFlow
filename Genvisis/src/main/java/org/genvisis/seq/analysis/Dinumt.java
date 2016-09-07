/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.analysis.genage.SRAPipeline;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BamOps.InsertSizeEstimate;
import org.genvisis.seq.qc.BamQC;
import org.genvisis.seq.qc.FilterNGS;

/**
 * @author Kitty
 *
 *         Wrapper for https://github.com/mills-lab/dinumt
 * 
 *         May become part of the {@link SRAPipeline}
 */
public class Dinumt {



  private static double getAvgCoverage(String bamQCOutput) {
    String[] dataToLoad = new String[] {};

    return 2;
  }

  private static void runQC(String inputBam, String outputDir, String targetLibraryFile,
                            int numThreads, Logger log) {
    String bamQCDir = outputDir + BamOps.getSampleName(inputBam, log) + "/";
    new File(bamQCDir).mkdirs();
    String bamList = bamQCDir + ".bamList.txt";
    String bamQCOutput = bamQCDir + ".bamQC.txt";
    Files.write(inputBam, bamList);
    double mappingQuality = 0;
    double phreadScore = 0.0;
    int[] readDepth = {0, 1, 2, 3, 4, 10, 20, 30, 40};
    FilterNGS filterNGS = new FilterNGS(mappingQuality, phreadScore, readDepth);

    BamQC.qcBams(null, outputDir, null, bamList, targetLibraryFile, null, 2, filterNGS, numThreads,
                 bamQCOutput, null, false, 0, log);
    InsertSizeEstimate insertSizeEstimate =
                                          BamOps.estimateInsertSize(inputBam,
                                                                    BamOps.NUM_READ_ESTIMATOR, log);


  }

  public static void main(String[] args) {
    String bam = "/Volumes/Beta/data/aric_sra/bams/SRR1654226.bam";
    String targetLibraryFile = "/Volumes/Beta/ref/VCRome_2_1_hg19_capture_targets.bed";
    String outputDir = "/Volumes/Beta/data/aric_sra/testDinumt/";

    runQC(bam, outputDir, targetLibraryFile, 1, new Logger());
  }

}
