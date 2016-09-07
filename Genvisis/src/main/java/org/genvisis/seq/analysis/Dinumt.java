/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
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



  private static double getAvgCoverage(String bamQCOutput, String bamFile, Logger log) {
    String[] dataToLoad = new String[] {"numOnTarget", "basepairsTargeted"};
    int[] indices = ext.indexFactors(dataToLoad, Files.getHeaderOfFile(bamQCOutput, log), true,
                                     false);


    String[] stats =
                   HashVec.loadFileToStringArray(bamQCOutput, true, indices, false)[0].split("\t");
    double bpCoverage = Double.parseDouble(stats[0]) * BamOps.estimateReadSize(bamFile, log);
    double targeted = Double.parseDouble(stats[1]);

    return bpCoverage / targeted;
  }

  private static class QCParams {
    private InsertSizeEstimate insertSizeEstimate;
    private double avgCoverage;

    public QCParams(InsertSizeEstimate insertSizeEstimate, double avgCoverage) {
      super();
      this.insertSizeEstimate = insertSizeEstimate;
      this.avgCoverage = avgCoverage;
    }


  }

  private static void runQC(String inputBam, String outputDir, String targetLibraryFile,
                            int numThreads, Logger log) {
    String sampleName = BamOps.getSampleName(inputBam, log);
    String bamQCDir = outputDir + sampleName + "/";
    new File(bamQCDir).mkdirs();
    String bamList = bamQCDir + sampleName + ".bamList.txt";
    String bamQCOutput = sampleName + ".bamQC.txt";
    if (!Files.exists(bamQCDir + bamQCOutput)) {
      Files.write(inputBam, bamList);
      double mappingQuality = 0;
      double phreadScore = 0.0;
      int[] readDepth = {0, 1, 2, 3, 4, 10, 20, 30, 40};
      FilterNGS filterNGS = new FilterNGS(mappingQuality, phreadScore, readDepth);

      BamQC.qcBams(null, outputDir, null, bamList, targetLibraryFile, null, 2, filterNGS,
                   numThreads, bamQCOutput, null, false, 0, log);
    }
    String insertSizeFile = bamQCDir + "insertSize.txt";
    InsertSizeEstimate insertSizeEstimate;
    if (!Files.exists(insertSizeFile)) {
      insertSizeEstimate = BamOps.estimateInsertSize(inputBam, BamOps.NUM_READ_ESTIMATOR, log);
      Files.write(insertSizeEstimate.getSummary(), insertSizeFile);
    } else {
      insertSizeEstimate = InsertSizeEstimate.load(insertSizeFile);
    }

  }

  public static void main(String[] args) {
    String bam = "/Volumes/Beta/data/aric_sra/bams/SRR1654226.bam";
    String targetLibraryFile = "/Volumes/Beta/ref/VCRome_2_1_hg19_capture_targets.bed";
    String outputDir = "/Volumes/Beta/data/aric_sra/testDinumt/";

    runQC(bam, outputDir, targetLibraryFile, 1, new Logger());
  }

}
