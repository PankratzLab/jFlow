/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.genage.SRAPipeline;
import org.genvisis.seq.manage.BamOps;
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



  private static QCParams getQCParams(String bamQCOutput, String bamFile, Logger log) {
    String[] dataToLoad = new String[] {"numOnTarget", "Total Base Pairs Targeted",
                                        "AverageOnTargetInsertSize", "OnTargetInsertSizeStdev"};
    int[] indices = ext.indexFactors(dataToLoad, Files.getHeaderOfFile(bamQCOutput, log), true,
                                     false);


    String[] stats =
                   HashVec.loadFileToStringArray(bamQCOutput, true, indices, false)[0].split("\t");
    double bpCoverage = Double.parseDouble(stats[0]) * BamOps.estimateReadSize(bamFile, log);
    double targeted = Double.parseDouble(stats[1]);
    double averageCoverage = bpCoverage / targeted;
    return new QCParams(Double.parseDouble(stats[2]), Double.parseDouble(stats[3]),
                        averageCoverage);
  }

  private static String[] developDinumtCommand(String inputBam, String fullPathToDinumt,
                                               String refNumts, String referenceGenome,
                                               String outputDir, Logger log) {
    ArrayList<String> command = new ArrayList<String>();
    String sampleName = BamOps.getSampleName(inputBam, log);
    String outputFile = outputDir + sampleName + ".vcf";
    command.add(fullPathToDinumt);
    command.add("--mask_filename=" + refNumts);
    command.add("--input_filename=" + inputBam);

    command.add("--reference=" + referenceGenome);
    command.add("--min_reads_cluster=1");
    command.add("-include_mask");
    command.add("--output_filename=" + outputFile);
    command.add("--prefix=" + sampleName);
    command.add("--output_filename=" + outputFile);
    command.add("--output_filename=" + outputFile);


    // --len_cluster_include=577 \
    // --len_cluster_link=1154 \
    // --insert_size=334.844984 \
    // --max_read_cov=29 \
    // --output_support \
    // --support_filename=sample1_support.sam

    return Array.toStringArray(command);
  }

  private static class QCParams {
    private double avgInsertSize;
    private double stDevInsertSize;

    private double avgCoverage;

    public QCParams(double avgInsertSize, double stDevInsertSize, double avgCoverage) {
      super();
      this.avgInsertSize = avgInsertSize;
      this.stDevInsertSize = stDevInsertSize;
      this.avgCoverage = avgCoverage;
    }


  }

  private static QCParams runQC(String inputBam, String outputDir, String targetLibraryFile,
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

      BamQC.qcBams(null, bamQCDir, null, bamList, targetLibraryFile, null, 2, filterNGS, numThreads,
                   bamQCOutput, null, false, 0, false, log);
    }

    return getQCParams(bamQCDir + bamQCOutput, inputBam, log);

  }

  public static void main(String[] args) {
    String bam = "/Volumes/Beta/data/aric_sra/bams/SRR1654226.bam";
    String targetLibraryFile = "/Volumes/Beta/ref/VCRome_2_1_hg19_capture_targets.bed";
    String outputDir = "/Volumes/Beta/data/aric_sra/testDinumt/";

    runQC(bam, outputDir, targetLibraryFile, 1, new Logger());
  }

}
