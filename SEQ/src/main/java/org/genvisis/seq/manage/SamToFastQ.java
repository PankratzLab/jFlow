/**
 * 
 */
package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.WorkerTrain.Producer;
import org.pankratzlab.core.CLI;
import org.pankratzlab.shared.qsub.Qsub;

/**
 * Prepares bam input for re-genotyping by handling bam-> fastq conversion.
 */
public class SamToFastQ {

  private SamToFastQ() {

  }

  /**
   * @param inputBam bam to convert
   * @param samToFastQLoc path to picard.jar, or path to SamToFastq.jar
   * @param r1 output fastq for
   * @param r2
   * @param log
   * @return
   */
  private static boolean convertToFasta(String inputBam, String samToFastQLoc, String r1, String r2,
                                        int memoryInMb, Logger log) {
    String[] inputs = new String[] {inputBam};
    String[] outputs = new String[] {r1, r2};
    ArrayList<String> command = new ArrayList<>();
    command.add("java");
    command.add("-Xmx" + memoryInMb + "m");
    command.add("-jar");

    command.add(samToFastQLoc);
    if (samToFastQLoc.endsWith("picard.jar")) {
      command.add("SamToFastq");
    }
    command.add("I=" + inputBam);
    command.add("F=" + r1);
    command.add("F2=" + r2);

    return CmdLine.runCommandWithFileChecks(ArrayUtils.toStringArray(command), "", inputs, outputs,
                                            true, false, false, log);
  }

  private static void prepBams(String bams, final String outDir, final String tag,
                               final String samToFastQ, int numThreads, int memoryInMb) {
    new File(outDir).mkdirs();
    final Logger log = new Logger(outDir + "samToFastq.log");
    final String[] bamFiles = getBams(bams);
    log.reportTimeInfo("Found " + bamFiles.length + " bams from input " + bams);

    Producer<Boolean> prepProducer = new Producer<Boolean>() {

      private int index = 0;

      @Override
      public boolean hasNext() {

        return index < bamFiles.length;
      }

      @Override
      public Callable<Boolean> next() {
        final String bamFile = bamFiles[index];

        Callable<Boolean> callable = new Callable<Boolean>() {

          @Override
          public Boolean call() throws Exception {
            String sampleName = null;
            try {
              if (BamOps.getHeader(bamFile, log).getReadGroups().size() != 1) {
                throw new IllegalArgumentException("This method currently supports bam to .fastq conversion for 1 and only 1 PE readgroups");
              }
              sampleName = BamOps.getSampleName(bamFile, log);
              String rootOut = outDir + (tag == null ? "" : tag + "-") + sampleName
                               + "_UUUUU-UUUUU_L001_.fastq";
              String r1 = ext.addToRoot(rootOut, "R1_001") + ".gz";
              String r2 = ext.addToRoot(rootOut, "R2_001") + ".gz";
              boolean success = convertToFasta(bamFile, samToFastQ, r1, r2, memoryInMb, log);

              if (!success) {
                log.reportError("Could not parse " + bamFile + ", removing any output");
                new File(r1).delete();
                new File(r2).delete();
              }
              return success;
            } catch (Exception e) {
              log.reportError("Could not process " + bamFile);
              log.reportException(e);
            }
            return false;
          }
        };
        index++;
        return callable;
      }

      @Override
      public void shutdown() {
        //
      }

      @Override
      public void remove() {
        //
      }
    };
    try (WorkerTrain<Boolean> train = new WorkerTrain<>(prepProducer, numThreads, 10, log)) {
      while (train.hasNext()) {
        train.next();
      }
    }
  }

  private static String[] getBams(String bams) {
    return Files.isDirectory(bams) ? Files.listFullPaths(bams, ".bam")
                                   : HashVec.loadFileToStringArray(bams, false, new int[] {0},
                                                                   true);
  }

  /**
   * 
   */
  private static final String GENVISIS = "genvisis";
  /**
   * 
   */
  private static final String BATCH = "batch";
  /**
   * 
   */
  private static final String TAG = "tag";
  /**
   * 
   */
  private static final String SAM_TO_FASTQ = "samToFastq";
  /**
   * 
   */
  private static final String BAMS = "bams";

  private static void batch(CLI c) {
    new File(c.get(CLI.ARG_OUTDIR)).mkdirs();
    final String[] bamFiles = getBams(c.get(BAMS));
    Logger log = new Logger(c.get(CLI.ARG_OUTDIR) + "batch.log");
    log.reportTimeInfo("Found " + bamFiles.length + " bams from input " + c.get(BAMS));
    String[][] batches = ArrayUtils.splitUpStringArray(bamFiles, c.getI(BATCH), log);
    List<String> cmd = new ArrayList<>();
    cmd.add("java");
    cmd.add("-jar");
    cmd.add(c.get(GENVISIS));
    cmd.add("seq.manage.SamToFastQ");
    cmd.add(SAM_TO_FASTQ + "=" + c.get(SAM_TO_FASTQ));
    if (c.get(TAG) != null) {
      cmd.add(TAG + "=" + c.get(TAG));
    }
    cmd.add(BATCH + "=-1");
    cmd.add(CLI.ARG_OUTDIR + "=" + c.get(CLI.ARG_OUTDIR));

    cmd.add(CLI.ARG_THREADS + "=" + c.getI(CLI.ARG_THREADS));

    for (int i = 0; i < batches.length; i++) {
      String batchFile = c.get(CLI.ARG_OUTDIR) + "samToFastQ" + i + ".txt";
      String pbs = c.get(CLI.ARG_OUTDIR) + "samToFastQ" + i + ".pbs";
      Files.writeArray(batches[i], batchFile);
      List<String> currentCmd = new ArrayList<>();
      currentCmd.addAll(cmd);
      currentCmd.add(BAMS + "=" + batchFile);

      Qsub.qsub(pbs, ArrayUtils.toStr(currentCmd, " "), c.getI(PSF.Ext.MEMORY_MB),
                c.getI(PSF.Ext.WALLTIME_HRS), c.getI(CLI.ARG_THREADS));

    }

  }

  public static void main(String[] args) {
    CLI c = new CLI(SamToFastQ.class);
    c.addArgWithDefault(BAMS,
                        "file listing bam files to analyze, one per line - or a directory of bams",
                        "bams.txt");
    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "out/");
    c.addArgWithDefault(SAM_TO_FASTQ, "full path to SamToFastq.jar", "SamToFastq.jar");
    c.addArgWithDefault(TAG, "custom ID tag to add to files", null);
    c.addArgWithDefault(BATCH,
                        "if batching is desired (set to >0 for no batching), number of batches",
                        -1);
    c.addArgWithDefault(GENVISIS, "location of genvisis.jar", "~/genvisis.jar");
    c.addArgWithDefault(PSF.Ext.MEMORY_MB, "memory in mb if batching", PSF.Ext.DEFAULT_MEMORY_MB);
    c.addArgWithDefault(PSF.Ext.WALLTIME_HRS, "walltime in hours, if batching", 48);

    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, 1);
    c.parseWithExit(args);

    if (c.getI(BATCH) > 0) {
      batch(c);
    } else {
      prepBams(c.get(BAMS), c.get(CLI.ARG_OUTDIR), c.get(TAG), c.get(SAM_TO_FASTQ),
               c.getI(CLI.ARG_THREADS), c.getI(PSF.Ext.MEMORY_MB));
    }

  }

}
