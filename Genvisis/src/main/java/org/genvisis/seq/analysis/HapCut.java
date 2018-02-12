/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.concurrent.Callable;
import org.genvisis.CLI;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Wrapper for https://github.com/vibansal/HapCUT2 and utility methods for creating single sample
 * VCFs Paper: http://genome.cshlp.org/content/early/2016/12/09/gr.213462.116.abstract Basic
 * process: <br>
 * 1. Create single sample .vcf (required by HapCut) <br>
 * 2. Run extractHairs<br>
 * 3. Run HapCut<br>
 * 4. Convert hapcut output to vcf<br>
 */
public class HapCut {

  private HapCut() {

  }

  /**
   * @param vcf input vcf
   * @param inputDir where bams are located
   * @param outDir where results will be stored
   * @param variantSet if vcf is a result of a previous merge, the variant set to use
   * @param hapCutLoc location of HAPCUT2
   * @param extractHairsLoc location of extractHAIRS
   * @param fgbioLoc location of fgbio jar
   * @param threads you know, for speed
   */
  public static void run(String vcf, String inputDir, String outDir, String variantSet,
                         String hapCutLoc, String extractHairsLoc, String fgbioLoc, int threads) {
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "hapCut.log");

    String[] bams = Files.listFullPaths(inputDir, ".bam");
    log.reportTimeInfo("found " + bams.length + " bams to HapCut");
    HapCutProducer producer = new HapCutProducer(bams, vcf, outDir, variantSet, hapCutLoc,
                                                 extractHairsLoc, fgbioLoc, log);
    WorkerTrain<HapCutResult> train = new WorkerTrain<>(producer, threads, 10, log);
    while (train.hasNext()) {
      train.next();
    }
  }

  /**
   * for threading the analysis
   */
  private static class HapCutProducer implements Producer<HapCutResult> {

    private final String[] bams;
    private final String vcf;
    private final String outDir;
    private final String variantSet;
    private final String hapCutLoc;
    private final String extractHairsLoc;
    private final String fgbioLoc;

    private final Logger log;
    private int index;

    private HapCutProducer(String[] bams, String vcf, String outDir, String variantSet,
                           String hapCutLoc, String extractHairsLoc, String fgbioLoc, Logger log) {
      super();
      this.bams = bams;
      this.vcf = vcf;
      this.outDir = outDir;
      this.variantSet = variantSet;
      this.hapCutLoc = hapCutLoc;
      this.extractHairsLoc = extractHairsLoc;
      this.fgbioLoc = fgbioLoc;
      this.log = log;
      this.index = 0;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext() {
      return index < bams.length;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    @Override
    public Callable<HapCutResult> next() {
      if (index >= bams.length) {
        throw new NoSuchElementException("Ran out of bams while iterating");
      }
      final String bam = bams[index];
      index++;
      return () -> processBam(vcf, outDir, variantSet, hapCutLoc, extractHairsLoc, fgbioLoc, log,
                              bam);
    }

    /*
     * (non-Javadoc)
     * @see org.genvisis.common.WorkerTrain.Producer#shutdown()
     */
    @Override
    public void shutdown() {
      // not needed
    }

    /**
     * Runs HapCut2 pipeline on a bam
     */
    private static HapCutResult processBam(String vcf, String outDir, String variantSet,
                                           String hapCutLoc, String extractHairsLoc,
                                           String fgBioLoc, Logger log, String bam) {
      HapCutResult result = extractSampleForBam(bam, vcf, outDir, variantSet, log);
      if (result.success) {
        result.success = runExtractHairs(result, bam, extractHairsLoc, log);
        if (result.success) {
          result.success = runHapCut2(result, hapCutLoc, log);
          if (!result.success) {
            log.reportTimeWarning("failed to generate haplotype blocks for " + bam);
          }
          if (result.success) {
            result.success = runFGBIOConversion(result, fgBioLoc, log);
            if (!result.success) {
              log.reportTimeWarning("failed to convert block results to vcf for " + bam);
            }
          }
        } else {
          log.reportTimeWarning("failed to extract hairs for " + bam);
        }
      }
      return result;
    }

    private static boolean runFGBIOConversion(HapCutResult result, String fgBioLoc, Logger log) {

      Set<String> necessaryInputFiles = new HashSet<>();
      necessaryInputFiles.add(result.extractedVCF);
      necessaryInputFiles.add(fgBioLoc);
      necessaryInputFiles.add(result.blockResults);

      Set<String> expectedOutputFiles = new HashSet<>();
      expectedOutputFiles.add(result.fgBioVCF);
      List<String> commandList = new ArrayList<>();
      commandList.add("java");
      commandList.add("-jar");
      commandList.add(fgBioLoc);
      commandList.add("HapCutToVcf");
      commandList.add("-v");
      commandList.add(result.extractedVCF);
      commandList.add("-i");
      commandList.add(result.blockResults);
      commandList.add("-o");
      commandList.add(result.fgBioVCF);
      return CmdLine.runCommandWithFileChecks(commandList, "", necessaryInputFiles,
                                              expectedOutputFiles, true, false, false, log);

    }

    /**
     * See
     * https://github.com/vibansal/HapCUT2/blob/eb3b64b14b40cec348916fd3f03222ca2b023ba0/README.md
     * 
     * @param result {@link HapCutResult}
     * @param hapCutLoc location of HapCut2
     * @param log
     * @return true if successful
     */
    private static boolean runHapCut2(HapCutResult result, String hapCutLoc, Logger log) {
      Set<String> necessaryInputFiles = new HashSet<>();
      necessaryInputFiles.add(result.extractedVCF);
      necessaryInputFiles.add(hapCutLoc);
      necessaryInputFiles.add(result.hairResults);

      Set<String> expectedOutputFiles = new HashSet<>();
      expectedOutputFiles.add(result.blockResults);
      List<String> commandList = new ArrayList<>();
      commandList.add(hapCutLoc);
      commandList.add("--fragments");
      commandList.add(result.hairResults);
      commandList.add("--VCF");
      commandList.add(result.extractedVCF);
      commandList.add("--output");
      commandList.add(result.blockResults);
      return CmdLine.runCommandWithFileChecks(commandList, "", necessaryInputFiles,
                                              expectedOutputFiles, true, false, false, log);

    }

    /**
     * See
     * https://github.com/vibansal/HapCUT2/blob/eb3b64b14b40cec348916fd3f03222ca2b023ba0/README.md
     * 
     * @param result {@link HapCutResult}
     * @param bam bam file to use
     * @param extractHairsLoc location of extractHairs
     * @param log
     * @return true if successful
     */
    private static boolean runExtractHairs(HapCutResult result, String bam, String extractHairsLoc,
                                           Logger log) {
      Set<String> necessaryInputFiles = new HashSet<>();
      necessaryInputFiles.add(result.extractedVCF);
      necessaryInputFiles.add(bam);
      necessaryInputFiles.add(extractHairsLoc);
      Set<String> expectedOutputFiles = new HashSet<>();
      expectedOutputFiles.add(result.hairResults);
      List<String> commandList = new ArrayList<>();
      commandList.add(extractHairsLoc);
      commandList.add("--bam");
      commandList.add(bam);
      commandList.add("--VCF");
      commandList.add(result.extractedVCF);
      commandList.add("--out");
      commandList.add(result.hairResults);

      return CmdLine.runCommandWithFileChecks(commandList, "", necessaryInputFiles,
                                              expectedOutputFiles, true, false, false, log);

    }

    /**
     * Create a single sample vcf for the sample-matched bam
     */
    private static HapCutResult extractSampleForBam(String bam, String vcf, String outDir,
                                                    String variantSet, Logger log) {
      String sample = BamOps.getSampleName(bam, log) + variantSet;
      String outputVcf = outDir + sample + "_" + ext.removeDirectoryInfo(vcf);
      if (!Files.exists(outputVcf)) {
        try (VCFFileReader reader = new VCFFileReader(new File(vcf), false)) {

          if (!reader.getFileHeader().getSampleNameToOffset().containsKey(sample)) {
            log.reportTimeWarning("Could not find sample " + sample + " in vcf " + vcf
                                  + ", skipping (verify variant set if this is unexpected");
            return new HapCutResult(false, outputVcf);
          }
          VariantContextWriter writer = VCFOps.initWriter(outputVcf, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                          VCFOps.getSequenceDictionary(reader));
          Set<String> limit = new HashSet<>();
          limit.add(sample);
          // add a new header to the output, using only the single sample of interest
          VCFOps.copyHeader(reader, writer, limit, HEADER_COPY_TYPE.SUBSET_STRICT, log);

          for (VariantContext vc : reader) {
            writer.add(VCOps.getSubset(vc, limit, VC_SUBSET_TYPE.SUBSET_STRICT, false));
          }
          writer.close();
        }
      }
      return new HapCutResult(true, outputVcf);
    }
  }

  /**
   * Store for HapCut analysis, mainly file name formats
   */
  private static class HapCutResult {

    private boolean success;
    private final String extractedVCF;
    private final String hairResults;
    private final String blockResults;
    private final String fgBioVCF;

    private HapCutResult(boolean success, String extractedVCF) {
      super();
      this.success = success;
      this.extractedVCF = extractedVCF;
      this.hairResults = VCFOps.getAppropriateRoot(extractedVCF, false) + ".hairs.txt";
      this.blockResults = ext.rootOf(hairResults, false) + ".blocks.txt";
      this.fgBioVCF = ext.rootOf(blockResults, false) + ".vcf.gz";
    }
  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(HapCut.class);
    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF, "variants.vcf", true);
    c.addArgWithDefault("extractHAIRS",
                        "full path to extractHAIRS, see https://github.com/vibansal/HapCUT2",
                        "~/git/HapCUT2/build/extractHAIRS");
    c.addArgWithDefault("HAPCUT2", "full path to HAPCUT2, see https://github.com/vibansal/HapCUT2",
                        "~/git/HapCUT2/build/HAPCUT2");
    c.addArgWithDefault(CLI.ARG_INDIR,
                        "Input directory containing bam file, each sample must be represented in the vcf",
                        "./bams/");

    c.addArgWithDefault("variantSet", "Variant set to use, if your vcf utilizes variant sets", "");
    c.addArgWithDefault("fgbio",
                        "full path to fgbio jar file, see https://github.com/fulcrumgenomics/fgbio",
                        "~/git/fgbio/");

    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "./hapCutRuns/");
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "1");

    c.parseWithExit(args);

    run(c.get(CLI.ARG_VCF), c.get(CLI.ARG_INDIR), c.get(CLI.ARG_OUTDIR), c.get("variantSet"),
        c.get("HAPCUT2"), c.get("extractHAIRS"), c.get("fgbio"), c.getI(CLI.ARG_THREADS));
  }
}
