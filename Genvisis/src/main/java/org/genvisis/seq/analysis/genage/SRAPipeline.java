package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.Map;
import java.util.concurrent.Callable;

import org.apache.commons.cli.Options;
import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.sra.SRARunTable;
import org.genvisis.sra.SRASample;
import org.genvisis.sra.SRAUtils;
import org.genvisis.sra.SRAUtils.SRABamWorker;
import org.genvisis.sra.SRAUtils.SRAConversionResult;

/**
 * more specific version of {@link Pipeline} that starts with a single SRA file
 *
 */
public class SRAPipeline implements Callable<Boolean> {

  private final SRASample sraSample;
  private final String inputSRA;
  private final String rootOutDir;
  private final String referenceGenome;
  private final String captureBed;
  private final int numThreads;
  private final Logger log;

  /**
   * @param sraSample an {@link SRASample} to analyze
   * @param inputSRA the input sra file in appropriate sra-toolkit directory
   * @param rootOutDir the output directory for the analysis
   * @param referenceGenome proper reference genome
   * @param captureBed the capture bed, only utilized with {@link ASSAY_TYPE#WXS}
   * @param atype {@link ASSAY_TYPE} of the sample
   * @param aName {@link ASSEMBLY_NAME} for the sample
   * @param numThreads number of threads for the pipeline branches
   * @param log
   */
  public SRAPipeline(SRASample sraSample, String inputSRA, String rootOutDir,
                     String referenceGenome, String captureBed, int numThreads, Logger log) {
    super();
    this.sraSample = sraSample;
    this.inputSRA = inputSRA;
    this.rootOutDir = rootOutDir;
    this.referenceGenome = referenceGenome;
    this.captureBed = captureBed;
    this.numThreads = numThreads;
    this.log = log;
  }

  @Override
  public Boolean call() throws Exception {
    String bamDir = rootOutDir + "bams/";
    new File(bamDir).mkdirs();
    String bam = bamDir + ext.rootOf(inputSRA) + ".bam";
    WorkerHive<SRAConversionResult> hive = new WorkerHive<SRAUtils.SRAConversionResult>(1, 10, log);
    hive.addCallable(new SRABamWorker(inputSRA, bam, log));
    hive.execute(true);
    Pipeline.pipeline(bam, rootOutDir, referenceGenome, captureBed, sraSample, numThreads, log);

    return true;
  }

  /**
   * This will be a bit "reversed" in the final pipeline version...This method is for processing
   * many pre-downloaded files
   */
  private static void runAll(String sraDir, String sraRunTableFile, String rootOutDir,
                             String referenceGenome, String captureBed, int numThreads) {
    Logger log = new Logger();
    String[] sraFiles = Files.list(sraDir, ".sra", false);
    SRARunTable srRunTable = SRARunTable.load(sraRunTableFile, log);
    log.reportTimeInfo("Found " + sraFiles.length + " sra files in " + sraDir);
    WorkerHive<Boolean> hive = new WorkerHive<Boolean>(numThreads, 10, log);
    for (String sraFile : sraFiles) {
      SRASample sample = srRunTable.get(ext.rootOf(sraFile));

      SRAPipeline pipeline = new SRAPipeline(sample, sraFile, rootOutDir, referenceGenome,
                                             captureBed, 1, log);
      hive.addCallable(pipeline);
    }

    hive.execute(true);
  }

  public static void main(String[] args) {

    Options options = CLI.defaultOptions();

    String sraDirDefault = "sra/";
    final String SRA_DRI = "sraDir";
    CLI.addArg(options, SRA_DRI, "directory with .sra files", sraDirDefault);

    String outDir = "out/";
    final String OUT_DIR = "outDir";
    CLI.addArg(options, OUT_DIR, "the output directory for results", outDir);

    String sraRunTableDefault = "sraRuntable.txt";
    final String SRA_RUN_TABLE = "sraRunTable";
    CLI.addArg(options, SRA_RUN_TABLE, "a sra run table providing sample information",
               sraRunTableDefault);

    int numThreads = 24;
    final String NUM_THREADS = "threads";
    CLI.addArg(options, NUM_THREADS, "a sra run table providing sample information",
               Integer.toString(numThreads));

    String refGenomeFasta = "hg19.canonical.fa";
    final String REFERENC_GENOME = "ref";
    CLI.addArg(options, REFERENC_GENOME, "appropriate reference genome file", refGenomeFasta);

    String captureBedFile = "VCRome_2_1_hg19_capture_targets.bed";
    final String CAPTURE_BED = "bed";
    CLI.addArg(options, CAPTURE_BED, "bed file of targeted capture", captureBedFile);

    Map<String, String> parsed = CLI.parseWithExit(SRAPipeline.class, options, args);

    runAll(parsed.get(SRA_DRI), parsed.get(SRA_RUN_TABLE), parsed.get(OUT_DIR),
           parsed.get(REFERENC_GENOME), parsed.get(CAPTURE_BED), numThreads);

  }

}
