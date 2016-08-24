package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BamImport;
import org.genvisis.seq.manage.ReferenceGenome;
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

  private SRASample sraSample;
  private String inputSRA;
  private String rootOutDir;
  private String referenceGenome;
  private String captureBed;
  private String binBed;
  private String vcf;
  private int numThreads;
  private Logger log;

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
                     String referenceGenome, String captureBed, String binBed, String vcf,
                     int numThreads, Logger log) {
    super();
    this.sraSample = sraSample;
    this.inputSRA = inputSRA;
    this.rootOutDir = rootOutDir;
    this.referenceGenome = referenceGenome;
    this.captureBed = captureBed;
    this.binBed = binBed;
    this.vcf = vcf;
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
    Pipeline.pipeline(bam, rootOutDir, referenceGenome, captureBed, binBed, vcf, sraSample,
                      numThreads, log);

    return true;
  }

  /**
   * This will be a bit "reversed" in the final pipeline version...This method is for processing
   * many pre-downloaded files
   */
  private static void runAll(String sraInput, String sraRunTableFile, String rootOutDir,
                             String referenceGenome, String captureBed, String binBed, String vcf,
                             int numThreads) {
    Logger log = new Logger();
    String[] sraFiles = null;
    if (Files.isDirectory(sraInput)) {
      log.reportTimeInfo("Gathering sra files from " + sraInput);
      sraFiles = Files.list(sraInput, ".sra", false);
    } else {
      log.reportTimeInfo("Reading sra files from " + sraInput);
      sraFiles = HashVec.loadFileToStringArray(sraInput, false, new int[] {0}, true);
    }
    SRARunTable srRunTable = SRARunTable.load(sraRunTableFile, log);
    log.reportTimeInfo("Found " + sraFiles.length + " sra files in " + sraInput);
    WorkerHive<Boolean> hive = new WorkerHive<Boolean>(numThreads, 10, log);
    boolean prelimGenvisisWGS = false;
    boolean prelimGenvisisWXS = false;
    ArrayList<String> sampleSummary = new ArrayList<String>();
    for (String sraFile : sraFiles) {
      SRASample sample = srRunTable.get(ext.rootOf(sraFile));
      sampleSummary.add(sraFile + "\t" + sample.toString());
      SRAPipeline pipeline = new SRAPipeline(sample, sraFile, rootOutDir, referenceGenome,
                                             captureBed, binBed, vcf, 1, log);
      switch (sample.getaType()) {// create the required markerSets for import...prior to threading
        case WGS:
          if (!prelimGenvisisWGS) {
            Project proj = Pipeline.getProjectFor(sample.getaType(), rootOutDir, referenceGenome);

            BamImport.generateAnalysisSet(proj, null, null, vcf, BamImport.CAPTURE_BUFFER,
                                          sample.getaType(), log,
                                          new ReferenceGenome(referenceGenome, log));
            prelimGenvisisWGS = true;
          }
          break;
        case WXS:
          if (!prelimGenvisisWXS) {
            Project proj = Pipeline.getProjectFor(sample.getaType(), rootOutDir, referenceGenome);
            BamImport.generateAnalysisSet(proj, binBed, captureBed, vcf, BamImport.CAPTURE_BUFFER,
                                          sample.getaType(), log,
                                          new ReferenceGenome(referenceGenome, log));
            prelimGenvisisWXS = true;
          }
          break;
        default:
          break;

      }
      hive.addCallable(pipeline);
    }
    Files.writeArrayList(sampleSummary, rootOutDir + "sampleAnalysis.summary.txt");
    System.exit(1);
    hive.execute(true);
  }



  public static void main(String[] args) {

    CLI c = new CLI();

    String sraDirDefault = "sra/";
    final String SRA_Input = "sraInput";
    c.addArg(SRA_Input, "directory or filename with .sra files", sraDirDefault);

    String outDir = "out/";
    final String OUT_DIR = "outDir";
    c.addArg(OUT_DIR, "the output directory for results", outDir);

    String sraRunTableDefault = "sraRuntable.txt";
    final String SRA_RUN_TABLE = "sraRunTable";
    c.addArg(SRA_RUN_TABLE, "a sra run table providing sample information", sraRunTableDefault);

    int numThreads = 24;
    final String NUM_THREADS = "threads";
    c.addArg(NUM_THREADS, "a sra run table providing sample information",
             Integer.toString(numThreads));

    String refGenomeFasta = "hg19.canonical.fa";
    final String REFERENCE_GENOME = "ref";
    c.addArg(REFERENCE_GENOME, "appropriate reference genome file", refGenomeFasta);

    String captureBedFile = "VCRome_2_1_hg19_capture_targets.bed";
    final String CAPTURE_BED = "bed";
    c.addArg(CAPTURE_BED, "bed file of targeted capture", captureBedFile);

    String binBed = "targetsOfInterest.bed";
    final String BIN_BED = "bin";
    c.addArg(BIN_BED, "bed file of targets of interests", binBed);

    String vcf = "vcf.vcf";
    final String VCF = "vcf";
    c.addArg(VCF, "vcf file of variants", vcf);


    c.parseWithExit(SRAPipeline.class, args);

    runAll(c.get(SRA_Input), c.get(SRA_RUN_TABLE), c.get(OUT_DIR), c.get(REFERENCE_GENOME),
           c.get(CAPTURE_BED), c.get(BIN_BED), c.get(VCF), numThreads);

  }

}
