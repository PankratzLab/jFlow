package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.analysis.genage.Pipeline.PipelinePart;
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
public class SRAPipeline implements Callable<List<PipelinePart>> {
  private static final String SRA_INPUT = "sraInput";
  private static final String OUT_DIR = "outDir";
  private static final String SRA_RUN_TABLE = "sraRunTable";
  private static final String NUM_THREADS = "threads";
  private static final String NUM_THREADS_PIPELINE = "threadsPipe";

  private static final String REFERENCE_GENOME = "ref";
  private static final String CAPTURE_BED = "bed";
  private static final String BIN_BED = "bin";
  private static final String VCF = "vcf";
  private static final String NUM_BATCHES = "batch";
  private static final String COMPILE = "compile";
  private static final String COMPUTEL = "computel";

  private SRASample sraSample;
  private String inputSRA;
  private String rootOutDir;
  private String referenceGenome;
  private String captureBed;
  private String binBed;
  private String vcfFile;
  private String computelLocation;
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
                     String computelLocation, int numThreads, Logger log) {
    super();
    this.sraSample = sraSample;
    this.inputSRA = inputSRA;
    this.rootOutDir = rootOutDir;
    this.referenceGenome = referenceGenome;
    this.captureBed = captureBed;
    this.binBed = binBed;
    this.vcfFile = vcf;
    this.computelLocation = computelLocation;
    this.numThreads = numThreads;
    this.log = log;
  }

  @Override
  public List<PipelinePart> call() throws Exception {
    String bamDir = getBamDirectory(rootOutDir);
    new File(bamDir).mkdirs();
    String bam = bamDir + ext.rootOf(inputSRA) + ".bam";
    WorkerHive<SRAConversionResult> hive = new WorkerHive<SRAUtils.SRAConversionResult>(1, 10, log);
    hive.addCallable(new SRABamWorker(inputSRA, bam, log));
    hive.execute(true);
    return Pipeline.pipeline(bam, rootOutDir, referenceGenome, captureBed, binBed, vcfFile,
                             sraSample, computelLocation, numThreads, log);
  }

  private static String getBamDirectory(String rootOutDir) {
    return rootOutDir + "bams/";
  }


  private static List<SRASample> loadSraSamples(String sraInput, String sraRunTableFile,
                                                Logger log) {
    ArrayList<SRASample> samples = new ArrayList<SRASample>();
    SRARunTable srRunTable = SRARunTable.load(sraRunTableFile, log);

    String[] sraFiles;
    if (Files.isDirectory(sraInput)) {
      log.reportTimeInfo("Gathering sra files from " + sraInput);
      sraFiles = Files.listFullPaths(sraInput, ".sra", false);
    } else {
      log.reportTimeInfo("Reading sra files from " + sraInput);
      sraFiles = HashVec.loadFileToStringArray(sraInput, false, new int[] {0}, true);
    }
    log.reportTimeInfo("Found " + sraFiles.length + " sra files in " + sraInput);
    for (String sraFile : sraFiles) {
      SRASample sample = srRunTable.get(ext.rootOf(sraFile));
      sample.setSraFile(sraFile);
      samples.add(sample);
    }


    return samples;
  }

  private static String[] getAssociatedBams(List<SRASample> samples, String rootOutDir,
                                            Logger log) {
    ArrayList<String> bams = new ArrayList<String>();
    for (SRASample sample : samples) {
      String bam = getBamDirectory(rootOutDir) + ext.rootOf(sample.getSraFile()) + ".bam";
      if (!Files.exists(bam)) {
        log.reportTimeWarning("Associated bam " + bam + " for " + sample.getSraFile()
                              + " did not exist, will not compile...run the import pipeline again");
      } else {
        bams.add(bam);
      }
    }
    return Array.toStringArray(bams);

  }

  private static void runCompile(List<SRASample> samples, String rootOutDir, String referenceGenome,
                                 String captureBed, String binBed, String vcf, int numThreads,
                                 Logger log) {

    String[] bams = getAssociatedBams(samples, rootOutDir, log);
    ASSAY_TYPE atType = samples.get(0).getaType();
    ASSEMBLY_NAME aName = samples.get(0).getaName();

    for (SRASample sample : samples) {
      if (sample.getaType() != atType) {
        throw new IllegalArgumentException("Mismatched assay types");
      }
      if (sample.getaName() != aName) {
        throw new IllegalArgumentException("Mismatched assembly name");
      }

    }
    Project proj = Pipeline.getProjectFor(atType, rootOutDir, referenceGenome);

    BamImport.importTheWholeBamProject(proj, binBed, captureBed, vcf, BamImport.CAPTURE_BUFFER, 4,
                                       true, atType, aName, bams, numThreads);
  }

  private static void compile(String sraInput, String sraRunTableFile, String rootOutDir,
                              String referenceGenome, String captureBed, String binBed, String vcf,
                              int numThreads) {
    Logger log = new Logger(rootOutDir + "compile.log");

    List<SRASample> samples = loadSraSamples(sraInput, sraRunTableFile, log);
    ArrayList<SRASample> wgsSamples = new ArrayList<SRASample>();
    ArrayList<SRASample> wxsSamples = new ArrayList<SRASample>();
    for (SRASample sample : samples) {
      switch (sample.getaType()) {
        case WGS:
          wgsSamples.add(sample);
          break;
        case WXS:
          wxsSamples.add(sample);
          break;
        default:
          throw new IllegalArgumentException("Invalid assay type " + sample.getaType());
      }
    }

    log.reportTimeInfo("Found " + wgsSamples.size() + " " + ASSAY_TYPE.WGS + " samples and "
                       + wxsSamples.size() + " " + ASSAY_TYPE.WXS + " samples");
    runCompile(wxsSamples, rootOutDir, referenceGenome, captureBed, binBed, vcf, numThreads, log);
    runCompile(wgsSamples, rootOutDir, referenceGenome, captureBed, binBed, vcf, numThreads, log);

  }


  /**
   * This will be a bit "reversed" in the final pipeline version...This method is for processing
   * many pre-downloaded files
   */
  private static void runAll(String sraInput, String sraRunTableFile, String rootOutDir,
                             String referenceGenome, String captureBed, String binBed, String vcf,
                             String computelLocation, int numThreads, int numThreadsPipeline,
                             int numBatches, CLI c) {
    Logger log = new Logger();

    WorkerHive<List<PipelinePart>> hive = new WorkerHive<List<PipelinePart>>(numThreads, 10, log);
    boolean prelimGenvisisWGS = false;
    boolean prelimGenvisisWXS = false;
    List<SRASample> samples = loadSraSamples(sraInput, sraRunTableFile, log);
    ArrayList<String> sampleSummary = new ArrayList<String>();
    ArrayList<String> sraFiles = new ArrayList<String>();

    for (SRASample sample : samples) {
      sraFiles.add(sample.getSraFile());
      sampleSummary.add(sample.getSraFile() + "\t" + sample.toString());
      SRAPipeline pipeline = new SRAPipeline(sample, sample.getSraFile(), rootOutDir,
                                             referenceGenome, captureBed, binBed, vcf,
                                             computelLocation, numThreadsPipeline, log);
      switch (sample.getaType()) {// create the required markerSets for import...prior to threading
        case WGS:
          if (!prelimGenvisisWGS) {
            Project proj = Pipeline.getProjectFor(sample.getaType(), rootOutDir, referenceGenome);
            if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {
              BamImport.generateAnalysisSet(proj, null, null, vcf, BamImport.CAPTURE_BUFFER,
                                            sample.getaType(), log,
                                            new ReferenceGenome(referenceGenome, log));
            }
            prelimGenvisisWGS = true;
          }

          break;
        case WXS:
          if (!prelimGenvisisWXS) {

            Project proj = Pipeline.getProjectFor(sample.getaType(), rootOutDir, referenceGenome);
            if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {

              BamImport.generateAnalysisSet(proj, binBed, captureBed, vcf, BamImport.CAPTURE_BUFFER,
                                            sample.getaType(), log,
                                            new ReferenceGenome(referenceGenome, log));
            }
            prelimGenvisisWXS = true;
          }
          break;
        default:
          throw new IllegalArgumentException("Invalid assay type " + sample.getaType());

      }
      hive.addCallable(pipeline);
    }
    Files.writeArrayList(sampleSummary, rootOutDir + "sampleAnalysis.summary.txt");
    if (numBatches > 0) {
      batch(Array.toStringArray(sraFiles), rootOutDir, c, log);
    } else {
      hive.execute(true);
    }
  }



  private static void batch(String[] sraFiles, String rootOutDir, CLI c, Logger log) {
    String[][] splits = Array.splitUpStringArray(sraFiles, c.getI(NUM_BATCHES), log);
    ArrayList<String> baseCommand = new ArrayList<String>();
    baseCommand.add("module load gcc/4.8.1\n");
    baseCommand.add("jcp seq.analysis.genage.SRAPipeline");
    baseCommand.add(OUT_DIR + "=" + c.get(OUT_DIR));
    baseCommand.add(SRA_RUN_TABLE + "=" + c.get(SRA_RUN_TABLE));
    baseCommand.add(NUM_THREADS + "=" + c.get(NUM_THREADS));
    baseCommand.add(NUM_THREADS_PIPELINE + "=" + c.get(NUM_THREADS_PIPELINE));
    baseCommand.add(REFERENCE_GENOME + "=" + c.get(REFERENCE_GENOME));
    baseCommand.add(CAPTURE_BED + "=" + c.get(CAPTURE_BED));
    baseCommand.add(BIN_BED + "=" + c.get(BIN_BED));
    baseCommand.add(VCF + "=" + c.get(VCF));
    baseCommand.add(COMPUTEL + "=" + c.get(COMPUTEL));

    String batchDir = rootOutDir + "batches/";
    new File(batchDir).mkdirs();
    for (int i = 0; i < splits.length; i++) {
      String batch = batchDir + "batch_" + i + ".txt";
      String qsub = batchDir + "batch_" + i + ".qsub";
      Files.writeList(splits[i], batch);
      ArrayList<String> currentCommand = new ArrayList<String>();
      currentCommand.addAll(baseCommand);
      currentCommand.add(SRA_INPUT + "=" + batch);
      Files.qsub(qsub, Array.toStr(Array.toStringArray(currentCommand), " "), 55000, 55,
                 c.getI(NUM_THREADS) * c.getI(NUM_THREADS_PIPELINE));
    }
  }



  public static void main(String[] args) {

    CLI c = new CLI();

    String sraDirDefault = "sra/";

    c.addArg(SRA_INPUT, "directory or filename with .sra files", sraDirDefault);

    String outDir = "out/";

    c.addArg(OUT_DIR, "the output directory for results", outDir);

    String sraRunTableDefault = "sraRuntable.txt";

    c.addArg(SRA_RUN_TABLE, "a sra run table providing sample information", sraRunTableDefault);

    int numThreads = 24;

    c.addArg(NUM_THREADS, "number of threads across samples", Integer.toString(numThreads));

    int numThreadsPipe = 1;
    c.addArg(NUM_THREADS_PIPELINE, "number of threads within samples",
             Integer.toString(numThreadsPipe));

    String refGenomeFasta = "hg19.canonical.fa";


    c.addArg(REFERENCE_GENOME, "appropriate reference genome file", refGenomeFasta);

    String captureBedFile = "VCRome_2_1_hg19_capture_targets.bed";

    c.addArg(CAPTURE_BED, "bed file of targeted capture", captureBedFile);

    String binBed = "targetsOfInterest.bed";

    c.addArg(BIN_BED, "bed file of targets of interests", binBed);

    String vcf = "vcf.vcf";
    c.addArg(VCF, "vcf file of variants", vcf);

    String computelLocation = null;
    c.addArg(COMPUTEL, "directory of computel", computelLocation);

    int batch = -1;
    c.addArg(NUM_BATCHES, "number of batches", Integer.toString(batch));



    c.addFlag(COMPILE, "Compile the genvisis portion of the pipeline");



    c.parseWithExit(SRAPipeline.class, args);


    if (c.has(COMPILE)) {
      compile(c.get(SRA_INPUT), c.get(SRA_RUN_TABLE), c.get(OUT_DIR), c.get(REFERENCE_GENOME),
              c.get(CAPTURE_BED), c.get(BIN_BED), c.get(VCF), c.getI(NUM_THREADS));
    } else {
      runAll(c.get(SRA_INPUT), c.get(SRA_RUN_TABLE), c.get(OUT_DIR), c.get(REFERENCE_GENOME),
             c.get(CAPTURE_BED), c.get(BIN_BED), c.get(VCF), c.get(COMPUTEL), c.getI(NUM_THREADS),
             c.getI(NUM_THREADS_PIPELINE), c.getI(NUM_BATCHES), c);
    }



  }

}
