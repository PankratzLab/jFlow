package org.pankratzlab.internal.utils.sra;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.seq.manage.BamImport;
import org.genvisis.cnv.seq.manage.BamSample.NORMALIZATON_METHOD;
import org.genvisis.seq.NGSSample;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.analysis.MitoSeqCN;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.telomere.Computel;
import org.genvisis.seq.telomere.TelSeq;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerHive;
import org.pankratzlab.common.ext;

/**
 * Going to be the pipeline of execution for a single input bam file
 */
public class Pipeline {

  private static final String MITO_DIR = "mtDNACN/";
  private static final String TELSEQ_DIR = "telseq/";
  private static final String COMPUTEL_DIR = "computel/";
  private static final String UNMAPPED_DIR = "unmapped/";

  private static final int[] TELOMERE_CAPTURE_BUFFER = new int[] {100, 0, 1000, 2000, 3000};

  private Pipeline() {

  }

  /**
   * abstract part of the pipeline
   */
  public abstract static class PipelinePart implements Callable<PipelinePart> {

    private List<String> input;
    private List<String> output;

    private PipelinePart() {}

    protected List<String> getInput() {
      return input;
    }

    protected void setInput(List<String> input) {
      this.input = input;
    }

    protected List<String> getOutput() {
      return output;
    }

    protected void setOutput(List<String> output) {
      this.output = output;
    }

  }

  /**
   * @author Kitty <br>
   *         The part of the pipeline that generates mtDNA CN estimates from a bam file
   */
  private static class MitoPipePart extends PipelinePart {

    private final String bamFile;
    private final String rootOutDir;
    private final String captureBed;
    private final String refGenome;
    private final NGSSample ngsSample;
    private final int numthreads;
    private final Logger log;

    private MitoPipePart(String bamFile, String rootOutDir, String captureBed, String refGenome,
                         NGSSample ngsSample, int numthreads, Logger log) {
      super();
      this.bamFile = bamFile;
      this.rootOutDir = rootOutDir;
      this.captureBed = captureBed;
      this.refGenome = refGenome;
      this.ngsSample = ngsSample;
      this.numthreads = numthreads;
      this.log = log;
    }

    @Override
    public PipelinePart call() throws Exception {
      String mitoDir = rootOutDir + MITO_DIR + ext.rootOf(bamFile) + "/";
      String bamList = mitoDir + "bam.list.txt";
      new File(mitoDir).mkdirs();
      Files.write(bamFile, bamList);
      String result = MitoSeqCN.run(bamList, mitoDir,
                                    ngsSample.getaType() == ASSAY_TYPE.WGS ? null : captureBed,
                                    refGenome, ngsSample.getaName(), ngsSample.getaType(),
                                    numthreads, log);

      ArrayList<String> input = new ArrayList<>();
      input.add(bamFile);
      setInput(input);
      ArrayList<String> output = new ArrayList<>();
      output.add(result);
      setOutput(output);
      return this;
    }

  }

  private static class TelSeqPart extends PipelinePart {

    private final String bamFile;
    private final String rootOutDir;
    private final String captureBed;
    private final NGSSample ngsSample;
    private final int numthreads;
    private final int captureBufferSize;
    private final Logger log;

    private TelSeqPart(String bam, String rootOutDir, String captureBed, NGSSample ngsSample,
                       int numthreads, int captureBufferSize, Logger log) {
      super();
      bamFile = bam;
      this.rootOutDir = rootOutDir;
      this.captureBed = captureBed;
      this.ngsSample = ngsSample;
      this.numthreads = numthreads;
      this.captureBufferSize = captureBufferSize;
      this.log = log;
    }

    @Override
    public PipelinePart call() throws Exception {
      String telSeqDir = rootOutDir + TELSEQ_DIR + ext.rootOf(bamFile) + "/";
      new File(telSeqDir).mkdir();
      String result = TelSeq.runTelSeq(new String[] {bamFile}, telSeqDir, captureBed, numthreads,
                                       ngsSample.getaType(), ngsSample.getaName(),
                                       captureBufferSize, log);
      ArrayList<String> input = new ArrayList<>();
      input.add(bamFile);
      setInput(input);
      ArrayList<String> output = new ArrayList<>();
      output.add(result);
      setOutput(output);

      return this;
    }

  }

  private static class GenvisisPart extends PipelinePart {

    private final String bamFile;
    private final String rootOutDir;
    private final String captureBed;
    private final NGSSample ngsSample;
    private final String binBed;
    private final int captureBufferSize;
    private final String vcf;
    private final String refGenome;

    private GenvisisPart(String bam, String rootOutDir, String captureBed, String binBed,
                         String vcf, NGSSample ngsSample, String refGenome, int captureBufferSize) {
      super();
      bamFile = bam;
      this.rootOutDir = rootOutDir;
      this.captureBed = captureBed;
      this.binBed = binBed;
      this.vcf = vcf;
      this.ngsSample = ngsSample;
      this.captureBufferSize = captureBufferSize;
      this.refGenome = refGenome;
    }

    /*
     * (non-Javadoc)
     * @see java.util.concurrent.Callable#call()
     */
    @Override
    public PipelinePart call() throws Exception {
      Project proj = getProjectFor(ngsSample.getaType(), rootOutDir);
      BamImport.importTheWholeBamProject(proj, binBed, captureBed, vcf, captureBufferSize, -1,
                                         false, ngsSample.getaType(), ngsSample.getaName(),
                                         NORMALIZATON_METHOD.GENOME, new String[] {bamFile},
                                         refGenome, false, true, 1);

      ArrayList<String> input = new ArrayList<>();
      input.add(bamFile);
      setInput(input);
      return this;
    }

  }

  private static class UnMappedPart extends PipelinePart {

    private final String inputBam;
    private final String rootOutDir;
    private final Logger log;

    public UnMappedPart(String inputBam, String rootOutDir, Logger log) {
      super();
      this.inputBam = inputBam;
      this.rootOutDir = rootOutDir;
      this.log = log;
    }

    /*
     * (non-Javadoc)
     * @see java.util.concurrent.Callable#call()
     */
    @Override
    public PipelinePart call() throws Exception {
      String unmappedDir = rootOutDir + UNMAPPED_DIR;
      new File(unmappedDir).mkdirs();
      String unMappedBam = unmappedDir + ext.rootOf(inputBam, true) + ".unmapped.bam";
      if (!Files.exists(unMappedBam)) {
        BamOps.dumpUnMappedReads(inputBam, unMappedBam, log);
      }
      return this;
    }

  }

  private static class ComputelPart extends PipelinePart {

    private final String bamFile;
    private final String rootOutDir;
    private final String computelLocation;
    private final Logger log;

    private ComputelPart(String bam, String rootOutDir, String computelLocation, String captureBed,
                         NGSSample ngsSample, int numthreads, int captureBufferSize, Logger log) {
      super();
      bamFile = bam;
      this.rootOutDir = rootOutDir;
      this.computelLocation = computelLocation;
      this.log = log;
    }

    @Override
    public PipelinePart call() throws Exception {
      String outputDir = rootOutDir + COMPUTEL_DIR + ext.rootOf(bamFile) + "/";
      new File(outputDir).mkdir();
      // TODO, remove capture bed prior to running
      Computel.runComputel(bamFile, outputDir, computelLocation, log);
      ArrayList<String> input = new ArrayList<>();
      input.add(bamFile);
      setInput(input);
      ArrayList<String> output = new ArrayList<>();
      output.add(null);
      setOutput(output);
      return this;
    }

  }

  /**
   * @param aType the {@link ASSAY_TYPE} to format the project for
   * @param rootOutDir where the results will be stored
   * @param genomeBuild the reference genome that will be set for the project
   * @return a project to use
   */
  public static Project getProjectFor(ASSAY_TYPE aType, String rootOutDir) {

    String projectName = aType.toString() + "_Genvisis_Project";
    String projectDir = rootOutDir + "genvisis/" + aType + "/";
    String projectFile = projectDir + projectName + ".properties";
    if (!Files.exists(projectFile)) {
      new File(projectDir).mkdirs();
      Files.writeArray(new String[] {"PROJECT_NAME=" + projectName,
                                     "PROJECT_DIRECTORY=" + projectDir},
                       projectFile);
    }
    Project proj = new Project(projectFile);
    proj.ARRAY_TYPE.setValue(ARRAY.NGS);
    proj.saveProperties();
    return proj;
  }

  /**
   * @author Kitty Used to specify which parts of the pipeline will be run
   */
  public enum PIPELINE_PARTS {
    /**
     * Create temporary source files for genvisis
     */
    GENVISIS,
    /**
     * Generate mtDNA CN estimates
     */
    MTDNACN,
    /**
     * Compute telomere length with TelSeq
     */
    TELSEQ,
    /**
     * Compute telomere length with Computel
     */
    COMPUTEL;
  }

  /**
   * @param inputBam The bam file to run through the pipeline
   * @param rootOutDir where output will be sent
   * @param genomeBuild
   * @param captureBed a bed file defining caputure regions
   * @param binBed a bed file defining input targets (intersection of this and capture bed will be
   *          imported)
   * @param vcf a vcf file defining variant sites
   * @param sample
   * @param parts the parts of the pipeline to run, see {@link PIPELINE_PARTS} for options
   * @param computelLocation the location of the computel directory, ideally freshly cloned
   * @param numThreads
   * @param log
   * @return
   */
  public static List<PipelinePart> pipeline(String inputBam, String rootOutDir, String captureBed,
                                            String binBed, String vcf, NGSSample sample,
                                            String refGenome, List<PIPELINE_PARTS> parts,
                                            String computelLocation, int numThreads, Logger log) {
    if (!Files.exists(inputBam)) {
      throw new IllegalArgumentException("Bam file " + inputBam + " must exist");
    }

    if (!Files.exists(refGenome)) {
      throw new IllegalArgumentException("Reference Genome for " + refGenome + " must exist");
    }
    if (sample.getaType() == ASSAY_TYPE.WXS && (!Files.exists(captureBed))) {
      throw new IllegalArgumentException(captureBed + " must exist");
    }

    WorkerHive<PipelinePart> hive = new WorkerHive<>(numThreads, 10, log);
    hive.addCallable(new UnMappedPart(inputBam, rootOutDir, log));// it's cheap to do
    for (PIPELINE_PARTS part : parts) {
      switch (part) {
        case COMPUTEL:
          hive.addCallable(new ComputelPart(inputBam, rootOutDir, computelLocation, captureBed,
                                            sample, 1, TELOMERE_CAPTURE_BUFFER[0], log));
          break;
        case GENVISIS:
          hive.addCallable(new GenvisisPart(inputBam, rootOutDir, captureBed, binBed, vcf, sample,
                                            refGenome, BamImport.CAPTURE_BUFFER));
          break;
        case MTDNACN:
          hive.addCallable(new MitoPipePart(inputBam, rootOutDir, captureBed, refGenome, sample, 1,
                                            log));
          break;
        case TELSEQ:
          for (int element : TELOMERE_CAPTURE_BUFFER) {
            hive.addCallable(new TelSeqPart(inputBam, rootOutDir, captureBed, sample, 1, element,
                                            log));
          }
          break;
        default:

          throw new IllegalArgumentException(part + " not implemented yet");
      }
    }
    hive.execute(true);
    return hive.getResults();
  }

}
