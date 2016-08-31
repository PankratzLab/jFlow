package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.NGSSample;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.analysis.MitoSeqCN;
import org.genvisis.seq.manage.BamImport;
import org.genvisis.seq.telomere.Computel;
import org.genvisis.seq.telomere.TelSeq;

/**
 * Going to be the pipeline of execution for a single input bam file
 *
 */
public class Pipeline {

  private static final String MITO_DIR = "mtDNACN/";
  private static final String TELSEQ_DIR = "telseq/";
  private static final String COMPUTEL_DIR = "computel/";


  private Pipeline() {

  }

  /**
   * abstract part of the pipeline
   *
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
   *
   */
  private static class MitoPipePart extends PipelinePart {
    private final String bamFile;
    private final String rootOutDir;
    private final String captureBed;
    private final String referenceGenomeFasta;
    private final NGSSample ngsSample;
    private final int numthreads;
    private final Logger log;

    private MitoPipePart(String bamFile, String rootOutDir, String captureBed,
                         String referenceGenomeFasta, NGSSample ngsSample, int numthreads,
                         Logger log) {
      super();
      this.bamFile = bamFile;
      this.rootOutDir = rootOutDir;
      this.captureBed = captureBed;
      this.referenceGenomeFasta = referenceGenomeFasta;
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
                                    referenceGenomeFasta, ngsSample.getaName(),
                                    ngsSample.getaType(), numthreads, log);

      ArrayList<String> input = new ArrayList<String>();
      input.add(bamFile);
      setInput(input);
      ArrayList<String> output = new ArrayList<String>();
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
      this.bamFile = bam;
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
      ArrayList<String> input = new ArrayList<String>();
      input.add(bamFile);
      setInput(input);
      ArrayList<String> output = new ArrayList<String>();
      output.add(result);
      setOutput(output);
      return this;
    }

  }

  private static class GenvisisPart extends PipelinePart {

    private final String bamFile;
    private final String rootOutDir;
    private final String referenceGenomeFasta;
    private final String captureBed;
    private final NGSSample ngsSample;
    private final String binBed;
    private final int captureBufferSize;
    private final String vcf;



    private GenvisisPart(String bam, String rootOutDir, String referenceGenomeFasta,
                         String captureBed, String binBed, String vcf, NGSSample ngsSample,
                         int captureBufferSize) {
      super();
      this.bamFile = bam;
      this.rootOutDir = rootOutDir;
      this.referenceGenomeFasta = referenceGenomeFasta;
      this.captureBed = captureBed;
      this.binBed = binBed;
      this.vcf = vcf;
      this.ngsSample = ngsSample;
      this.captureBufferSize = captureBufferSize;
    }



    /*
     * (non-Javadoc)
     * 
     * @see java.util.concurrent.Callable#call()
     */
    @Override
    public PipelinePart call() throws Exception {
      Project proj = getProjectFor(ngsSample.getaType(), rootOutDir, referenceGenomeFasta);
      BamImport.importTheWholeBamProject(proj, binBed, captureBed, vcf, captureBufferSize, -1,
                                         false, ngsSample.getaType(), ngsSample.getaName(),
                                         new String[] {bamFile}, 1);

      ArrayList<String> input = new ArrayList<String>();
      input.add(bamFile);
      setInput(input);
      return this;
    }

  }


  private static class ComputelPart extends PipelinePart {

    private final String bamFile;
    private final String rootOutDir;
    private final String computelLocation;
    private final String captureBed;
    private final NGSSample ngsSample;
    private final int numthreads;
    private final int captureBufferSize;
    private final Logger log;

    private ComputelPart(String bam, String rootOutDir, String computelLocation, String captureBed,
                         NGSSample ngsSample, int numthreads, int captureBufferSize, Logger log) {
      super();
      this.bamFile = bam;
      this.rootOutDir = rootOutDir;
      this.computelLocation = computelLocation;
      this.captureBed = captureBed;
      this.ngsSample = ngsSample;
      this.numthreads = numthreads;
      this.captureBufferSize = captureBufferSize;
      this.log = log;
    }

    @Override
    public PipelinePart call() throws Exception {
      String computelDirectory = rootOutDir + COMPUTEL_DIR + ext.rootOf(bamFile) + "/";
      new File(computelDirectory).mkdir();
      Computel.test(bamFile, computelLocation, computelDirectory);

      String result = TelSeq.runTelSeq(new String[] {bamFile}, computelDirectory, captureBed,
                                       numthreads, ngsSample.getaType(), ngsSample.getaName(),
                                       captureBufferSize, log);
      ArrayList<String> input = new ArrayList<String>();
      input.add(bamFile);
      setInput(input);
      ArrayList<String> output = new ArrayList<String>();
      output.add(result);
      setOutput(output);
      return this;
    }

  }

  /**
   * @param aType the {@link ASSAY_TYPE} to format the project for
   * @param rootOutDir where the results will be stored
   * @param referenceGenomeFasta the reference genome that will be set for the project
   * @return a project to use
   */
  public static Project getProjectFor(ASSAY_TYPE aType, String rootOutDir,
                                      String referenceGenomeFasta) {
    String projectName = aType.toString() + "_Genvisis_Project";
    String projectDir = rootOutDir + "genvisis/" + aType + "/";
    String projectFile = projectDir + projectName + ".properties";
    if (!Files.exists(projectFile)) {
      new File(projectDir).mkdirs();
      Files.writeList(new String[] {"PROJECT_NAME=" + projectName,
                                    "PROJECT_DIRECTORY=" + projectDir},
                      projectFile);
    }
    Project proj = new Project(projectFile, false);
    proj.ARRAY_TYPE.setValue(ARRAY.NGS);
    proj.REFERENCE_GENOME_FASTA_FILENAME.setValue(referenceGenomeFasta);
    proj.saveProperties();
    return proj;
  }


  /**
   * @param inputBam The bam file to run through the pipeline
   * @param rootOutDir where output will be sent
   * @param referenceGenome
   * @param captureBed a bed file defining caputure regions
   * @param binBed a bed file defining input targets (intersection of this and capture bed will be
   *        imported)
   * @param vcf a vcf file defining variant sites
   * @param sample
   * @param numThreads
   * @param log
   * @return
   */
  public static List<PipelinePart> pipeline(String inputBam, String rootOutDir,
                                            String referenceGenome, String captureBed,
                                            String binBed, String vcf, NGSSample sample,
                                            int numThreads, Logger log) {
    if (!Files.exists(inputBam)) {
      throw new IllegalArgumentException("Bam file " + inputBam + " must exist");
    }

    if (!Files.exists(referenceGenome)) {
      throw new IllegalArgumentException("Reference Genome " + referenceGenome + " must exist");
    } else {
      log.reportTimeWarning("Assuming " + referenceGenome + " matches assembly type "
                            + sample.getaName());
    }
    if (sample.getaType() == ASSAY_TYPE.WXS && (!Files.exists(captureBed))) {
      throw new IllegalArgumentException(captureBed + " must exist");
    }

    WorkerHive<PipelinePart> hive = new WorkerHive<Pipeline.PipelinePart>(numThreads, 10, log);

    // mtDNA CN
    hive.addCallable(new MitoPipePart(inputBam, rootOutDir, captureBed, referenceGenome, sample, 1,
                                      log));

    hive.addCallable(new TelSeqPart(inputBam, rootOutDir, captureBed, sample, 1, 100, log));

    hive.addCallable(new GenvisisPart(inputBam, rootOutDir, referenceGenome, captureBed, binBed,
                                      vcf, sample, BamImport.CAPTURE_BUFFER));

    hive.execute(true);

    return hive.getResults();
  }


}
