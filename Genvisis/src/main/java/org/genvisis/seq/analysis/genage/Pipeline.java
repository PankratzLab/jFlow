package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.NGSSample;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.analysis.MitoSeqCN;
import org.genvisis.seq.manage.BamImport;
import org.genvisis.seq.telomere.TelSeq;

/**
 * Going to be the pipeline of execution for a single input bam file
 *
 */
public class Pipeline {

  private static final String MITO_DIR = "mtDNACN/";
  private static final String TELSEQ_DIR = "telseq/";

  private Pipeline() {

  }

  private abstract static class PipelinePart implements Callable<PipelinePart> {
    private PipelinePart() {}

    protected void setOutput(List<String> output) {}

  }

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
      ArrayList<String> output = new ArrayList<String>();
      output.add(result);
      setOutput(output);
      return this;
    }

  }

  private static class TelSeqPart extends PipelinePart {

    private final String bam;
    private final String rootOutDir;
    private final String captureBed;
    private final NGSSample ngsSample;
    private final int numthreads;
    private final int captureBufferSize;
    private final Logger log;

    private TelSeqPart(String bam, String rootOutDir, String captureBed, NGSSample ngsSample,
                       int numthreads, int captureBufferSize, Logger log) {
      super();
      this.bam = bam;
      this.rootOutDir = rootOutDir;
      this.captureBed = captureBed;
      this.ngsSample = ngsSample;
      this.numthreads = numthreads;
      this.captureBufferSize = captureBufferSize;
      this.log = log;
    }

    @Override
    public PipelinePart call() throws Exception {
      String telSeqDir = rootOutDir + TELSEQ_DIR + ext.rootOf(bam) + "/";
      new File(telSeqDir).mkdir();
      String result = TelSeq.runTelSeq(new String[] {bam}, telSeqDir, captureBed, numthreads,
                                       ngsSample.getaType(), ngsSample.getaName(),
                                       captureBufferSize, log);
      ArrayList<String> output = new ArrayList<String>();
      output.add(result);
      setOutput(output);
      return this;
    }

  }

  private static class GenvisisPart extends PipelinePart {

    private final String bam;
    private final String rootOutDir;
    private final String captureBed;
    private final NGSSample ngsSample;
    private final String binBed;
    private final int numthreads;
    private final int captureBufferSize;
    private final String vcf;
    private final Logger log;



    public GenvisisPart(String bam, String rootOutDir, String captureBed, String binBed, String vcf,
                        NGSSample ngsSample, int numthreads, int captureBufferSize, Logger log) {
      super();
      this.bam = bam;
      this.rootOutDir = rootOutDir;
      this.captureBed = captureBed;
      this.binBed = binBed;
      this.vcf = vcf;
      this.ngsSample = ngsSample;
      this.numthreads = numthreads;
      this.captureBufferSize = captureBufferSize;
      this.log = log;
    }



    /*
     * (non-Javadoc)
     * 
     * @see java.util.concurrent.Callable#call()
     */
    @Override
    public PipelinePart call() throws Exception {

      BamImport.importTheWholeBamProject(null, binBed, captureBed, vcf, captureBufferSize, -1,
                                         false, ngsSample.getaType(), new String[] {bam},
                                         numthreads);

      return null;
    }

  }

  public static Project getProjectFor(ASSAY_TYPE aType, String rootOutDir) {


    return null;

  }


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

    WorkerHive<PipelinePart> hive = new WorkerHive<Pipeline.PipelinePart>(1, 10, log);
    // mtDNA CN
    hive.addCallable(new MitoPipePart(inputBam, rootOutDir, captureBed, referenceGenome, sample, 1,
                                      log));

    hive.addCallable(new TelSeqPart(inputBam, rootOutDir, captureBed, sample, 1, 100, log));

    hive.addCallable(new GenvisisPart(inputBam, rootOutDir, captureBed, binBed, vcf, sample, 1,
                                      BamImport.CAPTURE_BUFFER, log));

    hive.execute(true);

    return hive.getResults();
  }


}
