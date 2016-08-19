package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.seq.analysis.PlinkSeq.ANALYSIS_TYPES;
import org.genvisis.seq.analysis.PlinkSeq.PlinkSeqProducer;
import org.genvisis.seq.analysis.PlinkSeq.PlinkSeqWorker;
import org.genvisis.seq.analysis.PlinkSeqUtils.PlinkSeqBurdenSummary;
import org.genvisis.seq.analysis.PlinkSeqUtils.PseqProject;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.ChrSplitResults;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.pathway.GenomeRegions;
import org.genvisis.seq.pathway.Pathways;

import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212
 *
 *         Does a big plink-seq association with no filtering besides population mac<br>
 *
 *
 */
public class PlinkSeqMegs {

  public static void runBig(String vcf, String vpopFile, String resourceDirectory,
                            String geneTrackFile, String keggPathwayFile, double maf,
                            int numthreads, Logger log) {
    VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.CASE_CONTROL, log);
    vpop.report();

    GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);

    geneTrack.setGeneSetFilename(geneTrackFile);
    Pathways pathways = Pathways.load(keggPathwayFile);
    GenomeRegions gRegions = new GenomeRegions(geneTrack, pathways, log);

    VCFOps.verifyIndex(vcf, log);
    String locFile = resourceDirectory + ext.rootOf(gRegions.getGeneTrack().getGeneSetFilename())
                     + "_Gen.reg";
    PlinkSeqUtils.generatePlinkSeqLoc(gRegions, locFile, log);

    PlinkSeq plinkSeq = new PlinkSeq(false, true, log);

    PseqProject pseqProject = PlinkSeq.initialize(plinkSeq, ext.rootOf(vpop.getFileName()),
                                                  ext.parseDirectoryOfFile(vpop.getFileName())
                                                                                            + VCFOps.getAppropriateRoot(vcf,
                                                                                                                        true)
                                                                                            + "/",
                                                  vcf, vpop, resourceDirectory, true, true, log);
    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    // int macFilter = (int) Math.round((float) VCFOps.getSamplesInFile(reader).length * maf);
    reader.close();
    // System.exit(1);

    PlinkSeqWorker[] complete = plinkSeq.fullGamutAssoc(pseqProject,
                                                        new String[] {ext.rootOf(locFile)}, null,
                                                        -1, maf + "",
                                                        ext.rootOf(vpop.getFileName()), numthreads);
    PlinkSeqBurdenSummary[] summaries = new PlinkSeqBurdenSummary[1];
    int index = 0;
    for (PlinkSeqWorker element : complete) {
      ANALYSIS_TYPES type = element.getType();
      switch (type) {
        case BURDEN:
          String analysis = ext.rootOf(element.getOutputFiles()[0]);
          analysis = analysis.replaceAll(".*" + ext.rootOf(locFile) + ".", "");
          PlinkSeqBurdenSummary plinkSeqBurdenSummary = new PlinkSeqBurdenSummary(analysis,
                                                                                  element.getOutputFiles()[0],
                                                                                  log);
          plinkSeqBurdenSummary.load();
          plinkSeqBurdenSummary.correctPvalues();
          summaries[index] = plinkSeqBurdenSummary;
          index++;
          break;
        case I_SUMMARY:
          break;
        case V_ASSOC:
          break;
        case V_SUMMARY:
          break;
        default:
          log.reportTimeError("INVALID analysis type " + type);
          break;
      }
    }
  }

  private static class ImportProducer extends AbstractProducer<PlinkSeqWorker[]> {
    private final String[] vcfs;
    private final String vpopFile;
    private final String resourceDirectory;
    private final String geneTrackFile;
    private final String keggPathwayFile;
    private final double maf;
    private boolean loadLoc;
    private final Logger log;
    private int index;

    public ImportProducer(String[] vcfs, String vpopFile, String resourceDirectory,
                          String geneTrackFile, String keggPathwayFile, double maf, boolean loadLoc,
                          Logger log) {
      super();
      this.vcfs = vcfs;
      this.vpopFile = vpopFile;
      this.resourceDirectory = resourceDirectory;
      this.geneTrackFile = geneTrackFile;
      this.keggPathwayFile = keggPathwayFile;
      this.maf = maf;
      this.log = log;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < vcfs.length;
    }

    @Override
    public Callable<PlinkSeqWorker[]> next() {
      final String tmpVcf = vcfs[index];
      System.out.println("DSFDSFDS");
      Callable<PlinkSeqWorker[]> callable = new Callable<PlinkSeq.PlinkSeqWorker[]>() {

        @Override
        public PlinkSeqWorker[] call() throws Exception {
          return prepToImport(tmpVcf, vpopFile, resourceDirectory, ext.parseDirectoryOfFile(tmpVcf),
                              geneTrackFile, keggPathwayFile, maf, loadLoc && index == 0, 1, log);
        }

      };
      index++;
      return callable;
    }
  }

  public static void runVcfs(String vcfFile, String vpopFile, String resourceDirectory,
                             String geneTrackFile, String keggPathwayFile, double maf,
                             int numthreads, boolean loadLoc, Logger log) {

    String[] vcfs = HashVec.loadFileToStringArray(vcfFile, false, new int[] {0}, true);
    System.out.println(Array.toStr(vcfs));

    ArrayList<PlinkSeqWorker> workers = new ArrayList<PlinkSeq.PlinkSeqWorker>();
    ImportProducer importer = new ImportProducer(vcfs, vpopFile, resourceDirectory, geneTrackFile,
                                                 keggPathwayFile, maf, loadLoc, log);
    WorkerTrain<PlinkSeqWorker[]> importTrain = new WorkerTrain<PlinkSeq.PlinkSeqWorker[]>(importer,
                                                                                           numthreads,
                                                                                           numthreads,
                                                                                           log);
    while (importTrain.hasNext()) {
      PlinkSeqWorker[] tmp = importTrain.next();
      for (PlinkSeqWorker element : tmp) {
        workers.add(element);
      }
    }
    PlinkSeqProducer producer =
                              new PlinkSeqProducer(workers.toArray(new PlinkSeqWorker[workers.size()]),
                                                   log);
    WorkerTrain<PlinkSeqWorker> train =
                                      new WorkerTrain<PlinkSeq.PlinkSeqWorker>(producer, numthreads,
                                                                               numthreads, log);
    while (train.hasNext()) {
      train.next();
    }
  }

  private static PlinkSeqWorker[] prepToImport(String vcf, String vpopFile,
                                               String resourceDirectory, String directory,
                                               String geneTrackFile, String keggPathwayFile,
                                               double maf, boolean loadLoc, int numthreads,
                                               Logger log) {
    VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.CASE_CONTROL, log);
    vpop.report();

    GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);

    geneTrack.setGeneSetFilename(geneTrackFile);
    Pathways pathways = Pathways.load(keggPathwayFile);
    GenomeRegions gRegions = new GenomeRegions(geneTrack, pathways, log);

    VCFOps.verifyIndex(vcf, log);
    String locFile = resourceDirectory + ext.rootOf(gRegions.getGeneTrack().getGeneSetFilename())
                     + "_Gen.reg";
    PlinkSeqUtils.generatePlinkSeqLoc(gRegions, locFile, log);

    PlinkSeq plinkSeq = new PlinkSeq(false, true, log);
    // ext.parseDirectoryOfFile(vpop.getFileName())
    PseqProject pseqProject =
                            PlinkSeq.initialize(plinkSeq, ext.rootOf(vpop.getFileName()), directory,
                                                vcf, vpop, resourceDirectory, true, loadLoc, log);
    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    // int macFilter = (int) Math.round((float) VCFOps.getSamplesInFile(reader).length * maf);
    reader.close();
    return plinkSeq.fullGamutAssoc(pseqProject, new String[] {ext.rootOf(locFile)}, null, -1,
                                   maf + "", ext.rootOf(vpop.getFileName()), false, numthreads);
  }

  public static void prepareBatches(String vcf, String vpopFile, String fullPathTojarFile,
                                    int totalMemoryRequestedInMb, int walltimeRequestedInHours,
                                    int numthreads, int numBatches, Logger log) {
    log.reportTimeInfo("Utilizing " + (numthreads) + " total threads per batch");
    ChrSplitResults[] cSplitResults = VCFOps.splitByChrs(vcf, numthreads, true, log);
    ArrayList<ChrSplitResults[]> cSplitResultsBatched = Array.splitUpArray(cSplitResults,
                                                                           numBatches, log);
    String[] baseCommand =
                         Array.concatAll(PSF.Java.buildJavaCP(fullPathTojarFile),
                                         new String[] {"seq.analysis.PlinkSeqMegs",
                                                       "vpop=" + vpopFile, PSF.Ext.NUM_THREADS_COMMAND
                                                                           + numthreads + ""});
    ArrayList<String> batches = new ArrayList<String>();
    ArrayList<String> masterVCFS = new ArrayList<String>();
    String rootOut = ext.parseDirectoryOfFile(vpopFile);
    for (int i = 0; i < cSplitResultsBatched.size(); i++) {
      String batch = rootOut + "megs" + i + ".pbs";
      String vcfFile = rootOut + "vcfs." + i + ".txt";
      ArrayList<String> vcfs = new ArrayList<String>();
      for (int j = 0; j < cSplitResultsBatched.get(i).length; j++) {
        String tmpVCF = cSplitResultsBatched.get(i)[j].getOutputVCF();
        vcfs.add(tmpVCF);
        masterVCFS.add(tmpVCF);
      }
      Files.writeList(vcfs.toArray(new String[vcfs.size()]), vcfFile);
      String[] batchCommand = new String[] {"vcfs=" + vcfFile, (i == 0 ? "-loadLoc" : "")};

      batches.add("qsub -q " + batch);
      Files.qsub(batch, Array.toStr(Array.concatAll(baseCommand, batchCommand), " "),
                 totalMemoryRequestedInMb, walltimeRequestedInHours, numthreads);
    }
    Files.writeList(batches.toArray(new String[batches.size()]),
                    rootOut + ext.rootOf(vpopFile) + "master.sh");
    Files.writeList(masterVCFS.toArray(new String[masterVCFS.size()]),
                    rootOut + ext.rootOf(vpopFile) + "master.vcfs.txt");

    // Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir, batchRoot, maxJobsPerBatch,
    // forceMaxJobs, queueName, memoryPerProcRequestedInMb, totalMemoryRequestedInMb,
    // walltimeRequestedInHours);

    //
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String vcf =
    // "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
    // String vcf =
    // "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqTallyTest/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.errorRegions2.vcf";
    // String vcf =
    // "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqTallyTest/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vPopCaseControl.vcf";
    // String vcf =
    // "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.CUSHINGS.vcf.gz";
    String vcf =
               "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/vcf/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
    String fileOfVcfs = null;
    boolean batch = false;
    String vpopFile =
                    "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqProj_tsai_spector_joint_AgilentCaptureRecal/vPopCaseControl.txt";
    // String resourceDirectory = "/home/tsaim/public/bin/pseqRef/hg19/";
    String resourceDirectory = "/home/spectorl/public/bin/pseqRef/hg19/";
    Logger log = new Logger(ext.rootOf(vcf, false) + "tally.log");
    String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
    String keggPathwayFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/kegg.ser";
    String logfile = null;
    int numthreads = 4;
    int numBatches = 4;
    int memoryInMb = 62000;
    int wallTimeInHours = 24;
    String fullPathTojarFile = "/panfs/roc/groups/14/tsaim/lane0212/parkPseq.jar";
    boolean loadLoc = false;

    double maf = 0.05;
    String usage = "\n" + "seq.analysis.VCFTallyPSeq requires 0-1 arguments\n";
    usage += "   (1) vcf file (i.e. file=" + vcf + " (default))\n" + "";
    usage += "   (2) vpop file,(i.e. vpop=" + vpopFile + " (default))\n" + "";
    usage += "   (3) batch (i.e -batch (not the default))\n" + "";
    usage += PSF.Ext.getNumThreadsCommand(4, numthreads);
    usage += "   (5) number of batches,(i.e. numBatches=" + numBatches + " (default))\n" + "";
    usage += "   (6) file of vcfs to run(i.e. vcfs= (no default))\n" + "";
    usage += PSF.Ext.getMemoryMbCommand(7, memoryInMb);
    usage += PSF.Ext.getWallTimeCommand(8, wallTimeInHours);

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("vcf=")) {
        vcf = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("vcfs=")) {
        fileOfVcfs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("vpop=")) {
        vpopFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-batch")) {
        batch = true;
        numArgs--;
      } else if (arg.startsWith("-loadLoc")) {
        loadLoc = true;
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numthreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.MEMORY_MB)) {
        memoryInMb = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.WALLTIME_HRS)) {
        wallTimeInHours = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("numBatches=")) {
        numBatches = ext.parseIntArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      if (batch) {
        prepareBatches(vcf, vpopFile, fullPathTojarFile, memoryInMb, wallTimeInHours, numthreads,
                       numBatches, log);
      } else if (fileOfVcfs != null) {
        runVcfs(fileOfVcfs, vpopFile, resourceDirectory, geneTrackFile, keggPathwayFile, maf,
                numthreads, loadLoc, log);
      } else {
        runBig(vcf, vpopFile, resourceDirectory, geneTrackFile, keggPathwayFile, maf, numthreads,
               log);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
