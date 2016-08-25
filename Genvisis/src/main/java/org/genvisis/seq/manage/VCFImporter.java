package org.genvisis.seq.manage;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.CentroidCompute.CentroidBuilder;
import org.genvisis.cnv.analysis.PennCNV;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFSamplePrep.PREPPED_SAMPLE_TYPE;
import org.genvisis.seq.manage.VCFSamplePrep.VCFSamplePrepWorker;
import org.genvisis.seq.manage.VCOps.LocusID;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.contamination.MAF;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Class to import a vcf into Genvisis format
 *
 */
public class VCFImporter {
  private final VCFFileReader vcfFileReader;
  private final Logger log;

  public VCFImporter(VCFFileReader vcfFileReader, Logger log) {
    super();
    this.vcfFileReader = vcfFileReader;
    this.log = log;
  }

  public ConversionResults importVCF(Project proj, Set<String> samples, final String vcfFile,
                                     boolean createMarkerPositions) {
    HashSet<String> verifiedSamples = verifySamples(proj, vcfFile, samples);
    ArrayList<LocusID> markers = new ArrayList<VCOps.LocusID>(1000000);
    log.reportTimeWarning("Only un-ambigous and biallelic variants will be exported...until we figure these out");
    FilterNGS.VariantContextFilter niceAllele =
                                              new FilterNGS.VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {},
                                                                                 new VARIANT_FILTER_BOOLEAN[] {VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER,
                                                                                                               VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER},
                                                                                 null, null, log);
    SampleNGS[] vcSamples = SampleNGS.getSamples(verifiedSamples);
    int numSkipped = 0;
    int total = 0;

    for (VariantContext vc : vcfFileReader) {
      total++;
      if (total % 100000 == 0) {
        log.reportTimeInfo("Processed " + total + " variants");
      }
      if (niceAllele.filter(vc).passed()) {
        markers.add(new LocusID(vc));
        for (SampleNGS vcSample : vcSamples) {
          HashSet<String> tmp = new HashSet<String>();
          tmp.add(vcSample.getSampleName());
          VariantContext vcSub = VCOps.getSubset(vc, tmp, VC_SUBSET_TYPE.SUBSET_STRICT);
          vcSample.addGeno(vc, vcSub.getGenotype(0), log);
        }
      } else {
        numSkipped++;
      }
    }
    if (createMarkerPositions) {
      generateMarkerPositions(proj, markers);
    }
    long fingerprint = proj.getMarkerSet().getFingerprint();
    Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
    for (SampleNGS vcSample : vcSamples) {
      vcSample.dump(proj, allOutliers, fingerprint, log);
    }
    log.reportTimeWarning(numSkipped + " non-biallelic or ambigous varaints of " + total
                          + " total variants were skipped");
    return new ConversionResults(allOutliers);
  }

  /**
   * Make sure the samples are present in the current vcf
   */
  private HashSet<String> verifySamples(Project proj, String vcf, Set<String> samples) {

    HashSet<String> samplesHave = new HashSet<String>();

    ArrayList<String> samplesMissing = new ArrayList<String>();
    String[] samplesVCF = VCFOps.getSamplesInFile(new VCFFileReader(new File(vcf), true));
    if (samples == null) {
      for (String element : samplesVCF) {
        samplesHave.add(element);
      }
    } else {
      for (String desiredSample : samples) {
        if (ext.indexOfStr(desiredSample, samplesVCF, true, true) >= 0) {
          samplesHave.add(desiredSample);
        } else {
          samplesMissing.add(desiredSample);
          proj.getLog()
              .reportTimeWarning("Could not find sample " + desiredSample + " in vcf, skipping");
        }
      }
      if (samplesMissing.size() > 0) {
        try {
          PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
                                                              + "missing_samples.txt"));
          for (String missingSample : samplesMissing) {
            writer.print(missingSample);
          }
          writer.close();
        } catch (Exception e) {
          log.reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue()
                          + "missing_samples.txt");
          log.reportException(e);
        }
      }
    }

    return samplesHave;
  }

  /**
   * @param proj
   * @param markers
   */
  private void generateMarkerPositions(Project proj, ArrayList<LocusID> markers) {
    try {
      PrintWriter writer =
                         new PrintWriter(new FileWriter(proj.MARKER_POSITION_FILENAME.getValue(true,
                                                                                               false)));
      writer.println("Marker\tChr\tPosition");
      for (LocusID id : markers) {
        writer.println(id.getId() + "\t" + id.getChr() + "\t" + id.getStart());
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + proj.MARKER_POSITION_FILENAME.getValue());
      log.reportException(e);
    }
    new File(ext.parseDirectoryOfFile(proj.MARKERSET_FILENAME.getValue())).mkdirs();
    Markers.orderMarkers(null, proj.MARKER_POSITION_FILENAME.getValue(),
                         proj.MARKERSET_FILENAME.getValue(), proj.getLog());
  }

  /**
   * stores the outliers
   *
   */
  private static class ConversionResults {
    private final Hashtable<String, Float> allOutliers;

    public ConversionResults(Hashtable<String, Float> allOutliers) {
      super();
      this.allOutliers = allOutliers;
    }

    public Hashtable<String, Float> getAllOutliers() {
      return allOutliers;
    }

  }

  private static class VCFImporterWorker implements Callable<ConversionResults> {
    private final Project proj;
    private final String vcf;
    private final Set<String> samples;
    private final boolean createMarkerPositions;

    public VCFImporterWorker(Project proj, String vcf, Set<String> samples,
                             boolean createMarkerPositions) {
      super();
      this.proj = proj;
      this.vcf = vcf;
      this.samples = samples;
      this.createMarkerPositions = createMarkerPositions;
    }

    @Override
    public ConversionResults call() throws Exception {
      proj.getLog().reportTimeInfo("Exporting " + samples.size() + " samples on this round using "
                                   + Thread.currentThread().getName());
      VCFImporter importer = new VCFImporter(new VCFFileReader(new File(vcf), true), proj.getLog());
      ConversionResults results = importer.importVCF(proj, samples, vcf, createMarkerPositions);

      return results;
    }
  }

  private static HashSet<String> getset(String[] set) {
    HashSet<String> hset = new HashSet<String>();
    for (String element : set) {
      hset.add(element);
    }
    return hset;
  }

  public static void test2(Project proj, String vcf, int numRounds, int numThreads,
                           int numDecompressThreads) {
    PREPPED_SAMPLE_TYPE[] types = PREPPED_SAMPLE_TYPE.values();
    for (PREPPED_SAMPLE_TYPE type : types) {
      // if(types[i] ==PREPPED_SAMPLE_TYPE.NORMALIZED)
      try {
        test2(proj, vcf, numRounds, numThreads, numDecompressThreads, type);
      } catch (Exception e) {
        System.exit(1);
      }

    }

  }

  public static void test2(Project proj, String vcf, int numRounds, int numThreads,
                           int numDecompressThreads, PREPPED_SAMPLE_TYPE type) {
    String newProjectDir = proj.PROJECT_DIRECTORY.getValue() + type + "/";
    new File(newProjectDir).mkdirs();
    String newProjectFile = ext.addToRoot(proj.getPropertyFilename(), type + "");
    Files.copyFile(proj.getPropertyFilename(), newProjectFile);
    Project projNorm = new Project(newProjectFile, false);

    // TODO uncomment
    // projNorm.setProperty(Project.PROJECT_NAME, projNorm.getProperty(Project.PROJECT_NAME) + type
    // + "");
    // projNorm.setProperty(Project.PROJECT_DIRECTORY, newProjectDir);
    // projNorm.saveProperties();
    // projNorm.getDir(Project.DATA_DIRECTORY, true, true);
    // GcModel gcModel =
    // GcAdjustor.GcModel.populateFromFile(proj.getFilename(Project.GC_MODEL_FILENAME, false, true),
    // false, proj.getLog());
    // Files.copyFile(proj.getFilename(Project.MARKER_POSITION_FILENAME),
    // projNorm.getFilename(Project.MARKER_POSITION_FILENAME));
    // // Files.copyFile(proj.getFilename(Project.SAMPLELIST_FILENAME),
    // projNorm.getFilename(Project.SAMPLELIST_FILENAME));
    // Markers.orderMarkers(null, projNorm.getFilename(Project.MARKER_POSITION_FILENAME),
    // projNorm.getFilename(Project.MARKERSET_FILENAME, true, true), projNorm.getLog());
    // //
    // VCFSamplePrepWorker vPrepWorker = new VCFSamplePrepWorker(proj,
    // projNorm.getDir(Project.SAMPLE_DIRECTORY, true, true), type, gcModel);
    // Hashtable<String, Float> allNewOutliers = new Hashtable<String, Float>();
    // WorkerTrain<Hashtable<String, Float>> train = new WorkerTrain<Hashtable<String,
    // Float>>(vPrepWorker, numThreads, 0, proj.getLog());
    //
    // int index = 0;
    // while (train.hasNext()) {
    // index++;
    // if (index % 10 == 0) {
    // proj.getLog().reportTimeInfo(index + " of " + proj.getSamples().length);
    // }
    // Hashtable<String, Float> outliers = train.next();
    // if (outliers.size() > 0) {
    // allNewOutliers.putAll(outliers);
    // }
    // }
    // train.shutdown();
    // if (allNewOutliers.size() > 0) {
    // Files.writeSerial(allNewOutliers, projNorm.getDir(Project.SAMPLE_DIRECTORY) +
    // "outliers.ser");
    // }
    // proj.setProperty(Project.SAMPLE_SUBSET_FILENAME, proj.getProjectDir() + "samplesToUse.txt");
    // String[] samplesToUse = Array.subArray(proj.getSamples(), proj.getSamplesToInclude(null));
    // Files.writeList(samplesToUse, projNorm.getFilename(Project.SAMPLE_SUBSET_FILENAME));
    // TransposeData.transposeData(projNorm, 2000000000, false);
    //
    // Builder builder = new Builder();
    // String centFile = projNorm.getProjectDir() + ext.rootOf(vcf) + type + ".cent";
    // CentroidCompute.computeAndDumpCentroids(projNorm, centFile, builder, numThreads,
    // numDecompressThreads);
    // projNorm.saveProperties();
    // projNorm.setProperty(Project.CUSTOM_CENTROIDS_FILENAME, centFile);
    // Centroids.recompute(projNorm, centFile);
    // projNorm.setProperty(Project.INTENSITY_PC_FILENAME, "PCs" + type + ".PCs.extrapolated.txt");
    // projNorm.setProperty(Project.INTENSITY_PC_NUM_COMPONENTS, 10 + "");
    // //MitoPipeline.catAndCaboodle(projNorm, numDecompressThreads, "0.95",
    // proj.getFilename(Project.TARGET_MARKERS_FILENAME), 10, "PCs" + type, false, true, 0.98,
    // proj.getFilename(Project.SAMPLE_SUBSET_FILENAME), null, null, true, true, false, true);
    // TODO uncomment
    TransposeData.transposeData(projNorm, 2000000000, false);
    // String[] markers =
    // HashVec.loadFileToStringArray(proj.getFilename(Project.MARKER_POSITION_FILENAME), true, new
    // int[] { 0 }, true);

    MAF.computeMAF(projNorm);

    projNorm.saveProperties();
    // Files.writeList(markers, projNorm.getProjectDir() + "shadowMarkers.txt");
    // PennCNVPrep.exportSpecialPennCNV(projNorm, "shadowPrep" + type + "/", 10,
    // projNorm.getProjectDir() + "shadowMarkers.txt", 6, 1, false, false, false, 2);
    // PennCNVPrep.exportSpecialPennCNV(projNorm, "shadowPrep" + type + "/", 10,
    // projNorm.getProjectDir() + "shadowMarkers.txt", 6, 1, false, true, false, 2);

    // for (int i = 0; i < 10; i++) {
    // double qual = (double) i / 10;
    // builder.gcThreshold(qual);
    // CentroidCompute.computeAndDumpCentroids(projNorm, projNorm.getProjectDir() + ext.rootOf(vcf)
    // + type + "_" + qual + ".cent", builder, numThreads, numDecompressThreads);
    // projNorm.setProperty(Project.CUSTOM_CENTROIDS_FILENAME, ext.rootOf(vcf) + type + "_" + qual +
    // ".cent");
    // projNorm.saveProperties();
    // }
  }

  public static void test(Project proj, String vcf, String gc5Base, PREPPED_SAMPLE_TYPE type,
                          int numRounds, int numThreads) {
    String[] samples = VCFOps.getSamplesInFile(new VCFFileReader(new File(vcf), true));

    ArrayList<String[]> sampleChunks = Array.splitUpArray(samples, numRounds, proj.getLog());
    WorkerHive<ConversionResults> hive =
                                       new WorkerHive<VCFImporter.ConversionResults>(numThreads, 10,
                                                                                     proj.getLog());
    proj.SAMPLE_DIRECTORY.getValue(true, false);
    proj.MARKER_DATA_DIRECTORY.getValue(true, false);
    proj.DATA_DIRECTORY.getValue(true, false);
    Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
    if (proj.getSamples() == null || proj.getSamples().length == 0) {
      for (int i = 0; i < sampleChunks.size(); i++) {
        hive.addCallable(new VCFImporterWorker(proj, vcf, getset(sampleChunks.get(i)), i == 0));
      }
      hive.execute(true);
      hive.setReportEvery(1);
      ArrayList<ConversionResults> results = hive.getResults();
      for (int i = 0; i < results.size(); i++) {
        if (results.get(i).getAllOutliers().size() > 0) {
          allOutliers.putAll(results.get(i).getAllOutliers());
        }
      }
    }

    if (proj.getSamples() != null && proj.getSamples().length != samples.length) {
      proj.getLog().reportTimeError("A different number of samples appear to be parsed in "
                                    + proj.SAMPLE_DIRECTORY.getValue());
      proj.getLog().reportTimeError("Found " + samples.length + " samples in the vcf and "
                                    + proj.getSamples().length + " samples were parsed");
      return;
    }

    else {
      boolean haveAll = true;
      for (String sample : samples) {
        if (ext.indexOfStr(sample, proj.getSamples()) < 0) {
          proj.getLog()
              .reportTimeError("Did not find vcf sample " + sample + "in the sample directory");
          haveAll = false;
        }
      }
      if (!haveAll) {
        return;
      } else {
        proj.getLog()
            .reportTimeInfo("Detected all " + samples.length + " samples have been parsed");
        try {
          proj.getLog().reportTimeInfo("Constructing outlierHash");
          allOutliers = proj.loadOutliersFromSamples();
        } catch (Exception e) {
          proj.getLog().reportException(e);
          e.printStackTrace();
        }
      }
    }

    if (allOutliers.size() > 0) {
      SerializedFiles.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue() + "outliers.ser");
    }
    processExt(proj, gc5Base);

    processCentroids(proj, vcf, numThreads);

    String newProjectDir = proj.PROJECT_DIRECTORY.getValue() + type + "/";
    new File(newProjectDir).mkdirs();

    GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(), false,
                                                          proj.getLog());

    String newProjectFile = ext.addToRoot(proj.getPropertyFilename(), type + "");
    Files.copyFile(proj.getPropertyFilename(), newProjectFile);
    Project projNorm = new Project(newProjectFile, false);
    projNorm.PROJECT_NAME.setValue(projNorm.PROJECT_NAME.getValue() + type);
    projNorm.PROJECT_DIRECTORY.setValue(newProjectDir + type);
    projNorm.saveProperties();
    projNorm.DATA_DIRECTORY.getValue(true, false);

    if (projNorm.getSamples() == null || projNorm.getSamples().length == 0) {
      VCFSamplePrepWorker vPrepWorker = new VCFSamplePrepWorker(proj,
                                                                projNorm.SAMPLE_DIRECTORY.getValue(true,
                                                                                                   false),
                                                                PREPPED_SAMPLE_TYPE.NORMALIZED_GC_CORRECTED,
                                                                gcModel);
      Hashtable<String, Float> allNewOutliers = new Hashtable<String, Float>();
      WorkerTrain<Hashtable<String, Float>> train =
                                                  new WorkerTrain<Hashtable<String, Float>>(vPrepWorker,
                                                                                            numThreads,
                                                                                            0,
                                                                                            proj.getLog());

      int index = 0;
      while (train.hasNext()) {
        index++;
        if (index % 10 == 0) {
          proj.getLog().reportTimeInfo(index + " of " + proj.getSamples().length);
        }
        Hashtable<String, Float> outliers = train.next();
        if (outliers.size() > 0) {
          allNewOutliers.putAll(outliers);
        }
      }
      train.shutdown();
      if (allNewOutliers.size() > 0) {
        SerializedFiles.writeSerial(allNewOutliers,
                                    projNorm.SAMPLE_DIRECTORY.getValue() + "outliers.ser");
      }
      projNorm.getSamples();
    } else {
      projNorm.getLog()
              .reportTimeInfo("Assuming that all samples have been gc-corrected, skipping");
    }
    Files.copyFile(proj.MARKER_POSITION_FILENAME.getValue(),
                   projNorm.MARKER_POSITION_FILENAME.getValue());

    processExt(projNorm, gc5Base);
    projNorm.TARGET_MARKERS_FILENAMES.setValue(new String[] {ext.rootOf(vcf) + ".targetMarkers"});
    Files.writeList(projNorm.getAutosomalMarkers(),
                    projNorm.TARGET_MARKERS_FILENAMES.getValue()[0]);
    processCentroids(projNorm, vcf, numThreads);
    projNorm.LRRSD_CUTOFF.setValue(2.2);
    projNorm.SAMPLE_CALLRATE_THRESHOLD.setValue(0.98);
    projNorm.getLog().reportTimeError("Remember that the LRR_SD cutoff is set to 2.2");
    projNorm.saveProperties();
    String pretendMedian = projNorm.PROJECT_DIRECTORY.getValue() + "pretendMedian.txt";

    Files.writeList(Array.subArray(projNorm.getAutosomalMarkers(), 0, 100), pretendMedian);
    String useFile = projNorm.PROJECT_DIRECTORY.getValue() + "VCF_SAMPLES_TO_USE.txt";
    Files.writeList(Array.subArray(projNorm.getSamples(), projNorm.getSamplesToInclude(null, true)),
                    useFile);
    projNorm.getSamplesToInclude(null);
    MitoPipeline.catAndCaboodle(projNorm, numThreads, pretendMedian, 100,
                                projNorm.PROJECT_DIRECTORY.getValue() + "VCF_PCS", true, true, 0.98,
                                useFile, null, null, null, true, true, true, false, true, false,
                                null, -1, -1, GENOME_BUILD.HG19, MitoPipeline.DEFAULT_PVAL_OPTS,
                                null, false, true);
    SampleQC sampleQC = SampleQC.loadSampleQC(projNorm);
    sampleQC.addQCsToSampleData(5, true);
    sampleQC.addPCsToSampleData(5, 10, true);
  }

  private static void processExt(Project proj, String gc5Base) {
    if (!Files.exists(proj.SAMPLE_DATA_FILENAME.getValue())) {
      SampleData.createMinimalSampleData(proj);
    }
    if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {
      Markers.orderMarkers(null, proj.MARKER_POSITION_FILENAME.getValue(),
                           proj.MARKERSET_FILENAME.getValue(), proj.getLog());
    }
    if (!Files.exists(proj.MARKERLOOKUP_FILENAME.getValue())) {
      TransposeData.transposeData(proj, 2000000000, false);
    }
    if (!Files.exists(proj.GC_MODEL_FILENAME.getValue())) {
      PennCNV.gcModel(proj, gc5Base, proj.GC_MODEL_FILENAME.getValue(), 100);
    } else {
      proj.getLog().reportTimeInfo("Skipping transpose step");
    }
  }

  private static void processCentroids(Project proj, String vcf, int numThreads) {
    double qual = VARIANT_FILTER_DOUBLE.GQ.getDFilter() / 100;
    String cent = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(vcf) + "qual"
                  + ext.formDeci(qual, 2) + ".cent";
    if (!Files.exists(cent)) {
      CentroidBuilder builder = new CentroidBuilder();
      builder.gcThreshold(qual);
      builder.samplesToUse(proj.getSamplesToInclude(null));
      CentroidCompute.computeAndDumpCentroids(proj, cent, builder, numThreads, 3);
      proj.CUSTOM_CENTROIDS_FILENAME.setValue(cent);
      Centroids.recompute(proj, cent);
      TransposeData.transposeData(proj, 2000000000, false);
      proj.saveProperties();
    } else {
      proj.getLog().reportFileExists(proj.CUSTOM_CENTROIDS_FILENAME.getValue());
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "C:/workspace/Genvisis/projects/Project_Tsai_21_25_26_spector.properties";
    String vcf =
               "D:/data/Project_Tsai_21_25_26_spector/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf.gz";
    String gc5Base = "N:/statgen/NCBI/gc5Base.txt";
    int numRounds = 4;
    int numThreads = 4;
    String usage = "\n" + "seq.manage.VCFImporter requires 0-1 arguments\n";
    usage += "   (1) full path to project file name (i.e. file=" + filename + " (default))\n" + "";
    usage += "   (2) full path to a vcf file (i.e. vcf=" + vcf + " (default))\n" + "";
    usage += "   (3) number of rounds to process the vcf file in  (i.e. numRounds=" + numRounds
             + " (default))\n" + "";
    usage += PSF.Ext.getNumThreadsCommand(4, numThreads);

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("vcf=")) {
        vcf = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("numRounds")) {
        numRounds = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = ext.parseIntArg(arg);
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
      Project proj = new Project(filename, false);
      test(proj, vcf, gc5Base, PREPPED_SAMPLE_TYPE.NORMALIZED_GC_CORRECTED, numRounds, 1);
      // test2(proj, vcf, numRounds, numThreads, numDecompressThreads);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}

// VCFImporterProducer producer = new VCFImporterProducer(vcfFile, verifiedSamples);
// Producer<VariantContext> producer2 = new Producer<VariantContext>() {
// private final VCFFileReader reader = new VCFFileReader(vcfFile, true);
// private final CloseableIterator<VariantContext> iter = reader.iterator();
// private final VCFHeader headers = new VCFHeader(reader.getFileHeader());
//
// @Override
// public boolean hasNext() {
// return this.iter.hasNext();
// }
//
// @Override
// public Callable<VariantContext> next() {
//
// final VCFHeader header = new VCFHeader(this.reader.getFileHeader());
// final VariantContext vc = new VariantContextBuilder(iter.next()).make();
// //fullyDecode(reader.getFileHeader(), false);
//
// return new Callable<VariantContext>() {
// // private final VariantContext vch = new VariantContextBuilder(vc).make().fullyDecode(header,
// true);
//
//
// @Override
// public VariantContext call() throws Exception {
// //vch.getGenotypes();
// vc.getGenotypes().get(0);
// //LazyGenotypesContext lc =(LazyGenotypesContext) gc;
// // //System.out.println(lc.getUnparsedGenotypeData());
// // lc.getSampleNames();
// // lc.decode();
// // //System.out.println(lc.getUnparsedGenotypeData());
// // builder.genotypes(lc);
// //// gc.isLazyWithData();
// // TODO Auto-generated method stub
// // VariantContextBuilder vcBuilder = new VariantContextBuilder(vch);
// // vcBuilder.genotypes(vch.getGenotypes());
// // vcBuilder.alleles( vch.getAlleles());
// // //vcBuilder.genotypes(lc);
// // VariantContext decode = vcBuilder.make();
//
// // vcBuilder.
// // VariantContext decodeGeno =builder.make();
// // lc =(LazyGenotypesContext) decodeGeno.getGenotypes();
// return vc;
//
// }
// };
// }
//
// @Override
// public void remove() {
// // TODO Auto-generated method stub
//
// }
//
// @Override
// public void shutdown() {
// reader.close();
// // TODO Auto-generated method stub
//
// }
// };
// WorkerTrain<VariantContext> train = new WorkerTrain<VariantContext>(producer2, 2, 100, log);
//
// int num =100;
// int index2=0;
// VariantContext[] vcs = new VariantContext[num];
// for (VariantContext vc : vcfFileReader) {
// vcs[index2] =vc;
// index2++;
// if(index2>= num){
// break;
// }
// }
// vcfFileReader.close();
// vcfFileReader =null;
// // for (int i = 0; i < vcs.length; i++) {
// // VariantContext vc = vcs[i];
//
// private static class VCFDecodeWorker implements Callable<VariantContext> {
// private final VariantContext vc;
// private final VCFHeader vcfHeader;
// private final HashSet<String> verifiedSamples;
//
// public VCFDecodeWorker(VariantContext vc, VCFHeader vcfHeader, HashSet<String> verifiedSamples) {
// super();
// this.vc = vc;
// this.vcfHeader = vcfHeader;
// this.verifiedSamples = verifiedSamples;
// }
//
// @Override
// public VariantContext call() throws Exception {
// VariantContext vcDecode = null;
// try {
// vcDecode = vc.fullyDecode(vcfHeader, false);
// System.out.println(vcDecode.isFullyDecoded());
// GenotypesContext gc = vc.getGenotypes();
// gc.get(0).getDP();
// vc.getAlleles();
// vcDecode = vc;
// // vcDecode.fullyDecode(vcfHeader, false);
// gc = vcDecode.getGenotypes();
// vcDecode.getGenotypes().get(0).getSampleName();

// System.out.println(vc.getGenotypes().isLazyWithData());
// System.out.println(vc.getGenotypes().get(0).getSampleName());
// System.out.println(vc.getGenotypes().get(0).getGQ());
// System.exit(1);

// vcSub = VCOps.getSubset(builder.make(), verifiedSamples, VC_SUBSET_TYPE.SUBSET_STRICT);
// } catch (Exception e) {
// e.printStackTrace();
// System.exit(1);
// }
// return vcDecode;
//
// }
//
// }
//
// private static class VCFImporterProducer implements Producer<VariantContext> {
// private final VCFFileReader vcfFileReader;
// private final CloseableIterator<VariantContext> vcIterator;
// private final VCFHeader header;
// private final HashSet<String> verifiedSamples;
//
// public VCFImporterProducer(String vcfFile, HashSet<String> verifiedSamples) {
// super();
// this.vcfFileReader = new VCFFileReader(vcfFile, true);
// this.header =new VCFHeader( vcfFileReader.getFileHeader());
// this.vcIterator = vcfFileReader.iterator();
// this.verifiedSamples = verifiedSamples;
//
// }
//
// @Override
// public boolean hasNext() {
// // TODO Auto-generated method stub
// return vcIterator.hasNext();
// }
//
// @Override
// public Callable<VariantContext> next() {
// final VariantContext vc = vcIterator.next();
// final VCFHeader headerCopy = new VCFHeader(header);
// return new VCFDecodeWorker(vc, headerCopy, verifiedSamples);
// }
//
// @Override
// public void remove() {
// // TODO Auto-generated method stub
//
// }
//
// @Override
// public void shutdown() {
// vcfFileReader.close();
// }
// }

//
// private static class VCFImporterWorker implements Callable<VariantContext> {
//
// private HashSet<String> verifiedSamples;
// private VariantContext vc;
//
//
// public VCFImporterWorker(HashSet<String> verifiedSamples, VariantContext vc) {
// super();
// this.verifiedSamples = verifiedSamples;
// this.vc = vc;
// }
//
//
// @Override
// public VariantContext call() throws Exception {
// return VCOps.getSubset(vc, verifiedSamples, VC_SUBSET_TYPE.SUBSET_STRICT);
// }
// }

// Files.copyFile(proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false),
// proj.getFilename(Project.DISPLAY_MARKERS_FILENAME, false, false));
// String[] samplesToUse = Array.subArray(proj.getSamples(), proj.getSamplesToInclude(null));
// proj.setProperty(Project.SAMPLE_SUBSET_FILENAME, proj.getProjectDir() + "samplesToUse.txt");
// Files.writeList(samplesToUse, proj.getFilename(Project.SAMPLE_SUBSET_FILENAME));
// // PCA.computePrincipalComponents(proj, true, 100, false, true, true, true, true, true,
// proj.getFilename(Project.SAMPLE_SUBSET_FILENAME), "VCF_PCS_WOExclude_Centered");
// // PCA.computePrincipalComponents(proj, true, 100, false, false, true, true, true, true,
// proj.getFilename(Project.SAMPLE_SUBSET_FILENAME), "VCF_PCS_WOExclude_NO_Centered");
// TODO uncomment
// for (int i = 0; i < 99; i++) {
