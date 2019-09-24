/**
 * 
 */
package org.genvisis.cnv.hmm;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.Resource;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.hmm.CNVCaller.AdjustmentQC;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.cnv.seq.manage.BamImport.NGS_MARKER_TYPE;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.WorkerHive;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * Class to perform hmm optimization from an initial .hmm file (e.g as provided by PennCNV). Methods
 * currently mimic --train from detect_cnv.pl
 */

public class HmmOptimizer implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  /**
   * lrrs
   */
  private double[] o1;
  /**
   * bafs
   */
  private double[] o2;

  private double[] pfb;

  private int[] snpdist;

  private boolean[] copyNumberOnlyDef;

  private AdjustmentQC adjustmentQC;

  private static final int NEW_CHR_DIST = 100_000_000;

  /**
   * @param o1
   * @param o2
   * @param pfb
   * @param snpdist
   * @param copyNumberOnlyDef
   */
  public HmmOptimizer(double[] o1, double[] o2, double[] pfb, int[] snpdist,
                      boolean[] copyNumberOnlyDef) {
    super();
    this.o1 = o1;
    this.o2 = o2;
    this.pfb = pfb;
    this.snpdist = snpdist;
    this.copyNumberOnlyDef = copyNumberOnlyDef;
  }

  private HmmOptimizer(CNVCaller cnvCaller, Logger log) {
    this.o1 = cnvCaller.analysisLrrs;
    this.o2 = cnvCaller.analysisBafs;
    this.pfb = cnvCaller.analysisPfbs;
    this.copyNumberOnlyDef = cnvCaller.copyNumberDef;
    this.snpdist = new int[cnvCaller.analysisPositions.length];
    this.adjustmentQC = cnvCaller.adjustmentQC;
    adjustmentQC.setTotalMarkers(o1.length);
    HashMap<String, ArrayList<Integer>> chrIndices = new HashMap<>();
    boolean[] finalAnalysisSet = ArrayUtils.booleanArray(cnvCaller.markerSet.getMarkerNames().length,
                                                         false);
    cnvCaller.populateIndices(finalAnalysisSet, chrIndices, cnvCaller.markerSet.getChrs());

    for (String chr : chrIndices.keySet()) {
      int[] indices = Ints.toArray(chrIndices.get(chr));
      for (int i = 0; i < indices.length - 1; i++) {
        snpdist[indices[i]] = cnvCaller.analysisPositions[indices[i + 1]]
                              - cnvCaller.analysisPositions[indices[i]];
      }
      snpdist[indices[indices.length - 1]] = NEW_CHR_DIST;
    }
  }

  /**
   * THIS
   * https://github.com/WGLab/PennCNV/blob/122691c8178b1ac803c60c802dacec83dc2593a4/detect_cnv.pl#L675-L701
   * 
   * @param pennHmm
   * @return
   */
  private static HmmOptimizer generateRareEvents(PennHmm pennHmm, Logger log) {
    List<Double> o1 = new ArrayList<>();
    /**
     * bafs
     */
    List<Double> o2 = new ArrayList<>();

    List<Double> pfbVals = new ArrayList<>();

    List<Integer> snpdist = new ArrayList<>();

    // state 1

    for (int i = 0; i < 100; i++) {
      o1.add(pennHmm.getB1().getB_mean()[0]);
      o2.add(0.5);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);
    for (int i = 0; i < 1000; i++) {
      o1.add(0.0);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);

    for (int i = 0; i < 250; i++) {
      o2.add(0.001);
      o2.add(0.5);
      o2.add(0.5);
      o2.add(0.999);
    }
    for (int i = 0; i < 1100; i++) {
      pfbVals.add(0.5);
    }

    // state 2
    for (int i = 0; i < 100; i++) {
      o1.add(pennHmm.getB1().getB_mean()[1]);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);
    for (int i = 0; i < 50; i++) {
      o2.add(0.001);
      o2.add(0.999);
    }

    for (int i = 0; i < 1000; i++) {
      o1.add(0.0);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);

    for (int i = 0; i < 250; i++) {
      o2.add(0.001);
      o2.add(0.5);
      o2.add(0.5);
      o2.add(0.999);
    }
    for (int i = 0; i < 1100; i++) {
      pfbVals.add(0.5);
    }

    // state 4
    for (int i = 0; i < 100; i++) {
      o1.add(0.0);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);
    for (int i = 0; i < 50; i++) {
      o2.add(0.0);
      o2.add(1.0);
    }

    for (int i = 0; i < 1000; i++) {
      o1.add(0.0);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);

    for (int i = 0; i < 250; i++) {
      o2.add(0.001);
      o2.add(0.5);
      o2.add(0.5);
      o2.add(0.999);
    }
    for (int i = 0; i < 1100; i++) {
      pfbVals.add(0.5);
    }

    // state 5
    for (int i = 0; i < 100; i++) {
      o1.add(pennHmm.getB1().getB_mean()[4]);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);
    for (int i = 0; i < 25; i++) {
      o2.add(0.001);
      o2.add(0.33);
      o2.add(0.66);
      o2.add(0.999);

    }

    for (int i = 0; i < 1000; i++) {
      o1.add(0.0);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);

    for (int i = 0; i < 250; i++) {
      o2.add(0.001);
      o2.add(0.5);
      o2.add(0.5);
      o2.add(0.999);
    }
    for (int i = 0; i < 1100; i++) {
      pfbVals.add(0.5);
    }

    // state 5
    for (int i = 0; i < 100; i++) {
      o1.add(pennHmm.getB1().getB_mean()[5]);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);
    for (int i = 0; i < 25; i++) {
      o2.add(0.001);
      o2.add(0.75);
      o2.add(0.25);
      o2.add(0.5);

    }

    for (int i = 0; i < 1000; i++) {
      o1.add(0.0);
      snpdist.add(5000);
    }
    snpdist.set(snpdist.size() - 1, 10_000_000);

    for (int i = 0; i < 250; i++) {
      o2.add(0.001);
      o2.add(0.5);
      o2.add(0.5);
      o2.add(0.999);
    }
    for (int i = 0; i < 1100; i++) {
      pfbVals.add(0.5);
    }
    if (o1.size() != o2.size() || o1.size() != pfbVals.size() || o1.size() != snpdist.size()) {
      throw new IllegalStateException("Internal error when creating mock data");
    }
    log.reportTimeInfo("Adding mock data with " + o1.size() + "markers");

    return new HmmOptimizer(Doubles.toArray(o1), Doubles.toArray(o2), Doubles.toArray(pfbVals),
                            Ints.toArray(snpdist), ArrayUtils.booleanArray(o1.size(), false));
  }
  //
  // ArrayList<DATA_ADJUSTMENTS> da = new ArrayList<>();
  // da.add(DATA_ADJUSTMENTS.HANDLE_NAN);
  // da.add(DATA_ADJUSTMENTS.GC_ADJUST);
  // da.add(DATA_ADJUSTMENTS.SUBSET_TO_ANALYSIS_MARKERS);
  // da.add(DATA_ADJUSTMENTS.MEDIAN_ADJUST);
  // da.add(DATA_ADJUSTMENTS.ADJUST_HMM_SD);

  private static HmmOptimizer preprocessSample(Project proj, PennHmm pennHmm, String sampleName,
                                               double[] sampLrrs, double[] sampBafs,
                                               GcModel gcModel, PFB pfb,
                                               PreparedMarkerSet markerSet, boolean[] markersToUse,
                                               int[] chrsToCall, boolean callReverse,
                                               int minNumMarkers, double minConf, int numThreads,
                                               PFB_MANAGEMENT_TYPE pManagementType,
                                               boolean debugMode) {

    CNVCaller cnvCaller = CNVCaller.getCNVCaller(proj, pennHmm, sampleName, sampLrrs, sampBafs,
                                                 gcModel, pfb, markerSet, markersToUse, chrsToCall,
                                                 callReverse, minNumMarkers, minConf, numThreads,
                                                 pManagementType, debugMode);

    proj.getLog().reportTimeInfo("BAF drift:" + cnvCaller.adjustmentQC.baf_drift);
    proj.getLog().reportTimeInfo("LRR SD:" + cnvCaller.adjustmentQC.lrrSD);
    proj.getLog().reportTimeInfo("BAF median:" + cnvCaller.adjustmentQC.baf_median);
    return new HmmOptimizer(cnvCaller, proj.getLog());
  }

  static String getSer(String outDir, String sample) {
    return outDir + sample + "hmm.prep.ser.gz";
  }

  static AdjustmentQC process(Project proj, PreparedMarkerSet markerSet, PennHmm pennHmmOriginal,
                              PFB pfb, GcModel gcModel, String sample, String ser) {
    if (!Files.exists(ser)) {
      proj.getLog().reportTimeInfo("Processing sample " + sample);
      Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
      int[] autosomalMarkers = proj.getAutosomalMarkerIndices();
      boolean[] markersToUse = ArrayUtils.booleanArray(markerSet.getMarkerNames().length, false);
      for (int i = 0; i < autosomalMarkers.length; i++) {
        markersToUse[autosomalMarkers[i]] = true;
        if (proj.ARRAY_TYPE.getValue() == ARRAY.NGS_WES) {
          String name = markerSet.getMarkerNames()[autosomalMarkers[i]];
          markersToUse[autosomalMarkers[i]] = NGS_MARKER_TYPE.getType(name) != NGS_MARKER_TYPE.VARIANT_SITE;
        }
      }

      HmmOptimizer processed = preprocessSample(proj, pennHmmOriginal, sample,
                                                ArrayUtils.toDoubleArray(samp.getLRRs()),
                                                ArrayUtils.toDoubleArray(samp.getBAFs()), gcModel,
                                                pfb, markerSet, markersToUse, null, false, -1, -1,
                                                1, PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, true);
      proj.getLog().reportTimeInfo("Finished Processing sample " + sample);

      SerializedFiles.writeSerial(processed, ser, true);
      proj.getLog().reportTimeInfo("Finished saving progress for sample " + sample);
      return processed.adjustmentQC;
    } else {
      proj.getLog().reportFileExists(ser);
      HmmOptimizer tmp = (HmmOptimizer) SerializedFiles.readSerial(ser, proj.getLog(), false, true);
      return tmp.adjustmentQC;
    }
  }

  private static void run(CLI c) {
    Project proj = new Project(c.get(CLI.ARG_PROJ));
    String outDir = proj.PROJECT_DIRECTORY.getValue() + "hmmOpt/";
    proj.getLog().reportTimeInfo("output will be sent to " + outDir);
    new File(outDir).mkdirs();
    PreparedMarkerSet markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
    PennHmm pennHmmOriginal = PennHmm.loadPennHmm(proj.HMM_FILENAME.getValue(), new Logger());

    if (!Files.exists(proj.CUSTOM_PFB_FILENAME.getValue())) {
      proj.getLog().reportTimeInfo("Did not find " + proj.CUSTOM_PFB_FILENAME.getValue()
                                   + ", attempting to generate it now");
      PFB.populationBAF(proj);
    }
    PFB pfb = PFB.loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
    if (!Files.exists(proj.GC_MODEL_FILENAME.getValue(false, false))) {
      Resource gmodelBase = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
                                     .getModelBase();
      if (!Files.exists(proj.GC_MODEL_FILENAME.getValue()) && gmodelBase.isAvailable()) {
        proj.getLog()
            .reportTimeWarning("Generating gcModel for " + proj.GENOME_BUILD_VERSION.getValue()
                               + " at " + proj.GC_MODEL_FILENAME.getValue() + " from "
                               + gmodelBase.get());
        proj.getLog().setLevel(3);
        GcModel.gcModel(proj, gmodelBase.get(), proj.GC_MODEL_FILENAME.getValue(), 100);
      }
    }
    GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false,
                                                                                          false),
                                                          false, proj.getLog());

    // String[] samples = new String[] {"10287118"};
    String[] samples = proj.getSamples();

    String adjustmentFile = outDir + "qc.adjusts.ser.gz";
    if (!Files.exists(adjustmentFile)) {
      generateAdjustments(c, proj, outDir, markerSet, pennHmmOriginal, pfb, gcModel, samples,
                          adjustmentFile);
    }

    proj.getLog().reportTimeInfo("Loading qc " + adjustmentFile);
    @SuppressWarnings("unchecked")
    List<AdjustmentQC> adjustmentQcs = (List<AdjustmentQC>) SerializedFiles.readSerial(adjustmentFile,
                                                                                       proj.getLog(),
                                                                                       false, true);

    proj.getLog()
        .reportTimeInfo("Loaded " + adjustmentQcs.size() + " qc  results from" + adjustmentFile);

    List<AdjustmentQC> adjustmentQcFiltered = new ArrayList<>();

    int totalNumMarkers = 0;
    for (AdjustmentQC aQc : adjustmentQcs) {
      // pennHmmOriginal.getB1().
      if (Math.abs(aQc.lrrSD - pennHmmOriginal.getB1().getB_sd()[2]) < 0.05) {
        if (Math.abs(aQc.baf_median - 0.5) < 0.05) {
          if (aQc.baf_drift < pennHmmOriginal.getB1().getB_uf() / 10) {
            adjustmentQcFiltered.add(aQc);
            totalNumMarkers += aQc.totalMarkers;
          }
        }
      }
    }

    proj.getLog().reportTimeInfo(adjustmentQcFiltered.size() + " samples retained post QC");
    if (!adjustmentQcFiltered.isEmpty()) {

      HmmOptimizer knownEvents = generateRareEvents(pennHmmOriginal, proj.getLog());
      proj.getLog().reportTimeInfo("Generated known events across " + knownEvents.o1.length
                                   + " mock markers");
      totalNumMarkers += knownEvents.o1.length;
      proj.getLog()
          .reportTimeInfo("Preparing hmm optimization across " + adjustmentQcFiltered.size()
                          + " samples and " + totalNumMarkers + " markers (count includes "
                          + knownEvents.o1.length + " mock markers)");

      double[] o1 = new double[totalNumMarkers];
      double[] o2 = new double[totalNumMarkers];
      double[] pfbVals = new double[totalNumMarkers];
      int[] snpdist = new int[totalNumMarkers];
      boolean[] copyNumberOnlyDef = new boolean[totalNumMarkers];

      int index = 0;
      for (AdjustmentQC aQc : adjustmentQcFiltered) {
        String ser = getSer(outDir, aQc.sample);
        HmmOptimizer tmp = (HmmOptimizer) SerializedFiles.readSerial(ser, proj.getLog(), false,
                                                                     true);
        copyValues(o1, o2, pfbVals, snpdist, copyNumberOnlyDef, index, tmp);

        index += tmp.o1.length;
      }
      copyValues(o1, o2, pfbVals, snpdist, copyNumberOnlyDef, index, knownEvents);

      proj.getLog()
          .reportTimeInfo("Finished preparing hmm optimization across "
                          + adjustmentQcFiltered.size() + " samples and " + totalNumMarkers
                          + " markers");
      proj.getLog().reportTimeInfo("Running BaumWelchNP");
      BaumWelchNP.BaumWelchNP_CHMM(pennHmmOriginal, o1, o2, pfbVals, snpdist, copyNumberOnlyDef,
                                   proj.getLog());
      proj.getLog().reportTimeInfo("Finished running BaumWelchNP");

    } else {
      proj.getLog().reportError("No samples remained following qc filtering");
    }
  }

  static void generateAdjustments(CLI c, Project proj, String outDir, PreparedMarkerSet markerSet,
                                  PennHmm pennHmmOriginal, PFB pfb, GcModel gcModel,
                                  String[] samples, String adjustmentFile) {
    WorkerHive<AdjustmentQC> preprocessHive = new WorkerHive<>(c.getI(CLI.ARG_THREADS), 10,
                                                               proj.getLog());
    for (String sample : samples) {
      String ser = getSer(outDir, sample);
      preprocessHive.addCallable(() -> process(proj, markerSet, pennHmmOriginal, pfb, gcModel,
                                               sample, ser));
    }
    preprocessHive.execute(true);
    List<AdjustmentQC> adjustmentQc = preprocessHive.getResults();
    SerializedFiles.writeSerial(adjustmentQc, adjustmentFile, true);
  }

  static void copyValues(double[] o1, double[] o2, double[] pfbVals, int[] snpdist,
                         boolean[] copyNumberOnlyDef, int index, HmmOptimizer tmp) {
    System.arraycopy(tmp.o1, 0, o1, index, tmp.o1.length);
    System.arraycopy(tmp.o2, 0, o2, index, tmp.o2.length);
    System.arraycopy(tmp.pfb, 0, pfbVals, index, tmp.pfb.length);
    System.arraycopy(tmp.snpdist, 0, snpdist, index, tmp.snpdist.length);
    System.arraycopy(tmp.copyNumberOnlyDef, 0, copyNumberOnlyDef, index,
                     tmp.copyNumberOnlyDef.length);
  }

  public static void main(String[] args) {
    CLI c = new CLI(HmmOptimizer.class);
    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, true);
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "2");
    c.parseWithExit(args);
    run(c);
  }

}
