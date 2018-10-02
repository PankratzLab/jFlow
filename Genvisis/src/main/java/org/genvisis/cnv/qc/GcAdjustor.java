package org.genvisis.cnv.qc;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.stream.Stream;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.plots.ColorExt;
import org.genvisis.cnv.plots.ColorExt.ColorItem;
import org.genvisis.cnv.plots.ColorExt.ColorManager;
import org.genvisis.seq.manage.ReferenceGenome;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Positions;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;
import org.pankratzlab.core.CLI;
import org.pankratzlab.shared.filesys.Segment;
import org.pankratzlab.shared.stats.CrossValidation;
import org.pankratzlab.shared.stats.LeastSquares.LS_TYPE;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * Class for correcting intensity using gc content, see
 * http://nar.oxfordjournals.org/content/36/19/e126.full.pdf
 * <p>
 * We deviate slightly by determining median - binned gc content from a provided gc file, instead of
 * using the PennCNV reference below...however the refernence can be flagged
 * <p>
 * It is important to note that regression indices are determined by distance between markers, and
 * qc (WF,GCWF) are base-pair bins covering the genome
 * <p>
 * The current version takes ~ 0.5 seconds for a sample with ~600K markers. From PennCNV
 * <p>
 * #ref_median can be GC percentage in 135 1MB non-overlapping sliding windows in chr11 (2004 or
 * 2006 human genome assembly)
 * <p>
 * my @ref_median = qw/54.8207535282258 56.8381472081218 53.1218950320513 46.9484174679487
 * 39.9367227359694 38.3365384615385 41.9867788461538 40.4431401466837 44.5320512820513
 * 42.1979166666667 41.6984215561224 43.1598557692308 43.4388020833333 40.8104967948718
 * 39.8444475446429 41.5357572115385 38.7496995192308 45.0213249362245 42.3251201923077
 * 43.5287459935897 40.7440808354592 37.0492788461538 36.5006009615385 35.8518016581633
 * 35.2767427884615 35.1972155448718 36.5286192602041 39.4890825320513 36.5779246794872
 * 36.7275641025641 38.3256935586735 37.791266025641 41.1777844551282 41.950534119898
 * 42.3639823717949 41.9208733974359 41.2061543367347 35.4974959935897 35.2123397435897
 * 36.5101841517857 36.7135416666667 36.8268229166667 37.6945153061224 40.7453926282051
 * 47.7049278846154 47.3233173076923 44.7361288265306 46.6585536858974 39.1593549679487
 * 36.5684789540816 38.2718466806667 37.184425 37.184425 37.184425 37.184425 35.9227764423077
 * 41.1157852564103 41.6662348533163 39.7402844551282 40.0149238782051 46.6417211415816
 * 49.9136618589744 45.2016225961538 51.3019172512755 52.0818309294872 51.1320112179487
 * 49.9807185102302 49.9807185102302 49.5874473187766 50.547349024718 50.7186498397436
 * 45.6435347576531 46.3352363782051 42.4091546474359 46.6399274553571 43.7746394230769
 * 45.0160256410256 41.8526642628205 43.8899075255102 38.5112179487179 36.1038661858974
 * 36.1689851721939 39.8506610576923 37.0439703525641 36.8012595663265 40.2521033653846
 * 39.661858974359 37.5013769093564 35.5448717948718 36.9039979272959 35.2046274038462
 * 38.2195512820513 40.074537627551 40.7097355769231 40.5470753205128 38.4104380072343
 * 36.131109775641 35.3915264423077 34.9693080357143 36.2953725961538 37.9602363782051
 * 39.1942362882653 37.4464142628205 36.8879206730769 35.7242588141026 36.7556202168367
 * 37.0639022435897 40.6929086538462 38.385084502551 39.4121594551282 40.2410857371795
 * 42.0772879464286 43.2935697115385 43.2345753205128 40.9113919005102 44.9575320512821
 * 46.2513020833333 46.4753069196429 48.3886217948718 47.8520633012821 43.8001802884615
 * 39.808274872449 44.5042067307692 38.3835136217949 44.9097177933673 45.5366586538462
 * 41.7346754807692 39.2198461415816 41.9489182692308 44.3351362179487 42.7910754145408
 * 42.3190104166667 42.0425681089744 47.0514787946429 45.3482603740699/;
 */
public class GcAdjustor {

  public static final String[] GC_ADJUSTOR_TITLE = {"Adjust by GC content"};
  public static final String GC_BASE_FILE = "gc5base";

  /**
   * These defaults are also used by PennCNV - detect_cnv.pl, but might be worth adjusting sometime,
   * somewhere, in a galaxy called msi
   */
  public static final double DEFAULT_MINAUTOSOMALGC = 15d;
  public static final double DEFAULT_MAXAUTOSOMALGC = 80d;
  public static final double DEFAULT_MIN_DATA_VALUE = -1d;
  public static final double DEFAULT_MAX_DATA_VALUE = 1d;
  public static final int[] DEFAULT_REGRESSION_DISTANCE = {1000000};
  public static final int[] DEFAULT_SKIP_PER_CHR = {0};
  public static final int DEFAULT_NUM_SNP_MAD = 10;// not available as command line option currently
  public static final double[] DEFAULT_PENNCNV_CHR11_GC_BINS = {54.8207535282258, 56.8381472081218,
                                                                53.1218950320513, 46.9484174679487,
                                                                39.9367227359694, 38.3365384615385,
                                                                41.9867788461538, 40.4431401466837,
                                                                44.5320512820513, 42.1979166666667,
                                                                41.6984215561224, 43.1598557692308,
                                                                43.4388020833333, 40.8104967948718,
                                                                39.8444475446429, 41.5357572115385,
                                                                38.7496995192308, 45.0213249362245,
                                                                42.3251201923077, 43.5287459935897,
                                                                40.7440808354592, 37.0492788461538,
                                                                36.5006009615385, 35.8518016581633,
                                                                35.2767427884615, 35.1972155448718,
                                                                36.5286192602041, 39.4890825320513,
                                                                36.5779246794872, 36.7275641025641,
                                                                38.3256935586735, 37.791266025641,
                                                                41.1777844551282, 41.950534119898,
                                                                42.3639823717949, 41.9208733974359,
                                                                41.2061543367347, 35.4974959935897,
                                                                35.2123397435897, 36.5101841517857,
                                                                36.7135416666667, 36.8268229166667,
                                                                37.6945153061224, 40.7453926282051,
                                                                47.7049278846154, 47.3233173076923,
                                                                44.7361288265306, 46.6585536858974,
                                                                39.1593549679487, 36.5684789540816,
                                                                38.2718466806667, 37.184425,
                                                                37.184425, 37.184425, 37.184425,
                                                                35.9227764423077, 41.1157852564103,
                                                                41.6662348533163, 39.7402844551282,
                                                                40.0149238782051, 46.6417211415816,
                                                                49.9136618589744, 45.2016225961538,
                                                                51.3019172512755, 52.0818309294872,
                                                                51.1320112179487, 49.9807185102302,
                                                                49.9807185102302, 49.5874473187766,
                                                                50.547349024718, 50.7186498397436,
                                                                45.6435347576531, 46.3352363782051,
                                                                42.4091546474359, 46.6399274553571,
                                                                43.7746394230769, 45.0160256410256,
                                                                41.8526642628205, 43.8899075255102,
                                                                38.5112179487179, 36.1038661858974,
                                                                36.1689851721939, 39.8506610576923,
                                                                37.0439703525641, 36.8012595663265,
                                                                40.2521033653846, 39.661858974359,
                                                                37.5013769093564, 35.5448717948718,
                                                                36.9039979272959, 35.2046274038462,
                                                                38.2195512820513, 40.074537627551,
                                                                40.7097355769231, 40.5470753205128,
                                                                38.4104380072343, 36.131109775641,
                                                                35.3915264423077, 34.9693080357143,
                                                                36.2953725961538, 37.9602363782051,
                                                                39.1942362882653, 37.4464142628205,
                                                                36.8879206730769, 35.7242588141026,
                                                                36.7556202168367, 37.0639022435897,
                                                                40.6929086538462, 38.385084502551,
                                                                39.4121594551282, 40.2410857371795,
                                                                42.0772879464286, 43.2935697115385,
                                                                43.2345753205128, 40.9113919005102,
                                                                44.9575320512821, 46.2513020833333,
                                                                46.4753069196429, 48.3886217948718,
                                                                47.8520633012821, 43.8001802884615,
                                                                39.808274872449, 44.5042067307692,
                                                                38.3835136217949, 44.9097177933673,
                                                                45.5366586538462, 41.7346754807692,
                                                                39.2198461415816, 41.9489182692308,
                                                                44.3351362179487, 42.7910754145408,
                                                                42.3190104166667, 42.0425681089744,
                                                                47.0514787946429, 45.3482603740699};

  private final Project proj;
  private CrossValidation crossValidation;
  private List<Marker> correctedMarkers;// the indices we were able to correct, as ordered by the input
  // data
  private List<Double> lrrsToCorrect;
  private List<Double> gcsToCorrect;
  private List<? extends List<Marker>> qcMarkers;// the markers used to compute qc metrics, as ordered by the input-gc
  // matched arrays
  private List<? extends List<Marker>> chr11qcMarkers;// only used if defualy penncnv GCWF calculation is desired
  private final double minimumAutosomalGC, maximimumAutosomalGC, minIntensity, maxIntensity;

  private double wfPrior;

  private double wfPost;

  private double gcwfPrior;

  private double gcwfPost;
  private List<Double> regressionGcs;
  private List<Double> regressionIntensity;

  private final Map<Marker, Double> markerIntensities;
  private final Map<Marker, Double> markerGcs;

  private Map<Marker, Double> correctedIntensities;
  private final int regressionDistance, skipPerChr, numSnpMAD;
  private final Set<Marker> markerMask;
  private boolean fail;

  private final boolean verbose;
  private final GC_CORRECTION_METHOD correctionMethod;

  public enum GC_CORRECTION_METHOD {
    /**
     * Perform gc correction using the penncnv defaults (chr 11, etc)
     */
    PENNCNV_GC,
    /**
     * Perform gc correction using genvisis defaults (all available markers)
     */
    GENVISIS_GC;
  }

  // /**
  // * Constructor for PennCNV defaults, minus the chromosome 11 GC content bins for the GCWF
  // computation
  // */
  // public GcAdjustor(Project proj, PreparedMarkerSet preparedMarkerSet, GcModel gcModel, double[]
  // markerIntensities, boolean[] markerMask, GC_CORRECTION_METHOD correctionMethod, boolean
  // verbose) {
  // this(proj, preparedMarkerSet, gcModel, markerIntensities, DEFAULT_MINAUTOSOMALGC,
  // DEFAULT_MAXAUTOSOMALGC, DEFAULT_MIN_DATA_VALUE, DEFAULT_MAX_DATA_VALUE,
  // DEFAULT_REGRESSION_DISTANCE[0], DEFAULT_SKIP_PER_CHR[0], markerMask, correctionMethod,
  // verbose);
  // }
  //
  // /**
  // * @param proj
  // * current project
  // * @param gcModel
  // * valid {@link GcModel} object to use for correction
  // * @param markerIntensities
  // * array of intensities to correct
  // * @param minimumAutosomalGC
  // * gc contents less than this will not be included in the regression model
  // * @param maximimumAutosomalGC
  // * gc contents greater than this will not be included in the regression model
  // * @param minIntensity
  // * intensities less than this will not be included in the regression model
  // * @param maxIntensity
  // * intensities greater than this will not be included in the regression model
  // * @param regressionDistance
  // * the bin distance in base pairs, sliding window for regression, straight bins for qc
  // * @param skipPerChr
  // * skip this many markers at the beginning and end of each chromosome
  // * @param markerMask
  // * use at own peril, mainly a place holder currently
  // * @param PennCNVGCWF
  // * use the fixed PennCNV GC content array, based on chr 11
  // * @param verbose
  // * report
  // */
  // public GcAdjustor(Project proj, PreparedMarkerSet preparedMarkerSet, GcModel gcModel, double[]
  // markerIntensities, double minimumAutosomalGC, double maximimumAutosomalGC, double minIntensity,
  // double maxIntensity, int regressionDistance, int skipPerChr, boolean[] markerMask,
  // GC_CORRECTION_METHOD correctionMethod, boolean verbose) {
  // super();
  // this.proj = proj;
  // this.gcModel = gcModel;
  // this.markerIntensities = markerIntensities;
  // this.minimumAutosomalGC = minimumAutosomalGC;
  // this.maximimumAutosomalGC = maximimumAutosomalGC;
  // this.minIntensity = minIntensity;
  // this.maxIntensity = maxIntensity;
  // this.regressionDistance = regressionDistance;
  // this.skipPerChr = skipPerChr;
  // this.verbose = verbose;
  // this.fail = false;
  // this.wfPrior = Double.NaN;
  // this.wfPost = Double.NaN;
  // this.gcwfPrior = Double.NaN;
  // this.gcwfPost = Double.NaN;
  // this.correctionMethod = correctionMethod;
  // this.markerMask = markerMask;
  // this.preparedMarkerSet = preparedMarkerSet;
  // populateCurrentData();// Initialize everything needed
  // }

  public CrossValidation getCrossValidation() {
    return crossValidation;
  }

  public void correctIntensities() {
    correctIntensities(null, null);
  }

  /**
   * Builds a regression model from markers (usually ~3K) selected by
   * {@link GcAdjustor#populateCurrentData()}, the corrected intensities are the residuals after
   * applying the model to all intensities with a gc value
   */
  public void correctIntensities(String sample, GcAdjustorParameter gcParameters) {
    if (!fail) {
      double[] intensityArray = Doubles.toArray(lrrsToCorrect);
      double[][] gcArray = prepForRegression(gcsToCorrect.stream());
      if (gcParameters != null) {
        crossValidation = gcParameters.adjust(correctionMethod, sample, intensityArray, gcArray,
                                              verbose, proj.getLog());
      } else {
        crossValidation = new CrossValidation(Doubles.toArray(regressionIntensity),
                                              prepForRegression(regressionGcs.stream()),
                                              intensityArray, gcArray, true, LS_TYPE.REGULAR,
                                              proj.getLog());
        crossValidation.train();
        crossValidation.computePredictedValues();
        crossValidation.computeResiduals();
      }
      if (crossValidation.analysisFailed()) {
        if (verbose) {
          proj.getLog()
              .reportError("Error - the regression model has failed,  reverting to original intensity values");
        }
        assignOriginalIntensities();
      } else {
        if (verbose) {
          proj.getLog()
              .reportTimeInfo("regression model determined to be " + crossValidation.getBetas()[0]
                              + " + " + crossValidation.getBetas()[1] + " * GC content");
        }
        assignCorrectedIntensities();
      }
    } else {
      proj.getLog()
          .reportError("Error - cannot train gc model, reverting to original intensity values");
      assignOriginalIntensities();
    }
  }

  /**
   * @param computePrior compute wf and GCWF from original data
   * @param computePost compute wf and GCWF from corrected data, must call
   *          {@link GcAdjustor#correctIntensities()} first if this is flagged
   */

  public void computeQCMetrics(boolean computePrior, boolean computePost) {

    if (!fail) {
      if (qcMarkers == null || qcMarkers.isEmpty()) {
        proj.getLog()
            .reportError("Error - can not compute qc metrics, not enough qc markers were found");
        fail = true;
      }
      if (correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC
          && (chr11qcMarkers == null || chr11qcMarkers.isEmpty())) {
        proj.getLog()
            .reportError("Error - can not compute qc metrics using chromosome 11, not enough qc markers were found");
        fail = true;
      }
      if (!fail) {
        if (computePost && (crossValidation == null || crossValidation.analysisFailed())) {
          if (verbose) {
            proj.getLog()
                .report("Warning - intensity correction has not been performed or has failed,  corrected qc metrics will not be computed");
          }
        } else {
          if (computePost) {
            double[] wavesCorrected = getWave(correctedIntensities, markerGcs, qcMarkers,
                                              correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC ? chr11qcMarkers
                                                                                                  : null,
                                              verbose, proj.getLog());
            wfPost = wavesCorrected[0];
            gcwfPost = wavesCorrected[1];
          }
        }
        if (computePrior) {
          double[] wavesOriginal = getWave(markerIntensities, markerGcs, qcMarkers,
                                           correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC ? chr11qcMarkers
                                                                                               : null,
                                           verbose, proj.getLog());
          wfPrior = wavesOriginal[0];
          gcwfPrior = wavesOriginal[1];
        }
        if (verbose) {
          if (computePost || computePrior) {
            if (computePrior) {
              proj.getLog().reportTimeInfo("WF_PRIOR   -> " + wfPrior);
            }
            if (computePost) {
              proj.getLog().reportTimeInfo("WF_POST    -> " + wfPost);
            }
            if (computePrior) {
              proj.getLog().reportTimeInfo("GCWF_PRIOR -> " + gcwfPrior);
            }
            if (computePost) {
              proj.getLog().reportTimeInfo("GCWF_PRIOR -> " + gcwfPost);
            }
          }
        }
      } else {
        proj.getLog().reportError("Error - cannot compute qc metrics");
      }
    }
  }

  public boolean isFail() {
    return fail;
  }

  public String getQCString() {
    return wfPrior + "\t" + wfPost + "\t" + gcwfPrior + "\t" + gcwfPost;
  }

  public String getAnnotatedQCString() {
    return "WF_PRIOR: " + wfPrior + "\tWF_POST: " + wfPost + "\tGCWF_PRIOR" + gcwfPrior
           + "\tGCWF_POST" + gcwfPost;
  }

  public Map<Marker, Double> getCorrectedIntensities() {
    if (fail) {
      return markerIntensities;
    } else {
      return correctedIntensities;
    }
  }

  public double getWfPrior() {
    return wfPrior;
  }

  public double getWfPost() {
    return wfPost;
  }

  public double getGcwfPrior() {
    return gcwfPrior;
  }

  public double getGcwfPost() {
    return gcwfPost;
  }

  /**
   * We try to gather all we can in one pass, a long loop.... The loop iterates over markers by
   * chromosome, and gathers info for correcting and qc
   */
  private void populateCurrentData() {
    MarkerDetailSet markerSet = proj.getMarkerSet();
    if (markerIntensities.size() != markerSet.markersAsList().size()) {
      fail = true;
      if (verbose) {
        proj.getLog()
            .reportError("Error - the intensity data array must represent every marker in the project");
      }
    } else {
      List<Marker> markers = markerSet.markersAsList();

      ImmutableList.Builder<Double> regressGcsBuilder = ImmutableList.builderWithExpectedSize(3000);// usually around 3K, stores
      // info for regression model
      ImmutableList.Builder<Double> regressIntensityBuilder = ImmutableList.builderWithExpectedSize(3000);// usually around 3K,
      // stores info for
      // regression model
      ImmutableList.Builder<ImmutableList<Marker>> tmpQcMarkers = ImmutableList.builderWithExpectedSize(3000);// refers to the index matched
      // gc/LRR values
      ImmutableList.Builder<ImmutableList<Marker>> tmpchr11qcMarkers = ImmutableList.builderWithExpectedSize(3000);// refers to the index matched
      // gc/LRR values on chromosome
      // 11 only
      ImmutableList.Builder<Marker> correctedMarkersBuilder = ImmutableList.builderWithExpectedSize(markers.size());// refers to
      // all indices
      // we will
      // correct
      ImmutableList.Builder<Double> correctedLrrsBuilder = ImmutableList.builderWithExpectedSize(markers.size());
      ImmutableList.Builder<Double> correctedGcsBuilder = ImmutableList.builderWithExpectedSize(markers.size());

      int numPossibleBins = 0;
      for (Map.Entry<Byte, NavigableSet<Marker>> chrMarkers : markerSet.getChrMap().entrySet()) {
        int currentRegressDistance = 0;
        int currentBin = 0;
        List<Marker> tmpCurrentBin = new ArrayList<>(1000);
        List<Marker> tmpCurrentBinChr11 = new ArrayList<>(1000);// only used for
                                                                // chromosome 11...for
                                                                // default pennCNV
                                                                // behavior, and track
                                                                // even if not needed
        byte chr = chrMarkers.getKey();
        int markerCount = -1;
        for (Marker marker : chrMarkers.getValue()) {
          markerCount++;
          double gc = markerGcs.get(marker);// returns NaN if not present in
                                            // the gcModel file
          double intensity = markerIntensities.get(marker);
          if (!Double.isNaN(gc) && !Double.isNaN(intensity)) {// we take
                                                              // everything
                                                              // we can
                                                              // for
                                                              // correcting
                                                              // after
                                                              // the
                                                              // model
                                                              // is
                                                              // built
            correctedMarkersBuilder.add(marker);
            correctedLrrsBuilder.add(intensity);
            correctedGcsBuilder.add(gc);
          }
          if (chr > 0 && chr < 23) {// build from autosomal only
            if (markerCount >= skipPerChr
                && markerCount <= chrMarkers.getValue().size() - skipPerChr) {// allows to skip
              // first few and last
              // few markers per
              // chromosome
              if (markerMask == null || markerMask.contains(marker)) {// so we can skip poor
                                                                      // quality
                                                                      // markers/customize the
                                                                      // qc
                                                                      // note that this binning procedure is slightly different than the regression
                                                                      // indices
                int position = marker.getPosition();
                if (position >= regressionDistance * currentBin
                    && position < regressionDistance * (1 + currentBin)) {
                  if (!Double.isNaN(gc) && !Double.isNaN(intensity)) {
                    tmpCurrentBin.add(marker);
                    if (chr == 11) {
                      tmpCurrentBinChr11.add(marker);
                    }
                  }
                } else {// On to the next bin for qc, add indices if appropriate
                  numPossibleBins++;
                  currentBin++;
                  if (chr == 11) {
                    if (tmpCurrentBinChr11.size() > numSnpMAD) {
                      tmpchr11qcMarkers.add(ImmutableList.copyOf(tmpCurrentBinChr11));
                    } else {
                      tmpchr11qcMarkers.add(ImmutableList.of());
                    }
                  }
                  if (tmpCurrentBin.size() > numSnpMAD) {// skip small sized bins
                    tmpQcMarkers.add(ImmutableList.copyOf(tmpCurrentBin));
                  } else {
                    tmpQcMarkers.add(ImmutableList.of());
                  }
                  tmpCurrentBin = new ArrayList<>(1000);
                  tmpCurrentBinChr11 = new ArrayList<>(1000);
                  if (correctionMethod == GC_CORRECTION_METHOD.GENVISIS_GC) {// PennCNV skips this
                                                                             // marker, even though
                                                                             // it can be added to
                                                                             // a
                                                                             // bin, especially if
                                                                             // bins are small
                    if (regressionDistance * (currentBin + 1) <= position) {// this
                                                                            // is
                                                                            // in
                                                                            // order
                                                                            // to
                                                                            // handle
                                                                            // smaller
                                                                            // regression
                                                                            // distances
                      while (regressionDistance * (currentBin + 1) <= position) {
                        numPossibleBins++;
                        currentBin++;
                      }
                    }
                    if (position >= regressionDistance * currentBin
                        && position < regressionDistance * (1 + currentBin)) {
                      if (!Double.isNaN(gc) && !Double.isNaN(intensity)) {
                        tmpCurrentBin.add(marker);
                        if (chr == 11) {
                          tmpCurrentBinChr11.add(marker);
                        }
                      }
                    } else {
                      System.out.println(position + "\t" + (regressionDistance * currentBin) + "\t"
                                         + regressionDistance * (1 + currentBin));
                      proj.getLog().reportError("Should not happen");
                      System.exit(1);
                    }
                  }
                }
                if (!Double.isNaN(gc) && !Double.isNaN(intensity)) {// for
                                                                    // populating
                                                                    // regression
                                                                    // model
                  if (position - currentRegressDistance > regressionDistance) {

                    if (useMarker(intensity, gc, minimumAutosomalGC, maximimumAutosomalGC,
                                  minIntensity, maxIntensity)) {

                      // The gc and intensity will be used for the correction
                      regressGcsBuilder.add(gc);
                      regressIntensityBuilder.add(intensity);

                    }
                    currentRegressDistance = position;
                  }
                }
              }
            }
          }
        }
        if (tmpCurrentBin.size() > 0) {
          numPossibleBins++;

        }
        if (tmpCurrentBin.size() > numSnpMAD) {// add any leftovers
          tmpQcMarkers.add(ImmutableList.copyOf(tmpCurrentBin));
          if (chr == 11) {
            tmpchr11qcMarkers.add(ImmutableList.copyOf(tmpCurrentBinChr11));
          }
        } else {
          if (chrMarkers.getValue().size() > 0 && chr == 11) {
            tmpchr11qcMarkers.add(ImmutableList.of());
          }
          tmpQcMarkers.add(ImmutableList.of());
        }
        tmpCurrentBin = new ArrayList<>(1000);// reset the bins
        tmpCurrentBinChr11 = new ArrayList<>(1000);// reset the bins
      }
      regressionGcs = regressGcsBuilder.build();
      if (regressionGcs.size() == 0) {
        fail = true;
        regressionGcs = null;
        proj.getLog()
            .reportError("Error - could not find any autosomal markers to train the regression model");
        proj.getLog().reportError("Regression distance =" + regressionDistance);
        proj.getLog().reportError("Method =" + correctionMethod);

      } else {
        regressionIntensity = regressIntensityBuilder.build();
        correctedMarkers = correctedMarkersBuilder.build();
        lrrsToCorrect = correctedLrrsBuilder.build();
        gcsToCorrect = correctedGcsBuilder.build();
        if (verbose) {
          proj.getLog().report("Info - using " + regressionIntensity.size() + " of "
                               + markers.size() + " markers for regression model");
          proj.getLog().report("Info - " + correctedMarkers.size() + " of " + markers.size()
                               + " markers had a valid gc and a valid LRR for correction ");
        }
        if (regressionIntensity.size() == 0) {
          proj.getLog().reportError("Error - did not find enough valid markers for regression");
          fail = true;
        }

        qcMarkers = tmpQcMarkers.build();
        chr11qcMarkers = tmpchr11qcMarkers.build();
        if (verbose) {
          int totalMarkers = qcMarkers.stream().mapToInt(List::size).sum();
          proj.getLog()
              .report("Info - Correction method=" + correctionMethod + ", detected "
                      + qcMarkers.size() + " " + ext.prettyUpDistance(regressionDistance, 2)
                      + " sliding windows with >" + numSnpMAD + " markers out of a possible "
                      + numPossibleBins + " windows (" + totalMarkers
                      + " markers total) to compute WF"
                      + (correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC ? "" : " and GCWF"));
        }
      }
    }
  }

  /**
   * @param intensities gc matched intensities
   * @param gcs intensity matched gcs
   * @param wfBins bins for WF , and GCWF if pennCNVGCBins is not supplied should reflect all
   *          autosomal marker bins with valid gc/intensity
   * @param pennCNVGCBins (optional, can be null) if provided, this should reflect chromosome 11
   *          bins
   * @param verbose report some stuff
   * @param log
   * @return
   */
  private static double[] getWave(Map<Marker, Double> intensities, Map<Marker, Double> gcs,
                                  List<? extends List<Marker>> wfBins,
                                  List<? extends List<Marker>> pennCNVGCBins, boolean verbose,
                                  Logger log) {
    double[] waves = new double[2];// Organizes as WF, GCWF
    Arrays.fill(waves, Double.NaN);
    ArrayList<Double> medianIntensity = new ArrayList<>();
    ArrayList<Double> medianGc = new ArrayList<>();
    for (List<Marker> wFbin : wfBins) {
      if (!wFbin.isEmpty()) {
        medianIntensity.add(ArrayUtils.median(wFbin.stream().map(intensities::get), wFbin.size()));
        medianGc.add(ArrayUtils.median(wFbin.stream().map(gcs::get), wFbin.size()));
      }
    }

    double wf = ArrayUtils.mad(medianIntensity);
    if (pennCNVGCBins != null) {// Used for PennCNV bins, if not supplied we use what we found above
      if (pennCNVGCBins.size() != DEFAULT_PENNCNV_CHR11_GC_BINS.length) {
        log.reportError("Error - default PennCNV GC bins and current data do not match up, computing using full autosomal bins instead");
        log.reportError("Error - should have " + DEFAULT_PENNCNV_CHR11_GC_BINS.length
                        + " bins, but found " + pennCNVGCBins.size() + " bins instead");

      } else {
        medianIntensity = new ArrayList<>();
        medianGc = new ArrayList<>();
        for (int i = 0; i < pennCNVGCBins.size(); i++) {
          List<Marker> pennCNVGCBin = pennCNVGCBins.get(i);
          if (!pennCNVGCBin.isEmpty()) {
            medianGc.add(DEFAULT_PENNCNV_CHR11_GC_BINS[i]);
            medianIntensity.add(ArrayUtils.median(pennCNVGCBin.stream().map(intensities::get),
                                                  pennCNVGCBin.size()));
          }
        }
      }
      if (verbose) {
        log.report("Info - computing GCWF using " + medianGc.size() + " elements");
      }
    }
    double cc = org.pankratzlab.shared.stats.Correlation.Pearson(Doubles.toArray(medianIntensity),
                                                       Doubles.toArray(medianGc))[0];
    waves[0] = cc > 0 ? -1 * wf : wf;
    waves[1] = waves[0] * Math.abs(cc);
    return waves;
  }

  /**
   * When correction has failed, set to original
   */
  private void assignOriginalIntensities() {
    correctedIntensities = markerIntensities;
  }

  /**
   * Assign indices we have corrected to their new values
   */
  private void assignCorrectedIntensities() {
    Map<Marker, Double> correctedIntensitiesBuilder = Maps.newHashMap(markerIntensities);
    for (int i = 0; i < correctedMarkers.size(); i++) {
      Marker marker = correctedMarkers.get(i);
      correctedIntensitiesBuilder.put(marker, crossValidation.getResiduals()[i]);
    }
    correctedIntensities = Collections.unmodifiableMap(correctedIntensitiesBuilder);
  }

  private static double[][] prepForRegression(Stream<Double> gcs) {
    return gcs.map(g -> new double[] {g}).toArray(double[][]::new);
  }

  private static boolean useMarker(double intensity, double gc, double minimumAutosomalGC,
                                   double maximimumAutosomalGC, double minIntensity,
                                   double maxIntensity) {
    if (Double.isNaN(intensity)) {
      return false;
    }
    if (Double.isNaN(gc)) {
      return false;
    }
    if (gc <= minimumAutosomalGC || gc >= maximimumAutosomalGC) {
      return false;
    }
    if (intensity <= minIntensity || intensity >= maxIntensity) {
      return false;
    }
    return true;
  }

  /**
   * Put this in a separate class so it can be passed between samples
   */
  public static class GcModel implements Serializable {

    private static final long serialVersionUID = 1L;
    public static final String[] GC_HEADER = {"Name", "Chr", "Position", "GC"};
    public static final int DEFAULT_GC_MODEL_BIN_FASTA = 5120 * 100; // default bin generation for
                                                                    // reference genome-based gc
                                                                    // model
    private final String[] markers;
    private final byte[] chrs;
    private final int[] positions;
    private final double[] gcs;
    private Project mappedProj = null;
    private Map<Marker, Double> markerGcMap = null;
    private Map<String, Integer> index = new Hashtable<>();
    private ColorManager<String> colorManager;

    // private Logger log;

    public GcModel(GcModel gcmodel) {
      markers = gcmodel.markers;
      chrs = gcmodel.chrs;
      positions = gcmodel.positions;
      gcs = gcmodel.gcs;
      index = gcmodel.index;
      // this.log = gcmodel.log;
    }

    public GcModel(String[] markers, byte[] chrs, int[] positions, double[] gcs,
                   Map<String, Integer> index, Logger log) {
      super();
      this.markers = markers;
      this.chrs = chrs;
      this.positions = positions;
      this.gcs = gcs;
      this.index = index;
      // this.log = log;
    }

    private void developColorManager(int numBins, boolean redevelop) {
      if (redevelop || colorManager == null) {

        NormalDistribution nd = new NormalDistribution(ArrayUtils.mean(gcs, true),
                                                       ArrayUtils.stdev(getGcs(), true));
        Color[] colors = ColorExt.generatRGBScale(numBins); // bin gc to 100 bins
        Hashtable<String, String> lookup = new Hashtable<>();// items associated with
                                                             // category
                                                             // (marker->PoorQualityCategory)
        Hashtable<String, ColorItem<String>> manager = new Hashtable<>();
        for (int i = 0; i < gcs.length; i++) {
          int gcColorIndex = (int) Math.round(nd.cumulativeProbability(gcs[i]) * numBins - 1);
          gcColorIndex = Math.max(0, gcColorIndex);
          gcColorIndex = Math.min(numBins, gcColorIndex);
          lookup.put(markers[i], gcColorIndex + "");
          manager.put(gcColorIndex + "", new ColorItem<>(gcColorIndex + "", colors[gcColorIndex]));
        }
        colorManager = new ColorManager<String>(lookup, manager) {

          /**
           *
           */
          private static final long serialVersionUID = 1L;

        };
      }
    }

    public ColorManager<String> getColorManager() {
      if (colorManager == null) {
        developColorManager(100, false);
      }
      return colorManager;
    }

    public double[] getGcs() {
      return gcs;
    }

    public byte[] getChrs() {
      return chrs;
    }

    public int[] getPositions() {
      return positions;
    }

    // public Logger getLog() {
    // return log;
    // }

    public double[] getGcsFor(String[] markers) {
      double[] gcs = new double[markers.length];
      for (int i = 0; i < gcs.length; i++) {
        gcs[i] = getGcFor(markers[i]);
      }
      return gcs;
    }

    public double getGcFor(String marker) {
      if (hasMarker(marker)) {
        return gcs[index.get(marker)];
      } else {
        return Double.NaN;
      }
    }

    public double[] getGCsFor(String[] markers) {
      double[] gcs = new double[markers.length];
      for (int i = 0; i < gcs.length; i++) {
        gcs[i] = getGcFor(markers[i]);
      }
      return gcs;
    }

    public synchronized Map<Marker, Double> getMarkerGcMap(Project proj) {
      maybeSetProject(proj);
      if (markerGcMap == null) {
        markerGcMap = mappedProj.getMarkerSet().markersAsList().stream()
                                .collect(ImmutableMap.toImmutableMap(m -> m,
                                                                     m -> getGcFor(m.getName())));
      }
      return markerGcMap;
    }

    private void maybeSetProject(Project proj) {
      if (mappedProj == null) {
        mappedProj = proj;
      } else if (!mappedProj.equals(proj)) {
        throw new IllegalStateException("Same GC Model cannot be used with multiple projects");
      }
    }

    public boolean hasMarker(String marker) {
      return index.containsKey(marker);
    }

    public String[] getMarkers() {
      return markers;
    }

    public void Serialize(String fullPathToGcSer) {
      SerializedFiles.writeSerial(this, fullPathToGcSer);
    }

    public ColorManager<String> getGCColorManager() {

      return null;
    }

    /**
     * Generate the GCModel file needed by PennCNV software
     * (http://www.openbioinformatics.org/penncnv/).
     *
     * @param proj The project you are going to run PennCNV on.
     * @param inputGcBaseFullPath The user-supplied genome builds. Positions within each chromosome
     *          must be sorted by increasing order. For example,
     *          http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gc5Base.txt.gz
     * @param outputGcModelFullPath The name of the GCModel file.
     * @param numwindow For each SNP, GC Content is calculated for the range of numwindow*5120
     *          before its location through numwindow*5120 after. To use the default setting of 100,
     *          please enter 0. In order to be able to cross-reference with the same feature in
     *          PennCNV, this code intends to base on cal_gc_snp.pl in PennCNV package. But since
     *          the difference between Java and Perl, significant structural changes have been made.
     *          A sample gcBase file looks like below (There is no header line): 585 chr1 0 5120
     *          chr1.0 5 1024 0 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 59840 3942400 585 chr1 5120
     *          10240 chr1.1 5 1024 1024 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 59900 3904400 585
     *          chr1 10240 15360 chr1.2 5 1024 2048 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 55120
     *          3411200 585 chr1 15360 20480 chr1.3 5 1024 3072 /gbdb/hg18/wib/gc5Base.wib 0 100
     *          1024 49900 3078800 585 chr1 20480 25600 chr1.4 5 1024 4096
     *          /gbdb/hg18/wib/gc5Base.wib 0 100 1024 47600 2682400 A sample gcModel file (the
     *          output) look like this: Name Chr Position GC rs4961 4 2876505 48.6531211131841
     *          rs3025091 11 102219838 38.1080923507463 rs3092963 3 46371942 44.4687694340796
     *          rs3825776 15 56534122 40.7894123134328 rs17548390 6 3030512 45.3604050062189
     *          rs2276302 11 113355350 43.8200598569652 snpFile or pfb file (Population B Allele
     *          Frequency) is an additional data file that is required by the corresponding feature
     *          in PennCNV but not here. It looks like the this: Name Chr Position PFB rs1000113 5
     *          150220269 0.564615751221256 rs1000115 9 112834321 0.565931333264192 rs10001190 4
     *          6335534 0.5668604380025 rs10002186 4 38517993 0.57141752993563 rs10002743 4 6327482
     *          0.567557695424774
     */
    public static void gcModel(Project proj, String inputGcBaseFullPath,
                               String outputGcModelFullPath, int numwindow) {
      MarkerSetInfo markerSet;
      String[] markerNames;
      byte[] chrs;
      int[] positions;

      BufferedReader reader;
      String[] line;
      PrintWriter writer;
      byte curchr;
      int curstart, curend, curcount;
      boolean skip = false;
      double cursum;
      Hashtable<Byte, Byte> seen_chr = new Hashtable<>();

      int[] snp_count;
      double[] snp_sum;
      byte prechr = -1;
      int prestart = -1;
      int chr_index;
      Logger log;

      log = proj.getLog();

      // generate or load SnpFile (pbf file or Population B Allele Frequency)
      markerSet = proj.getMarkerSet();
      markerNames = markerSet.getMarkerNames();
      chrs = markerSet.getChrs(); // to be used only in the SnpFile (.pdf file) block. But then in the
                                 // outputFile.
      positions = markerSet.getPositions(); // to be used only in the SnpFile (.pdf file) block. But
                                           // then in the GcFile block.

      // How do we know whether "numwindow==null" ???
      if (numwindow == 0) {
        numwindow = 100;
      }

      snp_count = new int[markerNames.length];
      snp_sum = new double[markerNames.length];
      chr_index = 0;

      // load gcFile
      try {

        PSF.checkInterrupted();

        // If an invalid gc base path was given, try using the resources
        if (!new File(inputGcBaseFullPath).exists()) {
          Resource r = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), log).getModelBase();
          if (r.isAvailable()) {
            inputGcBaseFullPath = r.getAbsolute();
          }
        }
        reader = Files.getAppropriateReader(inputGcBaseFullPath);
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          curchr = Positions.chromosomeNumber(line[1]);
          curstart = Integer.parseInt(line[2]);
          curend = Integer.parseInt(line[3]);
          curcount = Integer.parseInt(line[11]);
          cursum = Double.parseDouble(line[12]);

          // For each record in gcFile, scan snpFile for the range of +/- numwindow*5120
          if (curchr > 0) {
            if (curchr == prechr) {
              if (curstart < prestart) {
                log.reportError("Error in gcFile: a record in chr" + curchr + " has position "
                                + curstart + ", less then the previous position $prestart");
                reader.close();
                return;
              }
            } else if (seen_chr.containsKey(curchr)) {
              log.reportError("Error in gcFile: rows of the same chromosome must be adjacent. But now chr"
                              + curchr + " occur multiple times in non-continuous segment of the "
                              + inputGcBaseFullPath + ": at " + curchr + ":" + curstart);
              reader.close();
              return;
            } else {
              seen_chr.put(curchr, (byte) 1);

              skip = false;
              if (chrs[chr_index] > curchr) {
                chr_index = 0;
              }
              while (chr_index < (chrs.length) && chrs[chr_index] != curchr) {
                chr_index++;
              }
            }

            if (!(curchr == prechr && skip)) {
              while (chr_index < chrs.length && chrs[chr_index] == curchr
                     && positions[chr_index] < Math.max(curstart - numwindow * 5120, 0)) {
                chr_index++;
              }
              if (chr_index >= chrs.length) {
                chr_index = 0;
                skip = true;
              } else if (chrs[chr_index] != curchr) {
                skip = true;
              } else {
                for (int i = chr_index; i < snp_count.length && chrs[i] == curchr
                                        && positions[i] <= (curend + numwindow * 5120); i++) {
                  snp_count[i] += curcount;
                  snp_sum[i] += cursum;
                }
              }
            }

            prestart = curstart;
            prechr = curchr;
          }
        }
        reader.close();
      } catch (Exception e) {
        log.reportError("Error reading from '" + inputGcBaseFullPath + "'");
        log.reportException(e);
        return;
      }

      PSF.checkInterrupted();
      // load pfb file or generate it

      // output the result
      try {
        writer = Files.openAppropriateWriter(outputGcModelFullPath);
        writer.println("Name\tChr\tPosition\tGC");
        for (int i = 0; i < markerNames.length; i++) {
          // writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]+"\t"+(snp_count[i]==0?(snp_sum[i]==0?0:"err"):(snp_sum[i]/snp_count[i])));
          writer.println(markerNames[i] + "\t"
                         + (chrs[i] < 23 ? chrs[i]
                                         : (chrs[i] == 23 ? "X"
                                                          : (chrs[i] == 24 ? "Y"
                                                                           : (chrs[i] == 25 ? "XY"
                                                                                            : (chrs[i] == 26 ? "M"
                                                                                                             : "Un")))))
                         + "\t" + positions[i] + "\t"
                         + (snp_count[i] == 0 ? (snp_sum[i] == 0 ? 0 : "err")
                                              : (snp_sum[i] / snp_count[i])));
        }
        writer.close();
        log.report("Generated population GC Model " + outputGcModelFullPath);
      } catch (Exception e) {
        log.reportError("Error writing to '" + outputGcModelFullPath + "'");
        log.reportException(e);
      }
    }

    public static GcModel loadSerial(String fullPathToGcSer) {
      return (GcModel) SerializedFiles.readSerial(fullPathToGcSer);
    }

    /**
     * This generates a gcModel file that takes the bin size to be explicitly the window around the
     * snp,
     */
    public static GcModel generateSnpWindowModel(Project proj, int snpWindow) {
      String refGenome = proj.getReferenceGenomeFASTAFilename();
      if (!Files.exists(refGenome)) {
        proj.getLog().reportFileNotFound(refGenome);
        return null;
      } else {
        ReferenceGenome referenceGenome = new ReferenceGenome(refGenome, proj.getLog());
        Map<String, Integer> indices = proj.getMarkerIndices();
        MarkerSetInfo markerSet = proj.getMarkerSet();
        double[] gcs = new double[markerSet.getMarkerNames().length];
        for (int i = 0; i < markerSet.getMarkerNames().length; i++) {
          if (i % 10000 == 0) {
            proj.getLog()
                .reportTimeInfo("Generating snp window gc model for window of " + snpWindow + " ("
                                + (i + 1) + " of " + markerSet.getMarkerNames().length + ")");
          }
          Segment seg = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i],
                                    markerSet.getPositions()[i]).getBufferedSegment(snpWindow);
          double gc = 100 * referenceGenome.getGCContentFor(seg, snpWindow > 100000);// not sure
                                                                                     // about the
                                                                                     // optimal
                                                                                     // query size
          gcs[i] = gc;
          if (Double.isNaN(gc) && seg.getChr() > 0) {
            proj.getLog()
                .reportError("Invalid gc content returned for query " + seg.getChromosomeUCSC()
                             + " and marker " + markerSet.getMarkerNames()[i]);
            return null;
          }
        }

        return new GcModel(markerSet.getMarkerNames(), markerSet.getChrs(),
                           markerSet.getPositions(), gcs, indices, proj.getLog());
      }

    }

    public static GcModel populateFromFile(String fullPathToGcModel, boolean verbose, Logger log) {
      ArrayList<String> markers = new ArrayList<>();
      ArrayList<Byte> chrs = new ArrayList<>();
      ArrayList<Integer> positions = new ArrayList<>();
      ArrayList<Double> gcs = new ArrayList<>();
      Hashtable<String, Integer> index = new Hashtable<>();
      String fullPathToGcSer = ext.rootOf(fullPathToGcModel, false) + ".gcmodel.ser";
      if (Files.exists(fullPathToGcSer)) {
        log.report("Info - loading gc model file " + fullPathToGcSer);
        GcModel model = null;
        try {
          model = loadSerial(fullPathToGcSer);
          model.getChrs();
          // having trouble with linux->windows SID transfer when gzipped
        } catch (Exception e) {}

        if (model != null) {
          return model;
        }
      }
      if (!Files.exists(fullPathToGcModel)) {
        log.reportError("Error - could not find gc model file " + fullPathToGcModel);
        return null;
      }
      log.report("Info - parsing gc model file " + fullPathToGcModel);

      int[] indices = ext.indexFactors(Files.getHeaderOfFile(fullPathToGcModel, log), GC_HEADER,
                                       true);
      if (ArrayUtils.countIf(indices, -1) > 0) {
        log.reportError("Error - could not find correct header for gc model file "
                        + fullPathToGcModel);
        log.reportError("		 - header must be:" + ArrayUtils.toStr(GC_HEADER));
        return null;
      } else {
        BufferedReader reader;
        try {
          reader = Files.getAppropriateReader(fullPathToGcModel);
          String[] line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);// header
          int lineNum = 0;
          while (reader.ready()) {
            line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
            try {
              String marker = line[indices[0]];
              byte chr = Positions.chromosomeNumber(line[indices[1]], log);
              int pos = Integer.parseInt(line[indices[2]]);
              double gc = Double.parseDouble(line[indices[3]]);
              markers.add(marker);
              chrs.add(chr);
              positions.add(pos);
              gcs.add(gc);
              index.put(marker, lineNum);
              lineNum++;
            } catch (NumberFormatException nfe) {
              if (verbose) {
                log.reportError("Error - found invalid number format on line "
                                + ArrayUtils.toStr(line) + " , skipping");
              }
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          log.reportError("Error: file \"" + fullPathToGcModel
                          + "\" not found in current directory");
          return null;
        } catch (IOException ioe) {
          log.reportError("Error reading file \"" + fullPathToGcModel + "\"");
          return null;
        }
        if (markers.size() == 0) {
          log.reportError("Error - did not find any valid markers in gc model file "
                          + fullPathToGcModel);
          return null;
        } else {
          log.report("Info - loaded " + markers.size() + " markers from gc model file "
                     + fullPathToGcModel);
          GcModel gcModel = new GcModel(ArrayUtils.toStringArray(markers), Bytes.toArray(chrs),
                                        Ints.toArray(positions), Doubles.toArray(gcs), index, log);
          gcModel.Serialize(fullPathToGcSer);
          return gcModel;
        }
      }
    }
  }

  /**
   * Similar to {@link GcAdjustor#getComputedAdjustor()} except it takes a {@link Sample} object as
   * input
   */
  public static GcAdjustor getComputedAdjustor(Project proj, Sample sample, GcModel gcModel,
                                               GC_CORRECTION_METHOD correctionMethod,
                                               boolean computePrior, boolean computePost,
                                               boolean verbose) {
    return getComputedAdjustor(proj, sample.markerLRRMap(proj.getMarkerSet()), gcModel,
                               correctionMethod, computePrior, computePost, verbose);
  }

  public static GcAdjustor getComputedAdjustor(Project proj, Map<Marker, Double> markerIntensities,
                                               GcModel gcModel,
                                               GC_CORRECTION_METHOD correctionMethod,
                                               boolean computePrior, boolean computePost,
                                               boolean verbose) {
    return getComputedAdjustor(proj, null, null, markerIntensities, gcModel, correctionMethod,
                               computePrior, computePost, verbose);
  }

  /**
   * @param proj
   * @param sample the current sample, only neccesary if a {@link GcAdjustorParameter} is provided
   * @param gcParameters if a valid {@link GcAdjustorParameter} is provided, the regression model
   *          will be skipped and the betas from this object will be used
   * @param markerIntensities float[] that is converted to double[] for correction
   * @param gcModel a valid model
   * @param computePrior compute WF and GCWF prior to gc correction
   * @param computePost compute WF and GCWF post gc correction
   * @param verbose reports things akin to PennCNV
   * @param pennCNVGCWF use the PennCNV default chr 11 bins for GCWF calculation
   * @return a {@link GcAdjustor} object with corrected intensities and qc metrics
   */
  public static GcAdjustor getComputedAdjustor(Project proj, String sample,
                                               GcAdjustorParameter gcParameters,
                                               Map<Marker, Double> markerIntensities,
                                               GcModel gcModel,
                                               GC_CORRECTION_METHOD correctionMethod,
                                               boolean computePrior, boolean computePost,
                                               boolean verbose) {
    GCAdjustorBuilder builder = new GCAdjustorBuilder();
    builder.correctionMethod(correctionMethod);
    return getComputedAdjustor(proj, builder, sample, gcParameters, markerIntensities, gcModel,
                               computePrior, computePost, verbose);
  }

  public static GcAdjustor getComputedAdjustor(Project proj, GCAdjustorBuilder builder,
                                               String sample, GcAdjustorParameter gcParameters,
                                               Map<Marker, Double> markerIntensities,
                                               GcModel gcModel, boolean computePrior,
                                               boolean computePost, boolean verbose) {
    GcAdjustor gcAdjustor = builder.build(proj, gcModel, markerIntensities);
    builder.verbose(verbose);
    gcAdjustor.correctIntensities(sample, gcParameters);
    gcAdjustor.computeQCMetrics(computePrior, computePost);
    return gcAdjustor;
  }

  public static void test(Project proj, String fullPathToGcModel,
                          String fullPathToFileOfTestSamples) {
    GcModel gcModel = GcModel.populateFromFile(fullPathToGcModel, false, proj.getLog());
    String[] samplesToTest = ArrayUtils.subArray(proj.getSamples(),
                                                 proj.getSamplesToInclude(fullPathToFileOfTestSamples));
    String fileTest = proj.PROJECT_DIRECTORY.getValue() + "testGCWF.txt";
    try {
      PrintWriter writer = Files.openAppropriateWriter(fileTest);
      for (String element : samplesToTest) {
        long time = System.currentTimeMillis();
        Sample samp = proj.getFullSampleFromRandomAccessFile(element);
        proj.getLog().report("Testing sample " + samp.getSampleName());
        GCAdjustorBuilder builder = new GCAdjustorBuilder();
        GcAdjustor gcAdjusterNew = builder.build(proj, gcModel,
                                                 samp.markerLRRMap(proj.getMarkerSet()));
        gcAdjusterNew.correctIntensities();
        gcAdjusterNew.computeQCMetrics(true, true);
        writer.println(samp.getSampleName() + "\t" + gcAdjusterNew.getQCString() + "\t"
                       + ext.getTimeElapsed(time));
      }
      writer.close();
    } catch (Exception e) {
      proj.getLog().reportError("Error writing to " + fileTest);
      proj.getLog().reportException(e);
    }
  }

  public static void main(String[] args) {
    CLI cli = new CLI(GcModel.class);

    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);
    cli.addArg(CLI.ARG_LOG, CLI.DESC_LOG, false);
    cli.addArg(GC_BASE_FILE, "GC5Base file");

    cli.parseWithExit(args);

    String filename = cli.get(CLI.ARG_PROJ);
    String logfile = cli.has(CLI.ARG_LOG) ? cli.get(CLI.ARG_LOG) : null;
    String gc5base = cli.get(GC_BASE_FILE);

    Project proj = new Project(filename, logfile);
    GcModel.gcModel(proj, gc5base, proj.GC_MODEL_FILENAME.getValue(), 100);
  }

  public static class GCAdjustorBuilder {

    private double minimumAutosomalGC = DEFAULT_MINAUTOSOMALGC;
    private double maximimumAutosomalGC = DEFAULT_MAXAUTOSOMALGC;
    private double minIntensity = DEFAULT_MIN_DATA_VALUE;
    private double maxIntensity = DEFAULT_MAX_DATA_VALUE;
    private int regressionDistance = DEFAULT_REGRESSION_DISTANCE[0];
    private int numSnpMAD = DEFAULT_NUM_SNP_MAD;
    private int skipPerChr = DEFAULT_SKIP_PER_CHR[0];
    private Set<Marker> markerMask = null;
    private boolean verbose = false;
    private GC_CORRECTION_METHOD correctionMethod = GC_CORRECTION_METHOD.GENVISIS_GC;

    public int getRegressionDistance() {
      return regressionDistance;
    }

    public int getNumSnpMAD() {
      return numSnpMAD;
    }

    public GCAdjustorBuilder minimumAutosomalGC(double minimumAutosomalGC) {
      this.minimumAutosomalGC = minimumAutosomalGC;
      return this;
    }

    public GCAdjustorBuilder numSnpMAD(int numSnpMAD) {
      this.numSnpMAD = numSnpMAD;
      return this;
    }

    public GCAdjustorBuilder maximimumAutosomalGC(double maximimumAutosomalGC) {
      this.maximimumAutosomalGC = maximimumAutosomalGC;
      return this;
    }

    public GCAdjustorBuilder minIntensity(double minIntensity) {
      this.minIntensity = minIntensity;
      return this;
    }

    public GCAdjustorBuilder maxIntensity(double maxIntensity) {
      this.maxIntensity = maxIntensity;
      return this;
    }

    public GCAdjustorBuilder regressionDistance(int regressionDistance) {
      this.regressionDistance = regressionDistance;
      return this;
    }

    public GCAdjustorBuilder skipPerChr(int skipPerChr) {
      this.skipPerChr = skipPerChr;
      return this;
    }

    public GCAdjustorBuilder markerMask(Set<Marker> markerMask) {
      this.markerMask = markerMask;
      return this;
    }

    public GCAdjustorBuilder verbose(boolean verbose) {
      this.verbose = verbose;
      return this;
    }

    public GCAdjustorBuilder correctionMethod(GC_CORRECTION_METHOD correctionMethod) {
      this.correctionMethod = correctionMethod;
      return this;
    }

    public GcAdjustor build(Project proj, GcModel gcModel, Map<Marker, Double> markerIntensities) {
      return new GcAdjustor(this, proj, gcModel, markerIntensities);
    }

    public GCAdjustorBuilder() {

    }

    public GCAdjustorBuilder(GCAdjustorBuilder builder) {
      numSnpMAD = builder.numSnpMAD;
      minimumAutosomalGC = builder.minimumAutosomalGC;
      maximimumAutosomalGC = builder.maximimumAutosomalGC;
      minIntensity = builder.minIntensity;
      maxIntensity = builder.maxIntensity;
      regressionDistance = builder.regressionDistance;
      skipPerChr = builder.skipPerChr;
      markerMask = builder.markerMask;
      verbose = builder.verbose;
      correctionMethod = builder.correctionMethod;
    }
  }

  private GcAdjustor(GCAdjustorBuilder builder, Project proj, GcModel gcModel,
                     Map<Marker, Double> markerIntensities) {
    numSnpMAD = builder.numSnpMAD;
    minimumAutosomalGC = builder.minimumAutosomalGC;
    maximimumAutosomalGC = builder.maximimumAutosomalGC;
    minIntensity = builder.minIntensity;
    maxIntensity = builder.maxIntensity;
    regressionDistance = builder.regressionDistance;
    skipPerChr = builder.skipPerChr;
    markerMask = builder.markerMask;
    verbose = builder.verbose;
    correctionMethod = builder.correctionMethod;
    this.proj = proj;
    this.markerIntensities = markerIntensities;
    if (gcModel == null) {
      fail = true;
      proj.getLog().reportError("Error - a gcModel was not supplied");
      markerGcs = null;
    } else {
      markerGcs = gcModel.getMarkerGcMap(proj);
    }
    fail = false;
    wfPrior = Double.NaN;
    wfPost = Double.NaN;
    gcwfPrior = Double.NaN;
    gcwfPost = Double.NaN;
    populateCurrentData();// Initialize everything needed
  }
}

// /**
// * @return an array with WF and GCWF prior to and post gc correction
// */
// public static double[] getQCMetrics(Project proj, Sample sample, GcModel gcModel, boolean
// pennCNVGCWF, boolean verbose) {
// double[] waves = new double[4];// stores WF prior, WF post, GCWF prior, GCWF post
// Arrays.fill(waves, Double.NaN);
// if (gcModel != null) {
// GcAdjustor gcAdjustor = getComputedAdjustor(proj, sample, gcModel, pennCNVGCWF, true, true,
// verbose);
// waves[0] = gcAdjustor.getWfPrior();
// waves[1] = gcAdjustor.getWfPost();
// waves[2] = gcAdjustor.getGcwfPrior();
// waves[3] = gcAdjustor.getGcwfPost();
// }
// return waves;
// }
