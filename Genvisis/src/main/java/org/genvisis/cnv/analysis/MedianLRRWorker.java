package org.genvisis.cnv.analysis;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.cnv.qc.CNVBDeviation;
import org.genvisis.cnv.qc.CNVBMAF;
import org.genvisis.cnv.qc.CNVBMAF.PoplulationBAFs;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.LeastSquares.LS_TYPE;

public class MedianLRRWorker extends SwingWorker<String, Integer> {

  // for storing region information
  private class MarkerRegion {

    private final String regionName;
    private final int regionType; // 0 = UCSC ,1 =markerOnly
    private final ArrayList<Integer> markerIndex;
    private final ArrayList<String> markersInRegion;
    private final Hashtable<String, Boolean> inRegion;
    private final Hashtable<String, Integer> index;

    private MarkerRegion(String regionName, int regionID, int regionType) {
      this.regionName = regionName;
      this.regionType = regionType;
      markerIndex = new ArrayList<Integer>();
      // stores markers in the region
      markersInRegion = new ArrayList<String>();
      // for assigning in region
      inRegion = new Hashtable<String, Boolean>();
      // for assigning index
      index = new Hashtable<String, Integer>();
    }

    private void addMarker(String marker) {
      markersInRegion.add(marker);
    }

    private void addMarkerIndex(int position) {
      markerIndex.add(position);
    }

    private void addMarkers(String[] input) {
      for (int i = MARKER_REGION_START_OF_MARKERS; i < input.length; i++) {
        if (hash.containsKey(input[i])) {
          computelog.report(input[i] + " was filtered out");
        } else {
          processTracker[0]++;
          addMarker(input[i]);
          markerInRegion(input[i]);
        }
      }
      checkNumMarkers();
    }

    private void addUCSCMarkers(String UCSCLine) {
      Segment seg = new Segment(UCSCLine);
      for (int i = 0; i < positions.length; i++) {
        if (chrs[i] == seg.getChr() && positions[i] >= seg.getStart()
            && positions[i] <= seg.getStop()) {
          if (hash.containsKey(markerNames[i])) {
            computelog.report(markerNames[i] + " was filtered out");
          } else {
            processTracker[0]++;
            addMarker(markerNames[i]);
            markerInRegion(markerNames[i]);
            addMarkerIndex(i);
            assignIndex(i, markerNames[i]);
          }
        }
      }
      checkNumMarkers();
    }

    private void assignIndex(int anindex, String markerName) {
      if (inRegion.containsKey(markerName)) {
        index.put(markerName, anindex);
      }
    }

    private void assignPositions() {
      for (int i = 0; i < markersInRegion.size(); i++) {
        if (index.containsKey(markersInRegion.get(i))) {
          markerIndex.add(index.get(markersInRegion.get(i)));
        } else {
          computelog.reportError("Error -could not find marker " + markersInRegion.get(i)
                                 + " in Marker Lookup");
        }
      }
    }

    private void checkNumMarkers() {
      if (markersInRegion.size() < 1) {
        String Error = "Error - All markers were filtered out of region " + regionName;
        computelog.reportError(Error);
        warnAndCancel(Error);
      }
    }

    private Hashtable<String, Integer> getIndex() {
      return index;
    }

    public ArrayList<Integer> getMarkerIndex() {
      return markerIndex;
    }

    private String getRegionName() {
      return regionName;
    }

    private int getRegionType() {
      return regionType;
    }

    public Segment getSegment(byte chr) {
      int start = positions[markerIndex.get(0)];
      int stop = positions[markerIndex.get(0)];
      for (int i = 0; i < markerIndex.size(); i++) {
        if (positions[markerIndex.get(i)] < start) {
          start = positions[markerIndex.get(i)];
        }
        if (positions[markerIndex.get(i)] > stop) {
          stop = positions[markerIndex.get(i)];
        }

      }
      return new Segment(chr, start, stop);
    }

    public boolean isOneChromosome() {
      boolean same = true;
      byte chr = chrs[markerIndex.get(0)];
      for (int i = 0; i < markerIndex.size(); i++) {
        if (chrs[markerIndex.get(i)] != chr) {
          same = false;
        }
      }
      return same;
    }

    private void markerInRegion(String marker) {
      inRegion.put(marker, true);
    }

    private String[] returnMarkers() {
      return markersInRegion.toArray(new String[markersInRegion.size()]);
    }
  }
  private class RegionResults {

    // stored as region, sample
    private final float[][] regionMedianValues;
    private final float[][] regionMADValues;
    private final float[][] regionBDeviationAll;
    private final float[][] regionBDeviationHet;
    private final float[][] regionBmafMetric;
    private final float[][] regionBmafMetricPenalty;
    private final float[][] regionPercentHet;

    // stored as region,sample,cnvFile
    private int[][][] cnvFileCNs;
    private String[] cnvFileNames;

    public RegionResults(float[][] regionMedianValues, float[][] regionMADValues,
                         float[][] regionBDeviationAll, float[][] regionBDeviationHet,
                         float[][] regionBmafMetric, float[][] regionBmafMetricPenalty,
                         float[][] regionPercentHet) {
      super();
      this.regionMedianValues = regionMedianValues;
      this.regionMADValues = regionMADValues;
      this.regionBDeviationAll = regionBDeviationAll;
      this.regionBDeviationHet = regionBDeviationHet;
      this.regionBmafMetric = regionBmafMetric;
      this.regionBmafMetricPenalty = regionBmafMetricPenalty;
      this.regionPercentHet = regionPercentHet;
      cnvFileCNs = null;
      cnvFileNames = null;
    }

    public float getBDeviationAllAt(int region, int sample) {
      if (regionBDeviationAll != null) {
        return regionBDeviationAll[region][sample];
      } else {
        return Float.NaN;
      }
    }

    public float getBDeviationHetAt(int region, int sample) {
      if (regionBDeviationHet != null) {
        return regionBDeviationHet[region][sample];
      } else {
        return Float.NaN;
      }
    }

    public float getBmafMetricAt(int region, int sample) {
      if (regionBmafMetric != null) {
        return regionBmafMetric[region][sample];
      } else {
        return Float.NaN;
      }
    }

    public float getBmafMetricPenaltyAt(int region, int sample) {
      if (regionBmafMetricPenalty != null) {
        return regionBmafMetricPenalty[region][sample];
      } else {
        return Float.NaN;
      }
    }

    public float getCNAt(int region, int sample, int cnvFile) {
      return cnvFileCNs[region][sample][cnvFile];
    }

    public String[] getCnvFileNames() {
      return cnvFileNames;
    }

    public float getMADAt(int region, int sample) {
      return regionMADValues[region][sample];
    }

    public float getMedianAt(int region, int sample) {
      return regionMedianValues[region][sample];
    }

    public float getPercentHetAt(int region, int sample) {
      if (regionPercentHet != null) {
        return regionPercentHet[region][sample];
      } else {
        return Float.NaN;
      }
    }

    public void setCnvFileCNs(int[][][] cnvFileCNs) {
      this.cnvFileCNs = cnvFileCNs;
    }

    public void setCnvFileNames(String[] cnvFileNames) {
      this.cnvFileNames = cnvFileNames;
    }

  }

  private static final String[] FILE_PREFIXES = {"LRR_MEDIAN_", "MarkersIn_"};
  private static final String[] FILE_EXT = {".xln"};
  private static final String MARKER_REGION_PREFIX = "probeset_id";
  private static final String MARKER_REGION_DELIMITER = ";";
  private static final String MARKER_REGION_REGEX = "[%0]";
  private static final int MARKER_REGION_START_OF_MARKERS = 2;
  private static final double[] QUANTILES = {0.5};
  private static final String[] MEDIAN_WORKER_JOBS = {"Parsing and intitializing Regions",
                                                      "Computing Median Log R Ratios for " + MARKER_REGION_REGEX,
                                                      "Waiting for data to Load for Region " + MARKER_REGION_REGEX,
                                                      "Creating Output Files",
                                                      "Assigning cnv copyNumber for "};
  public static final String[] CLASSES_TO_DUMP = {"FID", "IID"};

  private static final String[] MARKER_REGION_RESULTS_SUFFIX = {"MEDIAN", "MAD", "BDeviation_All",
                                                                "BDeviation_Het", "BMAF_Metric",
                                                                "BMAF_Metric_Penalty",
                                                                "Proportion_Het"};
  private static final String CNV_CLASSES = "COPYNUMBER" + MARKER_REGION_REGEX
                                            + ";5=CN0;1=CN1;2=CN2;3=CN3;4=CN4";

  private final Project proj;
  private final String[] input;
  private final String outputBase;
  private final int transformationType;
  private final int scope;
  private final Logger computelog; // considered using the project's default log or using the
                                  // proj.setLog() command, but nice to have the region specific
                                  // log named the same
  private final String[] samples, markerNames;
  private final byte[] chrs;

  private final int[] positions;

  private final Hashtable<String, String> hash;

  private final MarkerSetInfo markerSet;

  private final SampleList sampleList;

  private final int[] processTracker;

  // total,processed,updateat,tracker,scale;

  private final boolean[] transChrs;

  private final boolean recomputeLRR, correctXY, correctLRR;
  private final CHROMOSOME_X_STRATEGY strategy;
  private final boolean homozygousOnly;

  private final JProgressBar progressBar;

  public MedianLRRWorker(Project proj, String[] input, int transformationType, int scope,
                         String outputBase, JProgressBar jProgressBar, boolean recomputeLRR,
                         boolean correctLRR, boolean correctXY, boolean homozygousOnly,
                         CHROMOSOME_X_STRATEGY strategy, Logger log) {
    this.proj = proj;
    this.input = input;
    this.outputBase = outputBase;
    this.transformationType = transformationType;
    this.scope = scope;
    this.strategy = strategy;
    progressBar = jProgressBar == null ? new JProgressBar() : jProgressBar;
    addPropertyChangeListener(new PropertyChangeListener() {

      @Override
      public void propertyChange(PropertyChangeEvent evt) {
        if ("progress".equals(evt.getPropertyName())) {
          progressBar.setValue((Integer) evt.getNewValue());
        }
      }
    });

    computelog = log == null ? new Logger(proj.PROJECT_DIRECTORY.getValue() + outputBase + ".log")
                             : log;
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    sampleList = proj.getSampleList();
    hash = proj.getFilteredHash();
    samples = sampleList.getSamples();
    // total ,processed;
    processTracker = new int[2];
    transChrs = ArrayUtils.booleanArray(27, false);
    this.recomputeLRR = recomputeLRR;
    this.correctLRR = correctLRR;
    this.correctXY = correctXY;
    this.homozygousOnly = homozygousOnly;// only valid for marker-based approaches
  }

  // only if transforming by chromosome
  private void assignChr(MarkerRegion[] markerRegions) {
    for (int i = 0; i < markerRegions.length; i++) {
      String[] regionMarkers = markerRegions[i].returnMarkers();
      for (int k = 0; k < regionMarkers.length; k++) {
        if (!markerRegions[i].getIndex().containsKey(regionMarkers[k])) {
          computelog.reportError("Error - could not find " + regionMarkers[k]
                                 + " in markerLookup file , was this probeset analyzed?");
        } else {
          int index = markerRegions[i].getIndex().get(regionMarkers[k]);
          transChrs[chrs[index]] = true;
        }
      }
    }

  }

  // First Get Markers For Region if Region, else region =markers input;

  private void assignMarkerIndices(MarkerRegion[] markerRegions) {
    for (int i = 0; i < markerNames.length; i++) {
      for (MarkerRegion markerRegion : markerRegions) {
        // UCSC indices already Assigned;
        if (markerRegion.getRegionType() == 1) {
          markerRegion.assignIndex(i, markerNames[i]);
        }
      }
    }
    for (MarkerRegion markerRegion : markerRegions) {
      markerRegion.assignPositions();
    }
  }

  // For each region load lrrs from markerData and compute

  private void assignMarkerProgress() {
    progressBar.setMaximum(processTracker[0]);

  }

  private void assignSampleProgress() {
    progressBar.setMaximum(samples.length);
  }

  private boolean checkUCSCRegions(String region) {
    boolean valid = false;
    int[] newLocation = Positions.parseUCSClocation(region);
    if ((newLocation == null || newLocation.length != 3 || newLocation[0] < 0)
        || (newLocation[1] < 0) || (newLocation[2] < 0)) {
      valid = false;
    } else {
      valid = true;
    }
    return valid;
  }

  @Override
  protected String doInBackground() throws Exception {
    String fileNameToVisualize = "";
    System.out.println(proj.PROJECT_DIRECTORY.getValue() + outputBase);
    // FIXME It is NOT safe to make this call from an off-EDT thread. SwingWorkers are not intended
    // to directly interact with swing components. If a SwingWorker is being run on the EDT it is
    // being
    // used incorrectly.
    // Upcoming TaskManager API should make it easier to appropriately incorporate SwingWorkers
    // https://github.com/PankratzLab/Genvisis/issues/293
    progressBar.setValue(0);
    setProgress(0);
    newJob(MEDIAN_WORKER_JOBS[0]);
    MarkerRegion[] markerRegions = parseRegions();
    process(processTracker[0]);
    if (processTracker[0] < 1) {
      String Error = "Error - did not find any markers in the Region(s) of interest";
      computelog.reportError(Error);
      warnAndCancel(Error);
    } else {
      String job;
      process(0);
      RegionResults regionResults;
      if (transformationType == 0) {
        job = ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]",
                                 Transforms.TRANFORMATIONS[transformationType]);
        newJob(job);
        computelog.report(job);
        assignMarkerProgress();
        regionResults = getRawValueMedianForRegions(markerRegions);

        process(processTracker[0]);
        fileNameToVisualize = printResults(regionResults, markerRegions);
      } else {
        assignSampleProgress();
        job = ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]",
                                 Transforms.TRANFORMATIONS[transformationType] + " "
                                                                        + Transforms.SCOPES[scope]);
        computelog.report(job);
        newJob(job);
        regionResults = getNormalizedMedianForRegions(markerRegions);
        fileNameToVisualize = printResults(regionResults, markerRegions);
      }
      Thread.sleep(1000);
    }
    return fileNameToVisualize;
  }

  @Override
  protected void done() {
    try {
      setProgress(100);
      progressBar.setValue(100);
      progressBar.setStringPainted(false);
      progressBar.setVisible(false);
      // FIXME current understanding is that this method is a hook enqueued on the EDT after
      // doInBackground returns. This handling should be moved to whatever is relying on output from
      // this worker.
      get();
      JOptionPane.showMessageDialog(null, "Log R Ratio Summarization Complete");

    } catch (ExecutionException e) {
      computelog.reportError("Error - Could not Compute Median Log R Ratio Values");
      computelog.reportException(e);
    } catch (InterruptedException ie) {
      computelog.reportError("Error - Median Log R Ratio Computation was Interupted ");
      computelog.reportException(ie);
    } catch (CancellationException cce) {
      computelog.reportError("Error - Cancelling Median Log R Ratio Computation ");
      computelog.reportException(cce);
    }
  }

  private float getMAD(float[] regionLrrs, float median) {
    float[] diffs = new float[regionLrrs.length];
    for (int i = 0; i < regionLrrs.length; i++) {
      diffs[i] = Math.abs((regionLrrs[i] - median));
    }
    return ArrayUtils.quants(diffs, QUANTILES)[0];
  }

  private RegionResults getNormalizedMedianForRegions(MarkerRegion[] markerRegions) {
    int[][] indices;
    float[][] regionMedianValues = new float[markerRegions.length][samples.length];
    float[][] regionMADValues = new float[markerRegions.length][samples.length];
    assignMarkerIndices(markerRegions);
    computelog.report("Computing mean Log R ratios for:");
    // norm by genome
    if (scope == 1) {
      indices = new int[][] {ArrayUtils.arrayOfIndices(proj.getPartialSampleFromRandomAccessFile(samples[0])
                                                           .getLRRs().length)};
    } else {
      // only normalize chrs with markers in regions
      indices = markerSet.getIndicesByChr();
      assignChr(markerRegions);
    }
    long time = new Date().getTime();
    for (int i = 0; i < samples.length; i++) {
      newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", samples[i]));
      process(i + 1);
      if (i % 100 == 0) {
        time = new Date().getTime();
        computelog.report((i + 1) + " of " + samples.length + " (" + ext.getTimeElapsed(time)
                          + ")");
      }
      float[] lrrs = proj.getPartialSampleFromRandomAccessFile(samples[i]).getLRRs();
      // genome
      if (scope == 1) {
        lrrs = Transforms.transform(lrrs, transformationType, false, markerSet);
      }
      // default to norm by chromosome
      else {
        lrrs = Transforms.transform(lrrs, transformationType, indices, transChrs);
      }
      for (int j = 0; j < markerRegions.length; j++) {
        ArrayList<Integer> regionPositions = markerRegions[j].getMarkerIndex();
        ArrayList<Float> regionLrrs = new ArrayList<Float>();
        for (int k = 0; k < regionPositions.size(); k++) {
          if (Float.isNaN(lrrs[regionPositions.get(k)])) {
            continue;
          } else {
            regionLrrs.add(lrrs[regionPositions.get(k)]);
          }
        }

        regionMedianValues[j][i] = ArrayUtils.quants(toFloatArray(regionLrrs), QUANTILES)[0];
        regionMADValues[j][i] = getMAD(toFloatArray(regionLrrs), regionMedianValues[j][i]);
      }
    }
    return new RegionResults(regionMedianValues, regionMADValues, null, null, null, null, null);// TODO
                                                                                                // if
                                                                                                // b-deviation
                                                                                                // is
                                                                                                // effective
  }

  private RegionResults getRawValueMedianForRegions(MarkerRegion[] markerRegions) {
    float[][] rawMediansRegions = new float[markerRegions.length][samples.length];
    float[][] rawMADRegions = new float[markerRegions.length][samples.length];
    float[][] rawBDevRegionsAll = new float[markerRegions.length][samples.length];
    float[][] rawBDevRegionsHet = new float[markerRegions.length][samples.length];
    float[][] rawBDevBmafMetric = new float[markerRegions.length][samples.length];
    float[][] rawBDevBmafMetricPenalty = new float[markerRegions.length][samples.length];

    float[][] rawBDevPercentHet = new float[markerRegions.length][samples.length];

    assignMarkerIndices(markerRegions);
    for (int i = 0; i < markerRegions.length; i++) {
      // stores median ,MAD for each sample
      float[][] results = getRawValueResultsForRegion(markerRegions[i]);
      rawMediansRegions[i] = results[0];
      rawMADRegions[i] = results[1];
      rawBDevRegionsAll[i] = results[2];
      rawBDevRegionsHet[i] = results[3];
      rawBDevBmafMetric[i] = results[4];
      rawBDevBmafMetricPenalty[i] = results[5];
      rawBDevPercentHet[i] = results[6];
    }
    return new RegionResults(rawMediansRegions, rawMADRegions, rawBDevRegionsAll, rawBDevRegionsHet,
                             rawBDevBmafMetric, rawBDevBmafMetricPenalty, rawBDevPercentHet);
  }

  private float[][] getRawValueResultsForRegion(MarkerRegion markerRegion) {
    String[] regionMarkers = markerRegion.returnMarkers();
    float[][] sampleLrrs = new float[regionMarkers.length][samples.length];
    MarkerDataLoader markerDataLoader;
    boolean[] samplesToUse = null;
    PrincipalComponentsResiduals pcrs = null;
    Logger projLog;

    projLog = proj.getLog();
    if (correctLRR || correctXY || recomputeLRR) {// not truly necessary for recomputing LRR, but
                                                  // currently forcing it anyway
      if (!Files.exists(proj.INTENSITY_PC_FILENAME.getValue())) {
        String Error = "Error - could not find the file " + proj.INTENSITY_PC_FILENAME.getValue()
                       + ", cannot perform intensity correction";
        computelog.reportError(Error);
        warnAndCancel(Error);
      } else {
        pcrs = new PrincipalComponentsResiduals(proj, proj.INTENSITY_PC_FILENAME.getValue(), null,
                                                proj.getProperty(proj.INTENSITY_PC_NUM_COMPONENTS),
                                                false, 0, false, false, null);
        samplesToUse = proj.getSamplesToInclude(null);
        projLog.report("Info - using " + ArrayUtils.booleanArraySum(samplesToUse)
                       + " individuals that were not defined as excluded in "
                       + proj.SAMPLE_DATA_FILENAME.getValue() + " for correction clustering");
        projLog.reportTimeInfo("Correcting with " + proj.INTENSITY_PC_NUM_COMPONENTS.getValue()
                               + " components");
      }

    }
    regionMarkers = markerRegion.returnMarkers();
    sampleLrrs = new float[regionMarkers.length][samples.length];
    newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[2], "[%" + 0 + "]", markerRegion.getRegionName()));
    proj.setLog(computelog);
    markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, regionMarkers);
    proj.setLog(projLog);
    ClusterFilterCollection clusterFilterCollection = proj.getClusterFilterCollection();

    org.genvisis.cnv.qc.CNVBMAF.PoplulationBAFs pDeviation = new org.genvisis.cnv.qc.CNVBMAF.PoplulationBAFs(samples.length,
                                                                                                             CNVBDeviation.DEFAULT_INTENSITY_ONLY_FLAGS,
                                                                                                             CNVBDeviation.DEFAULT_GC_THRESHOLD);
    for (int i = 0; i < regionMarkers.length; i++) {
      if (i % 10 == 0 || correctXY) {
        newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]",
                                  markerRegion.getRegionName() + " (" + (i) + " of "
                                                                         + regionMarkers.length
                                                                         + ")"));
      }
      MarkerData markerData = markerDataLoader.requestMarkerData(i);
      processTracker[1]++;
      process(processTracker[1]);
      if (markerData.getFingerprint() != sampleList.getFingerprint()) {
        String Error = "Error - mismatched fingerprint for " + markerData.getMarkerName();
        computelog.reportError(Error);
        warnAndCancel(Error);
      }
      float[] lrrs = markerData.getLRRs();// default
      float[] bafs = markerData.getBAFs();
      float[] gcs = markerData.getGCs();
      byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection,
                                                               regionMarkers[i], 0, projLog);

      if (recomputeLRR || correctLRR || correctXY) {
        int numThreads = proj.getProperty(proj.NUM_THREADS);
        PrincipalComponentsIntensity pcIntensity = new PrincipalComponentsIntensity(pcrs,
                                                                                    markerData,
                                                                                    true, null,
                                                                                    samplesToUse, 1,
                                                                                    0,
                                                                                    clusterFilterCollection,
                                                                                    true,
                                                                                    LS_TYPE.REGULAR,
                                                                                    2, 5,
                                                                                    PrincipalComponentsIntensity.DEFAULT_RESID_STDV_FILTER,
                                                                                    PrincipalComponentsIntensity.DEFAULT_CORRECTION_RATIO,
                                                                                    numThreads,
                                                                                    false, null,
                                                                                    strategy);
        if (recomputeLRR && !correctLRR) {
          lrrs = pcIntensity.getCentroidCompute().getRecomputedLRR();
        } else if (correctLRR) {
          pcIntensity.correctLRRAt(proj.getProperty(proj.INTENSITY_PC_NUM_COMPONENTS));
          lrrs = pcIntensity.getCorrectedLRR();

          if (lrrs == null) {
            String Error = "Error - could not correct Log R Ratios for "
                           + markerData.getMarkerName();
            computelog.reportError(Error);
            warnAndCancel(Error);
          }
        } else if (correctXY) {
          pcIntensity.correctXYAt(proj.getProperty(proj.INTENSITY_PC_NUM_COMPONENTS));
          float[][] correctedLrrBafs = pcIntensity.getCorrectedIntensity(PrincipalComponentsIntensity.BAF_LRR_RETURN,
                                                                         true);
          lrrs = correctedLrrBafs[1];
          bafs = correctedLrrBafs[0];
        }
        if (lrrs == null) {
          String Error = "Error - could not correct Log R Ratios for " + markerData.getMarkerName();
          computelog.reportError(Error);
          warnAndCancel(Error);
        }
      }
      if (homozygousOnly) {
        for (int j = 0; j < genotypes.length; j++) {
          if (genotypes[j] != 0 && genotypes[j] != 2) {
            lrrs[j] = Float.NaN;// will no longer be included in the median
          }
        }
      }
      markerDataLoader.releaseIndex(i);
      pDeviation.add(regionMarkers[i], genotypes, bafs, gcs);
      for (int j = 0; j < samples.length; j++) {
        try {
          // marker;samples...
          sampleLrrs[i][j] = lrrs[j];
        } catch (ArrayIndexOutOfBoundsException aioobe) {
          computelog.report("" + i + "\t" + j);
          computelog.reportException(aioobe);
          System.exit(1);
        }
      }
    }
    process(processTracker[1]);
    return getSampleResults(sampleLrrs, pDeviation);
  }

  private int[][][] getRegionCNs(MarkerRegion[] markerRegions) {
    SampleData sampleData = proj.getSampleData(false);
    int[][][] cnvFileCNs = null;
    if (proj.CNV_FILENAMES.getValue() == null || proj.CNV_FILENAMES.getValue().length < 1) {
      computelog.report("Warning - no cnv file was found, not matching regions to cnvs");
    } else {
      // region,sample,cnvFile
      // CN are reported as 1=homzydel,2=hetrodel,3=none,4=hetrodup,5=homozydup;
      process(0);
      newJob("loading " + proj.CNV_FILENAMES.getValue().length + " cnv "
             + (proj.CNV_FILENAMES.getValue().length > 1 ? " files" : " file"));
      // TODO potential bug in old code - sets array length to length of String, not of value array
      cnvFileCNs = new int[markerRegions.length][samples.length][proj.CNV_FILENAMES.getValue().length];
      computelog.report("Info - assigning cnvs for " + proj.CNV_FILENAMES.getValue().length
                        + " cnv files");
      sampleData.loadCNVs(proj.CNV_FILENAMES.getValue(), proj.getLog());
      String[] cnvClasses = sampleData.getCnvClasses();
      computelog.report("Info - assigning cnvs for " + proj.CNV_FILENAMES.getValue().length);
      newJob(MEDIAN_WORKER_JOBS[4] + proj.CNV_FILENAMES.getValue().length + " cnv files");
      processTracker[1] = 0;
      for (int i = 0; i < markerRegions.length; i++) {
        assignSampleProgress();
        process(processTracker[1]);
        newJob(MEDIAN_WORKER_JOBS[4] + " region " + markerRegions[i].getRegionName());
        process(processTracker[1]);

        // Assuming that all markers are from the same chromosome...,but checking anyway
        byte chr = chrs[markerRegions[i].getMarkerIndex().get(0)];
        boolean oneChromosome = markerRegions[i].isOneChromosome();
        if (!oneChromosome) {
          computelog.report("Warning - region " + markerRegions[i].getRegionName()
                            + " spans more than one chromosome , setting all Copy number States to zero");
        }
        Segment region = markerRegions[i].getSegment(chr);
        for (int j = 0; j < samples.length; j++) {
          processTracker[1]++;
          if (j % 10 == 0) {
            newJob(MEDIAN_WORKER_JOBS[4] + " region " + markerRegions[i].getRegionName() + " (" + j
                   + " of " + samples.length + " samples)");
            process(processTracker[1]);
          }
          for (int k = 0; k < cnvClasses.length; k++) {
            // region,sample,cnvFile
            cnvFileCNs[i][j][k] = 2;
            CNVariant[] cnvs = sampleData.getIndiFromSampleHash(samples[j]).getCNVs(k, chr);
            if (cnvs != null && cnvs.length > 0 && oneChromosome) {
              for (CNVariant cnv : cnvs) {
                if (cnv.overlaps(region)) {
                  if (cnv.getCN() < 0) {
                    System.err.println("HAHAHA:");
                    System.exit(1);
                  }
                  cnvFileCNs[i][j][k] = cnv.getCN();
                }
              }
            }
          }
          process(j);
        }
      }
    }
    return cnvFileCNs;
  }

  private float[][] getSampleResults(float[][] lrrs, PoplulationBAFs pDeviation) {
    // store median, MAD,median B deviation all,median B deviation het
    float[][] results = new float[7][lrrs[0].length];
    ArrayList<ArrayList<Float>> sampleLRRS = new ArrayList<ArrayList<Float>>();
    for (int i = 0; i < lrrs[0].length; i++) {
      sampleLRRS.add(new ArrayList<Float>());
    }
    for (int i = 0; i < lrrs.length; i++) {
      for (int k = 0; k < lrrs[i].length; k++) {
        // Not Including NaNs
        // if marker i, sample k
        if (!Float.isNaN(lrrs[i][k])) {
          sampleLRRS.get(k).add((lrrs[i][k]));
        }
      }
    }
    pDeviation.summarize(CNVBMAF.DEFUALT_GENO_HET_PENALTY, CNVBMAF.DEFUALT_BAF_HET_PENALTY,
                         CNVBMAF.SUMMARY_BY_NUM_ALL_MARKERS, proj.getLog());
    double[] cnvbmafsNoPenalty = pDeviation.getCnvbmafsMetrics();
    pDeviation.summarize(CNVBMAF.DEFUALT_GENO_HET_PENALTY, CNVBMAF.DEFUALT_BAF_HET_PENALTY,
                         CNVBMAF.SUMMARY_BY_NUM_ALL_MARKERS_AVG_HET_PENALTY, proj.getLog());
    double[] cnvbmafsPenalty = pDeviation.getCnvbmafsMetrics();

    for (int i = 0; i < sampleLRRS.size(); i++) {
      if (sampleLRRS.get(i).size() > 0) {
        results[0][i] = ArrayUtils.quants(toFloatArray(sampleLRRS.get(i)), QUANTILES)[0];
        results[1][i] = getMAD(toFloatArray(sampleLRRS.get(i)), results[0][i]);
        results[2][i] = (float) pDeviation.getCnvbmafs()[i].getMedianBDeviationAll();
        results[3][i] = (float) pDeviation.getCnvbmafs()[i].getMedianBDeviationHet();
        results[4][i] = (float) cnvbmafsNoPenalty[i];
        results[5][i] = (float) cnvbmafsPenalty[i];
        results[6][i] = (float) pDeviation.getCnvbmafs()[i].getPercentHet();

      } else {
        proj.getLog().reportError("Warning - sample " + proj.getSamples()[i]
                                  + " did not have any data for a region, setting to NaN");
        results[0][i] = Float.NaN;
        results[1][i] = Float.NaN;
        results[2][i] = Float.NaN;
        results[3][i] = Float.NaN;
        results[4][i] = Float.NaN;
        results[5][i] = Float.NaN;
        results[6][i] = Float.NaN;

      }
    }
    return results;
  }

  // set ProgressBar Name
  protected void newJob(String job) {
    progressBar.setString(job);
  }

  private MarkerRegion[] parseRegions() {
    MarkerRegion[] regions = null;
    if (input != null) {
      regions = new MarkerRegion[input.length];
      for (int i = 0; i < input.length; i++) {
        process(Math.round((i + 1f) / input.length));
        if (input[i].split(MARKER_REGION_DELIMITER)[0].equals(MARKER_REGION_PREFIX)) {
          String[] regionMarkers = input[i].split(MARKER_REGION_DELIMITER);
          if (regionMarkers.length <= 1) {
            String Error = "Error - markers must be " + MARKER_REGION_DELIMITER
                           + " delimited for Marker Region " + input[i];
            computelog.reportError(Error);
            warnAndCancel(Error);
          }
          String regionName = regionMarkers[1];
          regions[i] = new MarkerRegion(regionName, i, 1);
          regions[i].addMarkers(regionMarkers);

        } else if (input[i].startsWith("chr") && checkUCSCRegions(input[i])) {
          regions[i] = new MarkerRegion(input[i], i, 0);
          regions[i].addUCSCMarkers(input[i]);
        } else {
          String Error = "Error - Improper Formatting for Input " + input[i]
                         + "\n Input must one per line of \n(1) UCSC Formatted\n or\n(2) "
                         + MARKER_REGION_PREFIX + MARKER_REGION_DELIMITER + "(Your Region Name)"
                         + MARKER_REGION_DELIMITER + "marker names(" + MARKER_REGION_DELIMITER
                         + ") delimited";
          computelog.reportError(Error);
          warnAndCancel(Error);
        }
      }
    } else {
      String Error = "Error - Zero Regions Were Entered!";
      computelog.reportError(Error);
      warnAndCancel(Error);
    }
    return regions;
  }

  // print median values
  private String printMedianLRRs(RegionResults regionResults, MarkerRegion[] markerRegions) {
    boolean reportCNVCN = false;
    String output = proj.PROJECT_DIRECTORY.getValue() + FILE_PREFIXES[0] + outputBase + FILE_EXT[0];
    Hashtable<String, String> hashSamps = HashVec.loadFileToHashString(proj.SAMPLE_DATA_FILENAME.getValue(),
                                                                       "DNA", CLASSES_TO_DUMP,
                                                                       "\t");
    int[][][] cnvFileCNs = getRegionCNs(markerRegions);
    String[] cnvFiles = proj.CNV_FILENAMES.getValue();
    if (cnvFileCNs != null) {
      computelog.report("Reporting cnvs");
      regionResults.setCnvFileCNs(cnvFileCNs);
      regionResults.setCnvFileNames(cnvFiles);
      reportCNVCN = true;
    }
    newJob(MEDIAN_WORKER_JOBS[3]);
    try {
      PrintWriter writer = Files.openAppropriateWriter(output);
      writer.print("Sample");
      for (String element : CLASSES_TO_DUMP) {
        writer.print("\t" + element.substring(element.lastIndexOf("=") + 1));
      }
      for (MarkerRegion markerRegion : markerRegions) {
        for (String element : MARKER_REGION_RESULTS_SUFFIX) {
          writer.print("\t" + element + "_" + markerRegion.getRegionName());
        }
        if (reportCNVCN) {
          for (int j = 0; j < cnvFiles.length; j++) {
            writer.print("\t"
                         + (ext.replaceAllWith(CNV_CLASSES, "[%" + 0
                                                            + "]",
                                               ext.rootOf(regionResults.getCnvFileNames()[j]) + "_"
                                                                   + markerRegion.getRegionName())));
          }
        }
      }
      writer.println();
      for (int i = 0; i < samples.length; i++) {
        writer.print(samples[i] + "\t"
                     + (hashSamps.containsKey(samples[i]) ? hashSamps.get(samples[i])
                                                          : ArrayUtils.stringArray(CLASSES_TO_DUMP.length,
                                                                                   ".")));
        for (int j = 0; j < markerRegions.length; j++) {
          writer.print("\t" + regionResults.getMedianAt(j, i) + "\t" + regionResults.getMADAt(j, i)
                       + "\t" + regionResults.getBDeviationAllAt(j, i) + "\t"
                       + regionResults.getBDeviationHetAt(j, i) + "\t"
                       + regionResults.getBmafMetricAt(j, i) + "\t"
                       + regionResults.getBmafMetricPenaltyAt(j, i) + "\t"
                       + regionResults.getPercentHetAt(j, i));
          if (reportCNVCN) {
            for (int k = 0; k < cnvFiles.length; k++) {
              writer.print("\t" + format(regionResults.getCNAt(j, i, k)));
            }
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      computelog.reportError("Error writing to Medain Log R Ratios to " + output);
      computelog.reportException(e);
    }
    return output;
  }

  // print markers in region
  private void printRegionMarkers(MarkerRegion[] markerRegions) {
    String output = proj.PROJECT_DIRECTORY.getValue() + FILE_PREFIXES[1] + outputBase + FILE_EXT[0];
    try {
      PrintWriter writer = Files.openAppropriateWriter(output);
      for (MarkerRegion markerRegion : markerRegions) {
        String[] regionMarkers = markerRegion.returnMarkers();
        if (markerRegion.getRegionType() == 0) {
          Segment seg = new Segment(markerRegion.getRegionName());
          writer.println(seg.getUCSClocation() + "\t" + regionMarkers.length + "\t" + seg.getChr()
                         + "\t" + seg.getStart() + "\t" + seg.getStop() + "\t"
                         + ArrayUtils.toStr(regionMarkers));
        } else if (markerRegion.getRegionType() == 1) {
          writer.println(markerRegion.getRegionName() + "\t" + regionMarkers.length
                         + "\tNA\tNA\tNA\t" + ArrayUtils.toStr(regionMarkers));
        } else {
          computelog.reportError("Error - Unknonwn Region Type " + markerRegion.getRegionType());
        }
      }
      writer.close();
    } catch (Exception e) {
      computelog.reportError("Error writing the list of marker names within the regions to "
                             + output);
      computelog.reportException(e);
    }
  }

  private String printResults(RegionResults regionResults, MarkerRegion[] markerRegions) {
    String output = "";
    printRegionMarkers(markerRegions);
    output = printMedianLRRs(regionResults, markerRegions);
    return output;
  }

  // update progressBar
  protected void process(Integer chunk) {
    progressBar.setValue(chunk);
  }

  protected void warnAndCancel(String message) {
    JOptionPane.showMessageDialog(null, message);
    cancel(true);
  }

  public static boolean checkExists(Project proj, String outputBase) {
    boolean exists = false;
    for (String element : FILE_PREFIXES) {
      for (String element2 : FILE_EXT) {
        if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + element + outputBase + element2)) {
          exists = true;
        }
      }
    }
    return exists;
  }

  // for access from command line
  public static void computeMedianLrrs(Project proj, String regionFileName, int transfromationType,
                                       int scope, CHROMOSOME_X_STRATEGY strategy,
                                       String outputBase) {
    Logger log = proj.getLog();
    String[] input = readToArray(proj.PROJECT_DIRECTORY.getValue() + regionFileName, log);
    MedianLRRWorker medianLRRWorker = new MedianLRRWorker(proj, input, transfromationType, scope,
                                                          outputBase, null, false, false, false,
                                                          false, strategy, log);
    log.report("Starting job for " + input.length + " regions");
    medianLRRWorker.execute();

  }

  public static String format(Number n) {
    NumberFormat format = DecimalFormat.getInstance();
    format.setRoundingMode(RoundingMode.FLOOR);
    format.setMinimumFractionDigits(0);
    format.setMaximumFractionDigits(2);
    return format.format(n);
  }

  private static String[] readToArray(String filename, Logger log) {

    ArrayList<String> lines = new ArrayList<String>();
    try {
      FileReader fileReader = new FileReader(filename);
      BufferedReader reader = new BufferedReader(fileReader);
      while (reader.ready()) {
        lines.add(reader.readLine().trim());
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error unable to find" + filename);

    } catch (IOException e) {
      log.reportError("Error reading file " + filename);
      log.reportException(e);
    }
    return lines.toArray(new String[lines.size()]);
  }

  // use this to eliminate NANs
  private static float[] toFloatArray(List<Float> al) {
    float[] d = new float[al.size()];
    for (int i = 0; i < al.size(); i++) {
      d[i] = al.get(i);
    }
    return d;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String regionFileName = "cnps.txt";
    String headless = "true";
    int transformationType = 0;
    int scope = 0;
    long time;
    boolean correctXY = false;
    String outputBase = Transforms.TRANFORMATIONS[transformationType];
    String logfile = outputBase + ".log";
    Logger log;
    Project proj;
    CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.BIOLOGICAL;

    String usage = "cnv.analysis.MedianLRRWorker requires 2 arguments\n";
    usage += "   (1) project properties filename (i.e. proj="
             + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n";
    usage += "   (2) filename of the regions (one per line) in UCSC format (chr8:25129632-25130278) \n";
    usage += "       OR:\n";
    usage += "       formatted as \"" + MARKER_REGION_PREFIX + MARKER_REGION_DELIMITER
             + "(Your Region Name)" + MARKER_REGION_DELIMITER + "marker name 1"
             + MARKER_REGION_DELIMITER + "marker name 2...\"" + MARKER_REGION_DELIMITER + "\n";
    usage += "       (i.e. regions=" + regionFileName + "(default))\n";
    usage += "       OPTIONAL:\n";
    usage += "   (3) transformation type (i.e. transform=0 (default, "
             + Transforms.TRANFORMATIONS[transformationType] + ")) \n";
    usage += "       transformations are: " + ArrayUtils.toStr(Transforms.TRANFORMATIONS) + "\n";
    usage += "   (4) scope of transformation (i.e. scope=0 (default))\n";
    usage += "       scopes are: " + ArrayUtils.toStr(Transforms.SCOPES) + "\n";
    usage += "   (5) base name of the output files (i.e out=" + outputBase + " (default))\n";
    usage += "   (6) name of the log file (i.e. log=" + logfile + "\n";
    usage += "   (7) run program in headless mode to quiet gui errors when X11 forwarding\n is un-available (i.e. headless=true (default));";
    usage += "   (8) correct data with principal components (must be defined by the properties file) (i.e. correctXY="
             + correctXY + " (default));";
    usage += "   (9) Chromosome X strategy.  Options include: ";
    usage += ArrayUtils.toStr(CHROMOSOME_X_STRATEGY.values(), ", ") + " (i.e. sexStrategy=";
    usage += strategy + " (default))\n";

    usage += "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("regions=")) {
        regionFileName = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("transform=")) {
        transformationType = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("scope=")) {
        scope = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outputBase = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("headless=")) {
        headless = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("correctXY=")) {
        correctXY = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("sexStrategy=")) {
        strategy = CHROMOSOME_X_STRATEGY.valueOf(ext.parseStringArg(arg, strategy.toString()));
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }
    try {
      time = new Date().getTime();
      proj = new Project(filename, logfile);
      log = proj.getLog();

      System.setProperty("java.awt.headless", headless);
      MedianLRRWorker medianLRRWorker = new MedianLRRWorker(proj,
                                                            readToArray(proj.PROJECT_DIRECTORY.getValue()
                                                                        + regionFileName, log),
                                                            transformationType, scope, outputBase,
                                                            null, false, false, correctXY, false,
                                                            strategy, log);
      medianLRRWorker.execute();
      while (!medianLRRWorker.isDone()) {
        Thread.sleep(100);
      }
      log.report("Finished in " + ext.getTimeElapsed(time));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
