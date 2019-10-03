package org.genvisis.cnv.filesys;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import org.genvisis.cnv.qc.GcAdjustorParameter;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.DoubleVector;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.WorkerHive;
import org.pankratzlab.common.ext;

import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Doubles;

public class Sample implements Serializable {

  public static final long serialVersionUID = 1L;
  public static final String[][] DATA_FIELDS = {{"GC Score", "GCscore", "confidence", "Confidence"},
                                                {"X Raw"}, {"Y Raw"},
                                                {"X", "Xvalue", "Log Ratio", "intensity_1",
                                                 "Signal A"},
                                                {"Y", "Yvalue", "Strength", "intensity_2",
                                                 "Signal B"},
                                                {"Theta"}, {"R"}, {"B Allele Freq"},
                                                {"Log R Ratio"}};
  public static final String[][] GENOTYPE_FIELDS = {{"Allele1 - Forward", "Allele1", "genotype1",
                                                     "Allele1 - Top", "Forward Strand Base Calls",
                                                     "Forced Call", "Forced Call Codes"},
                                                    {"Allele2 - Forward",
                                                     "Forward Strand Base Calls", "genotype2",
                                                     "Allele2 - Top", "Allele B", "Forced Call",
                                                     "Forced Call Codes"},
                                                    {"Allele1 - AB", "Call Codes", "Call"},
                                                    {"Allele2 - AB", "Call Codes", "Call"}}; // ,
                                                                                             // {"Forward
                                                                                             // Strand
                                                                                             // Base
                                                                                             // Calls"},
                                                                                             // {"Call
                                                                                             // Codes"}
  public static final String[] ALL_STANDARD_GENOTYPE_FIELDS = {"Allele1 - AB", "Allele2 - AB",
                                                               "Allele1 - Forward",
                                                               "Allele2 - Forward", "Allele1 - Top",
                                                               "Allele2 - Top", "Allele1 - Design",
                                                               "Allele2 - Design"};
  public static final String[] ALLELE_PAIRS = {"--", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                                               "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "DD",
                                               "DI", "II", "ID"};
  public static final String[] ALT_NULL = {"-", "0", ""};
  public static final String[] ALT_NULLS = {"--", "00", "---", "NoCall", "NC", "NN", "NA", ""};
  public static final String[] AB_PAIRS = {"AA", "AB", "BB"};
  public static final String SAMPLE_FILE_EXTENSION = ".sampRAF";
  // public static final byte PARAMETER_SECTION_BYTES = 13;
  public static final byte PARAMETER_SECTION_BYTES = 17;
  public static final byte PARAMETER_SECTION_NUMMARKERS_LOCATION = 0;
  public static final byte PARAMETER_SECTION_NUMMARKERS_LENGTH = 4;
  public static final byte PARAMETER_SECTION_NULLSTAT_LOCATION = 4;
  public static final byte PARAMETER_SECTION_NULLSTAT_LENGTH = 1;
  public static final byte PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION = 5;
  public static final byte PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LENGTH = 4;
  public static final byte PARAMETER_SECTION_FINGPRNT_LOCATION = 9;
  public static final byte PARAMETER_SECTION_FINGPRNT_LENGTH = 8;
  public static final int MAX_ROWS_PER_WRITE_OPERATION = 500;
  public static final byte NULLSTATUS_GC_LOCATION = 0;
  public static final byte NULLSTATUS_X_LOCATION = 1;
  public static final byte NULLSTATUS_Y_LOCATION = 2;
  public static final byte NULLSTATUS_BAF_LOCATION = 3;
  public static final byte NULLSTATUS_LRR_LOCATION = 4;
  public static final byte NULLSTATUS_ABGENOTYPE_LOCATION = 5;
  public static final byte NULLSTATUS_FOWARDGENOTYPE_LOCATION = 6;
  public static final byte NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION = 7;

  private byte[] abGenotypes;
  private final byte[] forwardGenotypes;
  private final float[] xs;
  private final float[] ys;
  private float[] thetas;
  private float[] rs;
  private final float[] lrrs;
  private final float[] bafs;
  private final float[] gcs;
  private final long fingerprint;
  private final String sampleName;
  private byte nullStatus;
  private final boolean canXYBeNegative;

  public Sample(String sampleName, long fingerprint, float[] gcs, float[] xs, float[] ys,
                float[] bafs, float[] lrrs, byte[] forwardGenotypes, byte[] abGenotypes,
                boolean canXYBeNegative) {
    this.sampleName = sampleName;
    this.fingerprint = fingerprint;
    this.gcs = gcs;
    this.xs = xs;
    this.ys = ys;
    this.bafs = bafs;
    this.lrrs = lrrs;
    this.forwardGenotypes = forwardGenotypes;
    this.abGenotypes = abGenotypes;
    this.canXYBeNegative = canXYBeNegative;
    updateNullStatus();
    // note for future expansion: 8th bit should indicate a second nullStatus byte will be used, to
    // ensure backward compatibility
  }

  /**
   * Data is expected to be:<br />
   * - data[0] => GCs<br />
   * - data[1] -> UNUSED<br />
   * - data[2] -> UNUSED<br />
   * - data[3] => Xs<br />
   * - data[4] => Ys<br />
   * - data[5] -> UNUSED<br />
   * - data[6] -> UNUSED<br />
   * - data[7] => BAFs<br />
   * - data[8] => LRRs<br />
   * 
   * @param sampleName
   * @param fingerprint
   * @param data
   * @param genotypes
   * @param canXYBeNegative
   */
  public Sample(String sampleName, long fingerprint, float[][] data, byte[][] genotypes,
                boolean canXYBeNegative) {
    this.sampleName = sampleName;
    this.fingerprint = fingerprint;
    gcs = data[0];
    xs = data[3];
    ys = data[4];
    thetas = null;
    rs = null;
    bafs = data[7];
    lrrs = data[8];
    forwardGenotypes = genotypes[0];
    abGenotypes = genotypes[1];
    this.canXYBeNegative = canXYBeNegative;
    updateNullStatus();
  }

  public void updateNullStatus() {
    nullStatus = 0;
    if (gcs == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_GC_LOCATION));
    }
    if (xs == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_X_LOCATION));
    }
    if (ys == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_Y_LOCATION));
    }
    if (bafs == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_BAF_LOCATION));
    }
    if (lrrs == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_LRR_LOCATION));
    }
    if (abGenotypes == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_ABGENOTYPE_LOCATION));
    }
    if (forwardGenotypes == null) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_FOWARDGENOTYPE_LOCATION));
    }
    if (canXYBeNegative) {
      nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION));
    }
  }

  public static byte computeNullStatus(float[] gcs, float[] xs, float[] ys, float[] bafs,
                                       float[] lrrs, byte[] abGenotypes, byte[] forwardGenotypes,
                                       boolean canXYBeNegative) {
    byte nullStatus;

    nullStatus = 0;
    if (gcs == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_GC_LOCATION));
    }
    if (xs == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_X_LOCATION));
    }
    if (ys == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_Y_LOCATION));
    }
    if (bafs == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_BAF_LOCATION));
    }
    if (lrrs == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_LRR_LOCATION));
    }
    if (abGenotypes == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_ABGENOTYPE_LOCATION));
    }
    if (forwardGenotypes == null) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_FOWARDGENOTYPE_LOCATION));
    }
    if (canXYBeNegative) {
      nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION));
    }

    return nullStatus;
  }

  public byte getNBytesPerSampleMarker() {
    byte nBytesPerSampleMarker;

    nBytesPerSampleMarker = (byte) Compression.BYTES_PER_SAMPLE_MARKER;
    if (gcs == null) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
    }
    if (xs == null) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
    }
    if (ys == null) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
    }
    if (bafs == null) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
    }
    if (lrrs == null) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
    }
    if (abGenotypes == null && forwardGenotypes == null) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
    }

    return nBytesPerSampleMarker;
  }

  public static byte getNBytesPerSampleMarker(byte nullStatus) {
    byte nBytesPerSampleMarker;

    nBytesPerSampleMarker = (byte) Compression.BYTES_PER_SAMPLE_MARKER;
    if (isGcNull(nullStatus)) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
    }
    if (isXOrYNull(nullStatus)) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
    }
    if (((nullStatus >> NULLSTATUS_Y_LOCATION) & 0x01) == 1) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
    }
    if (isBafOrLrrNull(nullStatus)) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
    }
    if (((nullStatus >> NULLSTATUS_LRR_LOCATION) & 0x01) == 1) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
    }
    if (isAbAndForwardGenotypeNull(nullStatus)) {
      nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
    }

    return nBytesPerSampleMarker;
  }

  public static boolean isNegativeXOrYAllowed(byte nullStatus) {
    return ((nullStatus >> NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION) & 0x01) == 1;
  }

  public static boolean isGcNull(byte nullStatus) {
    return ((nullStatus >> NULLSTATUS_GC_LOCATION) & 0x01) == 1;
  }

  public static boolean isXOrYNull(byte nullStatus) {
    return (((nullStatus >> NULLSTATUS_X_LOCATION) & 0x01) == 1
            || ((nullStatus >> NULLSTATUS_Y_LOCATION) & 0x01) == 1);
  }

  public static boolean isXNull(byte nullStatus) {
    return ((nullStatus >> NULLSTATUS_X_LOCATION) & 0x01) == 1;
  }

  public static boolean isYNull(byte nullStatus) {
    return ((nullStatus >> NULLSTATUS_Y_LOCATION) & 0x01) == 1;
  }

  public static boolean isBafOrLrrNull(byte nullStatus) {
    return (((nullStatus >> NULLSTATUS_BAF_LOCATION) & 0x01) == 1
            || ((nullStatus >> NULLSTATUS_LRR_LOCATION) & 0x01) == 1);
  }

  public static boolean isBafNull(byte nullStatus) {
    return ((nullStatus >> NULLSTATUS_BAF_LOCATION) & 0x01) == 1;
  }

  public static boolean isLrrNull(byte nullStatus) {
    return ((nullStatus >> NULLSTATUS_LRR_LOCATION) & 0x01) == 1;
  }

  public static boolean isAbAndForwardGenotypeNull(byte nullStatus) {
    return (((nullStatus >> NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) == 1
            && ((nullStatus >> NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) == 1);
  }

  public float[][] getAllData() {
    return new float[][] {gcs, xs, ys, thetas, rs, bafs, lrrs};
  }

  public byte[][] getAllGenotypes() {
    return new byte[][] {forwardGenotypes, abGenotypes};
  }

  public String getSampleName() {
    return sampleName;
  }

  public long getFingerprint() {
    return fingerprint;
  }

  public float[] getGCs() {
    return gcs;
  }

  public float[] getXs() {
    return xs;
  }

  public float[] getYs() {
    return ys;
  }

  int getDataLength() {
    int len = 0;
    for (float[] f : getAllData()) {
      if (f != null) {
        len = f.length;
      }
    }
    return len;
  }

  public float[] getThetas() {
    if (thetas == null && xs != null && ys != null) {
      thetas = new float[getDataLength()];
      for (int i = 0; i < thetas.length; i++) {
        thetas[i] = Centroids.calcTheta(xs[i], ys[i]);
      }
    }
    return thetas;
  }

  public float[] getRs() {
    if (rs == null && xs != null && ys != null) {
      rs = new float[getDataLength()];
      for (int i = 0; i < rs.length; i++) {
        rs[i] = Centroids.calcR(xs[i], ys[i]);
      }
    }
    return rs;
  }

  /**
   * @param markerSet
   * @return Map from {@link MarkerDetailSet.Marker} to BAF for this {@link Sample}
   */
  public Map<Marker, Double> markerBAFMap(MarkerDetailSet markerSet) {
    return floatsToDoubleMap(getBAFs(), markerSet);
  }

  /**
   * @param markerSet
   * @param centroids
   * @return Map from {@link MarkerDetailSet.Marker} to BAF for this {@link Sample}
   */
  public Map<Marker, Double> markerBAFMap(MarkerDetailSet markerSet, float[][][] centroids) {
    return floatsToDoubleMap(getBAFs(centroids), markerSet);
  }

  public float[] getBAFs() {
    return bafs;
  }

  public float[] getBAFs(float[][][] centroids) {
    float[] thetas, bafs;

    thetas = getThetas();
    if (thetas == null) {
      return null; // TODO should return null or orig BAFs? TODO log
    }
    bafs = new float[getDataLength()];
    for (int i = 0; i < bafs.length; i++) {
      if (centroids[i] == null) {
        bafs[i] = 1.2f;
        continue;
      }
      bafs[i] = Centroids.calcBAF(thetas[i], centroids[i]);
    }

    return bafs;
  }

  /**
   * @param markerSet
   * @return Map from {@link MarkerDetailSet.Marker} to LRR for this {@link Sample}
   */
  public Map<Marker, Double> markerLRRMap(MarkerDetailSet markerSet) {
    return floatsToDoubleMap(getLRRs(), markerSet);
  }

  public float[] getLRRs() {
    return lrrs;
  }

  /**
   * @param markerSet
   * @param centroids
   * @return Map from {@link MarkerDetailSet.Marker} to LRR for this {@link Sample}
   */
  public Map<Marker, Double> markerLRRMap(MarkerDetailSet markerSet, float[][][] centroids) {
    return floatsToDoubleMap(getLRRs(centroids), markerSet);
  }

  public float[] getLRRs(float[][][] centroids) {
    float[] thetas, rs, lrrs;

    thetas = getThetas();
    rs = getRs();
    if (thetas == null || rs == null) {
      return null; // TODO should return null or orig LRRs? TODO log
    }
    lrrs = new float[xs.length];
    for (int i = 0; i < xs.length; i++) {
      if (centroids[i] == null) {
        lrrs[i] = 1.2f;
        continue;
      }
      lrrs[i] = Centroids.calcLRR(thetas[i], rs[i], centroids[i]);
    }

    return lrrs;
  }

  public static Map<Marker, Double> floatsToDoubleMap(float[] lrrs, MarkerDetailSet markerSet) {
    List<Marker> projOrderMarkers = markerSet.markersAsList();
    ImmutableMap.Builder<Marker, Double> markerLRRMapBuilder = ImmutableMap.builder();
    for (int i = 0; i < lrrs.length; i++) {
      markerLRRMapBuilder.put(projOrderMarkers.get(i), (double) lrrs[i]);
    }
    return markerLRRMapBuilder.build();
  }

  /**
   * Use the centroids stored in the {@link GcAdjustorParameters} to recompute baf
   */
  public float[] getBAF(GcAdjustorParameters params, int sampleIndex, Logger log) {
    GcAdjustorParameter current = params.getGcAdjustorParameters()[sampleIndex];
    verifyParams(params, current);
    if (bafs == null) {
      return null;
    } else {
      float[] recompBaf = params.getCentroids() == null ? bafs.clone()
                                                        : getBAFs(params.getCentroids()
                                                                        .getCentroids());
      return recompBaf;
    }
  }

  /**
   * Use the centroids stored in the {@link GcAdjustorParameters} to recompute lrr, and then gc
   * correct the result
   */
  public float[] getGCCorrectedLRR(GcAdjustorParameters params, int sampleIndex, Logger log) {
    GcAdjustorParameter current = params.getGcAdjustorParameters()[sampleIndex];
    verifyParams(params, current);

    if (lrrs == null) {
      return null;
    } else {
      float[] recompLrrs = params.getCentroids() == null ? lrrs.clone()
                                                         : getLRRs(params.getCentroids()
                                                                         .getCentroids());
      for (int i = 0; i < recompLrrs.length; i++) {
        recompLrrs[i] = (float) current.adjust(GC_CORRECTION_METHOD.GENVISIS_GC, recompLrrs[i],
                                               params.getGcContent()[i]);
      }
      return recompLrrs;
    }
  }

  private void verifyParams(GcAdjustorParameters params, GcAdjustorParameter current) {
    if (!current.getSample().equals(sampleName)) {
      throw new IllegalArgumentException("Mismatched sample, was given " + current.getSample()
                                         + " and should have " + sampleName);
    }
    if (params.getMarkerFingerprint() != fingerprint) {
      throw new IllegalArgumentException("Mismatched sample, was given "
                                         + params.getMarkerFingerprint() + " and should have "
                                         + fingerprint);
    }
  }

  public byte[] getForwardGenotypes() {
    return forwardGenotypes;
  }

  public byte[] getForwardGenotypes(float gcThreshold) {
    byte[] result = new byte[forwardGenotypes.length];

    for (int i = 0; i < result.length; i++) {
      if (gcs[i] <= gcThreshold) {
        result[i] = (byte) 0;
      } else {
        result[i] = forwardGenotypes[i];
      }
    }

    return result;
  }

  public boolean getCanXYBeNegative() {
    return canXYBeNegative;
  }

  public byte[] getAB_GenotypesAfterFilters(String[] markerNames,
                                            ClusterFilterCollection clusterFilterCollection,
                                            float gcThreshold) {
    byte[] result;
    float realX;
    float realY;
    ClusterFilter clusterFilter;
    ArrayList<ClusterFilter> clusterFilterArray;

    getThetas();
    getRs();

    if (markerNames.length != abGenotypes.length) {
      System.err.println("Error - need to pass the full set of markerNames to getAB_GenotypesAfterFilters");
      return null;
    }

    result = new byte[abGenotypes.length];
    for (int i = 0; i < abGenotypes.length; i++) {
      if (gcs[i] < gcThreshold) {
        result[i] = (byte) -1;
      } else {
        result[i] = abGenotypes[i];
      }

      clusterFilterArray = clusterFilterCollection.getClusterFilters(markerNames[i]);
      if (clusterFilterArray != null) {
        for (int j = 0; j < clusterFilterArray.size(); j++) {
          clusterFilter = clusterFilterArray.get(j);
          switch (clusterFilter.getPlotType()) {
            case 0:
              // realX = -9;
              // realY = -9;
              // System.err.println("Error - PlotType cannot be 0 for ClusterFilter #"+(j+1)+" for
              // marker '"+markerNames[i]+"' as we've done away with raw Xs and Ys");
              // break;
              // case 1:
              realX = xs[i];
              realY = ys[i];
              break;
            case 1:
              realX = getThetas()[i];
              realY = getRs()[i];
              break;
            case 2:
              realX = bafs[i];
              realY = lrrs[i];
              break;
            default:
              System.err.println("Error - invalid PlotType in ClusterFilter #" + (j + 1)
                                 + " for marker '" + i + "'");
              realX = -9;
              realY = -9;
          }
          if (realX >= clusterFilter.getXMin() && realY >= clusterFilter.getYMin()
              && realX <= clusterFilter.getXMax() && realY <= clusterFilter.getYMax()) {
            result[i] = clusterFilter.getCluterGenotype();
          }
        }
      }
    }
    return result;
  }

  public byte[] getAB_GenotypesAfterFilters(String[] targetMarkers,
                                            int[] indicesOfTargetMarkersInProj,
                                            ClusterFilterCollection clusterFilterCollection,
                                            float gcThreshold) {
    byte[] result;
    float realX;
    float realY;
    ClusterFilter clusterFilter;
    ArrayList<ClusterFilter> clusterFilterArray;

    if (indicesOfTargetMarkersInProj == null) {
      return getAB_GenotypesAfterFilters(targetMarkers, clusterFilterCollection, gcThreshold);
    } else {
      result = new byte[targetMarkers.length];
      for (int i = 0; i < targetMarkers.length; i++) {
        if (gcs[indicesOfTargetMarkersInProj[i]] < gcThreshold) {
          result[i] = (byte) -1;
        } else {
          result[i] = abGenotypes[indicesOfTargetMarkersInProj[i]];
        }

        clusterFilterArray = clusterFilterCollection.getClusterFilters(targetMarkers[i]);
        if (clusterFilterArray != null) {
          for (int j = 0; j < clusterFilterArray.size(); j++) {
            clusterFilter = clusterFilterArray.get(j);
            switch (clusterFilter.getPlotType()) {
              case 0:
                // realX = -9;
                // realY = -9;
                // System.err.println("Error - PlotType cannot be 0 for ClusterFilter #" + (j + 1) +
                // " for marker '" + targetMarkers[i] + "' as we've done away with raw Xs and Ys");
                // break;
                // case 1:
                realX = xs[indicesOfTargetMarkersInProj[i]];
                realY = ys[indicesOfTargetMarkersInProj[i]];
                break;
              case 1:
                realX = getThetas()[indicesOfTargetMarkersInProj[i]];
                realY = getRs()[indicesOfTargetMarkersInProj[i]];
                break;
              case 2:
                realX = bafs[indicesOfTargetMarkersInProj[i]];
                realY = lrrs[indicesOfTargetMarkersInProj[i]];
                break;
              default:
                System.err.println("Error - invalid PlotType in ClusterFilter #" + (j + 1)
                                   + " for marker '" + targetMarkers[i] + "'");
                realX = -9;
                realY = -9;
            }
            if (realX >= clusterFilter.getXMin() && realY >= clusterFilter.getYMin()
                && realX <= clusterFilter.getXMax() && realY <= clusterFilter.getYMax()) {
              result[i] = clusterFilter.getCluterGenotype();
            }
          }
        }
      }
      return result;
    }
  }

  public byte[] getAB_Genotypes() {
    return abGenotypes;
  }

  public byte[] getAB_Genotypes(int[] indicesOfTargetMarkersInProj) {
    byte[] result;

    if (indicesOfTargetMarkersInProj == null) {
      return abGenotypes;

    } else {
      result = new byte[indicesOfTargetMarkersInProj.length];
      for (int i = 0; i < result.length; i++) {
        result[i] = abGenotypes[indicesOfTargetMarkersInProj[i]];
      }
      return result;
    }
  }

  public void setAB_Genotypes(byte[] abGenotypes) {
    this.abGenotypes = abGenotypes;
    updateNullStatus();
  }

  private String getMarkerString(int i, String markerName) {
    StringBuilder sb;
    sb = new StringBuilder(markerName == null ? "marker_" + i : markerName).append("\t");
    sb.append(gcs != null && gcs.length > i ? gcs[i] : ".").append("\t");
    sb.append(xs != null && xs.length > i ? xs[i] : ".").append("\t");
    sb.append(ys != null && ys.length > i ? ys[i] : ".").append("\t");
    sb.append(thetas != null && thetas.length > i ? thetas[i] : ".").append("\t");
    sb.append(rs != null && rs.length > i ? rs[i] : ".").append("\t");
    sb.append(lrrs != null && lrrs.length > i ? lrrs[i] : ".").append("\t");
    sb.append(bafs != null && bafs.length > i ? bafs[i] : ".").append("\t");
    sb.append(forwardGenotypes != null
              && forwardGenotypes.length > i ? ALLELE_PAIRS[forwardGenotypes[i]] : ".")
      .append("\t");
    sb.append(abGenotypes != null
              && abGenotypes.length > i ? (abGenotypes[i] == -1 ? "--" : AB_PAIRS[abGenotypes[i]])
                                        : ".");

    return sb.toString();
  }

  public void writeToOut(String[] markerNames) {
    if (markerNames != null && markerNames.length != lrrs.length) {
      System.err.println("Error - MarkerNames (n=" + markerNames.length
                         + ") do not match up with the number of LRRs/BAFs/etc (n=" + lrrs.length
                         + ")");
      System.exit(1);
    }
    System.out.println("SNP\tGC Score\tX\tY\tTheta\tR\tLRR\tBAF\tGenotypes\tAB_Genotypes");
    for (int i = 0; i < getDataLength(); i++) {
      System.out.println(getMarkerString(i, markerNames == null ? null : markerNames[i]));
    }

  }

  public void writeToFile(String[] markerNames, String filename) {
    PrintWriter writer;

    if (markerNames != null && markerNames.length != lrrs.length) {
      System.err.println("Error - MarkerNames (n=" + markerNames.length
                         + ") do not match up with the number of LRRs/BAFs/etc (n=" + lrrs.length
                         + ")");
      System.exit(1);
    }
    try {
      writer = Files.openAppropriateWriter(filename);
      writer.println("SNP\tGC Score\tX\tY\tTheta\tR\tLRR\tBAF\tGenotypes\tAB_Genotypes");
      for (int i = 0; i < getDataLength(); i++) {
        writer.println(getMarkerString(i, markerNames == null ? null : markerNames[i]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + filename);
      e.printStackTrace();
    }
  }

  public float[][][] loadCentroids(Project proj) {
    Centroids centroids;

    // centroids = Centroids.load(proj.getFilename(proj.CUSTOM_CENTROIDS_FILENAME),
    // proj.getJarStatus());
    centroids = Centroids.load(proj.CUSTOM_CENTROIDS_FILENAME.getValue());
    if (centroids.getFingerprint() != getFingerprint()) {
      System.err.println("Error - mismatched fingerprint for " + sampleName);
    }
    return centroids.getCentroids();
  }

  public void compareCalculationsFile(Project proj, String[] markerNames, String filename) {
    PrintWriter writer;
    float[][][] centroids;
    float[] compBAFs, compLRRs;

    centroids = loadCentroids(proj);
    compBAFs = getBAFs(centroids);
    compLRRs = getLRRs(centroids);

    if (markerNames.length != lrrs.length) {
      System.err.println("Error - MarkerNames (n=" + markerNames.length
                         + ") do not match up with the number of LRRs/BAFs/etc (n=" + lrrs.length
                         + ")");
      System.exit(1);
    }
    try {
      writer = Files.openAppropriateWriter(filename);
      writer.println("SNP\tX\tY\tTheta\tR\tcompTheta\tcompR\tLRR\tBAF\tcompLRR\tcompBAF");
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(markerNames[i] + "\t" + xs[i] + "\t" + ys[i] + "\t" + thetas[i] + "\t"
                       + rs[i] + "\t" + bafs[i] + "\t" + lrrs[i] + "\t" + compBAFs[i] + "\t"
                       + compLRRs[i]);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + filename);
      e.printStackTrace();
    }
  }

  public void compLRRs(Project proj) {
    float[] compLRRs, diffs;

    compLRRs = getLRRs(loadCentroids(proj));
    diffs = new float[lrrs.length];
    for (int i = 0; i < lrrs.length; i++) {
      diffs[i] = lrrs[i] - compLRRs[i];
    }
    SerializedFiles.writeSerial(diffs, proj.PROJECT_DIRECTORY.getValue() + "comps/" + sampleName
                                       + ".comp");
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  /**
   * Save the instance of Sample to hard drive in Random Access File format.
   *
   * @param filename
   */
  public void saveToRandomAccessFile(String filename) {
    saveToRandomAccessFile(filename, null, null);
  }

  public void saveToRandomAccessFile(String filename, Hashtable<String, Float> allOutliers,
                                     String sampleName) {
    File fileTmp;
    RandomAccessFile rafFile;
    Hashtable<String, Float> outOfRangeValuesEachSample;
    byte[] outOfRangeValuesWriteBuffer;
    int bytesRemained;
    byte[] writeBuffer = null;
    byte[] parameters;
    byte[] temp;
    int writeBufferIndex;
    byte bytesPerSampleMarker;
    boolean isOutlier;

    fileTmp = new File(filename);
    if (new File(ext.parseDirectoryOfFile(filename)).getFreeSpace() <= ((long) getDataLength()
                                                                        * Compression.BYTES_PER_SAMPLE_MARKER)) {
      System.err.println("Not enough space (available: "
                         + ext.prettyUpSize(new File(ext.parseDirectoryOfFile(filename)).getFreeSpace(),
                                            1)
                         + ") for all the new data to be created (required: "
                         + ext.prettyUpSize(((long) getDataLength()
                                             * Compression.BYTES_PER_SAMPLE_MARKER),
                                            1)
                         + ").");
      return;
    }
    if (fileTmp.exists()) {
      fileTmp.delete();
    }

    bytesPerSampleMarker = getNBytesPerSampleMarker();
    bytesRemained = getDataLength() * bytesPerSampleMarker;
    long time = new Date().getTime();
    outOfRangeValuesEachSample = new Hashtable<String, Float>();
    try {
      rafFile = new RandomAccessFile(filename, "rw");

      // rafFile.writeInt(getDataLength());
      // rafFile.writeByte(nullStatus);
      // rafFile.writeLong(fingerprint);
      parameters = new byte[PARAMETER_SECTION_BYTES];
      temp = Compression.intToBytes(getDataLength());
      for (int i = 0; i < temp.length; i++) {
        parameters[PARAMETER_SECTION_NUMMARKERS_LOCATION + i] = temp[i];
      }

      parameters[PARAMETER_SECTION_NULLSTAT_LOCATION] = nullStatus;

      temp = Compression.longToBytes(fingerprint);
      for (int i = 0; i < temp.length; i++) {
        parameters[PARAMETER_SECTION_FINGPRNT_LOCATION + i] = temp[i];
      }
      rafFile.write(parameters);

      writeBufferIndex = 0;
      for (int j = 0; j < getDataLength(); j++) {
        if (writeBufferIndex == 0) {
          writeBuffer = new byte[Math.min(Integer.MAX_VALUE, bytesRemained)];
        }
        if (gcs != null) {
          try {
            Compression.gcBafCompress(gcs[j], writeBuffer, writeBufferIndex);
          } catch (Elision e) {
            System.err.println("Error - problem saving sampleRAF file '" + filename
                               + "' for sample " + sampleName + " since the marker in index " + j
                               + " has a GC value of " + gcs[j]
                               + "; in the past this has happened if the final report file was truncated and did not contain all of the markers");
          }
          writeBufferIndex += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
        }
        if (xs != null) {
          if (canXYBeNegative) {
            isOutlier = !Compression.xyCompressAllowNegative(xs[j], writeBuffer, writeBufferIndex);
          } else {
            isOutlier = !Compression.xyCompressPositiveOnly(xs[j], writeBuffer, writeBufferIndex);
          }
          if (isOutlier) {
            outOfRangeValuesEachSample.put(j + "\tx", xs[j]);
            if (allOutliers != null) {
              allOutliers.put(j + "\t" + sampleName + "\tx", xs[j]);
            }
          }
          writeBufferIndex += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
        }
        if (ys != null) {
          if (canXYBeNegative) {
            isOutlier = !Compression.xyCompressAllowNegative(ys[j], writeBuffer, writeBufferIndex);
          } else {
            isOutlier = !Compression.xyCompressPositiveOnly(ys[j], writeBuffer, writeBufferIndex);
          }
          if (isOutlier) {
            outOfRangeValuesEachSample.put(j + "\ty", ys[j]);
            if (allOutliers != null) {
              allOutliers.put(j + "\t" + sampleName + "\ty", ys[j]);
            }
          }
          writeBufferIndex += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
        }
        if (bafs != null) {
          Compression.gcBafCompress(bafs[j], writeBuffer, writeBufferIndex);
          writeBufferIndex += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
        }
        if (lrrs != null) {
          if (Compression.lrrCompress(lrrs[j], writeBuffer, writeBufferIndex) == -1) {
            outOfRangeValuesEachSample.put(j + "\tlrr", lrrs[j]);
            if (allOutliers != null) {
              allOutliers.put(j + "\t" + sampleName + "\tlrr", lrrs[j]);
            }
          }
          writeBufferIndex += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
        }
        if (abGenotypes != null || forwardGenotypes != null) {
          writeBuffer[writeBufferIndex] = Compression.genotypeCompress(abGenotypes == null ? -1
                                                                                           : abGenotypes[j],
                                                                       forwardGenotypes == null ? 0
                                                                                                : forwardGenotypes[j]);
          writeBufferIndex += Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
        }
        bytesRemained -= bytesPerSampleMarker;
        if (writeBufferIndex >= writeBuffer.length) {
          rafFile.write(writeBuffer);
          writeBufferIndex = 0;
        }
      }

      if (outOfRangeValuesEachSample != null && outOfRangeValuesEachSample.size() > 0) {
        // TODO log here?
        outOfRangeValuesWriteBuffer = Compression.objToBytes(outOfRangeValuesEachSample);
        rafFile.write(outOfRangeValuesWriteBuffer);
        rafFile.seek(PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION);
        rafFile.writeInt(outOfRangeValuesWriteBuffer.length);
      }

      rafFile.close();
    } catch (IOException e) {
      e.printStackTrace();
    } catch (Elision e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    time = new Date().getTime() - time;
    // System.out.println("Random Access File. ---- Finished writing to all the samples in " +
    // (time/60000) + " min " + ((time%60000)/1000) + " sec");
  }

  public static Sample loadFromSerialized(String filename) {
    return (Sample) SerializedFiles.readSerial(filename, true);
  }

  public static Sample loadFromRandomAccessFile(String filename) {
    return loadFromRandomAccessFile(filename, true, true, true, true, true);
  }

  /**
   * Load data from Random Access Files organized by samples.
   */
  @SuppressWarnings("unchecked")
  public static Sample loadFromRandomAccessFile(String filename, boolean loadGC, boolean loadXY,
                                                boolean loadBAF, boolean loadLRR,
                                                boolean loadAbOrForwardGenotypes) {
    Sample result = null;
    int numMarkers;
    RandomAccessFile file;
    String sampleName;
    long fingerPrint;
    float[] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
    byte[] abGenotypes = null, fwdGenotypes = null, readBuffer;
    byte[] genoTypeTmp;
    byte[] temp;
    int index, indexStart;
    int numBytesOfOutOfRangeValues;
    Hashtable<String, Float> outOfRangeValues = null;
    int outlierSectionLocation;
    byte nullStatus;
    byte bytesPerSampleMarker;

    sampleName = ext.rootOf(filename);
    try {
      file = new RandomAccessFile(filename, "r");
      readBuffer = new byte[(int) file.length()]; // numMarkers * BYTES_PER_SAMPLE_MARKER
      file.read(readBuffer);
      file.close();

      temp = new byte[PARAMETER_SECTION_NUMMARKERS_LENGTH];
      for (int i = 0; i < temp.length; i++) {
        temp[i] = readBuffer[PARAMETER_SECTION_NUMMARKERS_LOCATION + i];
      }
      numMarkers = Compression.bytesToInt(temp);

      nullStatus = readBuffer[PARAMETER_SECTION_NULLSTAT_LOCATION];

      temp = new byte[PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LENGTH];
      for (int i = 0; i < temp.length; i++) {
        temp[i] = readBuffer[PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION + i];
      }
      numBytesOfOutOfRangeValues = Compression.bytesToInt(temp);

      temp = new byte[PARAMETER_SECTION_FINGPRNT_LENGTH];
      for (int i = 0; i < temp.length; i++) {
        temp[i] = readBuffer[PARAMETER_SECTION_FINGPRNT_LOCATION + i];
      }
      fingerPrint = Compression.bytesToLong(temp);

      bytesPerSampleMarker = getNBytesPerSampleMarker(nullStatus);

      if (numBytesOfOutOfRangeValues > 0) {
        outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampleMarker;
        outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer,
                                                                             outlierSectionLocation,
                                                                             numBytesOfOutOfRangeValues);
      } else {
        outOfRangeValues = new Hashtable<>();
      }

      indexStart = PARAMETER_SECTION_BYTES;
      index = indexStart;
      if (!isGcNull(nullStatus)) {
        if (loadGC) {
          gcs = new float[numMarkers];
          for (int j = 0; j < numMarkers; j++) {
            gcs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[index],
                                                             readBuffer[index + 1]});
            index += bytesPerSampleMarker;
          }
        }
        indexStart += 2;
      }
      index = indexStart;
      if (!isXNull(nullStatus)) {
        if (loadXY) {
          xs = new float[numMarkers];
          for (int j = 0; j < numMarkers; j++) {
            xs[j] = Compression.xyDecompressPositiveOnly(new byte[] {readBuffer[index],
                                                                     readBuffer[index + 1]});
            if (xs[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
              Float oor = outOfRangeValues.get(j + "\tx");
              if (oor == null) {
                System.err.println("Error - X value for sample " + sampleName + " and marker " + j
                                   + " was out of range, but couldn't be found.");
              } else {
                xs[j] = outOfRangeValues.get(j + "\tx");
              }
            }
            index += bytesPerSampleMarker;
          }
        }
        indexStart += 2;
      }
      index = indexStart;
      if (!isYNull(nullStatus)) {
        if (loadXY) {
          ys = new float[numMarkers];
          for (int j = 0; j < numMarkers; j++) {
            ys[j] = Compression.xyDecompressPositiveOnly(new byte[] {readBuffer[index],
                                                                     readBuffer[index + 1]});
            if (ys[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
              Float oor = outOfRangeValues.get(j + "\ty");
              if (oor == null) {
                System.err.println("Error - Y value for sample " + sampleName + " and marker " + j
                                   + " was out of range, but couldn't be found.");
              } else {
                ys[j] = oor.floatValue();
              }
            }
            index += bytesPerSampleMarker;
          }

        }
        indexStart += 2;
      }
      index = indexStart;
      if (!isBafNull(nullStatus)) {
        if (loadBAF) {
          bafs = new float[numMarkers];
          for (int j = 0; j < numMarkers; j++) {
            bafs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[index],
                                                              readBuffer[index + 1]});
            index += bytesPerSampleMarker;
          }
        }
        indexStart += 2;
      }
      index = indexStart;
      if (!isLrrNull(nullStatus)) {
        if (loadLRR) {
          lrrs = new float[numMarkers];
          for (int j = 0; j < numMarkers; j++) {
            lrrs[j] = Compression.lrrDecompress(new byte[] {readBuffer[index],
                                                            readBuffer[index + 1],
                                                            readBuffer[index + 2]});
            if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT) {
              Float oor = outOfRangeValues.get(j + "\tlrr");
              if (oor == null) {
                System.err.println("Error - LRR value for sample " + sampleName + " and marker " + j
                                   + " was out of range, but couldn't be found.");
              } else {
                lrrs[j] = oor.floatValue();
              }
            }
            index += bytesPerSampleMarker;
          }
        }
        indexStart += 3;
      }
      index = indexStart;
      if ((!isAbAndForwardGenotypeNull(nullStatus)) && loadAbOrForwardGenotypes) {
        abGenotypes = new byte[numMarkers];
        fwdGenotypes = new byte[numMarkers];
        for (int j = 0; j < numMarkers; j++) {
          genoTypeTmp = Compression.genotypeDecompress(readBuffer[index]);
          abGenotypes[j] = genoTypeTmp[0];
          fwdGenotypes[j] = genoTypeTmp[1];
          index += bytesPerSampleMarker;
        }
      }
      readBuffer = null;
      temp = null;
      result = new Sample(sampleName, fingerPrint, gcs, xs, ys, bafs, lrrs, fwdGenotypes,
                          abGenotypes, isNegativeXOrYAllowed(nullStatus));
    } catch (IOException e) {
      e.printStackTrace();
    } catch (ClassNotFoundException e) {
      e.printStackTrace();
    }
    return result;
  }

  // @SuppressWarnings("unchecked")
  // public static void loadFromRandomAccessFileWithoutDecompress(RandomAccessFile sampleFile,
  // byte[] readBuffer, int indexOfCurrentSample, int indexOfFirstMarkerToLoad, byte
  // bytesPerSampleMarker, int numMarkersInProj, Hashtable<String, Float> allOutliers) {
  // byte[] outliersBuffer;
  // Hashtable<String, Float> sampleOutlierHash;
  // Enumeration<String> keys;
  // String currentKey;
  // int outlierSectionSize = 0;
  //
  // try {
  // if (allOutliers != null) {
  // sampleFile.seek(Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOC);
  // outlierSectionSize = sampleFile.readInt();
  // }
  //
  // sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + indexOfFirstMarkerToLoad *
  // bytesPerSampleMarker);
  // sampleFile.read(readBuffer);
  //
  // if (outlierSectionSize > 0) {
  // sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker);
  // outliersBuffer = new byte[outlierSectionSize];
  // sampleFile.read(outliersBuffer);
  // sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(outliersBuffer);
  // keys = sampleOutlierHash.keys();
  // while (keys.hasMoreElements()) {
  // currentKey = keys.nextElement();
  // allOutliers.put(indexOfCurrentSample + "\t" + currentKey, sampleOutlierHash.get(currentKey));
  // }
  // }
  // } catch (IOException e) {
  // // TODO Auto-generated catch block
  // e.printStackTrace();
  // } catch (ClassNotFoundException e) {
  // // TODO Auto-generated catch block
  // e.printStackTrace();
  // }
  // }

  @SuppressWarnings("unchecked")
  public static void loadFromRandomAccessFileWithoutDecompress(RandomAccessFile sampleFile,
                                                               byte[] readBuffer,
                                                               boolean seekOrLoadWholeFile,
                                                               int indexOfCurrentSample,
                                                               int indexOfFirstMarkerToLoad,
                                                               byte bytesPerSampleMarker,
                                                               int numMarkersInProj,
                                                               Hashtable<String, Float> allOutliers,
                                                               Logger log) {
    byte[] outliersBuffer;
    Hashtable<String, Float> sampleOutlierHash;
    Enumeration<String> keys;
    String currentKey;
    int outlierSectionSize = 0;
    byte[] readBufferLocal;
    int pointer;
    long seekPointer;

    try {
      if (seekOrLoadWholeFile) {
        if (allOutliers != null) {
          sampleFile.seek(Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION);
          outlierSectionSize = sampleFile.readInt();
        }

        seekPointer = Sample.PARAMETER_SECTION_BYTES
                      + indexOfFirstMarkerToLoad * bytesPerSampleMarker;
        if (seekPointer != sampleFile.getFilePointer()) {
          sampleFile.seek(seekPointer);
        }
        sampleFile.read(readBuffer);

        if (outlierSectionSize > 0) {
          sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker);
          outliersBuffer = new byte[outlierSectionSize];
          sampleFile.read(outliersBuffer);
          sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(outliersBuffer);
          keys = sampleOutlierHash.keys();
          while (keys.hasMoreElements()) {
            currentKey = keys.nextElement();
            allOutliers.put(indexOfCurrentSample + "\t" + currentKey,
                            sampleOutlierHash.get(currentKey));
          }
        }

      } else {
        readBufferLocal = new byte[(int) sampleFile.length()];
        sampleFile.seek(0);
        sampleFile.read(readBufferLocal);
        pointer = Sample.PARAMETER_SECTION_BYTES + indexOfFirstMarkerToLoad * bytesPerSampleMarker;
        for (int i = 0; i < readBuffer.length; i++) {
          readBuffer[i] = readBufferLocal[pointer];
          pointer++;
        }

        if (allOutliers != null) {
          pointer = Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION;
          outlierSectionSize = Compression.bytesToInt(readBufferLocal, pointer);
          if (outlierSectionSize > 0) {
            pointer = Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker;
            sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(readBufferLocal,
                                                                                  pointer,
                                                                                  outlierSectionSize);
            keys = sampleOutlierHash.keys();
            while (keys.hasMoreElements()) {
              currentKey = keys.nextElement();
              allOutliers.put(indexOfCurrentSample + "\t" + currentKey,
                              sampleOutlierHash.get(currentKey));
            }
          }
        }
      }
    } catch (IOException ioe) {
      log.reportError("Error reading from a " + Sample.SAMPLE_FILE_EXTENSION + " file");
      log.reportException(ioe);
    } catch (ClassNotFoundException cnfe) {
      log.reportError("Error reading from a " + Sample.SAMPLE_FILE_EXTENSION + " file");
      log.reportException(cnfe);
    }
  }

  /*
   *
   */
  @SuppressWarnings("unchecked")
  public static void loadSampleFileWithoutDecompress(String sampleFileName, byte[] outputBuffer,
                                                     Hashtable<String, Float> allOutliers,
                                                     int indexOfCurrentSample,
                                                     int indexOfFirstMarkerToLoad,
                                                     byte bytesPerSampleMarker,
                                                     int numMarkersInProj) {
    RandomAccessFile sampleFile;
    byte[] outliersBuffer;
    Hashtable<String, Float> sampleOutlierHash;
    Enumeration<String> keys;
    String currentKey;
    int outlierSectionSize = 0;
    long seekPointer;

    try {
      sampleFile = new RandomAccessFile(sampleFileName, "r");
      if (allOutliers != null) {
        sampleFile.seek(Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION);
        outlierSectionSize = sampleFile.readInt();
      }

      seekPointer = Sample.PARAMETER_SECTION_BYTES
                    + indexOfFirstMarkerToLoad * bytesPerSampleMarker;
      if (seekPointer != sampleFile.getFilePointer()) {
        sampleFile.seek(seekPointer);
      }
      sampleFile.read(outputBuffer);

      if (outlierSectionSize > 0) {
        sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker);
        outliersBuffer = new byte[outlierSectionSize];
        sampleFile.read(outliersBuffer);
        sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(outliersBuffer);
        keys = sampleOutlierHash.keys();
        while (keys.hasMoreElements()) {
          currentKey = keys.nextElement();
          allOutliers.put(indexOfCurrentSample + "\t" + currentKey,
                          sampleOutlierHash.get(currentKey));
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    } catch (ClassNotFoundException e) {
      e.printStackTrace();
    } catch (ArrayIndexOutOfBoundsException e) {
      e.printStackTrace();
    }
  }

  public static byte getNullstatusFromRandomAccessFile(String filename) {
    byte nullStatusOfTheFile = Byte.MIN_VALUE;
    RandomAccessFile sampleFile;

    try {
      sampleFile = new RandomAccessFile(filename, "r");
      sampleFile.readInt();
      nullStatusOfTheFile = sampleFile.readByte();
      sampleFile.close();
    } catch (FileNotFoundException e1) {
      e1.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return nullStatusOfTheFile;
  }

  @SuppressWarnings("unchecked")
  public static Hashtable<String, Float> loadOutOfRangeValuesFromRandomAccessFile(String filename) throws Exception {
    int numMarkers;
    RandomAccessFile file;
    byte[] readBuffer;
    byte[] temp;
    int numBytesOfOutOfRangeValues;
    Hashtable<String, Float> outOfRangeValues = null;
    int outlierSectionLocation;
    byte nullStatus;
    byte bytesPerSampleMarker;

    file = new RandomAccessFile(filename, "r");
    readBuffer = new byte[(int) file.length()]; // numMarkers * BYTES_PER_SAMPLE_MARKER
    file.read(readBuffer);
    file.close();

    temp = new byte[PARAMETER_SECTION_NUMMARKERS_LENGTH];
    for (int i = 0; i < temp.length; i++) {
      temp[i] = readBuffer[PARAMETER_SECTION_NUMMARKERS_LOCATION + i];
    }
    numMarkers = Compression.bytesToInt(temp);

    nullStatus = readBuffer[PARAMETER_SECTION_NULLSTAT_LOCATION];

    temp = new byte[PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LENGTH];
    for (int i = 0; i < temp.length; i++) {
      temp[i] = readBuffer[PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION + i];
    }
    numBytesOfOutOfRangeValues = Compression.bytesToInt(temp);

    bytesPerSampleMarker = getNBytesPerSampleMarker(nullStatus);

    if (numBytesOfOutOfRangeValues > 0) {
      outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampleMarker;
      byte[] vals = new byte[numBytesOfOutOfRangeValues];
      System.arraycopy(readBuffer, outlierSectionLocation, vals, 0, numBytesOfOutOfRangeValues);
      outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(vals);
    }

    return outOfRangeValues;
  }

  @SuppressWarnings("unchecked")
  public static Hashtable<String, Float> loadOutOfRangeValues(RandomAccessFile file,
                                                              long outlierSectionLocation) throws Exception {
    Hashtable<String, Float> outOfRangeValues = null;
    int numBytesOfOutOfRangeValues;
    byte[] readBuffer;

    file.seek(outlierSectionLocation);
    numBytesOfOutOfRangeValues = file.readInt();
    if (numBytesOfOutOfRangeValues > 0) {
      readBuffer = new byte[numBytesOfOutOfRangeValues];
      file.read(readBuffer);
      outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer);
    }

    return outOfRangeValues;
  }

  // TODO if Outlier class is created, then move these two methods there
  public static Hashtable<String, Float> loadOutOfRangeValuesFromSerializable(Project proj) throws Exception {
    return loadOutOfRangeValuesFromSerializable(proj.SAMPLE_DIRECTORY.getValue(true, true)
                                                + "outliers.ser");
  }

  @SuppressWarnings("unchecked")
  public static Hashtable<String, Float> loadOutOfRangeValuesFromSerializable(String filenameOfSerializableOutOfRangeFile) throws Exception {
    return (Hashtable<String, Float>) SerializedFiles.readSerial(filenameOfSerializableOutOfRangeFile);
  }

  public static void dumpOutOfRangeValues(String outOfRangeFileName, String outFileName,
                                          boolean includeKeysOnly) {
    Hashtable<String, Float> outOfRangeValues;

    if (outOfRangeFileName == null || outFileName == null) {
      System.err.println("Error - outOfRangeFileName and outFileName cannot be null.");
    }
    try {
      outOfRangeValues = loadOutOfRangeValuesFromSerializable(outOfRangeFileName);
      dumpOutOfRangeValues(outOfRangeValues, outFileName, includeKeysOnly);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void dumpOutOfRangeValues(Hashtable<String, Float> outOfRangeValues,
                                          String outFileName, boolean includeKeysOnly) {
    String[] keys;
    PrintWriter outfile;

    try {
      keys = HashVec.getKeys(outOfRangeValues);
      outfile = new PrintWriter(new FileOutputStream(outFileName));
      if (!includeKeysOnly) {
        for (String key : keys) {
          outfile.println(key + "\t" + outOfRangeValues.get(key));
        }
      } else {
        for (String key : keys) {
          outfile.println(key);
        }
      }
      outfile.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

  @SuppressWarnings("unchecked")
  public static void verifyAndGenerateOutliers(Project proj, int numthreads, boolean verbose) {
    String outlierFileName = proj.SAMPLE_DIRECTORY.getValue() + "outliers.ser";
    Hashtable<String, Float> outliers = new Hashtable<String, Float>();

    if (Files.exists(outlierFileName)) {
      if (verbose) {
        proj.getLog().reportError(outlierFileName + " already exists, cancelling");
      }
    } else {
      proj.getLog().reportTimeInfo("Did not see " + outlierFileName + ", generating now");
      class HashLoadResult {

        private final Hashtable<String, Float> outlier;
        private final String sampleName;

        public HashLoadResult(Hashtable<String, Float> outlier, String sampleName) {
          super();
          this.outlier = outlier;
          this.sampleName = sampleName;
        }

        public Hashtable<String, Float> getOutlier() {
          return outlier;
        }

        public String getSampleName() {
          return sampleName;
        }

      }
      WorkerHive<HashLoadResult> hive = new WorkerHive<HashLoadResult>(numthreads, 10,
                                                                       proj.getLog());

      for (int i = 0; i < proj.getNumberOfParsedSamples(); i++) {
        final String currentSampleRAF = proj.SAMPLE_DIRECTORY.getValue() + proj.getSamples()[i]
                                        + Sample.SAMPLE_FILE_EXTENSION;
        final String sampleName = proj.getSamples()[i];

        final int index = i;
        final Logger log = proj.getLog();
        Callable<HashLoadResult> outCallable = new Callable<HashLoadResult>() {

          @Override
          public HashLoadResult call() throws Exception {
            // TODO Auto-generated method stub
            if (index % 100 == 0) {
              log.reportTimeInfo("loading sample " + (index + 1) + " on "
                                 + Thread.currentThread().getName());
            }

            Hashtable<String, Float> outs = Sample.loadOutOfRangeValuesFromRandomAccessFile(currentSampleRAF);
            return new HashLoadResult(outs, sampleName);
          }

        };
        hive.addCallable(outCallable);
      }
      hive.execute(true);
      ArrayList<HashLoadResult> hashLoadResults = hive.getResults();
      for (int i = 0; i < hashLoadResults.size(); i++) {
        if (hashLoadResults.get(i).getOutlier() != null
            && hashLoadResults.get(i).getOutlier().size() > 0) {
          for (String key : hashLoadResults.get(i).getOutlier().keySet()) {
            String[] newKey = key.split("\t");
            outliers.put(newKey[0] + "\t" + hashLoadResults.get(i).getSampleName() + "\t"
                         + newKey[1], hashLoadResults.get(i).getOutlier().get(key));
          }
        }
      }
      proj.getLog().reportTimeInfo("Writing outliers to " + outlierFileName);
      SerializedFiles.writeSerial(outliers, outlierFileName);
    }
    if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser")) {
      if (Files.exists(outlierFileName)) {
        outliers = (Hashtable<String, Float>) SerializedFiles.readSerial(outlierFileName);
      }
      proj.getLog()
          .reportTimeInfo("Writing outliers to " + proj.MARKER_DATA_DIRECTORY.getValue(false, true)
                          + "outliers.ser");
      SerializedFiles.writeSerial(outliers, proj.MARKER_DATA_DIRECTORY.getValue(false, true)
                                            + "outliers.ser");
    }
  }

  /**
   * Load data from Random Access Files organized by samples.
   */
  public static void testLoadTime(Project proj) {
    String[] samplesProj;
    RandomAccessFile file;
    byte[] readBuffer;
    long time1;

    samplesProj = proj.getSamples();
    time1 = new Date().getTime();
    try {
      // for (int i=0; i<samplesProj.length; i++) {
      for (int i = 0; i < 100; i++) {
        file = new RandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(true, true) + samplesProj[i]
                                    + SAMPLE_FILE_EXTENSION, "r");
        readBuffer = new byte[(int) file.length()]; // numMarkers * BYTES_PER_SAMPLE_MARKER
        file.read(readBuffer);
        file.close();
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    System.out.println("Finished loadToMarkerData from RAF in: " + (new Date().getTime() - time1)
                       + " ms");
  }

  public boolean hasBimodalBAF(byte chr, int startPosition, int endPosition) {
    DoubleVector dv = new DoubleVector();
    for (float baf : bafs) {
      if (baf > 0.15 && baf < 0.85) {
        dv.add((double) baf);
      }
    }
    return ArrayUtils.isBimodal(Doubles.toArray(dv));
  }

  public static void main(String[] args) {
    CLI cli = new CLI(Sample.class);

    cli.addArg("file", "SampRAF file");
    cli.addArg("proj", "Project properties file", false);
    cli.addFlag("out", "Write to standard out");

    cli.parseWithExit(args);

    Sample samp = Sample.loadFromRandomAccessFile(cli.get("file"), true, true, true, true, true);
    Project proj = cli.has("proj") ? new Project(cli.get("proj")) : null;

    if (cli.has("out")) {
      samp.writeToOut(proj == null ? null : proj.getMarkerNames());
    } else {
      samp.writeToFile(proj == null ? null : proj.getMarkerNames(),
                       (proj == null ? "./" : proj.PROJECT_DIRECTORY.getValue()) + cli.get("file")
                                                                    + ".txt");
    }
  }

}
