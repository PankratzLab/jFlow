
package org.genvisis.cnv.annotator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.CNVariant;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;
import org.genvisis.stats.ProbDist;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Ordering;
import com.google.common.collect.Table;

/**
 * {@link Annotator} that computes {@link BeastScore} annotations
 */
public class BeastScoreAnnotator extends AbstractAnnotator<CNVariant> {

  private static final String BEAST_SCORE_ANNO = "BEAST_SCORE";
  private static final String BEAST_HEIGHT_ANNO = "BEAST_HEIGHT";
  private static final String BEAST_LENGTH_ANNO = "BEAST_LENGTH";
  public static final double DEFAULT_MAD_SCALE_FACTOR = 0.134894516;
  public static final double DEFAULT_ALPHA = 0.5;
  private final SampleData sd;
  private final MarkerDetailSet md;
  private final Project proj;
  private final Map<String, Table<Byte, Marker, Double>> transformedBySample;
  private final double madScale;
  private final double alpha;

  public BeastScoreAnnotator(Project proj) {
    this(proj, DEFAULT_MAD_SCALE_FACTOR);
  }

  public BeastScoreAnnotator(Project proj, double scaleFactor) {
    this(proj, scaleFactor, DEFAULT_ALPHA);
  }

  public BeastScoreAnnotator(Project proj, double scaleFactor, double alpha) {
    // The project is used to load lrr data
    this.proj = proj;

    madScale = scaleFactor == 0 ? 1 : scaleFactor;
    this.alpha = alpha;
    transformedBySample = new HashMap<>();

    // We use the SampleData to look up sample IDs
    sd = proj.getSampleData(false);

    // And we need to use the markers
    md = proj.getMarkerSet();

  }

  @Override
  public Ordering<CNVariant> inputOrdering() {
    return Ordering.from(new SampleCNVariantComp());
  }

  @Override
  protected void calcAnnotations(CNVariant segmentToAnnotate, int start, int stop, String suffix) {
    String dna = sd.lookupDNA(segmentToAnnotate.getFamilyID() + "\t"
                              + segmentToAnnotate.getIndividualID());

    Byte targetChr = segmentToAnnotate.getChr();
    Table<Byte, Marker, Double> transformedMarkerData = null;

    // Check if we have a cached table for this individual
    if (transformedBySample.containsKey(dna)) {
      transformedMarkerData = transformedBySample.get(dna);
    } else {
      transformedBySample.clear();
      transformedMarkerData = HashBasedTable.create();
      transformedBySample.put(dna, transformedMarkerData);
    }

    // If we don't have transformed data for this chromosome, we do the transformation
    if (!transformedMarkerData.containsRow(targetChr)) {
      Sample samp = proj.getPartialSampleFromRandomAccessFile(dna, true, false, false, true, false);

      // Create a list of markers for this chromosome, sorted by their LRR value
      List<MarkerLRR> markersByLrr = md.getChrMap().get(targetChr).stream()
                                       .map(m -> new MarkerLRR(m, samp.getLRRs(),
                                                               md.getMarkerIndexMap()))
                                       .filter(MarkerLRR::notNaN)
                                       .sorted(BeastScoreAnnotator::compareLRR)
                                       .collect(Collectors.toList());

      // Inverse transform the quantiles of the sorted data
      List<MarkerTrans> markersTransformed = new ArrayList<>();
      double count = (double) (markersByLrr.size() + 1);

      for (int i = 0; i < markersByLrr.size(); i++) {
        MarkerLRR markerLrr = markersByLrr.get(i);
        double transformed = tDistReverse5df((i + 1) / count);
        if (!Double.isNaN(transformed)) {
          markersTransformed.add(new MarkerTrans(markerLrr.m, transformed));
        }
      }

      // Compute the median of the transformed data
      Collections.sort(markersTransformed, BeastScoreAnnotator::compareAbsTrans);
      double indicesMADScaled = ArrayUtils.medianSorted(markersTransformed) / madScale;

      // Record the scaled transformed data for each marker
      for (MarkerTrans trans : markersTransformed) {
        transformedMarkerData.put(targetChr, trans.m, trans.val / indicesMADScaled);
      }
    }

    // Get the transformed marker map for this chromosome
    final Map<Marker, Double> transformedRow = transformedMarkerData.row(targetChr);

    Stream<Marker> cnvMarkers = md.viewMarkersInSeg(segmentToAnnotate.getChr(), start, stop);

    // Get the transformed values for each marker in this CNV
    final boolean isNGS = ARRAY.NGS.equals(proj.ARRAY_TYPE.getValue());

    // TODO This NGS check seems overly complex
    List<Double> cnvMarkerVals = cnvMarkers.filter(m -> transformedRow.containsKey(m)
                                                        && !(isNGS
                                                             && NGS_MARKER_TYPE.VARIANT_SITE.equals(NGS_MARKER_TYPE.getType(m))))
                                           .map(transformedRow::get).sorted()
                                           .collect(Collectors.toList());

    // Score this CNV
    double beastHeight;
    double beastLength;
    double beastScore;
    if (cnvMarkerVals.isEmpty()) {
      beastHeight = 0.0;
      beastLength = 0.0;
      beastScore = 0.0;
    } else {
      beastHeight = ArrayUtils.medianSorted(cnvMarkerVals);
      beastLength = cnvMarkerVals.size();
      beastScore = Math.abs(Math.pow(beastLength, alpha) * beastHeight);
    }

    String beastSuffix = suffix + "_a" + alpha + "_s" + madScale;

    addAnnotation(segmentToAnnotate, BEAST_SCORE_ANNO, beastScore, beastSuffix);
    addAnnotation(segmentToAnnotate, BEAST_HEIGHT_ANNO, beastHeight, beastSuffix);
    addAnnotation(segmentToAnnotate, BEAST_LENGTH_ANNO, beastLength, beastSuffix);
  }

  @Override
  protected void cleanUp() {
    transformedBySample.clear();
  }

  private static double tDistReverse5df(double q) {
    int sign = 1;
    if (q < 0.5) {
      sign = -1;
    } else {
      q = 1 - q;
    }
    return ProbDist.TDistReverse(q * 2, 5) * sign;
  }

  /**
   * Comparator for sorting {@link MarkerLRR}s based on their LRRs
   */
  private static int compareLRR(MarkerLRR m1, MarkerLRR m2) {
    return Float.compare(m1.val, m2.val);
  }

  /**
   * Comparator for sorting {@link MarkerTrans}s based on the absolute value of their transformed
   * values
   */
  private static int compareAbsTrans(MarkerTrans m1, MarkerTrans m2) {
    return Double.compare(m1.doubleValue(), m2.doubleValue());
  }

  /**
   * Simple wrapper to keep markers and LRRs together
   */
  private static class MarkerLRR {

    private final Marker m;
    private float val;

    private MarkerLRR(Marker m, float[] lrrs, Map<Marker, Integer> markerIndexMap) {
      this.m = m;
      this.val = lrrs[markerIndexMap.get(m)];
    }

    private boolean notNaN() {
      return !Float.isNaN(val);
    }
  }

  /**
   * Simple wrapper to keep markers and transformed indices together. Extends {@link Number} for
   * purposes of computing median.
   */
  private static class MarkerTrans extends Number {

    private static final long serialVersionUID = 1L;

    private final Marker m;
    private final double val;

    private MarkerTrans(Marker m, double val) {
      this.m = m;
      this.val = val;
    }

    @Override
    public int intValue() {
      return (int) val;
    }

    @Override
    public long longValue() {
      return (long) val;
    }

    @Override
    public float floatValue() {
      return (float) val;
    }

    @Override
    public double doubleValue() {
      // NB: used for median computation
      return Math.abs(val);
    }
  }

  /**
   * Helper {@link Comparator} to collection {@link CNVariants} with the same DNA ids together.
   */
  private class SampleCNVariantComp implements Comparator<CNVariant> {

    @Override
    public int compare(CNVariant v1, CNVariant v2) {
      String dna1 = sd.lookupDNA(v1.getFamilyID() + "\t" + v1.getIndividualID());
      String dna2 = sd.lookupDNA(v2.getFamilyID() + "\t" + v2.getIndividualID());
      int c = dna1.compareTo(dna2);
      if (c == 0) {
        c = v1.compareTo(v2);
      }
      return c;
    }

  }
}
