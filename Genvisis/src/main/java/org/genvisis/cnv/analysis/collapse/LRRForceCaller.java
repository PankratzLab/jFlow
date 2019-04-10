/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.analysis.collapse.LRRForceCaller.LRRRegion;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.CNVariant.CNVBuilder;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Segment;

/**
 * Forced "calls" using LRR summary data
 */
public class LRRForceCaller extends ArrayBasedForcedCaller<LRRRegion> {

  private static final double MAD_FACTOR = 1.4826;
  private double[] lrrs;

  /**
   * @param proj
   * @param sample
   * @param indicesByChr
   * @param use
   */
  LRRForceCaller(Project proj, String sample, int[][] indicesByChr, double[] lrrs, boolean[] use) {
    super(proj, sample, indicesByChr, use);
    this.lrrs = lrrs;
  }

  static class LRRRegion extends CNVariant {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    private double lrrRegionMedian;
    private double lrrRegionMAD;
    private double autosomalLrrMedian;
    private double autosomalLrrMAD;

    LRRRegion(CNVBuilder builder, double lrrRegionMedian, double lrrRegionMAD,
              double autosomalLrrMedian, double autosomalLrrMAD) {
      super(builder);
      this.lrrRegionMedian = lrrRegionMedian;
      this.lrrRegionMAD = lrrRegionMAD;
      this.autosomalLrrMedian = autosomalLrrMedian;
      this.autosomalLrrMAD = autosomalLrrMAD;
    }

    @Override
    public String toAnalysisString() {
      return super.toAnalysisString() + "\t" + lrrRegionMedian + "\t" + lrrRegionMAD + "\t"
             + autosomalLrrMedian + "\t" + autosomalLrrMAD;
    }

    @Override
    public String[] getHeader() {
      return ArrayUtils.concatAll(super.getHeader(),
                                  new String[] {"LRR_REGION_MEDIAN", "LRR_REGION_MAD",
                                                "AUTOSOMAL_LRR_MEDIAN", "AUTOSOMAL_LRR_MAD"});
    }

  }

  /*
   * (non-Javadoc)
   * 
   * @see org.genvisis.cnv.analysis.collapse.ForcedCalling#forceCall(org.genvisis.filesys.LocusSet)
   */
  @Override
  public <V extends Segment> void forceCall(LocusSet<V> set) {
    int[] autosomalIndices = proj.getAutosomalMarkerIndices();
    List<Integer> tmp = new ArrayList<>();
    List<int[]> indicesToUse = getIndicesToCall(proj.getMarkerSet(), set, use);
    for (int i : autosomalIndices) {
      if (use == null || use[i]) {
        tmp.add(i);
      }
    }
    double[] sub = ArrayUtils.subArray(lrrs, tmp.stream().mapToInt(i -> i).toArray());
    double autosomalLrrMedian = ArrayUtils.median(sub);
    double autosomalLrrMAD = ArrayUtils.madFactor(sub, MAD_FACTOR);

    List<LRRRegion> lrrVariants = new ArrayList<>();
    for (int i = 0; i < set.getLoci().length; i++) {
      int[] current = indicesToUse.get(i);
      double[] currentLrrs = ArrayUtils.subArray(lrrs, current);
      double lrrRegionMedian = ArrayUtils.median(currentLrrs);
      double lrrRegionMAD = ArrayUtils.madFactor(currentLrrs, MAD_FACTOR);

      CNVBuilder builder = new CNVBuilder();
      V v = set.getLoci()[i];
      builder.chr(v.getChr());
      builder.start(v.getStart());
      builder.stop(v.getStop());
      builder.familyID(sample);
      builder.individualID(sample);
      builder.score(Double.NaN);
      lrrVariants.add(new LRRRegion(builder, lrrRegionMedian, lrrRegionMAD, autosomalLrrMedian,
                                    autosomalLrrMAD));
    }
    this.results = new LocusSet<>(lrrVariants, true, proj.getLog());
  }
}
