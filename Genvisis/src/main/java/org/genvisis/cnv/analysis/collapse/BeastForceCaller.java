/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.analysis.BeastScore.BeastVariant;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.CNVariant.CNVBuilder;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;

/**
 * 
 *
 */
class BeastForceCaller extends ArrayBasedForcedCaller<BeastVariant> {

  private double[] lrrs;

  public BeastForceCaller(Project proj, String sample, int[][] indicesByChr, double[] lrrs,
                          boolean[] use) {
    super(proj, sample, indicesByChr, use);
    this.lrrs = lrrs;
  }

  /*
   * (non-Javadoc)
   * @see org.genvisis.cnv.analysis.collapse.ForcedCalling#forceCall()
   */
  @Override
  public <V extends Segment> void forceCall(LocusSet<V> set) {

    List<int[]> tmp = getIndicesToCall(proj.getMarkerSet(), set, use);
    int[][] indicesForScores = new int[tmp.size()][];
    for (int i = 0; i < indicesForScores.length; i++) {
      indicesForScores[i] = tmp.get(i);
    }

    BeastScore beastScore = new BeastScore(ArrayUtils.toFloatArray(lrrs),
                                           getBackGroundIndices(use, indicesByChr, proj.getLog()),
                                           indicesForScores, proj.getLog());
    beastScore.computeBeastScores();
    if (set.getLoci().length != beastScore.getBeastHeights().length) {
      throw new IllegalStateException("Invalid beast score computation");
    }
    List<BeastVariant> beastVariants = new ArrayList<>();
    for (int i = 0; i < set.getLoci().length; i++) {
      CNVBuilder builder = new CNVBuilder();
      V v = set.getLoci()[i];
      builder.chr(v.getChr());
      builder.start(v.getStart());
      builder.stop(v.getStop());
      builder.familyID(sample);
      builder.individualID(sample);
      builder.score(Double.NaN);
      beastVariants.add(new BeastVariant(builder, beastScore.getBeastScores()[i],
                                         beastScore.getBeastLengths()[i],
                                         beastScore.getBeastHeights()[i]));
    }

    this.results = new LocusSet<>(beastVariants, true, proj.getLog());

  }

}
