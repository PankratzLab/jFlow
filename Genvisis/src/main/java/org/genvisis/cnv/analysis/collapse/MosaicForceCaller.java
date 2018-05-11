/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.analysis.MosaicismDetect;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicBuilder;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

/**
 * Force a mosaic call for a region
 */
class MosaicForceCaller extends ForcedCaller<MosaicRegion> {

  private final double[] bafs;
  private LocusSet<Segment> distributionalExcludes;

  MosaicForceCaller(Project proj, String sample, int[][] indicesByChr, double[] bafs,
                    boolean[] use) {
    super(proj, sample, indicesByChr, use);
    this.bafs = bafs;
    this.distributionalExcludes = null;
  }

  void setDistributionalExcludes(LocusSet<Segment> distributionalExcludes) {
    this.distributionalExcludes = distributionalExcludes;
  }

  @Override
  public <V extends Segment> void forceCall(LocusSet<V> set) {
    MosaicBuilder builderMosaic = new MosaicBuilder();
    builderMosaic.verbose(true);
    builderMosaic.use(use);
    boolean[] cns = proj.getCNMarkers();
    boolean[] tmp = ArrayUtils.booleanNegative(cns);
    proj.getLog()
        .reportTimeWarning("Starting with "
                           + (use == null ? 0
                                          : ArrayUtils.booleanArraySum(ArrayUtils.booleanNegative(use)))
                           + " pre-defined markers to remove of " + cns.length + " markers total");

    for (int i = 0; i < tmp.length; i++) {
      tmp[i] = tmp[i] && (use == null || use[i]);
    }
    builderMosaic.use(tmp);
    proj.getLog()
        .reportTimeWarning("Removed " + ArrayUtils.booleanArraySum(cns)
                           + " copy number only markers, leaving " + ArrayUtils.booleanArraySum(tmp)
                           + " markers for analysis");

    builderMosaic.indicesByChr(indicesByChr);
    MosaicismDetect md = builderMosaic.build(proj, sample, proj.getMarkerSet(), bafs);
    List<MosaicRegion> regions = new ArrayList<>();
    for (V v : set.getLoci()) {
      regions.add(md.callMosaic(v, true).getLoci()[0]);
    }
    this.results = new LocusSet<>(regions, true, proj.getLog());
  }

}
