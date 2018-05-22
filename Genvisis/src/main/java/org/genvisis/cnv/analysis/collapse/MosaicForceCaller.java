/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.analysis.MosaicismDetect;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicBuilder;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

/**
 * Force a mosaic call for a region
 */
class MosaicForceCaller extends MarkerBasedForcedCaller<MosaicRegion> {

  private final Map<Marker, Double> bafs;
  private LocusSet<Segment> distributionalExcludes;

  MosaicForceCaller(Project proj, String sample, Map<Marker, Double> bafs, Set<Marker> use) {
    super(proj, sample, use);
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
    List<Marker> allMarkers = proj.getMarkerSet().getMarkers();
    Set<Marker> nonCNMarkers = allMarkers.stream()
                                         .filter(m -> !proj.ARRAY_TYPE.getValue()
                                                                      .isCNOnly(m.getName()))
                                         .collect(ImmutableSet.toImmutableSet());
    proj.getLog()
        .reportTimeWarning("Starting with " + (use == null ? 0 : allMarkers.size() - use.size())
                           + " pre-defined markers to remove of " + nonCNMarkers.size()
                           + " markers total");

    Set<Marker> nonCNUseMarkers = use == null ? nonCNMarkers : Sets.intersection(use, nonCNMarkers);
    builderMosaic.use(nonCNUseMarkers);
    proj.getLog()
        .reportTimeWarning("Removed " + (allMarkers.size() - nonCNMarkers.size())
                           + " copy number only markers, leaving " + nonCNUseMarkers.size()
                           + " markers for analysis");

    MosaicismDetect md = builderMosaic.build(proj, sample, bafs);
    List<MosaicRegion> regions = new ArrayList<>();
    for (V v : set.getLoci()) {
      regions.add(md.callMosaic(v, true).getLoci()[0]);
    }
    this.results = new LocusSet<>(regions, true, proj.getLog());
  }

}
