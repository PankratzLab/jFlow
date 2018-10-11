package org.genvisis.cnv.analysis.collapse;

import java.util.Set;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.pankratzlab.common.filesys.Segment;
import org.genvisis.cnv.filesys.Project;

public abstract class MarkerBasedForcedCaller<T extends Segment> extends ForcedCaller<T> {

  protected final Set<Marker> use;

  MarkerBasedForcedCaller(Project proj, String sample, Set<Marker> use) {
    super(proj, sample);
    this.use = use;
  }

}
