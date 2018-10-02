package org.genvisis.cnv.analysis.collapse;

import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;

public abstract class ForcedCaller<T extends Segment> implements ForcedCalling<T> {

  protected final Project proj;
  protected final String sample;
  protected LocusSet<T> results;

  /**
   * @param proj
   * @param sample
   */
  protected ForcedCaller(Project proj, String sample) {
    super();
    this.proj = proj;
    this.sample = sample;
  }

  @Override
  public LocusSet<T> getResults() {
    return results;
  }

}
