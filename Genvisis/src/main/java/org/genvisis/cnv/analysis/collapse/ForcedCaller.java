/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

/**
 *
 */
abstract class ForcedCaller<T extends Segment> implements ForcedCalling<T> {

  protected Project proj;
  protected String sample;
  protected int[][] indicesByChr;
  protected LocusSet<T> results;
  protected boolean[] use;

  ForcedCaller(Project proj, String sample, int[][] indicesByChr, boolean[] use) {
    super();
    this.proj = proj;
    this.sample = sample;
    this.indicesByChr = indicesByChr;
    this.use = use;
  }

  @Override
  public LocusSet<T> getResults() {
    return results;
  }
}
