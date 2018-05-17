/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

/**
 * 
 *
 */
interface ForcedCalling<T extends Segment> {

  /**
   * @param builder the builder should at least contain FID/IID
   */
  <V extends Segment> void forceCall(LocusSet<V> set);

  /**
   * @return
   */
  LocusSet<T> getResults();

}
