/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;

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
