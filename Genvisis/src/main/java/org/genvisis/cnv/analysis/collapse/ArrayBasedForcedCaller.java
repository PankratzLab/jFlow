/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;

/**
 *
 */
abstract class ArrayBasedForcedCaller<T extends Segment> extends ForcedCaller<T> {

  protected final int[][] indicesByChr;
  protected final boolean[] use;

  ArrayBasedForcedCaller(Project proj, String sample, int[][] indicesByChr, boolean[] use) {
    super(proj, sample);
    this.indicesByChr = indicesByChr;
    this.use = use;
  }

  /**
   * Extract project indices for each of these regions
   * 
   * @param markerDetailSet
   * @param regions
   * @param use these project indices will not used, pass null if all markers are desired
   * @return
   */
  <E extends Segment> List<int[]> getIndicesToCall(MarkerDetailSet markerDetailSet,
                                                   LocusSet<E> regions, boolean[] use) {
    List<int[]> indicesInProjectToCall = new ArrayList<>();
    Map<Marker, Integer> map = markerDetailSet.getMarkerIndexMap();
    for (E t : regions.getLoci()) {
      Stream<Marker> markers = markerDetailSet.viewMarkersInSeg(t);
      List<Integer> useIndices = new ArrayList<>();

      if (use != null) {
        markers.mapToInt(map::get).filter(i -> use[i]).forEach(useIndices::add);
      }
      indicesInProjectToCall.add(useIndices.stream().mapToInt(i -> i).toArray());
    }
    return indicesInProjectToCall;
  }

  /**
   * Extract project indices for each of these regions
   * 
   * @param markerDetailSet
   * @param regions
   * @param use these project indices will not used, pass null if all markers are desired
   * @return
   */
  int[][] getBackGroundIndices(boolean[] use, int[][] originals, Logger log) {

    int old = 0;
    int news = 0;
    int[][] backGroundIndices = new int[originals.length][];
    for (int i = 0; i < originals.length; i++) {
      List<Integer> newIndices = new ArrayList<>();
      for (int index : originals[i]) {
        old++;
        if (use == null || use[index]) {
          newIndices.add(index);
          news++;
        }
      }
      backGroundIndices[i] = newIndices.stream().mapToInt(z -> z).toArray();
    }
    log.reportTimeInfo("Background indices went from n = " + old + " to n = " + news);
    return backGroundIndices;
  }
}
