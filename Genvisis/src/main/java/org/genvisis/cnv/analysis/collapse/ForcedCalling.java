/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.common.Logger;
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

  /**
   * Extract project indices for each of these regions
   * 
   * @param markerDetailSet
   * @param regions
   * @param use these project indices will not used, pass null if all markers are desired
   * @return
   */
  default <E extends Segment> List<int[]> getIndicesToCall(MarkerDetailSet markerDetailSet,
                                                           LocusSet<E> regions, boolean[] use) {
    List<int[]> indicesInProjectToCall = new ArrayList<>();
    Map<Marker, Integer> map = markerDetailSet.getMarkerIndexMap();
    for (E t : regions.getLoci()) {
      Iterable<Marker> markers = markerDetailSet.viewMarkersInSeg(t);
      List<Integer> useIndices = new ArrayList<>();

      for (Marker marker : markers) {
        int index = map.get(marker);
        if (use != null && use[index]) {
          useIndices.add(index);
        }
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
  default int[][] getBackGroundIndices(boolean[] use, int[][] originals, Logger log) {

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
