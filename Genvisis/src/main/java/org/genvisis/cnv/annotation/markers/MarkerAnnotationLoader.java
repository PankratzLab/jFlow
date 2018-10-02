package org.genvisis.cnv.annotation.markers;

import java.util.List;
import java.util.Map;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.filesys.Segment;

/**
 * @author lane0212 Class that concentrates on loading annotations for specific markers
 */
public class MarkerAnnotationLoader extends AnnotationFileLoader {

  private final MarkerDetailSet markerSet;
  private final byte[] chrs;
  private final int[] pos;
  private final Map<String, Integer> markerIndices;

  /**
   * @param annotationFilename
   * @param markerSet MarkerSet to load
   * @param indexRequired should always be true for now
   * @param log
   */
  public MarkerAnnotationLoader(AnalysisParams[] params, String annotationFilename,
                                MarkerDetailSet markerSet, boolean indexRequired, Logger log) {
    super(params, null, annotationFilename, indexRequired, log);
    this.markerSet = markerSet;
    chrs = markerSet.getChrs();
    pos = markerSet.getPositions();
    markerIndices = markerSet.getMarkerIndices();
  }

  public MarkerSetInfo getMarkerSet() {
    return markerSet;
  }

  public Map<String, Integer> getMarkerIndices() {
    return markerIndices;
  }

  /**
   * @param markers
   * @param parsersQueries typically each entry in the {@link AnnotationParser} array represents a
   *          single marker
   */
  public void fillAnnotations(final String[] markers,
                              List<Map<String, ? extends AnnotationParser>> parsersQueries,
                              QUERY_TYPE qOrder) {
    Segment[] markerSegments = null;
    if (markers != null || qOrder == QUERY_TYPE.DISCRETE_LIST) {
      markerSegments = getSegmentsForMarkers(markers);
    }
    query(markerSegments, parsersQueries, qOrder);
  }

  private Segment[] getSegmentsForMarkers(final String[] markers) {
    Segment[] segs = new Segment[markers.length];
    for (int i = 0; i < segs.length; i++) {
      int markerIndex = markerIndices.get(markers[i]);
      Segment markerSeg = new Segment(chrs[markerIndex], pos[markerIndex], pos[markerIndex]);
      segs[i] = markerSeg;
    }
    return segs;
  }

}
