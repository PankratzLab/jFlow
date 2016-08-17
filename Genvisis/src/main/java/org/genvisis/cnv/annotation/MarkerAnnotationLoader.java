package org.genvisis.cnv.annotation;

import java.util.Hashtable;
import java.util.List;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.filesys.Segment;

/**
 * @author lane0212 Class that concentrates on loading annotations for specific markers
 */
public class MarkerAnnotationLoader extends AnnotationFileLoader {
  private final MarkerSet markerSet;
  private final Hashtable<String, Integer> indices;
  private final byte[] chrs;
  private final int[] pos;
  private final Hashtable<String, Integer> markerIndices;

  /**
   * @param proj
   * @param annotationFilename
   * @param markerSet will be loaded if null, can save time if it is already available
   * @param indexRequired should always be true for now
   */
  public MarkerAnnotationLoader(Project proj, AnalysisParams[] params, String annotationFilename,
                                MarkerSet markerSet, boolean indexRequired) {
    super(proj, params, null, annotationFilename, indexRequired);
    if (markerSet == null) {
      markerSet = proj.getMarkerSet();
    }
    indices = proj.getMarkerIndices();
    this.markerSet = markerSet;
    chrs = markerSet.getChrs();
    pos = markerSet.getPositions();
    markerIndices = proj.getMarkerIndices();
  }

  /**
   * @param markers
   * @param parsersQueries typically each entry in the {@link AnnotationParser} array represents a
   *        single marker
   */
  public void fillAnnotations(final String[] markers, List<AnnotationParser[]> parsersQueries,
                              QUERY_ORDER qOrder) {
    Segment[] markerSegments = null;
    if (markers == null && qOrder != QUERY_ORDER.NO_ORDER) {
      proj.getLog()
          .reportTimeWarning("No marker names were provided, searching entire annotation file");
    } else {
      markerSegments = getSegmentsForMarkers(markers);
    }
    query(markerSegments, parsersQueries, qOrder);
  }

  public Hashtable<String, Integer> getIndices() {
    return indices;
  }

  public MarkerSet getMarkerSet() {
    return markerSet;
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
