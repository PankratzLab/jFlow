package org.genvisis.cnv.annotator;

import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.pankratzlab.common.filesys.Segment;
import com.google.common.collect.ImmutableSortedSet;

/**
 * Holds configuration for {@link Annotator}s to determine what region(s) to annotate
 */
public class AnnotatorConfig {

  private boolean doUpstream = false;
  private boolean doDownstream = false;
  private int basePairDistance;
  private int markerDistance;
  private boolean useSegmentSize = false;
  private MarkerDetailSet markers;

  /**
   * @return Whether or not to annotate region(s) upstream of the given segment
   */
  public boolean doUpstream() {
    return doUpstream;
  }

  /**
   * Flags this configuration to perform upstream annotation
   *
   * @return This {@link AnnotatorConfig} for chaining.
   */
  public AnnotatorConfig setDoUpstream() {
    this.doUpstream = true;
    // Return itself:
    return AnnotatorConfig.this;
  }

  /**
   * @return Whether or not to annotate region(s) downstream of the given segment
   */
  public boolean doDownstream() {
    return doDownstream;
  }

  /**
   * Flags this configuration to perform downstream annotation
   *
   * @return This {@link AnnotatorConfig} for chaining.
   */
  public AnnotatorConfig setDoDownstream() {
    this.doDownstream = true;
    // Return itself:
    return AnnotatorConfig.this;
  }

  /**
   * @return Fixed distance, in base pairs, to use to create up(down)stream segments
   */
  public int getBasePairDistance() {
    return basePairDistance;
  }

  /**
   * @param basePairDistance Desired distance to go up or downstream from a segment to create
   *          additional annotations
   * @return This {@link AnnotatorConfig} for chaining.
   */
  public AnnotatorConfig setBasePairDistance(int basePairDistance) {
    this.basePairDistance = basePairDistance;
    return AnnotatorConfig.this;
  }

  /**
   * @return True iff this configuration has marker data
   */
  public boolean hasMarkerData() {
    return markers != null;
  }

  /**
   * @param segment Reference segment
   * @param upstream Whether to go up or down stream in markers
   * @return The distance, in base pairs, from traveling a set number of markers up or downstream of
   *         the given {@link Segment}
   */
  public int getMarkerDistance(Segment segment, boolean upstream) {

    // If no marker data provided:
    if (!hasMarkerData()) {
      throw new IllegalStateException("Can't get marker distance: No MarkerDetailSet available.");
    }

    // The set of markers in a segment:
    ImmutableSortedSet<Marker> segmentMarkers = markers.getMarkersInSeg(segment);

    // Check for empty segment:
    if (segmentMarkers.isEmpty()) {
      System.out.println("Segment " + segment.toString() + " has no markers.");
      return -1;
    }

    int basePairOffset;
    int flankEndMarkerIndex;
    Marker referenceMarker = null;

    if (upstream) {
      // For upstream, get first marker:
      referenceMarker = segmentMarkers.iterator().next();

      flankEndMarkerIndex = markers.getMarkerIndexMap().get(referenceMarker);
      flankEndMarkerIndex -= markerDistance;
      // If you go over markers and there's none left in that direction (ensures this marker index
      // is valid):
      flankEndMarkerIndex = Math.max(flankEndMarkerIndex, 0);

    }
    // for downstream, get last marker (can alter markerDetailSet to avoid looping here):
    else {
      for (Marker next : segmentMarkers) {
        referenceMarker = next;
      }

      flankEndMarkerIndex = markers.getMarkerIndexMap().get(referenceMarker);
      flankEndMarkerIndex += markerDistance;
      // If you go over markers and there's none left in that direction (ensures this marker index
      // is valid):
      flankEndMarkerIndex = Math.min(flankEndMarkerIndex, markers.markersAsList().size() - 1);
    }

    basePairOffset = Math.abs(markers.markersAsList().get(flankEndMarkerIndex).getPosition()
                              - referenceMarker.getPosition());

    return basePairOffset;
  }

  /**
   * @param markerDistance How many markers to go up or downstream when annotating
   * @return This {@link AnnotatorConfig} for chaining.
   */
  public AnnotatorConfig setMarkerDistance(MarkerDetailSet markers, int markerDistance) {

    this.markers = markers;
    this.markerDistance = markerDistance;

    // Return itself:
    return AnnotatorConfig.this;

  }

  /**
   * @return the useSegmentSize
   */
  public boolean useSegmentSize() {
    return useSegmentSize;
  }

  /**
   * @param useSegmentSize the useSegmentSize to set
   */
  public AnnotatorConfig setUseSegmentSize() {
    this.useSegmentSize = true;
    // Return itself:
    return AnnotatorConfig.this;
  }
}
