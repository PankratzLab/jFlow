// 2)

package org.genvisis.cnv.annotator;

import java.util.ArrayList;
import java.util.List;

/*
 * 
 * 
 * 
 */
public abstract class AbstractAnnotator<A extends AnnotatedSegment> implements Annotator<A> {

  private AnnotatorConfig config;
  public static final String UPSTREAM = ".upstream";
  public static final String DOWNSTREAM = ".downstream";

  // Annotate collection:
  @Override
  public void annotate(AnnotatedCollection<? extends A> segments) {
    // Get configuration, to use when annotating segment:
    this.setConfig(segments.getConfig());
    // For each segment of type A in the segments collection...:
    for (A segment : segments) {
      annotate(segment);
    }
    segments.addAnnotator(this);
  }

  // Annotate individual segment:
  @Override
  public <T extends A> void annotate(T segment) {

    // List of labeled offsets:
    List<LabeledOffset> labeledOffsets = new ArrayList<>();

    labeledOffsets.add(new LabeledOffset(segment.getStart(), segment.getStop(), ""));

    // Use the configuration to figure out the actual ranges to annotate:
    // Helper method:
    addOffset(labeledOffsets, config.doUpstream(), config.doDownstream(),
              config.getBasePairDistance(), segment.getStart(), segment.getStop());

    //
    if (config.hasMarkerData()) {
      int offset = config.getMarkerDistance(segment, false);

      //
      addOffset(labeledOffsets, config.doUpstream(), config.doDownstream(), offset,
                segment.getStart(), segment.getStop());
    }

    //
    for (LabeledOffset o : labeledOffsets) {
      // Call below method: Compute the annotation and assigning to the segment (a field of the
      // segment object):
      calcAnnotations(segment, o.start, o.stop, o.suffix);
    }

  }

  // Helper method:
  private void addOffset(List<LabeledOffset> labeledOffsets, boolean doUpstream,
                         boolean doDownstream, int offset, int start, int stop) {

    // Making sure CNV distance is valid (>0):
    if (offset <= 0) {
      return;
    }
    // Upstream:
    if (doUpstream) {
      labeledOffsets.add(new LabeledOffset(start - offset, start, UPSTREAM));
    }
    // Downstream:
    if (doDownstream) {
      labeledOffsets.add(new LabeledOffset(stop, stop + offset, DOWNSTREAM));
    }
  }

  // Set configuration:
  private void setConfig(AnnotatorConfig config) {
    this.config = config;
  }

  /**
   * Helper method to attach a particular annotation to a segment. Call this from
   * {@link #calcAnnotations(AnnotatedSegment, int, int, String)}
   */
  protected void addAnnotation(A segmentToAnnotate, String annoName, Object annoValue,
                               String annoSuffix) {
    segmentToAnnotate.putAnnotation(this, annoName, annoValue.toString(), annoSuffix);
  }

  // Declare: Abstract method to give segment and offset to the annotator, to compute the annotation
  // value for that annotator:
  // suffix for renaming up/down stream:
  protected abstract void calcAnnotations(A segmentToAnnotate, int start, int stop, String suffix);

  // Helper class to keep offset and suffix together:
  private static class LabeledOffset {

    public LabeledOffset(int start, int stop, String suffix) {
      super();
      this.start = start;
      this.stop = stop;
      this.suffix = suffix;
    }

    private final int start;
    private final int stop;
    private final String suffix;
  }
}
