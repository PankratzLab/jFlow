package org.genvisis.cnv.annotator;

import java.util.ArrayList;
import java.util.List;
import org.pankratzlab.shared.filesys.Segment;

/**
 * Abstract {@link Annotator} implementation. Takes care of virtually all boilerplate so concrete
 * implementations only have to provide logic for actually computing a particular annotation.
 */
public abstract class AbstractAnnotator<A extends AnnotatedSegment> implements Annotator<A> {

  public static final String UPSTREAM = ".upstream";
  public static final String DOWNSTREAM = ".downstream";
  public static final String SELF = ".self";
  public static final String BP = ".bp";
  public static final String MARKER = ".marker";

  private AnnotatorConfig config;

  @Override
  public void annotate(AnnotatedCollection<? extends A> segments) {
    // Pull in the configuration from the AnnotatedCollection
    setConfig(segments.getConfig());

    // Annotate in sorted order
    segments.stream().sorted(inputOrdering()).forEach(this::annotate);

    // Record the annotation
    segments.addAnnotator(this);

    cleanUp();
  }

  @Override
  public <T extends A> void annotate(T segment) {
    List<LabeledOffset> labeledOffsets = new ArrayList<>();

    // Always annotate the segment itself
    labeledOffsets.add(new LabeledOffset(segment.getStart(), segment.getStop(), ""));

    // Use the configuration to figure out the actual ranges to annotate
    addOffsets(labeledOffsets, BP, config.doUpstream(), config.doDownstream(),
               config.getBasePairDistance(), segment.getStart(), segment.getStop());

    if (config.useSegmentSize()) {
      addOffsets(labeledOffsets, SELF, config.doUpstream(), config.doDownstream(),
                 segment.getLength(), segment.getStart(), segment.getStop());
    }

    if (config.hasMarkerData()) {
      int offset = config.getMarkerDistance(segment, false);

      addOffsets(labeledOffsets, MARKER, config.doUpstream(), config.doDownstream(), offset,
                 segment.getStart(), segment.getStop());
    }

    // Calculate all requested annotations for this segment
    for (LabeledOffset o : labeledOffsets) {
      calcAnnotations(segment, o.start, o.stop, o.suffix);
    }

  }

  /**
   * Helper method to generate {@link LabeledOffset}s
   */
  private void addOffsets(List<LabeledOffset> labeledOffsets, String type, boolean doUpstream,
                          boolean doDownstream, int offset, int start, int stop) {

    // Making sure CNV distance is valid
    if (offset <= 0) {
      return;
    }

    if (doUpstream) {
      labeledOffsets.add(new LabeledOffset(start - offset, start, type + UPSTREAM));
    }

    if (doDownstream) {
      labeledOffsets.add(new LabeledOffset(stop, stop + offset, type + DOWNSTREAM));
    }
  }

  /**
   * Set the {@link AnnotatorConfig} to use for annotation
   */
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

  /**
   * Optional hook for post-{@link Segment} annotation, providing the option to clean up cached
   * resources.
   */
  protected void cleanUp() {
    // NO-OP
  }

  /**
   * Calculate and apply annotation(s) to the specified range (via
   * {@link AnnotatedSegment#putAnnotation}
   */
  protected abstract void calcAnnotations(A segmentToAnnotate, int start, int stop, String suffix);

  /**
   * Helper class to keep offset metadata together
   */
  private static class LabeledOffset {

    public LabeledOffset(int start, int stop, String suffix) {
      this.start = start;
      this.stop = stop;
      this.suffix = suffix;
    }

    private final int start;
    private final int stop;
    private final String suffix;
  }
}
