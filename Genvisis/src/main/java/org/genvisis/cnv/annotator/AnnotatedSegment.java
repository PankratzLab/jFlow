// 3)

package org.genvisis.cnv.annotator;

import java.util.Collection;
import org.pankratzlab.common.filesys.Segment;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;

/**
 * A {@link Segment} that can hold any number of {@link Annotation}s
 */
public class AnnotatedSegment extends Segment {

  private static final long serialVersionUID = 1L;

  // Map of Annotators to a set of their annotations applied to this segment
  private Multimap<Annotator<? extends AnnotatedSegment>, Annotation> annotations = MultimapBuilder.hashKeys()
                                                                                                   .linkedHashSetValues()
                                                                                                   .build();

  /**
   * @see Segment#Segment(byte, int, int)
   */
  public AnnotatedSegment(byte chr, int start, int stop) {
    super(chr, start, stop);
  }

  /**
   * @see Segment#Segment(String)
   */
  public AnnotatedSegment(String ucscLocation) {
    super(ucscLocation);
  }

  /**
   * Convenience method for {@link Annotation#Annotation(String, String)}
   *
   * @see #putAnnotation(Annotator, Annotation)
   */
  public void putAnnotation(Annotator<? extends AnnotatedSegment> annotator, String annoName,
                            String annoValue) {
    putAnnotation(annotator, new Annotation(annoName, annoValue));
  }

  /**
   * Convenience method for {@link Annotation#Annotation(String, String, String)}
   *
   * @see #putAnnotation(Annotator, Annotation)
   */
  public void putAnnotation(Annotator<? extends AnnotatedSegment> annotator, String annoName,
                            String annoValue, String annoSuffix) {
    putAnnotation(annotator, new Annotation(annoName, annoValue, annoSuffix));
  }

  /**
   * Record an annotation as provided by a particular {@link Annotator}.
   */
  public void putAnnotation(Annotator<? extends AnnotatedSegment> annotator,
                            Annotation annotation) {
    this.annotations.put(annotator, annotation);
  }

  // FIXME change to annotator class?
  /**
   * @return All {@link Annotation}s on this segment created by the specified {@link Annotator}
   */
  public ImmutableList<Annotation> getAnnotations(Annotator<? extends AnnotatedSegment> annotator) {
    // Ensuring immutablility (anonymous object returned from get bc never had a pointer):
    return ImmutableList.copyOf(annotations.get(annotator));
  }

  /**
   * @return All {@link Annotation}s on this segment
   */
  public ImmutableList<Annotation> getAnnotations() {
    return ImmutableList.copyOf(annotations.values());
  }

  /**
   * @see #getAnnotation(Annotator, String)
   */
  public Annotation getAnnotation(Annotator<? extends AnnotatedSegment> annotator,
                                  Annotation annotation) {
    String annoName = annotation.getAnnotationName();
    //
    return getAnnotation(annotator, annoName);
  }

  /**
   * @return The {@link Annotation} on this segment created by the specified {@link Annotator}, with
   *         the given annotation name
   */
  public Annotation getAnnotation(Annotator<? extends AnnotatedSegment> annotator,
                                  String annoName) {

    // 1) Get the set of annotations for the given annotator
    Collection<Annotation> annoList = annotations.get(annotator);

    // 2) Then see if it contains an annotation with annoName:
    for (Annotation annoInstance : annoList) {

      // 3) If so return that annotation
      if (annoInstance.hasName(annoName)) {
        return annoInstance;
      }
    }

    return null;
  }
}
