
package org.genvisis.cnv.annotator;

import org.genvisis.filesys.Segment;

/**
 * Interface for classes that compute and record {@link Annotation}s for individual {@link Segment}s
 */
public interface Annotator<A extends AnnotatedSegment> {

  /**
   * @param segments Collection of {@link AnnotatedSegment}s to annotate
   */
  void annotate(AnnotatedCollection<? extends A> segments);

  /**
   * @param segment Individual {@link AnnotatedSegment} to annotate
   */
  <T extends A> void annotate(T segment);

}
