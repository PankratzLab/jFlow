
package org.genvisis.cnv.annotator;

import java.util.Comparator;
import org.genvisis.filesys.Segment;
import com.google.common.collect.Ordering;

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

  /**
   * Override this method to customize the order in which to annotate segments. By default, original
   * collection order is preserved.
   *
   * @return An optional {@link Ordering} to apply to input segments to determine the order they
   *         will be accessed by the {@link #annotate} method.
   */
  default Ordering<A> inputOrdering() {
    return Ordering.from(new Comparator<A>() {

      @Override
      public int compare(A o1, A o2) {
        // All-equal comparison to preserve existing order
        return 0;
      }
    });
  }
}
