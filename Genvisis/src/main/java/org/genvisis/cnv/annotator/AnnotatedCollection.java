package org.genvisis.cnv.annotator;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.StringJoiner;
import org.pankratzlab.common.Files;
import com.google.common.collect.ImmutableList;

/**
 * A collection of {@link AnnotatedSegment}s, to be consumed by {@link Annotator}s.
 */
public class AnnotatedCollection<A extends AnnotatedSegment> implements Collection<A> {

  private Collection<A> delegate;
  private AnnotatorConfig config;
  private List<Annotator<? super A>> annotators = new ArrayList<>();
  private ImmutableList<Annotator<? super A>> cachedAnnotators;

  /**
   * Create a collection with the default {@link AnnotatorConfig}
   *
   * @see #AnnotatedCollection(Collection, AnnotatorConfig)
   */
  public AnnotatedCollection(Collection<A> delegate) {
    this(delegate, new AnnotatorConfig());
  }

  /**
   * @param delegate Underlying collection of {@link AnnotatedSegment}s to annotate
   * @param config {@link AnnotatorConfig} to use when annotating this collection
   */
  public AnnotatedCollection(Collection<A> delegate, AnnotatorConfig config) {
    this.delegate = delegate;
    this.config = config;
  }

  /**
   * Record that an {@link Annotator} was used on this collection
   */
  public void addAnnotator(Annotator<? super A> annotator) {
    annotators.add(annotator);

    // Clear cached annotators
    if (cachedAnnotators != null) {
      cachedAnnotators = null;
    }
  }

  /**
   * @return List of all {@link Annotator}s used on this collection
   */
  public ImmutableList<Annotator<? super A>> getAnnotators() {
    // ImmutableList.copyOf will always copy the annotators list, but it is desirable to have the
    // strong return type, so we cache the copy until changes are made.
    if (cachedAnnotators == null) {
      cachedAnnotators = ImmutableList.copyOf(annotators);
    }
    return cachedAnnotators;
  }

  /**
   * @return the {@link AnnotatorConfig} to pass to any {@link Annotator}
   */
  public AnnotatorConfig getConfig() {
    return config;
  }

  /**
   * @param config The {@link AnnotatorConfig} to use when annotating this collection
   */
  public void setConfig(AnnotatorConfig config) {
    this.config = config;
  }

  /**
   * Writes a tab-delimited results file with all annotation columns, for all
   * {@link AnnotatedSegment}s in this collection.
   *
   * @param outputFile Path to desired output
   * @throws IOException
   */
  public void writeResults(String outputFile) throws IOException {
    if (isEmpty()) {
      return;
    }
    // FIXME not guaranteed to match toAnalysisString delim
    final String delim = "\t";
    PrintWriter writer = Files.openAppropriateWriter(outputFile);

    StringJoiner headerJoiner = new StringJoiner(delim);
    // Base segment used to build headers
    // FIXME this is annoying
    A ref = iterator().next();

    // Add base header
    for (String baseColumn : ref.getHeader()) {
      headerJoiner.add(baseColumn);
    }

    // Create annotation header
    for (Annotator<? super A> annotator : annotators) {
      for (Annotation annotation : ref.getAnnotations(annotator)) {
        headerJoiner.add(annotation.getAnnotationName());
      }
    }

    // Write header
    writer.println(headerJoiner.toString());

    // Add row values
    for (A segment : this) {
      StringJoiner rowJoiner = new StringJoiner(delim);
      // Core row values
      rowJoiner.add(segment.toAnalysisString());

      // Annotation row values
      for (Annotator<? super A> annotator : annotators) {
        for (Annotation annotation : segment.getAnnotations(annotator)) {
          rowJoiner.add(annotation.getAnnotationValue());
        }
      }

      writer.println(rowJoiner.toString());
    }

    writer.close();
  }

  // -- Delegation methods --

  /**
   * @return
   * @see java.util.Collection#size()
   */
  @Override
  public int size() {
    return delegate.size();
  }

  /**
   * @return
   * @see java.util.Collection#isEmpty()
   */
  @Override
  public boolean isEmpty() {
    return delegate.isEmpty();
  }

  /**
   * @param o
   * @return
   * @see java.util.Collection#contains(java.lang.Object)
   */
  @Override
  public boolean contains(Object o) {
    return delegate.contains(o);
  }

  /**
   * @return
   * @see java.util.Collection#iterator()
   */
  @Override
  public Iterator<A> iterator() {
    return delegate.iterator();
  }

  /**
   * @return
   * @see java.util.Collection#toArray()
   */
  @Override
  public Object[] toArray() {
    return delegate.toArray();
  }

  /**
   * @param a
   * @return
   * @see java.util.Collection#toArray(java.lang.Object[])
   */
  @Override
  public <T> T[] toArray(T[] a) {
    return delegate.toArray(a);
  }

  /**
   * @param e
   * @return
   * @see java.util.Collection#add(java.lang.Object)
   */
  @Override
  public boolean add(A e) {
    return delegate.add(e);
  }

  /**
   * @param o
   * @return
   * @see java.util.Collection#remove(java.lang.Object)
   */
  @Override
  public boolean remove(Object o) {
    return delegate.remove(o);
  }

  /**
   * @param c
   * @return
   * @see java.util.Collection#containsAll(java.util.Collection)
   */
  @Override
  public boolean containsAll(Collection<?> c) {
    return delegate.containsAll(c);
  }

  /**
   * @param c
   * @return
   * @see java.util.Collection#addAll(java.util.Collection)
   */
  @Override
  public boolean addAll(Collection<? extends A> c) {
    return delegate.addAll(c);
  }

  /**
   * @param c
   * @return
   * @see java.util.Collection#removeAll(java.util.Collection)
   */
  @Override
  public boolean removeAll(Collection<?> c) {
    return delegate.removeAll(c);
  }

  /**
   * @param c
   * @return
   * @see java.util.Collection#retainAll(java.util.Collection)
   */
  @Override
  public boolean retainAll(Collection<?> c) {
    return delegate.retainAll(c);
  }

  /**
   * @see java.util.Collection#clear()
   */
  @Override
  public void clear() {
    delegate.clear();
  }

  /**
   * @param o
   * @return
   * @see java.util.Collection#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    return delegate.equals(o);
  }

  /**
   * @return
   * @see java.util.Collection#hashCode()
   */
  @Override
  public int hashCode() {
    return delegate.hashCode();
  }

}
