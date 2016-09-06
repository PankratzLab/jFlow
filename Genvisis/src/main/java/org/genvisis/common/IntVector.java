package org.genvisis.common;

import java.util.Vector;

import com.google.common.primitives.Ints;

/**
 * {@link PrimitiveVector} of {@link Integer} types.
 */
public class IntVector extends Vector<Integer> implements PrimitiveVector {

  private static final long serialVersionUID = 6325652813358106316L;

  /**
   * Create an empty vector
   */
  public IntVector() {
    super();
  }

  /**
   * Create an empty vector with the specified initial size.
   */
  public IntVector(int initialSize) {
    super(initialSize);
  }

  /**
   * Create a vector with the given initial values.
   */
  public IntVector(int[] initialArray) {
    super(Ints.asList(initialArray));
  }

  /**
   * Add a given value to this vector if it is not already present.
   * <p>
   * FIXME: should use a set instead of creating this method.
   * </p>
   * <p>
   * FIXME: duplicate logic with {@link DoubleVector} due to primitives being terrible.
   * </p>
   *
   * @param v Value to add if not already present in this vector
   * @return {@code true} if the value was added.
   */
  public boolean addIfAbsent(int v) {
    boolean added = false;
    if (!Ints.contains(Ints.toArray(this), v)) {
      added = add(v);
    }
    return added;
  }

  @Override
  public IntVector clone() {
    return new IntVector(Ints.toArray(this));
  }
}
