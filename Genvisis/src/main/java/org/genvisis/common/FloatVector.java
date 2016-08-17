package org.genvisis.common;

import java.util.Vector;

import com.google.common.primitives.Floats;

/**
 * {@link PrimitiveVector} of {@link Float} types.
 */
public class FloatVector extends Vector<Float> implements PrimitiveVector {

  private static final long serialVersionUID = -4366233294864349627L;

  /**
   * Create an empty vector
   */
  public FloatVector() {
    super();
  }

  /**
   * Create a vector with the given initial values.
   */
  public FloatVector(float[] initialArray) {
    super(Floats.asList(initialArray));
  }

  /**
   * Create an empty vector with the specified initial size.
   */
  public FloatVector(int initialSize) {
    this(new float[initialSize]);
  }
}
