package org.genvisis.common;

import java.util.Vector;

import com.google.common.primitives.Bytes;

public class ByteVector extends Vector<Byte> implements PrimitiveVector {

  private static final long serialVersionUID = -7294186396263982757L;

  /**
   * Create an empty vector
   */
  public ByteVector() {
    super();
  }

  /**
   * Create a vector with the given initial values.
   */
  public ByteVector(byte[] initialArray) {
    super(Bytes.asList(initialArray));
  }

  /**
   * Create an empty vector with the specified initial size.
   */
  public ByteVector(int initialSize) {
    this(new byte[initialSize]);
  }
}
