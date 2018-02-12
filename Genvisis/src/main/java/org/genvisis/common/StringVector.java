package org.genvisis.common;

import java.util.Vector;

/**
 * {@link PrimitiveVector} of {@link String} types.
 */
public class StringVector extends Vector<String> implements PrimitiveVector {

  private static final long serialVersionUID = -3741807344249792544L;

  /**
   * Create an empty vector.
   */
  public StringVector() {
    super();
  }

  /**
   * Create a vector with the given initial values.
   */
  public StringVector(String... strings) {
    for (String s : strings) {
      add(s);
    }
  }
}
