package org.genvisis.common;

import com.google.common.base.Joiner;

public class IterableUtils {

  private IterableUtils() {
    // prevent instantiation of static utility class
  }
  
  /**
   * Returns a String of an {@link Iterable}'s elements, separated by the specified delimiter
   *
   * @param iterable an {@link Iterable} of Objects
   * @param delimiter String delimiter
   * @return String of printed objects
   */
  public static String toStr(Iterable<?> iterable, String delimiter) {
    return Joiner.on(delimiter).join(iterable);
  }

}
