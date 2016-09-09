package org.genvisis.common;

import com.google.common.base.Joiner;

public class IterableUtils {

  private IterableUtils() {
    // prevent instantiation of static utility class
  }
  
  /**
   * Prints an {@link Iterable}'s element's, separated by the specified delimiter
   *
   * @param array an iterable of objects
   * @param delimiter String delimiter
   * @return String of printed objects
   */
  public static String toStr(Iterable<?> iterable, String delimiter) {
    return Joiner.on(delimiter).join(iterable);
  }

}
