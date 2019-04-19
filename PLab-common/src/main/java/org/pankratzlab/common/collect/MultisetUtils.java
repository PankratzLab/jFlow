package org.pankratzlab.common.collect;

import java.util.Comparator;
import java.util.Optional;

import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;

/**
 * A static utility class that adds functionality missing from {@link Multisets}
 */
public class MultisetUtils {

  private MultisetUtils() {}

  /**
   * @return An {@link Optional} containing the {@link Multiset.Entry} in {@code multiset} with the
   *         maximum count
   */
  public static <E> Optional<Multiset.Entry<E>> maxCount(Multiset<E> multiset) {
    return multiset.entrySet().stream().max(Comparator.comparing(Multiset.Entry::getCount));
  }

  /**
   * @return An {@link Optional} containing the {@link Multiset.Entry} in {@code multiset} with the
   *         minimum count
   */
  public static <E> Optional<Multiset.Entry<E>> minCount(Multiset<E> multiset) {
    return multiset.entrySet().stream().min(Comparator.comparing(Multiset.Entry::getCount));
  }

}
