package org.pankratzlab.common.collect;

import java.util.Collection;
import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;

/**
 * An interface for a mapping from Ranges to multiple possible values, implementing both the
 * {@link RangeMap} and {@link RangeMultimap} interfaces
 *
 * @param <K>
 * @param <V>
 * @param <C>
 */
public interface RangeMultimap<K extends Comparable<?>, V, C extends Collection<V>> extends RangeMap<K, C>, Multimap<Range<K>, V> {

  /**
   * Removes any entries with empty collections as their value, such that {@link #asMapOfRanges()}
   * will contain only "legitimate" mappings to non-empty collections
   */
  void purgeEmptyValues();

  /**
   * Returns {@code asMapOfRanges().hashCode()} to match to equal {@link RangeMultimap}s, the same
   * caveats as with {@link #equals(Object)} apply
   */
  @Override
  int hashCode();

  /**
   * {@link RangeMultimap}s are only equal to other {@link RangeMultimap}s where their underlying
   * {@link RangeMap}s are equal. Unless {@link #putCoalescing(Range, ImmutableCollection)} was used
   * for every put in both {@link RangeMultimap}s and {@link #purgeEmptyValues()} is called prior to
   * equals, this method will not return true for two {@link RangeMultimap}s that would return the
   * same result for every possible call to {@link #get(Comparable)}
   */
  @Override
  boolean equals(Object obj);

  /**
   * @deprecated This method will only return the mapped values when {@link #asMapOfRanges()}
   *             contains the exact Range as a Key. Checking what value(s) a range is mapped to can
   *             better be performed using {@link #subRangeMap(Range)}
   */
  @Override
  @Deprecated
  C get(Range<K> key);

  /**
   * @deprecated This method will only return true when {@link #asMapOfRanges()} contains the exact
   *             Range as a Key. Checking if a range is mapped to any value(s) can better be
   *             performed using {@link #subRangeMap(Range)}
   */
  @Deprecated
  @Override
  boolean containsKey(Object key);

  /**
   * @deprecated This method will only return true when {@link #asMapOfRanges()} contains the exact
   *             entry. Checking if a range is mapped to any value(s) can better be performed using
   *             {@link #subRangeMap(Range)}
   */
  @Override
  @Deprecated
  boolean containsEntry(Object key, Object value);

  /**
   * Adds the items in value to the mappings for range. Specifically, after a call to put(range,
   * value), if range.contains(k), then get(k) will return a collection that contains everything in
   * value. To perform the replacing put that {@link RangeMap}s generally exhibit, use
   * {@link #replaceValues(Range, Iterable)} If range is empty, then this is a no-op.
   */
  @Override
  void put(Range<K> range, C value);

  /**
   * In keeping with the contract of the {@link Multimap} interface, this method will never return
   * null, instead returning an empty {@linkImmutableCollection} when the key is not mapped
   */
  @Override
  C get(K key);

}
