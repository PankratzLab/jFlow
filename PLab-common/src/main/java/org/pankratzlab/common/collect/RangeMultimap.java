package org.pankratzlab.common.collect;

import java.util.Collection;
import java.util.Map;
import java.util.Map.Entry;
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
 * @param <? extends Collection<V>>
 */
public interface RangeMultimap<K extends Comparable<?>, V> extends Multimap<Range<K>, V> {

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
  Collection<V> get(Range<K> key);

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
   * {@link #put(Object, Object)} without the (generally useless) return value based on size
   * changing
   */
  void putBlind(Range<K> range, V value);

  /**
   * In keeping with the contract of the {@link Multimap} interface, this method will never return
   * null, instead returning an empty {@link Collection} when the key is not mapped
   */
  Collection<V> get(K key);

  Map<Range<K>, ? extends Collection<V>> asDescendingMapOfRanges();

  Map<Range<K>, ? extends Collection<V>> asMapOfRanges();

  Entry<Range<K>, ? extends Collection<V>> getEntry(K key);

  void remove(Range<K> range);

  Range<K> span();

  RangeMultimap<K, V> subRangeMap(Range<K> range);

}
