package org.pankratzlab.common.collect;

import java.util.Collection;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

public interface RangeSetMultimap<K extends Comparable<?>, V> extends RangeMultimap<K, V> {

  /**
   * @deprecated This method will only return the mapped values when {@link #asMapOfRanges()}
   *             contains the exact Range as a Key. Checking what value(s) a range is mapped to can
   *             better be performed using {@link #subRangeMap(Range)}
   */
  @Override
  @Deprecated
  Set<V> get(Range<K> key);

  /**
   * In keeping with the contract of the {@link Multimap} interface, this method will never return
   * null, instead returning an empty {@link Collection} when the key is not mapped
   */
  @Override
  Set<V> get(K key);

  @Override
  Map<Range<K>, ? extends Set<V>> asDescendingMapOfRanges();

  @Override
  Map<Range<K>, ? extends Set<V>> asMapOfRanges();

  @Override
  Entry<Range<K>, ? extends Set<V>> getEntry(K key);

  @Override
  RangeSetMultimap<K, V> subRangeMap(Range<K> range);

}
