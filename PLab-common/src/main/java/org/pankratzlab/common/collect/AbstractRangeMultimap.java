package org.pankratzlab.common.collect;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import com.google.common.base.Functions;
import com.google.common.base.Suppliers;
import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

/**
 * A basic implementation of {@link RangeMultimap}
 *
 * @param <K>
 * @param <V>
 * @param <C>
 */
public abstract class AbstractRangeMultimap<K extends Comparable<?>, V, C extends ImmutableCollection<V>> implements RangeMultimap<K, V> {

  private final Map<Range<K>, C> mapOfRanges;
  private final RangeMap<K, C> rangeMap;

  protected AbstractRangeMultimap(RangeMap<K, C> underlyingRangeMap) {
    super();
    rangeMap = underlyingRangeMap;
    mapOfRanges = getRangeMap().asMapOfRanges();
  }

  @Override
  public Map<Range<K>, C> asDescendingMapOfRanges() {
    return getRangeMap().asDescendingMapOfRanges();
  }

  @Override
  public Map<Range<K>, Collection<V>> asMap() {
    // Necessary to return a Map<Range<K>, Collection<V>>
    return collectionMapViewMapOfRanges();
  }

  @Override
  public Map<Range<K>, C> asMapOfRanges() {
    return mapOfRanges;
  }

  @Override
  public void clear() {
    getRangeMap().clear();
  }

  /**
   * @deprecated This method will only return true when {@link #asMapOfRanges()} contains the exact
   *             entry. Checking if a range is mapped to any value(s) can better be performed using
   *             {@link #subRangeMap(Range)}
   */
  @Deprecated
  @Override
  public boolean containsEntry(Object key, Object value) {
    C mappedValue = mapOfRanges.get(key);
    return value.equals(mappedValue);
  }

  /**
   * @deprecated This method will only return true when {@link #asMapOfRanges()} contains the exact
   *             Range as a Key. Checking if a range is mapped to any value(s) can better be
   *             performed using {@link #subRangeMap(Range)}
   */
  @Deprecated
  @Override
  public boolean containsKey(Object key) {
    return mapOfRanges.containsKey(key);
  }

  @Deprecated
  @Override
  public boolean containsValue(Object value) {
    return mapOfRanges.values().stream().anyMatch(c -> c.contains(value));
  }

  /**
   * Unlike other implementations of {@link Multimap}, this implementation is not stored with
   * non-distinct keys and an easily accesible Collection of values. This method is therefore much
   * less efficient and returns an ImmutableCollection that cannot be used to modify the underlying
   * {@link AbstractRangeMultimap}
   */
  @Override
  public ImmutableCollection<Entry<Range<K>, V>> entries() {
    return mapOfRanges.entrySet().stream()
                      .collect(ImmutableSet.Builder<Entry<Range<K>, V>>::new,
                               (builder,
                                entry) -> entry.getValue().stream()
                                               .map(singleValue -> Maps.immutableEntry(entry.getKey(),
                                                                                       singleValue))
                                               .forEach(builder::add),
                               (builder1, builder2) -> builder1.addAll(builder2.build()))
                      .build();
  }

  @SuppressWarnings("rawtypes")
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof AbstractRangeMultimap)) return false;
    AbstractRangeMultimap other = (AbstractRangeMultimap) obj;
    if (getRangeMap() == null) {
      if (other.getRangeMap() != null) return false;
    } else if (!getRangeMap().equals(other.getRangeMap())) return false;
    return true;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((getRangeMap() == null) ? 0 : getRangeMap().hashCode());
    return result;
  }

  @Override
  public C get(K key) {
    C value = getRangeMap().get(key);
    return value != null ? value : emptyCollection();
  }

  /**
   * @deprecated This method will only return the mapped values when {@link #asMapOfRanges()}
   *             contains the exact Range as a Key. Checking what value(s) a range is mapped to can
   *             better be performed using {@link #subRangeMap(Range)}
   */
  @Override
  @Deprecated
  public C get(Range<K> key) {
    C values = mapOfRanges.get(key);
    return values != null ? values : emptyCollection();
  }

  @Override
  public Entry<Range<K>, C> getEntry(K key) {
    return getRangeMap().getEntry(key);
  }

  @Override
  public boolean isEmpty() {
    if (!mapOfRanges.isEmpty()) {
      Iterator<Entry<Range<K>, C>> entryIterator = mapOfRanges.entrySet().iterator();
      for (Map.Entry<Range<K>, C> entry = entryIterator.next(); entryIterator.hasNext(); entry = entryIterator.next()) {
        if (entry.getValue().isEmpty()) {
          entryIterator.remove();
        } else {
          break;
        }
      }
    }
    return mapOfRanges.isEmpty();
  }

  /**
   * Unlike other implementations of {@link Multimap}, this implementation is not stored with
   * non-distinct keys and an easily accessible Collection of values. As such, this method returns a
   * newly constructed {@link ImmutableMultiset} rather than a view that can be used to modify the
   * underlying {@link Multimap} {@link AbstractRangeMultimap}
   */
  @Override
  public ImmutableMultiset<Range<K>> keys() {
    return mapOfRanges.entrySet().stream()
                      .collect(ImmutableMultiset.toImmutableMultiset(Entry::getKey,
                                                                     e -> e.getValue().size()));
  }

  @Override
  public Set<Range<K>> keySet() {
    return mapOfRanges.keySet();
  }

  @Override
  public void purgeEmptyValues() {
    if (!mapOfRanges.isEmpty()) {
      Iterator<Entry<Range<K>, C>> entryIterator = mapOfRanges.entrySet().iterator();
      for (Map.Entry<Range<K>, C> entry = entryIterator.next(); entryIterator.hasNext(); entry = entryIterator.next())
        if (entry.getValue().isEmpty()) entryIterator.remove();
    }
  }

  @Override
  public void putBlind(Range<K> range, V value) {
    if (!range.isEmpty()) {
      RangeMap<K, C> existingSubRangeMap = getRangeMap().subRangeMap(range);
      Map<Range<K>, C> existingSubRangeMapOfRanges = existingSubRangeMap.asMapOfRanges();
      if (!existingSubRangeMapOfRanges.isEmpty()) {
        RangeMap<K, C> updatedMappings = newRangeMap();
        for (Map.Entry<Range<K>, C> existingEntry : existingSubRangeMapOfRanges.entrySet()) {
          updatedMappings.put(existingEntry.getKey(),
                              addToCollection(existingEntry.getValue(), value));
        }
        existingSubRangeMap.putAll(updatedMappings);
      }
      RangeSet<K> unmappedRangeSet = TreeRangeSet.create();
      unmappedRangeSet.add(range);
      unmappedRangeSet.removeAll(existingSubRangeMapOfRanges.keySet());
      Supplier<C> collectionSupplier = Suppliers.memoize(() -> createCollection(value));
      unmappedRangeSet.asRanges().forEach(r -> getRangeMap().put(r, collectionSupplier.get()));
    }

  }

  /**
   * The return value from this implementation should usually be ignored. Per contract, it returns
   * true when the size of the {@link AbstractRangeMultimap} changed but per the implementation of
   * {@link #size()} this does not necessarily indicate whether the put changed the mappings or not
   * Use {@link #putBlind(Range, Object)} to avoid the before and after size checking
   */
  @Override
  public boolean put(Range<K> key, V value) {
    int oldSize = size();
    putBlind(key, value);
    return oldSize != size();
  }

  /**
   * The return value from this implementation should usually be ignored. It returns true when the
   * size of the {@link AbstractRangeMultimap} changes but per the implementation of {@link #size()}
   * this does not necessarily indicate whether the put changed the mappings or not
   */
  @Override
  public boolean putAll(Multimap<? extends Range<K>, ? extends V> multimap) {
    int oldSize = size();
    for (Entry<? extends Range<K>, ? extends Collection<? extends V>> entry : multimap.asMap()
                                                                                      .entrySet()) {
      putAllBlind(entry.getKey(), entry.getValue());
    }
    return oldSize != size();
  }

  /**
   * The return value from this implementation should usually be ignored. It returns true when the
   * size of the {@link AbstractRangeMultimap} changes but per the implementation of {@link #size()}
   * this does not necessarily indicate whether the put changed the mappings or not
   */
  @Override
  public boolean putAll(Range<K> key, Iterable<? extends V> values) {
    int oldSize = size();
    putAllBlind(key, values);
    return oldSize != size();
  }

  private void putAllBlind(Range<K> key, Iterable<? extends V> values) {
    if (!key.isEmpty()) {
      RangeMap<K, C> existingSubRangeMap = getRangeMap().subRangeMap(key);
      Map<Range<K>, C> existingSubRangeMapOfRanges = existingSubRangeMap.asMapOfRanges();
      if (!existingSubRangeMapOfRanges.isEmpty()) {
        RangeMap<K, C> updatedMappings = newRangeMap();
        for (Map.Entry<Range<K>, C> existingEntry : existingSubRangeMapOfRanges.entrySet()) {
          updatedMappings.put(existingEntry.getKey(),
                              addAllToCollection(existingEntry.getValue(), values));
        }
        existingSubRangeMap.putAll(updatedMappings);
      }
      RangeSet<K> unmappedRangeSet = TreeRangeSet.create();
      unmappedRangeSet.add(key);
      unmappedRangeSet.removeAll(existingSubRangeMapOfRanges.keySet());
      if (!unmappedRangeSet.isEmpty()) {
        C valuesColl = createCollection(values);
        unmappedRangeSet.asRanges().forEach(r -> getRangeMap().put(r, valuesColl));
      }
    }
  }

  @Override
  public boolean remove(Object key, Object value) {
    Range<K> rangeKey = rangeKeyOrNull(key);
    if (rangeKey == null) return false;
    RangeMap<K, C> existingSubRangeMap = getRangeMap().subRangeMap(rangeKey);
    Map<Range<K>, C> existingSubRangeMapOfRanges = existingSubRangeMap.asMapOfRanges();
    if (existingSubRangeMapOfRanges.isEmpty()) return false;
    boolean changed = false;
    RangeMap<K, C> updatedMappings = newRangeMap();
    for (Map.Entry<Range<K>, C> existingEntry : existingSubRangeMapOfRanges.entrySet()) {
      C curValues = existingEntry.getValue();
      if (curValues.contains(value)) {
        changed = true;
        updatedMappings.put(existingEntry.getKey(),
                            removeFromCollection(existingEntry.getValue(), value));
      }
    }
    existingSubRangeMap.putAll(updatedMappings);
    return changed;
  }

  @Override
  public void remove(Range<K> range) {
    getRangeMap().remove(range);
  }

  /**
   * <b>The return value from this method should not be used.</b> While this method will remove all
   * values mapped with the key Range provided, previous mappings will only be returned as the prior
   * result of {@link #get(Range)} In many cases, an empty ImmutableCollection will be returned,
   * this <b>does</b> not indicate that no mappings were removed.
   */
  @Override
  public C removeAll(Object key) {
    Range<K> rangeKey = rangeKeyOrNull(key);
    if (rangeKey == null) return emptyCollection();
    C prevMapping = get(rangeKey);
    getRangeMap().remove(rangeKey);
    return prevMapping;
  }

  /**
   * <b>The return value from this method should not be used.</b> While this method will replace all
   * values mapped with the key Range provided, previous mappings will only be returned as the prior
   * result of {@link #get(Range)} In many cases, an empty ImmutableCollection will be returned,
   * this <b>does</b> not indicate that no mappings were replaced.
   */
  @Override
  public C replaceValues(Range<K> key, Iterable<? extends V> values) {
    C prevValues = get(key);
    getRangeMap().put(key, createCollection(values));
    return prevValues;
  }

  @Override
  public int size() {
    return mapOfRanges.values().parallelStream().collect(Collectors.summingInt(C::size));
  }

  @Override
  public Range<K> span() {
    return getRangeMap().span();
  }

  /**
   * Unlike other implementations of {@link Multimap}, this implementation is not stored with
   * non-distinct keys and an easily accesible Collection of values. This method is therefore much
   * less efficient and returns an ImmutableCollection that cannot be used to modify the underlying
   * {@link AbstractRangeMultimap}
   */
  @Override
  public ImmutableCollection<V> values() {
    return mapOfRanges.values().stream()
                      .collect(ImmutableList.Builder<V>::new,
                               (builder, values) -> builder.addAll(values),
                               (builder1, builder2) -> builder1.addAll(builder2.build()))
                      .build();
  }

  /**
   * @param currentCollection
   * @param additions
   * @return a new C containing the contents of currentCollection and additions
   */
  protected abstract C addAllToCollection(C currentCollection, Iterable<? extends V> additions);

  /**
   * @param currentCollection
   * @param addition
   * @return a new C containing the contents of currentCollection and addition
   */
  protected abstract C addToCollection(C currentCollection, V addition);

  /**
   * @param values
   * @return a new C containing values
   */
  protected abstract C createCollection(Iterable<? extends V> values);

  /**
   * @param value
   * @return a new C containing value
   */
  protected abstract C createCollection(V value);

  /**
   * @return an empty C
   */
  protected abstract C emptyCollection();

  /**
   * @return a new {@link RangeMap} of the correct underlying type
   */
  protected abstract RangeMap<K, C> newRangeMap();

  /**
   * @param underlyingRangeMap
   * @return a new {@link AbstractRangeMultimap} of the correct type with the specified
   *         underlyingRangeMap as its source
   */
  protected abstract RangeMultimap<K, V> newRangeMultiMap(RangeMap<K, C> underlyingRangeMap);

  /**
   * @param currentValues
   * @param removeValue
   * @return a new C containing the contents of currentCollection with removeValue removed
   */
  protected abstract C removeFromCollection(C currentValues, Object removeValue);

  /**
   * @return the underlying rangeMap
   */
  protected RangeMap<K, C> getRangeMap() {
    return rangeMap;
  }

  /**
   * A somewhat hacky way to "transform" the {@code Map<Range<K>, C>} to a
   * {@code Map<Range<K>, Collection<V>>}
   */
  private Map<Range<K>, Collection<V>> collectionMapViewMapOfRanges() {
    return Maps.transformValues(mapOfRanges, Functions.identity());
  }

  @SuppressWarnings("unchecked")
  private Range<K> rangeKeyOrNull(Object key) {
    if (key == null) return null;
    try {
      return (Range<K>) key;
    } catch (ClassCastException cce) {
      return null;
    }
  }

}
