package org.genvisis.common.collect;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import com.google.common.base.Functions;
import com.google.common.collect.BoundType;
import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableRangeSet;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterables;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

public abstract class RangeMultimap<K extends Comparable<?>, V, C extends ImmutableCollection<V>> implements RangeMap<K, C>, Multimap<Range<K>, V>

{

  private final Map<Range<K>, C> mapOfRanges;
  private final RangeMap<K, C> rangeMap;

  protected RangeMultimap(RangeMap<K, C> underlyingRangeMap) {
    super();
    rangeMap = underlyingRangeMap;
    mapOfRanges = rangeMap.asMapOfRanges();
  }

  @Override
  public Map<Range<K>, C> asDescendingMapOfRanges() {
    return rangeMap.asDescendingMapOfRanges();
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
    rangeMap.clear();
  }

  /**
   * @deprecated This method will only return true when {@link #asMapOfRanges()} contains the exact
   *             entry. Checking if a range is mapped to any value(s) can better be performed using
   *             {@link #subRangeMap(Range)}
   */
  @Override
  @Deprecated
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

  @Override
  public boolean containsValue(Object value) {
    return mapOfRanges.values().stream().anyMatch(c -> c.contains(value));
  }

  /**
   * Unlike other implementations of {@link Multimap}, this implementation is not stored with
   * non-distinct keys and an easily accesible Collection of values. This method is therefore much
   * less efficient and returns an ImmutableCollection that cannot be used to modify the underlying
   * {@link RangeMultimap}
   */
  @Override
  public ImmutableCollection<Entry<Range<K>, V>> entries() {
    return mapOfRanges.entrySet().stream()
                      .collect(ImmutableSet.Builder<Entry<Range<K>, V>>::new,
                               (builder,
                                entry) -> builder.addAll(entry.getValue().stream()
                                                              .collect(HashSet::new,
                                                                       (set,
                                                                        value) -> set.add(Maps.immutableEntry(entry.getKey(),
                                                                                                              value)),
                                                                       HashSet::addAll)),
                               (builder1, builder2) -> builder1.addAll(builder2.build()))
                      .build();
  }

  /**
   * {@link RangeMultimap}s are only equal to other {@link RangeMultimap}s where their underlying
   * {@link RangeMap}s are equal.
   */
  @SuppressWarnings("rawtypes")
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof RangeMultimap)) return false;
    RangeMultimap other = (RangeMultimap) obj;
    if (rangeMap == null) {
      if (other.rangeMap != null) return false;
    } else if (!rangeMap.equals(other.rangeMap)) return false;
    return true;
  }

  /**
   * Uses {@link RangeMap#hashCode()} to match to equal {@link RangeMultimap}s
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((rangeMap == null) ? 0 : rangeMap.hashCode());
    return result;
  }

  /**
   * In keeping with the contract of the {@link Multimap} interface, this method will never return
   * null, instead returning an empty {@linkImmutableCollection} when the key is not mapped
   */
  @Override
  public C get(K key) {
    C value = rangeMap.get(key);
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
    return rangeMap.getEntry(key);
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

  @Override
  public ImmutableMultiset<Range<K>> keys() {
    // TODO doc inefficiency and non-updateyness
    return mapOfRanges.entrySet().stream()
                      .collect(ImmutableMultiset.toImmutableMultiset(Entry::getKey,
                                                                     e -> e.getValue().size()));
  }

  @Override
  public Set<Range<K>> keySet() {
    return mapOfRanges.keySet();
  }

  /**
   * Removes any entries with empty collections as their value, such that
   * {@link RangeMultimap#size()} and {@link Map#size()} from {@link #asMapOfRanges()} will return
   * the size of legitimate mappings
   */
  public void purgeEmptyValues() {
    if (!mapOfRanges.isEmpty()) {
      Iterator<Entry<Range<K>, C>> entryIterator = mapOfRanges.entrySet().iterator();
      for (Map.Entry<Range<K>, C> entry = entryIterator.next(); entryIterator.hasNext(); entry = entryIterator.next())
        if (entry.getValue().isEmpty()) entryIterator.remove();
    }
  }

  /**
   * Adds the items in value to the mappings for range. Specifically, after a call to put(range,
   * value), if range.contains(k), then get(k) will return a collection that contains everything in
   * value. To perform the replacing put that {@link RangeMap}s generally exhibit, use
   * {@link #replaceValues(Range, Iterable)} If range is empty, then this is a no-op.
   */
  @Override
  public void put(Range<K> range, C value) {
    RangeMap<K, C> existingSubRangeMap = rangeMap.subRangeMap(range);
    Map<Range<K>, C> existingSubRangeMapOfRanges = existingSubRangeMap.asMapOfRanges();
    if (!existingSubRangeMapOfRanges.isEmpty()) {
      RangeMap<K, C> updatedMappings = newRangeMap();
      for (Map.Entry<Range<K>, C> existingEntry : existingSubRangeMapOfRanges.entrySet()) {
        updatedMappings.put(existingEntry.getKey(),
                            addAllToCollection(existingEntry.getValue(), value));
      }
      existingSubRangeMap.putAll(updatedMappings);
    }
    RangeSet<K> unmappedRangeSet = TreeRangeSet.create();
    unmappedRangeSet.add(range);
    unmappedRangeSet.removeAll(existingSubRangeMapOfRanges.keySet());
    unmappedRangeSet.asRanges().forEach(r -> rangeMap.put(r, value));
  }

  /**
   * The return value from this implementation should usually be ignored. Per contract, it returns
   * true when the size of the {@link RangeMultimap} changed but per the implementation of
   * {@link #size()} this does not necessarily indicate whether the put changed the mappings or not
   */
  @Override
  public boolean put(Range<K> key, V value) {
    int oldSize = size();
    put(key, createCollection(value));
    return oldSize != size();
  }

  /**
   * The return value from this implementation should usually be ignored. It returns true when the
   * size of the {@link RangeMultimap} changes but per the implementation of {@link #size()} this
   * does not necessarily indicate whether the put changed the mappings or not
   */
  @Override
  public boolean putAll(Multimap<? extends Range<K>, ? extends V> multimap) {
    boolean changed = false;
    for (Entry<? extends Range<K>, ? extends Collection<? extends V>> entry : multimap.asMap()
                                                                                      .entrySet()) {
      if (putAll(entry.getKey(), entry.getValue())) {
        changed = true;
      }
    }
    return changed;
  }

  /**
   * The return value from this implementation should usually be ignored. It returns true when the
   * size of the {@link RangeMultimap} changes but per the implementation of {@link #size()} this
   * does not necessarily indicate whether the put changed the mappings or not
   */
  @Override
  public boolean putAll(Range<K> key, Iterable<? extends V> values) {
    int oldSize = size();
    put(key, createCollection(values));
    return oldSize != size();
  }

  @Override
  public void putAll(RangeMap<K, C> rangeMap) {
    rangeMap.asMapOfRanges().entrySet().forEach(e -> put(e.getKey(), e.getValue()));
  }

  @Override
  public void putCoalescing(Range<K> range, C value) {
    if (mapOfRanges.isEmpty()) {
      put(range, value);
      return;
    }
    putCoalesceEnclosed(range, value);
    Range<K> prevRange;
    Range<K> newRange = range;
    do {
      prevRange = newRange;
      newRange = putCoalescedEdges(prevRange, value);
    } while (prevRange != newRange && !newRange.isEmpty());
    if (!newRange.isEmpty()) {
      put(newRange, value);
    }
  }

  @Override
  public boolean remove(Object key, Object value) {
    Range<K> rangeKey = rangeKeyOrNull(key);
    if (rangeKey == null) return false;
    RangeMap<K, C> existingSubRangeMap = rangeMap.subRangeMap(rangeKey);
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
    rangeMap.remove(range);
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
    rangeMap.remove(rangeKey);
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
    rangeMap.put(key, createCollection(values));
    return prevValues;
  }

  /**
   * @deprecated This method simply returns the size of {@link #asMapOfRanges()} and does not follow
   *             the contract of {@link Multimap}. Additionally, without calling
   *             {@link #purgeEmptyValues()}, the size may contain entries with empty collections as
   *             their value
   */
  @Override
  @Deprecated
  public int size() {
    return mapOfRanges.size();
  }

  @Override
  public Range<K> span() {
    return rangeMap.span();
  }

  @Override
  public RangeMultimap<K, V, C> subRangeMap(Range<K> range) {
    return newRangeMultiMap(rangeMap.subRangeMap(range));
  }

  /**
   * Unlike other implementations of {@link Multimap}, this implementation is not stored with
   * non-distinct keys and an easily accesible Collection of values. This method is therefore much
   * less efficient and returns an ImmutableCollection that cannot be used to modify the underlying
   * {@link RangeMultimap}
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
  protected abstract C addAllToCollection(C currentCollection, Iterable<V> additions);

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
   * @return a new {@link RangeMultimap} of the correct type with the specified underlyingRangeMap
   *         as its source
   */
  protected abstract RangeMultimap<K, V, C> newRangeMultiMap(RangeMap<K, C> underlyingRangeMap);

  /**
   * @param currentValues
   * @param removeValue
   * @return a new C containing the contents of currentCollection with removeValue removed
   */
  protected abstract C removeFromCollection(C currentValues, Object removeValue);

  /**
   * Coalesces the range enclosed by originalRange and puts value when possible
   */
  private void putCoalesceEnclosed(Range<K> originalRange, C value) {
    RangeMap<K, C> targetSubRangeMap = subRangeMap(originalRange);
    Set<C> distinctSubRangeVals = ImmutableSet.copyOf(targetSubRangeMap.asMapOfRanges().values());
    if (distinctSubRangeVals.size() == 1) {
      C currentVal = Iterables.getOnlyElement(distinctSubRangeVals);
      C newVal = addAllToCollection(currentVal, value);
      if ((currentVal.equals(newVal) && newVal.equals(value))
          || perfectCoalescedMatch(targetSubRangeMap, originalRange)) {
        rangeMap.put(originalRange, newVal);
      }
    }
  }

  /**
   * A somewhat hacky way to "transform" the {@code Map<Range<K>, C>} to a
   * {@code Map<Range<K>, Collection<V>>}
   */
  private Map<Range<K>, Collection<V>> collectionMapViewMapOfRanges() {
    return Maps.transformValues(mapOfRanges, Functions.identity());
  }

  private Range<K> emptyRange(K bound) {
    return Range.openClosed(bound, bound);
  }

  /**
   * Locates coalescable ranges at the edges of the originalRange and puts the value to the
   * coalesced ranges
   * 
   * @param originalRange
   * @param value
   * @return the remaining range to put (maybe empty)
   */
  private Range<K> putCoalescedEdges(Range<K> originalRange, C value) {
    Range<K> range = originalRange.intersection(Range.all());
    RangeMap<K, C> targetSubRangeMap = subRangeMap(range);
    Range<K> lowerEdgeRange = null;
    if (range.hasLowerBound()) {
      RangeMap<K, C> lowerSubRangeMap = subRangeMap(Range.upTo(range.lowerEndpoint(),
                                                               oppositeBound(range.lowerBoundType())));

      Range<K> leastTarget = Iterables.getFirst(targetSubRangeMap.asMapOfRanges().keySet(), null);
      Range<K> greatestLower = Iterables.getFirst(lowerSubRangeMap.asDescendingMapOfRanges()
                                                                  .keySet(),
                                                  null);
      C leastVal = (leastTarget != null && leastTarget.lowerEndpoint().equals(range.lowerEndpoint())
                    && leastTarget.lowerBoundType().equals(range.lowerBoundType()))
                                                                                    ? targetSubRangeMap.asMapOfRanges()
                                                                                                       .get(leastTarget)
                                                                                    : emptyCollection();
      C greatestLowerVal = (greatestLower != null
                            && greatestLower.isConnected(range)) ? lowerSubRangeMap.asMapOfRanges()
                                                                                   .get(greatestLower)
                                                                 : emptyCollection();
      C newVal = addAllToCollection(leastVal, value);
      if (greatestLowerVal.equals(newVal)) {
        if (greatestLower == null) greatestLower = Range.upTo(range.lowerEndpoint(),
                                                              oppositeBound(range.lowerBoundType()));
        else if (!greatestLower.isConnected(range)) greatestLower = Range.range(greatestLower.upperEndpoint(),
                                                                                oppositeBound(greatestLower.upperBoundType()),
                                                                                range.lowerEndpoint(),
                                                                                oppositeBound(range.lowerBoundType()));
        if (leastTarget == null) leastTarget = range;
        else if (!leastTarget.isConnected(greatestLower)) leastTarget = Range.range(range.lowerEndpoint(),
                                                                                    range.lowerBoundType(),
                                                                                    leastTarget.lowerEndpoint(),
                                                                                    oppositeBound(leastTarget.lowerBoundType()));
        if (!greatestLower.intersection(leastTarget).isEmpty()) {
          throw new IllegalStateException("Poorly formed coalescing ranges");
        }
        lowerEdgeRange = leastTarget.span(greatestLower);
        rangeMap.put(lowerEdgeRange, newVal);
        if (lowerEdgeRange.hasUpperBound()) {
          range = range.span(lowerEdgeRange);
        } else {
          return emptyRange(range.lowerEndpoint());
        }
      }
    }
    Range<K> upperEdgeRange = null;
    if (range.hasUpperBound()) {

      RangeMap<K, C> upperSubRangeMap = subRangeMap(Range.downTo(range.upperEndpoint(),
                                                                 oppositeBound(range.upperBoundType())));
      Range<K> greatestTarget = Iterables.getFirst(targetSubRangeMap.asDescendingMapOfRanges()
                                                                    .keySet(),
                                                   null);
      Range<K> leastUpper = Iterables.getFirst(upperSubRangeMap.asMapOfRanges().keySet(), null);

      C greatestVal = (greatestTarget != null
                       && greatestTarget.upperEndpoint().equals(range.upperEndpoint())
                       && greatestTarget.upperBoundType().equals(range.upperBoundType()))
                                                                                          ? targetSubRangeMap.asMapOfRanges()
                                                                                                             .get(greatestTarget)
                                                                                          : emptyCollection();
      C leastUpperVal = (leastUpper != null
                         && leastUpper.isConnected(range)) ? upperSubRangeMap.asMapOfRanges()
                                                                             .get(leastUpper)
                                                           : emptyCollection();
      C newVal = addAllToCollection(greatestVal, value);
      if (leastUpperVal.equals(newVal)) {

        if (leastUpper == null) leastUpper = Range.downTo(range.upperEndpoint(),
                                                          oppositeBound(range.upperBoundType()));
        else if (!leastUpper.isConnected(range)) leastUpper = Range.range(leastUpper.lowerEndpoint(),
                                                                          oppositeBound(leastUpper.lowerBoundType()),
                                                                          range.upperEndpoint(),
                                                                          oppositeBound(range.upperBoundType()));
        if (greatestTarget == null) greatestTarget = range;
        else if (!greatestTarget.isConnected(leastUpper)) greatestTarget = Range.range(greatestTarget.upperEndpoint(),
                                                                                       oppositeBound(greatestTarget.upperBoundType()),
                                                                                       range.upperEndpoint(),
                                                                                       range.upperBoundType());
        if (!leastUpper.intersection(greatestTarget).isEmpty()) {
          throw new IllegalStateException("Poorly formed coalescing ranges");
        }
        upperEdgeRange = greatestTarget.span(leastUpper);
        rangeMap.put(upperEdgeRange, newVal);
      }

    }
    if (upperEdgeRange != null && lowerEdgeRange != null) {
      if (upperEdgeRange.isConnected(lowerEdgeRange)) {
        rangeMap.put(upperEdgeRange.span(lowerEdgeRange), mapOfRanges.get(upperEdgeRange));
        return emptyRange(originalRange.lowerEndpoint());
      }
      return Range.range(lowerEdgeRange.upperEndpoint(),
                         oppositeBound(lowerEdgeRange.upperBoundType()),
                         upperEdgeRange.lowerEndpoint(),
                         oppositeBound(upperEdgeRange.lowerBoundType()));
    }
    if (upperEdgeRange != null) return originalRange.intersection(Range.upTo(upperEdgeRange.lowerEndpoint(),
                                                                             oppositeBound(upperEdgeRange.lowerBoundType())));
    if (lowerEdgeRange != null) return originalRange.intersection(Range.downTo(lowerEdgeRange.upperEndpoint(),
                                                                               oppositeBound(lowerEdgeRange.upperBoundType())));
    return originalRange;
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

  private static BoundType oppositeBound(BoundType type) {
    switch (type) {
      case OPEN:
        return BoundType.CLOSED;
      case CLOSED:
        return BoundType.OPEN;
      default:
        throw new IllegalStateException("Unknown " + BoundType.class.getName() + ": " + type);
    }
  }

  private static <K extends Comparable<?>, V> boolean perfectCoalescedMatch(RangeMap<K, V> rangeMap,
                                                                            Range<K> range) {
    RangeSet<K> rangeMapRanges = ImmutableRangeSet.copyOf(rangeMap.asMapOfRanges().keySet());
    return rangeMapRanges.asRanges().size() == 1 && rangeMapRanges.encloses(range);

  }
}
