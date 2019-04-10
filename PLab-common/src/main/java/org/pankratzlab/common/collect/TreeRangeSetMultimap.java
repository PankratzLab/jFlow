package org.pankratzlab.common.collect;

import java.util.Set;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeRangeMap;

/**
 * An implementation of {@link AbstractRangeMultimap} using a {@link TreeRangeMap} and
 * {@link ImmutableSet}
 *
 * @param <K>
 * @param <V>
 */
public class TreeRangeSetMultimap<K extends Comparable<?>, V>
                                 extends AbstractRangeMultimap<K, V, ImmutableSet<V>>
                                 implements RangeSetMultimap<K, V> {

  private TreeRangeSetMultimap() {
    this(TreeRangeMap.create());
  }

  private TreeRangeSetMultimap(RangeMap<K, ImmutableSet<V>> underlyingRangeMap) {
    super(underlyingRangeMap);
  }

  @Override
  protected ImmutableSet<V> addAllToCollection(ImmutableSet<V> currentCollection,
                                               Iterable<? extends V> additions) {
    return new ImmutableSet.Builder<V>().addAll(currentCollection).addAll(additions).build();
  }

  @Override
  protected ImmutableSet<V> addToCollection(ImmutableSet<V> currentCollection, V addition) {
    ImmutableSet.Builder<V> newSetBuilder = ImmutableSet.builderWithExpectedSize(currentCollection.size()
                                                                                 + 1);
    return newSetBuilder.addAll(currentCollection).add(addition).build();
  }

  @Override
  protected ImmutableSet<V> createCollection(Iterable<? extends V> values) {
    return new ImmutableSet.Builder<V>().addAll(values).build();
  }

  @Override
  protected ImmutableSet<V> createCollection(V value) {
    return ImmutableSet.of(value);
  }

  @Override
  protected ImmutableSet<V> emptyCollection() {
    return ImmutableSet.of();
  }

  @Override
  protected RangeMap<K, ImmutableSet<V>> newRangeMap() {
    return TreeRangeMap.create();
  }

  @Override
  protected TreeRangeSetMultimap<K, V> newRangeMultiMap(RangeMap<K, ImmutableSet<V>> underlyingRangeMap) {
    return new TreeRangeSetMultimap<>(underlyingRangeMap);
  }

  @Override
  protected ImmutableSet<V> removeFromCollection(ImmutableSet<V> currentValues,
                                                 Object removeValue) {
    Set<V> newSetItems = Sets.newHashSet(currentValues);
    newSetItems.remove(removeValue);
    return ImmutableSet.copyOf(newSetItems);
  }

  @Override
  public TreeRangeSetMultimap<K, V> subRangeMap(Range<K> range) {
    return newRangeMultiMap(getRangeMap().subRangeMap(range));
  }

  public static <K extends Comparable<?>, V> TreeRangeSetMultimap<K, V> create() {
    return new TreeRangeSetMultimap<>();
  }

}
