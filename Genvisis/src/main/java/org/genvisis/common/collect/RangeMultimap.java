package org.genvisis.common.collect;

import java.util.Collection;
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

}
