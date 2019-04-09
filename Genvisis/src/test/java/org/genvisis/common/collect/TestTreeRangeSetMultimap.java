package org.genvisis.common.collect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import org.genvisis.ExhaustiveUnitTests;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.pankratzlab.common.collect.RangeMultimap;
import org.pankratzlab.common.collect.TreeRangeSetMultimap;
import com.google.common.collect.BoundType;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

@Category(ExhaustiveUnitTests.class)
public class TestTreeRangeSetMultimap {

  public TestTreeRangeSetMultimap() {}

  @Test
  public void testAppendToArray() {}

  @Test
  public void basicTests() {
    TreeRangeSetMultimap<Integer, String> rangeMultimap = TreeRangeSetMultimap.create();
    rangeMultimap.put(Range.closed(1, 10), "foo");
    assertEquals(ImmutableSet.of("foo"), rangeMultimap.get(3));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(100));
    rangeMultimap.put(Range.closed(3, 7), "bar");
    assertEquals(ImmutableSet.of("foo"), rangeMultimap.get(1));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(100));
    assertEquals(ImmutableSet.of("foo", "bar"), rangeMultimap.get(5));
    rangeMultimap.remove(Range.closed(8, 10));
    assertEquals(ImmutableSet.of("foo"), rangeMultimap.get(1));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(8));
    assertEquals(ImmutableSet.of("foo", "bar"), rangeMultimap.get(5));
    rangeMultimap.remove(Range.closed(2, 5), "foo");
    assertEquals(ImmutableSet.of("foo"), rangeMultimap.get(1));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(2));
    assertEquals(ImmutableSet.of("bar"), rangeMultimap.get(3));
    assertEquals(ImmutableSet.of("foo", "bar"), rangeMultimap.get(6));
    rangeMultimap.replaceValues(Range.closed(0, 100),
                                ImmutableSet.of("junk", "more junk", "stuff"));
    assertEquals(ImmutableSet.of("junk", "more junk", "stuff"), rangeMultimap.get(1));
    assertEquals(ImmutableSet.of("junk", "more junk", "stuff"), rangeMultimap.get(2));
    assertEquals(ImmutableSet.of("junk", "more junk", "stuff"), rangeMultimap.get(3));
    assertEquals(ImmutableSet.of("junk", "more junk", "stuff"), rangeMultimap.get(6));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(101));
    rangeMultimap.clear();
    assertEquals(ImmutableSet.of(), rangeMultimap.get(1));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(2));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(3));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(6));
    assertEquals(ImmutableSet.of(), rangeMultimap.get(101));
    rangeMultimap.put(Range.closed(1, 10), "foo");
    rangeMultimap.put(Range.closed(3, 7), "bar");
    rangeMultimap.put(Range.closed(9, 15), "baz");
  }

  /**
   * Adapted from {@link HashMultimap}
   */
  @Test
  public void testCreate() {
    TreeRangeSetMultimap<String, Integer> multimap = TreeRangeSetMultimap.create();
    multimap.put(Range.singleton("foo"), 1);
    multimap.put(Range.singleton("bar"), 2);
    multimap.put(Range.singleton("foo"), 3);
    assertEquals(ImmutableSet.of(1, 3), multimap.get(Range.singleton("foo")));
  }

  /*
   * All tests below are adapted from {@link TreeRangeMultimap}
   */
  private static final ImmutableList<Range<Integer>> RANGES;
  private static final int MIN_BOUND = -2;
  private static final int MAX_BOUND = 2;

  static {
    ImmutableList.Builder<Range<Integer>> builder = ImmutableList.builder();

    builder.add(Range.<Integer>all());

    // Add one-ended ranges
    for (int i = MIN_BOUND; i <= MAX_BOUND; i++) {
      for (BoundType type : BoundType.values()) {
        builder.add(Range.upTo(i, type));
        builder.add(Range.downTo(i, type));
      }
    }

    // Add two-ended ranges
    for (int i = MIN_BOUND; i <= MAX_BOUND; i++) {
      for (int j = i; j <= MAX_BOUND; j++) {
        for (BoundType lowerType : BoundType.values()) {
          for (BoundType upperType : BoundType.values()) {
            if (i == j & lowerType == BoundType.OPEN & upperType == BoundType.OPEN) {
              continue;
            }
            builder.add(Range.range(i, lowerType, j, upperType));
          }
        }
      }
    }
    RANGES = builder.build();
  }

  @Test
  public void testSpanSingleRange() {
    for (Range<Integer> range : RANGES) {
      RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
      RangeMultimap.put(range, 1);

      try {
        assertEquals(range, RangeMultimap.span());
        assertFalse(range.isEmpty());
      } catch (NoSuchElementException e) {
        assertTrue(range.isEmpty());
      }
    }
  }

  @Test
  public void testSpanTwoRanges() {
    for (Range<Integer> range1 : RANGES) {
      for (Range<Integer> range2 : RANGES) {
        RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
        RangeMultimap.put(range1, 1);
        RangeMultimap.put(range2, 2);

        Range<Integer> expected;
        if (range1.isEmpty()) {
          if (range2.isEmpty()) {
            expected = null;
          } else {
            expected = range2;
          }
        } else {
          if (range2.isEmpty()) {
            expected = range1;
          } else {
            expected = range1.span(range2);
          }
        }

        try {
          assertEquals(expected, RangeMultimap.span());
          assertNotNull(expected);
        } catch (NoSuchElementException e) {
          assertNull(expected);
        }
      }
    }
  }

  @Test
  public void testAllRangesAlone() {
    for (Range<Integer> range : RANGES) {
      Multimap<Integer, Integer> model = HashMultimap.create();
      putModel(model, range, 1);
      RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test = TreeRangeSetMultimap.create();
      test.put(range, 1);
      verify(model, test);
    }
  }

  @Test
  public void testAllRangePairs() {
    for (Range<Integer> range1 : RANGES) {
      for (Range<Integer> range2 : RANGES) {
        Multimap<Integer, Integer> model = HashMultimap.create();
        putModel(model, range1, 1);
        putModel(model, range2, 2);
        RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test = TreeRangeSetMultimap.create();
        test.put(range1, 1);
        test.put(range2, 2);
        verify(model, test);
      }
    }
  }

  @Test
  public void testAllRangeTriples() {
    for (Range<Integer> range1 : RANGES) {
      for (Range<Integer> range2 : RANGES) {
        for (Range<Integer> range3 : RANGES) {
          Multimap<Integer, Integer> model = HashMultimap.create();
          putModel(model, range1, 1);
          putModel(model, range2, 2);
          putModel(model, range3, 3);
          RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test = TreeRangeSetMultimap.create();
          test.put(range1, 1);
          test.put(range2, 2);
          test.put(range3, 3);
          verify(model, test);
        }
      }
    }
  }

  @Test
  public void testPutAll() {
    for (Range<Integer> range1 : RANGES) {
      for (Range<Integer> range2 : RANGES) {
        for (Range<Integer> range3 : RANGES) {
          Multimap<Integer, Integer> model = HashMultimap.create();
          putModel(model, range1, 1);
          putModel(model, range2, 2);
          putModel(model, range3, 3);
          RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test = TreeRangeSetMultimap.create();
          RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test2 = TreeRangeSetMultimap.create();
          // put range2 and range3 into test2, and then put test2 into test
          test.put(range1, 1);
          test2.put(range2, 2);
          test2.put(range3, 3);
          test.putAll(test2);
          verify(model, test);
        }
      }
    }
  }

  @Test
  public void testPutAndRemove() {
    for (Range<Integer> rangeToPut : RANGES) {
      for (Range<Integer> rangeToRemove : RANGES) {
        Multimap<Integer, Integer> model = HashMultimap.create();
        putModel(model, rangeToPut, 1);
        removeModel(model, rangeToRemove);
        RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test = TreeRangeSetMultimap.create();
        test.put(rangeToPut, 1);
        test.remove(rangeToRemove);
        verify(model, test);
      }
    }
  }

  @Test
  public void testPutTwoAndRemove() {
    for (Range<Integer> rangeToPut1 : RANGES) {
      for (Range<Integer> rangeToPut2 : RANGES) {
        for (Range<Integer> rangeToRemove : RANGES) {
          Multimap<Integer, Integer> model = HashMultimap.create();
          putModel(model, rangeToPut1, 1);
          putModel(model, rangeToPut2, 2);
          removeModel(model, rangeToRemove);
          RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test = TreeRangeSetMultimap.create();
          test.put(rangeToPut1, 1);
          test.put(rangeToPut2, 2);
          test.remove(rangeToRemove);
          verify(model, test);
        }
      }
    }
  }

  @Test
  public void testsubRangeMapExhaustive() {
    for (Range<Integer> range1 : RANGES) {
      for (Range<Integer> range2 : RANGES) {
        RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
        RangeMultimap.put(range1, 1);
        RangeMultimap.put(range2, 2);

        for (Range<Integer> subRange : RANGES) {
          RangeMultimap<Integer, Integer, ImmutableSet<Integer>> expected = TreeRangeSetMultimap.create();
          for (Entry<Range<Integer>, ImmutableSet<Integer>> entry : RangeMultimap.asMapOfRanges()
                                                                                 .entrySet()) {
            if (entry.getKey().isConnected(subRange)) {
              expected.putAll(entry.getKey().intersection(subRange), entry.getValue());
            }
          }
          RangeMultimap<Integer, Integer, ImmutableSet<Integer>> subRangeMap = RangeMultimap.subRangeMap(subRange);
          assertEquals(expected, subRangeMap);
          assertEquals(expected.asMapOfRanges(), subRangeMap.asMapOfRanges());
          assertEquals(expected.asDescendingMapOfRanges(), subRangeMap.asDescendingMapOfRanges());
          assertEquals(ImmutableList.copyOf(subRangeMap.asMapOfRanges().entrySet()).reverse(),
                       ImmutableList.copyOf(subRangeMap.asDescendingMapOfRanges().entrySet()));

          if (!expected.asMapOfRanges().isEmpty()) {
            assertEquals(expected.span(), subRangeMap.span());
          }

          for (int i = MIN_BOUND; i <= MAX_BOUND; i++) {
            assertEquals(expected.get(i), subRangeMap.get(i));
          }

          for (Range<Integer> query : RANGES) {
            assertEquals(expected.asMapOfRanges().get(query),
                         subRangeMap.asMapOfRanges().get(query));
          }
        }
      }
    }
  }

  @Test
  public void testSubsubRangeMap() {
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
    RangeMultimap.put(Range.open(3, 7), 1);
    RangeMultimap.put(Range.closed(9, 10), 2);
    RangeMultimap.put(Range.closed(12, 16), 3);
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> sub1 = RangeMultimap.subRangeMap(Range.closed(5,
                                                                                                         11));
    assertEquals(ImmutableMap.of(Range.closedOpen(5, 7), ImmutableSet.of(1), Range.closed(9, 10),
                                 ImmutableSet.of(2)),
                 sub1.asMapOfRanges());
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> sub2 = sub1.subRangeMap(Range.open(6,
                                                                                              15));
    assertEquals(ImmutableMap.of(Range.open(6, 7), ImmutableSet.of(1), Range.closed(9, 10),
                                 ImmutableSet.of(2)),
                 sub2.asMapOfRanges());
  }

  @Test
  public void testsubRangeMapPut() {
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
    RangeMultimap.put(Range.open(3, 7), 1);
    RangeMultimap.put(Range.closed(9, 10), 2);
    RangeMultimap.put(Range.closed(12, 16), 3);
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> sub = RangeMultimap.subRangeMap(Range.closed(5,
                                                                                                        11));
    assertEquals(ImmutableMap.of(Range.closedOpen(5, 7), ImmutableSet.of(1), Range.closed(9, 10),
                                 ImmutableSet.of(2)),
                 sub.asMapOfRanges());
    sub.put(Range.closed(7, 9), 4);
    assertEquals(ImmutableMap.of(Range.closedOpen(5, 7), ImmutableSet.of(1), Range.closedOpen(7, 9),
                                 ImmutableSet.of(4), Range.singleton(9), ImmutableSet.of(4, 2),
                                 Range.openClosed(9, 10), ImmutableSet.of(2)),
                 sub.asMapOfRanges());
    assertEquals(ImmutableMap.of(Range.open(3, 7), ImmutableSet.of(1), Range.closedOpen(7, 9),
                                 ImmutableSet.of(4), Range.singleton(9), ImmutableSet.of(4, 2),
                                 Range.openClosed(9, 10), ImmutableSet.of(2), Range.closed(12, 16),
                                 ImmutableSet.of(3)),
                 RangeMultimap.asMapOfRanges());

    sub = sub.subRangeMap(Range.closedOpen(5, 5));
    sub.put(Range.closedOpen(5, 5), 6); // should be a no-op
    assertEquals(ImmutableMap.of(Range.open(3, 7), ImmutableSet.of(1), Range.closedOpen(7, 9),
                                 ImmutableSet.of(4), Range.singleton(9), ImmutableSet.of(4, 2),
                                 Range.openClosed(9, 10), ImmutableSet.of(2), Range.closed(12, 16),
                                 ImmutableSet.of(3)),
                 RangeMultimap.asMapOfRanges());
    try {
      sub.put(Range.open(9, 12), 5);
      fail("Expected IllegalArgumentException");
    } catch (IllegalArgumentException expected) {}
  }

  @Test
  public void testsubRangeMapRemove() {
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
    RangeMultimap.put(Range.open(3, 7), 1);
    RangeMultimap.put(Range.closed(9, 10), 2);
    RangeMultimap.put(Range.closed(12, 16), 3);
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> sub = RangeMultimap.subRangeMap(Range.closed(5,
                                                                                                        11));
    assertEquals(ImmutableMap.of(Range.closedOpen(5, 7), ImmutableSet.of(1), Range.closed(9, 10),
                                 ImmutableSet.of(2)),
                 sub.asMapOfRanges());
    sub.remove(Range.closed(7, 9));
    assertEquals(ImmutableMap.of(Range.closedOpen(5, 7), ImmutableSet.of(1),
                                 Range.openClosed(9, 10), ImmutableSet.of(2)),
                 sub.asMapOfRanges());
    assertEquals(ImmutableMap.of(Range.open(3, 7), ImmutableSet.of(1), Range.openClosed(9, 10),
                                 ImmutableSet.of(2), Range.closed(12, 16), ImmutableSet.of(3)),
                 RangeMultimap.asMapOfRanges());

    sub.remove(Range.closed(3, 9));
    assertEquals(ImmutableMap.of(Range.openClosed(9, 10), ImmutableSet.of(2)), sub.asMapOfRanges());
    assertEquals(ImmutableMap.of(Range.open(3, 5), ImmutableSet.of(1), Range.openClosed(9, 10),
                                 ImmutableSet.of(2), Range.closed(12, 16), ImmutableSet.of(3)),
                 RangeMultimap.asMapOfRanges());
  }

  @Test
  public void testsubRangeMapClear() {
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> RangeMultimap = TreeRangeSetMultimap.create();
    RangeMultimap.put(Range.open(3, 7), 1);
    RangeMultimap.put(Range.closed(9, 10), 2);
    RangeMultimap.put(Range.closed(12, 16), 3);
    RangeMultimap<Integer, Integer, ImmutableSet<Integer>> sub = RangeMultimap.subRangeMap(Range.closed(5,
                                                                                                        11));
    sub.clear();
    assertEquals(ImmutableMap.of(Range.open(3, 5), ImmutableSet.of(1), Range.closed(12, 16),
                                 ImmutableSet.of(3)),
                 RangeMultimap.asMapOfRanges());
  }

  private void verify(Multimap<Integer, Integer> model,
                      RangeMultimap<Integer, Integer, ImmutableSet<Integer>> test) {
    for (int i = MIN_BOUND - 1; i <= MAX_BOUND + 1; i++) {
      assertEquals(model.get(i), test.get(i));

      Entry<Range<Integer>, ImmutableSet<Integer>> entry = test.getEntry(i);
      assertEquals(model.containsKey(i), entry != null);
      if (entry != null) {
        assertTrue(test.asMapOfRanges().entrySet().contains(entry));
      }
    }
    for (Range<Integer> range : test.asMapOfRanges().keySet()) {
      assertFalse(range.isEmpty());
    }
  }

  private static void putModel(Multimap<Integer, Integer> model, Range<Integer> range,
                               Integer value) {
    for (int i = MIN_BOUND - 1; i <= MAX_BOUND + 1; i++) {
      if (range.contains(i)) {
        model.put(i, value);
      }
    }
  }

  private static void removeModel(Multimap<Integer, Integer> model, Range<Integer> range) {
    for (int i = MIN_BOUND - 1; i <= MAX_BOUND + 1; i++) {
      if (range.contains(i)) {
        model.removeAll(i);
      }
    }
  }
}
