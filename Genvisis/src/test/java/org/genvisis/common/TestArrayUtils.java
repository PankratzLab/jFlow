package org.genvisis.common;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import java.util.Arrays;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.junit.Test;
import org.pankratzlab.common.ArrayUtils;

/**
 * Tests for methods in {@link ArrayUtils}
 */
public class TestArrayUtils {

  @Test
  public void testAppendToArray() {
    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {"1", "2"}, new String[] {"3", "4"}));

    assertArrayEquals(new Object[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new Object[] {"1", "2"}, new Object[] {"3", "4"}));

    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {"1", "2"}, "3", "4"));

    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {"1", "2", "3"}, new String[] {"4"}));

    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {"1", "2", "3", "4"}, new String[] {}));

    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {"1", "2", "3"}, "4"));

    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {}, "1", "2", "3", "4"));

    assertArrayEquals(new String[] {"1", "2", "3", "4"},
                      ArrayUtils.appendToArray(new String[] {}, "1", "2", "3", "4"));
  }

  @Test
  public void testIntDeepCopy() {
    int[][][] source = {{{0, 0, 0}, {1, 2, 3}}, {{1}, {2}, {3}}};

    int[][][] copy = ArrayUtils.deepCopy(source);

    assertArrayEquals(source, copy);

    assertNotEquals(source, copy);

    assertNotEquals(source[1], copy[1]);

    assertNotEquals(source[0][0], copy[0][0]);

    assertEquals(source[0][0][0], copy[0][0][0]);

    source[0][0][0] = -1;

    assertNotEquals(source[0][0][0], copy[0][0][0]);

  }

  @Test
  public void testDoubleDeepCopy() {
    double[][][] source = {{{0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}}, {{1.0}, {2.0}, {3.0}}};

    double[][][] copy = ArrayUtils.deepCopy(source);

    assertArrayEquals(source, copy);

    assertNotEquals(source, copy);

    assertNotEquals(source[1], copy[1]);

    assertNotEquals(source[0][0], copy[0][0]);

    assertEquals(source[0][0][0], copy[0][0][0], 0.0);

    source[0][0][0] = -1.0;

    assertNotEquals(source[0][0][0], copy[0][0][0], 0.0);

  }

  @Test
  public void testMedian() {
    double[] test = new double[] {1, 2, 3};
    testMedianMethod(test);
    test = new double[] {1, 2, 3, 4};
    testMedianMethod(test);

    test = new double[] {1, 2, 3, 4};
    testMedianMethod(test);

    for (int i = 0; i < 1000; i++) {
      int size = (int) (Math.random() * 1000);
      test = new double[size];
      for (int j = 0; j < size; j++) {
        test[j] = Math.random() * i;
      }
      testMedianMethod(test);
    }
  }

  private static void testMedianMethod(double[] test) {
    double arrayVersion = new Median().evaluate(test);
    double streamVersion = ArrayUtils.median(Arrays.stream(test), test.length);
    assertEquals(arrayVersion, streamVersion, 0.000000001);
  }
}
