package org.genvisis.stats;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.pankratzlab.shared.stats.Maths;
import com.google.common.collect.ImmutableList;
import com.google.common.primitives.Doubles;

/**
 * Tests for the {@link Maths} class
 */
public class TestMaths {

  private static final List<Integer> POWERS_OF_10 = ImmutableList.of(1, 10, 100, 1000, 10000);
  private static final List<Integer> NOT_POWERS_OF_10 = ImmutableList.of(0, -1, -100, -2, 2, 15, 20,
                                                                         110, 111, 101);

  /**
   * Ensure that {@link Maths#isPowerOf10(int)} works as expected
   */
  @Test
  public void isPowerOf10Test() {
    for (int powerOf10 : POWERS_OF_10) {
      assertTrue(powerOf10 + " is a power of 10", Maths.isPowerOf10(powerOf10));
    }
    for (int notPowerOf10 : NOT_POWERS_OF_10) {
      assertFalse(notPowerOf10 + " is not a power of 10", Maths.isPowerOf10(notPowerOf10));
    }
  }

  @Test
  public void testMovingAverage() {

    double delta = 0.00001;

    double[] source = new double[] {1.0, 1.0, 1.0, 1.0};
    int window = 4;
    boolean skipNaN = true;
    assertArrayEquals(source,
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    skipNaN = false;
    assertArrayEquals(source,
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    source = new double[] {1.0, 2.0, 3.0, Double.NaN, 5.0};
    window = 3;
    skipNaN = true;
    assertArrayEquals(new double[] {1.0, 1.5, 2.0, Double.NaN, 3.3333333333333},
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    skipNaN = false;
    assertArrayEquals(new double[] {1.0, 1.5, 2.0, 2.5, 4.0},
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    // Test against apache commons
    for (int i = 0; i < 100; i++) {
      skipNaN = Math.random() > 0.5;
      int size = (int) (Math.random() * 1000.0) + 100;
      window = (int) (Math.random() * (size / 2)) + 1;
      source = new double[size];
      for (int j = 0; j < source.length; j++) {
        source[j] = Math.random() * 10000;
      }

      DescriptiveStatistics descStats = new DescriptiveStatistics(window);
      double[] apacheMethod = new double[size];
      for (int j = 0; j < apacheMethod.length; j++) {
        descStats.addValue(source[j]);
        apacheMethod[j] = descStats.getMean();
      }
      double[] mathsMethod = Doubles.toArray(Maths.movingAverageForward(window,
                                                                        Doubles.asList(source),
                                                                        skipNaN));

      assertArrayEquals(apacheMethod, mathsMethod, delta);
    }
  }

}
