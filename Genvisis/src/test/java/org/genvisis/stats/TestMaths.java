package org.genvisis.stats;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import java.util.List;
import org.genvisis.common.ArrayUtils;
import org.junit.Test;
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

    double[] source = new double[] {6745.765290861431, 399.446777449024, Double.NaN, Double.NaN,
                                    5978.602833264176, 235.2119804304731, 1873.0687036512327,
                                    Double.NaN, 8542.247047534842, 2854.7826390217024, Double.NaN,
                                    5111.187824725547, 936.6129342903928, 3245.1533930032074,
                                    897.6403968390579, 2164.463885403172, 4589.933012315601,
                                    Double.NaN, 1575.1297841465894, 6996.704815846775,
                                    1767.9816510459557, 7756.303914906478, 5089.036354014447,
                                    7732.033937265825, 5923.967440141675, 2970.062149117565,
                                    4899.9805919384935, 6903.090064092241, Double.NaN, Double.NaN,
                                    675.1036363732787, Double.NaN, 8614.938673936316, Double.NaN,
                                    197.43573951636395, 8685.195748202139, 4509.759740429229,
                                    2690.8294714841595, 4909.041916884391, Double.NaN,
                                    7086.009863537849, 1410.0541644855302, Double.NaN,
                                    7838.5165489508145, Double.NaN, 8937.968864151002,
                                    3100.524686903092, Double.NaN, 136.1349919315402,
                                    362.6193547799006, 595.4658224795595, 1050.6755698521142,
                                    9070.890221460402, 2361.5718585527934, 7391.875785027233,
                                    8023.05136196327, Double.NaN, 9489.074113928471,
                                    696.532062372317, 3473.0656516674053, 841.1489745596956,
                                    3878.337928865955, Double.NaN, 7835.927931789041,
                                    6220.122696464281, Double.NaN, 8984.5812581381,
                                    3017.5774265266473, 1374.3070229753673, 3778.0506847255856,
                                    8272.471980545262, 3925.085128660507, 7593.5379331399845,
                                    8046.635677183436, 2369.4887555508426, 417.0123826011185,
                                    Double.NaN, 1973.5990151879212, 200.50457015032384,
                                    4862.613222717878, Double.NaN, 9660.162614028728,
                                    2695.6736038006525, 5235.019560911128, 4578.1476673102425,
                                    6927.233461729101, 658.6581573905281, Double.NaN,
                                    8941.116544336961, 8721.544405273502, 7578.430275975322};
    int window = 29;
    boolean skipNaN = true;
    assertArrayEquals(ArrayUtils.movingAverageForward(window, source, skipNaN),
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    source = new double[] {1.0, 3.7, 4.2, Double.NaN, 78.8, 10008.0004, 32.0, 41.3, 50.34};
    window = 3;
    skipNaN = true;
    assertArrayEquals(ArrayUtils.movingAverageForward(window, source, skipNaN),
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    skipNaN = false;
    assertArrayEquals(ArrayUtils.movingAverageForward(window, source, skipNaN),
                      Doubles.toArray(Maths.movingAverageForward(window, Doubles.asList(source),
                                                                 skipNaN)),
                      delta);
    for (int i = 0; i < 100; i++) {
      skipNaN = Math.random() > 0.5;
      int size = (int) (Math.random() * 100000.0) + 1000;
      window = (int) (Math.random() * (size / 2)) + 1;
      source = new double[size];
      for (int j = 0; j < source.length; j++) {
        if (Math.random() > 0.75) source[j] = Double.NaN;
        else source[j] = Math.random() * 10000;
      }

      long time = System.currentTimeMillis();
      double[] oldMethod = ArrayUtils.movingAverageForward(window, source, skipNaN);
      double oldTime = (System.currentTimeMillis() - time) / 1000.0;
      time = System.currentTimeMillis();
      double[] newMethod = Doubles.toArray(Maths.movingAverageForward(window,
                                                                      Doubles.asList(source),
                                                                      skipNaN));
      double newTime = (System.currentTimeMillis() - time) / 1000.0;
      System.out.println("Moving Average Time Comparison: ");
      System.out.println("Old: " + oldTime);
      System.out.println("New: " + newTime);
      System.out.println("Diff: " + (newTime - oldTime));
      System.out.println();

      assertArrayEquals(oldMethod, newMethod, delta);
    }
  }

}
