package org.genvisis.stats;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;

import org.apache.commons.math3.stat.inference.TTest;
import org.genvisis.cnv.plots.stats.BatchEffects;
import org.genvisis.cnv.plots.stats.BatchEffectsBuilder;
import org.junit.Test;
import org.pankratzlab.common.Logger;

public class BatchEffectsTest {

  /**
   * Test of getNegLog10PValueMatrix method, of class BatchEffects.
   * 
   * @throws IOException
   * @throws ParseException
   * @throws FileNotFoundException
   */
  @Test
  public void getNegLog10PValueMatrixTEST() throws FileNotFoundException, ParseException,
                                            IOException {
    BatchEffectsBuilder builder = new BatchEffectsBuilder(new Logger());
    BatchEffects instance = builder.build("src/test/resources/BatchEffectsTEST_batch_input.txt",
                                          "src/test/resources/BatchEffectsTEST_factor_input.txt");
    String[][] negLog10PValueMatrix = instance.getNegLog10PValueMatrix(1E-300);
    double[] group1 = {0.010740732, 0.013343804, 0.025250002, 0.007692446, 0.006861333, 0.01084218};
    double[] group2 = {-0.015150839, -0.014317339, -0.005510375, -0.010304328, -0.007114341,
                       -0.004252832, 0.015766729, 0.016133045, 0.014680029, 0.007599692,
                       0.007581922, 0.00494472};
    TTest t = new TTest();
    double pValue = t.tTest(group1, group2);
    double expected = -1.0 * Math.log10(pValue);
    String s = negLog10PValueMatrix[1][1];
    double actual = Double.parseDouble(s);
    assertEquals(expected, actual, Double.MIN_VALUE);

    double[] group3 = {-0.003803057, -0.006277361, -0.007199143, 3.38E-04, -0.002828399,
                       -0.002692487};
    double[] group4 = {0.003361777, -0.006452373, 0.007192255, 0.00262634, 0.006921027,
                       -0.008258088, -0.005918212, 0.002103288, 0.010203016, -0.011983135,
                       -0.01139179, 0.0222021};
    pValue = t.tTest(group3, group4);
    expected = -1.0 * Math.log10(pValue);
    s = negLog10PValueMatrix[100][2];
    actual = Double.parseDouble(s);
    assertEquals(expected, actual, Double.MIN_VALUE);
  }
}
