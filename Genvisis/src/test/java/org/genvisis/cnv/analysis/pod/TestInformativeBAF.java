/**
 * 
 */
package org.genvisis.cnv.analysis.pod;

import org.genvisis.cnv.analysis.pod.InformativeBAF.BAF_STRATEGY;
import org.genvisis.cnv.analysis.pod.InformativeBAF.InformativeResult;
import org.genvisis.common.ArrayUtils;
import org.junit.Test;
import junit.framework.Assert;

/**
 * 
 *
 */
public class TestInformativeBAF {

  private static final double DELTA = 1e-15;

  @Test
  public void TestGetInformativeIndices() {
    byte[] genos = ArrayUtils.byteArray(10, (byte) 1);
    double[] bafs = new double[] {.5, .5, .5, .5, .5, .5, .5, .5, .5, .75};
    InformativeResult informativeResult = InformativeBAF.getInformativeIndices(bafs, genos,
                                                                               BAF_STRATEGY.HET_ONLY,
                                                                               InformativeBAF.CHEBYSHEV);
    Assert.assertEquals(1, informativeResult.getInformatives().size());
    Assert.assertEquals(9, (int) informativeResult.getInformatives().get(0));
    Assert.assertEquals(0.525, informativeResult.getNormalDistribution().getMean(), DELTA);
  }
}
