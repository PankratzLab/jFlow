package org.genvisis.stats;

import org.genvisis.stats.BinnedMovingStatistic.MovingStat;
import org.junit.Assert;
import org.junit.Test;

/**
 * Tests for {@link BinnedMovingStatistic}
 */
public class TestBinnedMovingStat {

  @Test
  public void testWindowEven() {
    BinnedMovingStatistic<Integer> bma = new BinnedMovingStatistic<Integer>(4, MovingStat.MEAN);
    bma.add(1, 1);
    bma.add(5, 2);
    // Need at least 2 values and the 2 bin should be under construction and not part of the stat
    // yet, so these should be -1
    Assert.assertEquals(-1, bma.mid());
    Assert.assertEquals(-1, bma.getValue(), 0.001);
    // The first and second bins should be combined to create the value for the first position
    bma.add(9, 3);
    Assert.assertEquals(1, bma.mid());
    Assert.assertEquals(3.0, bma.getValue(), 0.001);
    bma.add(9, 4);
    Assert.assertEquals(2, bma.mid());
    Assert.assertEquals(5.0, bma.getValue(), 0.001);
    // At this point we have a full window
    bma.add(1, 5);
    Assert.assertEquals(3, bma.mid());
    Assert.assertEquals(6.0, bma.getValue(), 0.001);
    bma.forceBinBreak();
    // The first bin should be ejected
    Assert.assertEquals(4, bma.mid());
    Assert.assertEquals(6.0, bma.getValue(), 0.001);
    // This pop should advance the index and pop off the oldest value, resulting in a size 3 mean
    bma.forceBinPop();
    Assert.assertEquals(5, bma.mid());
    Assert.assertEquals(6.33333, bma.getValue(), 0.001);
    // Because we have an even sized window, we need (window / 2 - 1) bins ahead and (window / 2)
    // bins behind to get a stat.
    // Advancing one more position should invalidate this stat.
    bma.forceBinPop();
    Assert.assertEquals(-1, bma.mid());
    Assert.assertEquals(-1, bma.getValue(), 0.001);
  }

  @Test
  public void testWindowOdd() {
    BinnedMovingStatistic<Integer> bma = new BinnedMovingStatistic<Integer>(3, MovingStat.MEAN);
    bma.add(5, 1);
    bma.add(9, 2);
    // Need at least 2 values and the 2 bin should be under construction and not part of the stat
    // yet, so these should be -1
    Assert.assertEquals(-1, bma.mid());
    Assert.assertEquals(-1, bma.getValue(), 0.001);
    // The first and second bins should be combined to create the value for the first position
    bma.add(9, 3);
    Assert.assertEquals(1, bma.mid());
    Assert.assertEquals(7.0, bma.getValue(), 0.001);
    // At this point we have a full window
    bma.add(1, 4);
    Assert.assertEquals(2, bma.mid());
    Assert.assertEquals(7.66666, bma.getValue(), 0.001);
    // The first bin should be ejected
    bma.forceBinBreak();
    Assert.assertEquals(3, bma.mid());
    Assert.assertEquals(6.33333, bma.getValue(), 0.001);
    bma.forceBinPop();
    Assert.assertEquals(4, bma.mid());
    Assert.assertEquals(5.0, bma.getValue(), 0.001);
    bma.forceBinPop();
    Assert.assertEquals(-1, bma.mid());
    Assert.assertEquals(-1, bma.getValue(), 0.001);
  }

}
