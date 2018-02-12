package org.genvisis.stats;

import static org.junit.Assert.assertEquals;
import org.junit.Test;

/**
 * Tests for {@link FishersExact2by2Calculator}
 */
public class TestFishersExact {

  /**
   * Test basic expected values from
   * <a href="https://en.wikipedia.org/wiki/Fisher%27s_exact_test">the wikipedia article</a>.
   */
  @Test
  public void testExpected() {
    double expected = 0.000033652;

    // Test array input
    double pVal = FishersExact2by2Calculator.getPvalue(new int[][] {{0, 10}, {12, 2}});
    assertEquals(expected, pVal, 0.000000001);

    // Should be the same for individual values
    pVal = FishersExact2by2Calculator.getPvalue(0, 10, 12, 2);
    assertEquals(expected, pVal, 0.000000001);
  }
}
