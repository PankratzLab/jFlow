package org.genvisis.common;

import static org.junit.Assert.assertEquals;
import org.junit.Test;

/**
 * Test for {@link ext} class
 */
public class Testext {

  @Test
  public void testIndexOfStartsWith() {
    assertEquals(1, ext.indexOfStartsWith("CAPS", new String[] {"CaPStone", "CAPSTONE", "CAP_NOPE"},
                                          false, true));
    assertEquals(0, ext.indexOfStartsWith("CAPS", new String[] {"CaPStone", "CAPSTONE", "CAP_NOPE"},
                                          false, false));
    assertEquals(-1,
                 ext.indexOfStartsWith("CAPS", new String[] {"CaPStone", "CAP=STONE", "CAP_NOPE"},
                                       false, true));

    assertEquals(1, ext.indexOfStartsWith("CAPSTONE", new String[] {"CaPS", "CAPS", "CAP_NOPE"},
                                          true, true));
    assertEquals(0, ext.indexOfStartsWith("CAPSTONE", new String[] {"CaPS", "CAPS", "CAP_NOPE"},
                                          true, false));
    assertEquals(-1, ext.indexOfStartsWith("CAPSTONE", new String[] {"CaPS", "CAP=S", "CAP_NOPE"},
                                           true, true));
  }

}
