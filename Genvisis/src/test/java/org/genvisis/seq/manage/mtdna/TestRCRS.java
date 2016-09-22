/**
 * 
 */
package org.genvisis.seq.manage.mtdna;

import org.genvisis.common.Logger;
import org.junit.Test;

import junit.framework.Assert;

/**
 * @author Kitty
 *
 */
public class TestRCRS {

  private static final String NAME = "chrM";// name
                                            // of
  private static final int LENGTH = 16569; // rcrs length
  private static final String BP_750 = "A";// diff between hg19 and rcrs
  private static final String BP_1018 = "A";// diff between hg19 and rcrs
  private static final String BP_13506 = "T";// diff between hg19 and rcrs

  /**
   * Ensure that we are likely loading the correct rcrs reference sequence
   * 
   */
  @Test
  public void rcrsTest() {
    RCRS rcrs = RCRS.getRCRS(new Logger());
    Assert.assertEquals(NAME, rcrs.getRcrsRef().getName());
    Assert.assertEquals(LENGTH, rcrs.getRcrsRef().length());
    Assert.assertEquals(LENGTH, rcrs.getBases().length);
    Assert.assertEquals(BP_750, rcrs.getBases()[750]);
    Assert.assertEquals(BP_1018, rcrs.getBases()[1018]);
    Assert.assertEquals(BP_13506, rcrs.getBases()[13506]);
  }

}
