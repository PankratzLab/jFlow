/**
 * 
 */
package org.genvisis.seq.manage;

import org.genvisis.seq.manage.SeqOps.GC_COMP_METHOD;
import org.junit.Test;

import junit.framework.Assert;

/**
 * 
 */
public class TestSeqOps {

  private static final double DELTA = 1e-15;

  /**
   * Test for {@link SeqOps#countMotif(String, String)}
   */
  @Test
  public void testCountMotif() {
    String motif = "TTAGGG";
    String seq = "TTAGGGdasfTTAGGGDSFTTAGGg";
    Assert.assertEquals(3, SeqOps.countMotif(seq, motif, false));
    Assert.assertEquals(2, SeqOps.countMotif(seq, motif, true));

  }

  @Test
  public void testGetProportionGC() {
    String[] seq = "AAAAAAAAAAA".split("");
    Assert.assertEquals(0, SeqOps.getProportionGC(seq, GC_COMP_METHOD.GCTA_ONLY), DELTA);
    Assert.assertEquals(0, SeqOps.getProportionGC(seq, GC_COMP_METHOD.N_COUNTS_FOR_TOTAL), DELTA);

    seq = "GNNgCNNc".split("");
    Assert.assertEquals(1, SeqOps.getProportionGC(seq, GC_COMP_METHOD.GCTA_ONLY), DELTA);
    Assert.assertEquals(.5, SeqOps.getProportionGC(seq, GC_COMP_METHOD.N_COUNTS_FOR_TOTAL), DELTA);

    seq = "GTagCAtc".split("");
    Assert.assertEquals(.5, SeqOps.getProportionGC(seq, GC_COMP_METHOD.GCTA_ONLY), DELTA);
    Assert.assertEquals(.5, SeqOps.getProportionGC(seq, GC_COMP_METHOD.N_COUNTS_FOR_TOTAL), DELTA);

  }

}
