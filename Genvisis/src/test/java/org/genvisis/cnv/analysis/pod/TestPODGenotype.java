/**
 * 
 */
package org.genvisis.cnv.analysis.pod;

import org.genvisis.cnv.analysis.pod.PODGenotype.BAF_EFFECT;
import org.genvisis.cnv.analysis.pod.PODGenotype.Genotype;
import org.genvisis.cnv.analysis.pod.PODGenotype.POD;
import org.junit.Test;

import junit.framework.Assert;

/**
 *
 *
 */
public class TestPODGenotype {

  @Test
  public void TestGetSharedAlleleCount() {

    // child has two alleles in common with parent
    Assert.assertEquals(2, Genotype.getSharedAlleleCount(Genotype.AA, Genotype.AA)
                                   .getNumAlleleMatch());
    Assert.assertEquals(2, Genotype.getSharedAlleleCount(Genotype.AA, Genotype.AB)
                                   .getNumAlleleMatch());

    Assert.assertEquals(2, Genotype.getSharedAlleleCount(Genotype.BB, Genotype.BB)
                                   .getNumAlleleMatch());
    Assert.assertEquals(2, Genotype.getSharedAlleleCount(Genotype.BB, Genotype.AB)
                                   .getNumAlleleMatch());
    Assert.assertEquals(2, Genotype.getSharedAlleleCount(Genotype.AB, Genotype.AB)
                                   .getNumAlleleMatch());
    // child has one allele in common with parent

    Assert.assertEquals(1, Genotype.getSharedAlleleCount(Genotype.AB, Genotype.AA)
                                   .getNumAlleleMatch());
    Assert.assertEquals(1, Genotype.getSharedAlleleCount(Genotype.AB, Genotype.BB)
                                   .getNumAlleleMatch());
    // child has zero alleles in common with parent
    Assert.assertEquals(0, Genotype.getSharedAlleleCount(Genotype.BB, Genotype.AA)
                                   .getNumAlleleMatch());
    Assert.assertEquals(0, Genotype.getSharedAlleleCount(Genotype.AA, Genotype.BB)
                                   .getNumAlleleMatch());

    Assert.assertEquals(0, Genotype.getSharedAlleleCount(Genotype.NC, Genotype.AA)
                                   .getNumAlleleMatch());
    Assert.assertEquals(0, Genotype.getSharedAlleleCount(Genotype.NC, Genotype.AB)
                                   .getNumAlleleMatch());
    Assert.assertEquals(0, Genotype.getSharedAlleleCount(Genotype.NC, Genotype.BB)
                                   .getNumAlleleMatch());
  }

  /**
   * Tests all mother and father genotype combinations for parental effect given a
   * {@link BAF_EFFECT}
   */
  @Test
  public void TestGetPodEffect() {
    for (Genotype g : PODGenotype.Genotype.values()) {
      Assert.assertEquals(true, PODGenotype.getPodEffect(g, g, BAF_EFFECT.INCREASE) == POD.NONE);
      Assert.assertEquals(true, PODGenotype.getPodEffect(g, g, BAF_EFFECT.DECREASE) == POD.NONE);
    }
    for (Genotype g : PODGenotype.Genotype.values()) {
      if (g != Genotype.BB) {
        Assert.assertEquals(true, PODGenotype.getPodEffect(Genotype.BB, g,
                                                           BAF_EFFECT.DECREASE) == POD.PATERNAL);

        Assert.assertEquals(true, PODGenotype.getPodEffect(g, Genotype.BB,
                                                           BAF_EFFECT.DECREASE) == POD.MATERNAL);
      }
      if (g != Genotype.AA) {

        Assert.assertEquals(true, PODGenotype.getPodEffect(Genotype.AA, g,
                                                           BAF_EFFECT.INCREASE) == POD.PATERNAL);
        Assert.assertEquals(true, PODGenotype.getPodEffect(g, Genotype.AA,
                                                           BAF_EFFECT.INCREASE) == POD.MATERNAL);
      }

    }
  }
}
