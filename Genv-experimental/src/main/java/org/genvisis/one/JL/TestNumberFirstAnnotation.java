package org.genvisis.one.JL;

import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.genvisis.seq.qc.FilterNGS.VcFilterJEXL;
import org.pankratzlab.common.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

/**
 * Demo of why we force info fields (typically from ANNOVAR) to start with a letter e.g
 * 1000g2014oct_afr -> g10002014oct_afr We filter our {@link VariantContext}s with the a
 * {@link VariantContextFilter} that uses {@link VcFilterJEXL}. JEXL magic depends on the
 * {@link VariantContextUtils#initializeMatchExps(String[], String[])}, and
 * {@link TestNumberFirstAnnotation#fail()} demonstrates failure
 */
class TestNumberFirstAnnotation {

  private TestNumberFirstAnnotation() {

  }

  private static void fail() {

    // works
    VariantContextUtils.initializeMatchExps(new String[] {"Test"},
                                            new String[] {"g10002014oct_afr=='.'"});
    new Logger().reportTimeInfo("Success1");

    // no go
    VariantContextUtils.initializeMatchExps(new String[] {"Test"},
                                            new String[] {"1000g2014oct_afr=='.'"});
    new Logger().reportTimeInfo("Success2");

  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    fail();
  }

}
