package org.genvisis.seq.qc;


import org.genvisis.common.Logger;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;

/**
 * @author lane0212 Class to mimic the variant filter reported for the CHARGE participants from
 *         ARIC, CHS and FHS,
 */
public class AricWesFilter {
  private static final String ARIC_FILTER_NAME = "ARIC_FILTER";
  private final Logger log;

  public AricWesFilter(Logger log) {

    this.log = log;
  }

  public VariantFilterSample getARICVariantContextFilters() {
    VariantFilterSample.Builder builder = new VariantFilterSample.Builder();

    // whole variant filters
    // bi-allelic

    VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
    // Missing <.20
    VARIANT_FILTER_DOUBLE callrate = VARIANT_FILTER_DOUBLE.CALL_RATE;

    callrate.setDFilter(.80);

    // Hardy-Weinberg equilibrium expectations (p<5x10-6)

    VARIANT_FILTER_DOUBLE hwe = VARIANT_FILTER_DOUBLE.HWE;
    hwe.setDFilter(0.000005);

    builder.vCFiltWholeVariant(new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {callrate, hwe},
                                                        new VARIANT_FILTER_BOOLEAN[] {biallelic},
                                                        null, null, log));

    // alt call filters....
    // low variant read count (<3),
    VARIANT_FILTER_DOUBLE altDP = VARIANT_FILTER_DOUBLE.ALT_ALLELE_DEPTH;
    altDP.setDFilter(3);

    // variant read ratio <0.25 or >0.75
    VARIANT_FILTER_DOUBLE altRatLow = VARIANT_FILTER_DOUBLE.HET_ALLELE_RATIO_LOW;
    altRatLow.setDFilter(.25);
    VARIANT_FILTER_DOUBLE altRatHigh = VARIANT_FILTER_DOUBLE.HET_ALLELE_RATIO_HIGH;
    altRatHigh.setDFilter(.75);

    builder.vCFiltAltGeno(new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {altDP, altRatLow,
                                                                                altRatHigh},
                                                   new VARIANT_FILTER_BOOLEAN[] {}, null, null,
                                                   log));

    // hom ref filters....

    // any genotype filter
    // coverage less than 10-fold
    VARIANT_FILTER_DOUBLE totalDepth = VARIANT_FILTER_DOUBLE.DP;
    totalDepth.setDFilter(10);

    VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_STRICT;
    gq.setDFilter(95);

    builder.vCFiltAllGeno(new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {gq, totalDepth},
                                                   new VARIANT_FILTER_BOOLEAN[] {}, null, null,
                                                   log));

    // indel filter
    // the same for indels except a total coverage less than 30-fold was used.
    VARIANT_FILTER_DOUBLE indelDepth = VARIANT_FILTER_DOUBLE.DP;
    totalDepth.setDFilter(30);

    builder.vCFiltIndel(new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {indelDepth},
                                                 new VARIANT_FILTER_BOOLEAN[] {}, null, null, log));

    return builder.build(ARIC_FILTER_NAME, log);
  }

}
