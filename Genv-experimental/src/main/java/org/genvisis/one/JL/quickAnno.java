package org.genvisis.one.JL;

import org.genvisis.seq.analysis.GATK_Genotyper;
import org.genvisis.seq.manage.CreateNonSiteOnlyVcf;
import org.genvisis.seq.manage.VCFOps;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.core.CLI;

/**
 * annotate a vcf using hard coded defualts...need I say this is dangerous? Need to sort, use
 * http://vcftools.sourceforge.net/perl_module.html vcf-sort file.vcf.gz
 */
public class quickAnno {

  private quickAnno() {

  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(quickAnno.class);

    c.addArgWithDefault("vcf", "vcf to annotate with default methods", "a.vcf");
    c.parseWithExit(args);
    String vcf = c.get("vcf");
    Logger log = new Logger(VCFOps.getAppropriateRoot(vcf, false) + ".anno.log");
    if (VCFOps.getSamplesInFile(vcf).length == 0) {
      vcf = CreateNonSiteOnlyVcf.createNonSiteOnlyVcf(vcf, false, false);
    }
    GATK_Genotyper.annotateOnlyWithDefualtLocations(vcf, null, PSF.Ext.DEFAULT_MEMORY_MB, true,
                                                    false, log);
  }

}
