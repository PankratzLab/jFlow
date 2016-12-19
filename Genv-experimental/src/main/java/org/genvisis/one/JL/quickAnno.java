package org.genvisis.one.JL;

import org.genvisis.CLI;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.seq.analysis.GATK_Genotyper;

/**
 * annotate a vcf using hard coded defualts...need I say this is dangerous?
 *
 */
public class quickAnno {

	public static void main(String[] args) {
		CLI c = new CLI(quickAnno.class);

		c.addArgWithDefault("vcf", "vcf to annotate with default methods", "a.vcf");
		c.parseWithExit(args);
		GATK_Genotyper.annotateOnlyWithDefualtLocations(c.get("vcf"), null, PSF.Ext.DEFAULT_MEMORY_MB,
																										true, false, new Logger());
	}

}
