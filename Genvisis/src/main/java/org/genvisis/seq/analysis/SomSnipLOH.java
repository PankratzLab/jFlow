package org.genvisis.seq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.GenotypeOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Extract variants called LOH by somatic sniper
 *
 */
public class SomSnipLOH {

	private static void extract(String vcf) {
		String root = VCFOps.getAppropriateRoot(vcf, false);
		Logger log = new Logger(root + "LOH.log");
		String out = root + "LOH.summary.txt";
		try {
			String[][] varAnno = VCFOps.getAnnotationKeys(vcf, log);
			String[][] geneAnno = GenotypeOps.getGenoFormatKeys(vcf, log);

			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("CHROM\tPOS\tID\tREF\tALT\tFILTER\tSAMPLE\t"	+ Array.toStr(geneAnno[1]) + "\t"
											+ Array.toStr(varAnno[1]));
			writer.println("CHROM\tPOS\tID\tREF\tALT\tFILTER\tSAMPLE\t"	+ Array.toStr(geneAnno[0]) + "\t"
											+ Array.toStr(varAnno[0]));

			VCFFileReader reader = new VCFFileReader(new File(vcf), true);
			for (VariantContext vc : reader) {
				String base = vc.getContig()	+ "\t" + vc.getStart() + "\t" + vc.getID() + "\t"
											+ vc.getReference().getBaseString() + "\t"
											+ vc.getAlternateAlleles().toString() + "\t" + vc.getFilters().toString();
				String[] vcAnnot = VCOps.getAnnotationsFor(varAnno[0], vc, ".");

				for (Genotype g : vc.getGenotypes()) {
					if (!g.isNoCall() && !g.isHomRef()) {
						writer.println(base	+ "\t" + g.getSampleName() + "\t"
														+ Array.toStr(GenotypeOps.getGenoAnnotationsFor(geneAnno[0], g, "."))
														+ "\t" + Array.toStr(vcAnnot));
					}
				}
			}
			writer.close();
			reader.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + out);
			log.reportException(e);
		}
	}

	public static void main(String[] args) {
		String vcf = args[0];
		extract(vcf);
	}

}
