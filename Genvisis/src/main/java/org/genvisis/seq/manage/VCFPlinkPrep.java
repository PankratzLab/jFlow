/**
 * 
 */
package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Class to prepare a vcf for plink import
 *
 */
public class VCFPlinkPrep {

	private VCFPlinkPrep() {

	}



	/**
	 * @param vcf input vcf
	 * @param outDir output directory
	 * @param gqs array of GQ thresholds
	 */
	private static void prep(String vcf, String outDir, int[] gqs) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "plinkvcfPrep.log");
		for (int gq : gqs) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);
			String root = outDir + VCFOps.getAppropriateRoot(vcf, true) + "_GQ_" + gq;
			String filtVcf = root + ".vcf.gz";
			log.reportTimeInfo("Filtering "+ vcf + ", setting genotypes with GQ < " + gq
													+ " to missing in " + filtVcf);
			VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, filtVcf,
																																VCFOps.DEFUALT_WRITER_OPTIONS,
																																new Logger());

			int numScanned = 0;
			for (VariantContext vc : reader) {
				numScanned++;
				if (numScanned % 10000 == 0) {
					log.reportTimeInfo(numScanned + " variants scanned for GQ" + gq);
				}
				VariantContextBuilder builder = new VariantContextBuilder(vc);
				ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
				for (Genotype g : vc.getGenotypes()) {
					if (g.hasGQ() && g.getGQ() < gq) {
						GenotypeBuilder gb = new GenotypeBuilder(g);
						gb.alleles(GenotypeOps.getNoCall());
						genotypes.add(gb.make());
					} else {
						genotypes.add(g);
					}
				}

				GenotypesContext bc = GenotypesContext.create(genotypes);
				builder.genotypes(bc);
				writer.add(builder.make());
			}
			reader.close();
			writer.close();
			Files.writeArray(VCFOps.getSamplesInFile(vcf), root + ".sampleList.txt");
		}
	}

	public static void main(String[] args) {
		CLI c = new CLI(VCFPlinkPrep.class);



		c.addArgWithDefault(CLI.ARG_VCF, CLI.DESC_VCF, "a.vcf");
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "outDir/");
		c.addArgWithDefault("gqs",
												"comma delimited GQ scores to filter genotypes by. NOTE: only genotypes with a GQ entry will be subjected to this filter",
												"20,50");


		c.parseWithExit(args);

		String vcf = c.get(CLI.ARG_VCF);
		String outDir = c.get(CLI.ARG_OUTDIR);
		int[] gqs = ArrayUtils.toIntArray(c.get("gqs").split(","));

		prep(vcf, outDir, gqs);
	}


}
