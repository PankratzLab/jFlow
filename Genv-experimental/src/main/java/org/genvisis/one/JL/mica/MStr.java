package org.genvisis.one.JL.mica;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.common.Files;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MStr {

	private static void run(String vcf, String outDir) {
		new File(outDir).mkdirs();
		Segment mSeg = new Segment("chr6:31380150-31380165");
		String outfile = outDir + "mica.matrix.txt";
		PrintWriter writer = Files.getAppropriateWriter(outfile);

		VCFFileReader reader = new VCFFileReader(new File(vcf), false);
		String[] samples = VCFOps.getSamplesInFile(reader);

		writer.println("CHR\tPOS\tREF\tSAMPLE\tA1\tA2");
		for (VariantContext vc : reader) {
			if (VCOps.getSegment(vc).overlaps(mSeg)) {
				if (vc.isIndel()) {
					int index = 0;
					StringBuilder builder = new StringBuilder(
							vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getReference().getBaseString());
					for (Genotype g : vc.getGenotypes()) {

						if (!g.getSampleName().equals(samples[index])) {
							writer.close();

							throw new IllegalStateException("Sample mismatch");
						}
						index++;
						StringBuilder builder2 = new StringBuilder(builder.toString());
						builder2.append("\t" + g.getSampleName());

						if (g.isCalled()) {

							ArrayList<String> alts = new ArrayList<String>();
							String ref = null;
							for (Allele a : g.getAlleles()) {
								if (a.isNonReference()) {
									alts.add(a.getBaseString());
								} else {
									ref = a.getBaseString();
								}
							}
							if (ref != null) {
								if (alts.isEmpty()) {
									builder2.append("\t" + ref + "\t" + ref);
								} else {
									builder2.append("\t" + ref + "\t" + alts.get(0));
									if (alts.size() != 1) {
										writer.close();

										throw new IllegalStateException("geno mismatch");
									}
								}
							} else {
								if (alts.size() != 2) {
									writer.close();

									throw new IllegalStateException("geno mismatch");

								}
								for (String al : alts) {
									builder2.append("\t" + al);
								}
							}
						} else {
							builder2.append("\t.\t.");
						}
						writer.println(builder2.toString());
					}
				}
			}
		}
		writer.close();
		reader.close();

	}

	public static void main(String[] args) {
		CLI c = new CLI(ABLookup.class);
		c.addArgWithDefault("vcf", "vcf file to explore", null);
		c.addArgWithDefault("outDir", "output directory", null);
		c.parseWithExit(args);
		run(c.get("vcf"), c.get("outDir"));
	}

}
