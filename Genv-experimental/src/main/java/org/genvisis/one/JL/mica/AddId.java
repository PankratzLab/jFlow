package org.genvisis.one.JL.mica;

import java.io.File;

import org.genvisis.CLI;
import org.genvisis.common.Logger;
import org.genvisis.one.JL.quickAnno;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class AddId {
	private AddId() {

	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		CLI c = new CLI(AddId.class);

		c.addArgWithDefault("vcf", "vcf to annotate with default methods", "a.vcf");
		c.parseWithExit(args);
		
		String vcfFile = c.get("vcf");
		VCFOps.verifyIndex(vcfFile, new Logger());
		String out = VCFOps.getAppropriateRoot(vcfFile, false) + ".ids.vcf";
		VCFFileReader reader = new VCFFileReader(new File(vcfFile), true);
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, out, VCFOps.DEFUALT_WRITER_OPTIONS,
				new Logger());
		for (VariantContext vc : reader) {
			String id = VCOps.getAnnotationsFor(new String[] { "snp138" }, vc, ".")[0];
			if (".".equals(id)) {
				id = new VCOps.LocusID(vc).getId();
			}
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			builder.id(id + "_" + VCOps.getAnnotationsFor(new String[] { "EFF" }, vc, ".")[0].split("\\(")[0]);
			writer.add(builder.make());
		}
		reader.close();
		writer.close();

	}
}
