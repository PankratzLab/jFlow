package org.genvisis.one.JL;

import java.io.File;

import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Test filtering of SV vcf using svtyper GQ of lumpy calls
 *
 */
public class LumpyFilter {

	private LumpyFilter() {

	}

	public static void main(String[] args) {
		String vcf = "/Volumes/Beta/data/aric_sra/SRAPipeline/private/testBams/lumpyTest/lumpyH_UM-Schiffman-129-SS-129.gt.sort.vcf";
		VCFOps.verifyIndex(vcf, new Logger());
		VCFFileReader reader = new VCFFileReader(new File(vcf));
		VariantContextWriter writer = VCFOps.initWriter(VCFOps.getAppropriateRoot(vcf, false) + ".filt.vcf",
				VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
		for (VariantContext vc : reader) {
			if (vc.getGenotype(0).getGQ() > 150) {
				writer.add(vc);
			}
		}
		reader.close();
		writer.close();
	}

}
