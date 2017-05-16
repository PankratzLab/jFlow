package org.genvisis.one.JL;

import java.io.File;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps;

import htsjdk.samtools.SAMSequenceDictionary;
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

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String vcf = "/Volumes/Beta/data/aric_sra/SRAPipeline/private/testBams/lumpyTest/H_UM-Schiffman-129-SS-129lumpy.gt.sort.vcf";
		Logger log = new Logger(ext.parseDirectoryOfFile(vcf) + "filt.log");
		SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome("/Volumes/Beta/ref/all_sequences.fa", log)
				.getIndexedFastaSequenceFile().getSequenceDictionary();
		VCFOps.verifyIndex(vcf, new Logger());
		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		VariantContextWriter writer = VCFOps.initBuilder(VCFOps.getAppropriateRoot(vcf, false) + ".filt.vcf",
				VCFOps.DEFUALT_WRITER_OPTIONS, samSequenceDictionary).build();
		for (VariantContext vc : reader) {
			if (vc.getGenotype(0).getGQ() > 150) {
				writer.add(vc);
			}
		}
		reader.close();
		writer.close();
	}

}
