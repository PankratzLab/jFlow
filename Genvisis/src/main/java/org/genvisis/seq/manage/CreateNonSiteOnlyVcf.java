/**
 * 
 */
package org.genvisis.seq.manage;

import java.io.File;
import java.util.HashSet;

import org.genvisis.CLI;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * @author takes a site only vcf, and appends a het genotype to all sites...
 *
 */
public class CreateNonSiteOnlyVcf {

	public static void main(String[] args) {
		CLI c = new CLI(CreateNonSiteOnlyVcf.class);



		c.addArgWithDefault("vcf", "vcf to generate non -site only output", "a.vcf");
		c.parseWithExit(args);

		String vcf = c.get("vcf");
		VCFFileReader reader = new VCFFileReader(new File(vcf), false);
		if (reader.getFileHeader().getNGenotypeSamples() > 0) {
			reader.close();
			throw new IllegalArgumentException("VCF "+ vcf
																					+ " did not appear to be site only, will not add mock sample");

		}
		Logger log = new Logger(ext.parseDirectoryOfFile(vcf) + "log.log");
		String randSample = "FAKE_SAMPLE_WITH_HET_ONLY_CALLS"; // will create all het calls for this
																														// sample
		HashSet<String> samps = new HashSet<String>();
		samps.add(randSample);

		String output = VCFOps.getAppropriateRoot(vcf, false) + randSample + ".vcf.gz";
		VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS,
																										reader.getFileHeader().getSequenceDictionary());
		VCFHeader header = reader.getFileHeader();

		VCFFormatHeaderLine format = new VCFFormatHeaderLine("GT", -1, VCFHeaderLineType.String, "GT");
		VCFInfoHeaderLine info = new VCFInfoHeaderLine(	"LEN", -1, VCFHeaderLineType.String,
																										"Don't know what this is, was not properly in the header -jlanej");
		VCFInfoHeaderLine info2 =
														new VCFInfoHeaderLine("TYPE", -1, VCFHeaderLineType.String,
				"Don't know what this is, was not properly in the header -jlanej");
		VCFInfoHeaderLine info3 =
														new VCFInfoHeaderLine("OLD_VARIANT", -1, VCFHeaderLineType.String,
																									"Don't know what this is, was not properly in the header -jlanej");



		header.addMetaDataLine(format);
		header.addMetaDataLine(info);
		header.addMetaDataLine(info2);
		header.addMetaDataLine(info3);


		writer.writeHeader(new VCFHeader(header.getMetaDataInInputOrder(), samps));

		int count = 0;
		for (VariantContext vc : reader) {
			count++;
			if (count % 100000 == 0) {
				log.reportTimeInfo("Converted " + count);
			}

			VariantContextBuilder builder = new VariantContextBuilder(vc);
			builder.genotypes(GenotypeBuilder.create(randSample, vc.getAlleles()));
			writer.add(builder.make());
		}
		reader.close();
		writer.close();


	}

}
