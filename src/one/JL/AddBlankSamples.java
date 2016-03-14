package one.JL;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import seq.manage.VCFOps;
import common.Logger;
import common.ext;

public class AddBlankSamples {

	public static void main(String[] args) {
		String vcf = "/panfs/roc/groups/5/pankrat2/shared/CHARGE/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";
		String segFile = "/panfs/roc/groups/5/pankrat2/shared/CHARGE/Blanks/segs.txt";
		int numBlanks = 5;
		String outDir = ext.parseDirectoryOfFile(vcf) + "Blanks/";
		new File(outDir).mkdirs();
		String outVcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".blanks.vcf.gz";
		if (!VCFOps.existsWithIndex(outVcf)) {
			VCFFileReader reader = new VCFFileReader(vcf, true);
			VariantContextWriter writer = VCFOps.initWriter(outVcf, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());

			Set<String> samps = new HashSet<String>();
			for (int i = 0; i < numBlanks; i++) {
				samps.add("Blank_" + i);
			}
			Set<VCFHeaderLine> newHeaderLines = new HashSet<VCFHeaderLine>();
			newHeaderLines.addAll(reader.getFileHeader().getInfoHeaderLines());
			newHeaderLines.addAll(reader.getFileHeader().getOtherHeaderLines());
			newHeaderLines.addAll(reader.getFileHeader().getContigLines());
			newHeaderLines.addAll(reader.getFileHeader().getFilterLines());
			VCFFormatHeaderLine newFormat = new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "A blank genotype");
			newHeaderLines.add(newFormat);

			final VCFHeader outHeader = new VCFHeader(newHeaderLines, samps);

			writer.writeHeader(outHeader);
			for (VariantContext vc : reader) {
				VariantContextBuilder builder = new VariantContextBuilder(vc);
				ArrayList<Genotype> genotypesBlanks = new ArrayList<Genotype>();
				for (int i = 0; i < numBlanks; i++) {
					GenotypeBuilder gb = new GenotypeBuilder("Blank_" + i, vc.getAlleles());
					genotypesBlanks.add(gb.make());
				}
				GenotypesContext gc = GenotypesContext.create(genotypesBlanks);
				builder.genotypes(gc);
				writer.add(builder.make());
			}
			reader.close();
			writer.close();
		}
		VCFOps.extractSegments(outVcf, segFile, 10, null, outDir, false, false, true, false, null, 1, new Logger());
	}

}
