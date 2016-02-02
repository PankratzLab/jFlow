package seq.manage;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import common.Logger;

/**
 * Consolidating some common operations for tumor normal calling, summarizing filtering, and etc
 *
 */
public class VCFTumorNormalOps {
	public static void renameTumorNormalVCF(String vcf, String tumorSamp, String normalSamp, String output, Logger log) {
		renameTumorNormalVCF(vcf, tumorSamp, "TUMOR", normalSamp, "NORMAL", output, log);
	}

	/**
	 * @param vcf
	 *            vcf to rename
	 * @param tumorSamp
	 *            the tumor sample
	 * @param tumorDef
	 *            tumorSamp will replace this sample
	 * @param normalSamp
	 *            the normal sample
	 * @param normalDef
	 *            normalSamp will replace this sample
	 * @param output
	 *            output vcf
	 * @param log
	 */
	public static void renameTumorNormalVCF(String vcf, String tumorSamp, String tumorDef, String normalSamp, String normalDef, String output, Logger log) {

		if (VCFOps.getSamplesInFile(vcf).length != 2) {
			throw new IllegalArgumentException("This method is only sdesignd for tumor normal renaming");
		}
		VCFFileReader reader = new VCFFileReader(output, false);
		VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
		Set<String> samps = new HashSet<String>();
		samps.add(normalSamp);
		samps.add(tumorSamp);
		final VCFHeader outHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), samps);
		writer.writeHeader(outHeader);
		for (VariantContext vc : reader) {

			VariantContextBuilder builder = new VariantContextBuilder(vc);
			ArrayList<Genotype> renamed = new ArrayList<Genotype>();
			Genotype normal = rename(vc.getGenotype(normalDef), normalSamp);
			Genotype tumor = rename(vc.getGenotype(tumorDef), tumorSamp);
			renamed.add(normal);
			renamed.add(tumor);
			builder.genotypes(renamed);
			if (!renamed.get(0).sameGenotype(vc.getGenotype(normalDef))) {
				reader.close();
				writer.close();
				throw new IllegalStateException("Improprer rename");
			}
			builder.genotypes(renamed);
			if (!renamed.get(1).sameGenotype(vc.getGenotype(tumorDef))) {
				reader.close();
				writer.close();
				throw new IllegalStateException("Improprer rename");
			}

			writer.add(builder.make());
		}
		log.reportTimeInfo("Re-named and indexed " + vcf + " to " + output);
		reader.close();
		writer.close();

	}

	private static Genotype rename(Genotype g, String newName) {
		GenotypeBuilder builder = new GenotypeBuilder(g);
		builder.name(newName);
		return builder.make();
	}

	/**
	 * GATK appends .variant## to each sample when merging individual tumor normal calls, this will re-name the samples to the original (removes .variant.*)
	 */
	public static void renameMergeVCF(String inputVCF, String outputVCF) {
		VCFFileReader reader = new VCFFileReader(inputVCF, true);
		Set<String> samps = new HashSet<String>();
		String[] sampIn = VCFOps.getSamplesInFile(inputVCF);
		for (int i = 0; i < sampIn.length; i++) {
			String fix = sampIn[i].replaceAll(".variant.*", "");
			samps.add(fix);
		}

		final VCFHeader outHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), samps);

		VariantContextWriter writer = VCFOps.initWriter(outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
		writer.writeHeader(outHeader);
		for (VariantContext vc : reader) {
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			ArrayList<Genotype> renameGeno = new ArrayList<Genotype>();
			for (Genotype g : vc.getGenotypes()) {
				GenotypeBuilder gBuilder = new GenotypeBuilder(g);
				gBuilder.name(g.getSampleName().replaceAll(".variant.*", ""));
				renameGeno.add(gBuilder.make());
			}
			builder.genotypes(renameGeno);
			VariantContext vcFilt = builder.make();
			if (!vcFilt.isMonomorphicInSamples()) {
				writer.add(vcFilt);
			}
		}
		reader.close();
		writer.close();
	}

	// double mapQ = 0;
	// double ssc = 0;
	// try {
	// if (g.hasAnyAttribute("MQ")) {
	// mapQ = Double.parseDouble(g.getAnyAttribute("MQ").toString());
	// }
	// if (g.hasAnyAttribute("SSC")) {
	// ssc = Double.parseDouble(g.getAnyAttribute("SSC").toString());
	// }
	// } catch (NumberFormatException nfe) {
	//
	// }
	// if (mapQ < 40 || ssc < 40) {
	// gBuilder.alleles(GenotypeOps.getNoCall());
	// }

}
