
package org.genvisis.seq.manage.mtdna;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.genvisis.CLI;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Class for managing mitochondrial VCFs
 *
 * http://haplogrep.uibk.ac.at/blog/rcrs-vs-rsrs-vs-hg19/
 *
 * Note that this class is not a completely general converter, complex indels and some specific
 * cases will require more logic
 */
public class VCFOpsMT {
	private static final String HG_19_POS = "HG19POS";
	private static final String HG_19_REF = "HG19REF";
	private static final String HG_19_ALTS = "HG19ALTS";

	private static final String VCF = "vcf";
	private static final String OUTPUT = "out";

	private static final int[] BP_DIFFS =
																			new int[] {	73, 150, 195, 263, 408, 750, 1438, 2352, 2483,
																									2706, 4769, 5580, 7028, 8701, 8860, 9377, 9540,
																									10398, 10819, 10873, 11017, 11719, 11722, 12705,
																									12850, 14212, 14580, 14766, 14905, 15301, 15326,
																									15932, 16172, 16183, 16189, 16223, 16320, 16519};

	/**
	 * Which mitochondriome is represented
	 *
	 */
	public enum MT_GENOME {
													HG19, RCRS;
	}



	private static class AlleleSwapper {
		private final HashSet<Integer> toSwaps;
		private final RCRS rcrs;

		private AlleleSwapper(Logger log) {
			toSwaps = new HashSet<Integer>();
			for (int element : BP_DIFFS) {
				toSwaps.add(element);
			}
			rcrs = RCRS.getRCRS(log);
		}

		private boolean shouldSwap(VariantContext vc) {
			if (toSwaps.contains(vc.getStart())) {
				return true;
			}
			if (vc.isIndel()) {
				for (int i = vc.getStart(); i <= vc.getEnd(); i++) {
					if (toSwaps.contains(i)) {
						return true;
					}
				}
			}
			return false;
		}



		private VariantContext swap(VariantContext vc, MT_GENOME posFrom, Logger log) {
			if (posFrom != MT_GENOME.RCRS) {
				throw new IllegalArgumentException("swap method assumes " + MT_GENOME.RCRS + " positions");
			}
			if (shouldSwap(vc)) {
				return swapAlleles(vc, log);
			} else {
				Allele rcrsRef = getRcrsAllles(vc);
				if (!rcrsRef.basesMatch(vc.getReference())) {
					throw new IllegalStateException("Reference alleles for "	+ vc.toStringWithoutGenotypes()
																					+ "\n" + "do not match rcrs " + rcrsRef.getBaseString());
				}

				return vc;
			}
		}

		private VariantContext swapAlleles(VariantContext vc, Logger log) {
			VariantContextBuilder builder = new VariantContextBuilder(vc);

			Allele rcrsRef = getRcrsAllles(vc);
			if (vc.getReference().basesMatch(rcrsRef)) {
				throw new IllegalStateException("Expected allele changes match already");
			}
			if (!vc.getAlleles().contains(Allele.create(rcrsRef, true))) {

				log.reportError(vc.toStringWithoutGenotypes()	+ "\n -->>" + rcrsRef.getBaseString() + "\t"
												+ vc.getAlleles().contains(rcrsRef));
				throw new IllegalStateException("not surprisingly, assumptions about allele swapping have been violated");
			}
			ArrayList<Allele> newAlleles = new ArrayList<Allele>();
			newAlleles.add(Allele.create(vc.getReference().getBases(), false));
			newAlleles.add(rcrsRef);
			for (Allele prevAlt : vc.getAlternateAlleles()) {
				if (!prevAlt.basesMatch(rcrsRef)) {// don't double add
					newAlleles.add(Allele.create(prevAlt.getBases(), false));
				}
			}
			reDeriveGenotypes(vc, builder, rcrsRef);// assign the re-derived genotypes
			builder.alleles(newAlleles);// assign the new alleles
			return builder.make();
		}

		private void reDeriveGenotypes(	VariantContext vc, VariantContextBuilder builder,
																		Allele rcrsRef) {
			GenotypesContext gc = vc.getGenotypes();// original genotypes
			ArrayList<Genotype> newGenos = new ArrayList<Genotype>();
			for (Genotype g : gc) {// handle the new alleles on a genotype basis
				if (g.isCalled()) {
					GenotypeBuilder b = new GenotypeBuilder(g);
					List<Allele> old = g.getAlleles();
					ArrayList<Allele> newGAlleles = new ArrayList<Allele>();
					for (Allele a : old) {
						swapGenotypes(rcrsRef, newGAlleles, a);
					}
					b.alleles(newGAlleles);
					newGenos.add(b.make());
				} else {
					newGenos.add(g);
				}
			}
			builder.genotypes(newGenos);
		}

		private void swapGenotypes(Allele rcrsRef, ArrayList<Allele> newGAlleles, Allele a) {
			if (a.basesMatch(rcrsRef)) {// If the allele matches the rcrsReference, it's the reference
				newGAlleles.add(rcrsRef);
			} else {
				newGAlleles.add(Allele.create(a.getBases(), false));// otherwise it's the new alt
			}
		}

		private Allele getRcrsAllles(VariantContext vc) {
			String rcrsS = Array.toStr(	Array.subArray(rcrs.getBases(), vc.getStart() - 1, vc.getEnd()),
																	"");
			return Allele.create(rcrsS, true);
		}
	}


	/**
	 * Shifts base pairs according to table from
	 * http://haplogrep.uibk.ac.at/blog/rcrs-vs-rsrs-vs-hg19/
	 *
	 */
	private static class BpShifter {

		private BpShifter() {

		}

		// Special cases to deal with as they matter.
		// 309.1C,
		// 315.1C,
		// 3107del,
		// 16193.1C,



		private VariantContext shift(VariantContext vc, MT_GENOME from) {
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			builder.start(shift(vc.getStart(), from));
			builder.stop(shift(vc.getEnd(), from));
			if (vc.getStart() < 309 && vc.getEnd() > 309) {// cushings had a deletion that spanned the 309
																											// del
				builder.start((long) vc.getStart() - 1);
				builder.stop((long) vc.getEnd() - 1);
			}
			return builder.make();
		}

		/**
		 * @param bp the base pair to shift
		 * @param from the source of this position
		 * @return
		 */
		private int shift(int bp, MT_GENOME from) {
			int shiftbp = bp;
			if (bp < 310) {
				shiftbp = bp;
			} else if ((bp >= 316 && bp <= 3108) || bp >= 16194) {
				shiftbp = shiftIt(from, shiftbp, -2);
			} else if ((bp >= 3108 && bp < 16194) || (bp >= 310 && bp < 316)) {
				shiftbp = shiftIt(from, shiftbp, -1);
			} else {
				throw new IllegalStateException("Unhandled base pair " + bp);
			}
			return shiftbp;
		}

		private int shiftIt(MT_GENOME from, int bp, int shift) {
			int newBp;
			switch (from) {
				case HG19:
					newBp = bp + shift;
					break;
				case RCRS:
				default:
					throw new IllegalArgumentException("mtDNA shift not implemented for " + from);
			}
			return newBp;
		}

	}



	private static HashMap<String, Integer> getAlleleCounts(VariantContext vc) {
		HashMap<String, Integer> alleleCounts = new HashMap<String, Integer>();

		for (Genotype g : vc.getGenotypes()) {
			for (Allele a : g.getAlleles()) {
				String base = a.getDisplayString();
				if (alleleCounts.containsKey(base)) {
					alleleCounts.put(base, alleleCounts.get(base) + 1);
				} else {
					alleleCounts.put(base, 1);
				}
			}
		}
		return alleleCounts;
	}

	private static boolean sameAlleleCounts(VariantContext vc1, VariantContext vc2) {
		HashMap<String, Integer> alleleCounts1 = getAlleleCounts(vc1);
		HashMap<String, Integer> alleleCounts2 = getAlleleCounts(vc2);

		for (String allele1 : alleleCounts1.keySet()) {
			if (!alleleCounts2.containsKey(allele1)
					|| !alleleCounts1.get(allele1).equals(alleleCounts2.get(allele1))) {
				return false;
			}
		}

		return true;
	}

	/**
	 * @param inputVCF a vcf with chrM only entries aligned to HG19
	 * @param outputVCF where the rcrs conversion will be written
	 * @param log
	 */
	public static void convertHg19ToRCRS(String inputVCF, String outputVCF, Logger log) {
		// TODO handle indels that span change sites and complex alleles in general, currently writing
		// to error file

		VCFFileReader reader = new VCFFileReader(new File(inputVCF), false);
		new File(ext.parseDirectoryOfFile(outputVCF)).mkdirs();
		String errorOut = VCFOps.getAppropriateRoot(outputVCF, false) + ".switch.swap.error.vcf";

		VCFHeader header = getRcrsHeader(reader, MT_GENOME.RCRS);

		VariantContextWriter writer = VCFOps.initWriter(outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
																										header.getSequenceDictionary());
		VariantContextWriter writerError = VCFOps.initWriter(	errorOut, VCFOps.DEFUALT_WRITER_OPTIONS,
																													header.getSequenceDictionary());
		writer.writeHeader(header);
		writerError.writeHeader(header);

		BpShifter shifter = new BpShifter();
		AlleleSwapper alleleSwapper = new AlleleSwapper(log);
		for (VariantContext vc : reader) {
			if (!"chrM".equals(vc.getContig())) {
				writer.close();
				writerError.close();
				throw new IllegalArgumentException("This method designed for HG19 to rCRS conversion only");
			}
			VariantContext vcShift = null;
			VariantContext vcSwap = null;
			try {
				vcShift = shifter.shift(vc, MT_GENOME.HG19);
				vcSwap = alleleSwapper.swap(vcShift, MT_GENOME.RCRS, log);
			} catch (Exception e) {
				log.reportError("Could not convert " + vc.toStringWithoutGenotypes());
				writerError.add(vc);
				log.reportException(e);
			}
			if (vcSwap != null) {
				if (sameAlleleCounts(vc, vcSwap)) {// allele counts should match despite any
																						// alt/ref position switches
					VariantContextBuilder builder = new VariantContextBuilder(vcSwap);
					builder.chr("MT");// GRCh mitochondrial contig
					builder.attribute(HG_19_POS, vc.getStart());
					builder.attribute(HG_19_REF, vc.getReference().toString());
					builder.attribute(HG_19_ALTS, vc.getAlternateAlleles().toString().replaceAll(" ", ""));

					writer.add(builder.make());
				} else {
					log.reportError("Mismatched Allele counts");
					writer.close();
					writerError.close();
					throw new IllegalStateException("Mismatched Allele counts");
				}
			}
		}
		reader.close();
		writer.close();
		writerError.close();
	}


	private static VCFHeader getRcrsHeader(VCFFileReader reader, MT_GENOME mGenome) {
		VCFHeader header = new VCFHeader(reader.getFileHeader());
		ArrayList<SAMSequenceRecord> records = new ArrayList<SAMSequenceRecord>();
		switch (mGenome) {
			case HG19:
				records.add(new SAMSequenceRecord("chrM", 16571));
				break;
			case RCRS:
				records.add(new SAMSequenceRecord("MT", 16569));
				break;
			default:
				throw new IllegalArgumentException("Invalid type " + mGenome);
		}
		VCFInfoHeaderLine hg19Pos = new VCFInfoHeaderLine(HG_19_POS, 1, VCFHeaderLineType.String,
																											"HG 19 variant position");
		VCFInfoHeaderLine hg19REF = new VCFInfoHeaderLine(HG_19_REF, 1, VCFHeaderLineType.String,
																											"HG 19 reference allele");
		VCFInfoHeaderLine hg19Alts = new VCFInfoHeaderLine(	HG_19_ALTS, 1, VCFHeaderLineType.String,
																												"HG 19 alternate alleles");
		header.addMetaDataLine(hg19Pos);
		header.addMetaDataLine(hg19REF);
		header.addMetaDataLine(hg19Alts);


		header.setSequenceDictionary(new SAMSequenceDictionary(records));
		return header;
	}

	/**
	 * Convert the vcf to the appropriate contigs
	 * 
	 * @param inputVcf
	 * @param outputVCF
	 * @param mtGenome
	 * @param log
	 */
	public static void convertContigs(String inputVcf, String outputVCF, MT_GENOME mtGenome,
																		Logger log) {
		VCFFileReader reader = new VCFFileReader(new File(inputVcf), false);
		new File(ext.parseDirectoryOfFile(outputVCF)).mkdirs();
		log.reportTimeInfo("converting contigs in " + inputVcf + " to " + mtGenome);
		VCFHeader header = getRcrsHeader(reader, mtGenome);

		VariantContextWriter writer = VCFOps.initWriter(outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
																										header.getSequenceDictionary());

		writer.writeHeader(header);
		for (VariantContext vc : reader) {
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			switch (mtGenome) {
				case HG19:
					builder.chr("chrM");
					break;
				case RCRS:
					builder.chr("MT");
					break;
				default:
					writer.close();
					throw new IllegalArgumentException("Invalid type " + mtGenome);
			}
			writer.add(builder.make());

		}
		reader.close();
		writer.close();


	}

	public static void main(String[] args) {
		String vcf = "chrM.vcf";
		String out = "chrM.rcrs.vcf";
		CLI c = new CLI(VCFOpsMT.class);
		c.addArgWithDefault(VCF, "a vcf file with chrM only entries to convert to rcrs", vcf);
		c.addArgWithDefault(OUTPUT, "the rcrs output vcf", out);
		convertHg19ToRCRS(vcf, ext.addToRoot(vcf, ".switch.swap"), new Logger());
	}
}
