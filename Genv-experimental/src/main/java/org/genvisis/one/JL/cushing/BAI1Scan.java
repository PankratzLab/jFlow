package org.genvisis.one.JL.cushing;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringJoiner;
import java.util.TreeSet;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

class BAI1Scan {

	private static final String[] BASIC_INFO = new String[] { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER" };
	private static final String[] SAMP_INFO = new String[] { "Allele1", "Allele2", "AlleleCombo", "GQ", "PID", "PS",
			"NE", "PR", "GT" };

	private static String[] getBasicInfo(VariantContext vc) {
		return new String[] { vc.getContig(), Integer.toString(vc.getStart()), vc.getID(),
				vc.getReference().getBaseString(), vc.getAlternateAlleles().toString(),
				Double.toString(vc.getPhredScaledQual()), vc.getFilters().toString() };

	}

	private static List<String> getSampHeader(Set<String> samps) {
		List<String> header = new ArrayList<>();
		for (String samp : samps) {
			for (String inf : SAMP_INFO) {
				header.add(inf + "_" + samp);
			}

		}
		return header;

	}

	private static class SampInf {
		private final List<String> infs;
		private final boolean hasPID;
		private final boolean hasPS;

		private SampInf(List<String> infs, boolean hasPID, boolean hasPS) {
			super();
			this.infs = infs;
			this.hasPID = hasPID;
			this.hasPS = hasPS;
		}

	}

	private static SampInf getSampInfo(Genotype g) {
		List<String> info = new ArrayList<>();
		Allele a1;
		Allele a2;
		if (g.getAlleles().size() == 2) {
			a1 = g.getAllele(0);
			a2 = g.getAllele(1);
			// if (!a1.isReference() && a2.isReference()) {
			// throw new IllegalArgumentException("Invalid alleles: invalid
			// order");
			//
			// }
		} else {
			throw new IllegalArgumentException("Invalid alleles");
		}
		info.add(a1.getBaseString());
		info.add(a2.getBaseString());
		info.add(a1.getBaseString() + a2.getBaseString());

		info.add(Integer.toString(g.getGQ()));

		info.add(g.hasAnyAttribute("PID") ? g.getAnyAttribute("PID").toString() : "NA");
		info.add(g.hasAnyAttribute("PS") ? g.getAnyAttribute("PS").toString() : "NA");
		info.add(g.hasAnyAttribute("NE") ? g.getAnyAttribute("NE").toString() : "NA");
		info.add(g.hasAnyAttribute("PR") ? g.getAnyAttribute("PR").toString() : "NA");

		info.add(g.toString());
		return new SampInf(info, g.hasAnyAttribute("PID"), g.hasAnyAttribute("PS"));

	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String outDir = "/Volumes/Beta/data/Cushings/BAI1LD/scans/";

		String[] initVCFs = new String[] { "/Volumes/Beta/data/Cushings/BAI1LD/Beagle/out.cushing.only.gt.phase.vcf.gz",
				"/Volumes/Beta/data/Cushings/BAI1LD/Beagle/bai1RecombinationNeutralCushingOnly.cr.0.95.vcf.gz",
				"/Volumes/Beta/data/Cushings/BAI1LD/Beagle/out.cushing.only.gt.noRef.phase.vcf.gz",
				"/Volumes/Beta/data/Cushings/BAI1LD/HapCut/bai1RecombinationNeutralCushingOnly.HapCut.phased.vcf.gz" };
		double gq = -1;
		for (String initVCF : initVCFs) {

			new File(outDir).mkdirs();
			Segment seg = new Segment("chr8:143249162-143750071");
			SortedSet<String> sampsInf = new TreeSet<>();
			sampsInf.add("PT352_03_TAAGGCGA-AGAGTAGA.variant");
			sampsInf.add("PT364_03_CGAGGCTG-TATCCTCT.variant");
			sampsInf.add("PT368_03_CGTACTAG-TAGATCGC.variant");
			// sampsInf.add("PT142_03_CAGAGAGG-TATCCTCT.variant");

			SortedSet<String> all = new TreeSet<>();
			String[] vcfsamps = VCFOps.getSamplesInFile(initVCF);
			for (String samp : vcfsamps) {
				all.add(samp);
			}
			String output = outDir + ext.replaceWithLinuxSafeCharacters(seg.getUCSClocation(), true) + "_"
					+ ext.rootOf(initVCF) + "_stats.txt";

			parse(initVCF, seg, all, output, sampsInf, gq);
			String outputLimit = outDir + ext.replaceWithLinuxSafeCharacters(seg.getUCSClocation(), true) + "_"
					+ ext.rootOf(initVCF) + "_limitstats.txt";
			parse(initVCF, seg, sampsInf, outputLimit, sampsInf, gq);

		}

	}

	private static void parse(String initVCF, Segment seg, SortedSet<String> sampsToScan, String output,
			SortedSet<String> sampsInf, double gq) {
		VCFFileReader reader = new VCFFileReader(new File(initVCF), false);
		String[] annos = VCFOps.getAnnotationKeys(initVCF, new Logger())[0];
		StringJoiner stats = new StringJoiner("\n");
		stats.add(ArrayUtils.toStr(BASIC_INFO)
				+ "\tMONOMORPHIC\tHOM_DIFF\tPhaseInfo\tALL_GQ\tSAMPS_INTEREST_MONOMORPHIC\t"
				+ ArrayUtils.toStr(getSampHeader(sampsToScan)) + "\t" + ArrayUtils.toStr(annos));

		for (VariantContext vc : reader) {
			if (seg.overlaps(VCOps.getSegment(vc))) {
				VariantContext vcSub = VCOps.getSubset(vc, sampsToScan, VC_SUBSET_TYPE.SUBSET_STRICT);

				if (vcSub.getNoCallCount() == 0) {
					StringJoiner sampStats = new StringJoiner("\t");
					boolean homRef = false;
					boolean homAlt = false;
					boolean phaseInfo = false;
					boolean allGQ = true;
					for (String samp : sampsToScan) {
						Genotype g = vcSub.getGenotype(samp);
						if (g.getGQ() < gq) {
							allGQ = false;
						}
						if (g.isHomRef() && sampsInf.contains(samp)) {
							homRef = true;

						}
						if (g.isHomVar() && sampsInf.contains(samp)) {
							homAlt = true;
						}
						SampInf info = getSampInfo(g);
						for (String inf : info.infs) {
							sampStats.add(inf);
						}
						if (info.hasPID || info.hasPS) {
							phaseInfo = true;
						}

					}

					stats.add(
							ArrayUtils.toStr(getBasicInfo(vcSub)) + "\t" + vcSub.isMonomorphicInSamples() + "\t"
									+ (homAlt && homRef) + "\t" + phaseInfo + "\t" + allGQ + "\t"
									+ VCOps.getSubset(vc, sampsInf, VC_SUBSET_TYPE.SUBSET_STRICT)
											.isMonomorphicInSamples()
									+ "\t" + sampStats.toString() + "\t"
									+ ArrayUtils.toStr(VCOps.getAnnotationsFor(annos, vcSub, ".")));

				}

			}
		}
		reader.close();
		Files.write(stats.toString(), output);
	}

}
