package seq.manage;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

import common.AlleleFreq;
import common.Array;
import common.Logger;
import common.Positions;
import filesys.Segment;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

/**
 * Class for common actions on {@link VariantContext} objects
 *
 */
public class VCOps {

	public enum GENOTYPE_INFO {
		GQ("GQ"), AD_REF("AD"), AD_ALT("AD"), DP("DP");
		private String flag;

		private GENOTYPE_INFO(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}
	}

	/**
	 * Fields shared across a variant context
	 *
	 */
	public enum COMMON_INFO {
		/**
		 * Read depth
		 */
		DP("DP"), /**
		 * Mapping quality
		 */
		MQ("MQ");

		/**
		 * flag represented in vcf file
		 */
		private String flag;

		private COMMON_INFO(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}

	}

	/**
	 * @param annosToGet
	 *            String keys to retrieve from the vc's common info
	 * @param vc
	 * @param defaultValue
	 *            if anno is not found, this will be returned
	 * @return
	 */
	public static String[] getAnnotationsFor(String[] annosToGet, VariantContext vc, String defaultValue) {
		String[] annos = new String[annosToGet.length];
		for (int i = 0; i < annos.length; i++) {
			annos[i] = vc.getCommonInfo().getAttributeAsString(annosToGet[i], defaultValue);
		}
		return annos;
	}

	public static double getMAF(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		int[] alleleCounts = getAlleleCounts(vcSub);
		double maf = AlleleFreq.calcMAF(alleleCounts[0], alleleCounts[1], alleleCounts[2]);
		return maf;
	}

	/**
	 * get the minor allele count for a variant context
	 */
	public static double getMAC(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		int[] alleleCounts = getAlleleCounts(vcSub);
		if (Array.sum(alleleCounts) == 0) {
			return Double.NaN;
		} else {
			double minor = Math.min(alleleCounts[0], alleleCounts[2]);
			return (minor * 2 + alleleCounts[1]);
		}
	}

	/**
	 * get the alternate allele count for a variant context
	 */
	public static double getAAC(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		int[] alleleCounts = getAlleleCounts(vcSub);
		return (alleleCounts[2] * 2 + alleleCounts[1]);

	}

	/**
	 * Tests whether the minor allele is the alternate allele
	 */
	public static boolean isMinorAlleleAlternate(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		double mac = getMAC(vcSub, null);
		return getAAC(vcSub, null) == mac;
	}

	/**
	 * @param names
	 *            the names corresponding the the array of jexl expressions
	 * @param expressions
	 *            the expressions to be evaluated
	 * @param log
	 * @return
	 */
	public static List<VariantContextUtils.JexlVCMatchExp> getJexlVCMathExp(String[] names, String[] expressions, Logger log) {
		List<VariantContextUtils.JexlVCMatchExp> jExps = null;
		try {
			jExps = VariantContextUtils.initializeMatchExps(names, expressions);
		} catch (IllegalArgumentException ile) {
			log.reportTimeError("Could not intitialize the jexl expressions:");
			log.reportTimeError("Names: " + Array.toStr(names));
			log.reportTimeError("JEXLs: " + Array.toStr(expressions));
			log.reportException(ile);
		}
		return jExps;
	}

	/**
	 * @param vc
	 * @param jExp
	 *            determine if the expression matches for the {@link VariantContext}
	 * @return
	 */
	public static boolean passesJexl(VariantContext vc, VariantContextUtils.JexlVCMatchExp jExp) {
		return VariantContextUtils.match(vc, jExp);
	}

	/**
	 * @param vc
	 * @param jExps
	 *            determine if all the expressions match for the {@link VariantContext}
	 * @return
	 */
	public static boolean passesJexls(VariantContext vc, List<VariantContextUtils.JexlVCMatchExp> jExps) {
		for (VariantContextUtils.JexlVCMatchExp jExp : jExps) {
			if (!passesJexl(vc, jExp)) {
				return false;
			}
		}
		return true;
	}

	public static Segment getSegment(VariantContext vc) {
		return new Segment(Positions.chromosomeNumber(vc.getChr()), vc.getStart(), vc.getStart());
	}

	/**
	 * Averages an info value across samples
	 */
	public static double getAverageCommonInfo(VariantContext vc, Set<String> sampleNames, COMMON_INFO info) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		double avgCI = Double.NaN;
		if (vcSub.getCommonInfo().hasAttribute(info.getFlag())) {
			int ci = vc.getCommonInfo().getAttributeAsInt(info.getFlag(), 0);
			if (vcSub.getNSamples() > 0) {
				avgCI = (double) ci / vcSub.getNSamples();
			}
		}
		return avgCI;
	}

	/**
	 * @author lane0212
	 * 
	 *         types of alternate contexts available
	 *
	 */
	public enum ALT_ALLELE_CONTEXT_TYPE {
		/**
		 * All genotypes with a one or two non-reference allele(s) will be returned
		 */
		HET_ONLY, /**
		 * All genotypes with a homozygous non-reference allele will be returned
		 */
		HOM_ONLY, /**
		 * All genotypes with a non-reference allele will be returned
		 */
		ALL;
	}

	public static VariantContext getAltAlleleContext(final VariantContext vc, int altAlleleDepth, Logger log) {
		return getAltAlleleContext(vc, altAlleleDepth, ALT_ALLELE_CONTEXT_TYPE.ALL, log);
	}

	/**
	 * @param vc
	 * @param altAlleleDepth
	 *            note that this depth will be applied to either the one non-ref alternate allele, or both non ref alternate alleles if applicable
	 * @param log
	 * @return
	 */
	public static VariantContext getAltAlleleContext(final VariantContext vc, int altAlleleDepth, ALT_ALLELE_CONTEXT_TYPE type, Logger log) {
		GenotypesContext gc = vc.getGenotypes();
		HashSet<String> samplesWithAlt = new HashSet<String>();
		if (altAlleleDepth >= 0) {
			log.reportTimeError("Alt Allele depth filter is currently not in here, JOHN");
			return null;
		}
		for (Genotype geno : gc) {
			boolean use = false;
			switch (type) {
			case ALL:
				use = !geno.isHomRef() && !geno.isNoCall();
				break;
			case HET_ONLY:
				use = geno.isHet();
				break;
			case HOM_ONLY:
				use = geno.isHomVar();
				break;
			default:
				log.reportTimeError("Invalid alt context type " + type);
				break;

			}
			if (use) {
				int[] AD = new int[] { 0, 0 };
				if (altAlleleDepth >= 0) {
					AD = getAppropriateAlleleDepths(vc, geno, log);
				}
				if (altAlleleDepth < 0 || AD[1] > altAlleleDepth) {
					if (altAlleleDepth < 0 || !geno.isHetNonRef() || AD[0] > altAlleleDepth) {// handles the case when both alleles are non-reference
						samplesWithAlt.add(geno.getSampleName());
					}
				}
			}
		}
		return getSubset(vc, samplesWithAlt);
	}

	/**
	 * Will return the average alt allele ratio for heterozygous calls in this context. Note that the ratio is not always representing the ratio to the reference genome (het non ref)<br>
	 * min/max will be used in the het-non ref case
	 */
	public static double getAverageHetAlleleRatio(VariantContext vc, int altAlleleDepth, Logger log) {
		double avg = 0;
		VariantContext vcAlts = getAltAlleleContext(vc, altAlleleDepth, ALT_ALLELE_CONTEXT_TYPE.HET_ONLY, log);
		GenotypesContext gc = vcAlts.getGenotypes();
		for (Genotype g : gc) {
			int[] AD = getAppropriateAlleleDepths(vc, g, log);
			if (g.isHetNonRef()) {
				int min = Math.min(AD[0], AD[1]);
				int max = Math.max(AD[0], AD[1]);
				avg += (double) min / max;
			} else if (AD[0] > 0) {
				avg += (double) AD[1] / AD[0];
			}
		}
		avg = (double) avg / gc.size();
		return avg;
	}

	private static int[] getAppropriateAlleleDepths(VariantContext vc, Genotype g, Logger log) {
		int[] AD = new int[2];
		Arrays.fill(AD, 0);
		List<Allele> gAlleles = g.getAlleles();

		if (gAlleles.size() > 2) {
			log.reportTimeError("Number of alleles should not be greater than 2");
			return null;
		} else if (gAlleles.size() == 0 || !g.hasAD()) {
			List<Allele> varAlleles = vc.getAlleles();
			int[] gAD = g.getAD();
			int indexVar = 0;
			int indexG = 0;
			for (Allele varAllele : varAlleles) {
				for (Allele gAllele : gAlleles) {
					if (gAllele.equals(varAllele)) {
						AD[indexG] = gAD[indexVar];
						indexG++;
					}
				}
				indexVar++;
			}
		}
		return AD;
		// //vc
		// g.ge
		// g.getAlleles();
	}

	public static double getAvgGenotypeInfo(VariantContext vc, Set<String> sampleNames, GENOTYPE_INFO info, Logger log) {
		double avgGI = 0;
		int numWith = 0;
		VariantContext vcSub = getSubset(vc, sampleNames);
		GenotypesContext gc = vcSub.getGenotypes();
		for (Genotype geno : gc) {
			if (geno.hasAnyAttribute(info.getFlag())) {
				numWith++;
				switch (info) {
				case AD_REF:
					avgGI += geno.getAD()[0];
					break;
				case AD_ALT:
					avgGI += geno.getAD()[1];
					break;
				case DP:
					avgGI += geno.getDP();
					break;
				case GQ:
					avgGI += geno.getGQ();
					break;
				default:
					log.reportTimeError("Invalid average type");
					break;
				}
			}
		}
		if (numWith > 0) {
			avgGI = (double) avgGI / numWith;
		} else {
			avgGI = Double.NaN;
		}
		return avgGI;
	}

	public enum VC_SUBSET_TYPE {
		/**
		 * Samples not contained in the vcf will be given missing genotypes
		 */
		SUBSET_LOOSE, /**
		 * A check will be performed and only samples present in the input set and the vcf file will be exported
		 */
		SUBSET_STRICT, /**
		 * Original variant context will be returned
		 */
		NO_SUBSET;

	}

	public static <E> Set<E> getOverlap(Set<E> subset, Set<E> superSet) {
		Hashtable<E, E> overlap = new Hashtable<E, E>();
		// Set<E> overlap = new Set<E>();
		for (E e : superSet) {
			if (subset.contains(e)) {
				overlap.put(e, e);
			}
		}
		return overlap.keySet();
	}

	public static Genotype getGenotypeFor(final VariantContext vc, final String sampleName, VC_SUBSET_TYPE type) {
		return getSubset(vc, sampleName, type).getGenotype(0);
	}

	/**
	 * Subsets to particular samples
	 */
	public static VariantContext getSubset(final VariantContext vc, final String sampleName, VC_SUBSET_TYPE type) {
		HashSet<String> tmp = new HashSet<String>();
		tmp.add(sampleName);
		return getSubset(vc, tmp, type);
	}

	/**
	 * Subsets to particular samples
	 */
	public static VariantContext getSubset(final VariantContext vc, final Set<String> sampleNames) {
		return getSubset(vc, sampleNames, sampleNames == null ? VC_SUBSET_TYPE.NO_SUBSET : VC_SUBSET_TYPE.SUBSET_STRICT);
	}

	/**
	 * Subsets to particular samples
	 */
	public static VariantContext getSubset(final VariantContext vc, final Set<String> sampleNames, VC_SUBSET_TYPE type) {
		VariantContext vcSub = null;
		switch (type) {
		case SUBSET_LOOSE:
			vcSub = vc.subContextFromSamples(sampleNames);
			break;
		case SUBSET_STRICT:
			vcSub = vc.subContextFromSamples(getOverlap(sampleNames, vc.getSampleNames()));
			break;
		case NO_SUBSET:
			vcSub = vc;
			break;
		default:
			break;

		}
		return vcSub;
	}

	public static boolean isInTheseSegments(VariantContext vc, Segment[] orderedSegs) {
		return getOverlappingSegments(vc, orderedSegs) != null;
	}

	public static int[] getOverlappingSegments(VariantContext vc, Segment[] orderedSegs) {
		Segment vcSeg = new Segment(Positions.chromosomeNumber(vc.getChr()), vc.getStart(), vc.getStart());
		int[] indices = Segment.binarySearchForAllOverLappingIndices(vcSeg, orderedSegs);
		return indices;
	}

	public static double getHWE(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		int[] alleleCounts = getAlleleCounts(vcSub);
		double hwe = AlleleFreq.HWE(alleleCounts[0], alleleCounts[1], alleleCounts[2]);
		return hwe;
	}

	public static int[] getAlleleCounts(VariantContext vc) {
		if ((vc.getHomRefCount() + vc.getHetCount() + vc.getHomVarCount() + vc.getNoCallCount()) != vc.getNSamples()) {
			System.err.println("Un accounted for genotypes...");
		}
		return new int[] { vc.getHomRefCount(), vc.getHetCount(), vc.getHomVarCount() };
	}

	public static double getCallRate(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		double callRate = Double.NaN;
		int noCalls = vcSub.getNoCallCount();
		int numSamps = vcSub.getNSamples();
		if (numSamps > 0) {
			int called = numSamps - noCalls;
			callRate = (double) called / numSamps;
		}
		return callRate;
	}

	public static boolean isBiallelic(VariantContext vc) {
		return vc.isBiallelic();
	}

	public static boolean isAmbiguous(VariantContext vc) {
		boolean ambiguous = false;
		String ref = vc.getReference().getDisplayString();
		List<Allele> alt = vc.getAlternateAlleles();
		String[][] ambiguousDefs = new String[][] { { "A", "T" }, { "T", "A" }, { "C", "G" }, { "G", "C" } };
		for (Allele b : alt) {
			for (int i = 0; i < ambiguousDefs.length; i++) {
				if (ref.equals(ambiguousDefs[i][0]) && b.getDisplayString().equals(ambiguousDefs[i][1])) {
					ambiguous = true;
					break;
				}
			}
		}
		return ambiguous;
	}

	public static class LocusID {
		private byte chr;
		private int start;
		private String ref;
		private String alt;
		private String id;

		public LocusID(VariantContext vc) {
			this.chr = Positions.chromosomeNumber(vc.getChr());
			this.start = vc.getStart();
			this.ref = vc.getReference().getDisplayString();
			this.alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();
			this.id = vc.getID().equals(".") ? chr + ":" + start + ":" + ref + ":" + alt : vc.getID();
		}

		public byte getChr() {
			return chr;
		}

		public int getStart() {
			return start;
		}

		public String getRef() {
			return ref;
		}

		public String getAlt() {
			return alt;
		}

		public String getId() {
			return id;
		}

	}

}