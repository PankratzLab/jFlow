package seq.manage;

import java.util.Hashtable;
import java.util.List;
import java.util.Set;

import common.AlleleFreq;
import common.Positions;
import filesys.Segment;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Class for common actions on {@link VariantContext} objects
 *
 */
public class VCOps {

	public enum INFO {
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

		private INFO(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}

	}

	public static double getMAF(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		int[] alleleCounts = getAlleleCounts(vcSub);
		double maf = AlleleFreq.calcMAF(alleleCounts[0], alleleCounts[1], alleleCounts[2]);
		return maf;
	}

	/**
	 * Averages an info value across samples
	 */
	public static double getAverageCommonInfo(VariantContext vc, Set<String> sampleNames, INFO info) {
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
		vcSub = sampleNames == null || sampleNames.size() < 1 ? vc : vc.subContextFromSamples(sampleNames);
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
		return new int[] { vc.getHomRefCount(), vc.getHetCount(), vc.getHomVarCount() };
	}

	public static double getCallRate(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = getSubset(vc, sampleNames);
		double callRate = Double.NaN;
		int noCalls = vcSub.getNoCallCount();
		int numSamps = vcSub.getNSamples();
		if (numSamps > 0) {
			callRate = (double) ((numSamps - noCalls) / numSamps);
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

}
