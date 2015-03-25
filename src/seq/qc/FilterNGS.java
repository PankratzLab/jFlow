package seq.qc;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import common.Logger;
import seq.manage.VCOps;
import seq.manage.VCOps.GENOTYPE_INFO;

public class FilterNGS implements Serializable {
	private static final long serialVersionUID = 1L;
	private double mappingQualityFilter;
	private double phreadScoreFilter;
	private int[] readDepthFilter;

	public FilterNGS() {

	}

	public FilterNGS(double mappingQualityFilter, double phreadScoreFilter, int[] readDepthFilter) {
		this.mappingQualityFilter = mappingQualityFilter;
		this.phreadScoreFilter = phreadScoreFilter;
		this.readDepthFilter = readDepthFilter;
	}

	public double getMappingQualityFilter() {
		return this.mappingQualityFilter;
	}

	public void setMappingQualityFilter(double mappingQualityFilter) {
		this.mappingQualityFilter = mappingQualityFilter;
	}

	public double getPhreadScoreFilter() {
		return this.phreadScoreFilter;
	}

	public void setPhreadScoreFilter(double phreadScoreFilter) {
		this.phreadScoreFilter = phreadScoreFilter;
	}

	public boolean passesPhreadScore(double phreadScore) {
		return phreadScore >= this.phreadScoreFilter;
	}

	public boolean passesMapQScore(double mapQScore) {
		return mapQScore >= this.mappingQualityFilter;
	}

	public int[] getReadDepthFilter() {
		return this.readDepthFilter;
	}

	public void setReadDepthFilter(int[] readDepthFilter) {
		this.readDepthFilter = readDepthFilter;
	}

	public enum FILTER_TYPE {
		/**
		 * Will always pass filter
		 */
		NO_FILTER, /**
		 * Check if greater than
		 */
		GT_FILTER, /**
		 * Check if less than
		 */
		LT_FILTER, /**
		 * Check if greater than or equal to
		 */
		GTE_FILTER, /**
		 * Check if less than or equal to
		 */
		LTE, /**
		 * Check if equal to
		 */
		ET_FILTER, /**
		 * Check if true
		 */
		TRUE_BOOL, /**
		 * Check if false
		 */
		FALSE_BOOL;
	}


	
	public enum VARIANT_FILTER_DOUBLE {

		MAC(2.0, FILTER_TYPE.GTE_FILTER),

		GQ(90, FILTER_TYPE.GTE_FILTER), DP(10, FILTER_TYPE.GTE_FILTER),
		/**
		 * Use for minor allele frequency filtering, maf must be greater
		 */
		MAF(0.05, FILTER_TYPE.GTE_FILTER), /**
		 * Hardy wienberg filtering
		 */
		HWE(0.0001, FILTER_TYPE.GTE_FILTER), /**
		 * Call rate Filtering
		 */
		CALL_RATE(0.95, FILTER_TYPE.GTE_FILTER), /**
		 * Call rate Filtering
		 */
		CALL_RATE_LOOSE(0.90, FILTER_TYPE.GTE_FILTER);

		private double dFilter;
		private FILTER_TYPE type;

		public static VARIANT_FILTER_DOUBLE[] getFiltersExcluding(VARIANT_FILTER_DOUBLE[] exclude) {
			ArrayList<VARIANT_FILTER_DOUBLE> remain = new ArrayList<VARIANT_FILTER_DOUBLE>();
			for (int i = 0; i < VARIANT_FILTER_DOUBLE.values().length; i++) {
				boolean use = true;
				for (int j = 0; j < exclude.length; j++) {
					if (VARIANT_FILTER_DOUBLE.values()[i] == exclude[j]) {
						use = false;
					}
				}
				if (use) {
					remain.add(VARIANT_FILTER_DOUBLE.values()[i]);
				}
			}
			return remain.toArray(new VARIANT_FILTER_DOUBLE[remain.size()]);
		}

		private VARIANT_FILTER_DOUBLE(double defaultD, FILTER_TYPE type) {
			this.dFilter = defaultD;
			this.type = type;
		}

		public double getDFilter() {
			return dFilter;
		}

		public void setDFilter(double dFilter) {
			this.dFilter = dFilter;
		}

		public FILTER_TYPE getType() {
			return type;
		}

		public void setType(FILTER_TYPE type) {
			this.type = type;
		}

	}

	public enum VARIANT_FILTER_BOOLEAN {
		/**
		 * will pass if the variant is biallelic
		 */
		BIALLELIC_FILTER(FILTER_TYPE.TRUE_BOOL), /**
		 * Will pass if the variant is un-ambiguous
		 */
		AMBIGUOUS_FILTER(FILTER_TYPE.FALSE_BOOL), /**
		 * Is true when the jexl formatted expression is evaluated
		 */
		JEXL(FILTER_TYPE.TRUE_BOOL);

		private FILTER_TYPE type;

		private VARIANT_FILTER_BOOLEAN(FILTER_TYPE type) {
			this.type = type;
		}

		public FILTER_TYPE getType() {
			return type;
		}

		public void setType(FILTER_TYPE type) {
			this.type = type;
		}

	}

	/**
	 * 
	 * Handles common variant filtering options
	 *
	 */
	private interface VcFilterI<T> {
		public boolean shouldFilter();

		public void setFilterType(FILTER_TYPE type);

		public VariantContextFilterPass filter(VariantContext vc, Logger log);

		public T getValue(VariantContext vc);

	}

	/**
	 * Handles command variant filtering when dealing with double data
	 *
	 */
	private abstract class VcFilterDouble implements VcFilterI<Double> {
		private VARIANT_FILTER_DOUBLE dfilter;
		private FILTER_TYPE type;
		private double filterThreshold;

		public VcFilterDouble(VARIANT_FILTER_DOUBLE dfilter) {
			super();
			this.dfilter = dfilter;
			this.type = dfilter.getType();
			this.filterThreshold = dfilter.getDFilter();
		}

		@Override
		public boolean shouldFilter() {
			return type != FILTER_TYPE.NO_FILTER;
		}

		public VARIANT_FILTER_DOUBLE getDfilter() {
			return dfilter;
		}

		@Override
		public void setFilterType(FILTER_TYPE type) {
			this.type = type;
		}

		@Override
		public VariantContextFilterPass filter(VariantContext vc, Logger log) {
			double value = getValue(vc);
			String testPerformed = "Type: " + dfilter + " Directon: " + type + " Threshold " + filterThreshold;
			boolean passes = false;
			switch (type) {
			case ET_FILTER:
				passes = value == filterThreshold;
				break;
			case GTE_FILTER:
				passes = value >= filterThreshold;
				break;
			case GT_FILTER:
				passes = value > filterThreshold;
				break;
			case LTE:
				passes = value <= filterThreshold;
				break;
			case LT_FILTER:
				passes = value < filterThreshold;
				break;
			case NO_FILTER:
				passes = true;
				break;
			default:
				passes = true;
				break;

			}

			return new VariantContextFilterPass(passes, testPerformed);
		}

		@Override
		public Double getValue(VariantContext vc) {
			return Double.NaN;
		}

	}

	/**
	 * Handles Filtering with {@link VariantContextUtils.JexlVCMatchExp}
	 *
	 */
	private class VcFilterJEXL extends VcFilterBoolean {
		private List<VariantContextUtils.JexlVCMatchExp> jExps;

		public VcFilterJEXL(VARIANT_FILTER_BOOLEAN bfilter, String[] names, String[] expressions, Logger log) {
			super(bfilter);
			this.jExps = VCOps.getJexlVCMathExp(names, expressions, log);
		}

		@Override
		public Boolean getValue(VariantContext vc) {
			return VCOps.passesJexls(vc, jExps);
		}
	}

	/**
	 * Handles boolean based filtering
	 *
	 */
	private abstract class VcFilterBoolean implements VcFilterI<Boolean> {
		private VARIANT_FILTER_BOOLEAN bfilter;
		private FILTER_TYPE type;

		public VcFilterBoolean(VARIANT_FILTER_BOOLEAN bfilter) {
			super();
			this.bfilter = bfilter;
			this.type = bfilter.getType();
		}

		public VARIANT_FILTER_BOOLEAN getBfilter() {
			return bfilter;
		}

		@Override
		public boolean shouldFilter() {
			return type != FILTER_TYPE.NO_FILTER;
		}

		@Override
		public void setFilterType(FILTER_TYPE type) {
			this.type = type;
		}

		@Override
		public VariantContextFilterPass filter(VariantContext vc, Logger log) {

			boolean value = getValue(vc);
			boolean passes = false;
			String testPerformed = "Type: " + bfilter + " Directon: " + type;

			switch (type) {
			case TRUE_BOOL:
				passes = value;
				break;
			case FALSE_BOOL:
				passes = !value;
				break;
			default:
				passes = true;
				break;
			}
			return new VariantContextFilterPass(passes, testPerformed);
		}

		@Override
		public Boolean getValue(VariantContext vc) {
			return true;
		}

	}

	private VcFilterDouble getMAFFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getMAF(vc, null);
			}
		};
	}

	private VcFilterDouble getHWEFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getHWE(vc, null);
			}
		};
	}

	private VcFilterDouble getCallRateFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getCallRate(vc, null);
			}
		};
	}

	private VcFilterDouble getMACFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getMAC(vc, null);
			}
		};
	}

	private VcFilterDouble getAvgGQFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.GQ, log);
			}
		};
	}

	private VcFilterDouble getAvgDPFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.DP, log);
			}
		};
	}

	private VcFilterBoolean getBiallelicFilter(VARIANT_FILTER_BOOLEAN bfilter) {
		return new VcFilterBoolean(bfilter) {
			@Override
			public Boolean getValue(VariantContext vc) {
				return VCOps.isBiallelic(vc);
			}
		};
	}

	private VcFilterBoolean getUnambiguousFilter(VARIANT_FILTER_BOOLEAN bfilter) {
		return new VcFilterBoolean(bfilter) {
			@Override
			public Boolean getValue(VariantContext vc) {
				return VCOps.isAmbiguous(vc);
			}
		};
	}

	private VcFilterJEXL getJEXLFilter(VARIANT_FILTER_BOOLEAN bfilter, String[] names, String[] expressions, Logger log) {
		return new VcFilterJEXL(bfilter, names, expressions, log);
	}

	private VcFilterBoolean[] getBooleanFilters(VARIANT_FILTER_BOOLEAN[] bFilters, Logger log) {
		VcFilterBoolean[] vBooleans = new VcFilterBoolean[bFilters.length];
		for (int i = 0; i < vBooleans.length; i++) {
			VARIANT_FILTER_BOOLEAN bfilter = bFilters[i];
			switch (bfilter) {
			case AMBIGUOUS_FILTER:
				vBooleans[i] = getUnambiguousFilter(bfilter);
				break;
			case BIALLELIC_FILTER:
				vBooleans[i] = getBiallelicFilter(bfilter);
				break;
			default:
				log.reportTimeError("Invalid boolean filter type " + bfilter);
				vBooleans[i] = null;
				break;
			}
		}
		return vBooleans;
	}

	private VcFilterDouble[] getDoubleFilters(VARIANT_FILTER_DOUBLE[] dFilters, Logger log) {
		VcFilterDouble[] vDoubles = new VcFilterDouble[dFilters.length];
		for (int i = 0; i < vDoubles.length; i++) {
			VARIANT_FILTER_DOUBLE dfilter = dFilters[i];
			switch (dfilter) {
			case MAC:
				vDoubles[i] = getMACFilter(dfilter);
				break;
			case GQ:
				vDoubles[i] = getAvgGQFilter(dfilter, log);
				break;
			case DP:
				vDoubles[i] = getAvgDPFilter(dfilter, log);
				break;
			case CALL_RATE:
				vDoubles[i] = getCallRateFilter(dfilter);
				break;
			case CALL_RATE_LOOSE:
				vDoubles[i] = getCallRateFilter(dfilter);
				break;
			case HWE:
				vDoubles[i] = getHWEFilter(dfilter);
				break;
			case MAF:
				vDoubles[i] = getMAFFilter(dfilter);
				break;
			default:
				log.reportTimeError("Invalid double filter type " + dfilter);
				vDoubles[i] = null;
				break;
			}
		}
		return vDoubles;
	}

	/**
	 * Combines the type specific filters
	 *
	 */
	public static class VariantContextFilter {
		private VcFilterDouble[] vDoubles;
		private VcFilterBoolean[] vBooleans;
		private VcFilterJEXL vFilterJEXL;
		private Logger log;

		public VariantContextFilter(VARIANT_FILTER_DOUBLE[] dFilters, VARIANT_FILTER_BOOLEAN[] bFilters, String[] jexlNames, String[] jexlExpression, Logger log) {
			this.log = log;
			this.vDoubles = new FilterNGS().getDoubleFilters(dFilters, log);
			this.vBooleans = new FilterNGS().getBooleanFilters(bFilters, log);
			if (jexlExpression != null) {
				vFilterJEXL = new FilterNGS().getJEXLFilter(VARIANT_FILTER_BOOLEAN.JEXL, jexlNames, jexlExpression, log);
			} else {
				vFilterJEXL = null;
			}
		}

		/**
		 * @param vc
		 * @return if the variant passed all filters Note that this method returns after the first filter it fails for speed
		 */
		public VariantContextFilterPass filter(VariantContext vc) {

			if (vFilterJEXL != null) {
				VariantContextFilterPass vcfp = vFilterJEXL.filter(vc, log);
				if (!vcfp.passed()) {
					return vcfp;
				}
			}
			for (int i = 0; i < vBooleans.length; i++) {
				VariantContextFilterPass vcfp = vBooleans[i].filter(vc, log);
				if (!vcfp.passed()) {
					// if(vBooleans[i].getBfilter()==VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER){
					// System.out.println(VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER);
					// System.exit(1);
					// }
					return vcfp;
				}
			}
			for (int i = 0; i < vDoubles.length; i++) {
				VariantContextFilterPass vcfp = vDoubles[i].filter(vc, log);

				if (!vcfp.passed()) {
					// if (vDoubles[i].getDfilter() == VARIANT_FILTER_DOUBLE.MAF) {
					// System.out.println(VARIANT_FILTER_DOUBLE.MAF);
					// System.exit(1);
					// }
					return vcfp;
				}
			}
			return new VariantContextFilterPass(true, "ALL");
		}
	}

	public static class VariantContextFilterPass {

		private boolean passed;
		private String testPerformed;

		public VariantContextFilterPass(boolean passed, String testPerformed) {
			super();
			this.passed = passed;
			this.testPerformed = testPerformed;
		}

		public boolean passed() {
			return passed;
		}

		public String getTestPerformed() {
			return passed ? "Passed :" : "Failed " + testPerformed;

		}

	}

	public static class RareVariantFilter {
		private Set<String> totalPop;
		private Set<String> casePop;
		private Set<String> refPop;
		private VARIANT_FILTER_DOUBLE mafRef = VARIANT_FILTER_DOUBLE.MAF;

		private VARIANT_FILTER_DOUBLE macCase = VARIANT_FILTER_DOUBLE.MAC;
		private VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		private VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
		private VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		private VariantContextFilter refFilters;
		private VariantContextFilter caseFilters;

		public RareVariantFilter(Set<String> casePop, Set<String> refPop) {
			super();
			this.casePop = casePop;
			this.refPop = refPop;
			mafRef.setType(FILTER_TYPE.LTE);
			this.totalPop = new HashSet<String>();
		}

		public void setMafRef(double maf) {
			mafRef.setDFilter(maf);
		}

		public void setMacCase(double mac) {
			macCase.setDFilter(mac);
		}

		public void setCallRate(double cr) {
			callRate.setDFilter(cr);
		}

		public void setDepth(double depth) {
			dp.setDFilter(depth);
		}

		public void setGQ(double genoQual) {
			gq.setDFilter(genoQual);
		}

		public void initFilters(Logger log) {
			VARIANT_FILTER_DOUBLE[] refF = new VARIANT_FILTER_DOUBLE[] { callRate, dp, gq };
			VARIANT_FILTER_DOUBLE[] caseF = new VARIANT_FILTER_DOUBLE[] { callRate, dp, gq };
			// VARIANT_FILTER_DOUBLE[] refF = new VARIANT_FILTER_DOUBLE[] { mafRef, callRate, dp, gq };
			// VARIANT_FILTER_DOUBLE[] caseF = new VARIANT_FILTER_DOUBLE[] { macCase, callRate, dp, gq };
			this.refFilters = new VariantContextFilter(refF, new VARIANT_FILTER_BOOLEAN[] {}, null, null, log);
			this.caseFilters = new VariantContextFilter(caseF, new VARIANT_FILTER_BOOLEAN[] {}, null, null, log);
		}

		public VariantContextFilterPass filter(final VariantContext vc, Logger log) {
			VariantContext vcRef = VCOps.getSubset(vc, refPop);
			VariantContextFilterPass pass = refFilters.filter(vcRef);
			if (pass.passed()) {
				VariantContext vcCase = VCOps.getSubset(vc, casePop);
				pass = caseFilters.filter(vcCase);
				if (pass.passed()) {
					VariantContext vcAlts = VCOps.getAltAlleleContext(vcCase);
					pass = caseFilters.filter(vcAlts);
					if (pass.passed()) {
					}
				}
			}
			return pass;
		}

	}

}
