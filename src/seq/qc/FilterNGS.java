package seq.qc;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import common.Array;
import common.Logger;
import seq.manage.VCOps;
import seq.manage.VCOps.GENOTYPE_INFO;

public class FilterNGS implements Serializable {
	private static final long serialVersionUID = 1L;
	private double mappingQualityFilter;
	private double phreadScoreFilter;
	private int[] readDepthFilter;
	private int[] altAlleleDepthFilter;
	private double[] altAlleleDepthRatioFilter;

	public FilterNGS() {

	}

	public FilterNGS(double mappingQualityFilter, double phreadScoreFilter, int[] readDepthFilter) {
		this.mappingQualityFilter = mappingQualityFilter;
		this.phreadScoreFilter = phreadScoreFilter;
		this.readDepthFilter = readDepthFilter;
	}

	public int[] getAltAlleleDepthFilter() {
		return altAlleleDepthFilter;
	}

	public double[] getAltAlleleDepthRatioFilter() {
		return altAlleleDepthRatioFilter;
	}

	public void setAltAlleleDepthRatioFilter(double[] altAlleleDepthRatioFilter) {
		this.altAlleleDepthRatioFilter = altAlleleDepthRatioFilter;
	}

	public void setAltAlleleDepthFilter(int[] altAlleleDepthFilter) {
		this.altAlleleDepthFilter = altAlleleDepthFilter;
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
		LTE_FILTER, /**
		 * Check if equal to
		 */
		ET_FILTER, /**
		 * Check if true
		 */
		TRUE_BOOL, /**
		 * Check if false
		 */
		FALSE_BOOL,

		;
	}

	public enum VARIANT_FILTER_DOUBLE {

		MAC(2.0, FILTER_TYPE.GTE_FILTER),

		/**
		 * Average GQ across a variant
		 */
		GQ_STRICT(90, FILTER_TYPE.GTE_FILTER),
		/**
		 * Average GQ across a variant
		 */
		GQ(50, FILTER_TYPE.GTE_FILTER),

		/**
		 * Average Depth across a variant
		 */
		DP(10, FILTER_TYPE.GTE_FILTER),
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
		CALL_RATE_LOOSE(0.90, FILTER_TYPE.GTE_FILTER),
		/**
		 * Used for strict filtering on VQSLOD score of variant
		 */
		VQSLOD_STRICT(3.0, FILTER_TYPE.GTE_FILTER),
		/**
		 * Used for loose filtering on VQSLOD score of variant
		 */
		VQSLOD_LOOSE(1.0, FILTER_TYPE.GTE_FILTER),
		/**
		 * Used for filtering by the lower bound for average allelic ratio across alternate calls
		 */
		HET_ALLELE_RATIO_LOW(.25, FILTER_TYPE.GTE_FILTER),

		/**
		 * Used for filtering by the lower bound for average allelic ratio across alternate calls
		 */
		HET_ALLELE_RATIO_HIGH(.75, FILTER_TYPE.LTE_FILTER),

		/**
		 * Filters by the alternate allele depth
		 */
		ALT_ALLELE_DEPTH(5, FILTER_TYPE.GTE_FILTER),

		/**
		 * Specific filters for tumor normal comparisons
		 */
		/**
		 * Total depth for the tumor
		 */
		AD_TUMOR(10, FILTER_TYPE.GTE_FILTER),

		/**
		 * Total depth for the normal
		 */
		AD_NORMAL(10, FILTER_TYPE.GTE_FILTER),

		/**
		 * Alt allele depth for tumor
		 */
		ALT_AD_TUMOR(4, FILTER_TYPE.GTE_FILTER),

		/**
		 * Alt allele depth for normal
		 */
		ALT_AD_NORMAL(0, FILTER_TYPE.LTE_FILTER),

		/**
		 * Allelic fraction for the tumor
		 */
		AF_TUMOR(0.2, FILTER_TYPE.GTE_FILTER),

		/**
		 * https://www.broadinstitute.org/cancer/cga/mutect for description of why 6.3 and 2.3
		 * 
		 */
		/**
		 * Set so errors are half the somatic mutation rate
		 */
		TLOD(6.3, FILTER_TYPE.GTE_FILTER),

		NLOD(2.3, FILTER_TYPE.GTE_FILTER)

		;

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

	// public enum GENOTYPE_FAIL_FILTER {
	//
	// MUTECT_FAIL_CONTAINS(FILTER_TYPE.TRUE_BOOL, new String[] { "PASS" }, GENOTYPE_INFO.MUTECT_FILTERS);
	//
	// private FILTER_TYPE type;
	// private String[] pass;
	// private String[] exceptions;
	// private GENOTYPE_INFO gInfo;
	//
	// private GENOTYPE_FAIL_FILTER(FILTER_TYPE type, String[] ops, GENOTYPE_INFO gInfo) {
	// this.type = type;
	// }
	//
	// public void setType(FILTER_TYPE type) {
	// this.type = type;
	// }
	//
	// public void setOps(String[] ops) {
	// this.pass = ops;
	// }
	// }

	public enum VARIANT_FILTER_BOOLEAN {
		/**
		 * will pass if the variant is biallelic
		 */
		BIALLELIC_FILTER(FILTER_TYPE.TRUE_BOOL), /**
		 * Will pass if the variant is un-ambiguous
		 */
		AMBIGUOUS_FILTER(FILTER_TYPE.FALSE_BOOL), /**
		 * Is true when the variant has not been filtered out previously
		 */
		FAILURE_FILTER(FILTER_TYPE.TRUE_BOOL), /**
		 * Is true when the jexl formatted expression is evaluated
		 */
		JEXL(FILTER_TYPE.TRUE_BOOL),
		/**
		 * When a genotype has a custom annotation denoting a mutect failure
		 */
		MUTECT_FAIL_FILTER(FILTER_TYPE.TRUE_BOOL);
		;

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
	public abstract class VcFilterDouble implements VcFilterI<Double>, Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private VARIANT_FILTER_DOUBLE dfilter;
		private FILTER_TYPE type;
		private double filterThreshold;

		public VcFilterDouble(VARIANT_FILTER_DOUBLE dfilter) {
			super();
			this.dfilter = dfilter;
			this.type = dfilter.getType();
			this.filterThreshold = dfilter.getDFilter();
		}

		public double getFilterThreshold() {
			return filterThreshold;
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
			// " Value :" + value + (vc.getSampleNames().size() == 1 ? vc.getGenotype(0).toString() : ""
			String testPerformed = "Type: " + dfilter + " Directon: " + type + " Threshold: " + filterThreshold + " Value :" + value;
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
			case LTE_FILTER:
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
	 *
	 */
	public class VcFilterString extends VcFilterBoolean implements Serializable {

		public VcFilterString(VARIANT_FILTER_BOOLEAN bfilter) {
			super(bfilter);
			// TODO Auto-generated constructor stub
		}

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

	}

	/**
	 * Handles Filtering with {@link VariantContextUtils.JexlVCMatchExp}
	 *
	 */
	public class VcFilterJEXL extends VcFilterBoolean implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private String[] names, expressions;
		private Logger log;

		// private List<VariantContextUtils.JexlVCMatchExp> jExps;

		public VcFilterJEXL(VARIANT_FILTER_BOOLEAN bfilter, String[] names, String[] expressions, Logger log) {
			super(bfilter);
			this.names = names;
			this.expressions = expressions;
		}

		public List<VariantContextUtils.JexlVCMatchExp> getjExps() {
			return VCOps.getJexlVCMathExp(names, expressions, log);
		}

		@Override
		public Boolean getValue(VariantContext vc) {
			return VCOps.passesJexls(vc, VCOps.getJexlVCMathExp(names, expressions, log));
		}
	}

	/**
	 * Handles boolean based filtering
	 *
	 */
	public abstract class VcFilterBoolean implements VcFilterI<Boolean>, Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
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
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getMAF(vc, null);
			}
		};
	}

	private VcFilterDouble getHWEFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				double hwe = VCOps.getHWE(vc, null);
				if (Double.isNaN(hwe)) {
					if (vc.getHomRefCount() + vc.getNoCallCount() == vc.getSampleNames().size()) {
						hwe = 1.0;
					} else if (vc.getHomVarCount() + vc.getNoCallCount() == vc.getSampleNames().size()) {
						hwe = 1.0;
					} else {
						String error = "HWE came back NaN but the variant context had valid data\n";
						error += "Total Samples: " + vc.getSampleNames().size();
						error += "\nHomRef: " + vc.getHomRefCount();
						error += "\nHet : " + vc.getHetCount();
						error += "\nHomVar : " + vc.getHomVarCount();
						error += "\nNoCall : " + vc.getNoCallCount();
						error += "\n" + vc.toStringWithoutGenotypes();
						throw new IllegalStateException(error);
					}
				}
				return hwe;
			}
		};
	}

	private VcFilterDouble getCallRateFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getCallRate(vc, null);
			}
		};
	}

	/**
	 * This will ONLY test Het Calls, use {@link FilterNGS#getAltAlleleDepthFilter()} for hom alt etc
	 */
	private VcFilterDouble getHetAltAlleleDepthRatioFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				if (vc.getSampleNames().size() > 1) {
					throw new IllegalArgumentException("Alt allele depth filter can only be applied on single sample variant contexts");
				} else if (vc.getGenotype(0).isHet()) {
					int[] AD = VCOps.getAppropriateAlleleDepths(vc, vc.getGenotype(0), false, new Logger());
					if (AD[0] == 0) {
						return 1.0;
					} else {
						double sum = Array.sum(AD);
						return (double) AD[1] / sum;
					}
				} else {// we only test het calls, and return a passing value if it is not
					switch (getDfilter().getType()) {
					case GTE_FILTER:
						return getDfilter().getDFilter() + .01;
					case GT_FILTER:
						return getDfilter().getDFilter() + .01;
					case LTE_FILTER:
						return getDfilter().getDFilter() - .01;
					case LT_FILTER:
						return getDfilter().getDFilter() - .01;
					case NO_FILTER:
						return 1.0;
					default:
						throw new IllegalArgumentException("Invalid Type filter " + getDfilter().getType() + " for double filter");

					}
				}
			}
		};
	}

	private VcFilterDouble getAltAlleleDepthFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				if (vc.getSampleNames().size() > 1) {
					throw new IllegalArgumentException("Alt allele depth filter can only be applied on single sample variant contexts");
				} else {
					double apad = 0;
					try {
						apad = (double) VCOps.getAppropriateAlleleDepths(vc, vc.getGenotype(0), true, new Logger())[1];
					} catch (IllegalStateException illegalStateException) {
						illegalStateException.printStackTrace();
					}
					return apad;
				}
			}
		};
	}

	private VcFilterDouble getVQSLODFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				if (!vc.hasAttribute("VQSLOD")) {
					return Double.NaN;
				} else {
					return vc.getCommonInfo().getAttributeAsDouble("VQSLOD", 0.0);
				}
			}
		};
	}

	// private VcFilterDouble getHetAlleleRatioFilter(VARIANT_FILTER_DOUBLE dfilter){
	// return new VcFilterDouble(dfilter) {
	// @Override
	// public Double getValue(VariantContext vc) {
	// if(VCOps.getAltAlleleContext(vc, altAlleleDepth, type, log))
	// }
	// };
	// }

	private VcFilterDouble getMACFilter(VARIANT_FILTER_DOUBLE dfilter) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getMAC(vc, null);
			}
		};
	}

	private VcFilterDouble getAvgGQFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.GQ, log);
			}
		};
	}

	// case AD_MUT:
	// vDoubles[i] =
	// break;
	// case AD_NORMAL:
	// vDoubles[i] =
	// break;
	// case ALT_AD_MUT:
	// vDoubles[i] =
	// break;
	// case AF_MUT:
	private VcFilterDouble getAvgAD_TUMORFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.AD_TUMOR, log);
			}
		};
	}

	private VcFilterDouble getAvgAD_NORMALFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.AD_NORMAL, log);
			}
		};
	}

	private VcFilterDouble getAvgALT_AD_TUMORFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.ALT_AD_TUMOR, log);
			}
		};
	}

	private VcFilterDouble getAvgALT_AD_NORMALFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.ALT_AD_NORMAL, log);
			}
		};
	}

	private VcFilterDouble getAvgAF_TUMORFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.AF_TUMOR, log);
			}
		};
	}

	private VcFilterDouble getAvgTLODFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.TLOD, log);
			}
		};
	}

	private VcFilterDouble getAvgNLODFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.NLOD, log);
			}
		};
	}

	private VcFilterDouble getAvgDPFilter(VARIANT_FILTER_DOUBLE dfilter, final Logger log) {
		return new VcFilterDouble(dfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Double getValue(VariantContext vc) {
				return VCOps.getAvgGenotypeInfo(vc, null, GENOTYPE_INFO.DP, log);
			}
		};
	}

	private VcFilterBoolean getBiallelicFilter(VARIANT_FILTER_BOOLEAN bfilter) {
		return new VcFilterBoolean(bfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Boolean getValue(VariantContext vc) {
				return VCOps.isBiallelic(vc);
			}
		};
	}

	private VcFilterBoolean getFailureFilter(VARIANT_FILTER_BOOLEAN bfilter) {
		return new VcFilterBoolean(bfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Boolean getValue(VariantContext vc) {
				return vc.isNotFiltered();
			}
		};
	}

	private VcFilterBoolean getMutectFailureFilter(VARIANT_FILTER_BOOLEAN bfilter) {
		return new VcFilterBoolean(bfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Boolean getValue(VariantContext vc) {
				if (vc.getGenotypes().size() != 1) {
					throw new IllegalArgumentException("Mutect failure checks can only be used on a single genotype VariantContex");
				}
				Genotype g = vc.getGenotype(0);
				boolean pass = true;
				if (g.hasAnyAttribute(GENOTYPE_INFO.MUTECT_FILTERS.getFlag())) {
					List<String> filts = Arrays.asList(g.getAnyAttribute(GENOTYPE_INFO.MUTECT_FILTERS.getFlag()).toString());
					pass = filts.size() == 0 || filts.get(0).equals("PASS");
				}

				return pass;
			}
		};
	}

	private VcFilterBoolean getUnambiguousFilter(VARIANT_FILTER_BOOLEAN bfilter) {
		return new VcFilterBoolean(bfilter) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

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
			case FAILURE_FILTER:
				vBooleans[i] = getFailureFilter(bfilter);
				break;
			case MUTECT_FAIL_FILTER:
				vBooleans[i] = getMutectFailureFilter(bfilter);
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
			case GQ_STRICT:
				vDoubles[i] = getAvgGQFilter(dfilter, log);
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
			case VQSLOD_STRICT:
				vDoubles[i] = getVQSLODFilter(dfilter);
				break;
			case VQSLOD_LOOSE:
				vDoubles[i] = getVQSLODFilter(dfilter);
				break;

			case ALT_ALLELE_DEPTH:
				vDoubles[i] = getAltAlleleDepthFilter(dfilter);
				break;
			case HET_ALLELE_RATIO_LOW:
				vDoubles[i] = getHetAltAlleleDepthRatioFilter(dfilter);
				break;
			case HET_ALLELE_RATIO_HIGH:
				vDoubles[i] = getHetAltAlleleDepthRatioFilter(dfilter);
				break;

			case AD_TUMOR:
				vDoubles[i] = getAvgAD_TUMORFilter(dfilter, log);
				break;
			case AD_NORMAL:
				vDoubles[i] = getAvgAD_NORMALFilter(dfilter, log);
				break;
			case ALT_AD_TUMOR:
				vDoubles[i] = getAvgALT_AD_TUMORFilter(dfilter, log);
				break;
			case AF_TUMOR:
				vDoubles[i] = getAvgAF_TUMORFilter(dfilter, log);
				break;

			case ALT_AD_NORMAL:
				vDoubles[i] = getAvgALT_AD_NORMALFilter(dfilter, log);
				break;

			case TLOD:
				vDoubles[i] = getAvgTLODFilter(dfilter, log);
				break;
			case NLOD:
				vDoubles[i] = getAvgNLODFilter(dfilter, log);
				break;
			default:
				log.reportTimeError("Invalid double filter type " + dfilter);
				vDoubles[i] = null;
				throw new IllegalArgumentException("Invalid double filter type " + dfilter);
				// break;
			}
		}
		return vDoubles;
	}

	/**
	 * Combines the type specific filters
	 *
	 */
	public static class VariantContextFilter implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
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

		public Logger getLog() {
			return log;
		}

		public VcFilterDouble[] getvDoubles() {
			return vDoubles;
		}

		public VcFilterBoolean[] getvBooleans() {
			return vBooleans;
		}

		public VcFilterJEXL getvFilterJEXL() {
			return vFilterJEXL;
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
			if (vBooleans != null) {
				for (int i = 0; i < vBooleans.length; i++) {
					VariantContextFilterPass vcfp = vBooleans[i].filter(vc, log);
					if (!vcfp.passed()) {
						return vcfp;
					}
				}
			}
			if (vDoubles != null) {
				for (int i = 0; i < vDoubles.length; i++) {
					VariantContextFilterPass vcfp = vDoubles[i].filter(vc, log);
					if (!vcfp.passed()) {
						return vcfp;
					}
				}
			}
			return new VariantContextFilterPass(true, "ALL");
		}

		public static class VariantContextFilterBuilder {
			private VcFilterDouble[] vDoubles = new VcFilterDouble[] {};
			private VcFilterBoolean[] vBooleans = new VcFilterBoolean[] {};
			private VcFilterJEXL vFilterJEXL = new FilterNGS().getJEXLFilter(VARIANT_FILTER_BOOLEAN.JEXL, new String[] {}, new String[] {}, new Logger());

			public VariantContextFilterBuilder altAlleleDepthRatioFilter(double[] altAlleleDepthRatioFilter) {
				return this;
			}

			public VariantContextFilter build(Logger log) {
				return new VariantContextFilter(this, log);
			}
		}

		private VariantContextFilter(VariantContextFilterBuilder builder, Logger log) {
			this.vDoubles = builder.vDoubles;
			this.vBooleans = builder.vBooleans;
			this.vFilterJEXL = builder.vFilterJEXL;
			this.log = log;
		}
	}

	public static AggregateFilter initializeFilters(FilterNGS filterNGS, SAM_FILTER_TYPE filterType, Logger log) {
		if (filterNGS == null) {
			filterNGS = new FilterNGS();
		}
		ArrayList<SamRecordFilter> filters = filterNGS.getStandardSAMRecordFilters(filterType, log);
		filters.add(filterNGS.getSamRecordMapQFilter(filterNGS.getMappingQualityFilter()));
		AggregateFilter filter = new AggregateFilter(filters);
		return filter;
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
		// private Set<String> totalPop;
		private Set<String> casePop;
		private Set<String> refPop;
		private VARIANT_FILTER_DOUBLE mafRef = VARIANT_FILTER_DOUBLE.MAF;

		private VARIANT_FILTER_DOUBLE macCase = VARIANT_FILTER_DOUBLE.MAC;
		private VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		private VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_STRICT;
		private VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		private VariantContextFilter refFilters;
		private VariantContextFilter caseFilters;

		public RareVariantFilter(Set<String> casePop, Set<String> refPop) {
			super();
			this.casePop = casePop;
			this.refPop = refPop;
			mafRef.setType(FILTER_TYPE.LTE_FILTER);
			// this.totalPop = new HashSet<String>();
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
					VariantContext vcAlts = VCOps.getAltAlleleContext(vcCase, null, log);
					pass = caseFilters.filter(vcAlts);
					if (pass.passed()) {
					}
				}
			}
			return pass;
		}

	}

	// Some Same Record Filters

	public abstract class SamRecordExtFilter implements SamRecordFilter {
		protected double doubleFilter;

		public SamRecordExtFilter(double doubleFilter) {
			super();
			this.doubleFilter = doubleFilter;
		}
	}

	/**
	 * @param mapQ
	 * @return filter that will remove records with mapping qualities less than mapQ
	 */
	public FilterNGS.SamRecordExtFilter getSamRecordMapQFilter(double mapQ) {
		return new SamRecordExtFilter(mapQ) {

			@Override
			public boolean filterOut(SAMRecord first, SAMRecord second) {
				// TODO Auto-generated method stub
				return false;
			}

			@Override
			public boolean filterOut(SAMRecord record) {
				int mapQual = record.getMappingQuality();
				if (mapQual == 255 || (double) mapQual < doubleFilter) {
					// TODO Auto-generated method stub
					return true;
				} else {
					return false;
				}
			}
		};
	}

	/**
	 * @return filter that will remove invalid records
	 */
	public SamRecordFilter getValidRecordFilter() {
		return new SamRecordFilter() {

			@Override
			public boolean filterOut(SAMRecord first, SAMRecord second) {
				// TODO Auto-generated method stub
				return false;
			}

			@Override
			public boolean filterOut(SAMRecord record) {
				// TODO Auto-generated method stub
				return record.isValid() != null;
			}
		};
	}

	/**
	 * @return filter that removes reads with reference of star
	 */
	public SamRecordFilter getValidReferenceFilter() {
		return new SamRecordFilter() {

			@Override
			public boolean filterOut(SAMRecord first, SAMRecord second) {
				// TODO Auto-generated method stub
				return false;
			}

			@Override
			public boolean filterOut(SAMRecord record) {
				// TODO Auto-generated method stub
				return record.getReferenceName().equals("*");
			}
		};
	}

	public SamRecordFilter getProperlyPairedFilter() {
		return new SamRecordFilter() {

			@Override
			public boolean filterOut(SAMRecord first, SAMRecord second) {
				// TODO Auto-generated method stub
				return false;
			}

			@Override
			public boolean filterOut(SAMRecord record) {
				return !record.getProperPairFlag();
			}
		};
	}

	public enum SAM_FILTER_TYPE {
		GENOTYPE, COPY_NUMBER;
	}

	/**
	 * @return pretty standard filteres for samrecords
	 */
	public ArrayList<SamRecordFilter> getStandardSAMRecordFilters(SAM_FILTER_TYPE type, Logger log) {
		ArrayList<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
		filters.add(new DuplicateReadFilter());
		filters.add(new AlignedFilter(true));
		filters.add(new SecondaryAlignmentFilter());
		filters.add(getValidRecordFilter());
		filters.add(getValidReferenceFilter());
		switch (type) {
		case COPY_NUMBER:
			break;
		case GENOTYPE:
			filters.add(getProperlyPairedFilter());// This can be evidence of a CNV
			break;
		default:
			log.reportTimeError("Invalid filter type" + type);
			filters = null;
			break;

		}

		return filters;
	}

	private static final String ESP_FILTER = "(esp6500si_all=='.'||esp6500si_all <=";
	private static final String ESPV2_FILTER = "(esp6500siv2_all=='.'||esp6500siv2_all <=";
	private static final String G10002014_FILTER = "(g10002014oct_all=='.'||g10002014oct_all <=";
	private static final String G10002015_FILTER = "(g10002015aug_all=='.'||g10002015aug_all <=";
	private static final String POPFREQ_MAXFILTER = "(PopFreqMax=='.'||PopFreqMax <=";
	private static final String AND = "&&";

	public static String getPopFreqFilterString(double maf) {
		String freq = ESP_FILTER + maf + ")" + AND + G10002014_FILTER + maf + ")" + AND + G10002015_FILTER + maf + ")" + AND + ESPV2_FILTER + maf + ")" + AND + POPFREQ_MAXFILTER + maf + ")";
		return freq;
	}

	/**
	 * Currently paramaterized to be similar to http://www.nature.com/cr/journal/v25/n3/extref/cr201520x14.pdf <br>
	 * The standards to reduce false positives were as follows: (1) a minimum depth of 10x in both tumors and normal pairs; (2) read depths of variant alleles in tumors should be more than 4x; (3) allelic fractions in tumors should be more than 20%." <br>
	 * NOTE: we also add a tlod, nlod, and normal alt count filter
	 */
	public static VariantContextFilter getTumorNormalFilter(double maf, Logger log) {

		VARIANT_FILTER_DOUBLE dpMut = VARIANT_FILTER_DOUBLE.AD_TUMOR;
		dpMut.setDFilter(10);

		VARIANT_FILTER_DOUBLE dpNormal = VARIANT_FILTER_DOUBLE.AD_NORMAL;
		dpNormal.setDFilter(10);

		VARIANT_FILTER_DOUBLE altADTumor = VARIANT_FILTER_DOUBLE.ALT_AD_TUMOR;
		altADTumor.setDFilter(4);

		VARIANT_FILTER_DOUBLE altAdNormal = VARIANT_FILTER_DOUBLE.ALT_AD_NORMAL;
		altAdNormal.setDFilter(0);

		VARIANT_FILTER_DOUBLE mutAF = VARIANT_FILTER_DOUBLE.AF_TUMOR;
		mutAF.setDFilter(.2);

		VARIANT_FILTER_DOUBLE tlod = VARIANT_FILTER_DOUBLE.TLOD;
		VARIANT_FILTER_DOUBLE nlod = VARIANT_FILTER_DOUBLE.NLOD;

		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		// VARIANT_FILTER_BOOLEAN mutectFail = VARIANT_FILTER_BOOLEAN.MUTECT_FAIL_FILTER;

		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { tlod, nlod, dpMut, dpNormal, altADTumor, altAdNormal, mutAF };

		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, new String[] { "G1000Freq" }, new String[] { FilterNGS.getPopFreqFilterString(maf) }, log);
		return vContextFilter;

	}

}
