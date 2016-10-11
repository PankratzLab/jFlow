package org.genvisis.seq.manage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

import org.genvisis.common.AlleleFreq;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;

/**
 * Class for common actions on {@link VariantContext} objects
 *
 */
public class VCOps {
	private static final String SNPEFF_GENE_NAME = "SNPEFF_GENE_NAME";
	private static final String SNPEFF_IMPACT = "SNPEFF_IMPACT";
	private static final String[] SNPEFF_IMPACT_IMPACTS = new String[] {"HIGH", "MODERATE", "LOW"};

	public enum GENOTYPE_INFO {
															GQ("GQ"), AD_REF("AD"), AD_ALT("AD"), DP("DP"), AD_TUMOR("AD_TUMOR"), AD_NORMAL("AD_NORMAL"), ALT_AD_TUMOR("AD_TUMOR"), AF_TUMOR("AF"), ALT_AD_NORMAL("AD_NORMAL"), MUTECT_FILTERS("MUTF"), TLOD("TLOD"), NLOD("NLOD");


		private final String flag;

		private GENOTYPE_INFO(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}
	}

	public enum GENOTYPE_FLAG_INFO {
																	HQ_DNM("HQ_DNM"), EHQ_DNM("EHQ_DNM");

		private final String flag;

		private GENOTYPE_FLAG_INFO(String flag) {
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
														DP("DP"),
														/**
														 * Mapping quality
														 */
														MQ("MQ"),

														/**
														 * Mapping quality
														 */
														GQ("GQ");

		/**
		 * flag represented in vcf file
		 */
		private final String flag;

		private COMMON_INFO(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}

	}

	/**
	 * @param annosToGet String keys to retrieve from the vc's common info
	 * @param vc
	 * @param defaultValue if anno is not found, this will be returned
	 * @return
	 */
	public static String[] getAnnotationsFor(	String[] annosToGet, VariantContext vc,
																						String defaultValue) {
		String[] annos = new String[annosToGet.length];
		for (int i = 0; i < annos.length; i++) {
			annos[i] = vc.getCommonInfo().getAttributeAsString(annosToGet[i], defaultValue);
		}
		return annos;
	}

	public static String getSNP_EFFGeneName(VariantContext vc) {
		String geneName = ".";
		return getAnnotationsFor(new String[] {SNPEFF_GENE_NAME}, vc, geneName)[0];
	}

	public static String getSNP_EFFImpact(VariantContext vc) {
		String impact = ".";
		return getAnnotationsFor(new String[] {SNPEFF_IMPACT}, vc, impact)[0];
	}

	public static boolean isHighModLowSNP_EFFImpact(VariantContext vc) {
		return ext.indexOfStr(getSNP_EFFImpact(vc), SNPEFF_IMPACT_IMPACTS) >= 0;
	}

	public static double getMAF(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = sampleNames == null ? vc : getSubset(vc, sampleNames);
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
		VariantContext vcSub = sampleNames == null ? vc : getSubset(vc, sampleNames);
		int[] alleleCounts = getAlleleCounts(vcSub);
		return (alleleCounts[2] * 2 + alleleCounts[1]);

	}

	/**
	 * Tests whether the minor allele is the alternate allele
	 */
	public static boolean isMinorAlleleAlternate(VariantContext vc, Set<String> sampleNames) {
		VariantContext vcSub = sampleNames == null ? vc : getSubset(vc, sampleNames);
		double mac = getMAC(vcSub, null);
		double aac = getAAC(vcSub, null);
		return aac == mac;
	}

	/**
	 * @param names the names corresponding the the array of jexl expressions
	 * @param expressions the expressions to be evaluated
	 * @param log
	 * @return
	 */
	public static List<VariantContextUtils.JexlVCMatchExp> getJexlVCMathExp(String[] names,
																																					String[] expressions,
																																					Logger log) {
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
	 * @param jExp determine if the expression matches for the {@link VariantContext}
	 * @return
	 */
	public static boolean passesJexl(VariantContext vc, VariantContextUtils.JexlVCMatchExp jExp) {
		return VariantContextUtils.match(vc, jExp);
	}

	/**
	 * @param vc
	 * @param jExps determine if all the expressions match for the {@link VariantContext}
	 * @return
	 */
	public static boolean passesJexls(VariantContext vc,
																		List<VariantContextUtils.JexlVCMatchExp> jExps) {
		for (VariantContextUtils.JexlVCMatchExp jExp : jExps) {
			if (!passesJexl(vc, jExp)) {
				return false;
			}
		}
		return true;
	}

	public static Segment getSegment(VariantContext vc) {
		return new Segment(	Positions.chromosomeNumber(vc.getContig()), vc.getStart(),
												vc.getEnd() <= vc.getStart() ? vc.getStart() : vc.getEnd());
	}

	/**
	 * Averages an info value across samples
	 */
	public static double getAverageCommonInfo(VariantContext vc, Set<String> sampleNames,
																						COMMON_INFO info) {
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
																				 * All genotypes with a one or two non-reference allele(s)
																				 * will be returned
																				 */
																				HET_ONLY,
																				/**
																				 * All genotypes with a homozygous non-reference allele will
																				 * be returned
																				 */
																				HOM_ONLY,
																				/**
																				 * All genotypes with a non-reference allele will be
																				 * returned
																				 */
																				ALL;
	}

	public static VariantContext getAltAlleleContext(	final VariantContext vc, FilterNGS filterNGS,
																										Logger log) {
		return getAltAlleleContext(vc, filterNGS, null, ALT_ALLELE_CONTEXT_TYPE.ALL, log);
	}

	public static VariantContext getIndividualPassingContext(	final VariantContext vc,
																														final VariantContextFilter variantContextFilter,
																														Logger log) {
		HashSet<String> passing = new HashSet<String>();
		Set<String> curSamps = vc.getSampleNames();
		for (String samp : curSamps) {
			HashSet<String> tmp = new HashSet<String>();
			tmp.add(samp);
			if (variantContextFilter.filter(getSubset(vc, tmp)).passed()) {
				passing.add(samp);
			}
		}
		return getSubset(vc, passing);
	}

	/**
	 * @param vc
	 * @param filterNGS used for alternate allele depth,note that this depth will be applied to either
	 *        the one non-ref alternate allele, or both non ref alternate alleles if applicable
	 * @param variantContextFilter will be applied to each sample with an alternate call
	 * @param type see {@link ALT_ALLELE_CONTEXT_TYPE}
	 * @param log
	 * @return
	 */
	public static VariantContext getAltAlleleContext(	final VariantContext vc, FilterNGS filterNGS,
																										final VariantContextFilter variantContextFilter,
																										final ALT_ALLELE_CONTEXT_TYPE type,
																										final Logger log) {
		return getAltAlleleContext(vc, filterNGS, variantContextFilter, type, true, log);
	}

	public static VariantContext getAltAlleleContext(	final VariantContext vc, FilterNGS filterNGS,
																										final VariantContextFilter variantContextFilter,
																										final ALT_ALLELE_CONTEXT_TYPE type,
																										boolean verbose, final Logger log) {
		GenotypesContext gc = vc.getGenotypes();
		HashSet<String> samplesWithAlt = new HashSet<String>();
		int altAlleleDepth = -1;
		double altAlleleDepthRatio = -1;
		if (filterNGS != null	&& filterNGS.getAltAlleleDepthFilter() != null
				&& filterNGS.getAltAlleleDepthFilter()[0] > 0) {
			// log.reportTimeError("Alt Allele depth filter is currently not in here, JOHN");
			altAlleleDepth = filterNGS.getAltAlleleDepthFilter()[0];
			// return null;
		}
		if (filterNGS != null	&& filterNGS.getAltAlleleDepthRatioFilter() != null
				&& filterNGS.getAltAlleleDepthRatioFilter()[0] > 0) {
			// log.reportTimeError("Alt Allele depth filter is currently not in here, JOHN");
			altAlleleDepthRatio = filterNGS.getAltAlleleDepthRatioFilter()[0];
			// return null;
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
				int[] AD = new int[] {0, 0};
				if (altAlleleDepth >= 0 || altAlleleDepthRatio >= 0) {
					try {
						if (geno.hasAD()) {
							AD = getAppropriateAlleleDepths(vc, geno, verbose, log);
						} else {
							AD = new int[] {altAlleleDepth + 1, altAlleleDepth + 1};
							if (verbose) {
								log.reportTimeWarning(geno.toString()
																				+ " did not have allele depths, setting depths to "
																			+ Array.toStr(AD));
							}
						}
					} catch (IllegalStateException ise) {
						if (verbose) {
							// TODO, report?
							log.reportTimeError("Could not compute appropriate allele Depths");
							log.reportException(ise);
						}
					}
				}
				// TODO gte vs gt
				if (altAlleleDepth <= 0 || AD[1] >= altAlleleDepth) {
					if (altAlleleDepth <= 0 || !geno.isHetNonRef() || AD[0] >= altAlleleDepth) {// handles the
																																											// case when
																																											// both
																																											// alleles are
																																											// non-reference
						double ratio = (double) AD[1] / geno.getDP();
						if (altAlleleDepthRatio <= 0 || ratio >= altAlleleDepthRatio) {
							HashSet<String> tmp = new HashSet<String>();
							tmp.add(geno.getSampleName());
							if (variantContextFilter == null
									|| variantContextFilter.filter(getSubset(vc, tmp)).passed()) {
								samplesWithAlt.add(geno.getSampleName());
							}
						}
					}
				}
			}
		}
		return getSubset(vc, samplesWithAlt);
	}

	/**
	 * Will return the average alt allele ratio for heterozygous calls in this context. Note that the
	 * ratio is not always representing the ratio to the reference genome (het non ref)<br>
	 * min/max will be used in the het-non ref case
	 */
	public static double getAverageHetAlleleRatio(VariantContext vc, FilterNGS filterNGS,
																								VariantContextFilter variantContextFilter,
																								Logger log) {
		double avg = 0;
		VariantContext vcAlts = getAltAlleleContext(vc, filterNGS, variantContextFilter,
																								ALT_ALLELE_CONTEXT_TYPE.HET_ONLY, log);
		GenotypesContext gc = vcAlts.getGenotypes();
		for (Genotype g : gc) {
			int[] AD = getAppropriateAlleleDepths(vc, g, false, log);
			if (g.isHetNonRef()) {
				int min = Math.min(AD[0], AD[1]);
				int max = Math.max(AD[0], AD[1]);
				avg += (double) min / max;
			} else if (AD[0] > 0) {
				avg += (double) AD[1] / AD[0];
			}
		}
		avg = avg / gc.size();
		return avg;
	}

	public static int[] getAppropriateAlleleDepths(	VariantContext vc, Genotype g, boolean verbose,
																									Logger log) {
		int[] AD = new int[2];
		if (!vc.isBiallelic()) {
			log.reportTimeWarning("JOHN REMEMBER THE BIALLELIC ISSUE!");
		}
		Arrays.fill(AD, 0);
		if (g.isNoCall()) {
			return AD;
		}
		List<Allele> gAlleles = g.getAlleles();

		List<Allele> varAlleles = vc.getAlleles();

		if (gAlleles.size() != 2 && g.hasAD()) {
			log.reportTimeError("Number of alleles must equal 2, return AD[0,0]");
			return AD;
		} else if (gAlleles.size() == 0 || !g.hasAD()) {
			String error = "Invalid Allele retrieval";
			error += g.toString();
			error += vc.toStringWithoutGenotypes();
			throw new IllegalStateException(error);

		} else {
			int[] gAD = g.getAD();
			if (gAlleles.size() == 1) {
				AD[0] = gAD[0];
				AD[1] = gAD[1];
			} else {
				if (gAD.length != varAlleles.size()) {
					if (verbose) {
						String error = "an allelic depth must be present for every allele in the variant";
						log.reportTimeError(error);
						log.reportTimeError(varAlleles.toString());
						log.reportTimeError(gAlleles.toString());
						log.reportTimeError(Array.toStr(gAD) + " > " + Array.toStr(AD));
						log.reportTimeError(vc.toStringWithoutGenotypes());
						log.reportTimeError(g.toString());
					}
					return null;
				}
				ArrayList<Allele> uniqs = new ArrayList<Allele>();
				for (int i = 0; i < gAlleles.size(); i++) {
					boolean have = false;
					for (int j = 0; j < uniqs.size(); j++) {
						if (gAlleles.get(i).equals(uniqs.get(j), false)) {
							have = true;
						}
					}
					if (!have) {
						uniqs.add(gAlleles.get(i));
					}
				}
				if (varAlleles.size() < uniqs.size()) {
					throw new IllegalStateException("Variant does not capture genotyped alleles\t"
																					+ varAlleles.toString() + "\t" + uniqs.toString());
				}
				ArrayList<Integer> gAlleleIndices = new ArrayList<Integer>();
				int index = 0;
				for (Allele varAllele : varAlleles) {
					for (Allele gAllele : uniqs) {
						if (varAllele.equals(gAllele)) {
							gAlleleIndices.add(index);
						}
					}
					index++;
				}

				if (gAlleleIndices.size() > 2 || gAlleleIndices.size() == 0) {
					if (verbose) {

						log.reportTimeError("Invalid Allelic depth extraction");
						log.reportTimeError(varAlleles.toString());
						log.reportTimeError(gAlleles.toString());
						log.reportTimeError(Array.toStr(gAD) + " > " + Array.toStr(AD));
						log.reportTimeError(vc.toStringWithoutGenotypes());
						log.reportTimeError(g.toString());
					}
					throw new IllegalStateException("Invalid Allelic depth extraction");
				}
				if (gAlleleIndices.size() == 1) {
					if (g.isHomRef()) {
						AD[0] = gAD[gAlleleIndices.get(0)];
					} else if (g.isHomVar()) {
						AD[1] = gAD[gAlleleIndices.get(0)];
					} else {
						if (verbose) {

							log.reportTimeError("Invalid Allelic depth extraction for 1 index");
							log.reportTimeError(varAlleles.toString());
							log.reportTimeError(gAlleles.toString());
							log.reportTimeError(Array.toStr(gAD) + " > " + Array.toStr(AD));
							log.reportTimeError(vc.toStringWithoutGenotypes());
							log.reportTimeError(g.toString());
						}
						throw new IllegalStateException("Invalid Allelic depth extraction for 1 index");
					}
				} else {
					for (int i = 0; i < gAlleleIndices.size(); i++) {
						AD[i] = gAD[gAlleleIndices.get(i)];
					}
				}

				if (g.isHet()	&& (AD[1] == 0 || AD[0] == 0) && Array.sum(gAD) != AD[0]
						&& Array.sum(gAD) != AD[1]) {// there can actually be het calls with 0 ref or 0 alt, or
																					// both...apparently, I would'nt do that but whatever. So
																					// anyways we do not test AD[0]
					if (verbose) {

						log.reportTimeError("Invalid Het allele depth, Het non-ref " + g.isHetNonRef());
						log.reportTimeError(varAlleles.toString());
						log.reportTimeError(gAlleles.toString());
						log.reportTimeError(Array.toStr(gAD) + " > " + Array.toStr(AD));
						log.reportTimeError(vc.toStringWithoutGenotypes());
						log.reportTimeError(g.toString());
					}
					throw new IllegalStateException("Invalid Het allele depth");
				} else if (g.isHomVar()	&& AD[1] == 0 && Array.sum(gAD) > 0
										&& Array.sum(gAD) != Array.sum(AD)) {
					if (verbose) {

						log.reportTimeError("Invalid Hom Var allele depth");
						log.reportTimeError(varAlleles.toString());
						log.reportTimeError(gAlleles.toString());
						log.reportTimeError(Array.toStr(gAD) + " > " + Array.toStr(AD));
						log.reportTimeError(vc.toStringWithoutGenotypes());
						log.reportTimeError(g.toString());
					}
					throw new IllegalStateException("Invalid Hom Var  allele depth");
				} else if (g.isHomRef() && AD[0] == 0 && Array.sum(gAD) > 0) {
					if (verbose) {

						log.reportTimeError("Invalid Hom Ref allele depth");
						log.reportTimeError(varAlleles.toString());
						log.reportTimeError(gAlleles.toString());
						log.reportTimeError(Array.toStr(gAD) + " > " + Array.toStr(AD));
						log.reportTimeError(vc.toStringWithoutGenotypes());
						log.reportTimeError(g.toString());
					}
					throw new IllegalStateException("Invalid Hom Ref allele depth");
				}
			}
		}

		return AD;
		// //vc
		// g.ge
		// g.getAlleles();
	}

	public static boolean getFlagEqualsInfo(VariantContext vc, Set<String> sampleNames,
																					GENOTYPE_FLAG_INFO info, String equal, Logger log) {
		VariantContext vcSub = sampleNames == null ? vc : getSubset(vc, sampleNames);
		GenotypesContext gc = vcSub.getGenotypes();
		for (Genotype geno : gc) {
			if (!geno.hasAnyAttribute(info.getFlag())
					|| !geno.getAnyAttribute(info.getFlag()).toString().equals(equal)) {
				return false;
			}
		}
		return true;
	}

	public static double getAvgGenotypeInfo(VariantContext vc, Set<String> sampleNames,
																					GENOTYPE_INFO info, Logger log) {
		double avgGI = 0;
		int numWith = 0;
		VariantContext vcSub = sampleNames == null ? vc : getSubset(vc, sampleNames);
		GenotypesContext gc = vcSub.getGenotypes();
		for (Genotype geno : gc) {
			if (geno.hasAnyAttribute(info.getFlag())) {
				numWith++;
				switch (info) {
					case AD_REF:
						avgGI += getAppropriateAlleleDepths(vcSub, geno, false, log)[0];
						break;
					case AD_ALT:
						avgGI += getAppropriateAlleleDepths(vcSub, geno, false, log)[1];
						break;
					case DP:
						avgGI += geno.getDP();
						break;
					case GQ:
						avgGI += geno.getGQ();
						break;
					case AD_TUMOR:
					case AD_NORMAL:
						double[] adTotal = Array.toDoubleArray(geno	.getAnyAttribute(info.getFlag()).toString()
																												.split(","));
						avgGI += Array.sum(adTotal);
						break;
					case ALT_AD_TUMOR:
					case ALT_AD_NORMAL:
						avgGI += Array.toDoubleArray(geno	.getAnyAttribute(info.getFlag()).toString()
																							.split(","))[1];
						break;
					case AF_TUMOR:
					case NLOD:
					case TLOD:
						avgGI += Array.toDoubleArray(geno	.getAnyAttribute(info.getFlag()).toString()
																							.split(","))[0];
						break;
					default:
						throw new IllegalArgumentException("Invalid genotype flag " + info);
				}
			}
		}
		if (numWith > 0) {
			avgGI = avgGI / numWith;
		} else {
			avgGI = Double.NaN;
		}
		return avgGI;
	}

	public enum VC_SUBSET_TYPE {
															/**
															 * Samples not contained in the vcf will be given missing genotypes
															 */
															SUBSET_LOOSE,
															/**
															 * A check will be performed and only samples present in the input set
															 * and the vcf file will be exported
															 */
															SUBSET_STRICT,
															/**
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

	public static Genotype getGenotypeFor(final VariantContext vc, final String sampleName,
																				VC_SUBSET_TYPE type) {
		return getSubset(vc, sampleName, type).getGenotype(0);
	}

	/**
	 * Subsets to particular samples
	 */
	public static VariantContext getSubset(	final VariantContext vc, final String sampleName,
																					VC_SUBSET_TYPE type) {
		HashSet<String> tmp = new HashSet<String>();
		tmp.add(sampleName);
		return getSubset(vc, tmp, type);
	}

	/**
	 * Subsets to particular samples
	 */
	public static VariantContext getSubset(final VariantContext vc, final Set<String> sampleNames) {
		return getSubset(vc, sampleNames, sampleNames == null	? VC_SUBSET_TYPE.NO_SUBSET
																													: VC_SUBSET_TYPE.SUBSET_STRICT);
	}

	public static VariantContext getSubset(	final VariantContext vc, final Set<String> sampleNames,
																					VC_SUBSET_TYPE type) {
		return getSubset(vc, sampleNames, type, true);
	}

	/**
	 * @param vc
	 * @return the number of samples that do not have homozygous reference calls or no calls
	 */
	public static int getNumWithAlt(VariantContext vc) {
		int numSamps = vc.getSampleNames().size();
		return (numSamps - vc.getHomRefCount() - vc.getNoCallCount());
	}

	/**
	 * Subsets to particular samples
	 */
	public static VariantContext getSubset(	final VariantContext vc, final Set<String> sampleNames,
																					VC_SUBSET_TYPE type,
																					boolean rederiveAllelesFromGenotypes) {
		VariantContext vcSub = null;
		switch (type) {
			case SUBSET_LOOSE:
				vcSub = vc.subContextFromSamples(sampleNames, rederiveAllelesFromGenotypes);
				break;
			case SUBSET_STRICT:
				vcSub = vc.subContextFromSamples(	getOverlap(sampleNames, vc.getSampleNames()),
																					rederiveAllelesFromGenotypes);
				break;
			case NO_SUBSET:
				vcSub = vc;
				break;
			default:
				break;

		}
		return vcSub;
	}

	/**
	 * @param vc
	 * @param samplesToSetMissing sample names in this hashset will be set to missing in the returned
	 *        {@link VariantContext}
	 * @return
	 */
	public static VariantContext setTheseSamplesToMissing(VariantContext vc,
																												HashSet<String> samplesToSetMissing) {
		VariantContextBuilder builder = new VariantContextBuilder(vc);
		if (samplesToSetMissing == null) {
			return builder.make();
		} else {
			Set<String> samples = vc.getSampleNames();
			ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
			for (String sample : samples) {

				GenotypeBuilder gbuilder = new GenotypeBuilder(vc.getGenotype(sample));
				if (samplesToSetMissing.contains(sample)) {
					gbuilder.alleles(GenotypeOps.getNoCall());
				}
				genotypes.add(gbuilder.make());
			}
			builder.genotypes(genotypes);
			return builder.make();
		}
	}

	public static GeneData[] getGenesThatOverlap(VariantContext vc, GeneTrack geneTrack, Logger log) {
		Segment vcSeg = getSegment(vc);
		return geneTrack.getOverlappingGenes(vcSeg);
	}

	public static boolean isInTheseSegments(VariantContext vc, Segment[] orderedSegs) {
		return getOverlappingSegments(vc, orderedSegs) != null;
	}

	public static int[] getOverlappingSegments(VariantContext vc, Segment[] orderedSegs) {
		Segment vcSeg = getSegment(vc);
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
		if ((vc.getHomRefCount()	+ vc.getHetCount() + vc.getHomVarCount()
					+ vc.getNoCallCount()) != vc.getNSamples()) {
			System.err.println("Un accounted for genotypes...");
		}
		return new int[] {vc.getHomRefCount(), vc.getHetCount(), vc.getHomVarCount()};
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
		String[][] ambiguousDefs = new String[][] {{"A", "T"}, {"T", "A"}, {"C", "G"}, {"G", "C"}};
		for (Allele b : alt) {
			for (String[] ambiguousDef : ambiguousDefs) {
				if (ref.equals(ambiguousDef[0]) && b.getDisplayString().equals(ambiguousDef[1])) {
					ambiguous = true;
					break;
				}
			}
		}
		return ambiguous;
	}

	public static List<Allele> getAllelesFor(VariantContext vc, String sampleName) {
		HashSet<String> tmp = new HashSet<String>();
		tmp.add(sampleName);
		return getAllelesFor(vc, tmp);
	}

	public static List<Allele> getAllelesFor(VariantContext vc, Set<String> sampleNames) {
		return getSubset(vc, sampleNames).getAlleles();
	}

	public static String allelesToString(List<Allele> alleles) {
		String tmp = "";
		boolean first = true;
		if (alleles.size() == 1) {
			tmp = alleles.get(0).getDisplayString() + "/" + alleles.get(0).getDisplayString();
		} else {
			for (Allele a : alleles) {
				tmp += first ? a.getDisplayString() : "/" + a.getDisplayString();
			}
		}
		return tmp;
	}

	public static class Transmission {
		private final String child;
		private final String p1;
		private final String p2;

		private List<Allele> childAlleles;
		private List<Allele> p1Alleles;
		private List<Allele> p2Alleles;
		private Allele ref;

		public Transmission(String child, String p1, String p2) {
			super();
			this.child = child;
			this.p1 = p1;
			this.p2 = p2;
		}

		public void parseAlleles(VariantContext vc) {
			childAlleles = getAllelesFor(vc, child);
			p1Alleles = getAllelesFor(vc, p1);
			p2Alleles = getAllelesFor(vc, p2);
			for (Allele a : vc.getAlleles()) {
				if (a.isReference()) {
					ref = a;
					break;
				}
			}
		}

		public String getSummary() {
			String summary = "";
			summary += "\t" + child + ": " + allelesToString(childAlleles);
			summary += "\t" + p1 + ": " + allelesToString(p1Alleles);
			summary += "\t" + p2 + ": " + allelesToString(p2Alleles);
			summary += "\t" + "ref: " + ref.getDisplayString();
			return summary;
		}
	}

	public static class LocusID {
		private final byte chr;
		private final int start;
		private final String ref;
		private final String alt;
		private final String id;

		public LocusID(VariantContext vc) {
			chr = Positions.chromosomeNumber(vc.getContig());
			start = vc.getStart();
			ref = vc.getReference().getDisplayString();
			alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();
			id = vc.getID().equals(".") ? chr + ":" + start + ":" + ref + ":" + alt : vc.getID();
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
