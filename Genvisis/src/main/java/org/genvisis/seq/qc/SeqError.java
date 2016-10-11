package org.genvisis.seq.qc;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * @author lane0212 Class to compute error rates between samples akin to
 *         http://onlinelibrary.wiley.com/doi/10.1002/gepi.21881/epdf
 *
 */
public class SeqError {
	public static final String[] OUTPUT_HEADER = new String[] {	"Samp1", "Samp2", "total", "missing",
																															"matched", "proportionMissing",
																															"proportionAgree", "averageGC"};
	private final String vcfFile;
	private DuplicateETwo[] dETwos;
	private final Logger log;

	public SeqError(String vcfFile, DuplicateETwo[] dETwos, Logger log) {
		super();
		this.vcfFile = vcfFile;
		this.dETwos = dETwos;
		this.log = log;
	}

	public DuplicateETwo[] getdETwos() {
		return dETwos;
	}

	/**
	 * @param vContextFilterVariant filter applied to the entire variant context
	 * @param vContextFilterSample filter applied to the pairwise comparison
	 * @param filterNGS
	 * @param numthreads if a large number of pairwise comparisons are performed, the number of
	 *        threads can be increased
	 */
	public void populateError(VariantContextFilter setFilter, ReferenceGenome referenceGenome,
														int numVariantsToTest, int numthreads) {

		VCFFileReader reader = new VCFFileReader(new File(vcfFile), true);
		VCFHeader header = reader.getFileHeader();
		HashSet<String> allDEs = DuplicateETwo.getUniqSamples(dETwos);

		log.reportTimeInfo("Computing concordance for (" + dETwos.length + ") comparisons per variant");
		int numTotal = 0;
		int numSetPass = 0;
		WorkerTrain<DuplicateETwo> train = new WorkerTrain<SeqError.DuplicateETwo>(	null, numthreads,
																																								numthreads, log);
		train.setAutoShutDown(false);
		long time = System.currentTimeMillis();
		for (VariantContext vcTmp : reader) {
			numTotal++;
			if (numTotal % 10000 == 0) {
				log.reportTimeInfo(numTotal	+ " variants processed...with " + numSetPass
														+ " passing the set filter " + ext.getTimeElapsed(time));
				time = System.currentTimeMillis();
			}
			if (numVariantsToTest >= 0 && numSetPass == numVariantsToTest) {
				log.reportTimeInfo(numTotal	+ " variants processed...," + numVariantsToTest
														+ " variants to test reached " + ext.getTimeElapsed(time));
				reader.close();
				train.shutdown();
				return;
			}
			if (setFilter == null || setFilter.filter(vcTmp).passed()) {
				// System.out.println(vcTmp.getCommonInfo().getAttribute("esp6500si_all"));
				if (VCOps.getAAC(vcTmp, allDEs) > 0) {
					VariantContext vc = vcTmp.fullyDecode(header, false);
					numSetPass++;
					DuplicateProducer producer = new DuplicateProducer(vc, dETwos, referenceGenome);

					train.setProducer(producer);
					int tmpI = 0;
					dETwos = new DuplicateETwo[dETwos.length];
					while (train.hasNext()) {
						dETwos[tmpI] = train.next();
						tmpI++;
					}
				}
			}
		}
		reader.close();
		train.shutdown();
	}

	public void summarize(String fullPathToOutput) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.println(Array.toStr(OUTPUT_HEADER));
			for (DuplicateETwo dETwo : dETwos) {
				writer.println(dETwo.getSummary());
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathToOutput);
			log.reportException(e);
		}
	}

	private static class DuplicateProducer extends AbstractProducer<DuplicateETwo> {
		private final VariantContext vc;
		private final DuplicateETwo[] dETwos;
		private final ReferenceGenome referenceGenome;
		private int index;

		private DuplicateProducer(final VariantContext vc, final DuplicateETwo[] dETwos,
															final ReferenceGenome referenceGenome) {
			super();
			this.vc = vc;
			this.dETwos = dETwos;
			index = 0;
			this.referenceGenome = referenceGenome;
		}

		@Override
		public boolean hasNext() {
			return index < dETwos.length;
		}

		@Override
		public Callable<DuplicateETwo> next() {
			DuplicateWorker worker = new DuplicateWorker(vc, dETwos[index], referenceGenome);
			index++;
			return worker;
		}
	}

	private static class DuplicateWorker implements Callable<DuplicateETwo> {
		private final VariantContext vc;
		private final DuplicateETwo deTwo;
		private final ReferenceGenome referenceGenome;

		private DuplicateWorker(VariantContext vc, DuplicateETwo deTwo,
														ReferenceGenome referenceGenome) {
			super();
			this.vc = vc;
			this.deTwo = deTwo;
			this.referenceGenome = referenceGenome;
		}

		@Override
		public DuplicateETwo call() throws Exception {
			deTwo.addVC(vc, referenceGenome);
			return deTwo;
		}
	}

	public static DuplicateETwo[] getDups(String[] allSamps, DUPLICATE_COMP_TYPE type, MODE mode,
																				VariantContextFilter vContextFilterSample,
																				FilterNGS filterNGS, Logger log) {
		DuplicateETwo[] dETwos = new DuplicateETwo[(allSamps.length * (allSamps.length - 1)) / 2];
		int index = 0;
		for (int i = 0; i < allSamps.length; i++) {
			for (int j = i + 1; j < allSamps.length; j++) {
				HashSet<String> curDups = new HashSet<String>();
				curDups.add(allSamps[i]);
				curDups.add(allSamps[j]);
				dETwos[index] =
											new DuplicateETwo(curDups, type, mode, vContextFilterSample, filterNGS, log);
				index++;
			}
		}
		return dETwos;
	}

	public enum DUPLICATE_COMP_TYPE {
																		/**
																		 * Variant must past the sample filter for both
																		 */
																		ALL_PASS,
																		/**
																		 * Variant must past the sample filter for one
																		 */
																		ONE_PASS,
		/**
		 * // * // * Variant must past the sample filter on the average of the two //
		 */
		// AVERAGE_PASS
		;
	}

	public enum MODE {
										/**
										 * Will only compute concordance if both samples are called and pass according
										 * to {@link DUPLICATE_COMP_TYPE}
										 */
										BOTH_MUST_BE_CALLED,
										/**
										*
										*/
										ONE_MUST_BE_CALLED;
	}

	/**
	 * Class for tracking the e2 error rate across all comparisons
	 *
	 */
	public static class DuplicateETwo {
		private final Set<String> dups;
		private final DUPLICATE_COMP_TYPE type;
		private final MODE mode;
		private final VariantContextFilter vContextFilterSample;
		private final FilterNGS filterNGS;
		private int missing;
		private int matched;
		private int total;
		private double averageGC;
		private final Logger log;

		// private Logger log;

		public DuplicateETwo(	Set<String> dups, DUPLICATE_COMP_TYPE type, MODE mode,
													VariantContextFilter vContextFilterSample, FilterNGS filterNGS,
													Logger log) {
			super();
			this.dups = dups;
			this.type = type;
			this.mode = mode;
			this.vContextFilterSample = vContextFilterSample;
			this.filterNGS = filterNGS;
			this.log = log;
		}

		public VariantContextFilter getvContextFilterSample() {
			return vContextFilterSample;
		}

		public FilterNGS getFilterNGS() {
			return filterNGS;
		}

		public Set<String> getDups() {
			return dups;
		}

		public String getSummary() {
			String summary = "";
			int index = 0;
			for (String dup : dups) {
				summary += (index == 0 ? "" : "\t") + dup;
				index++;
			}
			double percentMissing = (double) missing / total;
			double percentMatched = (double) matched / total;
			averageGC = averageGC / total;
			// percentMissing = 1 - percentMissing;
			// percentMatched = 1 - percentMatched;
			summary += "\t" + total;
			summary += "\t" + missing;
			summary += "\t" + matched;
			summary += "\t" + percentMissing;
			summary += "\t" + percentMatched;
			summary += "\t" + averageGC;
			return summary;
		}

		private void addVC(VariantContext vc, ReferenceGenome referenceGenome) {
			VariantContext vcSub = VCOps.getSubset(vc, dups, VC_SUBSET_TYPE.SUBSET_LOOSE);
			VariantContext vcAlts = VCOps.getAltAlleleContext(vcSub, null, null,
																												ALT_ALLELE_CONTEXT_TYPE.ALL, false, log);// start
																																																	// with
																																																	// unfiltered
																																																	// easy
																																																	// test;
			if (vcAlts.getSampleNames().size() > 0) {// no variant calls, we do not care
				boolean tally = true;
				VariantContext vcFilteredAlts = null;
				switch (type) {
					case ALL_PASS:
						vcFilteredAlts = VCOps.getAltAlleleContext(	vcSub, filterNGS, vContextFilterSample,
																												ALT_ALLELE_CONTEXT_TYPE.ALL, false, log);
						tally = vcFilteredAlts.getSampleNames().size() == vcAlts.getSampleNames().size();// all
																																															// dup
																																															// variants
																																															// pass
						if (tally) {
							tally = vContextFilterSample == null	? true
																										: VCOps	.getIndividualPassingContext(vcSub,
																																												vContextFilterSample,
																																												log)
																														.getSampleNames().size() == dups.size();// all
																																																		// dups
																																																		// pass

						}
						break;
					case ONE_PASS:
						vcFilteredAlts = VCOps.getAltAlleleContext(	vcSub, filterNGS, vContextFilterSample,
																												ALT_ALLELE_CONTEXT_TYPE.ALL, false, log);
						tally = vcFilteredAlts.getSampleNames().size() > 0;// one of the dup variants pass
						if (tally) {
							tally = vContextFilterSample == null	? true
																										: VCOps	.getIndividualPassingContext(vcSub,
																																												vContextFilterSample,
																																												log)
																														.getSampleNames().size() > 0;
						}
						break;
					default:
						log.reportTimeError("Invalid Comparison type " + type);
						break;
				}
				if (tally) {
					if (vcFilteredAlts.getSampleNames().size() < dups.size()
							&& vcFilteredAlts.getSampleNames().size() > 0) {
						tallyMissingAlts();
						tally = false;
					} else if (vcFilteredAlts.getSampleNames().size() > 0) {// if ==0, all were filtered
																																	// out,so no countem
						GenotypesContext gc = vcSub.getGenotypes();

						boolean allCalled = true;
						boolean allMatched = true;
						for (Genotype g : gc) {
							for (Genotype g2 : gc) {
								if (g.isNoCall()) {
									allCalled = false;
								}
								if (!g.sameGenotype(g2)) {
									allMatched = false;
								}
							}
						}
						if (allCalled && allMatched) {
							matched++;
							total++;
							if (referenceGenome != null) {
								averageGC += referenceGenome.getGCContentFor(vcSub);
							}
						}
						if (!allMatched) {
							tallyMissingAlts();
						}
					}
				}
			}

		}

		private void tallyMissingAlts() {
			switch (mode) {
				case BOTH_MUST_BE_CALLED:
					break;
				case ONE_MUST_BE_CALLED:
					missing++;
					total++;
					break;
				default:
					log.reportTimeError("Invalid Mode type " + type);
					break;

			}
		}

		public static HashSet<String> getUniqSamples(DuplicateETwo[] dETwos) {
			HashSet<String> uniq = new HashSet<String>();
			for (DuplicateETwo dETwo : dETwos) {
				uniq.addAll(dETwo.getDups());
			}
			return uniq;
		}
	}
}

// private void addVC(VariantContext vc) {
// VariantContext vcSub = VCOps.getSubset(vc, dups, VC_SUBSET_TYPE.SUBSET_STRICT);
// VariantContext vcAlts = VCOps.getAltAlleleContext(vcSub, null, null, ALT_ALLELE_CONTEXT_TYPE.ALL,
// log);// start with unfiltered easy test;
// if (vcAlts.getSampleNames().size() > 0) {// no variant calls, we do not care
//
// boolean tally = true;
// VariantContext vcFilteredAlts = null;
// if (vContextFilterSample == null) {
// switch (type) {
// // case AVERAGE_PASS:
// // vcFilteredAlts = VCOps.getAltAlleleContext(vcSub, null, null, ALT_ALLELE_CONTEXT_TYPE.ALL,
// log);
// // tally = vContextFilterSample.filter(vcSub).passed();
// // log.reportTimeWarning("Not sure the best method here any more");
// // break;
// case ALL_PASS:
// vcFilteredAlts = VCOps.getAltAlleleContext(vcSub, filterNGS, vContextFilterSample,
// ALT_ALLELE_CONTEXT_TYPE.ALL, log);
// tally = vcFilteredAlts.getSampleNames().size() == vcAlts.getSampleNames().size();
// break;
// case ONE_PASS:
// vcFilteredAlts = VCOps.getAltAlleleContext(vcSub, filterNGS, vContextFilterSample,
// ALT_ALLELE_CONTEXT_TYPE.ALL, log);
// tally = vcFilteredAlts.getSampleNames().size() > 0;
// break;
// default:
// log.reportTimeError("Invalid Comparison type " + type);
// break;
// }
// }
// if (tally) {
// tally = vcFilteredAlts.getSampleNames().size() > 0;// all alts were filtered out
// }
// if (tally) {
// if (vcFilteredAlts.getSampleNames().size() < dups.size()) {// one of the alts was filtered out
// tallyMissingAlts();
// } else {
// GenotypesContext gc = vcSub.getGenotypes();
// if (vcSub.getNoCallCount() < dups.size()) {
// total++;
// boolean allCalled = true;
// boolean allMatched = true;
// for (Genotype g : gc) {
// for (Genotype g2 : gc) {
// if (g.isNoCall()) {
// allCalled = false;
// }
// if (!g.sameGenotype(g2)) {
// allMatched = false;
// }
// }
// }
// if (allCalled && allMatched) {
// matched++;
// }
// if (!allMatched) {
// tallyMissingAlts();
// }
// }
// }
// }
// }
// }
