package org.genvisis.seq.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.genvisis.seq.qc.FilterNGS.VcFilterDouble;
import org.genvisis.seq.qc.SeqError.DUPLICATE_COMP_TYPE;
import org.genvisis.seq.qc.SeqError.DuplicateETwo;
import org.genvisis.seq.qc.SeqError.MODE;

public class SeqQCValidation {
	private static String[] SNP_SETS = new String[] {	"(esp6500si_all=='.'||esp6500si_all <= 0.01)&&(g10002014oct_all=='.'||g10002014oct_all <= 0.01)",
																										"snp138=='.'", "snp138!='.'"};
	// , "(esp6500si_all=='.'||esp6500si_all <= 0.01)"

	private final SeqError seqError;
	private final VariantContextFilter setFilter;
	private final ReferenceGenome referenceGenome;

	public SeqQCValidation(	SeqError seqError, VariantContextFilter setFilter,
													ReferenceGenome referenceGenome) {
		super();
		this.seqError = seqError;
		this.setFilter = setFilter;
		this.referenceGenome = referenceGenome;

	}

	public SeqError getSeqError() {
		return seqError;
	}

	public VariantContextFilter getSetFilter() {
		return setFilter;
	}

	public void validate(int numVariantsToTest, int numthreads) {
		seqError.populateError(setFilter, referenceGenome, numVariantsToTest, numthreads);
	}

	private static class SeqQCValidationWorker implements Callable<SeqQCValidation> {
		private final SeqQCValidation seqQCValidation;
		private final int numthreads;
		private final int numVariantsToTest;

		public SeqQCValidationWorker(	SeqQCValidation seqQCValidation, int numVariantsToTest,
																	int numthreads) {
			super();
			this.seqQCValidation = seqQCValidation;
			this.numthreads = numthreads;
			this.numVariantsToTest = numVariantsToTest;
		}

		@Override
		public SeqQCValidation call() throws Exception {
			seqQCValidation.validate(numVariantsToTest, numthreads);
			return seqQCValidation;
		}

	}

	private static class SeqQCValidationProducer extends AbstractProducer<SeqQCValidation> {

		private final SeqQCValidation[] validations;
		private int index;
		private final int numVariantsToTest;
		int numInternalThreads;
		// private Logger log;

		public SeqQCValidationProducer(	SeqQCValidation[] validations, int numVariantsToTest,
																		int numInternalThreads, Logger log) {
			super();
			this.validations = validations;
			this.numVariantsToTest = numVariantsToTest;
			this.numInternalThreads = numInternalThreads;
			// this.log = log;
			index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < validations.length;
		}

		@Override
		public Callable<SeqQCValidation> next() {
			SeqQCValidationWorker worker = new SeqQCValidationWorker(	validations[index],
																																numVariantsToTest,
																																numInternalThreads);
			index++;
			return worker;
		}
	}

	private static VariantContextFilter[] getSampleVariantFilters(Logger log) {
		double[] GQ_Values = new double[] {	-1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 90,
																				98};
		// double[] GQ_Values = new double[] { -1 };

		double[] VQSLOD_Values = new double[] {-1};
		double[] DEPTH = new double[] {-1, 3, 5, 8, 10, 12, 15, 20, 25};
		// double[] DEPTH = new double[] { -1 };

		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		VariantContextFilter[] vContextFilterSample = new VariantContextFilter[VQSLOD_Values.length
																																							* GQ_Values.length
																																						* DEPTH.length];
		int index = 0;
		for (double vqslod_Value : VQSLOD_Values) {
			for (double gq_Value : GQ_Values) {
				for (double element : DEPTH) {
					VARIANT_FILTER_DOUBLE vq = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
					vq.setDFilter(vqslod_Value);
					VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_STRICT;
					gq.setDFilter(gq_Value);
					VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
					dp.setDFilter(element);
					vContextFilterSample[index] =
																			new VariantContextFilter(	new VARIANT_FILTER_DOUBLE[] {vq, gq,
																																														dp},
																																new VARIANT_FILTER_BOOLEAN[] {fail},
																																null, null, log);
					index++;
				}
			}
		}
		return vContextFilterSample;
	}

	private static FilterNGS[] getAltAlleleDepthFilters(int startDepth, int stopDepth, int jumpDepth,
																											double startRatio, double stopRatio,
																											double jumpRatio) {
		ArrayList<FilterNGS> filterNGSs = new ArrayList<FilterNGS>();
		for (int i = startDepth; i <= stopDepth; i += jumpDepth) {
			for (double j = startRatio; j <= stopRatio; j += jumpRatio) {
				FilterNGS filterNGS = new FilterNGS(0, 0, null);
				filterNGS.setAltAlleleDepthFilter(new int[] {i});
				filterNGS.setAltAlleleDepthRatioFilter(new double[] {j});
				filterNGSs.add(filterNGS);
			}

		}
		return filterNGSs.toArray(new FilterNGS[filterNGSs.size()]);
	}

	private static String[][] loadPairWise(String fileOFSamplesForPairWise, Logger log) {
		ArrayList<String[]> comps = new ArrayList<String[]>();
		if (Files.getHeaderOfFile(fileOFSamplesForPairWise, log).length == 2) {
			log.reportTimeInfo("Assuming pair-wise comparisions have been pre-defined");
			try {
				BufferedReader reader = Files.getAppropriateReader(fileOFSamplesForPairWise);

				while (reader.ready()) {
					comps.add(reader.readLine().trim().split("[\t]+"));
				}

			} catch (FileNotFoundException e) {
				log.reportFileNotFound(fileOFSamplesForPairWise);
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				log.reportException(e);
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} else {
			log.reportTimeError("Not implemented yet");
		}
		return comps.toArray(new String[comps.size()][]);
	}

	public static void validate(String vcf, String fileOFSamplesForPairWise,
															String referenceGenomeFasta, int numVariantsToTest, int numthreads,
															Logger log) {
		String[][] sampsForPairWise = loadPairWise(fileOFSamplesForPairWise, log);
		String output = ext.addToRoot(fileOFSamplesForPairWise, ".validation.summary");
		// FilterNGS[] altDepthFilterNGSs = getAltAlleleDepthFilters(-1, 0, 1, 0.0, 0.1, 0.1);
		FilterNGS[] altDepthFilterNGSs = getAltAlleleDepthFilters(-1, 12, 1, 0.0, 0.7, 0.1);
		log.reportTimeInfo(altDepthFilterNGSs.length + " alt depth filters");
		VariantContextFilter[] vContextFilterSamples = getSampleVariantFilters(log);
		int numComp =
								altDepthFilterNGSs.length * vContextFilterSamples.length * sampsForPairWise.length;
		log.reportTimeInfo("Running "	+ (altDepthFilterNGSs.length * vContextFilterSamples.length)
												+ " validations over " + sampsForPairWise.length
												+ " pairwise comparisons for a total of " + numComp
												+ " tests on each variant");
		ReferenceGenome referenceGenome = referenceGenomeFasta == null	? null
																																		: new ReferenceGenome(referenceGenomeFasta,
																																													log);
		if (referenceGenome == null) {
			log.reportTimeWarning("Was unable to find a reference genome, skipping gc computations");
		} else {
			referenceGenome.setDefaultBuffer(50);
		}
		ArrayList<DuplicateETwo> deETwos = new ArrayList<SeqError.DuplicateETwo>(numComp);
		for (VariantContextFilter vContextFilterSample : vContextFilterSamples) {
			for (FilterNGS altDepthFilterNGS : altDepthFilterNGSs) {
				for (String[] element : sampsForPairWise) {
					HashSet<String> curDups = new HashSet<String>();
					curDups.add(element[0]);
					curDups.add(element[1]);
					deETwos.add(new DuplicateETwo(curDups, DUPLICATE_COMP_TYPE.ALL_PASS,
																				MODE.ONE_MUST_BE_CALLED, vContextFilterSample,
																				altDepthFilterNGS, log));
				}
			}
		}
		ArrayList<DuplicateETwo[]> duplicateETwos =
																							Array.splitUpArray(	deETwos.toArray(new DuplicateETwo[deETwos.size()]),
																																	numthreads, log);
		SeqQCValidation[] setSeqQCValidations = new SeqQCValidation[SNP_SETS.length
																																* duplicateETwos.size()];

		int index = 0;
		for (String element : SNP_SETS) {
			for (int j = 0; j < duplicateETwos.size(); j++) {
				SeqError seqError = new SeqError(vcf, duplicateETwos.get(j), log);
				VariantContextFilter setFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {},
																																	new VARIANT_FILTER_BOOLEAN[] {VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER},
																																	new String[] {element},
																																	new String[] {element}, log);
				setSeqQCValidations[index] = new SeqQCValidation(seqError, setFilter, referenceGenome);
				index++;
			}
		}
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println(Array.toStr(SeqError.OUTPUT_HEADER)	+ "\t"
											+ "SET\tAltDepth\tAltDepthRatio\tGQ\tVQSLOD\tDepth");
			SeqQCValidationProducer producer = new SeqQCValidationProducer(	setSeqQCValidations,
																																			numVariantsToTest, 1, log);
			WorkerTrain<SeqQCValidation> train = new WorkerTrain<SeqQCValidation>(producer, numthreads, 1,
																																						log);
			while (train.hasNext()) {
				SeqQCValidation tmp = train.next();
				DuplicateETwo[] dETwos = tmp.getSeqError().getdETwos();
				for (DuplicateETwo dETwo : dETwos) {
					writer.print(dETwo.getSummary()	+ "\t"
												+ tmp.getSetFilter().getvFilterJEXL().getjExps().get(0).name + "\t"
												+ dETwo.getFilterNGS().getAltAlleleDepthFilter()[0] + "\t"
												+ dETwo.getFilterNGS().getAltAlleleDepthRatioFilter()[0]);
					VcFilterDouble[] doubles = dETwo.getvContextFilterSample().getvDoubles();
					int GQIndex = -1;
					int VQSLODIndex = -1;
					int dpIndex = -1;
					for (int j = 0; j < doubles.length; j++) {
						if (doubles[j].getDfilter() == VARIANT_FILTER_DOUBLE.GQ_STRICT) {
							GQIndex = j;
						}
						if (doubles[j].getDfilter() == VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE) {
							VQSLODIndex = j;
						}
						if (doubles[j].getDfilter() == VARIANT_FILTER_DOUBLE.DP) {
							dpIndex = j;
						}
					}
					if (GQIndex >= 0) {
						writer.print("\t" + doubles[GQIndex].getFilterThreshold());
					} else {
						writer.print("\tNA");
					}
					if (VQSLODIndex >= 0) {
						writer.print("\t" + doubles[VQSLODIndex].getFilterThreshold());
					} else {
						writer.print("\tNA");
					}
					if (dpIndex >= 0) {
						writer.print("\t" + doubles[dpIndex].getFilterThreshold());
					} else {
						writer.print("\tNA");
					}
					writer.println();
				}
				writer.flush();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "SeqErrorIterator.vcf";
		String fileOFSamplesForPairWise = null;
		String referenceGenomeFasta = null;
		int numVariantsToTest = -1;
		// String fileOFSamplesForPairWise = "D:/data/Project_Tsai_Spector_Joint/ErrorRates/pwise.txt";
		// String vcf =
		// "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";

		int numThreads = 8;
		// String logfile = null;
		Logger log;

		String usage = "\n" + "seq.qc.SeqErrorIterator requires 0-1 arguments\n";
		usage += "   (1) vcf filename (i.e. vcf=" + vcf + " (default))\n" + "";
		usage +=
					"   (2) a file of samples for pairwise comparisions (i.e. samples= (no default))\n" + "";
		usage += "   (3) number of variants to test per comparision (i.e. numVar="	+ numVariantsToTest
							+ " (defaults to all variants))\n" + "";
		usage += "   (4) full path to a reference genome file (i.e. ref=(no default))\n" + "";

		usage += PSF.Ext.getNumThreadsCommand(3, numThreads);
		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("vcf=")) {
				vcf = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("ref=")) {
				referenceGenomeFasta = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("samples=")) {
				fileOFSamplesForPairWise = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("numVar=")) {
				numVariantsToTest = ext.parseIntArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(Files.exists(fileOFSamplesForPairWise)	? fileOFSamplesForPairWise + "er.log"
																															: null);
			validate(	vcf, fileOFSamplesForPairWise, referenceGenomeFasta, numVariantsToTest, numThreads,
								log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

//
// public static void test(String vcf, String fileOFSamplesForPairWise, Logger log) {
//
// String[] sampsForPairWise = fileOFSamplesForPairWise == null ? null :
// HashVec.loadFileToStringArray(fileOFSamplesForPairWise, false, new int[] { 0 }, true);
// int numthreads = 8;
// VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
//
// VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
// VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
// VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
// // VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
// // VARIANT_FILTER_BOOLEAN amb = VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER;
// VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;
//
// // VARIANT_FILTER_BOOLEAN[] bQualFilts = new VARIANT_FILTER_BOOLEAN[] { amb };
// VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, dp, gq, vqslod };
// VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new
// VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
//
// VariantContextFilter vContextFilterAnno = null;
// log.reportTimeInfo("Intializing snp138 filter");
// vContextFilterAnno = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new
// VARIANT_FILTER_BOOLEAN[] {}, new String[] { "SNP138" }, new String[] { "snp138=='.'" }, log);
//
// FilterNGS filterNGS = new FilterNGS(0, 0, new int[] { 6 });
// SeqError seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, vContextFilter, filterNGS, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138FilterSummary.e2"));
//
// seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, null, filterNGS, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise,
// ".pwiseSnp138FilterJustAltSummary.e2"));
//
// seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, null, null, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138Summary.e2"));
//
// seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, vContextFilter, filterNGS, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseFilterSummary.e2"));
//
// seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, null, filterNGS, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseFilterJustAltSummary.e2"));
//
// seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, null, null, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSummary.e2"));
//
// }
//
// VARIANT_FILTER_DOUBLE vq = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
// vq.setDFilter(VQSLOD_Values[i]);
// vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { vq }, new
// VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
// index++;
// }
// for (int j = 0; j < GQ_Values.length; j++) {
// VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
// gq.setDFilter(GQ_Values[j]);
// vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { gq }, new
// VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
// index++;
// }
// for (int k = 0; k < DEPTH.length; k++) {
// VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
// dp.setDFilter(DEPTH[k]);
// vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { dp }, new
// VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
// index++;
// }
