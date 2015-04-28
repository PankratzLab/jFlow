package seq.qc;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.concurrent.Callable;

import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import seq.qc.FilterNGS.VcFilterDouble;
import seq.qc.SeqError.DuplicateETwo;
import common.Array;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;

public class SeqQCValidation {
	private static String[] DB_SNP_SETS = new String[] { "snp138=='.'", "snp138!='.'", "(snp138=='.'||snp138!='.')" };

	private FilterNGS filterNGS;
	private VariantContextFilter vContextFilterSample;
	private VariantContextFilter vContextFilterVariant;
	private SeqError seqError;

	public SeqQCValidation(FilterNGS filterNGS, VariantContextFilter vContextFilterSample, VariantContextFilter vContextFilterVariant, SeqError seqError) {
		super();
		this.filterNGS = filterNGS;
		this.vContextFilterSample = vContextFilterSample;
		this.vContextFilterVariant = vContextFilterVariant;
		this.seqError = seqError;
	}

	public SeqError getSeqError() {
		return seqError;
	}

	public FilterNGS getFilterNGS() {
		return filterNGS;
	}

	public VariantContextFilter getvContextFilterSample() {
		return vContextFilterSample;
	}

	public VariantContextFilter getvContextFilterVariant() {
		return vContextFilterVariant;
	}

	public void validate(int numthreads) {
		seqError.populateError(vContextFilterVariant, vContextFilterSample, filterNGS, numthreads);
	}

	private static class SeqQCValidationWorker implements Callable<SeqQCValidation> {
		private SeqQCValidation seqQCValidation;
		private int numthreads;

		public SeqQCValidationWorker(SeqQCValidation seqQCValidation, int numthreads) {
			super();
			this.seqQCValidation = seqQCValidation;
			this.numthreads = numthreads;
		}

		@Override
		public SeqQCValidation call() throws Exception {
			seqQCValidation.validate(numthreads);
			return seqQCValidation;
		}

	}

	private static class SeqQCValidationProducer implements Producer<SeqQCValidation> {

		private SeqQCValidation[] validations;
		private int index;
		private Logger log;

		public SeqQCValidationProducer(SeqQCValidation[] validations, Logger log) {
			super();
			this.validations = validations;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < validations.length;
		}

		@Override
		public Callable<SeqQCValidation> next() {
			SeqQCValidationWorker worker = new SeqQCValidationWorker(validations[index], 1);
			index++;
			return worker;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static VariantContextFilter[] getWholeVariantFilters(Logger log) {
		VariantContextFilter[] vContextFilterVariant = new VariantContextFilter[DB_SNP_SETS.length];
		for (int i = 0; i < vContextFilterVariant.length; i++) {
			vContextFilterVariant[i] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] {}, new String[] { DB_SNP_SETS[i] }, new String[] { DB_SNP_SETS[i] }, log);
		}
		return vContextFilterVariant;
	}

	private static VariantContextFilter[] getSampleVariantFilters(Logger log) {
		double[] GQ_Values = new double[] { -1, 20, 50, 80, 90, 98 };
		double[] VQSLOD_Values = new double[] { -1, 0, 1, 2, 3 };
		double[] DEPTH = new double[] { -1, 0, 10, 20 };
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		VariantContextFilter[] vContextFilterSample = new VariantContextFilter[VQSLOD_Values.length * GQ_Values.length * DEPTH.length];
		int index = 0;
		for (int i = 0; i < VQSLOD_Values.length; i++) {
			for (int j = 0; j < GQ_Values.length; j++) {
				for (int k = 0; k < DEPTH.length; k++) {
					VARIANT_FILTER_DOUBLE vq = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
					vq.setDFilter(VQSLOD_Values[i]);
					VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
					gq.setDFilter(GQ_Values[j]);
					VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
					dp.setDFilter(DEPTH[k]);
					vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { vq, gq, dp }, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
					index++;
				}
			}
		}
		//
		// VARIANT_FILTER_DOUBLE vq = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		// vq.setDFilter(VQSLOD_Values[i]);
		// vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { vq }, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
		// index++;
		// }
		// for (int j = 0; j < GQ_Values.length; j++) {
		// VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
		// gq.setDFilter(GQ_Values[j]);
		// vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { gq }, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
		// index++;
		// }
		// for (int k = 0; k < DEPTH.length; k++) {
		// VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		// dp.setDFilter(DEPTH[k]);
		// vContextFilterSample[index] = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { dp }, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
		// index++;
		// }
		return vContextFilterSample;
	}

	private static FilterNGS[] getAltAlleleDepthFilters(int start, int stop) {
		FilterNGS[] filterNGSs = new FilterNGS[stop + 1 - start];
		int index = 0;
		for (int i = start; i <= stop; i++) {
			filterNGSs[index] = new FilterNGS(0, 0, new int[] { i });
			index++;
		}
		return filterNGSs;
	}

	public static void validate(String vcf, String fileOFSamplesForPairWise, int numthreads, Logger log) {
		String[] sampsForPairWise = fileOFSamplesForPairWise == null ? null : HashVec.loadFileToStringArray(fileOFSamplesForPairWise, false, new int[] { 0 }, true);
		String output = ext.addToRoot(fileOFSamplesForPairWise, ".validation.summary");
		FilterNGS[] altDepthFilterNGSs = getAltAlleleDepthFilters(0, 10);

		VariantContextFilter[] vContextFilterVariants = getWholeVariantFilters(log);
		VariantContextFilter[] vContextFilterSamples = getSampleVariantFilters(log);

		SeqQCValidation[] validations = new SeqQCValidation[altDepthFilterNGSs.length * vContextFilterVariants.length * vContextFilterSamples.length];
		log.reportTimeInfo("Running " + validations.length + " validations ");
		int index = 0;
		for (int i = 0; i < altDepthFilterNGSs.length; i++) {
			for (int j = 0; j < vContextFilterVariants.length; j++) {
				for (int j2 = 0; j2 < vContextFilterSamples.length; j2++) {
					validations[index] = new SeqQCValidation(altDepthFilterNGSs[i], vContextFilterSamples[j2], vContextFilterVariants[j], new SeqError(vcf, sampsForPairWise, log));
					index++;
				}

			}
		}
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println(Array.toStr(SeqError.OUTPUT_HEADER) + "\t" + "SET\tAltDepth\tGQ\tVQSLOD\tDepth");
			SeqQCValidationProducer producer = new SeqQCValidationProducer(validations, log);
			WorkerTrain<SeqQCValidation> train = new WorkerTrain<SeqQCValidation>(producer, numthreads, numthreads, log);
			while (train.hasNext()) {
				SeqQCValidation tmp = train.next();
				DuplicateETwo[] dETwos = tmp.getSeqError().getdETwos();
				for (int i = 0; i < dETwos.length; i++) {
					writer.print(dETwos[i].getSummary() + "\t" + tmp.getvContextFilterVariant().getvFilterJEXL().getjExps().get(0).name + "\t" + tmp.getFilterNGS().getReadDepthFilter()[0]);
					VcFilterDouble[] doubles = tmp.getvContextFilterSample().getvDoubles();
					int GQIndex = -1;
					int VQSLODIndex = -1;
					int dpIndex = -1;
					for (int j = 0; j < doubles.length; j++) {
						if (doubles[j].getDfilter() == VARIANT_FILTER_DOUBLE.GQ) {
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
		// String fileOFSamplesForPairWise = "D:/data/Project_Tsai_Spector_Joint/ErrorRates/pwise.txt";
		// String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";

		int numThreads = 7;
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.qc.SeqErrorIterator requires 0-1 arguments\n";
		usage += "   (1) vcf filename (i.e. vcf=" + vcf + " (default))\n" + "";
		usage += "   (2) a file of samples for pairwise comparisions (i.e. samples= (no default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(3, numThreads);
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("samples=")) {
				fileOFSamplesForPairWise = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			validate(vcf, fileOFSamplesForPairWise, numThreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

//
// public static void test(String vcf, String fileOFSamplesForPairWise, Logger log) {
//
// String[] sampsForPairWise = fileOFSamplesForPairWise == null ? null : HashVec.loadFileToStringArray(fileOFSamplesForPairWise, false, new int[] { 0 }, true);
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
// VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
//
// VariantContextFilter vContextFilterAnno = null;
// log.reportTimeInfo("Intializing snp138 filter");
// vContextFilterAnno = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] {}, new String[] { "SNP138" }, new String[] { "snp138=='.'" }, log);
//
// FilterNGS filterNGS = new FilterNGS(0, 0, new int[] { 6 });
// SeqError seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, vContextFilter, filterNGS, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138FilterSummary.e2"));
//
// seqError = new SeqError(vcf, sampsForPairWise, log);
// seqError.populateError(vContextFilterAnno, null, filterNGS, numthreads);
// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138FilterJustAltSummary.e2"));
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
