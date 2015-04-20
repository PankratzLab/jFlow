package seq.qc;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Callable;

import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import common.Array;
import common.HashVec;
import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;

/**
 * @author lane0212 Class to compute error rates between samples akin to http://onlinelibrary.wiley.com/doi/10.1002/gepi.21881/epdf
 *
 */
public class SeqError {
	private static final String[] OUTPUT_HEADER = new String[] { "Samp1", "Samp2", "total", "allMismatched", "allNonMissingMismatched", "proportionAgreeTotal", "proportionAgreeNonMissing" };
	private static final String SNP138 = "snp138";
	private String vcfFile;
	private String[] sampsForPairWise;
	private DuplicateETwo[] dETwos;
	private Logger log;

	public SeqError(String vcfFile, String[] sampsForPairWise, Logger log) {
		super();
		this.vcfFile = vcfFile;
		this.sampsForPairWise = sampsForPairWise == null ? VCFOps.getSamplesInFile(new VCFFileReader(vcfFile, true)) : sampsForPairWise;
		this.dETwos = getDups(this.sampsForPairWise, log);
		this.log = log;
	}

	public void populateError(boolean skipSNP138, VariantContextFilter vContextFilter2, FilterNGS filterNGS, int numthreads) {
		VCFFileReader reader = new VCFFileReader(vcfFile, true);
		VCFHeader header = reader.getFileHeader();

		log.reportTimeInfo("Computing pairwise comparisons for " + sampsForPairWise.length + " samples (" + dETwos.length + ") comparisons per variant");
		VariantContextFilter vContextFilter = null;
		if (skipSNP138) {
			log.reportTimeInfo("Intializing snp138 filter");
			vContextFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] {}, new String[] { "SNP138" }, new String[] { "snp138=='.'" }, log);
		}
		if (skipSNP138 && !VCFOps.hasInfoLine(reader, SNP138)) {
			log.reportTimeError("Cannot find annotation for " + SNP138 + " and skipping known sites was flagged");
		} else {
			int index = 0;
			WorkerTrain<DuplicateETwo> train = new WorkerTrain<SeqError.DuplicateETwo>(null, numthreads, numthreads, log);
			train.setAutoShutDown(false);
			long time = System.currentTimeMillis();
			for (VariantContext vc : reader) {

				if (vContextFilter == null || vContextFilter.filter(vc).passed()) {
					index++;

					if (index % 400000 == 0) {
						log.reportTimeInfo(index + " variants processed..." + ext.getTimeElapsed(time));
						time = System.currentTimeMillis();
						reader.close();
						train.shutdown();
						return;
					}
					vc.fullyDecode(header, false);
					DuplicateProducer producer = new DuplicateProducer(vc, vContextFilter2, dETwos, filterNGS);
					train.setProducer(producer);
					int tmpI = 0;
					while (train.hasNext()) {
						dETwos[tmpI] = train.next();
						tmpI++;
					}
				}

			}
			reader.close();
			train.shutdown();
		}
	}

	public void summarize(String fullPathToOutput) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.println(Array.toStr(OUTPUT_HEADER));
			for (int i = 0; i < dETwos.length; i++) {
				writer.println(dETwos[i].getSummary());
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathToOutput);
			log.reportException(e);
		}
	}

	private static class DuplicateProducer implements Producer<DuplicateETwo> {
		private VariantContext vc;
		private DuplicateETwo[] dETwos;
		private int index;
		private VariantContextFilter variantContextFilter;
		private FilterNGS filterNGS;

		public DuplicateProducer(final VariantContext vc, final VariantContextFilter variantContextFilter, final DuplicateETwo[] dETwos, FilterNGS filterNGS) {
			super();
			this.vc = vc;
			this.dETwos = dETwos;
			this.index = 0;
			this.variantContextFilter = variantContextFilter;
			this.filterNGS = filterNGS;
		}

		@Override
		public boolean hasNext() {
			return index < dETwos.length;
		}

		@Override
		public Callable<DuplicateETwo> next() {
			DuplicateWorker worker = new DuplicateWorker(vc, dETwos[index], variantContextFilter, filterNGS);
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

	private static class DuplicateWorker implements Callable<DuplicateETwo> {
		private VariantContext vc;
		private DuplicateETwo deTwo;
		private VariantContextFilter variantContextFilter;
		private FilterNGS filterNGS;

		public DuplicateWorker(VariantContext vc, DuplicateETwo deTwo, VariantContextFilter variantContextFilter, FilterNGS filterNGS) {
			super();
			this.vc = vc;
			this.deTwo = deTwo;
			this.variantContextFilter = variantContextFilter;
			this.filterNGS = filterNGS;
		}

		@Override
		public DuplicateETwo call() throws Exception {
			deTwo.addVC(vc, variantContextFilter, filterNGS);
			return deTwo;
		}
	}

	private static DuplicateETwo[] getDups(String[] allSamps, Logger log) {
		DuplicateETwo[] dETwos = new DuplicateETwo[(allSamps.length * (allSamps.length - 1)) / 2];
		int index = 0;
		for (int i = 0; i < allSamps.length; i++) {
			for (int j = i + 1; j < allSamps.length; j++) {
				if (j == i) {
				}
				HashSet<String> curDups = new HashSet<String>();
				curDups.add(allSamps[i]);
				curDups.add(allSamps[j]);

				dETwos[index] = new DuplicateETwo(curDups, log);
				index++;

			}
		}
		return dETwos;
	}

	private static class DuplicateETwo {
		private Set<String> dups;
		private int allMismatched;
		private int nonMissingMismatched;
		private int total;
		private Logger log;

		public DuplicateETwo(Set<String> dups, Logger log) {
			super();
			this.dups = dups;
		}

		public String getSummary() {
			String summary = "";
			int index = 0;
			for (String dup : dups) {
				summary += (index == 0 ? "" : "\t") + dup;
				index++;
			}

			double mm6 = (double) allMismatched / total;
			double mm7 = (double) nonMissingMismatched / total;
			mm6 = 1 - mm6;
			mm7 = 1 - mm7;

			summary += "\t" + total;
			summary += "\t" + allMismatched;
			summary += "\t" + nonMissingMismatched;
			summary += "\t" + mm6;
			summary += "\t" + mm7;
			return summary;
		}

		public void addVC(VariantContext vc, VariantContextFilter variantContextFilter, FilterNGS filterNGS) {
			VariantContext vcSub = VCOps.getSubset(vc, dups, VC_SUBSET_TYPE.SUBSET_STRICT);
			if (variantContextFilter == null || variantContextFilter.filter(vcSub).passed()) {
				if (VCOps.getAltAlleleContext(vcSub, filterNGS == null || filterNGS.getReadDepthFilter() == null ? 0 : filterNGS.getReadDepthFilter()[0]).getSampleNames().size() > 0) {// has alt call
					GenotypesContext gc = vcSub.getGenotypes();
					if (vcSub.getNoCallCount() < dups.size()) {
						total++;
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
						if (allCalled && !allMatched) {
							nonMissingMismatched++;
						}
						if (!allMatched) {
							allMismatched++;
						}
					}
				}
			}
		}
	}

	public static void test(String vcf, Logger log) {
		String fileOFSamplesForPairWise = "D:/data/Project_Tsai_Spector_Joint/ErrorRates/pwise.txt";
		String[] sampsForPairWise = fileOFSamplesForPairWise == null ? null : HashVec.loadFileToStringArray(fileOFSamplesForPairWise, false, new int[] { 0 }, true);
		int numthreads = 8;
		VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		// VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
		// VARIANT_FILTER_BOOLEAN amb = VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER;
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		// VARIANT_FILTER_BOOLEAN[] bQualFilts = new VARIANT_FILTER_BOOLEAN[] { amb };
		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, dp, gq, vqslod };
		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);

		FilterNGS filterNGS = new FilterNGS(0, 0, new int[] { 6 });
		SeqError seqError = new SeqError(vcf, sampsForPairWise, log);
		seqError.populateError(true, vContextFilter, filterNGS, numthreads);
		seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138FilterSummary.e2"));

		seqError = new SeqError(vcf, sampsForPairWise, log);
		seqError.populateError(true, null, filterNGS, numthreads);
		seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138FilterJustAltSummary.e2"));

		seqError = new SeqError(vcf, sampsForPairWise, log);
		seqError.populateError(true, null, null, numthreads);
		seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSnp138Summary.e2"));

		seqError = new SeqError(vcf, sampsForPairWise, log);
		seqError.populateError(false, vContextFilter, filterNGS, numthreads);
		seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseFilterSummary.e2"));

		seqError = new SeqError(vcf, sampsForPairWise, log);
		seqError.populateError(false, null, filterNGS, numthreads);
		seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseFilterJustAltSummary.e2"));

		seqError = new SeqError(vcf, sampsForPairWise, log);
		seqError.populateError(false, null, null, numthreads);
		seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSummary.e2"));

		// seqError = new SeqError(vcf, sampsForPairWise, log);
		// seqError.populateError(false, vContextFilter, filterNGS, 5);
		// seqError.summarize(ext.addToRoot(fileOFSamplesForPairWise, ".pwiseSummary"));

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.qc.SeqError requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + vcf + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			log = new Logger(ext.rootOf(vcf) + ".log");
			test(vcf, log);
			// parse(filename, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
