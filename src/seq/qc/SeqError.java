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
import seq.qc.FilterNGS.VariantContextFilter;
import common.Array;
import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;

/**
 * @author lane0212 Class to compute error rates between samples akin to http://onlinelibrary.wiley.com/doi/10.1002/gepi.21881/epdf
 *
 */
public class SeqError {
	public static final String[] OUTPUT_HEADER = new String[] { "Samp1", "Samp2", "total", "allMismatched", "allNonMissingMismatched", "proportionAgreeTotal", "proportionAgreeNonMissing" };
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

	public DuplicateETwo[] getdETwos() {
		return dETwos;
	}

	/**
	 * @param vContextFilterVariant
	 *            filter applied to the entire variant context
	 * @param vContextFilterSample
	 *            filter applied to the pairwise comparison
	 * @param filterNGS
	 * @param numthreads
	 *            if a large number of pairwise comparisons are performed, the number of threads can be increased
	 */
	public void populateError(VariantContextFilter vContextFilterVariant, VariantContextFilter vContextFilterSample, FilterNGS filterNGS, int numVariantsToTest, int numthreads) {
		VCFFileReader reader = new VCFFileReader(vcfFile, true);
		VCFHeader header = reader.getFileHeader();

		log.reportTimeInfo("Computing pairwise comparisons for " + sampsForPairWise.length + " samples (" + dETwos.length + ") comparisons per variant");
		int index = 0;
		WorkerTrain<DuplicateETwo> train = new WorkerTrain<SeqError.DuplicateETwo>(null, numthreads, numthreads, log);
		train.setAutoShutDown(false);
		long time = System.currentTimeMillis();
		for (VariantContext vc : reader) {

			if (vContextFilterVariant == null || vContextFilterVariant.filter(vc).passed()) {
				index++;

				if (index % 100000 == 0) {
					log.reportTimeInfo(index + " variants processed..." + ext.getTimeElapsed(time));
					time = System.currentTimeMillis();
					// reader.close();
					// train.shutdown();
					// return;
				}
				if (numVariantsToTest >= 0 && index == numVariantsToTest) {
					log.reportTimeInfo(index + " variants processed...," + numVariantsToTest + " variants to test reached " + ext.getTimeElapsed(time));
					reader.close();
					train.shutdown();
					return;
				}
				vc.fullyDecode(header, false);
				DuplicateProducer producer = new DuplicateProducer(vc, vContextFilterSample, dETwos, filterNGS);
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

		private DuplicateProducer(final VariantContext vc, final VariantContextFilter variantContextFilter, final DuplicateETwo[] dETwos, FilterNGS filterNGS) {
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

		private DuplicateWorker(VariantContext vc, DuplicateETwo deTwo, VariantContextFilter variantContextFilter, FilterNGS filterNGS) {
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
				HashSet<String> curDups = new HashSet<String>();
				curDups.add(allSamps[i]);
				curDups.add(allSamps[j]);
				dETwos[index] = new DuplicateETwo(curDups, log);
				index++;
			}
		}
		return dETwos;
	}

	/**
	 * Class for tracking the e2 error rate across all comparisons
	 *
	 */
	public static class DuplicateETwo {
		private Set<String> dups;
		private int allMismatched;
		private int nonMissingMismatched;
		private int total;

		// private Logger log;

		private DuplicateETwo(Set<String> dups, Logger log) {
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

		private void addVC(VariantContext vc, VariantContextFilter variantContextFilter, FilterNGS filterNGS) {
			VariantContext vcSub = VCOps.getSubset(vc, dups, VC_SUBSET_TYPE.SUBSET_STRICT);
			if (variantContextFilter == null || variantContextFilter.filter(vcSub).passed()) {
				int altAlleleDepth = filterNGS == null || filterNGS.getReadDepthFilter() == null ? 0 : filterNGS.getReadDepthFilter()[0];
				if (VCOps.getAltAlleleContext(vcSub, altAlleleDepth,variantContextFilter.getLog()).getSampleNames().size() > 0) {// has alt call
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
}
