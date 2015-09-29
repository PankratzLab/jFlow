package seq.analysis;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.ChrSplitResults;
import seq.manage.VCFOps.HEADER_COPY_TYPE;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import common.Array;
import common.ArraySpecialList;
import common.Files;
import common.Logger;
import common.WorkerHive;
import common.ext;

/**
 *
 *
 */
public class VCFSimpleTally {
	private static final String[] EFF = { "HIGH", "MODERATE", "LOW" };
	private static final String[][] EFF_DEFS = new String[][] { new String[] { EFF[0] }, new String[] { EFF[0], EFF[1] }, new String[] { EFF[0], EFF[1], EFF[2] } };
	private static final String ESP_FILTER = "(esp6500si_all=='.'||esp6500si_all <= 0.01)";
	private static final String G1000_FILTER = "(g10002014oct_all=='.'||g10002014oct_all <= 0.01)";
	private static final String AND = "&&";
	private static final String SNPEFF_IMPACTS = "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE'||SNPEFF_IMPACT=='LOW')";
	private static final String SNPEFF_NAMES = "G1000_esp_charge_aricFreq_SNPEFF_HIGH_MODERATE_LOW";
	private static final String CHARGE_B_FILTER = "(MAF_blacks=='.'||MAF_blacks <= 0.01)";
	private static final String CHARGE_W_FILTER = "(MAF_whites=='.'||MAF_whites <= 0.01)";
	private static final String[] ANNO_BASE = new String[] { "CHROM", "POS", "ID", "REF", "ALT" };
	private static final String[] ANNO_ADD = new String[] { "_NUM_WITH_CALLS", "_NUM_WITH_ALT", "_AAC", "_HQ_NUM_WITH_ALT", "_HQ_AAC", };

	private static void filter(String vcf, String output, String casePop, VcfPopulation vpop, double controlFreq, Logger log) {
		if (!Files.exists(output)) {
			Set<String> cases = vpop.getSuperPop().get(casePop);
			Hashtable<String, Set<String>> controls = vpop.getSuperPop();

			log.reportTimeInfo("CASE :" + casePop + " n: " + cases.size());
			controls.remove(casePop);
			controls.remove(VcfPopulation.EXCLUDE);
			for (String control : controls.keySet()) {
				log.reportTimeInfo("Control: " + control + " n: " + controls.get(control).size());
			}
			VCFFileReader reader = new VCFFileReader(vcf, true);
			VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
			VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
			VariantContextFilter freqFilter = getFreqFilter(log);
			int numScanned = 0;
			int numPass = 0;

			for (VariantContext vc : reader) {
				numScanned++;
				if (numScanned % 10000 == 0) {
					log.reportTimeInfo(numScanned + " variants scanned, " + numPass + " variants passed");
				}
				if (!vc.isFiltered() && vc.isBiallelic()) {// no tranche
					VariantContext vcCase = VCOps.getSubset(vc, cases, VC_SUBSET_TYPE.SUBSET_STRICT, false);
					if (vcCase.getSampleNames().size() != cases.size()) {
						throw new IllegalArgumentException("could not find all cases for " + casePop);
					}
					if (vcCase.getHomRefCount() != cases.size() && vcCase.getNoCallCount() != cases.size() && freqFilter.filter(vcCase).passed()) {// as alts in rare esp/1000g
						boolean controlPass = true;
						for (String controlPop : controls.keySet()) {
							VariantContext vcControl = VCOps.getSubset(vc, controls.get(controlPop), VC_SUBSET_TYPE.SUBSET_STRICT, false);
							if (vcControl.getSampleNames().size() != controls.get(controlPop).size()) {
								throw new IllegalArgumentException("could not find all controls for " + controlPop);
							}
							double maf = VCOps.getMAF(vcControl, null);
							if (maf >= controlFreq && vcControl.getNoCallCount() != controls.get(controlPop).size()) {// rare in control
								controlPass = false;
								break;
							}
						}
						if (controlPass) {
							numPass++;
							writer.add(vc);
						}
					}
				}
			}
			reader.close();
			writer.close();
		} else {
			log.reportFileExists(output);
		}
	}

	private static VariantContextFilter getFreqFilter(Logger log) {
		VariantContextFilter vContextFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] {}, new String[] { "RARE_" + SNPEFF_NAMES }, new String[] { ESP_FILTER + AND + G1000_FILTER + AND + SNPEFF_IMPACTS + AND + CHARGE_B_FILTER + AND + CHARGE_W_FILTER }, log);
		return vContextFilter;
	}

	public static void test() {
		String vcf = "/home/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/aric_merge/vcf/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.aric.chargeMaf.vcf";
		String popDir = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/aric_merge/vcf/Freq/";
		String[] vpops = new String[] { popDir + "ANIRDIA.vpop", popDir + "ANOTIA.vpop", popDir + "CUSHING.vpop" };

		double maf = 0.01;
		int numThreads = 24;
		for (int i = 0; i < vpops.length; i++) {
			String outDir = ext.parseDirectoryOfFile(vpops[i]);
			new File(outDir).mkdirs();
			Logger log = new Logger(ext.rootOf(vpops[i], false) + ".log");
			runSimpleTally(vcf, vpops[i], maf, numThreads, outDir, log);
		}
	}

	private static void runSimpleTally(String vcf, String vpop, double maf, int numThreads, String outDir, Logger log) {
		VcfPopulation vpopAc = VcfPopulation.load(vpop, POPULATION_TYPE.ANY, log);
		vpopAc.report();
		String caseDef = ext.rootOf(vpop);
		ChrSplitResults[] chrSplitResults = VCFOps.splitByChrs(vcf, outDir, numThreads, false, log);
		WorkerHive<String> hive = new WorkerHive<String>(numThreads, 10, log);
		ArrayList<String> filtVcfs = new ArrayList<String>();

		for (int i = 0; i < chrSplitResults.length; i++) {// filter each chromosome
			String outputVcf = outDir + ext.rootOf(vpop) + chrSplitResults[i].getChr() + VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral();
			filtVcfs.add(outputVcf);
			FilterWorker worker = new FilterWorker(chrSplitResults[i].getOutputVCF(), outputVcf, ext.rootOf(vpop), VcfPopulation.load(vpop, POPULATION_TYPE.ANY, log), maf, log);
			hive.addCallable(worker);
		}
		hive.execute(true);

		String finalOut = outDir + ext.rootOf(vpop) + ".final";
		String finalOutVCF = finalOut + VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral();
		String finalAnnot = finalOut + ".summary";
		VCFFileReader tmp = new VCFFileReader(filtVcfs.get(0), true);

		VariantContextWriter writer = VCFOps.initWriter(finalOutVCF, VCFOps.DEFUALT_WRITER_OPTIONS, VCFOps.getSequenceDictionary(tmp));
		VCFOps.copyHeader(tmp, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
		tmp.close();

		PrintWriter annoWriter = Files.getAppropriateWriter(finalAnnot);
		annoWriter.print(Array.toStr(ANNO_BASE));
		Set<String> cases = vpopAc.getSuperPop().get(caseDef);
		annoWriter.print("\t" + Array.tagOn(ANNO_ADD, caseDef, null));
		Hashtable<String, Set<String>> controls = vpopAc.getSuperPop();
		controls.remove(caseDef);
		controls.remove(VcfPopulation.EXCLUDE);
		ArrayList<String> controlsOrdered = new ArrayList<String>();
		controlsOrdered.addAll(controls.keySet());
		for (int i = 0; i < controlsOrdered.size(); i++) {
			annoWriter.print("\t" + Array.tagOn(ANNO_ADD, controlsOrdered.get(i), null));
		}
		String[][] annotations = VCFOps.getAnnotationKeys(vcf, log);
		annoWriter.print("\t" + Array.toStr(annotations[0]));

		VariantContextFilter qual = getQualityFilterwkggseq(log);
		Hashtable<String, ArrayList<GeneSummary[]>> geneSummaries = new Hashtable<String, ArrayList<GeneSummary[]>>();
		for (int i = 0; i < filtVcfs.size(); i++) {
			VCFFileReader result = new VCFFileReader(filtVcfs.get(i), true);
			for (VariantContext vc : result) {
				writer.add(vc);
				String geneName = VCOps.getSNP_EFFGeneName(vc);
				if (!geneSummaries.containsKey(geneName)) {
					addEntries(caseDef, controlsOrdered, geneSummaries, geneName);
				}

				VcGroupSummary vcCaseGroup = new VcGroupSummary(caseDef, cases, vc, qual, log);
				for (int j = 0; j < geneSummaries.get(geneName).get(0).length; j++) {
					geneSummaries.get(geneName).get(0)[j].add(vcCaseGroup);
				}
				annoWriter.print(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getReference().getBaseString() + "\t" + vc.getAlternateAlleles().toString());
				annoWriter.print("\t" + Array.toStr(vcCaseGroup.getSummary()));
				for (int j = 0; j < controlsOrdered.size(); j++) {
					VcGroupSummary vcControlGroup = new VcGroupSummary(controlsOrdered.get(i), controls.get(controlsOrdered.get(i)), vc, qual, log);
					for (int k = 0; k < geneSummaries.get(geneName).get(0).length; k++) {
						geneSummaries.get(geneName).get(0)[k].add(vcControlGroup);
					}
					annoWriter.print("\t" + Array.toStr(vcControlGroup.getSummary()));
				}
			}
			result.close();
		}
		annoWriter.close();
		writer.close();
	}

	private static void addEntries(String caseDef, ArrayList<String> controlsOrdered, Hashtable<String, ArrayList<GeneSummary[]>> geneSummaries, String geneName) {
		geneSummaries.put(geneName, new ArrayList<GeneSummary[]>());
		GeneSummary[] caseGeneSummaries = new GeneSummary[EFF_DEFS.length];
		for (int j = 0; j < caseGeneSummaries.length; j++) {
			caseGeneSummaries[j] = new GeneSummary(geneName, caseDef, EFF_DEFS[j]);
		}
		geneSummaries.get(geneName).add(caseGeneSummaries);
		for (int j = 0; j < controlsOrdered.size(); j++) {
			GeneSummary[] controlGeneSummaries = new GeneSummary[EFF_DEFS.length];
			for (int k = 0; k < controlGeneSummaries.length; k++) {
				controlGeneSummaries[k] = new GeneSummary(geneName, controlsOrdered.get(j), EFF_DEFS[k]);
			}
		}
	}

	private static class GeneSummary {
		private static final String SNPEFF_IMPACT = "SNPEFF_IMPACT";
		private String geneName;
		private String group;
		private String[] effects;
		private int numVar;
		private HashSet<String> uniqueIndsWithVar;
		private int hqNumVar;
		private HashSet<String> uniqueHqIndsWithVar;

		private GeneSummary(String geneName, String group, String[] effects) {
			super();
			this.group = group;
			this.geneName = geneName;
			this.numVar = 0;
			this.hqNumVar = 0;
			this.uniqueHqIndsWithVar = new HashSet<String>();
			this.uniqueIndsWithVar = new HashSet<String>();
		}

		private void add(VcGroupSummary vcGroupSummary) {
			if (!vcGroupSummary.getGroupName().equals(group)) {
				throw new IllegalArgumentException("Mismatched group names");
			}
			if (!VCOps.getSNP_EFFGeneName(vcGroupSummary.getVcOriginal()).equals(geneName)) {
				throw new IllegalArgumentException("Mismatched gene names");
			} else {
				String impact = VCOps.getSNP_EFFImpact(vcGroupSummary.getVcOriginal());
				if (ext.indexOfStr(impact, effects) >= 0) {
					if (vcGroupSummary.getIndsWithAlt().size() > 0) {
						numVar++;
						uniqueIndsWithVar.addAll(vcGroupSummary.getIndsWithAlt());
						if (vcGroupSummary.getHqIndsWithAlt().size() > 0) {
							hqNumVar++;
							uniqueIndsWithVar.addAll(uniqueHqIndsWithVar);
						}
					}
				}
			}
		}
	}

	private static class VcGroupSummary {
		private String groupName;
		private Set<String> group;
		private VariantContext vcOriginal;
		private VariantContext vcAlt;
		private VariantContext vcAltHq;

		private int numIndsAlt;
		private int numWithCalls;
		private int aac;
		private int numHqIndsAlt;
		private int hqAac;
		private Set<String> indsWithAlt;
		private Set<String> hqIndsWithAlt;
		private VariantContextFilter filter;
		private Logger log;

		public VcGroupSummary(String groupName, Set<String> group, VariantContext vcOriginal, VariantContextFilter filter, Logger log) {
			super();
			this.groupName = groupName;
			this.group = group;
			this.vcOriginal = vcOriginal;
			this.filter = filter;
			this.log = log;
			summarize();
		}

		public String getGroupName() {
			return groupName;
		}

		public VariantContext getVcOriginal() {
			return vcOriginal;
		}

		public Set<String> getIndsWithAlt() {
			return indsWithAlt;
		}

		public Set<String> getHqIndsWithAlt() {
			return hqIndsWithAlt;
		}

		private void summarize() {
			VariantContext vcSub = VCOps.getSubset(vcOriginal, group, VC_SUBSET_TYPE.SUBSET_STRICT, false);
			this.vcAlt = VCOps.getAltAlleleContext(vcOriginal, null, filter, ALT_ALLELE_CONTEXT_TYPE.ALL, log);
			this.aac = (int) VCOps.getAAC(vcSub, null);
			this.indsWithAlt = vcAlt.getSampleNames();
			this.numIndsAlt = indsWithAlt.size();
			this.numWithCalls = vcSub.getSampleNames().size() - vcSub.getNoCallCount();
			this.vcAltHq = VCOps.getAltAlleleContext(vcOriginal, null, filter, ALT_ALLELE_CONTEXT_TYPE.ALL, log);
			this.hqAac = (int) VCOps.getAAC(vcAltHq, null);
			this.hqIndsWithAlt = vcAltHq.getSampleNames();
			this.numHqIndsAlt = hqIndsWithAlt.size();
		}

		private String[] getSummary() {
			ArrayList<String> summary = new ArrayList<String>();
			summary.add(numWithCalls + "");
			summary.add(numIndsAlt + "");
			summary.add(aac + "");
			summary.add(numHqIndsAlt + "");
			summary.add(hqAac + "");
			return Array.toStringArray(summary);
		}
	}

	private static class FilterWorker implements Callable<String> {
		private String vcf;
		private String outputVcf;
		private String casePop;
		private VcfPopulation vpop;
		private double maf;
		private Logger log;

		public FilterWorker(String vcf, String outputVcf, String casePop, VcfPopulation vpop, double maf, Logger log) {
			super();
			this.vcf = vcf;
			this.outputVcf = outputVcf;
			this.casePop = casePop;
			this.vpop = vpop;
			this.maf = maf;
			this.log = log;
		}

		@Override
		public String call() throws Exception {
			filter(vcf, outputVcf, casePop, vpop, maf, log);

			return outputVcf;
		}

	}

	private static VariantContextFilter getQualityFilterwkggseq(Logger log) {
		VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_LOOSE;
		gq.setDFilter(20);
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		dp.setDFilter(4);
		VARIANT_FILTER_DOUBLE altD = VARIANT_FILTER_DOUBLE.ALT_ALLELE_DEPTH;
		altD.setDFilter(1);
		VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		// VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
		// VARIANT_FILTER_BOOLEAN amb = VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER;
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		// VARIANT_FILTER_BOOLEAN[] bQualFilts = new VARIANT_FILTER_BOOLEAN[] { amb };
		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, dp, altD, gq, vqslod };
		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, new String[] { "G1000Freq" }, new String[] { ESP_FILTER + AND + G1000_FILTER }, log);
		return vContextFilter;
	}

	public static void main(String[] args) {
		test();
	}
}
