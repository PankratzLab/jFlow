package seq.analysis;

import filesys.Segment;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import bioinformatics.OMIM;
import bioinformatics.OMIM.OMIMGene;
import seq.analysis.VCFSimpleTally.GeneVariantPositionSummary.ADD_TYPE;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.ChrSplitResults;
import seq.manage.VCFOps.HEADER_COPY_TYPE;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import seq.manage.VCOps.GENOTYPE_INFO;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import common.Array;
import common.Files;
import common.Logger;
import common.Sort;
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
	// private static final String CHARGE_B_FILTER = "(charge.MAF_blacks=='.'||charge.MAF_blacks <= 0.01)";
	// private static final String CHARGE_W_FILTER = "(charge.MAF_whites=='.'||charge.MAF_whites <= 0.01)";
	private static final String[] ANNO_BASE = new String[] { "CHROM", "POS", "ID", "REF", "ALT" };
	private static final String[] ANNO_BASE_SAMPLE = Array.concatAll(ANNO_BASE, new String[] { "SAMPLE", "GENOTYPE" });

	private static final String[] ANNO_ADD = new String[] { "_AVG_GQ", "_AVG_DP", "_NUM_WITH_CALLS", "_NUM_WITH_ALT", "_AAC", "_HQ_NUM_WITH_ALT", "_HQ_AAC", };
	private static final String[] GENE_BASE = new String[] { "GENE/GENE_SET", "FUNCTIONAL_TYPE" };
	private static final String[] GENE_ADD = new String[] { "numVar", "uniqInds", "numCompoundHets", "numCompoundHetsDiffHaplotype", "hqNumVar", "hqUniqInds", "numHQCompoundHets", "numHQCompoundHetsDiffHaplotype" };

	private static boolean filterCHARGE(VariantContext vc, double maf) {
		boolean pass = true;
		if (vc.hasAttribute("charge.MAF_whites")) {
			pass = vc.getCommonInfo().getAttributeAsDouble("charge.MAF_whites", 0) < maf;
		}
		if (pass && vc.hasAttribute("charge.MAF_blacks")) {
			pass = vc.getCommonInfo().getAttributeAsDouble("charge.MAF_blacks", 0) < maf;
		}
		return pass;
	}

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

				// System.out.println(Array.toStr(VCOps.getAnnotationsFor(new String[] { "charge.MAF_whites", "charge.MAF_blacks" }, vc, ".")));

				if (!vc.isFiltered() && vc.isBiallelic()) {// no tranche
					VariantContext vcCase = VCOps.getSubset(vc, cases, VC_SUBSET_TYPE.SUBSET_STRICT, false);
					if (vcCase.getSampleNames().size() != cases.size()) {
						throw new IllegalArgumentException("could not find all cases for " + casePop);
					}
					if (!vcCase.isMonomorphicInSamples() && vcCase.getNoCallCount() != cases.size() && (!vc.hasAttribute("esp6500si_all") || !vc.hasAttribute("g10002014oct_all"))) {
						// String error = "Expected annotations esp6500si_all, g10002014oct_all were not present";
						// error += "\n" + vc.toStringWithoutGenotypes();
						// throw new IllegalStateException(error);
					} else if (vcCase.getHomRefCount() + vcCase.getNoCallCount() != cases.size() && vcCase.getNoCallCount() != cases.size() && freqFilter.filter(vcCase).passed() && filterCHARGE(vcCase, 0.01)) {// as alts in rare esp/1000g
						boolean controlPass = true;
						for (String controlPop : controls.keySet()) {
							VariantContext vcControl = VCOps.getSubset(vc, controls.get(controlPop), VC_SUBSET_TYPE.SUBSET_STRICT, false);
							if (vcControl.getSampleNames().size() != controls.get(controlPop).size()) {
								throw new IllegalArgumentException("could not find all controls for " + controlPop);
							}
							double maf = VCOps.getMAF(vcControl, null);
							if ((!VCOps.isMinorAlleleAlternate(vcControl, null) || maf >= controlFreq) && vcControl.getNoCallCount() != controls.get(controlPop).size()) {// rare in control
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
		// CHARGE_B_FILTER + AND + CHARGE_W_FILTER
		VariantContextFilter vContextFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] {}, new String[] { "RARE_" + SNPEFF_NAMES }, new String[] { ESP_FILTER + AND + G1000_FILTER + AND + SNPEFF_IMPACTS }, log);
		return vContextFilter;
	}

	private static Hashtable<String, String[]> loadToGeneFuncHash(String geneFile, Logger log) {
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		try {
			BufferedReader reader = Files.getAppropriateReader(geneFile);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				String gene = line[0];
				String func = line[1];
				String[] rest = Array.subArray(line, 2);
				hash.put(gene + "_" + func, rest);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + geneFile + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + geneFile + "\"");
			return null;
		}
		return hash;

	}

	private static class GeneSet {
		private String fileName;
		private Logger log;
		private Hashtable<String, String> genes;
		private String tag;

		public GeneSet(String fileName, String tag, Logger log) {
			super();
			this.fileName = fileName;
			this.tag = tag;
			this.log = log;
		}

		public String getTag() {
			return tag;
		}

		private static GeneSet[] load(String[] setFiles, String vpop, Logger log) {
			GeneSet[] sets = new GeneSet[setFiles.length];
			for (int i = 0; i < sets.length; i++) {
				sets[i] = new GeneSet(setFiles[i], ext.rootOf(setFiles[i]).replaceAll("_" + vpop, ""), log);
				sets[i].load();
			}
			return sets;
		}

		public Hashtable<String, String> getGenes() {
			return genes;
		}

		private void load() {
			this.genes = new Hashtable<String, String>();
			try {
				BufferedReader reader = Files.getAppropriateReader(fileName);

				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					String gene = line[0];
					if (genes.contains(line[0])) {
						log.reportTimeWarning(gene + " was sen twice, using first annotation");

						// throw new IllegalAccessError("One entry per gene per file");

					} else {
						genes.put(gene.toUpperCase(), line[0]);
					}

				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + fileName + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + fileName + "\"");
				return;
			}
		}

	}

	private static class OtherGeneInfo {
		private static final String GENE_TAG = "GENE";
		private String fileName;
		private Logger log;
		private Hashtable<String, String[]> anno;
		private String[] header;

		private OtherGeneInfo(String fileName, Logger log) {
			super();

			this.fileName = fileName;
			this.log = log;
		}

		public Hashtable<String, String[]> getAnno() {
			return anno;
		}

		public String[] getHeader() {
			return header;
		}

		private void load() {
			this.anno = new Hashtable<String, String[]>();
			try {
				BufferedReader reader = Files.getAppropriateReader(fileName);
				String[] head = reader.readLine().trim().split("\t");
				if (!head[0].equals(GENE_TAG)) {
					throw new IllegalArgumentException("File " + fileName + " must have " + GENE_TAG + " to be used for extra gene annotation, found " + head[0] + " instead");
				}

				this.header = Array.tagOn(head, ext.removeDirectoryInfo(fileName), null);
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					String gene = line[0];
					if (anno.contains(gene)) {
						throw new IllegalAccessError("One entry per gene per file");

					} else {
						anno.put(gene.toUpperCase(), line);
					}

				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + fileName + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + fileName + "\"");
				return;
			}
		}

	}

	private static Hashtable<String, GeneVariantPositionSummary> parseAvailable(GeneVariantPositionSummary[] summaries) {
		Hashtable<String, GeneVariantPositionSummary> parse = new Hashtable<String, GeneVariantPositionSummary>();
		for (int i = 0; i < summaries.length; i++) {
			String key = summaries[i].getKey();

			if (parse.containsKey(key)) {
				throw new IllegalArgumentException("Dup keys");
			} else {
				parse.put(key, summaries[i]);
			}
		}
		return parse;
	}

	private static ArrayList<Integer> removeNeg(ArrayList<Integer> al) {
		ArrayList<Integer> pos = new ArrayList<Integer>();
		for (int i = 0; i < al.size(); i++) {
			if (al.get(i) >= 0) {
				pos.add(al.get(i));
			}
		}
		return pos;
	}

	private static double centroid(ArrayList<Integer> al) {
		return Array.mean(Array.toIntArray(al));
	}

	// private static double distanceSED(double cent, ArrayList<Integer> al) {
	// double sum = 0;
	// for (int i = 0; i < al.size(); i++) {
	// sum += Math.pow((double) al.get(i) - cent, 2);
	// }
	// return Math.sqrt(sum);
	// }
	private static class DistanceResult {
		private double distance;
		private int iter;

		public DistanceResult(double distance, int iter) {
			super();
			this.distance = distance;
			this.iter = iter;
		}

		public double getDistance() {
			return distance;
		}

		public int getIter() {
			return iter;
		}

	}

	private static DistanceResult distance(ArrayList<Integer> from, ArrayList<Integer> to, boolean identical) {
		double sum = 0;
		int num = 0;
		for (int i = 0; i < from.size(); i++) {
			for (int j = identical ? i + 1 : 0; j < to.size(); j++) {
				sum += (double) Math.abs(from.get(i) - to.get(j));
				num++;
			}
		}
		double distance = (double) sum / num;
		return new DistanceResult(distance, num);
	}

	// private static double distance(double cent, ArrayList<Integer> al) {
	// double sum = 0;
	// for (int i = 0; i < al.size(); i++) {
	// sum += Math.abs((double) (al.get(i) - cent));
	// }
	// return sum / al.size();
	// }

	private static class PosCluster {
		private static final String[] BASE = new String[] { "_AA_CENTROID", "_AVG_DISTANCE", "_N_COMP" };
		private static final String[] OTHER = new String[] { "_DIS_FROM_", "_N_COMP_" };
		private static final String[] NN = new String[] { "_NumNearestTo_" };

		// private static final String[] FINAL_METRIC = new String[] { "FINAL_METRIC" };
		private static final String MEMBER = "member";
		private static final String NONMEMBER = "non";
		private double centriodBelong;
		private DistanceResult distanceBelong;
		private double centriodOut;
		private DistanceResult distanceOut;
		private DistanceResult distanceBelongOut;
		private DistanceResult distanceOutBelong;
		private int nnCase;
		private int nnControl;

		private ArrayList<Integer> al;

		public PosCluster(ArrayList<Integer> al) {
			super();
			this.al = removeNeg(al);
			this.distanceOut = new DistanceResult(Double.NaN, 0);
			this.distanceOutBelong = new DistanceResult(Double.NaN, 0);
			this.distanceBelongOut = new DistanceResult(Double.NaN, 0);
			init();
		}

		private void init() {
			if (al.size() > 0) {
				this.centriodBelong = centroid(al);
				this.distanceBelong = distance(al, al, true);
			} else {
				this.centriodBelong = Double.NaN;
				this.distanceBelong = new DistanceResult(Double.NaN, 0);
			}
		}

		private void computeBelongs(ArrayList<Integer> others) throws IllegalStateException {
			others = removeNeg(others);
			this.centriodOut = centroid(others);
			this.distanceOut = distance(others, others, true);
			this.distanceOutBelong = distance(others, al, false);
			this.distanceBelongOut = distance(al, others, false);
			Hashtable<String, Integer> tmp = computeNN(al, others);
			this.nnCase = tmp.containsKey(MEMBER) ? tmp.get(MEMBER) : 0;
			this.nnControl = tmp.containsKey(NONMEMBER) ? tmp.get(NONMEMBER) : 0;
		}

		private static Hashtable<String, Integer> computeNN(ArrayList<Integer> al, ArrayList<Integer> others) throws IllegalStateException {
			ArrayList<Integer> combined = new ArrayList<Integer>();
			ArrayList<String> member = new ArrayList<String>();
			Hashtable<String, Integer> totals = new Hashtable<String, Integer>();
			for (int i = 0; i < al.size(); i++) {
				combined.add(al.get(i));
				member.add(MEMBER);

			}
			for (int i = 0; i < others.size(); i++) {
				combined.add(others.get(i));
				member.add(NONMEMBER);
			}
			for (int i = 0; i < al.size(); i++) {
				ArrayList<Double> distances = new ArrayList<Double>();
				ArrayList<String> currentMemberDistance = new ArrayList<String>();

				for (int j = 0; j < combined.size(); j++) {
					distances.add((double) Math.abs(al.get(i) - combined.get(j)));
					currentMemberDistance.add(member.get(j));
				}

				int[] sort = Sort.trickSort(Array.toDoubleArray(distances), Array.toStringArray(currentMemberDistance));// so that member is favored when there are ties
				// for (int j = 0; j < sort.length; j++) {
				// System.out.println(j + "\t" + currentMemberDistance.get(sort[j]) + "\t" + distances.get(sort[j])+"\t"+al.size());
				// }
				// try {
				// Thread.sleep(1000);
				// } catch (InterruptedException ie) {
				// }
				int rank = 0;
				while (rank < Math.min(al.size(), others.size())) {
					String key = member.get(sort[rank]);
					if (rank == 0 && key != MEMBER) {
						throw new IllegalStateException("ALL v ALL did not return self when sorted");
					} else if (rank > 0) {
						if (totals.containsKey(key)) {
							totals.put(key, totals.get(key) + 1);
						} else {
							totals.put(key, 1);
						}
					}
					rank++;
				}
			}
			return totals;
		}

		private static String[] getHeader(String cases, String controls) {
			String[] header = Array.tagOn(BASE, cases, null);
			header = Array.concatAll(header, Array.tagOn(OTHER, cases, controls));
			// header = Array.concatAll(header, Array.tagOn(FINAL_METRIC, cases, null));
			header = Array.concatAll(header, Array.tagOn(BASE, controls, null));
			header = Array.concatAll(header, Array.tagOn(OTHER, controls, cases));

			header = Array.concatAll(header, Array.tagOn(NN, cases, cases));
			header = Array.concatAll(header, Array.tagOn(NN, cases, controls));

			// header = Array.concatAll(header, Array.tagOn(FINAL_METRIC, controls, null));
			return header;
		}

		private String[] getData() {
			ArrayList<String> data = new ArrayList<String>();
			data.add(centriodBelong + "");
			data.add(distanceBelong.getDistance() + "");
			data.add(distanceBelong.getIter() + "");
			data.add(distanceBelongOut.getDistance() + "");
			data.add(distanceBelongOut.getIter() + "");
			data.add(centriodOut + "");
			data.add(distanceOut.getDistance() + "");
			data.add(distanceOut.getIter() + "");

			data.add(distanceOutBelong.getDistance() + "");
			data.add(distanceOutBelong.getIter() + "");
			data.add(nnCase + "");
			data.add(nnControl + "");

			return Array.toStringArray(data);
		}
	}

	private static Hashtable<String, PosCluster[]> densityEnrichment(SimpleTallyResult cases, SimpleTallyResult controls, Logger log) throws IllegalStateException {
		Hashtable<String, PosCluster[]> cluster = new Hashtable<String, PosCluster[]>();
		System.out.println(cases.getFinalGeneVariantPositions());
		Hashtable<String, GeneVariantPositionSummary> casesSummaries = parseAvailable(GeneVariantPositionSummary.readSerial(cases.getFinalGeneVariantPositions(), log));
		System.out.println(controls.getFinalGeneVariantPositions());

		Hashtable<String, GeneVariantPositionSummary> controlSummaries = parseAvailable(GeneVariantPositionSummary.readSerial(controls.getFinalGeneVariantPositions(), log));
		for (String key : casesSummaries.keySet()) {
			GeneVariantPositionSummary currentCase = casesSummaries.get(key);
			PosCluster currentCaseClusReg = new PosCluster(currentCase.variantAAPositions);
			PosCluster currentCaseClusHQ = new PosCluster(currentCase.hqVariantAAPositions);

			GeneVariantPositionSummary currentControl = null;
			if (controlSummaries.containsKey(key)) {
				currentControl = controlSummaries.get(key);
				currentCaseClusReg.computeBelongs(currentControl.variantAAPositions);
				currentCaseClusHQ.computeBelongs(currentControl.hqVariantAAPositions);
			}
			cluster.put(key, new PosCluster[] { currentCaseClusReg, currentCaseClusHQ });
		}
		return cluster;
	}

	public static void test() {
		String vcf = "/home/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/aric_merge/vcf/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.aric.chargeMaf.vcf.gz";
		String popDir = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/aric_merge/vcf/Freq/";
		String[] vpopsCase = new String[] { popDir + "OSTEO_OFF_INHERIT.vpop", popDir + "ANIRIDIA.vpop", popDir + "ANOTIA.vpop" };
		// popDir + "CUSHING_FREQ.vpop", popDir + "EPP.vpop" };
		String omimDir = "/home/pankrat2/public/bin/ref/OMIM/";
		String[] otherGenesOfInterest = new String[] { popDir + "SB_T1.txt", popDir + "SB_T2.txt", "/home/pankrat2/public/bin/ref/COSMIC/cancer_gene_census.txt" };
		// ,popDir + "ALL_CONTROL_EPP.vpop", popDir + "ANIRIDIA.vpop", popDir + "ANOTIA.vpop" };
		int numThreads = 24;
		for (int i = 0; i < vpopsCase.length; i++) {

			double maf = 0.01;
			ArrayList<String> filesToWrite = new ArrayList<String>();
			ArrayList<String> names = new ArrayList<String>();
			if (vpopsCase[i].endsWith("OSTEO_OFF.vpop")) {
				maf = 0.001;
			}
			String outDir = ext.parseDirectoryOfFile(vpopsCase[i]);
			new File(outDir).mkdirs();
			Logger log = new Logger(ext.rootOf(vpopsCase[i], false) + ".log");
			if (i > 0) {
				log.reportTimeWarning("JOHN remember break");
				break;
			}
			GeneSet[] currentSets = GeneSet.load(Files.listFullPaths(ext.parseDirectoryOfFile(vpopsCase[i]), ext.rootOf(vpopsCase[i]) + ".geneset", false), ext.rootOf(vpopsCase[i]), log);
			log.reportTimeInfo("Found " + currentSets.length + " gene sets for " + vpopsCase[i]);
			OMIM omim = new OMIM(omimDir, log);
			OtherGeneInfo[] otherGeneInfos = null;
			if (otherGenesOfInterest != null) {
				otherGeneInfos = new OtherGeneInfo[otherGenesOfInterest.length];
				for (int j = 0; j < otherGeneInfos.length; j++) {
					otherGeneInfos[j] = new OtherGeneInfo(otherGenesOfInterest[j], log);
					otherGeneInfos[j].load();
				}
			}
			SimpleTallyResult caseResult = runSimpleTally(vcf, vpopsCase[i], maf, numThreads, outDir, currentSets, log);
			filesToWrite.add(caseResult.getFinalsampSummary());
			names.add("AnalysisInfo");
			filesToWrite.add(caseResult.getFinalAnnot());
			names.add("VARIANT_LEVEL");
			filesToWrite.add(caseResult.getFinalAnnotSample());
			names.add("VARIANT_SAMP_LEVEL");
			summarizeVariantsBySample(caseResult, log);
			VcfPopulation controls = caseResult.getControls();
			String controlFile = ext.parseDirectoryOfFile(vpopsCase[i]) + controls.getUniqSuperPop().get(0) + ".vpop";
			controls.report();
			controls.dump(controlFile);
			SimpleTallyResult controlResult = runSimpleTally(vcf, controlFile, maf, numThreads, outDir, currentSets, log);

			String geneFileCase = caseResult.getFinalAnnotGene();
			String geneFileControl = controlResult.getFinalAnnotGene();
			String caseWithControls = ext.addToRoot(geneFileCase, "_" + ext.rootOf(controlFile));
			filesToWrite.add(caseWithControls);
			names.add("GENE");
			Hashtable<String, PosCluster[]> cluster;
			try {
				cluster = densityEnrichment(caseResult, controlResult, log);
			} catch (IllegalStateException e) {
				log.reportTimeError("Could not enrich");
				e.printStackTrace();
				return;
			}
			Hashtable<String, String[]> controlsFuncHash = loadToGeneFuncHash(geneFileControl, log);
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(caseWithControls));

				BufferedReader reader = Files.getAppropriateReader(geneFileCase);
				String[] blanksEnrichment = new String[PosCluster.getHeader("BLANK", "BLANK").length * 2];
				Arrays.fill(blanksEnrichment, "0");

				String[] blanks = new String[GENE_ADD.length];
				Arrays.fill(blanks, "0");
				int line = 0;
				while (reader.ready()) {
					line++;
					String[] caseLine = reader.readLine().trim().split("\t");
					String key = caseLine[0] + "_" + caseLine[1];
					ArrayList<OMIMGene> oGene = omim.getOmimGene(caseLine[0]);
					if (controlsFuncHash.containsKey(key)) {
						writer.print(Array.toStr(caseLine) + "\t" + Array.toStr(Array.subArray(controlsFuncHash.get(key), 0, GENE_ADD.length)));// skip non-count data (like geneset denote)
					} else {
						writer.print(Array.toStr(caseLine) + "\t" + Array.toStr(blanks));
					}
					if (line == 1) {
						writer.print("\tOMIM_DISORDER(S)\tOMIM_GENE_STATUS\tOMIM_GENE_Symbol(s)\tOMIM_GENE_TITLE(s)");
						if (otherGeneInfos != null) {
							for (int j = 0; j < otherGeneInfos.length; j++) {
								writer.print("\t" + Array.toStr(otherGeneInfos[j].getHeader()));
							}
						}
						String[] baseHeader = PosCluster.getHeader(ext.rootOf(vpopsCase[i]), ext.rootOf(controlFile));
						writer.print("\t" + Array.toStr(baseHeader) + "\t" + Array.toStr(Array.tagOn(baseHeader, "HQ_", null)));

					} else {
						writer.print("\t");
						for (int j = 0; j < oGene.size(); j++) {
							writer.print((j == 0 ? "" : "|AdditionalOMIM") + (oGene.get(j).getDisorders().equals("") ? "NA" : oGene.get(j).getDisorders()) + (j == 0 ? "\t" : "|") + oGene.get(j).getStatus() + (j == 0 ? "\t" : "|") + Array.toStr(oGene.get(j).getGeneSymbols(), "/") + (j == 0 ? "\t" : "|") + oGene.get(j).getTitle());
						}
						if (otherGeneInfos != null) {
							for (int j = 0; j < otherGeneInfos.length; j++) {
								String gene = caseLine[0].toUpperCase();
								if (otherGeneInfos[j].getAnno().containsKey(gene)) {
									writer.print("\t" + Array.toStr(otherGeneInfos[j].getAnno().get(gene)));
								} else {
									String[] blank = Array.stringArray(otherGeneInfos[j].getHeader().length, "NA");
									writer.print("\t" + Array.toStr(blank));
								}
							}
						}
						if (cluster.containsKey(key)) {
							PosCluster[] tmp = cluster.get(key);
							writer.print("\t" + Array.toStr(tmp[0].getData()) + "\t" + Array.toStr(tmp[1].getData()));
						} else {
							writer.print("\t" + Array.toStr(blanksEnrichment));
						}
					}
					writer.println();
				}
				reader.close();
				writer.close();
				// String outXl = caseWithControls + ".xls";
				// ExcelWriter writerxl = new ExcelWriter(Array.toStringArray(filesToWrite), Array.toStringArray(names), log);
				// writerxl.write(outXl);

			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + geneFileCase + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + geneFileCase + "\"");
				return;
			}

		}

		// String[] vpopsControl = new String[] { popDir + "EPP.vpop", popDir + "ALL_CONTROL_EPP.vpop", popDir + "ALL_CONTROL_ANIRIDIA.vpop", popDir + "ALL_CONTROL_ANOTIA.vpop", popDir + "ANIRIDIA.vpop", popDir + "ANOTIA.vpop" };
		// for (int i = 0; i < vpopsControl.length; i++) {
		// String outDir = ext.parseDirectoryOfFile(vpopsControl[i]);
		// new File(outDir).mkdirs();
		// Logger log = new Logger(ext.rootOf(vpopsControl[i], false) + ".log");
		// runSimpleTally(vcf, vpopsControl[i], maf, numThreads, outDir, log);
		// }
	}

	private static void summarizeVariantsBySample(SimpleTallyResult sr, Logger log) {

		try {
			BufferedReader reader = Files.getAppropriateReader(sr.getFinalAnnotSample());
			int snpEFFIndex = ext.indexOfStr("SNPEFF_IMPACT", Files.getHeaderOfFile(sr.getFinalAnnotSample(), log));
			int sampIndex = ext.indexOfStr("SAMPLE", Files.getHeaderOfFile(sr.getFinalAnnotSample(), log));

			Hashtable<String, Integer> counts = new Hashtable<String, Integer>();
			reader.readLine();
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");

				String anyKey = line[sampIndex] + "\tANY_EFF";
				String funcKey = line[sampIndex] + "\t" + line[snpEFFIndex];
				if (!counts.containsKey(anyKey)) {
					counts.put(anyKey, 1);
				} else {
					counts.put(anyKey, counts.get(anyKey) + 1);
				}
				if (!counts.containsKey(funcKey)) {
					counts.put(funcKey, 1);
				} else {
					counts.put(funcKey, counts.get(funcKey) + 1);
				}
			}
			reader.close();
			String out = sr.getFinalAnnotSample() + ".varCounts";
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				writer.println("Sample\tSNPEFF_IMPACT\tCOUNTS");
				for (String key : counts.keySet()) {
					writer.println(key + "\t" + counts.get(key));
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + out);
				log.reportException(e);
			}
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + sr.getFinalAnnotSample() + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + sr.getFinalAnnotSample() + "\"");
			return;
		}

	}

	private static class SimpleTallyResult {
		private VcfPopulation controls;
		// private String finalOut;
		// private String finalOutVCF;
		private String finalAnnot;
		private String finalAnnotSample;
		private String finalsampSummary;
		private String finalAnnotGene;
		private String finalGeneVariantPositions;

		public SimpleTallyResult(VcfPopulation controls, String finalOut, String finalOutVCF, String finalsampSummary, String finalAnnot, String finalAnnotSample, String finalAnnotGene, String finalGeneVariantPositions) {
			super();
			this.controls = controls;
			// this.finalOut = finalOut;
			// this.finalOutVCF = finalOutVCF;
			this.finalsampSummary = finalsampSummary;
			this.finalAnnot = finalAnnot;
			this.finalAnnotGene = finalAnnotGene;
			this.finalAnnotSample = finalAnnotSample;
			this.finalGeneVariantPositions = finalGeneVariantPositions;
		}

		public String getFinalAnnot() {
			return finalAnnot;
		}

		public String getFinalGeneVariantPositions() {
			return finalGeneVariantPositions;
		}

		public VcfPopulation getControls() {
			return controls;
		}

		public String getFinalAnnotGene() {
			return finalAnnotGene;
		}

		public String getFinalAnnotSample() {
			return finalAnnotSample;
		}

		public String getFinalsampSummary() {
			return finalsampSummary;
		}

	}

	private static SimpleTallyResult runSimpleTally(String vcf, String vpop, double maf, int numThreads, String outDir, GeneSet[] geneSets, Logger log) {
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
		String finalAnnotSample = finalOut + ".summary.sample";
		String finalsampSummary = finalOut + ".sampSummary.txt";
		String finalAnnotGene = finalOut + ".gene";
		String finalAnnotGeneBed = finalOut + ".gene.bed";
		String finalGeneSetSummary = finalOut + ".geneset.summmary";
		String finalGeneVariantPositions = finalOut + ".gene.position.counts.ser";
		// String finalAnnotGeneSample = finalOut + ".gene.sample";

		Set<String> cases = vpopAc.getSuperPop().get(caseDef);
		vpopAc.report();
		Hashtable<String, Set<String>> controls = vpopAc.getSuperPop();
		controls.remove(caseDef);
		controls.remove(VcfPopulation.EXCLUDE);

		Set<String> tmpSet = new HashSet<String>();
		Hashtable<String, Set<String>> controlPop = new Hashtable<String, Set<String>>();
		controlPop.put(caseDef + "_" + VcfPopulation.CONTROL, tmpSet);

		for (String acontrolPop : controls.keySet()) {
			controlPop.get(caseDef + "_" + VcfPopulation.CONTROL).addAll(controls.get(acontrolPop));
		}
		VcfPopulation controlVcfPopulation = new VcfPopulation(controlPop, controlPop, POPULATION_TYPE.CASE_CONTROL, new Logger());
		SimpleTallyResult simpleTallyResult = new SimpleTallyResult(controlVcfPopulation, finalOut, finalOutVCF, finalsampSummary, finalAnnot, finalAnnotSample, finalAnnotGene, finalGeneVariantPositions);
		simpleTallyResult.getFinalGeneVariantPositions();

		summarizeAnalysisParams(finalsampSummary, caseDef, cases, controls, maf, log);
		summarizeGeneSets(geneSets, finalGeneSetSummary, log);
		if (!Files.exists(finalAnnotGene) || !Files.exists(finalAnnotGeneBed) || !Files.exists(finalAnnotSample) || !Files.exists(finalGeneVariantPositions)) {
			VCFFileReader tmp = new VCFFileReader(filtVcfs.get(0), true);

			VariantContextWriter writer = VCFOps.initWriter(finalOutVCF, VCFOps.DEFUALT_WRITER_OPTIONS, VCFOps.getSequenceDictionary(tmp));
			VCFOps.copyHeader(tmp, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
			tmp.close();

			PrintWriter annoGeneWriter = Files.getAppropriateWriter(finalAnnotGene);
			PrintWriter annoGeneBedWriter = Files.getAppropriateWriter(finalAnnotGeneBed);
			annoGeneBedWriter.println("CHR\tSTART\tSTOP\tGENENAME_FUNCTION");
			annoGeneWriter.print(Array.toStr(GENE_BASE));

			PrintWriter annoWriterSample = Files.getAppropriateWriter(finalAnnotSample);

			PrintWriter annoWriter = Files.getAppropriateWriter(finalAnnot);

			annoWriter.print(Array.toStr(ANNO_BASE));
			annoWriterSample.print(Array.toStr(ANNO_BASE_SAMPLE));
			annoWriter.print("\t" + Array.toStr(Array.tagOn(ANNO_ADD, caseDef + "_N_" + cases.size(), null)));
			annoWriterSample.print("\t" + Array.toStr(Array.tagOn(ANNO_ADD, caseDef + "_N_" + cases.size(), null)));

			annoGeneWriter.print("\t" + Array.toStr(Array.tagOn(GENE_ADD, caseDef + "_N_" + cases.size(), null)));

			ArrayList<String> controlsOrdered = new ArrayList<String>();
			controlsOrdered.addAll(controls.keySet());
			for (int i = 0; i < controlsOrdered.size(); i++) {
				String annotLine = "\t" + Array.toStr(Array.tagOn(ANNO_ADD, controlsOrdered.get(i) + "_N_" + controls.get(controlsOrdered.get(i)).size(), null));
				annoWriter.print(annotLine);
				annoWriterSample.print(annotLine);
				annoGeneWriter.print("\t" + Array.toStr(Array.tagOn(GENE_ADD, controlsOrdered.get(i) + "_N_" + controls.get(controlsOrdered.get(i)).size(), null)));
			}
			annoGeneWriter.print("\tIS_GENE_SET");

			String[][] annotations = VCFOps.getAnnotationKeys(vcf, log);
			annoWriter.print("\t" + Array.toStr(annotations[0]));
			annoWriterSample.print("\t" + Array.toStr(annotations[0]));
			for (int i = 0; i < geneSets.length; i++) {
				annoGeneWriter.print("\t" + geneSets[i].getTag() + "_Membership");
				annoWriter.print("\t" + geneSets[i].getTag() + "_Membership");
				annoWriterSample.print("\t" + geneSets[i].getTag() + "_Membership");
			}
			annoGeneWriter.println();
			annoWriter.println();
			annoWriterSample.println();

			VariantContextFilter qual = getQualityFilterwkggseq(log);
			Hashtable<String, ArrayList<GeneSummary[]>> geneSummaries = new Hashtable<String, ArrayList<GeneSummary[]>>();
			for (int i = 0; i < filtVcfs.size(); i++) {
				// if (i > 1) {
				// log.reportTimeWarning("JOHN remember break");
				// break;
				// }
				log.reportTimeInfo("Summarizing " + filtVcfs.get(i));
				VCFFileReader result = new VCFFileReader(filtVcfs.get(i), true);
				for (VariantContext vc : result) {
					writer.add(vc);
					Segment seg = VCOps.getSegment(vc);

					String geneName = VCOps.getSNP_EFFGeneName(vc);
					String func = VCOps.getSNP_EFFImpact(vc);
					annoGeneBedWriter.println(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop() + "\t" + geneName + ":" + func);
					if (!geneSummaries.containsKey(geneName)) {
						addEntries(caseDef, controlsOrdered, geneSummaries, geneName);
					}
					for (int j = 0; j < geneSets.length; j++) {
						if (geneSets[j].getGenes().containsKey(geneName) && !geneSummaries.containsKey(geneSets[j].getTag())) {
							addEntries(caseDef, controlsOrdered, geneSummaries, geneSets[j].getTag());
						}
					}
					VcGroupSummary vcCaseGroup = new VcGroupSummary(caseDef, cases, vc, qual, log);
					// TODO, check
					for (int j = 0; j < geneSummaries.get(geneName).get(0).length; j++) {
						geneSummaries.get(geneName).get(0)[j].add(vcCaseGroup, null);
						for (int j2 = 0; j2 < geneSets.length; j2++) {
							if (geneSets[j2].getGenes().containsKey(geneName)) {

								geneSummaries.get(geneSets[j2].getTag()).get(0)[j].add(vcCaseGroup, geneSets[j2].getTag());
							}
						}
					}

					annoWriter.print(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getID() + "\t" + vc.getReference().getBaseString() + "\t" + vc.getAlternateAlleles().toString());
					annoWriter.print("\t" + Array.toStr(vcCaseGroup.getSummary()));

					GenotypesContext gc = vcCaseGroup.getVcAlt().getGenotypes();

					ArrayList<VcGroupSummary> controlGroupSummaries = new ArrayList<VCFSimpleTally.VcGroupSummary>();
					for (int j = 0; j < controlsOrdered.size(); j++) {
						VcGroupSummary vcControlGroup = new VcGroupSummary(controlsOrdered.get(j), controls.get(controlsOrdered.get(j)), vc, qual, log);
						controlGroupSummaries.add(vcControlGroup);
					}

					for (Genotype g : gc) {
						annoWriterSample.print(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getID() + "\t" + vc.getReference().getBaseString() + "\t" + vc.getAlternateAlleles().toString() + "\t" + g.getSampleName() + "\t" + g.toString());
						annoWriterSample.print("\t" + Array.toStr(vcCaseGroup.getSummary()));
						for (int j = 0; j < controlsOrdered.size(); j++) {
							annoWriterSample.print("\t" + Array.toStr(controlGroupSummaries.get(j).getSummary()));
						}
						annoWriterSample.print("\t" + Array.toStr(VCOps.getAnnotationsFor(annotations[0], vc, ".")));
						for (int j = 0; j < geneSets.length; j++) {
							annoWriterSample.print("\t" + geneSets[j].getGenes().containsKey(geneName));
						}
						annoWriterSample.println();
					}

					for (int j = 0; j < controlsOrdered.size(); j++) {
						VcGroupSummary vcControlGroup = controlGroupSummaries.get(j);
						for (int k = 0; k < geneSummaries.get(geneName).get(0).length; k++) {
							geneSummaries.get(geneName).get(j + 1)[k].add(vcControlGroup, null);
							for (int j2 = 0; j2 < geneSets.length; j2++) {
								if (geneSets[j2].getGenes().containsKey(geneName)) {
									geneSummaries.get(geneSets[j2].getTag()).get(j + 1)[k].add(vcControlGroup, geneSets[j2].getTag());
								}
							}
						}
						annoWriter.print("\t" + Array.toStr(vcControlGroup.getSummary()));
					}
					annoWriter.print("\t" + Array.toStr(VCOps.getAnnotationsFor(annotations[0], vc, ".")));
					for (int j = 0; j < geneSets.length; j++) {
						annoWriter.print("\t" + geneSets[j].getGenes().containsKey(geneName));
					}
					annoWriter.println();
				}
				result.close();
			}
			annoWriterSample.close();
			annoGeneBedWriter.close();
			annoWriter.close();
			writer.close();
			ArrayList<GeneVariantPositionSummary> controlPos = new ArrayList<GeneVariantPositionSummary>();
			for (String gene : geneSummaries.keySet()) {
				ArrayList<GeneSummary[]> geneSummariesCurrent = geneSummaries.get(gene);

				for (int i = 0; i < geneSummariesCurrent.get(0).length; i++) {
					controlPos.add(geneSummariesCurrent.get(0)[i].getGeneVariantPositionSummary());// case only
					annoGeneWriter.print(gene + "\t" + Array.toStr(geneSummariesCurrent.get(0)[i].getEffects(), "||"));
					for (int j = 0; j < geneSummariesCurrent.size(); j++) {
						annoGeneWriter.print("\t" + Array.toStr(geneSummariesCurrent.get(j)[i].getSummary()));
					}

					annoGeneWriter.print("\t" + isGeneSet(geneSets, gene));
					for (int j = 0; j < geneSets.length; j++) {
						annoGeneWriter.print("\t" + geneSets[j].getGenes().containsKey(gene));
					}
					annoGeneWriter.println();
				}
			}
			annoGeneWriter.close();
			GeneVariantPositionSummary.writeSerial(controlPos.toArray(new GeneVariantPositionSummary[controlPos.size()]), finalGeneVariantPositions, log);
			VCFOps.VcfPopulation.splitVcfByPopulation(finalOutVCF, vpop, true, true, log);
		} else {
			log.reportTimeWarning(finalAnnotGene + " exists so skipping summarize");
		}
		return simpleTallyResult;
	}

	private static boolean isGeneSet(GeneSet[] geneSets, String tag) {
		for (int i = 0; i < geneSets.length; i++) {
			if (geneSets[i].getTag().equals(tag)) {
				return true;
			}
		}
		return false;
	}

	private static void summarizeGeneSets(GeneSet[] geneSets, String sumFile, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(sumFile));
			writer.println("GENE_SET\tN_GENES\tGENES->");
			for (int i = 0; i < geneSets.length; i++) {
				writer.print(geneSets[i].getTag() + "\t" + geneSets[i].getGenes().size());
				for (String gene : geneSets[i].getGenes().keySet()) {
					writer.print("\t" + gene);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + sumFile);
			log.reportException(e);
		}
	}

	private static void summarizeAnalysisParams(String sumFile, String caseDef, Set<String> cases, Hashtable<String, Set<String>> controls, double maf, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(sumFile));
			writer.print("#CASE\t" + caseDef + "\tn=" + cases.size());
			for (String acase : cases) {
				writer.print("\t" + acase);
			}
			writer.println();
			for (String control : controls.keySet()) {
				writer.print("#CONTROL\t" + control + "\tn=" + controls.get(control).size());
				for (String acontrol : controls.get(control)) {
					writer.print("\t" + acontrol);
				}
				writer.println();
			}
			// writer.println("#FILTERS:");
			// writer.println(ESP_FILTER);
			// writer.println(G1000_FILTER);
			// writer.println(CHARGE_B_FILTER);
			// writer.println(CHARGE_W_FILTER);
			// for (String control : controls.keySet()) {
			// writer.print(control + ": maf < " + maf);
			// writer.println();
			// }
			// writer.println("#HQ_GENO_FILTERS:");

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + sumFile);
			log.reportException(e);
		}
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
			geneSummaries.get(geneName).add(controlGeneSummaries);
		}
	}

	public static class GeneVariantPositionSummary implements Serializable {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		private String name;
		private String group;
		private String key;
		private ArrayList<Integer> variantNPositions;
		private ArrayList<Integer> hqVariantNPositions;
		private ArrayList<Integer> variantAAPositions;
		private ArrayList<Integer> hqVariantAAPositions;
		private String[] effects;

		public GeneVariantPositionSummary(String name, String group, String[] effects) {
			super();
			this.name = name;
			this.group = group;
			this.effects = effects;
			this.variantNPositions = new ArrayList<Integer>();
			this.variantAAPositions = new ArrayList<Integer>();
			this.hqVariantNPositions = new ArrayList<Integer>();
			this.hqVariantAAPositions = new ArrayList<Integer>();
			this.key = name + "_" + Array.toStr(effects, "||");
		}

		public String getKey() {
			return key;
		}

		public String[] getEffects() {
			return effects;
		}

		public static void writeSerial(GeneVariantPositionSummary[] geneVariantSummaries, String filename, Logger log) {
			Files.writeSerial(geneVariantSummaries, filename, true);
		}

		public static GeneVariantPositionSummary[] readSerial(String filename, Logger log) {
			return (GeneVariantPositionSummary[]) Files.readSerial(filename, false, log, false, true);
		}

		public static long getSerialversionuid() {
			return serialVersionUID;
		}

		public String getName() {
			return name;
		}

		public String getGroup() {
			return group;
		}

		public enum ADD_TYPE {
			REGULAR, HQ;
		}

		public void add(VariantContext vc, ADD_TYPE type) {
			String aaChange = VCOps.getAnnotationsFor(new String[] { "SNPEFF_AMINO_ACID_CHANGE" }, vc, ".")[0];
			int aapos = -1;
			if (!aaChange.equals(".")) {
				aapos = Integer.parseInt(aaChange.replaceAll("[^\\d.]", ""));
			}
			int nucPos = vc.getStart();

			switch (type) {
			case HQ:
				for (int i = 0; i < vc.getSampleNames().size(); i++) {
					hqVariantAAPositions.add(aapos);
					hqVariantNPositions.add(nucPos);
				}

				break;
			case REGULAR:
				for (int i = 0; i < vc.getSampleNames().size(); i++) {
					variantAAPositions.add(aapos);
					variantNPositions.add(nucPos);
				}
				break;
			default:
				break;
			}
		}
	}

	private static class GeneSummary {
		private GeneVariantPositionSummary geneVariantPositionSummary;
		private String geneName;
		private String group;
		private String[] effects;
		private int numVar;
		private HashSet<String> uniqueIndsWithVar;
		private int hqNumVar;
		private HashSet<String> uniqueHqIndsWithVar;
		private Hashtable<String, Integer> numVarsPerInd;
		private Hashtable<String, Integer> numVarsDiffHaploPerInd;

		private Hashtable<String, Integer> numHQVarsPerInd;
		private Hashtable<String, Integer> numHQDiffHaploVarsPerInd;

		private Hashtable<String, ArrayList<String>> sampleHQHaplotype;
		private Hashtable<String, ArrayList<String>> sampleHaplotype;

		private GeneSummary(String geneName, String group, String[] effects) {
			super();
			this.group = group;
			this.geneName = geneName;
			this.numVar = 0;
			this.hqNumVar = 0;
			this.effects = effects;
			this.uniqueHqIndsWithVar = new HashSet<String>();
			this.uniqueIndsWithVar = new HashSet<String>();
			this.numVarsPerInd = new Hashtable<String, Integer>();
			this.numHQVarsPerInd = new Hashtable<String, Integer>();
			this.numHQDiffHaploVarsPerInd = new Hashtable<String, Integer>();
			this.numVarsDiffHaploPerInd = new Hashtable<String, Integer>();
			this.sampleHaplotype = new Hashtable<String, ArrayList<String>>();
			this.sampleHQHaplotype = new Hashtable<String, ArrayList<String>>();
			this.geneVariantPositionSummary = new GeneVariantPositionSummary(geneName, group, effects);
		}

		public GeneVariantPositionSummary getGeneVariantPositionSummary() {
			return geneVariantPositionSummary;
		}

		public String[] getEffects() {
			return effects;
		}

		private static int numGreaterThan(Hashtable<String, Integer> numHash, int numReq) {
			int num = 0;
			for (String key : numHash.keySet()) {
				if (numHash.get(key) >= numReq) {
					num++;
				}
			}
			return num;
		}

		private String[] getSummary() {
			ArrayList<String> summary = new ArrayList<String>();
			summary.add(numVar + "");
			summary.add(uniqueIndsWithVar.size() + "");
			summary.add(numGreaterThan(numVarsPerInd, 2) + "");
			summary.add(numGreaterThan(numVarsDiffHaploPerInd, 2) + "");

			summary.add(hqNumVar + "");
			summary.add(uniqueHqIndsWithVar.size() + "");
			summary.add(numGreaterThan(numHQVarsPerInd, 2) + "");
			summary.add(numGreaterThan(numHQDiffHaploVarsPerInd, 2) + "");

			return Array.toStringArray(summary);
		}

		private static void addHash(Hashtable<String, Integer> toAdd, Set<String> from) {
			for (String key : from) {
				if (!toAdd.containsKey(key)) {
					toAdd.put(key, 1);
				} else {
					int num = toAdd.get(key) + 1;
					toAdd.put(key, num);
				}
			}
		}

		private void add(VcGroupSummary vcGroupSummary, String setTag) {
			if (!vcGroupSummary.getGroupName().equals(group)) {

				throw new IllegalArgumentException("Mismatched group names: should be " + group + " but actually " + vcGroupSummary.getGroupName());

			}
			if (!VCOps.getSNP_EFFGeneName(vcGroupSummary.getVcOriginal()).equals(geneName) && (setTag != null && !setTag.equals(geneName))) {
				throw new IllegalArgumentException("Mismatched gene/geneSet names");
			} else {
				String impact = VCOps.getSNP_EFFImpact(vcGroupSummary.getVcOriginal());
				if (ext.indexOfStr(impact, effects) >= 0) {
					if (vcGroupSummary.getIndsWithAlt().size() > 0) {
						numVar++;
						uniqueIndsWithVar.addAll(vcGroupSummary.getIndsWithAlt());
						addHash(numVarsPerInd, vcGroupSummary.getIndsWithAlt());

						GenotypesContext gc = vcGroupSummary.getVcAlt().getGenotypes();
						HashSet<String> newHaplos = new HashSet<String>();
						for (Genotype g : gc) {
							if (!sampleHaplotype.containsKey(g.getSampleName())) {
								sampleHaplotype.put(g.getSampleName(), new ArrayList<String>());
							}
							String pid = (String) g.getAnyAttribute("PID");
							if (pid == null || !sampleHaplotype.get(g.getSampleName()).contains(pid)) {
								newHaplos.add(g.getSampleName());
								sampleHaplotype.get(g.getSampleName()).add(pid);
							}
						}
						addHash(numVarsDiffHaploPerInd, newHaplos);
						geneVariantPositionSummary.add(vcGroupSummary.getVcAlt(), ADD_TYPE.REGULAR);
						if (vcGroupSummary.getHqIndsWithAlt().size() > 0) {
							hqNumVar++;
							uniqueHqIndsWithVar.addAll(vcGroupSummary.getHqIndsWithAlt());
							addHash(numHQVarsPerInd, vcGroupSummary.getHqIndsWithAlt());

							GenotypesContext gcHQ = vcGroupSummary.getVcAltHq().getGenotypes();
							HashSet<String> newHaplosHQ = new HashSet<String>();
							for (Genotype g : gcHQ) {
								if (!sampleHQHaplotype.containsKey(g.getSampleName())) {
									sampleHQHaplotype.put(g.getSampleName(), new ArrayList<String>());
								}
								String pid = (String) g.getAnyAttribute("PID");
								if (pid == null || !sampleHQHaplotype.get(g.getSampleName()).contains(pid)) {
									newHaplosHQ.add(g.getSampleName());
									sampleHQHaplotype.get(g.getSampleName()).add(pid);
								}
							}
							geneVariantPositionSummary.add(vcGroupSummary.getVcAltHq(), ADD_TYPE.HQ);
							addHash(numHQDiffHaploVarsPerInd, newHaplosHQ);
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
		private double avgGq, avgDp;
		private int numIndsAlt;
		private int numWithCalls;
		private int aac;
		private int numHqIndsAlt;
		private int hqAac;
		private Set<String> indsWithAlt;
		private Set<String> hqIndsWithAlt;
		private VariantContextFilter filter;
		private VariantContext vcSub;
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

		// public VariantContext getVcSub() {
		// return vcSub;
		// }

		public VariantContext getVcAlt() {
			return vcAlt;
		}

		public VariantContext getVcAltHq() {
			return vcAltHq;
		}

		private void summarize() {
			this.vcSub = VCOps.getSubset(vcOriginal, group, VC_SUBSET_TYPE.SUBSET_STRICT, false);
			this.avgGq = VCOps.getAvgGenotypeInfo(vcSub, null, GENOTYPE_INFO.GQ, log);
			this.avgDp = VCOps.getAvgGenotypeInfo(vcSub, null, GENOTYPE_INFO.DP, log);
			this.vcAlt = VCOps.getAltAlleleContext(vcSub, null, null, ALT_ALLELE_CONTEXT_TYPE.ALL, false, log);
			this.aac = (int) VCOps.getAAC(vcSub, null);
			this.indsWithAlt = vcAlt.getSampleNames();
			this.numIndsAlt = indsWithAlt.size();
			this.numWithCalls = vcSub.getSampleNames().size() - vcSub.getNoCallCount();
			this.vcAltHq = VCOps.getAltAlleleContext(vcSub, null, filter, ALT_ALLELE_CONTEXT_TYPE.ALL, false, log);
			this.hqAac = (int) VCOps.getAAC(vcAltHq, null);
			this.hqIndsWithAlt = vcAltHq.getSampleNames();
			this.numHqIndsAlt = hqIndsWithAlt.size();
		}

		private String[] getSummary() {
			ArrayList<String> summary = new ArrayList<String>();
			summary.add(avgGq + "");
			summary.add(avgDp + "");
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

	/**
	 * @param log
	 * @return the {@link htsjdk.variant.variantcontext.filter.VariantContextFilter} used for genotype qc
	 */
	private static VariantContextFilter getQualityFilterwkggseq(Logger log) {
		VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_LOOSE;
		gq.setDFilter(50);
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		dp.setDFilter(10);
		VARIANT_FILTER_DOUBLE altD = VARIANT_FILTER_DOUBLE.ALT_ALLELE_DEPTH;
		altD.setDFilter(3);

		// VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		// VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
		// VARIANT_FILTER_BOOLEAN amb = VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER;
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		// VARIANT_FILTER_BOOLEAN[] bQualFilts = new VARIANT_FILTER_BOOLEAN[] { amb };
		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, dp, altD, gq };
		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, new String[] { "G1000Freq" }, new String[] { ESP_FILTER + AND + G1000_FILTER }, log);
		return vContextFilter;
	}

	public static void main(String[] args) {
		test();
	}

	// boolean setMember = false;
	// for (int j = 0; j < geneSets.length; j++) {
	// if (geneSets[j].getGenes().containsKey(gene)) {
	// if (!setMember) {
	// annoGeneWriter.print("\t" + geneSets[j].getTag());
	// } else {
	// annoGeneWriter.print(";" + geneSets[j].getTag());
	// }
	// setMember = true;
	//
	// }
	// }
	// if (!setMember) {
	// annoGeneWriter.print("\t.");
	// }
}
