package org.genvisis.seq.analysis;

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
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import org.genvisis.bioinformatics.OMIM;
import org.genvisis.bioinformatics.OMIM.OMIMGene;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ExcelConverter;
import org.genvisis.common.ExcelConverter.ExcelConversionParams;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.VCFSimpleTally.GeneVariantPositionSummary.ADD_TYPE;
import org.genvisis.seq.manage.GenotypeOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.ChrSplitResults;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import org.genvisis.seq.manage.VCOps.GENOTYPE_INFO;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilterPass;
import org.genvisis.seq.qc.FilterNGS.VcFilterBoolean;
import org.genvisis.seq.qc.FilterNGS.VcFilterDouble;

import com.google.common.primitives.Ints;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 *
 *
 */
public class VCFSimpleTally {
	public static final String[] EFF = {"HIGH", "MODERATE", "LOW"};
	public static final String ANY_EFF = "OTHER";
	private static final String[][] EFF_DEFS = new String[][] {new String[] {EFF[0]},
																														 new String[] {EFF[0], EFF[1]},
																														 new String[] {EFF[0], EFF[1], EFF[2]},
																														 new String[] {EFF[0], EFF[1], EFF[2],
																																					 ANY_EFF}};

	private static final String ANY_GENE_SET = "*";
	// private static final String SNPEFF_IMPACTS =
	// "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE'||SNPEFF_IMPACT=='LOW')";
	private static final String SNPEFF_NAMES = "G1000_esp_charge_aricFreq_SNPEFF_HIGH_MODERATE_LOW";
	// private static final String CHARGE_B_FILTER = "(charge.MAF_blacks=='.'||charge.MAF_blacks <=
	// 0.01)";
	// private static final String CHARGE_W_FILTER = "(charge.MAF_whites=='.'||charge.MAF_whites <=
	// 0.01)";
	private static final String[] ANNO_BASE = new String[] {"CHROM", "POS", "ID", "REF", "ALT",
																													"BIALLELIC", "FILTERS",
																													"HIGH||MODERATE||LOW"};
	private static final String[] ANNO_BASE_SAMPLE = ArrayUtils.concatAll(ANNO_BASE,
																																				new String[] {"SAMPLE",
																																											"GENOTYPE",
																																											"HQ"});

	private static final String[] ANNO_ADD = new String[] {"_AVG_GQ", "_AVG_DP", "_NUM_WITH_CALLS",
																												 "_NUM_WITH_ALT", "_AAC",
																												 "_HQ_NUM_WITH_ALT", "_HQ_AAC",};
	private static final String[] GENE_BASE = new String[] {"GENE/GENE_SET", "FUNCTIONAL_TYPE"};
	private static final String[] GENE_ADD = new String[] {"numVar", "numMut", "uniqInds",
																												 "numCompoundHets",
																												 "numCompoundHetsDiffHaplotype", "hqNumVar",
																												 "hqNumMut", "hqUniqInds",
																												 "numHQCompoundHets",
																												 "numHQCompoundHetsDiffHaplotype"};

	private static boolean filterCHARGE(VariantContext vc, double maf) {
		boolean pass = true;
		if (vc.hasAttribute("charge.MAF_whites")) {
			pass = vc.getCommonInfo().getAttributeAsDouble("charge.MAF_whites", 0) <= maf;
		}
		if (pass && vc.hasAttribute("charge.MAF_blacks")) {
			pass = vc.getCommonInfo().getAttributeAsDouble("charge.MAF_blacks", 0) <= maf;
		}
		return pass;
	}

	private static void filter(String vcf, String output, String casePop, VcfPopulation vpop,
														 double controlFreq, Logger log) {
		if (!Files.exists(output)) {
			Set<String> cases = vpop.getSuperPop().get(casePop);
			Hashtable<String, Set<String>> controls = vpop.getSuperPop();

			log.reportTimeInfo("CASE :" + casePop + " n: " + cases.size());
			controls.remove(casePop);
			controls.remove(VcfPopulation.EXCLUDE);
			for (String control : controls.keySet()) {
				log.reportTimeInfo("Control: " + control + " n: " + controls.get(control).size());
			}
			VCFFileReader reader = new VCFFileReader(new File(vcf), true);
			VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS,
																											reader.getFileHeader()
																														.getSequenceDictionary());
			VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
			VariantContextFilter freqFilter = getFreqFilter(controlFreq, log);
			int numScanned = 0;
			int numPass = 0;
			log.reportTimeWarning("Note, bialliec variant frequencies will be inaccurate");
			for (VariantContext vc : reader) {
				numScanned++;
				if (numScanned % 10000 == 0) {
					log.reportTimeInfo(numScanned + " variants scanned, " + numPass + " variants passed");
				}

				// System.out.println(Array.toStr(VCOps.getAnnotationsFor(new String[] {
				// "charge.MAF_whites", "charge.MAF_blacks" }, vc, ".")));

				// if (!vc.isFiltered()) {// no tranche
				VariantContext vcCase = VCOps.getSubset(vc, cases, VC_SUBSET_TYPE.SUBSET_STRICT, false);
				if (vcCase.getSampleNames().size() != cases.size()) {
					String[] allSamps = VCFOps.getSamplesInFile(vcf);
					for (String aCase : cases) {
						log.reportTimeWarning("CASE: " + aCase + " in vcf "
																	+ (ext.indexOfStr(aCase, allSamps) >= 0));
					}
					Files.writeArray(allSamps, ext.rootOf(casePop, false) + "samplesToPickFrom.txt");
					throw new IllegalArgumentException("could not find all cases for " + casePop);
				}
				if (!vcCase.isMonomorphicInSamples() && vcCase.getNoCallCount() != cases.size()
						&& (!vc.hasAttribute("esp6500si_all") || !vc.hasAttribute("g10002014oct_all")
								|| !vc.hasAttribute("g10002015aug_all") || !vc.hasAttribute("esp6500siv2_all"))) {
					// String error = "Expected annotations esp6500si_all, g10002014oct_all were not present";
					// error += "\n" + vc.toStringWithoutGenotypes();
					// throw new IllegalStateException(error);
				} else if (vcCase.getHomRefCount() + vcCase.getNoCallCount() != cases.size()
									 && vcCase.getNoCallCount() != cases.size() && freqFilter.filter(vcCase).passed()
									 && filterCHARGE(vcCase, controlFreq)) {// as alts in rare esp/1000g
					boolean controlPass = true;
					for (String controlPop : controls.keySet()) {
						VariantContext vcControl = VCOps.getSubset(vc, controls.get(controlPop),
																											 VC_SUBSET_TYPE.SUBSET_STRICT, false);
						if (vcControl.getSampleNames().size() != controls.get(controlPop).size()) {
							throw new IllegalArgumentException("could not find all controls for " + controlPop);
						}
						if (controlFreq < 1) {
							double maf = VCOps.getMAF(vcControl, null);
							if ((!VCOps.isMinorAlleleAlternate(vcControl, null) || maf > controlFreq)
									&& vcControl.getNoCallCount() != controls.get(controlPop).size()) {// rare in
																																										 // control
								controlPass = false;
								break;
							}
						}
					}
					if (controlPass) {
						numPass++;
						writer.add(vc);
					}
				}
			}
			// }
			reader.close();
			writer.close();
		} else {
			log.reportFileExists(output);
		}
	}

	private static VariantContextFilter getFreqFilter(double maf, Logger log) {
		VariantContextFilter vContextFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {},
																																	 new VARIANT_FILTER_BOOLEAN[] {},
																																	 new String[] {"RARE_"
																																								 + SNPEFF_NAMES},
																																	 new String[] {FilterNGS.getPopFreqFilterString(maf)},
																																	 log);
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
				String[] rest = ArrayUtils.subArray(line, 2);
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
		private final String fileName;
		private final Logger log;
		private Hashtable<String, String> genes;
		private final String tag;

		public GeneSet(String fileName, Logger log, Hashtable<String, String> genes, String tag) {
			super();
			this.fileName = fileName;
			this.log = log;
			this.genes = genes;
			this.tag = tag;
		}

		private static GeneSet getAnySet(Logger log) {
			Hashtable<String, String> permiscuousSet = new Hashtable<String, String>();
			permiscuousSet.put(ANY_GENE_SET, ANY_GENE_SET);
			return new GeneSet(null, log, permiscuousSet, ANY_GENE_SET);
		}

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
			sets = ArrayUtils.concatAll(sets, new GeneSet[] {getAnySet(log)});
			return sets;
		}

		public Hashtable<String, String> getGenes() {
			return genes;
		}

		private void load() {
			genes = new Hashtable<String, String>();
			try {
				BufferedReader reader = Files.getAppropriateReader(fileName);

				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					String gene = line[0].trim();
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
		private final String fileName;
		private final Logger log;
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
			anno = new Hashtable<String, String[]>();
			try {
				BufferedReader reader = Files.getAppropriateReader(fileName);
				String[] head = reader.readLine().trim().split("\t");
				if (!head[0].equals(GENE_TAG)) {
					throw new IllegalArgumentException("File " + fileName + " must have " + GENE_TAG
																						 + " to be used for extra gene annotation, found "
																						 + head[0] + " instead");
				}

				header = ArrayUtils.tagOn(head, ext.removeDirectoryInfo(fileName), null);
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					String gene = line[0].toUpperCase();
					if (anno.containsKey(gene)) {
						throw new IllegalAccessError("One entry per gene per file, " + gene
																				 + " was seen twice");

					} else {
						if (line.length > header.length) {
							throw new IllegalArgumentException("Error, lines cannot be greater length than the header in"
																								 + fileName);
						}
						ArrayList<String> filled = new ArrayList<String>();
						for (int i = 0; i < header.length; i++) {
							if (i < line.length) {
								filled.add(line[i]);
							} else {
								filled.add("NA");
							}
						}
						if (filled.size() != header.length) {
							throw new IllegalArgumentException("Error, could not fill line in" + fileName);

						}
						anno.put(gene, ArrayUtils.toStringArray(filled));
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
		for (GeneVariantPositionSummary summarie : summaries) {
			String key = summarie.getKey();

			if (parse.containsKey(key)) {
				throw new IllegalArgumentException("Dup keys");
			} else {
				parse.put(key, summarie);
			}
		}
		return parse;
	}

	private static ArrayList<Integer> removeNeg(List<Integer> others) {
		ArrayList<Integer> pos = new ArrayList<Integer>();
		for (int i = 0; i < others.size(); i++) {
			if (others.get(i) >= 0) {
				pos.add(others.get(i));
			}
		}
		return pos;
	}

	private static double centroid(List<Integer> al) {
		return ArrayUtils.mean(Ints.toArray(al));
	}

	// private static double distanceSED(double cent, ArrayList<Integer> al) {
	// double sum = 0;
	// for (int i = 0; i < al.size(); i++) {
	// sum += Math.pow((double) al.get(i) - cent, 2);
	// }
	// return Math.sqrt(sum);
	// }
	private static class DistanceResult {
		private final double distance;
		private final int iter;

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

	private static DistanceResult distance(List<Integer> from, List<Integer> to,
																				 boolean identical) {
		double sum = 0;
		int num = 0;
		for (int i = 0; i < from.size(); i++) {
			for (int j = identical ? i + 1 : 0; j < to.size(); j++) {
				sum += Math.abs(from.get(i) - to.get(j));
				num++;
			}
		}
		double distance = sum / num;
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
		private static final String[] BASE = new String[] {"_AA_CENTROID", "_AVG_DISTANCE", "_N_COMP"};
		private static final String[] OTHER = new String[] {"_DIS_FROM_", "_N_COMP_"};
		private static final String[] NN = new String[] {"_NumNearestTo_"};

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

		private final ArrayList<Integer> al;

		public PosCluster(List<Integer> al) {
			super();
			this.al = removeNeg(al);
			distanceOut = new DistanceResult(Double.NaN, 0);
			distanceOutBelong = new DistanceResult(Double.NaN, 0);
			distanceBelongOut = new DistanceResult(Double.NaN, 0);
			init();
		}

		private void init() {
			if (al.size() > 0) {
				centriodBelong = centroid(al);
				distanceBelong = distance(al, al, true);
			} else {
				centriodBelong = Double.NaN;
				distanceBelong = new DistanceResult(Double.NaN, 0);
			}
		}

		private void computeBelongs(List<Integer> others) throws IllegalStateException {
			others = removeNeg(others);
			centriodOut = centroid(others);
			distanceOut = distance(others, others, true);
			distanceOutBelong = distance(others, al, false);
			distanceBelongOut = distance(al, others, false);
			Hashtable<String, Integer> tmp = computeNN(al, others);
			nnCase = tmp.containsKey(MEMBER) ? tmp.get(MEMBER) : 0;
			nnControl = tmp.containsKey(NONMEMBER) ? tmp.get(NONMEMBER) : 0;
		}

		private static Hashtable<String, Integer> computeNN(List<Integer> al,
																												List<Integer> others) throws IllegalStateException {
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

				int[] sort = Sort.getSort2DIndices(distances, currentMemberDistance);
				// for (int j = 0; j < sort.length; j++) {
				// System.out.println(j + "\t" + currentMemberDistance.get(sort[j]) + "\t" +
				// distances.get(sort[j])+"\t"+al.size());
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
			String[] header = ArrayUtils.tagOn(BASE, cases, null);
			header = ArrayUtils.concatAll(header, ArrayUtils.tagOn(OTHER, cases, controls));
			// header = Array.concatAll(header, Array.tagOn(FINAL_METRIC, cases, null));
			header = ArrayUtils.concatAll(header, ArrayUtils.tagOn(BASE, controls, null));
			header = ArrayUtils.concatAll(header, ArrayUtils.tagOn(OTHER, controls, cases));

			header = ArrayUtils.concatAll(header, ArrayUtils.tagOn(NN, cases, cases));
			header = ArrayUtils.concatAll(header, ArrayUtils.tagOn(NN, cases, controls));

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

			return ArrayUtils.toStringArray(data);
		}
	}

	private static Hashtable<String, PosCluster[]> densityEnrichment(SimpleTallyResult cases,
																																	 SimpleTallyResult controls,
																																	 Logger log) throws IllegalStateException {
		Hashtable<String, PosCluster[]> cluster = new Hashtable<String, PosCluster[]>();
		Hashtable<String, GeneVariantPositionSummary> casesSummaries = parseAvailable(GeneVariantPositionSummary.readSerial(cases.getFinalGeneVariantPositions(),
																																																												log));
		log.reportTimeInfo("Computing density enrichments  for comp "
											 + ext.rootOf(cases.getFinalGeneVariantPositions()) + " Vs "
											 + controls.getFinalGeneVariantPositions());

		Hashtable<String, GeneVariantPositionSummary> controlSummaries = parseAvailable(GeneVariantPositionSummary.readSerial(controls.getFinalGeneVariantPositions(),
																																																													log));
		int index = 0;
		for (String key : casesSummaries.keySet()) {
			if (index % 1000 == 0) {
				log.reportTimeInfo("On index " + index + ", gene " + key + " of "
													 + casesSummaries.keySet().size() + " for comp "
													 + ext.rootOf(cases.getFinalGeneVariantPositions()) + " Vs "
													 + controls.getFinalGeneVariantPositions());
			}
			index++;
			GeneVariantPositionSummary currentCase = casesSummaries.get(key);
			PosCluster currentCaseClusReg = new PosCluster(currentCase.variantAAPositions);
			PosCluster currentCaseClusHQ = new PosCluster(currentCase.hqVariantAAPositions);

			GeneVariantPositionSummary currentControl = null;
			if (controlSummaries.containsKey(key)) {
				currentControl = controlSummaries.get(key);
				currentCaseClusReg.computeBelongs(currentControl.variantAAPositions);
				currentCaseClusHQ.computeBelongs(currentControl.hqVariantAAPositions);
			}
			cluster.put(key, new PosCluster[] {currentCaseClusReg, currentCaseClusHQ});
		}
		return cluster;
	}

	private static void summarizeVariantsBySample(SimpleTallyResult sr, Set<String> lqs, Logger log) {

		try {

			BufferedReader reader = Files.getAppropriateReader(sr.getFinalAnnotSample());
			int snpEFFIndex = ext.indexOfStr("SNPEFF_IMPACT",
																			 Files.getHeaderOfFile(sr.getFinalAnnotSample(), log));
			int sampIndex = ext.indexOfStr("SAMPLE",
																		 Files.getHeaderOfFile(sr.getFinalAnnotSample(), log));
			int geneIndex = ext.indexOfStr("SNPEFF_GENE_NAME",
																		 Files.getHeaderOfFile(sr.getFinalAnnotSample(), log));
			int hqIndex = ext.indexOfStr("HQ", Files.getHeaderOfFile(sr.getFinalAnnotSample(), log));

			Hashtable<String, Integer> counts = new Hashtable<String, Integer>();
			Hashtable<String, Integer> countsGene = new Hashtable<String, Integer>();

			reader.readLine();
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				// writer.println("Sample\tSNPEFF_IMPACT\tQUALITY\tGENE\tCOUNTS");

				String anyKey = line[sampIndex] + "\tANY_EFF\tANY";
				String funcKey = line[sampIndex] + "\t" + line[snpEFFIndex] + "\tANY";
				String geneAnyKey = line[sampIndex] + "\tANY_EFF\tANY\t" + line[geneIndex];
				String geneFuncKey = line[sampIndex] + "\t" + line[snpEFFIndex] + "\tANY\t"
														 + line[geneIndex];

				String anyKeyHQ = line[sampIndex] + "\tANY_EFF\tHQ";
				String funcKeyHQ = line[sampIndex] + "\t" + line[snpEFFIndex] + "\tHQ";
				String geneAnyKeyHQ = line[sampIndex] + "\tANY_EFF\tHQ\t" + line[geneIndex];
				String geneFuncKeyHQ = line[sampIndex] + "\t" + line[snpEFFIndex] + "\tHQ\t"
															 + line[geneIndex];
				String qual = line[hqIndex];
				if (!counts.containsKey(anyKey)) {
					counts.put(anyKey, 1);
				} else {
					counts.put(anyKey, counts.get(anyKey) + 1);
				}
				if (!countsGene.containsKey(geneAnyKey)) {
					countsGene.put(geneAnyKey, 1);
				} else {
					countsGene.put(geneAnyKey, countsGene.get(geneAnyKey) + 1);
				}

				if (!counts.containsKey(funcKey)) {
					counts.put(funcKey, 1);
				} else {
					counts.put(funcKey, counts.get(funcKey) + 1);
				}

				if (!countsGene.containsKey(geneFuncKey)) {
					countsGene.put(geneFuncKey, 1);
				} else {
					countsGene.put(geneFuncKey, countsGene.get(geneFuncKey) + 1);
				}
				if (qual.contains("Passed")) {
					if (!counts.containsKey(anyKeyHQ)) {
						counts.put(anyKeyHQ, 1);
					} else {
						counts.put(anyKeyHQ, counts.get(anyKeyHQ) + 1);
					}
					if (!countsGene.containsKey(geneAnyKeyHQ)) {
						countsGene.put(geneAnyKeyHQ, 1);
					} else {
						countsGene.put(geneAnyKeyHQ, countsGene.get(geneAnyKeyHQ) + 1);
					}
					if (!counts.containsKey(funcKeyHQ)) {
						counts.put(funcKeyHQ, 1);
					} else {
						counts.put(funcKeyHQ, counts.get(funcKeyHQ) + 1);
					}

					if (!countsGene.containsKey(geneFuncKeyHQ)) {
						countsGene.put(geneFuncKeyHQ, 1);
					} else {
						countsGene.put(geneFuncKeyHQ, countsGene.get(geneFuncKeyHQ) + 1);
					}
				}

			}
			reader.close();
			String out = sr.getSampleVarCounts();

			try {
				PrintWriter writer = Files.openAppropriateWriter(out);
				writer.println("Sample\tSNPEFF_IMPACT\tQUALITY\tGENE\tCOUNTS\tHQ_Sample");
				for (String key : counts.keySet()) {
					writer.println(key + "\tANY\t" + counts.get(key) + "\t"
												 + !lqs.contains(key.split("\t")[0]));
				}
				for (String key : countsGene.keySet()) {
					writer.println(key + "\t" + countsGene.get(key) + "\t"
												 + !lqs.contains(key.split("\t")[0]));
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + out);
				log.reportException(e);
			}
			Files.copyFile(out, sr.getBundleDir() + ext.removeDirectoryInfo(out));

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + sr.getFinalAnnotSample()
											+ "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + sr.getFinalAnnotSample() + "\"");
			return;
		}

	}

	private static class SimpleTallyResult {
		private final VcfPopulation controls;
		// private String finalOut;
		private final String finalOutVCF;
		private final String finalAnnot;
		private final String finalAnnotSample;
		private final String finalsampSummary;
		private final String finalAnnotGene;
		private final String finalGeneVariantPositions;
		private final String bundleDir;
		private final String annotKeysFile;
		private final String filters;
		private final String sampleVarCounts;
		private final String finalGeneSetSummary;

		public SimpleTallyResult(VcfPopulation controls, String finalOut, String finalOutVCF,
														 String finalsampSummary, String finalAnnot, String finalAnnotSample,
														 String finalAnnotGene, String finalGeneVariantPositions,
														 String bundleDir, String annotKeysFile, String filters,
														 String sampleVarCounts, String finalGeneSetSummary) {
			super();
			this.controls = controls;
			// this.finalOut = finalOut;
			this.finalOutVCF = finalOutVCF;
			this.finalsampSummary = finalsampSummary;
			this.finalAnnot = finalAnnot;
			this.finalAnnotGene = finalAnnotGene;
			this.finalAnnotSample = finalAnnotSample;
			this.finalGeneVariantPositions = finalGeneVariantPositions;
			this.bundleDir = bundleDir;
			this.annotKeysFile = annotKeysFile;
			this.filters = filters;
			this.sampleVarCounts = sampleVarCounts;
			this.finalGeneSetSummary = finalGeneSetSummary;
		}

		public String getFinalOutVCF() {
			return finalOutVCF;
		}

		public String getSampleVarCounts() {
			return sampleVarCounts;
		}

		public String getFilters() {
			return filters;
		}

		public String getAnnotKeysFile() {
			return annotKeysFile;
		}

		public String getBundleDir() {
			return bundleDir;
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

	private static SimpleTallyResult runSimpleTally(String vcf, String vpop, double maf,
																									int numThreads, String outDir, GeneSet[] geneSets,
																									VariantContextFilter qualCase,
																									HashSet<String> lqs, HashSet<String> genelimit,
																									boolean splitPopVcf,
																									Logger log) {
		VcfPopulation vpopAc = VcfPopulation.load(vpop, POPULATION_TYPE.ANY, log);
		vpopAc.report();
		String caseDef = ext.rootOf(vpop);
		ChrSplitResults[] chrSplitResults = VCFOps.splitByChrs(vcf, outDir, numThreads, false, log);
		WorkerHive<String> hive = new WorkerHive<String>(numThreads, 10, log);
		ArrayList<String> filtVcfs = new ArrayList<String>();

		for (ChrSplitResults chrSplitResult : chrSplitResults) {// filter each chromosome
			String outputVcf = outDir + ext.rootOf(vpop) + ".maf_" + maf + "." + chrSplitResult.getChr()
												 + VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral();
			filtVcfs.add(outputVcf);
			FilterWorker worker = new FilterWorker(chrSplitResult.getOutputVCF(), outputVcf,
																						 ext.rootOf(vpop),
																						 VcfPopulation.load(vpop, POPULATION_TYPE.ANY, log),
																						 maf, log);
			hive.addCallable(worker);
		}
		hive.execute(true);

		String finalOut = outDir + ext.rootOf(vpop) + ".maf_" + maf + ".final";
		String finalOutVCF = finalOut + VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral();
		String finalAnnot = finalOut + ".summary";
		String finalAnnotSample = finalOut + ".summary.sample";
		String finalsampSummary = finalOut + ".sampSummary.txt";
		String finalAnnotGene = finalOut + ".gene";
		String finalAnnotGeneBed = finalOut + ".gene.bed";
		String finalGeneSetSummary = finalOut + ".geneset.summmary";
		String finalGeneVariantPositions = finalOut + ".gene.position.counts.ser";
		String annotKeysFile = finalOut + ".annotation.keys.txt";
		String filterFile = finalOut + "filters.applied.txt";
		String varCounts = finalAnnotSample + ".varCounts";
		// String finalAnnotGeneSample = finalOut + ".gene.sample";

		Set<String> cases = vpopAc.getSuperPop().get(caseDef);
		Set<String> hqCases = new HashSet<String>();

		for (String removeCase : lqs) {
			if (!cases.contains(removeCase)) {
				throw new IllegalArgumentException("Invalid case to remove " + removeCase
																					 + ", must be present in actual case file");
			}
		}
		for (String acase : cases) {
			if (!lqs.contains(acase)) {
				hqCases.add(acase);
			}
		}
		String hqCaseDef = hqCases.size() == cases.size() ? caseDef : "HIGH_Q_SAMPLE_" + caseDef;
		log.reportTimeInfo(cases.size() + " total cases, " + hqCases.size() + " HQ cases");
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
		Hashtable<String, Set<String>> controlSubPop = new Hashtable<String, Set<String>>();

		for (String aControl : controlPop.get(caseDef + "_" + VcfPopulation.CONTROL)) {
			if (!controlSubPop.containsKey(caseDef + "_"
																		 + vpopAc.getPopulationForInd(aControl,
																																	RETRIEVE_TYPE.SUB)[0])) {
				Set<String> tmpSetSub = new HashSet<String>();
				controlSubPop.put(caseDef + "_"
													+ vpopAc.getPopulationForInd(aControl, RETRIEVE_TYPE.SUB)[0], tmpSetSub);
			}
			controlSubPop.get(caseDef + "_" + vpopAc.getPopulationForInd(aControl, RETRIEVE_TYPE.SUB)[0])
									 .add(aControl);
		}
		VcfPopulation controlVcfPopulation = new VcfPopulation(controlSubPop, controlPop,
																													 POPULATION_TYPE.CASE_CONTROL,
																													 new Logger());
		String bundleDir = outDir + "bundle_" + caseDef + "maf_" + maf + "/";
		new File(bundleDir).mkdirs();
		SimpleTallyResult simpleTallyResult = new SimpleTallyResult(controlVcfPopulation, finalOut,
																																finalOutVCF, finalsampSummary,
																																finalAnnot, finalAnnotSample,
																																finalAnnotGene,
																																finalGeneVariantPositions,
																																bundleDir, annotKeysFile,
																																filterFile, varCounts,
																																finalGeneSetSummary);
		simpleTallyResult.getFinalGeneVariantPositions();

		String[][] genotypeAnnotations = GenotypeOps.getGenoFormatKeys(vcf, log);
		String[][] variantAnnotations = VCFOps.getAnnotationKeys(vcf, log);
		String annos = "CONTEXT\tID\tDescription\n";
		for (int j = 0; j < genotypeAnnotations[0].length; j++) {
			annos += "GENOTYPE\t" + genotypeAnnotations[0][j] + "\t" + genotypeAnnotations[1][j] + "\n";
		}
		for (int j = 0; j < variantAnnotations[0].length; j++) {
			annos += "VARIANT\t" + variantAnnotations[0][j] + "\t" + variantAnnotations[1][j] + "\n";
		}

		Files.write(annos, annotKeysFile);

		summarizeAnalysisParams(finalsampSummary, caseDef, cases, lqs, controls, maf, log);
		summarizeGeneSets(geneSets, finalGeneSetSummary, log);
		if (qualCase != null) {
			log.reportTimeWarning("Using different qc filter for cases, this is appropriate for Tumor Normal analysis, but not much else currently");
		} else {
			qualCase = getQualityFilterwkggseq(maf, log);
		}
		VariantContextFilter qualControl = getQualityFilterwkggseq(maf, log);
		if (!Files.exists(finalAnnotGene) || !Files.exists(finalAnnotGeneBed)
				|| !Files.exists(finalAnnotSample) || !Files.exists(finalGeneVariantPositions)
				|| !Files.exists(filterFile)) {
			VCFFileReader tmp = new VCFFileReader(new File(filtVcfs.get(0)), true);
			summarizeQC(caseDef, filterFile, qualCase, qualControl);
			VariantContextWriter writer = VCFOps.initWriter(finalOutVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
																											VCFOps.getSequenceDictionary(tmp));
			VCFOps.copyHeader(tmp, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
			tmp.close();

			PrintWriter annoGeneWriter = Files.getAppropriateWriter(finalAnnotGene);
			PrintWriter annoGeneBedWriter = Files.getAppropriateWriter(finalAnnotGeneBed);
			annoGeneBedWriter.println("CHR\tSTART\tSTOP\tGENENAME_FUNCTION");
			annoGeneWriter.print(ArrayUtils.toStr(GENE_BASE));

			PrintWriter annoWriterSample = Files.getAppropriateWriter(finalAnnotSample);

			PrintWriter annoWriter = Files.getAppropriateWriter(finalAnnot);

			annoWriter.print(ArrayUtils.toStr(ANNO_BASE));
			annoWriterSample.print(ArrayUtils.toStr(ANNO_BASE_SAMPLE));
			annoWriterSample.print("\t" + ArrayUtils.toStr(genotypeAnnotations[0]));

			annoWriter.print("\t" + ArrayUtils.toStr(ArrayUtils.tagOn(ANNO_ADD,
																																hqCaseDef + "_N_" + hqCases.size(),
																																null)));
			// annoWriterSample.print("\t" + Array.toStr(Array.tagOn(ANNO_ADD, hqCaseDef + "_N_" +
			// hqCases.size(), null)));
			annoGeneWriter.print("\t"
													 + ArrayUtils.toStr(ArrayUtils.tagOn(GENE_ADD,
																															 hqCaseDef + "_N_" + hqCases.size(),
																															 null)));

			annoWriter.print("\t"
											 + ArrayUtils.toStr(ArrayUtils.tagOn(ANNO_ADD, caseDef + "_N_" + cases.size(),
																													 null)));
			annoWriterSample.print("\t"
														 + ArrayUtils.toStr(ArrayUtils.tagOn(ANNO_ADD,
																																 caseDef + "_N_" + cases.size(),
																																 null)));
			annoGeneWriter.print("\t" + ArrayUtils.toStr(ArrayUtils.tagOn(GENE_ADD,
																																		caseDef + "_N_" + cases.size(),
																																		null)));

			ArrayList<String> controlsOrdered = new ArrayList<String>();
			// general note, this is no longer a simple tally and would recommend re-doing this whthing.
			controlsOrdered.addAll(controls.keySet());
			for (int i = 0; i < controlsOrdered.size(); i++) {
				String annotLine = "\t" + ArrayUtils.toStr(ArrayUtils.tagOn(ANNO_ADD,
																																		controlsOrdered.get(i) + "_N_"
																																							+ controls.get(controlsOrdered.get(i))
																																												.size(),
																																		null));
				annoWriter.print(annotLine);
				annoWriterSample.print(annotLine);
				annoGeneWriter.print("\t" + ArrayUtils.toStr(ArrayUtils.tagOn(GENE_ADD,
																																			controlsOrdered.get(i) + "_N_"
																																								+ controls.get(controlsOrdered.get(i))
																																													.size(),
																																			null)));
			}
			annoGeneWriter.print("\tIS_GENE_SET");

			annoWriter.print("\t" + ArrayUtils.toStr(variantAnnotations[0]));
			annoWriterSample.print("\t" + ArrayUtils.toStr(variantAnnotations[0]));
			for (GeneSet geneSet : geneSets) {
				annoGeneWriter.print("\t" + geneSet.getTag() + "_Membership");
				annoWriter.print("\t" + geneSet.getTag() + "_Membership");
				annoWriterSample.print("\t" + geneSet.getTag() + "_Membership");
			}
			annoWriterSample.print("\tHQ_SAMPLE");

			annoGeneWriter.println();
			annoWriter.println();
			annoWriterSample.println();

			Hashtable<String, ArrayList<GeneSummary[]>> geneSummaries = new Hashtable<String, ArrayList<GeneSummary[]>>();
			for (int i = 0; i < filtVcfs.size(); i++) {
				log.reportTimeInfo("Summarizing " + filtVcfs.get(i));
				VCFFileReader result = new VCFFileReader(new File(filtVcfs.get(i)), true);
				for (VariantContext vc : result) {
					String geneName = VCOps.getSNP_EFFGeneName(vc);
					if (genelimit == null || genelimit.isEmpty() || genelimit.contains(geneName)) {
						writer.add(vc);
						Segment seg = VCOps.getSegment(vc);

						String func = VCOps.getSNP_EFFImpact(vc);
						annoGeneBedWriter.println(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop()
																			+ "\t" + geneName + ":" + func);
						if (!geneSummaries.containsKey(geneName)) {
							addEntries(caseDef, hqCaseDef, controlsOrdered, geneSummaries, geneName);
						}
						for (int j = 0; j < geneSets.length; j++) {
							if ((geneSets[j].getGenes().containsKey(geneName)
									 || geneSets[j].getGenes().containsKey(ANY_GENE_SET))
									&& !geneSummaries.containsKey(geneSets[j].getTag())) {
								addEntries(caseDef, hqCaseDef, controlsOrdered, geneSummaries,
													 geneSets[j].getTag());
							}
						}
						VcGroupSummary vcCaseGroup = new VcGroupSummary(caseDef, cases, vc, qualCase, log);
						VcGroupSummary vcHqCaseGroup = vcCaseGroup;
						if (hqCases.size() < cases.size()) {
							vcHqCaseGroup = new VcGroupSummary(hqCaseDef, hqCases, vc, qualCase, log);
						}

						// TODO, check
						for (int j = 0; j < geneSummaries.get(geneName).get(0).length; j++) {
							geneSummaries.get(geneName).get(0)[j].add(vcHqCaseGroup, null);
							geneSummaries.get(geneName).get(1)[j].add(vcCaseGroup, null);

							for (GeneSet geneSet : geneSets) {
								if (geneSet.getGenes().containsKey(geneName)
										|| geneSet.getGenes().containsKey(ANY_GENE_SET)) {
									geneSummaries.get(geneSet.getTag()).get(0)[j].add(vcHqCaseGroup,
																																		geneSet.getTag());
									geneSummaries.get(geneSet.getTag()).get(1)[j].add(vcCaseGroup, geneSet.getTag());
								}
							}
						}
						boolean highModLow = ext.indexOfStr(func, EFF) >= 0;
						annoWriter.print(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getID() + "\t"
														 + vc.getReference().getBaseString() + "\t"
														 + vc.getAlternateAlleles().toString() + "\t" + vc.isBiallelic() + "\t"
														 + vc.getFilters().toString() + "\t" + highModLow);
						annoWriter.print("\t" + ArrayUtils.toStr(vcHqCaseGroup.getSummary()) + "\t"
														 + ArrayUtils.toStr(vcCaseGroup.getSummary()));

						GenotypesContext gc = vcCaseGroup.getVcAlt().getGenotypes();

						ArrayList<VcGroupSummary> controlGroupSummaries = new ArrayList<VCFSimpleTally.VcGroupSummary>();
						for (int j = 0; j < controlsOrdered.size(); j++) {
							VcGroupSummary vcControlGroup = new VcGroupSummary(controlsOrdered.get(j),
																																 controls.get(controlsOrdered.get(j)),
																																 vc, qualControl, log);
							controlGroupSummaries.add(vcControlGroup);
						}

						for (Genotype g : gc) {
							HashSet<String> tmpVCSub = new HashSet<String>();
							tmpVCSub.add(g.getSampleName());
							VariantContextFilterPass pass = qualCase.filter(VCOps.getSubset(vc, tmpVCSub));
							annoWriterSample.print(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getID()
																		 + "\t"
																		 + vc.getReference().getBaseString() + "\t"
																		 + vc.getAlternateAlleles().toString() + "\t" + vc.isBiallelic()
																		 + "\t" + vc.getFilters().toString() + "\t" + highModLow + "\t"
																		 + g.getSampleName() + "\t" + g.toString() + "\t"
																		 + pass.getTestPerformed() + "\t"
																		 + ArrayUtils.toStr(GenotypeOps.getGenoAnnotationsFor(genotypeAnnotations[0],
																																													g, ".")));
							annoWriterSample.print("\t" + ArrayUtils.toStr(vcCaseGroup.getSummary()));
							for (int j = 0; j < controlsOrdered.size(); j++) {
								annoWriterSample.print("\t"
																			 + ArrayUtils.toStr(controlGroupSummaries.get(j)
																																							 .getSummary()));
							}
							annoWriterSample.print("\t"
																		 + ArrayUtils.toStr(VCOps.getAnnotationsFor(variantAnnotations[0],
																																								vc, ".")));
							for (GeneSet geneSet : geneSets) {
								annoWriterSample.print("\t" + (geneSet.getGenes().containsKey(geneName)
																							 || geneSet.getGenes().containsKey(ANY_GENE_SET)));
							}
							annoWriterSample.print("\t" + !lqs.contains(g.getSampleName()));
							annoWriterSample.println();
						}

						for (int j = 0; j < controlsOrdered.size(); j++) {
							VcGroupSummary vcControlGroup = controlGroupSummaries.get(j);
							for (int k = 0; k < geneSummaries.get(geneName).get(0).length; k++) {
								geneSummaries.get(geneName).get(j + 2)[k].add(vcControlGroup, null);// all, hq
								for (GeneSet geneSet : geneSets) {
									if (geneSet.getGenes().containsKey(geneName)
											|| geneSet.getGenes().containsKey(ANY_GENE_SET)) {
										geneSummaries.get(geneSet.getTag()).get(j + 2)[k].add(vcControlGroup,
																																					geneSet.getTag());
									}
								}
							}
							annoWriter.print("\t" + ArrayUtils.toStr(vcControlGroup.getSummary()));
						}
						annoWriter.print("\t" + ArrayUtils.toStr(VCOps.getAnnotationsFor(variantAnnotations[0],
																																						 vc, ".")));
						for (GeneSet geneSet : geneSets) {
							annoWriter.print("\t" + (geneSet.getGenes().containsKey(geneName)
																			 || geneSet.getGenes().containsKey(ANY_GENE_SET)));
						}
						annoWriter.println();
					}
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
					controlPos.add(geneSummariesCurrent.get(0)[i].getGeneVariantPositionSummary());// case
																																												 // only
					annoGeneWriter.print(gene + "\t"
															 + ArrayUtils.toStr(geneSummariesCurrent.get(0)[i].getEffects(),
																									"||"));
					for (int j = 0; j < geneSummariesCurrent.size(); j++) {
						annoGeneWriter.print("\t"
																 + ArrayUtils.toStr(geneSummariesCurrent.get(j)[i].getSummary()));
					}

					annoGeneWriter.print("\t" + isGeneSet(geneSets, gene));
					for (GeneSet geneSet : geneSets) {
						annoGeneWriter.print("\t" + (geneSet.getGenes().containsKey(gene)
																				 || geneSet.getGenes().containsKey(ANY_GENE_SET)));
					}
					annoGeneWriter.println();
				}
			}
			annoGeneWriter.close();
			GeneVariantPositionSummary.writeSerial(controlPos.toArray(new GeneVariantPositionSummary[controlPos.size()]),
																						 finalGeneVariantPositions, log);
			if (splitPopVcf) {
				VCFOps.VcfPopulation.splitVcfByPopulation(finalOutVCF, vpop, true, true, false, log);
			}
		} else

		{
			log.reportTimeWarning(finalAnnotGene + " exists so skipping summarize");
		}
		Files.copyFile(finalAnnot, bundleDir + ext.removeDirectoryInfo(finalAnnot));
		Files.copyFile(filterFile, bundleDir + ext.removeDirectoryInfo(filterFile));
		Files.copyFile(annotKeysFile, bundleDir + ext.removeDirectoryInfo(annotKeysFile));
		Files.copyFile(finalAnnotGene, bundleDir + ext.removeDirectoryInfo(finalAnnotGene));
		Files.copyFile(finalsampSummary, bundleDir + ext.removeDirectoryInfo(finalsampSummary));
		Files.copyFile(finalGeneSetSummary, bundleDir + ext.removeDirectoryInfo(finalGeneSetSummary));
		Files.copyFile(finalAnnotSample, bundleDir + ext.removeDirectoryInfo(finalAnnotSample));
		return simpleTallyResult;

	}

	private static void summarizeQC(String controlGroup, String output,
																	VariantContextFilter caseFilter,
																	VariantContextFilter controlFilter) {
		ArrayList<String> filterSummary = new ArrayList<String>();
		filterSummary.add("GROUP\tTYPE\tID\tVALUE\tDIRECTION");
		for (int i = 0; i < caseFilter.getvBooleans().length; i++) {
			VcFilterBoolean vb = caseFilter.getvBooleans()[i];
			filterSummary.add(controlGroup + "\tBOOLEAN\t" + vb.getBfilter() + "\t"
												+ vb.getBfilter().getType() + "\tNA");
		}
		for (int i = 0; i < caseFilter.getvDoubles().length; i++) {
			VcFilterDouble vd = caseFilter.getvDoubles()[i];
			String add = controlGroup + "\tNUMERIC\t" + vd.getDfilter() + "\t" + vd.getFilterThreshold()
									 + "\t" + vd.getDfilter().getType();
			filterSummary.add(add);
		}

		if (caseFilter.getvFilterJEXL().getjExps().size() > 0) {
			filterSummary.add(controlGroup + "\tSTRING\t"
												+ caseFilter.getvFilterJEXL().getjExps().get(0).name + "\t"
												+ caseFilter.getvFilterJEXL().getjExps().get(0).exp + "\tNA");
		}

		for (int i = 0; i < controlFilter.getvBooleans().length; i++) {
			VcFilterBoolean vb = controlFilter.getvBooleans()[i];
			filterSummary.add("CONTROL\tBOOLEAN\t" + vb.getBfilter() + "\t" + vb.getBfilter().getType()
												+ "\tNA");
		}
		for (int i = 0; i < controlFilter.getvDoubles().length; i++) {
			VcFilterDouble vd = controlFilter.getvDoubles()[i];
			filterSummary.add("CONTROL\tNUMERIC\t" + vd.getDfilter() + "\t" + vd.getFilterThreshold()
												+ "\t" + vd.getDfilter().getType());
		}
		if (controlFilter.getvFilterJEXL().getjExps().size() > 0) {
			filterSummary.add("CONTROL\tSTRING\t" + controlFilter.getvFilterJEXL().getjExps().get(0).name
												+ "\t" + controlFilter.getvFilterJEXL().getjExps().get(0).exp + "\tNA");
		}
		Files.writeArray(ArrayUtils.toStringArray(filterSummary), output);
	}

	private static boolean isGeneSet(GeneSet[] geneSets, String tag) {
		for (GeneSet geneSet : geneSets) {
			if (geneSet.getTag().equals(tag)) {
				return true;
			}
		}
		return false;
	}

	private static void summarizeGeneSets(GeneSet[] geneSets, String sumFile, Logger log) {
		try {
			PrintWriter writer = Files.openAppropriateWriter(sumFile);
			writer.println("GENE_SET\tN_GENES\tGENES->");
			for (GeneSet geneSet : geneSets) {
				writer.print(geneSet.getTag() + "\t" + geneSet.getGenes().size());
				for (String gene : geneSet.getGenes().keySet()) {
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

	private static void summarizeAnalysisParams(String sumFile, String caseDef, Set<String> cases,
																							Set<String> lowerQualityCases,
																							Hashtable<String, Set<String>> controls, double maf,
																							Logger log) {
		try {
			PrintWriter writer = Files.openAppropriateWriter(sumFile);
			writer.print("#CASE\t" + caseDef + "\tn=" + cases.size());
			for (String acase : cases) {
				writer.print("\t" + acase);
			}
			writer.println();
			writer.print("#CASE_FLAGGED_AS_LOWER_QUALITY\t" + caseDef + "\tn="
									 + lowerQualityCases.size());
			for (String lqc : lowerQualityCases) {
				writer.print("\t" + lqc);
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

	private static void addEntries(String caseDef, String hqCaseDef,
																 List<String> controlsOrdered,
																 Hashtable<String, ArrayList<GeneSummary[]>> geneSummaries,
																 String geneName) {
		geneSummaries.put(geneName, new ArrayList<GeneSummary[]>());

		GeneSummary[] hqcaseGeneSummaries = new GeneSummary[EFF_DEFS.length];
		for (int j = 0; j < hqcaseGeneSummaries.length; j++) {
			hqcaseGeneSummaries[j] = new GeneSummary(geneName, hqCaseDef, EFF_DEFS[j]);
		}
		geneSummaries.get(geneName).add(hqcaseGeneSummaries);

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

		private final String name;
		private final String group;
		private final String key;
		private final ArrayList<Integer> variantNPositions;
		private final ArrayList<Integer> hqVariantNPositions;
		private final ArrayList<Integer> variantAAPositions;
		private final ArrayList<Integer> hqVariantAAPositions;
		private final String[] effects;

		public GeneVariantPositionSummary(String name, String group, String[] effects) {
			super();
			this.name = name;
			this.group = group;
			this.effects = effects;
			variantNPositions = new ArrayList<Integer>();
			variantAAPositions = new ArrayList<Integer>();
			hqVariantNPositions = new ArrayList<Integer>();
			hqVariantAAPositions = new ArrayList<Integer>();
			key = name + "_" + ArrayUtils.toStr(effects, "||");
		}

		public String getKey() {
			return key;
		}

		public String[] getEffects() {
			return effects;
		}

		public static void writeSerial(GeneVariantPositionSummary[] geneVariantSummaries,
																	 String filename, Logger log) {
			SerializedFiles.writeSerial(geneVariantSummaries, filename, true);
		}

		public static GeneVariantPositionSummary[] readSerial(String filename, Logger log) {
			return (GeneVariantPositionSummary[]) SerializedFiles.readSerial(filename, false, log, false,
																																			 true);
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
			String aaChange = VCOps.getAnnotationsFor(new String[] {"SNPEFF_AMINO_ACID_CHANGE"}, vc,
																								".")[0];
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
		private final GeneVariantPositionSummary geneVariantPositionSummary;
		private final String geneName;
		private final String group;
		private final String[] effects;
		private int numVar;
		private int numMut;
		private final HashSet<String> uniqueIndsWithVar;
		private int hqNumVar;
		private int hqNumMut;
		private final HashSet<String> uniqueHqIndsWithVar;
		private final Hashtable<String, Integer> numVarsPerInd;
		private final Hashtable<String, Integer> numVarsDiffHaploPerInd;

		private final Hashtable<String, Integer> numHQVarsPerInd;
		private final Hashtable<String, Integer> numHQDiffHaploVarsPerInd;

		private final Hashtable<String, ArrayList<String>> sampleHQHaplotype;
		private final Hashtable<String, ArrayList<String>> sampleHaplotype;

		private GeneSummary(String geneName, String group, String[] effects) {
			super();
			this.group = group;
			this.geneName = geneName;
			numVar = 0;
			numMut = 0;
			hqNumVar = 0;
			hqNumMut = 0;

			this.effects = effects;
			uniqueHqIndsWithVar = new HashSet<String>();
			uniqueIndsWithVar = new HashSet<String>();
			numVarsPerInd = new Hashtable<String, Integer>();
			numHQVarsPerInd = new Hashtable<String, Integer>();
			numHQDiffHaploVarsPerInd = new Hashtable<String, Integer>();
			numVarsDiffHaploPerInd = new Hashtable<String, Integer>();
			sampleHaplotype = new Hashtable<String, ArrayList<String>>();
			sampleHQHaplotype = new Hashtable<String, ArrayList<String>>();
			geneVariantPositionSummary = new GeneVariantPositionSummary(geneName, group, effects);
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
			summary.add(numMut + "");

			summary.add(uniqueIndsWithVar.size() + "");
			summary.add(numGreaterThan(numVarsPerInd, 2) + "");
			summary.add(numGreaterThan(numVarsDiffHaploPerInd, 2) + "");

			summary.add(hqNumVar + "");
			summary.add(hqNumMut + "");

			summary.add(uniqueHqIndsWithVar.size() + "");
			summary.add(numGreaterThan(numHQVarsPerInd, 2) + "");
			summary.add(numGreaterThan(numHQDiffHaploVarsPerInd, 2) + "");

			return ArrayUtils.toStringArray(summary);
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

				throw new IllegalArgumentException("Mismatched group names: should be " + group
																					 + " but actually " + vcGroupSummary.getGroupName());

			}
			if (!VCOps.getSNP_EFFGeneName(vcGroupSummary.getVcOriginal()).equals(geneName)
					&& (setTag != null && !setTag.equals(geneName))) {
				throw new IllegalArgumentException("Mismatched gene/geneSet names");
			} else {
				String impact = VCOps.getSNP_EFFImpact(vcGroupSummary.getVcOriginal());
				if (ext.indexOfStr(impact, effects) >= 0 || effects[effects.length - 1].equals(ANY_EFF)) {
					if (vcGroupSummary.getIndsWithAlt().size() > 0) {
						numVar++;
						numMut += vcGroupSummary.getIndsWithAlt().size();

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
							hqNumMut += vcGroupSummary.getHqIndsWithAlt().size();
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
		private final String groupName;
		private final Set<String> group;
		private final VariantContext vcOriginal;
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
		private final VariantContextFilter filter;
		private VariantContext vcSub;
		private final Logger log;

		public VcGroupSummary(String groupName, Set<String> group, VariantContext vcOriginal,
													VariantContextFilter filter, Logger log) {
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
			vcSub = VCOps.getSubset(vcOriginal, group, VC_SUBSET_TYPE.SUBSET_STRICT, false);
			avgGq = VCOps.getAvgGenotypeInfo(vcSub, null, GENOTYPE_INFO.GQ, log);
			avgDp = VCOps.getAvgGenotypeInfo(vcSub, null, GENOTYPE_INFO.DP, log);
			vcAlt = VCOps.getAltAlleleContext(vcSub, null, null, ALT_ALLELE_CONTEXT_TYPE.ALL, false, log);
			aac = (int) VCOps.getAAC(vcSub, null);
			indsWithAlt = vcAlt.getSampleNames();
			numIndsAlt = indsWithAlt.size();
			numWithCalls = vcSub.getSampleNames().size() - vcSub.getNoCallCount();
			vcAltHq = VCOps.getAltAlleleContext(vcSub, null, filter, ALT_ALLELE_CONTEXT_TYPE.ALL, false,
																					log);
			hqAac = (int) VCOps.getAAC(vcAltHq, null);
			hqIndsWithAlt = vcAltHq.getSampleNames();
			numHqIndsAlt = hqIndsWithAlt.size();
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
			return ArrayUtils.toStringArray(summary);
		}
	}

	private static class FilterWorker implements Callable<String> {
		private final String vcf;
		private final String outputVcf;
		private final String casePop;
		private final VcfPopulation vpop;
		private final double maf;
		private final Logger log;

		public FilterWorker(String vcf, String outputVcf, String casePop, VcfPopulation vpop,
												double maf, Logger log) {
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
	 * @return the {@link htsjdk.variant.variantcontext.filter.VariantContextFilter} used for genotype
	 *         qc
	 */
	private static VariantContextFilter getQualityFilterwkggseq(double maf, Logger log) {
		// VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
		gq.setDFilter(50);
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		dp.setDFilter(10);
		// VARIANT_FILTER_DOUBLE altD = VARIANT_FILTER_DOUBLE.ALT_ALLELE_DEPTH;
		// altD.setDFilter(3);

		// VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		// VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
		// VARIANT_FILTER_BOOLEAN amb = VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER;
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		// VARIANT_FILTER_BOOLEAN[] bQualFilts = new VARIANT_FILTER_BOOLEAN[] { amb };

		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] {dp, gq};

		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts,
																																	 new VARIANT_FILTER_BOOLEAN[] {fail},
																																	 new String[] {"G1000Freq"},
																																	 new String[] {FilterNGS.getPopFreqFilterString(maf)},
																																	 log);
		return vContextFilter;
	}

	public static void test(String vcf, String[] vpopsCase, String omimDir,
													String[] otherGenesOfInterest, String genesetDir, double maf,
													boolean controlSpecifiComp, boolean controlSpecifiEnrich,
													VariantContextFilter caseQualFilter, boolean doclustering) {
		// popDir + "CUSHING_FREQ.vpop", popDir + "EPP.vpop" };
		// ,popDir + "ALL_CONTROL_EPP.vpop", popDir + "ANIRIDIA.vpop", popDir + "ANOTIA.vpop" };
		int numThreads = 24;
		for (int i = 0; i < vpopsCase.length; i++) {

			ArrayList<ExcelConversionParams> excelFile = new ArrayList<ExcelConversionParams>();

			String outDir = ext.parseDirectoryOfFile(vpopsCase[i]);
			new File(outDir).mkdirs();
			Logger log = new Logger(ext.rootOf(vpopsCase[i], false) + ".log");
			// if (i > 0) {
			// log.reportTimeWarning("JOHN remember break");
			// break;
			// }
			log.reportTimeWarning("Now loading all .geneset");
			GeneSet[] currentSets = GeneSet.load(Files.listFullPaths(genesetDir == null ? ext.parseDirectoryOfFile(vpopsCase[i])
																																									: genesetDir,
																															 ".geneset", false),
																					 ext.rootOf(vpopsCase[i]), log);
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
			String lowerQualitySamples = ext.rootOf(vpopsCase[i], false) + ".lowerQualitySamples.txt";
			HashSet<String> lqs = new HashSet<String>();
			if (Files.exists(lowerQualitySamples)) {
				log.reportTimeInfo("Loading lower quality samples from " + lowerQualitySamples);
				String[] lq = HashVec.loadFileToStringArray(lowerQualitySamples, false, new int[] {0},
																										true);
				for (String element : lq) {
					lqs.add(element);
				}
			} else {
				log.reportError("JOHN you could probably remove this");
			}

			String geneLimiters = ext.rootOf(vpopsCase[i], false) + ".geneLimiter.txt";
			HashSet<String> genelimit = new HashSet<String>();
			if (Files.exists(geneLimiters)) {
				log.reportTimeInfo("Loading gene limiters " + geneLimiters);
				genelimit = HashVec.loadFileToHashSet(geneLimiters, false);

			} else {
				log.reportError("JOHN you could probably remove this");
			}
			log.reportTimeInfo("Loaded " + genelimit.size() + " gene limiters");
			SimpleTallyResult caseResult = runSimpleTally(vcf, vpopsCase[i], maf, numThreads, outDir,
																										currentSets, caseQualFilter, lqs, genelimit,
																										false, log);

			summarizeVariantsBySample(caseResult, lqs, log);
			VcfPopulation controls = caseResult.getControls();
			String controlFile = ext.parseDirectoryOfFile(vpopsCase[i])
													 + controls.getUniqSuperPop().get(0) + ".vpop";
			controls.report();
			controls.dump(controlFile);
			SimpleTallyResult controlResult = runSimpleTally(vcf, controlFile, maf, numThreads, outDir,
																											 currentSets, null, new HashSet<String>(),
																											 genelimit,
																											 false, log);
			VCFOps.VcfPopulation.splitVcfByPopulation(controlResult.getFinalOutVCF(), vpopsCase[i], true,
																								true, false, log);
			String geneFileCase = caseResult.getFinalAnnotGene();
			String geneFileControl = controlResult.getFinalAnnotGene();
			String caseWithControls = ext.addToRoot(geneFileCase, "_" + ext.rootOf(controlFile));

			ArrayList<Hashtable<String, PosCluster[]>> clusters = new ArrayList<Hashtable<String, PosCluster[]>>();
			ArrayList<Hashtable<String, String[]>> controlFuncHashes = new ArrayList<Hashtable<String, String[]>>();
			ArrayList<String> allControlFiles = new ArrayList<String>();
			Hashtable<String, PosCluster[]> clusterAllControls = new Hashtable<String, VCFSimpleTally.PosCluster[]>();
			if (doclustering) {
				try {
					clusterAllControls = densityEnrichment(caseResult, controlResult, log);
				} catch (IllegalStateException e) {
					log.reportError("Could not enrich");
					e.printStackTrace();
					return;
				}
			}
			clusters.add(clusterAllControls);

			Hashtable<String, String[]> controlsAllFuncHash = loadToGeneFuncHash(geneFileControl, log);
			controlFuncHashes.add(controlsAllFuncHash);
			allControlFiles.add(controlFile);

			if (controlSpecifiComp) {
				log.reportTimeWarning("Not doing individual control comparisions");

				for (String controlGroup : controls.getSubPop().keySet()) {
					log.reportTimeInfo("Generating control specific comparison for " + controlGroup);
					Hashtable<String, Set<String>> specificControls = new Hashtable<String, Set<String>>();
					specificControls.put(controlGroup, controls.getSubPop().get(controlGroup));


					VcfPopulation tmpPop = new VcfPopulation(specificControls, specificControls,
																									 POPULATION_TYPE.CASE_CONTROL, log);
					String out = ext.parseDirectoryOfFile(vpopsCase[i]) + controlGroup + ".vpop";
					tmpPop.dump(out);
					SimpleTallyResult controlSpecificResult = runSimpleTally(vcf, out, maf, numThreads,
																																	 outDir, currentSets, null,
																																	 new HashSet<String>(), genelimit,
																																	 false,
																																	 log);
					controlFuncHashes.add(loadToGeneFuncHash(controlSpecificResult.getFinalAnnotGene(), log));
					if (controlSpecifiEnrich) {
						Hashtable<String, PosCluster[]> clusterSpecific;
						try {
							clusterSpecific = densityEnrichment(caseResult, controlSpecificResult, log);
						} catch (IllegalStateException e) {
							log.reportError("Could not enrich");
							e.printStackTrace();
							return;
						}
						clusters.add(clusterSpecific);
					}
					allControlFiles.add(out);
				}
			} else {
				log.reportTimeWarning("Not doing individual control comparisions");
			}
			try {
				PrintWriter writer = Files.openAppropriateWriter(caseWithControls);

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
					writer.print(ArrayUtils.toStr(caseLine));

					for (int j = 0; j < controlFuncHashes.size(); j++) {
						if (controlFuncHashes.get(j).containsKey(key)) {
							writer.print("\t"
													 + ArrayUtils.toStr(ArrayUtils.subArray(controlFuncHashes.get(j).get(key),
																																	0, GENE_ADD.length)));// skip
																																												// non-count
																																												// data
																																												// (like
																																												// geneset
																																												// denote)
						} else {
							writer.print("\t" + ArrayUtils.toStr(blanks));
						}
					}

					if (line == 1) {
						for (int j = 0; j < clusters.size(); j++) {
							String[] baseHeader = PosCluster.getHeader(ext.rootOf(vpopsCase[i]),
																												 ext.rootOf(allControlFiles.get(j)));
							writer.print("\t" + ArrayUtils.toStr(baseHeader) + "\t"
													 + ArrayUtils.toStr(ArrayUtils.tagOn(baseHeader, "HQ_", null)));
						}
						writer.print("\tOMIM_DISORDER(S)\tOMIM_GENE_STATUS\tOMIM_GENE_Symbol(s)\tOMIM_GENE_TITLE(s)");
						if (otherGeneInfos != null) {
							for (OtherGeneInfo otherGeneInfo : otherGeneInfos) {
								writer.print("\t" + ArrayUtils.toStr(otherGeneInfo.getHeader()));
							}
						}

					} else {
						for (int j = 0; j < clusters.size(); j++) {
							if (clusters.get(j).containsKey(key)) {
								PosCluster[] tmp = clusters.get(j).get(key);
								writer.print("\t" + ArrayUtils.toStr(tmp[0].getData()) + "\t"
														 + ArrayUtils.toStr(tmp[1].getData()));
							} else {
								writer.print("\t" + ArrayUtils.toStr(blanksEnrichment));
							}
						}

						writer.print("\t");
						for (int j = 0; j < oGene.size(); j++) {
							writer.print((j == 0 ? "" : "|AdditionalOMIM")
													 + (oGene.get(j).getDisorders().equals("") ? "NA"
																																		 : oGene.get(j).getDisorders())
													 + (j == 0 ? "\t" : "|") + oGene.get(j).getStatus()
													 + (j == 0 ? "\t" : "|")
													 + ArrayUtils.toStr(oGene.get(j).getGeneSymbols(), "/")
													 + (j == 0 ? "\t" : "|") + oGene.get(j).getTitle());
						}
						if (otherGeneInfos != null) {
							for (OtherGeneInfo otherGeneInfo : otherGeneInfos) {
								String gene = caseLine[0].toUpperCase();
								if (otherGeneInfo.getAnno().containsKey(gene)) {
									writer.print("\t" + ArrayUtils.toStr(otherGeneInfo.getAnno().get(gene)));
								} else {
									String[] blank = ArrayUtils.stringArray(otherGeneInfo.getHeader().length, "NA");
									writer.print("\t" + ArrayUtils.toStr(blank));
								}
							}
						}

					}
					writer.println();
				}
				reader.close();
				writer.close();
				Files.copyFile(caseWithControls,
											 caseResult.getBundleDir() + ext.removeDirectoryInfo(caseWithControls));
				// String outXl = caseWithControls + ".xls";
				// ExcelWriter writerxl = new ExcelWriter(Array.toStringArray(filesToWrite),
				// Array.toStringArray(names), log);
				// writerxl.write(outXl);

			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + geneFileCase + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + geneFileCase + "\"");
				return;
			}
			ExcelConversionParams samples = new ExcelConversionParams(caseResult.getFinalsampSummary(),
																																"\t", "samples");
			excelFile.add(samples);
			ExcelConversionParams annotations = new ExcelConversionParams(caseResult.getAnnotKeysFile(),
																																		"\t", "Annotations");
			excelFile.add(annotations);
			ExcelConversionParams filters = new ExcelConversionParams(caseResult.getFilters(), "\t",
																																"HQ_filters");
			excelFile.add(filters);
			ExcelConversionParams varCount = new ExcelConversionParams(caseResult.getSampleVarCounts(),
																																 "\t", "samplesVarCount");
			excelFile.add(varCount);

			ExcelConversionParams var = new ExcelConversionParams(caseResult.getFinalAnnot(), "\t",
																														"Summary_Var");
			excelFile.add(var);
			ExcelConversionParams varSample = new ExcelConversionParams(caseResult.getFinalAnnotSample(),
																																	"\t", "Summary_Var_Sample");
			excelFile.add(varSample);
			ExcelConversionParams gene = new ExcelConversionParams(caseWithControls, "\t", "Gene");
			excelFile.add(gene);

			ExcelConversionParams geneSetSummary = new ExcelConversionParams(caseResult.finalGeneSetSummary,
																																			 "\t", "GeneSetSummary");
			excelFile.add(geneSetSummary);
			String output = caseResult.getBundleDir() + ext.rootOf(vpopsCase[i]) + "_summary.xlsx";
			new ExcelConverter(excelFile, output, log).convert(true);

		}

		// String[] vpopsControl = new String[] { popDir + "EPP.vpop", popDir + "ALL_CONTROL_EPP.vpop",
		// popDir + "ALL_CONTROL_ANIRIDIA.vpop", popDir + "ALL_CONTROL_ANOTIA.vpop", popDir +
		// "ANIRIDIA.vpop", popDir + "ANOTIA.vpop" };
		// for (int i = 0; i < vpopsControl.length; i++) {
		// String outDir = ext.parseDirectoryOfFile(vpopsControl[i]);
		// new File(outDir).mkdirs();
		// Logger log = new Logger(ext.rootOf(vpopsControl[i], false) + ".log");
		// runSimpleTally(vcf, vpopsControl[i], maf, numThreads, outDir, log);
		// }
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "main.vcf.gz";
		// String[] vpopsCase = new String[] { popDir + "OSTEO_OFF_INHERIT.vpop", popDir +
		// "ANIRIDIA.vpop", popDir + "ANOTIA.vpop" };
		// =
		// "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/aric_merge/vcf/Freq/";
		String popDir = "~/freq/";
		String[] vpopsCase = new String[] {popDir + "apop.vpop"};
		String omimDir = "~/Omim";
		// String omimDir = "/home/pankrat2/public/bin/ref/OMIM/";
		// String[] otherGenesOfInterest = new String[] { popDir + "SB_T1.txt", popDir + "SB_T2.txt",
		// "/home/pankrat2/public/bin/ref/COSMIC/cancer_gene_census.txt" };
		String[] otherGenesOfInterest = null;
		boolean controlSpecifiComp = false;
		boolean controlSpecifiEnrich = false;
		boolean doclustering = false;
		double[] mafs = new double[] {0.01};

		String usage = "\n" + "seq.analysis.VCFSimpleTally requires 0-1 arguments\n"
									 + "   (1) vcf (i.e. vcf=" + vcf + " (default))\n"
									 + "   (2) population directory (i.e. popDir=" + popDir + " (default))\n"
									 + "   (3) comma-delimited .vpop files in the popDir (i.e. vpops="
									 + ArrayUtils.toStr(vpopsCase, ",") + " (default))\n"
									 + "   (4) omim directory (i.e. omim=" + omimDir + " (default))\n"
									 + "   (5) comma-delimited extra-info files (i.e. extra= (no default))\n"
									 + "   (6) comma separated mafs  (i.e. maf="
									 + ArrayUtils.toStr(ArrayUtils.toStringArray(mafs), ",") + " (default))\n"
									 + "   (7) for each control group, run a specific tally  (i.e. controlSpecifiComp="
									 + controlSpecifiComp + " (default))\n"
									 + "   (8) for each control group, do enrichment (i.e. controlSpecifiEnrich="
									 + controlSpecifiEnrich + " (default))\n" +

									 "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("vcf=")) {
				vcf = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("popDir=")) {
				popDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("omim=")) {
				omimDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("extra=")) {
				otherGenesOfInterest = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("controlSpecifiComp=")) {
				controlSpecifiComp = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("controlSpecifiEnrich=")) {
				controlSpecifiEnrich = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("maf=")) {
				mafs = ArrayUtils.toDoubleArray(ext.parseStringArg(arg, "NA").split(","));
				numArgs--;
			} else if (arg.startsWith("vpops=")) {
				vpopsCase = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("-doclustering")) {
				doclustering = true;
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
			vpopsCase = ArrayUtils.tagOn(vpopsCase, popDir, null);
			for (String vpops : vpopsCase) {
				for (double maf : mafs) {
					test(vcf, new String[] {vpops}, omimDir, otherGenesOfInterest, null, maf,
							 controlSpecifiComp,
							 controlSpecifiEnrich, null, doclustering);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
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
