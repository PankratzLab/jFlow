package seq.analysis;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.pathway.GenomeRegions;
import seq.pathway.Pathway;
import seq.qc.FilterNGS;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

/**
 * @author lane0212 A somewhat specific class that sums gene counts in a case/control fashion and compares to charge mafs
 */
public class VCFTally implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String vcf;
	protected GenomeRegions genomeRegions;
	protected VcfPopulation vpop;
	protected TallyTracker[] trackersCase, trackersControl, trackersCharge;
	protected Logger log;
	private static final String AND = "&&";
	private static final String[] SNPEFF_IMPACTS = { "(SNPEFF_EFFECT=='STOP_LOST'||SNPEFF_EFFECT=='STOP_GAINED')", "SNPEFF_IMPACT=='HIGH'", "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE')", "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE'||SNPEFF_IMPACT=='LOW')", "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE'||SNPEFF_IMPACT=='LOW'||SNPEFF_IMPACT=='MODIFIER')" };
	protected static final String[] SNPEFF_NAMES = { "SNPEFF_STOPGAIN_LOSS", "SNPEFF_HIGH", "SNPEFF_HIGH_MODERATE", "SNPEFF_HIGH_MODERATE_LOW", "SNPEFF_HIGH_MODERATE_LOW_MODIFIER", "NO_SNPEFF" };
	private static final String ESP_FILTER = "(esp6500si_all=='.'||esp6500si_all <= 0.01)";
	private static final String G1000_FILTER = "(g10002014oct_all=='.'||g10002014oct_all <= 0.01)";
	protected CASE_CONTROL_TYPE type;

	public enum CASE_CONTROL_TYPE {
		BOTH_PASS, ONE_PASS;
	}

	public VCFTally(String vcf, GenomeRegions genomeRegions, VcfPopulation vpop, CASE_CONTROL_TYPE type, Logger log) {
		super();
		this.vcf = vcf;
		this.genomeRegions = genomeRegions;
		this.vpop = vpop;
		this.log = log;
		this.type = type;
	}

	public VCFTally(String vcf, GenomeRegions genomeRegions, VcfPopulation vpop, TallyTracker[] trackersCase, TallyTracker[] trackersControl, TallyTracker[] trackersCharge, Logger log, CASE_CONTROL_TYPE type) {
		super();
		this.vcf = vcf;
		this.genomeRegions = genomeRegions;
		this.vpop = vpop;
		this.trackersCase = trackersCase;
		this.trackersControl = trackersControl;
		this.trackersCharge = trackersCharge;
		this.log = log;
		this.type = type;
	}

	public void tallyCaseControlVCF(FilterNGS filterNGS, String outputList) {

		if (vpop != null) {
			VariantContextFilter totalQuality = initializeTrackers();

			int index = 0;

			VCFFileReader reader = new VCFFileReader(vcf, true);
			Set<String> cases = vpop.getSubPop().get(VCFOps.VcfPopulation.CASE);
			Set<String> controls = vpop.getSubPop().get(VCFOps.VcfPopulation.CONTROL);
			log.reportTimeInfo(cases.size() + " Cases and " + controls.size() + " Controls");
			Set<String> all = new HashSet<String>();
			all.addAll(cases);
			all.addAll(controls);
			if (cases.size() >= controls.size()) {
				log.reportTimeError("Since the case/control ratio is gte 0.5, a filtering procedure to remove minor alleles that are not alternate alleles may exclude informative variants, please add more controls");
			}
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(outputList, false));
				for (VariantContext vc : reader) {
					index++;
					VariantContext vcPop = VCOps.getSubset(vc, all, VC_SUBSET_TYPE.SUBSET_STRICT, false);
					if (VCOps.isMinorAlleleAlternate(vcPop, null) && VCOps.getAAC(vcPop, null) > 0) {// minor allele is alt

						VariantContext vcCase = VCOps.getSubset(vc, cases, VC_SUBSET_TYPE.SUBSET_STRICT, false);
						VariantContext vcControl = VCOps.getSubset(vc, controls, VC_SUBSET_TYPE.SUBSET_STRICT, false);

						if ((type == CASE_CONTROL_TYPE.BOTH_PASS && totalQuality.filter(vcCase).passed() && totalQuality.filter(vcControl).passed()) || (type == CASE_CONTROL_TYPE.ONE_PASS && (totalQuality.filter(vcCase).passed() || totalQuality.filter(vcControl).passed()))) {

							VariantContext vcAltCase = VCOps.getAltAlleleContext(vcCase, filterNGS, log);
							VariantContext vcAltControl = VCOps.getAltAlleleContext(vcControl, filterNGS, log);
							boolean hasCaseAlt = vcAltCase.getSampleNames().size() > 0;
							boolean hasAltControl = vcAltControl.getSampleNames().size() > 0;

							if ((type == CASE_CONTROL_TYPE.BOTH_PASS && (!hasCaseAlt || totalQuality.filter(vcAltCase).passed()) && (!hasAltControl || totalQuality.filter(vcAltControl).passed())) || (type == CASE_CONTROL_TYPE.ONE_PASS && (totalQuality.filter(vcAltCase).passed() || totalQuality.filter(vcAltControl).passed()))) {
								GeneData[] genesThatOverlap = VCOps.getGenesThatOverlap(vc, genomeRegions.getGeneTrack(), log);

								for (int i = 0; i < trackersCase.length; i++) {
									int caseAdded = trackersCase[i].addIfPasses(vcCase, genesThatOverlap, filterNGS);
									int controlAdded = trackersControl[i].addIfPasses(vcControl, genesThatOverlap, filterNGS);
									if (caseAdded + controlAdded > 0) {
										Segment vcSeg = VCOps.getSegment(vc);
										writer.println("VAR\t" + SNPEFF_NAMES[i] + "_" + type + "\t" + Positions.getChromosomeUCSC((int) vcSeg.getChr(), true) + ":" + vcSeg.getStart() + (vcSeg.getStop() == vcSeg.getStart() ? "" : ".." + vcSeg.getStop()));
									}
								}
							}
						}
					}
					// }
					if (index % 10000 == 0) {
						log.reportTimeInfo("Scanned " + index + " total variants from " + ext.removeDirectoryInfo(vcf));
						log.reportTimeInfo("Currently on " + VCOps.getSegment(vc).getUCSClocation());
						// log.reportTimeError("JOHN REMEMBER THIS CUTOFF");
						// writer.close();
						//
						// return;
					}
				}
				reader.close();
				writer.close();

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			log.reportTimeError("THIS is the case/control method, sorry");
		}
	}

	protected VariantContextFilter initializeTrackers() {
		this.trackersCase = null;
		this.trackersControl = new TallyTracker[SNPEFF_NAMES.length];
		log.reportTimeInfo("Using vpopulation to report");
		if (!vpop.getSubPop().containsKey(VCFOps.VcfPopulation.CASE)) {
			log.reportTimeError("Vpopulation must contain " + VCFOps.VcfPopulation.CASE);
		}
		if (!vpop.getSubPop().containsKey(VCFOps.VcfPopulation.CONTROL)) {
			log.reportTimeError("Vpopulation must contain " + VCFOps.VcfPopulation.CONTROL);
		}
		trackersCase = new TallyTracker[SNPEFF_NAMES.length];

		VariantContextFilter totalQuality = getQualityFilter(log);
		for (int i = 0; i < SNPEFF_IMPACTS.length; i++) {
			VariantContextFilter[] vContextFilters = new VariantContextFilter[2];
			vContextFilters[0] = getQualityFilter(log);
			vContextFilters[1] = getJEXLAt(i, log);
			trackersCase[i] = new TallyTracker("CASE_" + SNPEFF_NAMES[i], genomeRegions, vContextFilters, log);
			trackersControl[i] = new TallyTracker("CONTROL_" + SNPEFF_NAMES[i], genomeRegions, vContextFilters, log);
		}
		VariantContextFilter[] vContextFilters = new VariantContextFilter[2];
		vContextFilters[0] = getQualityFilter(log);
		vContextFilters[1] = getJEXLAt(SNPEFF_NAMES.length - 1, log);
		trackersCase[SNPEFF_NAMES.length - 1] = new TallyTracker("CASE_" + SNPEFF_NAMES[SNPEFF_NAMES.length - 1], genomeRegions, vContextFilters, log);
		trackersControl[SNPEFF_NAMES.length - 1] = new TallyTracker("CONTROL_" + SNPEFF_NAMES[SNPEFF_NAMES.length - 1], genomeRegions, vContextFilters, log);
		return totalQuality;
	}

	public void tallyCharge(String fullpathToChargeVCF) {
		this.trackersCharge = new TallyTracker[SNPEFF_NAMES.length];

		for (int i = 0; i < SNPEFF_IMPACTS.length; i++) {
			VariantContextFilter[] vContextFilters = new VariantContextFilter[1];
			vContextFilters[0] = getJEXLAt(i, log);
			trackersCharge[i] = new TallyTracker("CHARGE_" + SNPEFF_NAMES[i], genomeRegions, vContextFilters, log);
			trackersCharge[i].setCharge(true);
		}
		VariantContextFilter[] vContextFilters = new VariantContextFilter[1];
		vContextFilters[0] = getJEXLAt(SNPEFF_NAMES.length - 1, log);
		trackersCharge[SNPEFF_NAMES.length - 1] = new TallyTracker("CHARGE_" + SNPEFF_NAMES[SNPEFF_NAMES.length - 1], genomeRegions, vContextFilters, log);
		trackersCharge[SNPEFF_NAMES.length - 1].setCharge(true);
		int count = 0;
		VCFFileReader reader = new VCFFileReader(fullpathToChargeVCF, false);

		for (VariantContext vc : reader) {
			GeneData[] genesThatOverlap = VCOps.getGenesThatOverlap(vc, genomeRegions.getGeneTrack(), log);
			for (int i = 0; i < trackersCharge.length; i++) {
				trackersCharge[i].addIfPasses(vc, genesThatOverlap, null);
			}
			count++;
			if (count % 10000 == 0) {
				log.reportTimeInfo("Scanned " + count + " total variants from " + ext.removeDirectoryInfo(fullpathToChargeVCF));
				// log.reportTimeError("JOHN REMEMBER THIS CUTOFF");
				// reader.close();
				// return;
			}
		}
		reader.close();
	}

	public void summarize(String fullPathToOutput) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.print("SET\tCHR\tSTART\tSTOP\tUCSC\tTotal_length\tMrna_length");
			for (int i = 0; i < trackersCase.length; i++) {
				writer.print("\t" + trackersCase[i].getTallyName() + "_NUM_VAR" + "\t" + trackersCase[i].getTallyName() + "_NUM_INDS" + "\t" + trackersControl[i].getTallyName() + "_NUM_VAR" + "\t" + trackersControl[i].getTallyName() + "_NUM_INDS" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B_CASE_Expect");
			}
			for (int i = 0; i < trackersCase.length; i++) {
				writer.print("\t" + trackersCase[i].getTallyName() + "_ALT_ALLELE_COUNT" + "\t" + trackersCase[i].getTallyName() + "_NUM_INDS" + "\t" + trackersControl[i].getTallyName() + "_ALT_ALLELE_COUNT" + "\t" + trackersControl[i].getTallyName() + "_NUM_INDS" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B_CASE_Expect");
			}
			writer.println();
			Set<String> sets = trackersCase[0].getTally().keySet();
			for (String set : sets) {
				if (trackersCase[0].getGene(set) != null && trackersCase[0].getGene(set).length > 0) {
					String chr = Positions.getChromosomeUCSC(trackersCase[0].getGene(set)[0].getChr(), true);
					String start = trackersCase[0].getGene(set)[0].getStart() + "";
					String stop = trackersCase[0].getGene(set)[0].getStop() + "";

					writer.print(set + "\t" + chr + "\t" + start + "\t" + stop + "\t" + trackersCase[0].getGene(set)[0].getUCSCLink("hg19") + "\t" + trackersCase[0].getGeneTotalLength(set) + "\t" + trackersCase[0].getGeneTotalMrnaLength(set));
					printGene(writer, set);
				} else {
					writer.print(set + "\tNA\tNA\tNA\tNA\tNA\tNA");
					printGene(writer, set);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathToOutput);
			log.reportException(e);
		}
		String bed = fullPathToOutput + ".CASE.bed";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(bed));
			for (int i = 0; i < trackersCase.length; i++) {
				String[] tbed = trackersCase[i].getBed();
				for (int j = 0; j < tbed.length; j++) {
					writer.println(tbed[j]);
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + bed);
			log.reportException(e);
		}
	}

	// public void serialize(String filename) {
	// Files.writeSerial(this, filename);
	// }

	public static VCFTally load(String filename) {
		return (VCFTally) Files.readSerial(filename, false, false);
	}

	private void printGene(PrintWriter writer, String set) {
		for (int i = 0; i < trackersCase.length; i++) {
			int numCases = trackersCase[i].getUniqs().get(set).size();
			int numControls = trackersControl[i].getUniqs().get(set).size();
			writer.print("\t" + trackersCase[i].getTally().get(set) + "\t" + numCases + "\t" + trackersControl[i].getTally().get(set) + "\t" + numControls + "\t" + trackersCharge[i].getTallyMAC().get(set) + "\t" + ((double) vpop.getSubPop().get("CASE").size() * trackersCharge[i].getTallyMAC().get(set)));
		}
		for (int i = 0; i < trackersCase.length; i++) {
			int numCases = trackersCase[i].getUniqs().get(set).size();
			int numControls = trackersControl[i].getUniqs().get(set).size();
			writer.print("\t" + trackersCase[i].getTallyMAC().get(set) + "\t" + numCases + "\t" + trackersControl[i].getTallyMAC().get(set) + "\t" + numControls + "\t" + trackersCharge[i].getTallyMAC().get(set) + "\t" + ((double) vpop.getSubPop().get("CASE").size() * trackersCharge[i].getTallyMAC().get(set)));
		}
	}

	protected static class TallyTracker implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private Hashtable<String, Float> tally;
		private Hashtable<String, Float> tallyMAC;

		private GenomeRegions gRegions;
		private String tallyName;
		private Hashtable<Integer, int[]> geneSegTrack;
		private VariantContextFilter[] vContextFilters;
		private ArrayList<GenePass> passingLocs;
		private Hashtable<String, HashSet<String>> uniqs;
		private Hashtable<String, ArrayList<String>> all;
		private Logger log;
		private boolean charge;

		public TallyTracker(String tallyName, GenomeRegions geneTrack, VariantContextFilter[] vContextFilters, Logger log) {
			super();
			this.log = log;

			this.tallyName = tallyName;
			this.tally = getTally(geneTrack, log);
			this.tallyMAC = getTally(geneTrack, log);
			this.vContextFilters = vContextFilters;
			this.gRegions = geneTrack;
			this.charge = false;
			this.passingLocs = new ArrayList<VCFTally.GenePass>();
			this.uniqs = new Hashtable<String, HashSet<String>>();
			this.all = new Hashtable<String, ArrayList<String>>();
			populateQuery();
		}

		public TallyTracker(Hashtable<String, Float> tally, Hashtable<String, Float> tallyMAC, GenomeRegions gRegions, String tallyName, Segment[] geneSegs, Hashtable<Integer, int[]> geneSegTrack, VariantContextFilter[] vContextFilters, ArrayList<GenePass> passingLocs, Hashtable<String, HashSet<String>> uniqs, Hashtable<String, ArrayList<String>> all, Logger log, boolean charge) {
			super();
			this.tally = tally;
			this.tallyMAC = tallyMAC;
			this.gRegions = gRegions;
			this.tallyName = tallyName;
			this.geneSegTrack = geneSegTrack;
			this.vContextFilters = vContextFilters;
			this.passingLocs = passingLocs;
			this.uniqs = uniqs;
			this.all = all;
			this.log = log;
			this.charge = charge;
		}

		public void setCharge(boolean charge) {
			this.charge = charge;
		}

		private void populateQuery() {
			this.geneSegTrack = new Hashtable<Integer, int[]>();
			GeneData[][] gDatas = gRegions.getGeneTrack().getGenes();
			int index = 0;
			int curChr = 0;
			int curPos = 0;
			for (int chr = 0; chr < gDatas.length; chr++) {
				curPos = 0;
				for (int chrIndex = 0; chrIndex < gDatas[chr].length; chrIndex++) {
					if (gDatas[chr][chrIndex].getChr() < curChr || gDatas[chr][chrIndex].getStart() < curPos) {
						log.reportTimeError("Unsorted Genetrack!");
						return;
					} else {
						curChr = gDatas[chr][chrIndex].getChr();
						curPos = gDatas[chr][chrIndex].getStart();
					}
					all.put(gDatas[chr][chrIndex].getGeneName(), new ArrayList<String>());
					uniqs.put(gDatas[chr][chrIndex].getGeneName(), new HashSet<String>());
					geneSegTrack.put(index, new int[] { chr, chrIndex });
					index++;
				}
			}
			Pathway[] ways = gRegions.getPathways().getPathways();
			for (int i = 0; i < ways.length; i++) {
				uniqs.put(ways[i].getPathwayName(), new HashSet<String>());
				all.put(ways[i].getPathwayName(), new ArrayList<String>());
			}
		}

		private static int total(GeneTrack geneTrack) {
			int total = 0;
			GeneData[][] gDatas = geneTrack.getGenes();
			for (int i = 0; i < gDatas.length; i++) {
				total += gDatas[i].length;
			}
			return total;
		}

		private static Hashtable<String, Float> getTally(GenomeRegions genomeRegions, Logger log) {
			Hashtable<String, Float> tally = new Hashtable<String, Float>();
			GeneData[][] gDatas = genomeRegions.getGeneTrack().getGenes();
			for (int i = 0; i < gDatas.length; i++) {
				for (int j = 0; j < gDatas[i].length; j++) {
					tally.put(gDatas[i][j].getGeneName(), (float) 0.0);
				}
			}
			Pathway[] ways = genomeRegions.getPathways().getPathways();
			for (int i = 0; i < ways.length; i++) {
				tally.put(ways[i].getPathwayName(), (float) 0.0);
			}

			return tally;
		}

		public String getTallyName() {
			return tallyName;
		}

		public Hashtable<String, Float> getTally() {
			return tally;
		}

		public Hashtable<String, Float> getTallyMAC() {
			return tallyMAC;
		}

		public int addIfPasses(VariantContext vc, GeneData[] geneDatas, FilterNGS filterNGS) {
			int numAdded = 0;
			for (int i = 0; i < vContextFilters.length; i++) {
				if (!vContextFilters[i].filter(vc).passed()) {
					return numAdded;
				}
			}
			if (!charge && VCOps.getAAC(vc, null) <= 0) {
				return numAdded;
			}

			double mac = Double.NaN;
			if (!charge) {
				mac = VCOps.getAAC(vc, null);
			} else {
				mac = vc.getCommonInfo().getAttributeAsDouble("MAF_whites", 0.0) + vc.getCommonInfo().getAttributeAsDouble("MAF_blacks", 0.0);
			}

			String snpEffGeneName = vc.getCommonInfo().getAttributeAsString("SNPEFF_GENE_NAME", ".");
			if (mac > 0) {
				if (geneDatas.length > 0) {
					boolean foundMatch = false;
					for (int i = 0; i < geneDatas.length; i++) {
						if (geneDatas[i].getGeneName().equals(snpEffGeneName)) {
							foundMatch = true;
							Pathway[] ways = gRegions.getPathways().getPathwaysFor(geneDatas[i]);
							addVC(vc, geneDatas, filterNGS, mac, i, ways);
							numAdded++;
						}
					}
					if (!foundMatch&&!tallyName.contains("NO_SNPEFF")&&!tallyName.contains("SNPEFF_HIGH_MODERATE_LOW_MODIFIER")) {
						log.reportTimeError("Error - could not find matching SNPEFF gene for " + vc.toStringWithoutGenotypes());
						log.reportTimeError("Available overlapping genes are");
						for (int i = 0; i < geneDatas.length; i++) {
							log.reportTimeError(geneDatas[i].getGeneName());
						}
						System.exit(1);
					}
				} else if (!snpEffGeneName.equals(".") && !tallyName.contains("NO_SNPEFF") && !tallyName.contains("SNPEFF_HIGH_MODERATE_LOW_MODIFIER")) {

					log.reportTimeError("Error - could not find matching GENVISIS gene for " + vc.toStringWithoutGenotypes());
					System.exit(1);
				}
			}
			return numAdded;
		}

		private void addVC(VariantContext vc, GeneData[] geneDatas, FilterNGS filterNGS, double mac, int i, Pathway[] ways) {
			if (!charge) {
				passingLocs.add(new GenePass(VCOps.getSegment(vc), geneDatas[i].getGeneName()));

				addTally(geneDatas[i].getGeneName(), 1, (float) mac);
				// TODO alt allele context
				Set<String> samplesWithAlts = VCOps.getAltAlleleContext(vc, filterNGS, log).getSampleNames();
				uniqs.get(geneDatas[i].getGeneName()).addAll(samplesWithAlts);
				all.get(geneDatas[i].getGeneName()).addAll(samplesWithAlts);
				for (int j = 0; j < ways.length; j++) {
					uniqs.get(ways[j].getPathwayName()).addAll(samplesWithAlts);
					all.get(ways[j].getPathwayName()).addAll(samplesWithAlts);
					addTally(ways[j].getPathwayName(), 1, (float) mac);
				}
			} else {
				addTally(geneDatas[i].getGeneName(), 1, (float) mac);
				for (int j = 0; j < ways.length; j++) {
					addTally(ways[j].getPathwayName(), 1, (float) mac);
				}
			}
		}

		public String[] getBed() {
			String[] bed = new String[passingLocs.size()];
			for (int i = 0; i < bed.length; i++) {
				GenePass pass = passingLocs.get(i);
				bed[i] = Positions.getChromosomeUCSC(pass.getVarSeg().getChr(), true) + "\t" + pass.getVarSeg().getStart() + "\t" + pass.getVarSeg().getStop() + "\t" + pass.getGeneName() + ":" + tallyName;
			}
			return bed;
		}

		public void addTally(String gene, float num, float mac) {
			if (!tally.containsKey(gene)) {
				log.reportTimeError("Could not find " + gene + " in the initialized tally");
			} else {
				tally.put(gene, tally.get(gene) + num);
				tallyMAC.put(gene, tallyMAC.get(gene) + mac);
			}
		}

		public Hashtable<String, HashSet<String>> getUniqs() {
			return uniqs;
		}

		public Hashtable<String, ArrayList<String>> getAll() {
			return all;
		}

		public GeneData[] getGene(String gene) {
			return gRegions.getGeneTrack().lookupAllGeneDatas(gene);
		}

		public int getGeneTotalLength(String agene) {
			GeneData[] geneDatas = getGene(agene);
			int total = 0;
			for (int i = 0; i < geneDatas.length; i++) {
				total += geneDatas[i].getSize();
			}
			return total;
		}

		public int getGeneTotalMrnaLength(String gene) {
			int total = 0;
			GeneData[] geneDatas = getGene(gene);
			for (int i = 0; i < geneDatas.length; i++) {
				int[][] exonBoundaries = geneDatas[i].getExonBoundaries();
				for (int j = 0; j < exonBoundaries.length; j++) {
					total += exonBoundaries[j][1] - exonBoundaries[j][0];
				}
			}
			return total;
		}
	}

	private static VariantContextFilter getJEXLAt(int index, Logger log) {
		if (index < SNPEFF_IMPACTS.length) {
			return getJEXLFilter(new String[] { "RARE_" + SNPEFF_NAMES[index] }, new String[] { SNPEFF_IMPACTS[index] + AND + ESP_FILTER + AND + G1000_FILTER }, log);

		} else {
			return getJEXLFilter(new String[] { "RARE_" + SNPEFF_NAMES[index] }, new String[] { ESP_FILTER + AND + G1000_FILTER }, log);
		}
	}

	private static VariantContextFilter getJEXLFilter(String[] jexlNames, String[] jexls, Logger log) {
		VariantContextFilter vContextFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] {}, jexlNames, jexls, log);
		return vContextFilter;
	}

	private static VariantContextFilter getQualityFilter(Logger log) {
		VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_LOOSE;
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		// VARIANT_FILTER_BOOLEAN biallelic = VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER;
		// VARIANT_FILTER_BOOLEAN amb = VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER;
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;

		// VARIANT_FILTER_BOOLEAN[] bQualFilts = new VARIANT_FILTER_BOOLEAN[] { amb };
		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, dp, gq, vqslod };
		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
		return vContextFilter;
	}

	private static class GenePass implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private Segment varSeg;
		private String geneName;

		public GenePass(Segment varSeg, String geneName) {
			super();
			this.varSeg = varSeg;
			this.geneName = geneName;
		}

		public Segment getVarSeg() {
			return varSeg;
		}

		public String getGeneName() {
			return geneName;
		}
	}

	public static void test() {
		String vcf = "D:/data/Project_Tsai_Project_021/JointAnalysis/joint_genotypes.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.AgilentCaptureRegions.vcf.gz";
		String[] vpopFiles = new String[] { "D:/data/Project_Tsai_Project_021/JointAnalysis/vPopCaseControl.txt", "D:/data/Project_Tsai_Project_021/JointAnalysis/vPopCaseControlAllRaces.txt" };
		String fullpathToChargeVCF = "D:/data/CHARGE/CHARGE_MAFS/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";

		Logger log = new Logger(ext.rootOf(vcf, false) + "tally.log");
		for (int i = 0; i < vpopFiles.length; i++) {
			String outputList = ext.rootOf(vpopFiles[i], false) + ".plinkSeqVar";
			if (Files.exists(outputList)) {
				new File(outputList).delete();
			}
			for (int j = 0; j < CASE_CONTROL_TYPE.values().length; j++) {
				VcfPopulation vpop = VcfPopulation.load(vpopFiles[i], POPULATION_TYPE.CASE_CONTROL, log);
				vpop.report();
				String geneTrackFile = "N:/statgen/NCBI/RefSeq_hg19.gtrack";
				int altAlleleDepth = 0;
				VCFTally tally = new VCFTally(vcf, GenomeRegions.load(geneTrackFile), vpop, CASE_CONTROL_TYPE.values()[j], log);
				tally.tallyCaseControlVCF(null, outputList);
				tally.tallyCharge(fullpathToChargeVCF);
				tally.summarize(ext.parseDirectoryOfFile(vcf) + ext.rootOf(vpopFiles[i]) + "tallyCounts." + CASE_CONTROL_TYPE.values()[j] + ".txt");
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "VCFTally.dat";

		String usage = "\n" + "seq.analysis.VCFTally requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			test();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public String getVcf() {
		return vcf;
	}

	public GenomeRegions getGenomeRegions() {
		return genomeRegions;
	}

	public VcfPopulation getVpop() {
		return vpop;
	}

	public TallyTracker[] getTrackersCase() {
		return trackersCase;
	}

	public TallyTracker[] getTrackersControl() {
		return trackersControl;
	}

	public TallyTracker[] getTrackersCharge() {
		return trackersCharge;
	}

	public Logger getLog() {
		return log;
	}

	public static String[] getSnpeffNames() {
		return SNPEFF_NAMES;
	}

	public CASE_CONTROL_TYPE getType() {
		return type;
	}

}
