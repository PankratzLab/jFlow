package seq.analysis;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import common.Logger;
import common.Positions;
import common.ext;
import seq.manage.VCOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

/**
 * @author lane0212
 *A somewhat specific class that tallies gene counts in a case/control fashion and compares to charge mafs
 */
public class VCFTally {
	private String vcf;
	private GeneTrack geneTrack;
	private VcfPopulation vpop;
	private TallyTracker[] trackersCase, trackersControl, trackersCharge;
	private Logger log;
	private static final String AND = "&&";
	private static final String[] SNPEFF_IMPACTS = { "SNPEFF_IMPACT=='HIGH'", "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE')", "(SNPEFF_IMPACT=='HIGH'||SNPEFF_IMPACT=='MODERATE'||SNPEFF_IMPACT=='LOW')" };
	private static final String[] SNPEFF_NAMES = { "SNPEFF_HIGH", "SNPEFF_HIGH_MODERATE", "SNPEFF_HIGH_MODERATE_LOW", "NO_SNPEFF" };
	private static final String[] CHARGdE_JEXL = { "vc.getNoCallCount()<40" };
	private static final String ESP_FILTER = "(esp6500si_all=='.'||esp6500si_all <= 0.01)";
	private static final String G1000_FILTER = "(g10002014oct_all=='.'||g10002014oct_all <= 0.01)";
	private static final String CHARGE_MAF = "MAF_whites";

	public VCFTally(String vcf, GeneTrack geneTrack, VcfPopulation vpop, Logger log) {
		super();
		this.vcf = vcf;
		this.geneTrack = geneTrack;
		this.vpop = vpop;
		this.log = log;
	}

	public void tallyCaseControlVCF() {

		if (vpop != null) {
			VariantContextFilter totalQuality = initializeTrackers();

			int index = 0;

			VCFFileReader reader = new VCFFileReader(vcf, true);
			Set<String> cases = vpop.getSubPop().get("CASE");
			Set<String> controls = vpop.getSubPop().get("CONTROL");
			Set<String> all = new HashSet<String>();
			all.addAll(cases);
			all.addAll(controls);
			for (VariantContext vc : reader) {
				index++;
				VariantContext vcCase = VCOps.getSubset(vc, cases);
				VariantContext vcControl = VCOps.getSubset(vc, controls);
				if (totalQuality.filter(vcCase).passed() && totalQuality.filter(vcControl).passed()) {
					for (int i = 0; i < trackersCase.length; i++) {
						trackersCase[i].addIfPasses(vcCase);
					}
					for (int i = 0; i < trackersControl.length; i++) {
						trackersControl[i].addIfPasses(vcControl);
					}
				}
				if (index % 100000 == 0) {
					log.reportTimeInfo("Scanned " + index + " total variants from " + ext.removeDirectoryInfo(vcf));
				}
			}
			reader.close();
		} else {
			log.reportTimeError("THIS is the case/control method, sorry");
		}
	}

	private VariantContextFilter initializeTrackers() {
		this.trackersCase = null;
		this.trackersControl = new TallyTracker[SNPEFF_NAMES.length];
		log.reportTimeInfo("Using vpopulation to report");
		if (!vpop.getSubPop().containsKey("CASE")) {
			log.reportTimeError("Vpopulation must contain CASE");
		}
		if (!vpop.getSubPop().containsKey("CONTROL")) {
			log.reportTimeError("Vpopulation must contain CONTROL");
		}
		trackersCase = new TallyTracker[SNPEFF_NAMES.length];

		VariantContextFilter totalQuality = getQualityFilter(log);
		for (int i = 0; i < SNPEFF_IMPACTS.length; i++) {
			VariantContextFilter[] vContextFilters = new VariantContextFilter[2];
			vContextFilters[0] = getQualityFilter(log);
			vContextFilters[1] = getJEXLAt(i, log);
			trackersCase[i] = new TallyTracker("CASE_" + SNPEFF_NAMES[i], geneTrack, vContextFilters, log);
			trackersControl[i] = new TallyTracker("CONTROL_" + SNPEFF_NAMES[i], geneTrack, vContextFilters, log);
		}
		VariantContextFilter[] vContextFilters = new VariantContextFilter[2];
		vContextFilters[0] = getQualityFilter(log);
		vContextFilters[1] = getJEXLAt(SNPEFF_NAMES.length - 1, log);
		trackersCase[SNPEFF_NAMES.length - 1] = new TallyTracker("CASE_" + SNPEFF_NAMES[SNPEFF_NAMES.length - 1], geneTrack, vContextFilters, log);
		trackersControl[SNPEFF_NAMES.length - 1] = new TallyTracker("CONTROL_" + SNPEFF_NAMES[SNPEFF_NAMES.length - 1], geneTrack, vContextFilters, log);
		return totalQuality;
	}

	public void tallyCharge(String fullpathToChargeVCF) {
		this.trackersCharge = new TallyTracker[SNPEFF_NAMES.length];

		for (int i = 0; i < SNPEFF_IMPACTS.length; i++) {
			VariantContextFilter[] vContextFilters = new VariantContextFilter[1];
			vContextFilters[0] = getJEXLAt(i, log);
			trackersCharge[i] = new TallyTracker("CHARGE_" + SNPEFF_NAMES[i], geneTrack, vContextFilters, log);
			trackersCharge[i].setCharge(true);
		}
		VariantContextFilter[] vContextFilters = new VariantContextFilter[1];
		vContextFilters[0] = getJEXLAt(SNPEFF_NAMES.length - 1, log);
		trackersCharge[SNPEFF_NAMES.length - 1] = new TallyTracker("CHARGE_" + SNPEFF_NAMES[SNPEFF_NAMES.length - 1], geneTrack, vContextFilters, log);
		trackersCharge[SNPEFF_NAMES.length - 1].setCharge(true);
		int count = 0;
		VCFFileReader reader = new VCFFileReader(fullpathToChargeVCF, false);
		for (VariantContext vc : reader) {
			for (int i = 0; i < trackersCharge.length; i++) {
				trackersCharge[i].addIfPasses(vc);
			}
			count++;
			if (count % 100000 == 0) {
				log.reportTimeInfo("Scanned " + count + " total variants from " + ext.removeDirectoryInfo(fullpathToChargeVCF));
			}
		}
		reader.close();

	}

	public void summarize(String fullPathToOutput) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.print("GENE\tCHR\tSTART\tSTOP\tUCSC\tTotal_length\tMrna_length");
			for (int i = 0; i < trackersCase.length; i++) {
				writer.print("\t" + trackersCase[i].getTallyName() + "_NUM_VAR" + "\t" + trackersControl[i].getTallyName() + "_NUM_VAR" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W");
			}
			for (int i = 0; i < trackersCase.length; i++) {
				writer.print("\t" + trackersCase[i].getTallyName() + "_MAC" + "\t" + trackersControl[i].getTallyName() + "_MAC" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W");
			}
			writer.println();
			Set<String> genes = trackersCase[0].getTally().keySet();
			for (String gene : genes) {
				String chr = Positions.getChromosomeUCSC(trackersCase[0].getGene(gene)[0].getChr(), true);
				String start = trackersCase[0].getGene(gene)[0].getStart() + "";
				String stop = trackersCase[0].getGene(gene)[0].getStop() + "";

				writer.print(gene + "\t" + chr + "\t" + start + "\t" + stop + "\t" + trackersCase[0].getGene(gene)[0].getUCSCLink("hg19") + "\t" + trackersCase[0].getGeneTotalLength(gene) + "\t" + trackersCase[0].getGeneTotalMrnaLength(gene));
				for (int i = 0; i < trackersCase.length; i++) {
					writer.print("\t" + trackersCase[i].getTally().get(gene) + "\t" + trackersControl[i].getTally().get(gene) + "\t" + trackersCharge[i].getTallyMAC().get(gene));
				}
				for (int i = 0; i < trackersCase.length; i++) {
					writer.print("\t" + trackersCase[i].getTallyMAC().get(gene) + "\t" + trackersControl[i].getTallyMAC().get(gene) + "\t" + trackersCharge[i].getTallyMAC().get(gene));
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

	private static class TallyTracker {
		private Hashtable<String, Float> tally;
		private Hashtable<String, Float> tallyMAC;

		private GeneTrack geneTrack;
		private String tallyName;
		private Segment[] geneSegs;
		private Hashtable<Integer, int[]> geneSegTrack;
		private VariantContextFilter[] vContextFilters;
		private ArrayList<GenePass> passingLocs;
		private Logger log;
		private boolean charge;
		private int count = 0;

		public TallyTracker(String tallyName, GeneTrack geneTrack, VariantContextFilter[] vContextFilters, Logger log) {
			super();
			this.log = log;

			this.tallyName = tallyName;
			this.tally = getTally(geneTrack, log);
			this.tallyMAC = getTally(geneTrack, log);
			this.vContextFilters = vContextFilters;
			this.geneTrack = geneTrack;
			this.charge = false;
			this.passingLocs = new ArrayList<VCFTally.GenePass>();
			populateQuery();
		}

		public void setCharge(boolean charge) {
			this.charge = charge;
		}

		private void populateQuery() {
			this.geneSegTrack = new Hashtable<Integer, int[]>();
			this.geneSegs = new Segment[total(geneTrack)];
			GeneData[][] gDatas = geneTrack.getGenes();
			int index = 0;
			int curChr = 0;
			int curPos = 0;
			int total = 0;
			for (int chr = 0; chr < gDatas.length; chr++) {
				curPos = 0;
				total += gDatas[chr].length;
				for (int chrIndex = 0; chrIndex < gDatas[chr].length; chrIndex++) {
					if (gDatas[chr][chrIndex].getChr() < curChr || gDatas[chr][chrIndex].getStart() < curPos) {
						log.reportTimeError("Unsorted Genetrack!");
						return;
					} else {
						curChr = gDatas[chr][chrIndex].getChr();
						curPos = gDatas[chr][chrIndex].getStart();
					}
					geneSegs[index] = new Segment(gDatas[chr][chrIndex].getChr(), gDatas[chr][chrIndex].getStart(), gDatas[chr][chrIndex].getStop());
					geneSegTrack.put(index, new int[] { chr, chrIndex });
					index++;
				}
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

		private static Hashtable<String, Float> getTally(GeneTrack geneTrack, Logger log) {
			Hashtable<String, Float> tally = new Hashtable<String, Float>();
			GeneData[][] gDatas = geneTrack.getGenes();
			for (int i = 0; i < gDatas.length; i++) {
				for (int j = 0; j < gDatas[i].length; j++) {
					tally.put(gDatas[i][j].getGeneName(), (float) 0.0);
				}
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

		public int addIfPasses(VariantContext vc) {
			int numAdded = 0;
			for (int i = 0; i < vContextFilters.length; i++) {
				if (!vContextFilters[i].filter(vc).passed()) {

					return numAdded;
				}
			}
			double mac = Double.NaN;
			if (!charge) {
				VCOps.getMAC(vc, null);
			} else {
				mac = vc.getCommonInfo().getAttributeAsDouble("MAF_whites", 0.0);
			}

			if (mac > 0) {
				int[] matches = Segment.binarySearchForAllOverLappingIndices(new Segment(Positions.chromosomeNumber(vc.getChr()), vc.getStart(), vc.getEnd()), geneSegs);
				if (matches != null && matches.length > 0) {
					GeneData[] geneDatas = new GeneData[matches.length];
					for (int i = 0; i < geneDatas.length; i++) {
						int[] indices = geneSegTrack.get(matches[i]);
						geneDatas[i] = geneTrack.getGenes()[indices[0]][indices[1]];
					}
					if (geneDatas.length > 0) {

						for (int i = 0; i < geneDatas.length; i++) {
							if (!charge) {
								passingLocs.add(new GenePass(VCOps.getSegment(vc), geneDatas[i].getGeneName()));
								addTally(geneDatas[i].getGeneName(), 1, (float) mac);
							} else {
								addTally(geneDatas[i].getGeneName(), 1, (float) mac);
							}
							numAdded++;
						}
					}
				}
			}
			return numAdded;
		}

		public String[] getBed() {
			String[] bed = new String[passingLocs.size()];
			for (int i = 0; i < bed.length; i++) {
				GenePass pass = passingLocs.get(i);
				bed[i] = pass.getVarSeg().getChr() + "\t" + pass.getVarSeg().getStart() + "\t" + pass.getVarSeg().getStop() + "\t" + pass.getGeneName() + ":" + tallyName;
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

		public GeneData[] getGene(String gene) {
			return geneTrack.lookupAllGeneDatas(gene);
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
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, dp, gq };
		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] {}, null, null, log);
		return vContextFilter;
	}

	private static class GenePass {
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
		String vpopFile = "D:/data/Project_Tsai_Project_021/JointAnalysis/vPopCaseControl.txt";
		String fullpathToChargeVCF = "D:/data/CHARGE/CHARGE_MAFS/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.sed.vcf";
		VcfPopulation vpop = VcfPopulation.load(vpopFile, new Logger());
		String geneTrackFile = "N:/statgen/NCBI/RefSeq_hg19.gtrack";
		VCFTally tally = new VCFTally(vcf, GeneTrack.load(geneTrackFile, false), vpop, new Logger());
		tally.tallyCaseControlVCF();
		tally.tallyCharge(fullpathToChargeVCF);
		tally.summarize(ext.parseDirectoryOfFile(vcf) + "tallyCounts.txt");
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "VCFTally.dat";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.analysis.VCFTally requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			log = new Logger(logfile);
			test();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
