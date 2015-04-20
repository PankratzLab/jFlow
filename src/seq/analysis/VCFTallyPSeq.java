package seq.analysis;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Set;

import seq.analysis.PlinkSeq.ANALYSIS_TYPES;
import seq.analysis.PlinkSeq.BURDEN_Tests;
import seq.analysis.PlinkSeq.PlinkSeqWorker;
import seq.analysis.PlinkSeqUtils.PlinkSeqBurdenSummary;
import seq.analysis.PlinkSeqUtils.PlinkSeqTestSummary;
import seq.analysis.PlinkSeqUtils.PseqProject;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.pathway.GenomeRegions;
import seq.pathway.Pathways;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.GeneTrack;

/**
 * Does a vcf tally and incorporates plink seq results
 *
 */
public class VCFTallyPSeq extends VCFTally {

	private PlinkSeq plinkSeq;
	private PseqProject pseqProject;
	private String varList;
	private String locFile;
	private PlinkSeqBurdenSummary[] summaries;

	public VCFTallyPSeq(String vcf, GenomeRegions geneTrack, VcfPopulation vpop, CASE_CONTROL_TYPE type, String plinkSeqResourceDirectory, String plinkSeqProjName, Logger log) {
		super(vcf, geneTrack, vpop, type, log);
		VCFOps.verifyIndexRegular(vcf, log);
		this.locFile = plinkSeqResourceDirectory + ext.rootOf(geneTrack.getGeneTrack().getGeneSetFilename()) + ".reg";
		PlinkSeqUtils.generatePlinkSeqLoc(geneTrack, locFile, log);
		// VCFOps.gzipAndIndex(vcf, log);
		this.plinkSeq = new PlinkSeq(false, true, log);
		this.pseqProject = PlinkSeq.initialize(plinkSeq, plinkSeqProjName, vcf, vpop, plinkSeqResourceDirectory, log);

	}

	private void tallyCaseControlVCF(int altAlleleDepth) {
		this.varList = pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName(), true) + ".varList";
		// if (!Files.exists(varList)) {
		tallyCaseControlVCF(altAlleleDepth, varList);
		plinkSeq.eraseAndLoadVarSet(pseqProject, varList);
		//
		// } else {
		// log.reportFileExists(varList);
		// log.reportTimeWarning("Assuming Variant sets have not changed");
		// }
	}

	public void fullGamutAssoc(int numPerm, String mac, int altAlleleDepth, String fullpathToChargeVCF, int numThreads) {
		tallyCaseControlVCF(altAlleleDepth);
		tallyCharge(fullpathToChargeVCF);

		String fullPathToOutput = pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName()) + "_" + type + ".summary";

		String[] varMasks = Array.unique(HashVec.loadFileToStringArray(varList, false, new int[] { 1 }, true));
		summarize(fullPathToOutput);
		String locFile = pseqProject.getResourceDirectory() + ext.rootOf(genomeRegions.getGeneTrack().getGeneSetFilename() + ".reg");
		PlinkSeqUtils.generatePlinkSeqLoc(genomeRegions, locFile, log);

		PlinkSeqWorker[] complete = plinkSeq.fullGamutAssoc(pseqProject, new String[] { ext.rootOf(locFile) }, varMasks, numPerm, mac, ext.rootOf(vpop.getFileName()), numThreads);
		this.summaries = new PlinkSeqBurdenSummary[SNPEFF_NAMES.length];
		int index = 0;
		for (int i = 0; i < complete.length; i++) {
			ANALYSIS_TYPES type = complete[i].getType();
			switch (type) {
			case BURDEN:
				String analysis = ext.rootOf(complete[i].getOutputFiles()[0]);
				analysis = analysis.replaceAll(".*" + ext.rootOf(locFile) + ".", "");
				analysis = analysis.replaceAll("_" + this.type + ".*", "");
				PlinkSeqBurdenSummary plinkSeqBurdenSummary = new PlinkSeqBurdenSummary(analysis, complete[i].getOutputFiles()[0], log);
				plinkSeqBurdenSummary.load();
				summaries[index] = plinkSeqBurdenSummary;
				index++;
				break;
			case I_SUMMARY:
				break;
			case V_ASSOC:
				break;
			case V_SUMMARY:
				break;
			default:
				log.reportTimeError("INVALID analysis type " + type);
				break;

			}
		}
	}

	public PseqProject getPseqProject() {
		return pseqProject;
	}

	public void summarize() {
		String fullPathToOutput = pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName()) + "_" + type + ".hit.summary";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.print("GENE\tCHR\tSTART\tSTOP\tUCSC\tTotal_length\tMrna_length");
			for (int i = 0; i < trackersCase.length; i++) {
				writer.print("\t" + trackersCase[i].getTallyName() + "_NUM_VAR" + "\t" + trackersCase[i].getTallyName() + "_NUM_INDS" + "\t" + trackersControl[i].getTallyName() + "_NUM_VAR" + "\t" + trackersControl[i].getTallyName() + "_NUM_INDS" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B_CASE_Expect");

			}
			for (int i = 0; i < trackersCase.length; i++) {
				writer.print("\t" + trackersCase[i].getTallyName() + "_ALT_ALLELE_COUNT" + "\t" + trackersCase[i].getTallyName() + "_NUM_INDS" + "\t" + trackersControl[i].getTallyName() + "_ALT_ALLELE_COUNT" + "\t" + trackersControl[i].getTallyName() + "_NUM_INDS" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B" + "\t" + trackersCharge[i].getTallyName() + "_CMAF_W_B_CASE_Expect");
				for (int j = 0; j < BURDEN_Tests.values().length; j++) {
					BURDEN_Tests test = BURDEN_Tests.values()[j];
					writer.print("\t" + trackersCase[i].getTallyName() + "_" + test + "_P" + "\t" + trackersCase[i].getTallyName() + "_" + test + "_I" + "\t" + trackersCase[i].getTallyName() + "_" + test + "_DESC");
				}
			}
			writer.println();
			writer.flush();
			Set<String> sets = trackersCase[0].getTally().keySet();
			for (String set : sets) {
				if (trackersCase[0].getGene(set) != null&&trackersCase[0].getGene(set).length>0) {
					String chr = Positions.getChromosomeUCSC(trackersCase[0].getGene(set)[0].getChr(), true);
					String start = trackersCase[0].getGene(set)[0].getStart() + "";
					String stop = trackersCase[0].getGene(set)[0].getStop() + "";

					writer.print(set + "\t" + chr + "\t" + start + "\t" + stop + "\t" + trackersCase[0].getGene(set)[0].getUCSCLink("hg19") + "\t" + trackersCase[0].getGeneTotalLength(set) + "\t" + trackersCase[0].getGeneTotalMrnaLength(set));
				} else {
					writer.print(set + "\tNA\tNA\tNA\tNA\tNA\tNA");

				}

				for (int i = 0; i < trackersCase.length; i++) {
					int numCases = trackersCase[i].getUniqs().get(set).size();
					int numControls = trackersControl[i].getUniqs().get(set).size();
					writer.print("\t" + trackersCase[i].getTally().get(set) + "\t" + numCases + "\t" + trackersControl[i].getTally().get(set) + "\t" + numControls + "\t" + trackersCharge[i].getTallyMAC().get(set) + "\t" + ((double) vpop.getSubPop().get("CASE").size() * trackersCharge[i].getTallyMAC().get(set)));

				}
				for (int i = 0; i < trackersCase.length; i++) {
					int numCases = trackersCase[i].getUniqs().get(set).size();
					int numControls = trackersControl[i].getUniqs().get(set).size();
					writer.print("\t" + trackersCase[i].getTallyMAC().get(set) + "\t" + numCases + "\t" + trackersControl[i].getTallyMAC().get(set) + "\t" + numControls + "\t" + trackersCharge[i].getTallyMAC().get(set) + "\t" + ((double) vpop.getSubPop().get("CASE").size() * trackersCharge[i].getTallyMAC().get(set)));
					for (int j = 0; j < summaries.length; j++) {
						// System.out.println(trackersCase[1].getTallyName());
						// System.out.println(i+"\t"+j+"\t"+trackersCase[i].getTallyName()+"\t"+trackersCase.length+"\t"+summaries[j].getAnalysis());
						if (trackersCase[i].getTallyName().endsWith(summaries[j].getAnalysis())) {
							for (int j2 = 0; j2 < BURDEN_Tests.values().length; j2++) {
								BURDEN_Tests test = BURDEN_Tests.values()[j2];
								if (summaries[j].hasSummaryFor(set) && summaries[j].getPlinkSeqLocSummaryFor(set).getSummaries()[j2] != null) {
									PlinkSeqTestSummary pstSummary = summaries[j].getPlinkSeqLocSummaryFor(set).getSummaries()[j2];
									if (pstSummary.getType() != test) {
										System.out.println(pstSummary.getType() + "\t" + test);
										log.reportTimeError("Mismatched parsing error, halting...");
									}
									writer.print("\t" + pstSummary.getP() + "\t" + pstSummary.getI() + "\t" + "\"" + pstSummary.getDesc().replaceAll("/", "-") + "\"");
								} else {
									writer.print("\t1\t1\tNO_TEST");
								}
							}
						}
					}
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

	public static void test() {
		String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
		String[] vpopFiles = new String[] { "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqProj_tsai_spector_joint_AgilentCaptureRecal/vPopCaseControl.txt" };
		String fullpathToChargeVCF = "/panfs/roc/groups/14/tsaim/shared/bin/CHARGE/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";
		String resourceDirectory = "/home/tsaim/public/bin/pseqRef/hg19/";
		Logger log = new Logger(ext.rootOf(vcf, false) + "tally.log");
		String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
		String keggPathwayFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/kegg.ser";
		int altAlleleDepth = 6;

		VcfPopulation vpop = VcfPopulation.load(vpopFiles[0], POPULATION_TYPE.CASE_CONTROL, log);
		vpop.report();
		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);
		geneTrack.setGeneSetFilename(geneTrackFile);
		Pathways pathways = Pathways.load(keggPathwayFile);
		GenomeRegions gRegions = new GenomeRegions(geneTrack, pathways, log);
		VCFTallyPSeq vcfTallyPSeq = new VCFTallyPSeq(vcf, gRegions, vpop, CASE_CONTROL_TYPE.BOTH_PASS, resourceDirectory, null, log);
		vcfTallyPSeq.fullGamutAssoc(-1, "0", altAlleleDepth, fullpathToChargeVCF, 4);
		vcfTallyPSeq.summarize();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "VCFTallyPSeq.dat";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.analysis.VCFTallyPSeq requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

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

//
// for (int i = 0; i < vpopFiles.length; i++) {
// String outputList = ext.rootOf(vpopFiles[i], false) + ".plinkSeqVar";
// if (Files.exists(outputList)) {
// new File(outputList).delete();
// }
// for (int j = 0; j < VCFTally.CASE_CONTROL_TYPE.values().length; j++) {
// vpop.report();
// VCFTally tally = new VCFTally(vcf, GeneTrack.load(geneTrackFile, false), vpop, CASE_CONTROL_TYPE.values()[j], log);
// tally.tallyCaseControlVCF(altAlleleDepth, outputList);
// tally.tallyCharge(fullpathToChargeVCF);
// tally.summarize(ext.parseDirectoryOfFile(vcf) + ext.rootOf(vpopFiles[i]) + "tallyCounts." + CASE_CONTROL_TYPE.values()[j] + ".txt");
// }
// }
