package org.genvisis.seq.analysis;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Hashtable;
import java.util.Set;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.DumpMultiLoc;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.seq.analysis.PlinkSeq.ANALYSIS_TYPES;
import org.genvisis.seq.analysis.PlinkSeq.BURDEN_Tests;
import org.genvisis.seq.analysis.PlinkSeq.PlinkSeqWorker;
import org.genvisis.seq.analysis.PlinkSeqUtils.PlinkSeqBurdenSummary;
import org.genvisis.seq.analysis.PlinkSeqUtils.PlinkSeqTestSummary;
import org.genvisis.seq.analysis.PlinkSeqUtils.PseqProject;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.pathway.GenomeRegions;
import org.genvisis.seq.pathway.Pathways;
import org.genvisis.seq.qc.FilterNGS;

/**
 * Does a vcf tally and incorporates plink seq results
 *
 */
public class VCFTallyPSeq extends VCFTally implements Serializable {

	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	private final PlinkSeq plinkSeq;
	private final PseqProject pseqProject;
	private final String varList;
	private final String locFile;
	private PlinkSeqBurdenSummary[] summaries;

	public VCFTallyPSeq(String vcf, GenomeRegions geneTrack, VcfPopulation vpop,
											CASE_CONTROL_TYPE type, String plinkSeqResourceDirectory,
											String plinkSeqProjName, Logger log) {
		super(vcf, geneTrack, vpop, type, log);
		VCFOps.verifyIndex(vcf, log);
		locFile = plinkSeqResourceDirectory + ext.rootOf(geneTrack.getGeneTrack().getGeneSetFilename())
							+ "_Gen.reg";
		PlinkSeqUtils.generatePlinkSeqLoc(geneTrack, locFile, log);
		plinkSeq = new PlinkSeq(false, true, log);
		pseqProject = PlinkSeq.initialize(plinkSeq, plinkSeqProjName,
																			ext.parseDirectoryOfFile(vpop.getFileName()), vcf, vpop,
																			plinkSeqResourceDirectory, false, false, log);
		varList = pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName(), true) + ".varList";
	}

	public VCFTallyPSeq(String vcf, GenomeRegions genomeRegions, VcfPopulation vpop,
											CASE_CONTROL_TYPE type, Logger log, PlinkSeq plinkSeq,
											PseqProject pseqProject, String varList, String locFile,
											PlinkSeqBurdenSummary[] summaries) {
		super(vcf, genomeRegions, vpop, type, log);
		this.plinkSeq = plinkSeq;
		this.pseqProject = pseqProject;
		this.varList = varList;
		this.locFile = locFile;
		this.summaries = summaries;
	}

	public VCFTallyPSeq(String vcf, GenomeRegions genomeRegions, VcfPopulation vpop,
											TallyTracker[] trackersCase, TallyTracker[] trackersControl,
											TallyTracker[] trackersCharge, Logger log, CASE_CONTROL_TYPE type,
											PlinkSeq plinkSeq, PseqProject pseqProject, String varList, String locFile,
											PlinkSeqBurdenSummary[] summaries) {
		super(vcf, genomeRegions, vpop, trackersCase, trackersControl, trackersCharge, log, type);
		this.plinkSeq = plinkSeq;
		this.pseqProject = pseqProject;
		this.varList = varList;
		this.locFile = locFile;
		this.summaries = summaries;
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static VCFTallyPSeq load(String filename) {
		return (VCFTallyPSeq) SerializedFiles.readSerial(filename, false);
	}

	public String getOutput() {
		return pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName()) + "_" + type
					 + ".summary";
	}

	public String getSerFile() {
		return getOutput() + ".ser";
	}

	private void tally(FilterNGS filterNGS, String fullpathToChargeVCF) {
		String serFile = getSerFile();
		if (Files.exists(serFile) && load(serFile) != null) {
			log.reportTimeWarning(serFile + " loading pre-tallied values instead");
		} else {
			tallyCaseControlVCF(filterNGS, varList);
			tallyCharge(fullpathToChargeVCF);
			// summarize(fullPathToOutput);
			serialize(serFile);
			plinkSeq.eraseAndLoadVarSet(pseqProject, varList);
		}
	}

	public void fullGamutAssoc(int numPerm, String mac, FilterNGS filterNGS,
														 String fullpathToChargeVCF, int numThreads) {
		tally(filterNGS, fullpathToChargeVCF);
		String[] varMasks = ArrayUtils.unique(HashVec.loadFileToStringArray(varList, false,
																																				new int[] {1}, true));
		String locFile = pseqProject.getResourceDirectory()
										 + ext.rootOf(genomeRegions.getGeneTrack().getGeneSetFilename() + ".reg");

		PlinkSeqWorker[] complete = plinkSeq.fullGamutAssoc(pseqProject,
																												new String[] {ext.rootOf(locFile)},
																												varMasks, numPerm, mac,
																												ext.rootOf(vpop.getFileName()), numThreads);
		summaries = new PlinkSeqBurdenSummary[SNPEFF_NAMES.length];
		int index = 0;
		for (PlinkSeqWorker element : complete) {
			ANALYSIS_TYPES type = element.getType();
			switch (type) {
				case BURDEN:
					String analysis = ext.rootOf(element.getOutputFiles()[0]);
					analysis = analysis.replaceAll(".*" + ext.rootOf(locFile) + ".", "");
					analysis = analysis.replaceAll("_" + this.type + ".*", "");
					PlinkSeqBurdenSummary plinkSeqBurdenSummary = new PlinkSeqBurdenSummary(analysis,
																																									element.getOutputFiles()[0],
																																									log);
					plinkSeqBurdenSummary.load();
					plinkSeqBurdenSummary.correctPvalues();
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
					log.reportError("INVALID analysis type " + type);
					break;
			}
		}
	}

	public PseqProject getPseqProject() {
		return pseqProject;
	}

	public void summarize() {
		String fullPathToOutput = pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName())
															+ "_" + type + ".hit.summary5_snpEffGeneProper";
		try {
			if (!Files.exists(fullPathToOutput)) {
				PrintWriter writer = Files.openAppropriateWriter(fullPathToOutput);
				writer.print("GENE\tCHR\tSTART\tSTOP\tUCSC\tTotal_length\tMrna_length\tmultiLoc");
				for (int i = 0; i < trackersCase.length; i++) {
					writer.print("\t" + trackersCase[i].getTallyName() + "_NUM_VAR" + "\t"
											 + trackersCase[i].getTallyName() + "_NUM_UNIQ_INDS" + "\t"
											 + trackersCase[i].getTallyName() + "_TotalWithAlt" + "\t"
											 + trackersControl[i].getTallyName() + "_NUM_VAR" + "\t"
											 + trackersControl[i].getTallyName() + "_NUM_UNIQ_INDS" + "\t"
											 + trackersControl[i].getTallyName() + "_TotalWithAlt" + "\t"
											 + trackersCharge[i].getTallyName() + "_CMAF_W_B" + "\t"
											 + trackersCharge[i].getTallyName() + "_CMAF_W_B_CASE_Expect");
				}
				for (int i = 0; i < trackersCase.length; i++) {
					writer.print("\t" + trackersCase[i].getTallyName() + "_ALT_ALLELE_COUNT" + "\t"
											 + trackersCase[i].getTallyName() + "_NUM_UNIQ_INDS" + "\t"
											 + trackersCase[i].getTallyName() + "_TotalWithAlt" + "\t"
											 + trackersControl[i].getTallyName() + "_ALT_ALLELE_COUNT" + "\t"
											 + trackersControl[i].getTallyName() + "_NUM_UNIQ_INDS" + "\t"
											 + trackersControl[i].getTallyName() + "_TotalWithAlt" + "\t"
											 + trackersCharge[i].getTallyName() + "_CMAF_W_B" + "\t"
											 + trackersCharge[i].getTallyName() + "_CMAF_W_B_CASE_Expect");
					for (int j = 0; j < BURDEN_Tests.values().length; j++) {
						BURDEN_Tests test = BURDEN_Tests.values()[j];
						writer.print("\t" + trackersCase[i].getTallyName() + "_" + test + "_P" + "\t"
												 + trackersCase[i].getTallyName() + "_" + test + "_I" + "\t"
												 + trackersCase[i].getTallyName() + "_" + test + "_DESC"
												 + (test == BURDEN_Tests.BURDEN ? "\t" + trackersCase[i].getTallyName()
																													+ "_" + test + "_NUMCASE" + "\t"
																													+ trackersCase[i].getTallyName() + "_"
																													+ test + "_NUMCONTROLS"
																												: ""));
						writer.print("\t" + trackersCase[i].getTallyName() + "_" + test + "_BONF_FULL_P");
						for (double element : PlinkSeqUtils.PlinkSeqBurdenSummary.I_THRESHOLDS) {
							writer.print("\t" + trackersCase[i].getTallyName() + "_" + test + "_BONF_I_" + element
													 + "_P");
						}
					}
				}
				writer.println();
				Set<String> sets = trackersCase[0].getTally().keySet();

				for (String set : sets) {
					long time = System.currentTimeMillis();
					log.reportTimeInfo("1" + ext.getTimeElapsed(time));

					String curString = "";
					if (trackersCase[0].getGene(set) != null && trackersCase[0].getGene(set).length > 0) {
						String chr = Positions.getChromosomeUCSC(trackersCase[0].getGene(set)[0].getChr(),
																										 true);
						String start = trackersCase[0].getGene(set)[0].getStart() + "";
						String stop = trackersCase[0].getGene(set)[0].getStop() + "";
						curString += set + "\t" + chr + "\t" + start + "\t" + stop + "\t"
												 + trackersCase[0].getGene(set)[0].getUCSCLink("hg19") + "\t"
												 + trackersCase[0].getGeneTotalLength(set) + "\t"
												 + trackersCase[0].getGeneTotalMrnaLength(set) + "\t"
												 + (trackersCase[0].getGene(set)[0].getMultiLoc() > 0);
					} else {
						curString += set + "\tNA\tNA\tNA\tNA\tNA\tNA\tTRUE";
					}
					for (int i = 0; i < trackersCase.length; i++) {
						int numCases = trackersCase[i].getUniqs().get(set).size();
						int numControls = trackersControl[i].getUniqs().get(set).size();
						int totalCasesWithAlt = trackersCase[i].getAll().get(set).size();
						int totalControlsWithAlt = trackersControl[i].getAll().get(set).size();
						curString += "\t" + trackersCase[i].getTally().get(set) + "\t" + numCases + "\t"
												 + totalCasesWithAlt + "\t" + trackersControl[i].getTally().get(set) + "\t"
												 + numControls + "\t" + totalControlsWithAlt + "\t"
												 + trackersCharge[i].getTallyMAC().get(set) + "\t"
												 + ((double) vpop.getSubPop().get("CASE").size()
														* trackersCharge[i].getTallyMAC().get(set));
					}

					for (int i = 0; i < trackersCase.length; i++) {
						log.reportTimeInfo("3" + ext.getTimeElapsed(time));

						int numCases = trackersCase[i].getUniqs().get(set).size();
						int numControls = trackersControl[i].getUniqs().get(set).size();
						int totalCasesWithAlt = trackersCase[i].getAll().get(set).size();
						int totalControlsWithAlt = trackersControl[i].getAll().get(set).size();
						log.reportTimeInfo("4" + ext.getTimeElapsed(time));

						Hashtable<String, Float> caseMac = trackersCase[i].getTallyMAC();
						Hashtable<String, Float> controlMac = trackersControl[i].getTallyMAC();
						Hashtable<String, Float> chargeMac = trackersCharge[i].getTallyMAC();

						curString += "\t" + caseMac.get(set) + "\t" + numCases + "\t" + totalCasesWithAlt + "\t"
												 + controlMac.get(set) + "\t" + numControls + "\t" + totalControlsWithAlt
												 + "\t" + chargeMac.get(set) + "\t"
												 + ((double) vpop.getSubPop().get("CASE").size() * chargeMac.get(set));
						for (PlinkSeqBurdenSummary curSummary : summaries) {
							if (trackersCase[i].getTallyName().endsWith(curSummary.getAnalysis())) {
								for (int j2 = 0; j2 < BURDEN_Tests.values().length; j2++) {
									BURDEN_Tests test = BURDEN_Tests.values()[j2];
									String key = set;
									boolean hasGenvisisTest = false;
									if (curSummary.hasSummaryFor(key + PlinkSeqUtils.GENVISIS_GENE)) {
										key = set + PlinkSeqUtils.GENVISIS_GENE;
										hasGenvisisTest = true;
									} else if (curSummary.hasSummaryFor(key + PlinkSeqUtils.GENVISIS_PATHWAY)) {
										key = set + PlinkSeqUtils.GENVISIS_PATHWAY;
										hasGenvisisTest = true;
									} else if (curSummary.hasSummaryFor(key + PlinkSeqUtils.GENVISIS)) {
										key = set + PlinkSeqUtils.GENVISIS;
										hasGenvisisTest = true;
									}
									if (hasGenvisisTest && curSummary.hasSummaryFor(key)
											&& curSummary.getPlinkSeqLocSummaryFor(key).getSummaries()[j2] != null) {
										PlinkSeqTestSummary pstSummary = curSummary.getPlinkSeqLocSummaryFor(key)
																															 .getSummaries()[j2];
										if (pstSummary.getType() != test) {
											System.out.println(pstSummary.getType() + "\t" + test);
											log.reportError("Mismatched parsing error, halting...");
										}
										String desc[] = pstSummary.getDesc().split("/");
										if (desc.length != 2 && test == BURDEN_Tests.BURDEN) {
											log.reportError("Did not find two counts for burden test");
										}
										curString += "\t" + pstSummary.getP() + "\t" + pstSummary.getI() + "\t" + ""
																 + pstSummary.getDesc().replaceAll("/", "::").replaceAll("/", "::")
																 + (test == BURDEN_Tests.BURDEN ? "\t" + ArrayUtils.toStr(desc)
																																: "");
										curString += "\t" + pstSummary.getBonfFull();
										for (int k = 0; k < PlinkSeqUtils.PlinkSeqBurdenSummary.I_THRESHOLDS.length; k++) {
											curString += "\t" + pstSummary.getBonfsI()[k];
										}
									} else {
										curString += "\t1\t1\tNO_TEST" + (test == BURDEN_Tests.BURDEN ? "\t0\t0" : "");
										curString += "\t1";
										for (double element : PlinkSeqUtils.PlinkSeqBurdenSummary.I_THRESHOLDS) {
											curString += "\t1";
										}
									}
								}
							}
						}
					}
					log.reportTimeInfo(ext.getTimeElapsed(time));
					writer.println(curString);
				}
				writer.close();
			} else {
				log.reportFileExists(fullPathToOutput);
			}
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathToOutput);
			log.reportException(e);
		}

		String bed = fullPathToOutput + ".CASE.bed";
		try {
			PrintWriter writer = Files.openAppropriateWriter(bed);
			for (TallyTracker element : trackersCase) {
				String[] tbed = element.getBed();
				for (String element2 : tbed) {
					writer.println(element2);
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + bed);
			log.reportException(e);
		}

	}

	// private void loadFromFile(String fullPathToOutput) {
	// initializeTrackers();
	// try {
	// BufferedReader reader = Files.getAppropriateReader(fullPathToOutput);
	// String[] header = reader.readLine().trim().split("\t");
	//
	// while (reader.ready()) {
	//
	// }
	// reader.close();
	// } catch (FileNotFoundException fnfe) {
	// log.reportError("Error: file \"" + fullPathToOutput + "\" not found in current directory");
	// return;
	// } catch (IOException ioe) {
	// log.reportError("Error reading file \"" + fullPathToOutput + "\"");
	// return;
	// }
	// }

	public static void test(String vcf, String[] vpopFiles) {

		for (String vpopFile : vpopFiles) {
			String fullpathToChargeVCF = "/panfs/roc/groups/14/tsaim/shared/bin/CHARGE/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";
			String resourceDirectory = "/home/tsaim/public/bin/pseqRef/hg19/";
			Logger log = new Logger(ext.rootOf(vcf, false) + "tally.log");
			String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
			String keggPathwayFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/kegg.ser";
			int altAlleleDepth = -1;
			FilterNGS filterNGS = new FilterNGS();
			filterNGS.setAltAlleleDepthFilter(new int[] {altAlleleDepth});
			VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.CASE_CONTROL, log);
			vpop.report();
			GeneTrack geneTrack = GeneTrack.load(geneTrackFile);

			DumpMultiLoc.dumpMultiLoc(geneTrackFile, ext.rootOf(vpopFile) + ".multiLoc", log);
			// GeneData[][] genes = geneTrack.getGenes();
			// for (int i = 0; i < genes.length; i++) {
			// for (int j = 0; j < genes[i].length; j++) {
			// if (genes[i][j].getGeneName().equals("ULK4P3")) {
			// System.out.println(genes[i][j].getChr() + "\t" + genes[i][j].getStart() + "\t" +
			// genes[i][j].getStop() + "\t" + genes[i][j].getMultiLoc());
			// }
			// }
			// }
			geneTrack.setGeneSetFilename(geneTrackFile);
			Pathways pathways = Pathways.load(keggPathwayFile);
			GenomeRegions gRegions = new GenomeRegions(geneTrack, pathways, log);
			VCFTallyPSeq vcfTallyPSeq = new VCFTallyPSeq(vcf, gRegions, vpop, CASE_CONTROL_TYPE.BOTH_PASS,
																									 resourceDirectory,
																									 ext.rootOf(vpop.getFileName()), log);
			if (Files.exists(vcfTallyPSeq.getSerFile())) {
				vcfTallyPSeq = load(vcfTallyPSeq.getSerFile());
			}
			vcfTallyPSeq.fullGamutAssoc(-1, "0", filterNGS, fullpathToChargeVCF, 24);
			vcfTallyPSeq.summarize();
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String vcf =
		// "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
		// String vcf =
		// "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqTallyTest/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.errorRegions2.vcf";
		// String vcf =
		// "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqTallyTest/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vPopCaseControl.vcf";
		// String vcf =
		// "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.CUSHINGS.vcf.gz";
		String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/vcf/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
		String[] vpopFiles = new String[] {"/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqProj_tsai_spector_joint_AgilentCaptureRecal/vPopCaseControl.txt",
																			 "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqProj_tsai_spector_joint_AgilentCaptureRecal/vPopCaseControlAllRaces.txt"};

		// String logfile = null;
		// Logger log;
		String usage = "\n" + "seq.analysis.VCFTallyPSeq requires 0-1 arguments\n";
		usage += "   (1) vcf file (i.e. vcf=" + vcf + " (default))\n" + "";
		usage += "   (2) vpop files, comma delimited (i.e. vpop=" + ArrayUtils.toStr(vpopFiles)
						 + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("vcf=")) {
				vcf = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("vpop=")) {
				vpopFiles = arg.split("=")[1].split(",");
				numArgs--;
			}
			// else if (args[i].startsWith("log=")) {
			// logfile = args[i].split("=")[1];
			// numArgs--;
			// }
			else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			// log = new Logger(logfile);
			test(vcf, vpopFiles);
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
// VCFTally tally = new VCFTally(vcf, GeneTrack.load(geneTrackFile, false), vpop,
// CASE_CONTROL_TYPE.values()[j], log);
// tally.tallyCaseControlVCF(altAlleleDepth, outputList);
// tally.tallyCharge(fullpathToChargeVCF);
// tally.summarize(ext.parseDirectoryOfFile(vcf) + ext.rootOf(vpopFiles[i]) + "tallyCounts." +
// CASE_CONTROL_TYPE.values()[j] + ".txt");
// }
// }
