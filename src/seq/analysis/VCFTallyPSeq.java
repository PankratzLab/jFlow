package seq.analysis;

import seq.analysis.PlinkSeq.PlinkSeqWorker;
import seq.analysis.PlinkSeqUtils.PseqProject;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
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

	public VCFTallyPSeq(String vcf, GeneTrack geneTrack, VcfPopulation vpop, CASE_CONTROL_TYPE type, String plinkSeqResourceDirectory, String plinkSeqProjName, Logger log) {
		super(vcf, geneTrack, vpop, type, log);
		VCFOps.verifyIndexRegular(vcf, log);
		this.locFile = plinkSeqResourceDirectory + ext.rootOf(geneTrack.getGeneSetFilename()) + ".reg";
		PlinkSeqUtils.generatePlinkSeqLoc(geneTrack, locFile, log);
		// VCFOps.gzipAndIndex(vcf, log);
		this.plinkSeq = new PlinkSeq(true, true, log);
		this.pseqProject = PlinkSeq.initialize(plinkSeq, plinkSeqProjName, vcf, vpop, plinkSeqResourceDirectory, log);

	}

	private void tallyCaseControlVCF(int altAlleleDepth) {
		this.varList = pseqProject.getProjectDirectory() + ext.rootOf(vpop.getFileName(), true) + ".varList";
		if (!Files.exists(varList)) {
			tallyCaseControlVCF(altAlleleDepth, varList);
		}
	}

	public void fullGamutAssoc(int numPerm, String mac, int altAlleleDepth, int numThreads) {
		tallyCaseControlVCF(altAlleleDepth);
		String[] varMasks = Array.unique(HashVec.loadFileToStringArray(varList, false, new int[] { 1 }, true));
		System.out.println(Array.toStr(varMasks));
		String locFile = pseqProject.getResourceDirectory() + ext.rootOf(geneTrack.getGeneSetFilename() + ".reg");
		PlinkSeqUtils.generatePlinkSeqLoc(geneTrack, locFile, log);
		PlinkSeqWorker[] complete = plinkSeq.fullGamutAssoc(pseqProject, new String[] { ext.rootOf(locFile) }, varMasks, numPerm, mac, ext.rootOf(vpop.getFileName()), numThreads);

	}

	public PseqProject getPseqProject() {
		return pseqProject;
	}

	public static void test() {
		String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
		String[] vpopFiles = new String[] { "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqProj_tsai_spector_joint_AgilentCaptureRecal/vPopCaseControl.txt" };
		String fullpathToChargeVCF = "/panfs/roc/groups/14/tsaim/shared/bin/CHARGE/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";
		String resourceDirectory = "/home/tsaim/public/bin/pseqRef/hg19/";
		Logger log = new Logger(ext.rootOf(vcf, false) + "tally.log");
		String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
		int altAlleleDepth = 6;
		VcfPopulation vpop = VcfPopulation.load(vpopFiles[0], POPULATION_TYPE.CASE_CONTROL, log);
		vpop.report();
		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);
		geneTrack.setGeneSetFilename(geneTrackFile);
		VCFTallyPSeq vcfTallyPSeq = new VCFTallyPSeq(vcf, geneTrack, vpop, CASE_CONTROL_TYPE.BOTH_PASS, resourceDirectory, null, log);
		vcfTallyPSeq.tallyCaseControlVCF(altAlleleDepth);

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
