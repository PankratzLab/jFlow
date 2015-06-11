package seq.analysis;

import htsjdk.variant.vcf.VCFFileReader;
import seq.analysis.PlinkSeq.ANALYSIS_TYPES;
import seq.analysis.PlinkSeq.PlinkSeqWorker;
import seq.analysis.PlinkSeqUtils.PlinkSeqBurdenSummary;
import seq.analysis.PlinkSeqUtils.PseqProject;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.pathway.GenomeRegions;
import seq.pathway.Pathways;
import common.Logger;
import common.ext;
import filesys.GeneTrack;

/**
 * @author lane0212
 *
 *         Does a big plink-seq association with no filtering besides population mac<br>
 * 
 * 
 */
public class PlinkSeqMegs {

	public static void runBig(String vcf, String vpopFile, String resourceDirectory, String geneTrackFile, String keggPathwayFile, double maf, int numthreads, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.CASE_CONTROL, log);
		vpop.report();

		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);

		geneTrack.setGeneSetFilename(geneTrackFile);
		Pathways pathways = Pathways.load(keggPathwayFile);
		GenomeRegions gRegions = new GenomeRegions(geneTrack, pathways, log);

		VCFOps.verifyIndex(vcf, log);
		String locFile = resourceDirectory + ext.rootOf(gRegions.getGeneTrack().getGeneSetFilename()) + "_Gen.reg";
		PlinkSeqUtils.generatePlinkSeqLoc(gRegions, locFile, log);

		PlinkSeq plinkSeq = new PlinkSeq(false, true, log);

		PseqProject pseqProject = PlinkSeq.initialize(plinkSeq, ext.rootOf(vpop.getFileName()), ext.parseDirectoryOfFile(vpop.getFileName()), vcf, vpop, resourceDirectory, false, false, log);
		VCFFileReader reader = new VCFFileReader(vcf, true);
		int macFilter = (int) Math.round((float) VCFOps.getSamplesInFile(reader).length * maf);
		reader.close();
		//System.exit(1);

		PlinkSeqWorker[] complete = plinkSeq.fullGamutAssoc(pseqProject, new String[] { ext.rootOf(locFile) }, null, -1, macFilter + "", ext.rootOf(vpop.getFileName()), numthreads);
		PlinkSeqBurdenSummary[] summaries = new PlinkSeqBurdenSummary[1];
		int index = 0;
		for (int i = 0; i < complete.length; i++) {
			ANALYSIS_TYPES type = complete[i].getType();
			switch (type) {
			case BURDEN:
				String analysis = ext.rootOf(complete[i].getOutputFiles()[0]);
				analysis = analysis.replaceAll(".*" + ext.rootOf(locFile) + ".", "");
				PlinkSeqBurdenSummary plinkSeqBurdenSummary = new PlinkSeqBurdenSummary(analysis, complete[i].getOutputFiles()[0], log);
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
				log.reportTimeError("INVALID analysis type " + type);
				break;
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";
		// String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqTallyTest/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.errorRegions2.vcf";
		// String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqTallyTest/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vPopCaseControl.vcf";
		// String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/joint_genotypes.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.CUSHINGS.vcf.gz";
		String vcf = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/vcf/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf";


		
		String vpopFile = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Spector_Joint/vcf/pseqProj_tsai_spector_joint_AgilentCaptureRecal/vPopCaseControl.txt";
		String resourceDirectory = "/home/tsaim/public/bin/pseqRef/hg19/";
		Logger log = new Logger(ext.rootOf(vcf, false) + "tally.log");
		String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
		String keggPathwayFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/kegg.ser";
		String logfile = null;
		int numthreads = 24;
		double maf = 0.05;
		String usage = "\n" + "seq.analysis.VCFTallyPSeq requires 0-1 arguments\n";
		usage += "   (1) vcf file (i.e. file=" + vcf + " (default))\n" + "";
		usage += "   (2) vpop files, comma delimited (i.e. vpop=" + vpopFile + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vpop=")) {
				vpopFile = args[i].split("=")[1];
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
			runBig(vcf, vpopFile, resourceDirectory, geneTrackFile, keggPathwayFile, maf, numthreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
