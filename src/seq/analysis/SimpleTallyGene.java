package seq.analysis;

import java.io.File;
import java.util.concurrent.Callable;

import seq.manage.VCFOps;
import seq.qc.FilterNGS;
import seq.qc.FilterNGS.FILTER_GENERATION_TYPE;
import seq.qc.FilterNGS.VariantContextFilter;
import common.Files;
import common.Logger;
import common.WorkerHive;
import common.ext;
import filesys.Segment;

public class SimpleTallyGene {

	private static class Params implements Callable<Params> {
		private String vcf;
		private Segment seg;
		private String name;
		private String vpopFile;
		private String omimDir;
		private VariantContextFilter filter;

		public Params(String vcf, Segment seg, String name, String vpopFile, String omimDir, VariantContextFilter filter) {
			super();
			this.vcf = vcf;
			this.seg = seg;
			this.name = name;
			this.vpopFile = vpopFile;
			this.omimDir = omimDir;
			this.filter = filter;

		}

		@Override
		public Params call() throws Exception {
			String dir = ext.rootOf(vpopFile, false) + "_" + name + "/";
			new File(dir).mkdirs();
			Logger log = new Logger(dir + "log.log");

			String newVpop = dir + ext.removeDirectoryInfo(vpopFile);
			String segFile = dir + name + ".segs";
			Files.write(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop(), segFile);
			String subVcf = VCFOps.extractSegments(vcf, segFile, 100, null, ext.rootOf(vpopFile, false) + "_extractedVcfs/", false, true, true, false, null, 1, log);
			Files.copyFile(vpopFile, newVpop);
			VCFSimpleTally.test(subVcf, new String[] { newVpop }, omimDir, null, null, 1.2, true, true, filter);
			return this;
		}

	}

	private static void run() {
		String vcfGermline = "D:/data/Project_Tsai_21_25_26_28_spector/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz";

		Segment[] segs = new Segment[] { new Segment("chr6:31,371,371-31,383,090"), new Segment("chr17:7,571,720-7,590,868") };
		String[] names = new String[] { "MICA","TP53" };
		String vpopFileGermlineOsteo = "D:/data/logan/OSv2_seq/SRGAP2/OSTEO_OFF_INHERIT.vpop";
		String vpopFileGermlineCushing = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/CUSHING_FREQ.vpop";
		String vpopFileGermlineEPP = "D:/data/Project_Tsai_21_25_26_28_spector/Freq/EPP.vpop";
		String newDir = "D:/data/Project_Tsai_21_25_26_28_spector/Look_Freq/";
		new File(newDir).mkdirs();
		Files.copyFile(vpopFileGermlineCushing, newDir + ext.removeDirectoryInfo(vpopFileGermlineCushing));
		Files.copyFile(vpopFileGermlineOsteo, newDir + ext.removeDirectoryInfo(vpopFileGermlineOsteo));
		Files.copyFile(vpopFileGermlineEPP, newDir + ext.removeDirectoryInfo(vpopFileGermlineEPP));

		vpopFileGermlineOsteo = newDir + ext.removeDirectoryInfo(vpopFileGermlineOsteo);
		vpopFileGermlineCushing = newDir + ext.removeDirectoryInfo(vpopFileGermlineCushing);
		vpopFileGermlineEPP = newDir + ext.removeDirectoryInfo(vpopFileGermlineEPP);

		WorkerHive<Params> hive = new WorkerHive<SimpleTallyGene.Params>(2, 1, new Logger());
		String omimDir = "C:/bin/ref/OMIM/";
		for (int i = 0; i < names.length; i++) {
			hive.addCallable(new Params(vcfGermline, segs[i], names[i], vpopFileGermlineOsteo, omimDir, null));
			hive.addCallable(new Params(vcfGermline, segs[i], names[i], vpopFileGermlineCushing, omimDir, null));
			hive.addCallable(new Params(vcfGermline, segs[i], names[i], vpopFileGermlineEPP, omimDir, null));

		}
		hive.execute(true);
		// String vcfTumor = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/tnMatch.merged.renamed.hg19_multianno.eff.gatk.anno_charge.sed1000g.finalMerge.vcf.gz";

		// String vpopFileGermline = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/CUSHING_FREQ.vpop";
		// String vpopFileTumor = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/CUSHINGS_TUMOR.vpop";
	}

	public static void main(String[] args) {
		run();// TODO, cmdline if use this again

	}

}
