package seq.analysis;

import java.io.File;
import java.util.concurrent.Callable;

import seq.manage.VCFOps;
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

		public Params(String vcf, Segment seg, String name, String vpopFile, String omimDir) {
			super();
			this.vcf = vcf;
			this.seg = seg;
			this.name = name;
			this.vpopFile = vpopFile;
			this.omimDir = omimDir;
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
			VCFSimpleTally.test(subVcf, new String[] { newVpop }, omimDir, null, null, 1.2, true, true, null);
			return this;
		}

	}

	private static void run() {
		String vcfGermline = "D:/data/Project_Tsai_21_25_26_28_spector/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz";
		String vcfTumor = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/tnMatch.merged.renamed.hg19_multianno.eff.gatk.anno_charge.sed1000g.finalMerge.vcf.gz";

		Segment[] segs = new Segment[] { new Segment("chr18:20,714,528-20,840,434"), new Segment("chr15:50714579-50795277") };
		String[] names = new String[] { "CABLES1", "USP8", };

		String vpopFileGermline = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/CUSHING_FREQ.vpop";
		String vpopFileTumor = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/CUSHINGS_TUMOR.vpop";

		WorkerHive<Params> hive = new WorkerHive<SimpleTallyGene.Params>(8, 1, new Logger());
		String omimDir = "C:/bin/ref/OMIM/";
		for (int i = 0; i < names.length; i++) {
			hive.addCallable(new Params(vcfGermline, segs[i], names[i], vpopFileGermline, omimDir));
			hive.addCallable(new Params(vcfTumor, segs[i], names[i], vpopFileTumor, omimDir));
		}
		hive.execute(true);
	}

	public static void main(String[] args) {
		run();

	}

}
