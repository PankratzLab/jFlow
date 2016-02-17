package seq.analysis;

import java.io.File;

import seq.manage.VCFOps;
import common.Files;
import common.Logger;
import common.ext;
import filesys.Segment;

public class SimpleTallyGene {

	private static void run() {
		String vcf = "D:/data/Project_Tsai_21_25_26_28_spector/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz";

		Segment[] segs = new Segment[] { new Segment("chr15:50714579-50795277") };

		String[] names = new String[] { "USP8" };
		String vpopFile = "D:/data/Project_Tsai_21_25_26_28_spector/Cushings/candidateGenes/CUSHING_FREQ.vpop";
		String omimDir = "C:/bin/ref/OMIM/";
		double[] mafs = new double[] { 1.2 };
		for (int i = 0; i < names.length; i++) {

			String dir = ext.parseDirectoryOfFile(vpopFile) + names[i] + "/";
			Logger log = new Logger(dir + "log.log");
			String newVpop = dir + ext.removeDirectoryInfo(vpopFile);
			String segFile = dir + names[i] + ".segs";
			Files.write(segs[i].getChr() + "\t" + segs[i].getStart() + "\t" + segs[i].getStop(), segFile);
			String subVcf = VCFOps.extractSegments(vcf, segFile, 100, null, dir, false, true, true, false, null, 1, log);
			new File(dir).mkdirs();
			Files.copyFile(vpopFile, newVpop);
			for (int j = 0; j < mafs.length; j++) {
				VCFSimpleTally.test(subVcf, new String[] { newVpop }, omimDir, null, null, mafs[j], true, true, null);
			}
		}

	}

	public static void main(String[] args) {
		run();

	}

}
