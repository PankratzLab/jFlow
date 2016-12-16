package org.genvisis.seq.analysis;

import java.io.File;
import java.util.concurrent.Callable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;

public class SimpleTallyGene {

	private static class Params implements Callable<Params> {
		private final String vcf;
		private final Segment seg;
		private final String name;
		private final String vpopFile;
		private final String omimDir;
		private final VariantContextFilter filter;
		private final double maf;

		public Params(String vcf, Segment seg, String name, String vpopFile, String omimDir,
									VariantContextFilter filter, double maf) {
			super();
			this.vcf = vcf;
			this.seg = seg;
			this.name = name;
			this.vpopFile = vpopFile;
			this.omimDir = omimDir;
			this.filter = filter;
			this.maf = maf;

		}

		@Override
		public Params call() throws Exception {
			String dir = ext.rootOf(vpopFile, false) + "_" + name + "/";
			new File(dir).mkdirs();
			Logger log = new Logger(dir + "log.log");

			String newVpop = dir + ext.removeDirectoryInfo(vpopFile);
			String segFile = dir + name + ".segs";
			Files.write(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop(), segFile);
			String subVcf = VCFOps.extractSegments(	vcf, segFile, 100, null,
																							ext.rootOf(vpopFile, false) + "_extractedVcfs/",
																							false, true, true, false, null, 1, log);
			Files.copyFile(vpopFile, newVpop);
			if (Files.exists(ext.rootOf(vpopFile, false) + ".lowerQualitySamples.txt")) {
				Files.copyFile(ext.rootOf(vpopFile, false)	+ ".lowerQualitySamples.txt",
												dir + ext.rootOf(vpopFile, true) + ".lowerQualitySamples.txt");
			}
			VCFSimpleTally.test(subVcf, new String[] {newVpop}, omimDir, null, null, maf, true, true,
													filter, false);
			return this;
		}

	}

	public static void run(String vcf, String vpop, double[] mafs) {

		String[] names = new String[] {"Mito"};
		Segment[] segs = new Segment[] {new Segment("chrM:1-20000")};
		mafs = mafs == null ? new double[] {1.2} : mafs;
		String tnVpop = "/Volumes/Beta/data/Cushings/mito/CUSHINGS_TUMOR.vpop";

		WorkerHive<Params> hive = new WorkerHive<SimpleTallyGene.Params>(1, 1, new Logger());
		String omimDir = "/Volumes/Beta/ref/OMIM/";
		for (int i = 0; i < names.length; i++) {
			for (double maf : mafs) {
				hive.addCallable(new Params(vcf, segs[i], names[i], vpop, omimDir, null, maf));
				if (maf == 0) {
					hive.addCallable(new Params(vcf, segs[i], names[i], tnVpop, omimDir, null, maf));
				}
			}
		}
		hive.execute(true);


	}

	public static void main(String[] args) {


		run(null, null, null);// TODO, cmdline if use this again

	}



}


