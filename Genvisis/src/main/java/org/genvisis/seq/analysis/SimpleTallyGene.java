package org.genvisis.seq.analysis;

import java.io.File;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
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

		/**
		 * @param vcf vcf file to use
		 * @param seg seqment to tally
		 * @param name name of this segment
		 * @param vpopFile
		 * @param omimDir
		 * @param filter
		 * @param maf max maf
		 */
		private Params(	String vcf, Segment seg, String name, String vpopFile, String omimDir,
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
				Files.copyFile(ext.rootOf(vpopFile, false)+ ".lowerQualitySamples.txt",
												dir + ext.rootOf(vpopFile, true) + ".lowerQualitySamples.txt");
			}
			VCFSimpleTally.test(subVcf, new String[] {newVpop}, omimDir, null, null, maf, true, true,
													filter, false);
			return this;
		}

	}

	/**
	 * @param vcf
	 * @param vpop
	 * @param maf
	 * @param seg
	 * @param name
	 * @param omimDir
	 */
	public static void run(	String vcf, String vpop, double maf, Segment seg, String name,
													String omimDir) {

		WorkerHive<Params> hive = new WorkerHive<SimpleTallyGene.Params>(1, 1, new Logger());
		hive.addCallable(new Params(vcf, seg, name, vpop, omimDir, null, maf));
		hive.execute(true);


	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CLI c = new CLI(SimpleTallyGene.class);



		c.addArgWithDefault("vcf", "vcf to tally", "a.vcf");
		c.addArgWithDefault("vpop", "vpop to use", "a.vpop");
		c.addArgWithDefault("maf", "maf to use", "1.2");
		c.addArgWithDefault("segment", "UCSC segments , ; delimited",
												"chr18:20714428-20840534;chr17:7571720-7590868");
		c.addArgWithDefault("name", "typically gene name, comma delimited", "CABLES1,TP53");
		c.addArgWithDefault("omimDir", "omim directory", "/Volumes/Beta/ref/OMIM/");


		c.parseWithExit(args);

		String vcf = c.get("vcf");
		String vpop = c.get("vpop");
		double maf = c.getD("maf");

		Segment[] segs = Segment.getSegments(c.get("segment").split(";"));
		String[] names = c.get("name").split(",");
		String omimDir = c.get("omimDir");
		if (names.length != segs.length) {
			throw new IllegalArgumentException("Must have a name for each segment");
		}
		for (int i = 0; i < names.length; i++) {
			run(vcf, vpop, maf, segs[i], names[i], omimDir);
		}

	}



}


