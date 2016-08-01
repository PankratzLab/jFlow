package org.genvisis.cnv.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;

/**
 * Class for running GCTA (http://cnsgenomics.com/software/gcta/)
 * 
 *
 */
public class GCTA {

	// gcta64 --mgrm grm_chrs.txt --make-grm --out test
	private static void mergeGRMs(ArrayList<GRM> grms, Logger log) {

	}

	/**
	 * @param plinkRoot
	 *            e.g. test.bed, test.bim and test.fam.
	 * @param output
	 *            output root
	 * @param maf
	 *            minor allele freq threshold
	 * @param numChrthreads
	 *            number of threads to use within a chromosome
	 * @param numBetweenThreads
	 *            number of chromosomes to analyze at once
	 * @param log
	 */
	private static ArrayList<GRM> splitRunGCTA(final String plinkRoot, String output, final double maf,
			final int numChrthreads, int numBetweenThreads, final Logger log) {
		WorkerHive<GRM> hive = new WorkerHive<GCTA.GRM>(numBetweenThreads, 10, log);

		for (int i = 1; i < 23; i++) {
			final String chrOut = output + "_chr" + i;
			final byte chr = (byte) i;
			hive.addCallable(new Callable<GCTA.GRM>() {

				@Override
				public GRM call() throws Exception {

					return runGCTA(plinkRoot, chrOut, maf, chr, numChrthreads, log);
				}
			});
		}
		hive.execute(true);
		ArrayList<GRM> results = hive.getResults();
		return results;
	}

	private static class GRM {
		private boolean success;
		private String[] output;

		private GRM(boolean success, String[] output) {
			super();
			this.success = success;
			this.output = output;
		}

	}

	// gcta64 --bfile test --chr 1 --maf 0.01 --make-grm --out test_chr1
	// --thread-num 10
	/**
	 * Focusing on generating the grm for a single autosomal chromosome
	 */
	private static GRM runGCTA(String plinkRoot, String output, double maf, byte chr, int numthreads, Logger log) {

		String[] inputs = PSF.Plink.getPlinkBedBimFam(plinkRoot);
		ArrayList<String> command = new ArrayList<String>();
		command.add("gcta64");
		command.add("--bfile");
		command.add(plinkRoot);
		if (chr > 0) {
			command.add("--chr");
			command.add(chr + "");
		} else {
			log.reportTimeWarning("GCTA running on full data set, consider breaking the analyis up by chromosome");
		}
		command.add("--maf");
		command.add(maf + "");
		command.add("--out");
		command.add(output);
		command.add("--thread-num");
		command.add(numthreads + "");
		command.add("--make-grm");
		boolean success = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, null, true, true,
				false, log);

		return new GRM(success, null);
	}

	private static void test(Project proj, String[] samples, String[] markers) {
		String outDir = proj.PROJECT_DIRECTORY.getValue()+"gcta/";
		new File(outDir).mkdirs();
		String fakePed = outDir+"gctaPedigree.dat";
		proj.PEDIGREE_FILENAME.setValue(fakePed);
		Pedigree.build(proj, null, samples, false);

	}

}

// private Project proj;
// /**
// * Samples to use, should already have been qc'ed
// */
// private String[] samples;
// /**
// * Markers to use, also should already have been qc'ed
// */
// private String[] markers;
//

// private String plinkRoot;
//
// /**
// * Maf filter for any markers used
// */
// private double maf;
// /**
// * Even remotely related pairs of individuals (genetic similarity greater
// * than 0.025, which represents fifth-degree relatives) are excluded so
// that
// * chance genetic similarity is used as a random effect in a linear mixed
// * model.
// */
// private double crypticRelatedParam;
// /**
// * Full path to plink root (i.e /home/usr/data/plinkFileRoot
// */
//
// /**
// * Set to -1 to run all chrs
// */
// private byte chr;
// private int numthreads;
//
// private String output;
