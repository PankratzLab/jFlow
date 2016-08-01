package org.genvisis.cnv.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.gwas.Qc;

/**
 * Class for running GCTA (http://cnsgenomics.com/software/gcta/)
 * 
 *
 */
public class GCTA {

	
	//TODO
	// gcta64 --mgrm grm_chrs.txt --make-grm --out test
	private static void mergeGRMs(ArrayList<GRM> grms, String output, Logger log) {
		
	}
	
	//TODO
	private static void generatePCACovars(){
		
	}
	
	
	//TODO
//	gcta64 --grm test --grm-cutoff 0.025 --make-grm --out test_rm025
//
	private static void removeCrypticRelated(){
		
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

	private static void run(Project proj, String sampFile, String markerFile, int numthreads) {
		String[] samples = sampFile == null ? null
				: HashVec.loadFileToStringArray(sampFile, false, false, new int[] { 0 }, false, true, "\t");

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "gcta/";
		new File(outDir).mkdirs();
		String fakePed = outDir + "gctaPedigree.dat";
		proj.PEDIGREE_FILENAME.setValue(fakePed);

		String[] plinks = PSF.Plink.getPlinkBedBimFam(outDir + "gcta");
		if (!Files.exists("", plinks)) {
			Pedigree.build(proj, null, samples, false);
			PlinkData.saveGenvisisToPlinkBedSet(proj, "gcta/gcta", null, markerFile, -1, true);
		}
		Qc.fullGamut(outDir, "gcta", false, proj.getLog());

		splitRunGCTA("FILLIN", outDir + "gcta", 0.01, 1, numthreads, proj.getLog());

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String sampFile = null;
		String markerFile = null;

		String out = null;
		boolean overwrite = false;
		Project proj;
		int numthreads = 24;

		String usage = "\n" + "cnv.analysis.GCTA requires 1-3 arguments\n"
				+ "   (1) project properties filename (i.e. proj="
				+ org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
				+ "   (2) samples to use (i.e. samps= (defaults to all samples))\n"
				+ PSF.Ext.getNumThreadsCommand(3, numthreads) + "\n"

				// + " (3) markers to use (i.e. markers= (defaults to all
				// markers))\n"

				+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("samps=")) {
				sampFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				out = args[i].split("=")[1];
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
			proj = new Project(filename, false);
			run(proj, sampFile, markerFile, numthreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
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
