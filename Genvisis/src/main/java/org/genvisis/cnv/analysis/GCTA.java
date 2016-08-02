package org.genvisis.cnv.analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
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

	// Lil note about whether to remove related samples
	/**
	 * Can I use GCTA to estimate the variance explained by a subset of SNP in
	 * family data? Yes, you can. GCTA does not assume that the individuals
	 * should be unrelated. The reason for excluding close-relatives in Yang et
	 * al. (Nat. Genet. 2010 and 2011) is because we do not want our estimates
	 * to be confounded with some possible shared environment effects and the
	 * effects of some possible causal variants that are not tagged by the SNPs
	 * but captured by pedigree information. If you are interested in the
	 * variance explained by a subset of SNPs in family data, you could fit the
	 * genetic relationship matrix (GRM) estimated from these SNPs along with a
	 * matrix of pedigree structure using the option --mgrm when running the
	 * REML analysis (--reml). Alternatively, we could fit the GRM of the subset
	 * of SNPs together with another GRM estimated from the SNPs in the rest of
	 * the genome.
	 *
	 */

	// gcta64 --mgrm grm_chrs.txt --make-grm --out test
	private static boolean mergeGRMs(ArrayList<GRM> grms, String output, int numthreads, Logger log) {
		ArrayList<String> chrGRMs = new ArrayList<String>();
		for (GRM grm : grms) {
			if (!grm.success) {
				log.reportTimeError("GRM generation has failed (" + grm.grmFile + ")");
				return false;
			}
			chrGRMs.add(grm.grmFile);
		}
		String chrListFile = output + "_chrsGRM.txt";
		Files.writeArrayList(chrGRMs, chrListFile);

		// String[] inputs = Array.toStringArray(chrGRMs);

		String[] outputs = new String[] { output + ".grm.N.bin", output + ".grm.id" };
		ArrayList<String> command = new ArrayList<String>();
		command.add("gcta64");
		command.add("--mgrm");
		command.add(chrListFile);
		command.add("--make-grm");
		command.add("--out");
		command.add(output);
		command.add("--thread-num");
		command.add(numthreads + "");
		boolean success = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", null, outputs, true, false,
				false, log);
		return success;

	}

	// gcta64 --grm test --keep test.indi.list --pca 20 --out test
	private static boolean generatePCACovars(String inputGrm, String output, int numPCs, int numthreads, Logger log) {
		String[] inputs = new String[] { inputGrm + ".grm.bin", inputGrm + ".grm.N.bin", inputGrm + ".grm.id" };

		String[] outputs = new String[] { output + ".eigenval", output + ".eigenvec" };
		ArrayList<String> command = new ArrayList<String>();
		command.add("gcta64");
		command.add("--grm");
		command.add(inputGrm);

		command.add("--pca");
		command.add(numPCs + "");

		command.add("--out");
		command.add(output);
		command.add("--thread-num");
		command.add(numthreads + "");
		boolean success = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
				false, false, log);
		return success;
	}

	// // TODO
	// // gcta64 --grm test --grm-cutoff 0.025 --make-grm --out test_rm025
	// //
	// private static boolean removeCrypticRelated(String inputGrm, String
	// outputGrm, double grmCutoff, Logger log) {
	// String[] inputs = new String[] { inputGrm };
	//
	// String[] outputs = new String[] { outputGrm };
	// ArrayList<String> command = new ArrayList<String>();
	// command.add("gcta64");
	// command.add("--grm");
	// command.add(inputGrm);
	// command.add("--grm-cutoff");
	// command.add(grmCutoff + "");
	// command.add("--make-grm");
	// command.add("--out");
	// command.add(outputGrm);
	// boolean success =
	// CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "",
	// inputs, outputs, true,
	// false, false, log);
	// return success;
	// }

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
		private String grmFile;

		private GRM(boolean success, String grmFile) {
			super();
			this.success = success;
			this.grmFile = grmFile;
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

		String[] outputs = new String[] { output + ".grm.N.bin", output + ".grm.id" };
		boolean success = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
				false, false, log);

		return new GRM(success, output);
	}

	private static class VarianceResult {
		private String summaryFile;
		private boolean success;

		public VarianceResult(String summaryFile, boolean success) {
			super();
			this.summaryFile = summaryFile;
			this.success = success;
		}

	}

	// gcta64 --grm test --pheno test.phen --reml --qcovar test_10PCs.txt --out
	// test --thread-num 10
	private static VarianceResult determineVarianceExplained(String inputGrm, String output, String phenoFile,
			String covarFile, int numthreads, Logger log) {
		String[] inputs = new String[] { inputGrm + ".grm.bin", inputGrm + ".grm.N.bin", inputGrm + ".grm.id" };
		String[] outputs = new String[] { "FDDD" };
		ArrayList<String> command = new ArrayList<String>();
		command.add("gcta64");
		command.add("--grm");
		command.add(inputGrm);

		command.add("--pheno");
		command.add(phenoFile + "");

		command.add("--reml");
		if (covarFile != null) {
			command.add("--qcovar");
			command.add(covarFile);
		}
		command.add("--out");
		command.add(output);
		command.add("--thread-num");
		command.add(numthreads + "");
		boolean success = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
				false, false, log);
		return new VarianceResult(null, success);
	}

	public enum PHENO_TYPE {
		PRE_PROCESSED, MITO_FILE;
	}

	private static void processMitoFile(Project proj, String mitoFile, String mergedGRM, String covarFile, Logger log) {
		ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
		builder.sampleBased(true);
		builder.requireAll(true);
		builder.dataKeyColumnIndex(0);
		builder.treatAllNumeric(true);
		try {
			ExtProjectDataParser parser = builder.build(proj, mitoFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			String resultsDir = ext.parseDirectoryOfFile(mergedGRM) + ext.removeDirectoryInfo(mitoFile);
			new File(resultsDir).mkdirs();
			SampleData sampleData = proj.getSampleData(0, false);
			String[] samples = proj.getSamples();
			ArrayList<String> fidIID = new ArrayList<String>();
			for (int i = 0; i < samples.length; i++) {
				fidIID.add(sampleData.lookup(samples[i])[1]);
			}
			for (int i = 0; i < parser.getNumericDataTitles().length; i++) {
				String current = parser.getNumericDataTitles()[i];
				ArrayList<String> pheno = new ArrayList<String>();
				double[] data = parser.getNumericDataForTitle(current);
				for (int j = 0; j < data.length; j++) {
					pheno.add(fidIID.get(j) + "\t" + data[j]);
				}
				String phenoFile = resultsDir + current + ".txt";
				Files.writeArrayList(pheno, phenoFile);
			}
		} catch (FileNotFoundException nfe) {
			log.reportException(nfe);
		}
	}

	private static void run(Project proj, String sampFile, String phenoFile, PHENO_TYPE pType, int pcCovars,
			int numthreads) {
		String[] samples = sampFile == null ? null
				: HashVec.loadFileToStringArray(sampFile, false, false, new int[] { 0 }, false, true, "\t");

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "gcta/";
		String plinkRoot = outDir + "gcta";
		new File(outDir).mkdirs();
		String fakePed = outDir + "gctaPedigree.dat";
		proj.PEDIGREE_FILENAME.setValue(fakePed);

		String[] plinks = PSF.Plink.getPlinkBedBimFam(plinkRoot);
		if (!Files.exists("", plinks)) {
			String nonCNFile = outDir + "markersToQC.txt";
			Files.writeList(proj.getNonCNMarkers(), nonCNFile);
			Pedigree.build(proj, null, samples, false);
			PlinkData.saveGenvisisToPlinkBedSet(proj, "gcta/gcta", null, nonCNFile, -1, true);
		}
		Qc.fullGamut(outDir, "gcta", false, proj.getLog());

		String plinkRootQC = outDir + "gcta_qc";

		ArrayList<GRM> grms = splitRunGCTA(outDir + "quality_control/genome/gcta", plinkRootQC, 0.01, 1, numthreads,
				proj.getLog());
		String mergedGRM = plinkRootQC + "_merge";
		boolean success = mergeGRMs(grms, mergedGRM, numthreads, proj.getLog());

		String covarFile = null;
		if (success) {
			if (pcCovars > 0) {
				success = generatePCACovars(mergedGRM, mergedGRM, pcCovars, numthreads, proj.getLog());
				if (success) {
					covarFile = mergedGRM + ".eigenvec";
				} else {
					throw new IllegalStateException("Failed to generate grms");
				}
			}
			switch (pType) {
			case MITO_FILE:

				break;
			case PRE_PROCESSED:
				// determineVarianceExplained(mergedGRM, mergedGRM, phenoFile,
				// covarFile, numthreads, proj.getLog());
			default:
				throw new IllegalArgumentException(pType + " is not implemented yet");

			}
		} else {
			throw new IllegalStateException("Failed to generate grms");
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String sampFile = null;
		String phenoFile = null;
		String out = null;
		Project proj;
		int numthreads = 24;
		int pcCovars = 10;

		String usage = "\n" + "cnv.analysis.GCTA requires 1-3 arguments\n"
				+ "   (1) project properties filename (i.e. proj="
				+ org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
				+ "   (2) samples to use (i.e. samps= (defaults to all samples))\n"
				+ PSF.Ext.getNumThreadsCommand(3, numthreads) + "\n"
				+ "   (4) phenotype file (i.e. pheno= (no default))\n"
				+ "   (5) number of Principal components to compute for covariate use, set to -1 to skip (i.e. pcCovars="
				+ pcCovars + " (default))\n" + "";

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
			} else if (args[i].startsWith("pheno=")) {
				phenoFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pcCovars=")) {
				pcCovars = ext.parseIntArg(args[i]);
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
			run(proj, sampFile, phenoFile, PHENO_TYPE.MITO_FILE, pcCovars, numthreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}