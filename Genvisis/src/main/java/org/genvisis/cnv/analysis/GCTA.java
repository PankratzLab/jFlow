package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
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

	// QC NOTES
	// We chose the PLINK25 compact binary file format (∗.bed, ∗.bim, and ∗.fam)
	// as the input data format for GCTA because of its popularity in the
	// genetics community and its efficiency of data storage. For the imputed
	// dosage data, we use the output files of the imputation program MACH32
	// (∗.mldose.gz and ∗.mlinfo.gz) as the inputs for GCTA. For the convenience
	// of analysis, we provide options to extract a subset of individuals and/or
	// SNPs and to filter SNPs based on certain criteria, such as chromosome
	// position, minor allele frequency (MAF), and imputation R2 (for the
	// imputed data). However, we do not provide functions for a thorough
	// quality control (QC) of the data, such as Hardy-Weinberg equilibrium test
	// and missingness, because these functions have been well developed in many
	// other genetic analysis packages, e.g., PLINK, GenABEL,33 and SNPTEST.34
	// We assume that the data have been cleaned by a standard QC process before
	// entering into GCTA.

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

	// gcta64 --grm test --grm-cutoff 0.025 --make-grm --out test_rm025
	//

	// We provide a function to iteratively exclude one individual of a pair
	// whose relationship is greater than a specified cutoff value, e.g., 0.025,
	// while retaining the maximum number of individuals in the data. For data
	// collected from family or twin studies, we recommend that users estimate
	// the genetic relationships with all of the autosomal SNPs and then use
	// this option to exclude close relatives. The reason for exclusion is that
	// the objective of the analysis is to estimate genetic variation captured
	// by all the SNPs, just as GWAS does for single SNPs. Including close
	// relatives, such as parent-offspring pairs and siblings, would result in
	// the estimate of genetic variance being driven by the phenotypic
	// correlations for these pairs (just as in pedigree analysis), and this
	// estimate could be a biased estimate of total genetic variance, for
	// example because of common environmental effects. Even if the estimate is
	// not biased, its interpretation is different from the estimate from
	// “unrelated” individuals: a pedigree-based estimator captures the
	// contribution from all causal variants (across the entire allele frequency
	// spectrum), whereas our method captures the contribution from causal
	// variants that are in LD with the genotyped SNPs.

	// Nother source

	// https://espace.library.uq.edu.au/view/UQ:342517/UQ342517_OA.pdf
	private static boolean removeCrypticRelated(String inputGrm, String outputGrm, double grmCutoff, Logger log) {
		// String[] inputs = new String[] { inputGrm };

		String[] outputs = new String[] { outputGrm + ".grm.N.bin", outputGrm + ".grm.id" };
		ArrayList<String> command = new ArrayList<String>();
		command.add("gcta64");
		command.add("--grm");
		command.add(inputGrm);
		command.add("--grm-cutoff");
		command.add(grmCutoff + "");
		command.add("--make-grm");
		command.add("--out");
		command.add(outputGrm);
		boolean success = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", null, outputs, true, false,
				false, log);
		return success;
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
		private ArrayList<String> source;
		private ArrayList<Double> variance;
		private ArrayList<Double> se;

		private VarianceResult(String summaryFile, boolean success, Logger log) {
			super();
			this.summaryFile = summaryFile;
			this.success = success;
			loadResults(log);
		}

		// Output file format
		// test.hsq (rows are
		// header line;
		// name of genetic variance, estimate and standard error (SE);
		// residual variance, estimate and SE;
		// phenotypic variance, estimate and SE;
		// ratio of genetic variance to phenotypic variance, estimate and SE;
		// log-likelihood;
		// sample size). If there are multiple GRMs included in the REML
		// analysis, there will be multiple rows for the genetic variance (as
		// well as their ratios to phenotypic variance) with the names of V(1),
		// V(2), … .
		// Source Variance SE
		// V(1) 0.389350 0.161719
		// V(e) 0.582633 0.160044
		// Vp 0.971984 0.031341
		// V(1)/Vp 0.400573 0.164937
		// The estimate of variance explained on the observed scale is
		// transformed to that on the underlying scale:
		// (Proportion of cases in the sample = 0.5; User-specified disease
		// prevalence = 0.1)
		// V(1)/Vp_L 0.657621 0.189123
		// logL -945.65
		// logL0 -940.12
		// LRT 11.06
		// Pval 4.41e-4
		// n 2000
		//
		// Source Variance SE
		// V(G) 0.006697 0.000593
		// V(e) 0.014934 0.000567
		// Vp 0.021631 0.000291
		// V(G)/Vp 0.309599 0.026340
		// logL 16299.208
		// logL0 16219.785
		// LRT 158.847
		// df 1
		// Pval 0
		// n 11443
		private void loadResults(Logger log) {
			if (success) {
				try {
					BufferedReader reader = Files.getAppropriateReader(summaryFile);
					reader.readLine();
					this.source = new ArrayList<String>();
					this.variance = new ArrayList<Double>();
					this.se = new ArrayList<Double>();

					while (reader.ready()) {
						String[] line = reader.readLine().trim().split("\t");
						source.add(line[0]);
						variance.add(Double.parseDouble(line[1]));
						if (line.length > 2) {
							se.add(Double.parseDouble(line[2]));
						} else {
							se.add(Double.NaN);
						}
					}
					reader.close();
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
		}

	}

	// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4295936/
	// gcta64 --grm test --pheno test.phen --reml --qcovar test_10PCs.txt --out
	// test --thread-num 10
	private static VarianceResult determineVarianceExplained(String inputGrm, String output, String phenoFile,
			String covarFile, int numthreads, Logger log) {
		String[] inputs = new String[] { inputGrm + ".grm.bin", inputGrm + ".grm.N.bin", inputGrm + ".grm.id" };
		String summaryFile = output + ".hsq";
		// if (!Files.exists(summaryFile)) {
		// System.out.println(summaryFile);
		// System.exit(1);
		// }
		String[] outputs = new String[] { summaryFile };
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
		return new VarianceResult(summaryFile, success, log);
	}

	public enum PHENO_TYPE {
		PRE_PROCESSED, MITO_FILE;
	}

	private static ArrayList<VarianceResult> processMitoFile(final Project proj, String mitoFile,
			final String mergedGRM, final String covarFile, int numThreads, Logger log) {
		ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
		builder.sampleBased(true);
		builder.requireAll(true);
		builder.dataKeyColumnIndex(0);
		builder.treatAllNumeric(true);
		try {
			ExtProjectDataParser parser = builder.build(proj, mitoFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			final String resultsDir = ext.parseDirectoryOfFile(mergedGRM) + ext.rootOf(mitoFile) + "/";
			new File(resultsDir).mkdirs();
			SampleData sampleData = proj.getSampleData(0, false);
			String[] samples = proj.getSamples();
			ArrayList<String> fidIID = new ArrayList<String>();
			for (int i = 0; i < samples.length; i++) {
				fidIID.add(sampleData.lookup(samples[i])[1]);
			}

			WorkerHive<VarianceResult> hive = new WorkerHive<GCTA.VarianceResult>(numThreads, 10, proj.getLog());
			for (int i = 0; i < Math.min(125, parser.getNumericDataTitles().length); i++) {
				final String current = parser.getNumericDataTitles()[i];
				ArrayList<String> pheno = new ArrayList<String>();
				double[] data = parser.getNumericDataForTitle(current);
				for (int j = 0; j < data.length; j++) {
					pheno.add(fidIID.get(j) + "\t" + data[j]);
				}
				final String phenoFile = resultsDir + current + ".txt";
				Files.writeArrayList(pheno, phenoFile);
				hive.addCallable(new Callable<GCTA.VarianceResult>() {

					@Override
					public VarianceResult call() throws Exception {
						return determineVarianceExplained(mergedGRM, resultsDir + current, phenoFile, covarFile, 1,
								proj.getLog());
					}
				});

			}
			hive.execute(true);
			return hive.getResults();

		} catch (FileNotFoundException nfe) {
			log.reportException(nfe);
		}
		return null;
	}

	private static void run(Project proj, String sampFile, String[] phenoFiles, PHENO_TYPE pType, int pcCovars,
			double grmCutoff, int numthreads) {
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

		ArrayList<GRM> grms = splitRunGCTA(outDir + "quality_control/ld_pruning/gcta", plinkRootQC, 0.01, 1, numthreads,
				proj.getLog());
		String mergedGRM = plinkRootQC + "_merge";
		boolean success = mergeGRMs(grms, mergedGRM, numthreads, proj.getLog());

		String covarFile = null;
		if (success) {
			String mergeRmGRM = mergedGRM + "_rm_" + grmCutoff;
			success = removeCrypticRelated(mergedGRM, mergeRmGRM, grmCutoff, proj.getLog());

			if (pcCovars > 0) {
				success = generatePCACovars(mergeRmGRM, mergeRmGRM, pcCovars, numthreads, proj.getLog());
				if (success) {
					covarFile = mergeRmGRM + ".eigenvec";
				} else {
					throw new IllegalStateException("Failed to generate grms");
				}
			}
			switch (pType) {
			case MITO_FILE:
				String output = mergeRmGRM + "_mitoResults.txt";
				PrintWriter writer = Files.getAppropriateWriter(output);

				for (int i = 0; i < phenoFiles.length; i++) {
					ArrayList<VarianceResult> results = processMitoFile(proj, phenoFiles[i], mergeRmGRM, covarFile,
							numthreads, proj.getLog());
					if (i == 0) {
						writer.println("PhenoFile\tresultFile\tPC\t"
								+ Array.toStr(Array.toStringArray(results.get(0).source)));
					}
					for (VarianceResult varianceResult : results) {
						writer.println(ext.removeDirectoryInfo(phenoFiles[i]) + "\t" + varianceResult.summaryFile + "\t"
								+ ext.removeDirectoryInfo(varianceResult.summaryFile).replaceAll("\\.hsq", "")
										.replaceAll("PC", "")
								+ "\t" + Array.toStr(Array.toDoubleArray(varianceResult.variance)));
					}
				}
				writer.close();
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
		String[] phenoFiles = null;
		// String out = null;
		double grmCutoff = 0.025;
		Project proj;
		int numthreads = 24;
		int pcCovars = 10;

		String usage = "\n" + "cnv.analysis.GCTA requires 1-3 arguments\n"
				+ "   (1) project properties filename (i.e. proj="
				+ org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
				+ "   (2) samples to use (i.e. samps= (defaults to all samples))\n"
				+ PSF.Ext.getNumThreadsCommand(3, numthreads) + "\n"
				+ "   (4) phenotype file(s) (i.e. pheno= (no default))\n"
				+ "   (5) number of Principal components to compute for covariate use, set to -1 to skip (i.e. pcCovars="
				+ pcCovars + " (default))\n" + "" + "   (6) grm relatedness cutoff (i.e. grmCutoff=" + grmCutoff
				+ " (default))\n" + "";

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
				phenoFiles = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("pcCovars=")) {
				pcCovars = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("grmCutoff=")) {
				grmCutoff = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			}
			// else if (args[i].startsWith("out=")) {
			// out = args[i].split("=")[1];
			// numArgs--;
			// }
			else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj = new Project(filename, false);
			run(proj, sampFile, phenoFiles, PHENO_TYPE.MITO_FILE, pcCovars, grmCutoff, numthreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}