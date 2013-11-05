package gwas;

import java.io.*;
import java.util.*;

import stats.Rscript;
import common.*;

public class SkatMetaPrimary {

	public static void batch(String cohort, String genos, String phenoFilename, String snpInfo) {
		String phenoDir;
		String phenoRoot;
		String resultDir;
		String currentGeno;
		String currentSnpInfo;
		String rCode;
		String batchDir;
		PrintWriter out;
		Vector<String> v;
		String filename, commands;
		String[][] iterations;
		boolean foundGenos;
		boolean foundSnpInfo;

		if (! new File(phenoFilename).exists()) {
			System.err.println("Error - File not found " + phenoFilename);
			return;
		}

		// create directory called dir+root+"/"; example: "c:/diffpath/pheno_F7_studyIDs/"
		phenoDir = ext.parseDirectoryOfFile(new File(phenoFilename).getAbsolutePath());
		phenoRoot = ext.rootOf(phenoFilename);
		if(! new File(phenoDir + phenoRoot + "/").exists() || ! new File(phenoDir + phenoRoot + "/").isDirectory()) {
			new File(phenoDir + phenoRoot + "/").mkdirs();
		}

		resultDir = phenoDir + phenoRoot + "/results/";
		if(! new File(resultDir).exists() || ! new File(resultDir).isDirectory()) {
			new File(resultDir).mkdirs();
		}

		batchDir = phenoDir + phenoRoot + "/batchFiles/";
		if(! new File(batchDir).exists() || ! new File(batchDir).isDirectory()) {
			new File(batchDir).mkdirs();
		}

		try {

			v = new Vector<String>();
			// generate batch files in dir+root+"/batchFiles/"; example: "c:/diffpath/pheno_F7_studyIDs/batchFiles/chr1.R"
			foundGenos = false;
			foundSnpInfo = false;
			for (int i = 1; i <= 26; i++) {
				currentGeno = ext.replaceAllWith(genos, "#", i+""); //Files.getAppropriateWriter();
				if (snpInfo.contains("#")) {
					currentSnpInfo = ext.replaceAllWith(snpInfo, "#", i+"");
				} else {
					currentSnpInfo = snpInfo;
				}
				if (new File(currentGeno).exists() && new File(currentSnpInfo).exists()) {
					foundGenos = true;
					foundSnpInfo = true;
					rCode = "library(\"skatMeta\")\n"
							+ "setwd(\"" + resultDir + "\")\n"
							+ "\n"
							+(currentSnpInfo.toLowerCase().endsWith(".rdata")
									?"obj_name <- load(\"" + currentSnpInfo + "\")\n"
									+"SNPInfo <- get(obj_name)\n"
											+"rm(list=obj_name)\n"
											+"rm(obj_name)\n"
									:"SNPInfo <- read.csv(\"" + currentSnpInfo + "\", header=T, as.is=T)\n"
									)
							+ "\n"
							
							+(genos.toLowerCase().endsWith(".rdata")
									?"genoName <- load(\"" + currentGeno + "\")\n"
									+"Z <- get(genoName)\n"
									+ "names <- colnames(Z)\n"
									+ "for (i in 1:ncol(Z)) {\n"
									+ "    names[i] <- paste(\"chr\", names[i], sep=\"\")\n"
									+ "}\n"
									:"Z <- t(read.csv(\"" + currentGeno + "\", header=T, as.is=T, row.names=1))\n"
									+ "names <- colnames(Z)\n"
									+ "for (i in 1:ncol(Z)) {\n"
									+ "    tmp <- names[i]\n"
									+ "    if (\"_\" == substr(tmp, start=nchar(tmp)-1, stop=nchar(tmp)-1)) {\n"
									+ "        names[i] = substr(tmp, start=1, stop=nchar(tmp)-2);\n"
									+ "    }\n"
									+ "}\n"
									)
							+ "colnames(Z) <- names\n"
							+ "\n"
							+ "pheno <- read.csv(\"" + phenoFilename + "\", header=T, as.is=T, row.names=1)\n"
							+ "xphen<-na.omit(pheno)\n"
							+ "merged <- merge(xphen, Z, by=\"row.names\")\n"
							+ "mPheno <- merged[,1:ncol(pheno)+1]\n"
							+ "names <- colnames(mPheno)\n"
							+ "formu <- paste(names[1], \"~\", names[2])\n"
							+ "if (length(names)>2) {\n"
							+ "    for (i in 3:length(names)) {\n"
							+ "        formu <- paste(formu, \"+\", names[i])\n"
							+ "    }\n"
							+ "}\n"
							+ "\n"
							+ "offset <- 1+ncol(pheno)\n"
							+ "mGeno <- merged[,1:ncol(Z)+offset]\n"
							+ "\n"
							+ cohort + "_chr" + i + " <- skatCohort(Z=mGeno, formula(formu), SNPInfo=SNPInfo, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)\n"
//							+ "results <- singlesnpMeta(" + cohort + "_chr" + i + ", SNPInfo=SNPInfo, snpNames = \"Name\", cohortBetas = TRUE)\n"
							+ "results <- burdenMeta(" + cohort + "_chr" + i + ", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), SNPInfo=subset(SNPInfo, sc_nonsynSplice==TRUE), snpNames=\"SNP\", wts = 1)\n"
							+ "write.table(results, \"" + phenoRoot + "_chr" + i + "_beforeSave_results.csv\", sep=\",\", row.names = F)\n"
							+ "save(" + cohort + "_chr" + i + ", file=\"" + cohort + "_chr" + i + ".RData\")";

					filename = batchDir + cohort+ "_chr" + i + ".R";
					out = new PrintWriter(new FileOutputStream(filename));
					out.println(rCode);
					out.close();

					// then run the R code, if it only takes a few minutes then on your machine
//					Runtime.getRuntime().exec("C:/Progra~1/R/R-3.0.1/bin/Rscript " + batchDir + "chr" + i + ".R");
					v.add(filename);
					
				}
			}
			if (! foundGenos) {
				System.err.println("Error - Files not found " + genos);
			}
			if (! foundSnpInfo) {
				System.err.println("Error - Files not found " + snpInfo);
			}
			
			iterations = Matrix.toMatrix(Array.toStringArray(v));
			if (System.getProperty("os.name").startsWith("Windows")) {
				commands = "Rscript --no-save [%0]";
				Files.batchIt(batchDir + "run", "", 5, commands, iterations);
			} else {
//				commands = "/soft/R/3.0.1/bin/Rscript --no-save [%0]";
				commands = Rscript.getRscriptExecutable(new Logger())+" --no-save [%0]";
//				Files.qsub("checkObject", dir, -1, commands, iterations, 30000, 24);
				Files.qsub(batchDir + "run_" + cohort, batchDir, -1, commands, iterations, 60000, 24);
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		String cohort;
		String genos;
		String pheno;
		String snpInfo;

		cohort="aric";
		genos="D:/SkatMeta/genotypes_blacks_AA/AA_ARIC_noJHS_chr#t.csv";
		pheno="D:/SkatMeta/results_hemostasis/pheno_F7_studyIDs.csv";
		snpInfo="N:/statgen/skatMeta/fullExample/SNPInfo_HumanExome-12v1_rev5_justGeneSpliced.csv";

		String usage = "\n"+
				"gwas.SkatMetaPrimary requires 4 arguments\n"+
				"   (1) cohort name (i.e. cohort=" + cohort + " (not the default))\n"+
				"   (2) genotype file name (i.e. geno=" + genos + " (not the default))\n"+
				"   (3) phenotype file name (i.e. pheno=" + pheno + " (not the default))\n"+
				"   (4) snpInfo file name (i.e. snpInfo=" + snpInfo + " (not the default))\n"+
				"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("cohort=")) {
				cohort = args[i].split("=")[1];
			} else if (args[i].startsWith("geno=")) {
				genos = args[i].split("=")[1];
			} else if (args[i].startsWith("pheno=")) {
				pheno = args[i].split("=")[1];
			} else if (args[i].startsWith("snpInfo=")) {
				snpInfo = args[i].split("=")[1];
			} else {
				System.err.println("Error - invalid argument: "+args[i]);
			}
		}
		
		batch(cohort, genos, pheno, snpInfo);
	}

}
