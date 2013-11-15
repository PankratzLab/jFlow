package gwas;

import java.io.*;
import java.util.*;

import stats.Rscript;
import common.*;

public class SkatMetaPrimary {

	public static void batch(String cohort, String genos, String phenoFilename, String snpInfo, int qsubMem, double qsubWalltime) {
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
		String consolidate;
		Vector<String> jobNamesWithAbsolutePaths;
		IntVector jobSizes;

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
							+ "xphen <- na.omit(pheno)\n"
							+ "merged <- merge(xphen, Z, by=\"row.names\")\n"
							+ "mPheno <- merged[,1:ncol(pheno)+1]\n"
							+ "names <- colnames(pheno)\n"
							+ "if (length(names)>1) {\n"
							+ "    formu <- paste(names[1], \"~\", names[2])\n"
							+ "    for (i in 3:length(names)) {\n"
							+ "        formu <- paste(formu, \"+\", names[i])\n"
							+ "    }\n"
							+ "} else {\n"
							+ "    len <- length(mPheno)\n"
							+ "    mPheno <- c(mPheno, rep(1, len))\n"
							+ "    dim(mPheno) <- c(len, 2)\n"
							+ "    mPheno <- as.data.frame(mPheno)\n"
							+ "    colnames(mPheno) <- c(names[1], \"dummy\")\n"
							+ "    formu <- paste(names[1], \"~ 1\")\n"
							+ "}\n"
							+ "\n"
							+ "offset <- 1+ncol(pheno)\n"
							+ "mGeno <- merged[,1:ncol(Z)+offset]\n"
							+ "\n"
							+ cohort + "_chr" + i + " <- skatCohort(Z=mGeno, formula(formu), SNPInfo=SNPInfo, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)\n"
//							+ "results <- singlesnpMeta(" + cohort + "_chr" + i + ", SNPInfo=SNPInfo, snpNames = \"Name\", cohortBetas = TRUE)\n"
							+ "results <- burdenMeta(" + cohort + "_chr" + i + ", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), SNPInfo=subset(SNPInfo, sc_nonsynSplice==TRUE), snpNames=\"SNP\", wts = 1)\n"
							+ "write.table(results, \"" + cohort + "_chr" + i + "_beforeSave_results.csv\", sep=\",\", row.names = F)\n"
							+ "save(" + cohort + "_chr" + i + ", file=\"" + cohort + "_chr" + i + ".RData\", compress=\"bzip2\")";

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
//				Files.qsub("checkObject", dir, -1, commands, iterations, qsubMem, qsubWalltime);
				Files.qsub(batchDir + "run_" + cohort, batchDir, -1, commands, iterations, qsubMem, qsubWalltime);
			}

			v = new Vector<String>();
			jobNamesWithAbsolutePaths = new Vector<String>();
			jobSizes = new IntVector();
			v.add("setwd(\"" + resultDir + "\")");
			consolidate = cohort + "<- c(";
			for (int i = 0; i < iterations.length; i++) {
				v.add("load(\"" + ext.rootOf(iterations[i][0])+ ".RData\")");
				consolidate += (i==0?"":", ")+ext.rootOf(iterations[i][0]);
				jobNamesWithAbsolutePaths.add(batchDir + "run_" + cohort+"_"+(i+1)+".qsub");
				jobSizes.add((int)new File(ext.replaceAllWith(genos, "#", (i+1)+"")).length());
			}
			consolidate += ")";
			v.add(consolidate);
			v.add("class("+cohort + ") <- \"skatCohort\"");
			v.add("save(" + cohort + ", file=\"" + cohort + ".RData\", compress=\"bzip2\")");
			Files.writeList(Array.toStringArray(v), batchDir + "mergeRdataFiles.R");
			commands = Rscript.getRscriptExecutable(new Logger())+" --no-save "+batchDir + "mergeRdataFiles.R";
			Files.qsub(batchDir + "run_mergeRdataFiles_" + cohort, commands, qsubMem*4, qsubWalltime);
			Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir, "chunks", 8, true, 2000, 1);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void batchMany(String cohort, String genos, String phenosCommaDelimited, String racesCommaDelimited, String snpInfo, int qsubMem, double qsubWalltime) {
		String[] phenos, races;
		Vector<String> v;
		
		phenos = phenosCommaDelimited.split(",");
		races = racesCommaDelimited.split(",");
		
		v = new Vector<String>();
		for (int i = 0; i < phenos.length; i++) {
			for (int j = 0; j < races.length; j++) {
				try {
					batch(cohort+"_"+races[j]+"_"+phenos[i], ext.replaceAllWith(genos, "[%race]", races[j]), ext.pwd()+cohort+"_"+races[j]+"_"+phenos[i]+".csv", snpInfo, qsubMem, qsubWalltime);
				} catch (Exception e) {
					System.err.println("Error - failed to script up "+phenos[i]+"/"+races[j]);
				}
				v.add("cd "+cohort+"_"+races[j]+"_"+phenos[i]+"/batchFiles/");
				v.add("./master.run_"+cohort+"_"+races[j]+"_"+phenos[i]);
				v.add("cd ../../");
				v.add("");
			}
		}
		Files.writeList(Array.toStringArray(v), "scriptAll");

		v = new Vector<String>();
		for (int i = 0; i < phenos.length; i++) {
			for (int j = 0; j < races.length; j++) {
				v.add("cd "+cohort+"_"+races[j]+"_"+phenos[i]+"/batchFiles/");
				v.add(Rscript.getRscriptExecutable(new Logger())+" --no-save mergeRdataFiles.R");
				v.add("mv ../results/"+cohort+"_"+races[j]+"_"+phenos[i]+".RData ../../");
				v.add("cd ../../");
				v.add("");
			}
		}
		Files.writeList(Array.toStringArray(v), "mergeAll");
	}

	public static void main(String[] args) {
		String cohort;
		String genos;
		String pheno;
		String snpInfo;
		String phenos = null;
		String races = "EA,AA";
		int qsubMem;
		double qsubWalltime;

		cohort="aric";
		genos="D:/SkatMeta/genotypes_blacks_AA/AA_ARIC_noJHS_chr#t.csv";
		pheno="D:/SkatMeta/results_hemostasis/pheno_F7_studyIDs.csv";
		snpInfo="N:/statgen/skatMeta/fullExample/SNPInfo_HumanExome-12v1_rev5_justGeneSpliced.csv";
		qsubMem = 15000;
		qsubWalltime = 2;

		String usage = "\n"+
				"gwas.SkatMetaPrimary requires 4 arguments\n"+
				"   (1) cohort name (i.e. cohort=" + cohort + " (not the default))\n"+
				"   (2) genotype file name (i.e. geno=" + genos + " (not the default))\n"+
				"   (3) phenotype file name (i.e. pheno=" + pheno + " (not the default))\n"+
				"   (4) snpInfo file name (i.e. snpInfo=" + snpInfo + " (not the default))\n"+
				"   (5) qsub memory in megabytes (i.e. qsubmem=" + qsubMem + " (default))\n"+
				"   (6) qsub walltime in hours (i.e. qsubwalltime=" + qsubWalltime + " (default))\n"+
				"   (7) (optional) phenotypes to run (i.e. phenos=CRP,ICAM1,TNFA (not the default))\n"+
				"   (8) (optional) races to run (i.e. races="+races+" (default))\n"+
				" ( the final optional parameters require a fixed phenotype format of [cohort name]_[race]_[phenotype]_.csv )\n"+
				" ( as well as a genotype file with the format of regular_filename_[%race]_chr#.csv )\n"+
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
			} else if (args[i].startsWith("qsubmem=")) {
				qsubMem = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("qsubwalltime=")) {
				qsubWalltime = Double.parseDouble(args[i].split("=")[1]);
			} else if (args[i].startsWith("phenos=")) {
				phenos = args[i].split("=")[1];
			} else if (args[i].startsWith("races=")) {
				races = args[i].split("=")[1];
			} else {
				System.err.println("Error - invalid argument: "+args[i]);
			}
		}
		
		if (phenos != null) {
			batchMany(cohort, genos, phenos, races, snpInfo, qsubMem, qsubWalltime);
		} else {
			batch(cohort, genos, pheno, snpInfo, qsubMem, qsubWalltime);
		}
	}

}
