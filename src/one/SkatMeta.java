package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import parse.LookupTable;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.Sort;
import common.ext;

public class SkatMeta {
	public static final String[] GENE_RESULT_COLUMNS = new String[] {"gene", "p"};
	public static final String[] SNPINFO_COLUMNS = new String[] {"Name", "Chr", "MapInfo"};
	public static final String FILENAME_ETHNIC_SEGMENT = "%ethnic%";
	public static final String FILENAME_PHENO_SEGMENT = "%pheno%";
	public static final String FILENAME_CONDITION_SEGMENT = "%cond%";
	public static final String FILENAME_CHROMOSOME_SEGMENT = "%chr%";
	public static final String FILENAME_ANALYSIS_SEGMENT = "%analysis%";
	public static final String[] OPERATORS = new String[] {"<=", "==", ">=", "!=", "<", ">"};	//NOTE: For current parser, longer operators should be placed before shorter ones.
	public static final String[][] SPECIAL_CHROMOSOMES = new String[][] {{"X", "Y", "XY", "MT"}, {"23", "24", "25", "26"}};


	/* This is a working copy before the 2nd draft of the method with the same name
	public static void generateSkatMetaRScriptConditional(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
		String[] files;
		PrintWriter writer;
		String line;
		String[] temp1, temp2;
		String phenoAndCondition;
		String condFile;
		byte count;

		files = Files.list(sourceRDataFilesDir, ".RData", false);
		try {
			phenoAndCondition = "";
			temp1 = sourceRDataFilesDir.split("/");
			for (int i = 0; i < 3; i++) {
				phenoAndCondition = temp1[temp1.length - i - 1] + "_" + phenoAndCondition;
			}

			phenoAndCondition = phenoAndCondition.substring(0, phenoAndCondition.length()-1);
			writer = new PrintWriter(new FileOutputStream(rScriptDir + phenoAndCondition + ".R"));
			writer.println("library(seqMeta)\n"
							+ "temp <- load(\"" + snpInfoFile + "\")\n"
							+ "SNPInfo <- get(temp)\n"
							+ "rm(list=temp)\n"
							+ "rm(temp)\n"
//							+ "names(SNPInfo) <- c(\"Name\", names(SNPInfo)[2:length(SNPInfo)])\n"
							+ "\nsetwd(\"" + sourceRDataFilesDir + "\")");
			line = new String();
			for (int i = 0; i < files.length; i++) {
				writer.println(   "temp <- load(\"" + files[i] + "\")\n"
								+ "Cohort" + i + " <- get(temp)\n"
								+ "rm(list=temp)\n"
								+ "rm(temp)");
				line += "Cohort" + i + ", ";
			}

			condFile = "";
			files = Files.list(condFileDir, ".txt", false);
			for (int i = 0; i < files.length; i++) {
				count = 0;
				temp2 = files[i].split("\\.")[0].split("_");
				for (int j = 0; j < temp2.length; j++) {
					for (int k = 0; k < 3; k++) {
						if(temp2[j].equalsIgnoreCase(temp1[temp1.length - k - 1])) {
							count ++;
							break;
						}
					}
				}
				if(count == (temp2.length - (temp2[0].equalsIgnoreCase("f8")?1:0))) {
					condFile = files[i];
					break;
				}
			}
			writer.println("\ntemp <- read.table(\"" + condFileDir + condFile + "\", header=TRUE)\n"
					+ "genes <- SNPInfo$SKATgene %in% temp$SKATgene\n"
					+ "results <- singlesnpMeta(" + line + "SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE)\n"
					+ "write.table(results, \"" + resultsDir + phenoAndCondition + "_SingleSNP.csv\", sep=\",\", row.names = F)\n\n"
					+ "results <- burdenMeta(" + line + "SNPInfo=subset(SNPInfo[genes,], sc_functional==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1)\n"
					+ "write.table(results, \"" + resultsDir + phenoAndCondition + "_T5Count.csv\", sep=\",\", row.names = F)\n\n"
					+ "results <- skatMeta(" + line + "SNPInfo=subset(SNPInfo[genes,], sc_functional==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)})\n"
					+ "write.table( results, \"" + resultsDir + phenoAndCondition + "_SKAT_T5.csv\", sep=\",\", row.names = F)\n");
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	*/

	public static void generateSkatMetaRScriptConditional(String sourceRDataFilesDir, String snpInfoFile, String condFileDirFilenameTemplate, String rScriptDir, String resultsDir) {
		PrintWriter writer;
		String[] temp1;
		String[] phenoCondStratum;
		String condFile;
		String root;
		
		Files.list(sourceRDataFilesDir, ".RData", false);
		phenoCondStratum = new String[3];
		try {
			temp1 = sourceRDataFilesDir.split("/");
			for (int i = 0; i < 3; i++) {
				phenoCondStratum[i] = temp1[temp1.length - 3 + i];
			}
			
			condFile = condFileDirFilenameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenoCondStratum[0]).replaceAll(FILENAME_CONDITION_SEGMENT, phenoCondStratum[1]);

			root = Array.toStr(phenoCondStratum, "_");
			writer = new PrintWriter(new FileOutputStream(rScriptDir + root + ".R"));
			writer.println(getRScriptForConditional(sourceRDataFilesDir, snpInfoFile, condFile, rScriptDir, resultsDir));
			writer.close();

			writer = new PrintWriter(new FileOutputStream(rScriptDir + "runAll.bat", true));
			writer.println("RScript --no-save "+root + ".R > "+root + ".log");
			writer.close();
			
			if (!Files.exists(rScriptDir+"splitRuns.crf")) {
				writer = new PrintWriter(new FileOutputStream(rScriptDir+"splitRuns.crf"));
				writer.println("split");
				writer.println("runAll.bat numFiles=4 sizeOfHeader=0 blockSize=1 root=run ext=.bat");
				writer.close();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public static void summaryTable(String[] phenos, String[] races, String condResultsDir, String originalResultsDir, String summaryDir, String snpInfoFile, Logger log) {
		PrintWriter writer;
		int numConds;
		String minPstatement;
		
		summaryDir = ext.verifyDirFormat(summaryDir);
		if (!Files.exists(summaryDir)) {
			new File(summaryDir).mkdirs();
		}

		for (int i = 0; i < phenos.length; i++) {
			numConds = 1;
			while (Files.exists(condResultsDir + phenos[i]+"_cond"+numConds+"_"+races[0]+"_SingleSNP.csv")) {
				numConds++;
			}
			numConds -= 1;
			
			if (numConds == 0) {
				log.reportError("Error - zero conditions detected for "+ phenos[i]+"/"+races[0]+"_in "+condResultsDir);
				continue;
			}
			
			// If anything else is added to each row, then make sure to change the +7 to conform to the first index of a p-value
			minPstatement="=MIN(";
			for (int j = 0; j < races.length*(numConds+1); j++) {
				minPstatement += (j==0?"":",")+ext.getExcelColumn(j*2+7)+"1";
			}
			minPstatement += ")";
			
			try {
				writer = new PrintWriter(new FileOutputStream(summaryDir + "parse"+phenos[i]+".crf"));
				writer.println("lookup");
				writer.println(condResultsDir + phenos[i]+"_cond1_"+races[0]+"_SingleSNP.csv 1 header out="+summaryDir+phenos[i]+"_conditionalResults.xln head=MarkerName");
				
				if (Files.exists(ext.rootOf(snpInfoFile, false)+".csv")) {
					writer.println(ext.rootOf(snpInfoFile, false)+".csv 'SNP' 'SKATgene'=Gene 'func_region'=Function $@"+minPstatement);
				} else {
					log.report("Warning - could not find text dump of snpInfo file at '"+ext.rootOf(snpInfoFile, false)+".csv"+"'");
				}

				for (int k = 1; k < races.length; k++) {
					writer.println(originalResultsDir + phenos[i]+"_"+races[k]+"_SingleSNP.csv 'Name' 'maf'="+races[k]+"_maf");
				}

				for (int k = 0; k < races.length; k++) {
					writer.println(originalResultsDir + phenos[i]+"_"+races[k]+"_SingleSNP.csv 'Name' 'beta'=beta_uncond_"+races[k]+" 'p'=p_uncond_"+races[k]);
					for (int cond = 1; cond <= numConds; cond++) {
						writer.println(condResultsDir + phenos[i]+"_cond"+cond+"_"+races[k]+"_SingleSNP.csv 'Name' 'beta'=beta_cond"+cond+"_"+races[k]+" 'p'=p_cond"+cond+"_"+races[k]);
					}
				}
				writer.close();

				log.report("Building results for "+phenos[i]);
				LookupTable.fromParameters(summaryDir + "parse"+phenos[i]+".crf", log);
				
			} catch (Exception e) {
				log.reportError("Error writing to " + "parse"+phenos[i]+".crf");
				log.reportException(e);
			}
		}
	}

	public static void generateSkatMetaRScript(String phenoDirAndNameTemplate, String genoDirAndNameTemplate, String condFileDirAndNameTemplate, String snpInfoDirAndNameTemplate, String[] phenos, String[] ethnics, String startCondition, String rScriptDir, String resultsDir, String rcommand) {
		String condFile, phenoFile, genoFile, snpInfoFile = null, outputFile, rScript;
		Hashtable<String, String[]> qsubIterations;
//		String[][] qsubIterationsMatrix;
		String[] chrs, tmp;
		String[][] tmp2;
		qsubIterations = new Hashtable<String, String[]>();
		for (int i = 0; i < phenos.length; i++) {
			condFile = condFileDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
			condFile = condFile.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + startCondition);
			tmp2 = new String[ethnics.length][1];
			rScript = "";
			for (int j = 0; j < ethnics.length; j++) {
				condFile = condFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
				chrs = loadFile(condFile, null, new String[] {"CHR"}, null, null, null).keySet().toArray(new String[0]);
				if (chrs != null && chrs.length > 0) {
					qsubIterations.put(phenos[i] + "\t" + ethnics[j] + "\t" + startCondition, chrs);
					tmp = new String[chrs.length];
					for (int k = 0; k < chrs.length; k++) {
						snpInfoFile = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chrs[k]);
						genoFile = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chrs[k]);
						genoFile = genoFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
						phenoFile = phenoDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
						phenoFile = phenoFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
						outputFile = resultsDir + "ARIC_" + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_chr" + chrs[k] + "_prepCondScores_seqMeta.RData";
						tmp[k] = outputFile;
//						tmp[k] = rScriptDir + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_cohort_chr" + chrs[?] + ".R";
//						Files.write(getRScriptForCohort(snpInfoFile, genoFile, phenoFile, condFile, outputFile), tmp[l]);
						rScript += (getRScriptForCohort(snpInfoFile, genoFile, phenoFile, condFile, outputFile) + "\n");
					}
					tmp2[j][0] = resultsDir + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_cohort_merged.RData";
//					Files.write(getRScriptForMergingChromosomes(tmp, tmp2[k][0]), resultsDir + phenos[i] + "_" + ethnics[k] + "_cond" + startCondition + "_cohort_merge.R");
					rScript += (getRScriptForMergingChromosomes(tmp, tmp2[j][0]) + "\n");
				}
			}
//			Files.write(getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFileDirAndNameTemplate, rScriptDir + phenos[i] + "_cond" + startCondition + "_", null), rScriptDir + phenos[i] + "_cond" + startCondition + "_meta.R");
			rScript += (getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFile, resultsDir + phenos[i] + "_cond" + startCondition, null));
			Files.write(rScript, rScriptDir + phenos[i] + "_cond" + startCondition + "_all.R");
			Files.qsub(rScriptDir + "[%0]_cond[%1]_all", rcommand + " " + rScriptDir + "[%0]_cond[%1]_all.R", new String[][] {{phenos[i], startCondition}}, -1, 6, -1);
			new File(rScriptDir + "master.qsub").renameTo(new File(rScriptDir + phenos[i] + "_cond" + startCondition + "_all.sh"));
		}
	}

	/* This is a working copy before the 2nd draft of the method with the same name
	public static void generateSkatMetaRScript(String phenoDirAndNameTemplate, String genoDirAndNameTemplate, String condFileDirAndNameTemplate, String snpInfoDirAndNameTemplate, String[] phenos, String[] ethnics, String startCondition, String rScriptDir, String resultsDir, String rcommand) {
		String condFile, phenoFile, genoFile, snpInfoFile = null, outputFile, rScript;
		Hashtable<String, String[]> qsubIterations;
		String[][] qsubIterationsMatrix;
		String[] chrs, tmp1;
		String[][] tmp2;
		String tmp3;
		int counter;

		counter = 0;
		qsubIterations = new Hashtable<String, String[]>();
		qsubIterationsMatrix = new String[phenos.length][2];
		for (int i = 0; i < phenos.length; i++) {
			condFile = condFileDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
			condFile = condFile.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + startCondition);
			tmp2 = new String[ethnics.length][1];
			rScript = "";
			for (int j = 0; j < ethnics.length; j++) {
				condFile = condFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
				chrs = loadFile(condFile, null, new String[] {"CHR"}, null, null, null).keySet().toArray(new String[0]);
				if (chrs != null && chrs.length > 0) {
					qsubIterations.put(phenos[i] + "\t" + ethnics[j] + "\t" + startCondition, chrs);
					tmp1 = new String[chrs.length];
					for (int k = 0; k < chrs.length; k++) {
						snpInfoFile = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chrs[k]);
						genoFile = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chrs[k]);
						genoFile = genoFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
						phenoFile = phenoDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
						phenoFile = phenoFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
						outputFile = resultsDir + "ARIC_" + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_chr" + chrs[k] + "_prepCondScores_seqMeta.RData";
						tmp1[k] = outputFile;
//						tmp[k] = rScriptDir + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_cohort_chr" + chrs[?] + ".R";
//						Files.write(getRScriptForCohort(snpInfoFile, genoFile, phenoFile, condFile, outputFile), tmp[l]);
						rScript += (getRScriptForCohort(snpInfoFile, genoFile, phenoFile, condFile, outputFile) + "\n");
						counter ++;
					}
					tmp2[j][0] = resultsDir + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_cohort_merged.RData";
//					Files.write(getRScriptForMergingChromosomes(tmp, tmp2[k][0]), resultsDir + phenos[i] + "_" + ethnics[k] + "_cond" + startCondition + "_cohort_merge.R");
					rScript += (getRScriptForMergingChromosomes(tmp1, tmp2[j][0]) + "\n");
				}
			}
//			Files.write(getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFileDirAndNameTemplate, rScriptDir + phenos[i] + "_cond" + startCondition + "_", null), rScriptDir + phenos[i] + "_cond" + startCondition + "_meta.R");
			rScript += (getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFile, resultsDir + phenos[i] + "_cond" + startCondition, null));
			Files.write(rScript, rScriptDir + phenos[i] + "_cond" + startCondition + "_all.R");
			qsubIterationsMatrix[i] = new String[] {phenos[i], startCondition};
		}
		Files.qsub(rScriptDir + "[%0]_cond[%1]_all", rcommand + " " + rScriptDir + "[%0]_cond[%1]_all.R", qsubIterationsMatrix, -1, 6, -1);
		tmp3 = "";
		for (int i = 0; i < phenos.length; i++) {
			tmp3 += ("_" + phenos[i]);
		}
		new File(rScriptDir + "master.qsub").renameTo(new File(rScriptDir + "master" + tmp3 + "_cond" + startCondition + ".sh"));

//		if (qsubIterations.size() > 0) {
////			if (new File(rScriptDir, "master.qsub").exists()) {
////				new File(rScriptDir, "master.qsub").renameTo(new File(rScriptDir, "master" + Files.list(rScriptDir, "master", ".qsub", false, false).length + ".qsub"));
////			}
//			if (rcommand == null || rcommand.equals("")) {
//				rcommand = "Rscript";
//			}
//
//			qsubIterationsMatrix = new String[counter][4];
//			counter = 0;
//			for (String phenoEthnicCondition : qsubIterations.keySet()) {
//				tmp = phenoEthnicCondition.split("\t");
//				chrs = qsubIterations.get(phenoEthnicCondition);
//				for (int i = 0; i < chrs.length; i++) {
//					for (int j = 0; j < tmp.length; j++) {
//						qsubIterationsMatrix[counter][j] = tmp[j];
//					}
//					qsubIterationsMatrix[counter][3] = chrs[i];
//					counter ++;
//				}
//			}
//			Files.qsub(rScriptDir + "[%0]_[%1]_cond[%2]_cohort_chr[%3]", rcommand + " " + rScriptDir + "[%0]_[%1]_cond[%2]_cohort_chr[%3].R", qsubIterationsMatrix, -1, 6, -1);
//			outputFile = rScriptDir + "masterqsub";
//			for (int i = 0; i < phenos.length; i++) {
//				outputFile += ("_" + phenos[i]);
//			}
//			new File(rScriptDir + "master.qsub").renameTo(new File(outputFile + "_cohort.sh"));
//
//			qsubIterationsMatrix = new String[qsubIterations.size()][];
//			counter = 0;
//			for (String phenoEthnicCondition : qsubIterations.keySet()) {
//				qsubIterationsMatrix[counter] = phenoEthnicCondition.split("\t");
//				counter ++;
//			}
//			Files.qsub(rScriptDir + "[%0]_[%1]_cond[%2]_cohort_merge", rcommand + " " + rScriptDir + "[%0]_[%1]_cond[%2]_cohort_merge.R", qsubIterationsMatrix, -1, 6, -1);
//			new File(rScriptDir + "master.qsub").renameTo(new File(outputFile + "_cohort_merge.sh"));
//
//			
//		}
	}
	*/

	public static String getRScriptForConditional(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
		String[] files;
		String rscript;
		String line;
		String[] temp1, temp2;
		String condFile;
		byte count;
		String phenoAndCondition;

		phenoAndCondition = "";
		temp1 = sourceRDataFilesDir.split("/");
		for (int i = 0; i < 3; i++) {
			phenoAndCondition = temp1[temp1.length - i - 1] + "_" + phenoAndCondition;
		}

		files = Files.list(sourceRDataFilesDir, ".RData", false);
		rscript = ("library(seqMeta)\n"
						+ "library(\"methods\")\n"
						+ "temp <- load(\"" + snpInfoFile + "\")\n"
						+ "SNPInfo <- get(temp)\n"
						+ "rm(list=temp)\n"
						+ "rm(temp)\n"
//							+ "names(SNPInfo) <- c(\"Name\", names(SNPInfo)[2:length(SNPInfo)])\n"
						+ "\nsetwd(\"" + sourceRDataFilesDir + "\")\n");
		line = new String();
		for (int i = 0; i < files.length; i++) {
			rscript += (   "temp <- load(\"" + files[i] + "\")\n"
							+ "Cohort" + i + " <- get(temp)\n"
							+ "rm(list=temp)\n"
							+ "rm(temp)\n");
			line += "Cohort" + i + ", ";
		}

		condFile = "";
		files = Files.list(condFileDir, ".txt", false);
		for (int i = 0; i < files.length; i++) {
			count = 0;
			temp2 = files[i].split("\\.")[0].split("_");
			for (int j = 0; j < temp2.length; j++) {
				for (int k = 0; k < 3; k++) {
					if(temp2[j].equalsIgnoreCase(temp1[temp1.length - k - 1])) {
						count ++;
						break;
					}
				}
			}
			if(count == (temp2.length - (temp2[0].equalsIgnoreCase("f8")?1:0))) {
				condFile = files[i];
				break;
			}
		}
		rscript += ("\ntemp <- read.table(\"" + condFileDir + condFile + "\", header=TRUE)\n"
				+ "genes <- SNPInfo$SKATgene %in% temp$SKATgene\n"
				+ "results <- singlesnpMeta(" + line + "SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE)\n"
				+ "write.table(results, \"" + resultsDir + phenoAndCondition + "_SingleSNP.csv\", sep=\",\", row.names = F)\n\n"
				+ "results <- burdenMeta(" + line + "SNPInfo=subset(SNPInfo[genes,], sc_functional==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1)\n"
				+ "write.table(results, \"" + resultsDir + phenoAndCondition + "_T5Count.csv\", sep=\",\", row.names = F)\n\n"
				+ "results <- skatMeta(" + line + "SNPInfo=subset(SNPInfo[genes,], sc_functional==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)})\n"
				+ "write.table( results, \"" + resultsDir + phenoAndCondition + "_SKAT_T5.csv\", sep=\",\", row.names = F)\n");
		return rscript;
	}

	public static String getRScriptForCohort (String snpInfoFile, String genoFile, String phenoFile, String condFile, String resultFileName) {
		return ("library(seqMeta)\nlibrary(\"methods\")\n"
					+ "\ntmp <- load(\"" + snpInfoFile + "\")\nSNPInfo <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n"
					+ (genoFile.endsWith(".RData")?
							  "\ntmp <- load(\"" + genoFile + "\")\nZ <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n"
							: "\nZ <- t(read.csv(\"" + genoFile + "\", header=T, as.is=T, row.names=1));\n")
					+ "names <- colnames(Z);\n"
					+ "for (i in 1:ncol(Z)) {\n"
					+ "	tmp <- names[i];\n"
					+ "	if (\"_\" == substr(tmp, start=nchar(tmp)-1, stop=nchar(tmp)-1)) {\n"
					+ "		names[i] = substr(tmp, start=1, stop=nchar(tmp)-2);\n"
					+ "	}\n"
					+ "}\n"
					+ "colnames(Z) <- names;\n"
					+ "\npheno <- na.omit(read.csv(\"" + phenoFile + "\", header=T, as.is=T, row.names=1));\n"
					+ "pheno_Z_merged <- merge(pheno, Z, by=\"row.names\")\n"
					+ "mPheno <- pheno_Z_merged[, (1 : ncol(pheno)) + 1]\n"
					+ "names <- colnames(pheno)\n"
					+ "if (length(names)>1) {\n"
					+ "	formu <- paste(names[1], \"~\", names[2])\n"
					+ "	for (i in 3:length(names)) {\n"
					+ "		formu <- paste(formu, \"+\", names[i])\n"
					+ "	}\n"
					+ "} else {\n"
					+ "	len <- length(mPheno)\n"
					+ "	mPheno <- c(mPheno, rep(1, len))\n"
					+ "	dim(mPheno) <- c(len, 2)\n"
					+ "	mPheno <- as.data.frame(mPheno)\n"
					+ "	colnames(mPheno) <- c(names[1], \"dummy\")\n"
					+ "	formu <- paste(names[1], \"~ 1\")\n"
					+ "}\n"
					+ "\nmGeno <- pheno_Z_merged[,(1 : ncol(Z)) + (1 + ncol(pheno))]\n"
					+ "\ncond <- read.table(\"" + condFile + "\", header=TRUE)\n"
//					+ "genes <- SNPInfo$SKATgene %in% cond$SKATgene\n"
					+ "result <- prepCondScores(Z=mGeno, formula(formu), SNPInfo=SNPInfo, adjustments=cond, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)\n"
					+ "save(result, file=\"" + resultFileName + "\", compress=\"bzip2\");\n"
					);
	}

	public static String getRScriptForCohort1 (String snpInfoFile, String genoFile, String phenoFile, String condFile, String resultFileName) {
		String resultId;
		resultId = ext.rootOf(resultFileName);

		return ("library(seqMeta)\nlibrary(\"methods\")\n"
					+ "\ntmp <- load(\"" + snpInfoFile + "\")\nSNPInfo <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n"
					+ (genoFile.endsWith(".RData")?
							  "\ntmp <- load(\"" + genoFile + "\")\nZ <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n"
							: "\nZ <- t(read.csv(\"" + genoFile + "\", header=T, as.is=T, row.names=1));\n")
					+ "names <- colnames(Z);\n"
					+ "for (i in 1:ncol(Z)) {\n"
					+ "	tmp <- names[i];\n"
					+ "	if (\"_\" == substr(tmp, start=nchar(tmp)-1, stop=nchar(tmp)-1)) {\n"
					+ "		names[i] = substr(tmp, start=1, stop=nchar(tmp)-2);\n"
					+ "	}\n"
					+ "}\n"
					+ "colnames(Z) <- names;\n"
					+ "\npheno <- na.omit(read.csv(\"" + phenoFile + "\", header=T, as.is=T, row.names=1));\n"
					+ "pheno_Z_merged <- merge(pheno, Z, by=\"row.names\")\n"
					+ "mPheno <- pheno_Z_merged[, (1 : ncol(pheno)) + 1]\n"
					+ "names <- colnames(pheno)\n"
					+ "if (length(names)>1) {\n"
					+ "	formu <- paste(names[1], \"~\", names[2])\n"
					+ "	for (i in 3:length(names)) {\n"
					+ "		formu <- paste(formu, \"+\", names[i])\n"
					+ "	}\n"
					+ "} else {\n"
					+ "	len <- length(mPheno)\n"
					+ "	mPheno <- c(mPheno, rep(1, len))\n"
					+ "	dim(mPheno) <- c(len, 2)\n"
					+ "	mPheno <- as.data.frame(mPheno)\n"
					+ "	colnames(mPheno) <- c(names[1], \"dummy\")\n"
					+ "	formu <- paste(names[1], \"~ 1\")\n"
					+ "}\n"
					+ "\nmGeno <- pheno_Z_merged[,(1 : ncol(Z)) + (1 + ncol(pheno))]\n"
					+ "\ncond <- read.table(\"" + condFile + "\", header=TRUE)\n"
//					+ "genes <- SNPInfo$SKATgene %in% cond$SKATgene\n"
					+ resultId + " <- prepCondScores(Z=mGeno, formula(formu), SNPInfo=SNPInfo, adjustments=cond, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)\n"
					+ "save(" + resultId + ", file=\"" + resultFileName + "\", compress=\"bzip2\");\n"
					);
	}

	public static String getRScriptForMergingChromosomes(String[] filenames, String outputFileName) {
		String rScript;

		rScript = "";
		for (int i = 0; i < filenames.length; i++) {
			rScript += "load(\"" + filenames[i] + "\");\ntmp" + i + " <- result;\nrm(result);\n\n";
		}
		rScript += "combined <- c(";
		for (int i = 0; i < filenames.length; i++) {
			rScript += ("tmp" + i + (i == (filenames.length - 1)? "" : ","));
		}
		rScript += (");\n"
					+ "class(combined) <- \"seqMeta\"\n"
					+ "save(combined, file=\"" + outputFileName + "\", compress=\"bzip2\");"
					);

		return rScript;
	}

	public static String getRScriptForMetaAnalysis (String snpInfoFile, String[] ethnics, String[][] cohortResultFilesnames, String conditionFile, String outputDirAndRoot, Logger log) {
		String rScript;
		String[] allDataByEthnic;
		String allEthnics, allData;

		if (log == null) {
			log = new Logger();
		}
		if (ethnics.length != cohortResultFilesnames.length) {
			log.reportError("Error - array 'ethnics' size (" + ethnics.length + ") is different from array 'filenames' size (" + cohortResultFilesnames.length + ")");
			System.exit(0);
		}

		allDataByEthnic = new String[ethnics.length];
		rScript = "library(seqMeta)\nlibrary(\"methods\")\n\ntmp <- load(\"" + snpInfoFile + "\")\nSNPInfo <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n\n";
		for (int i = 0; i < cohortResultFilesnames.length; i++) {
			allDataByEthnic[i] = "";
			for (int j = 0; j < cohortResultFilesnames[i].length; j++) {
				if (cohortResultFilesnames[i][j] != null) {
					allDataByEthnic[i] += (ethnics[i] + j + ",");
					rScript += "tmp <- load(\"" + cohortResultFilesnames[i][j] + "\");\n" + ethnics[i] + j + " <- get(tmp);\nrm(list=tmp);\nrm(tmp);\n\n";
				}
			}
		}
		rScript += "temp <- read.table(\"" + conditionFile + "\", header=TRUE);\ngenes <- SNPInfo$gene %in% temp$SKATgene;\n\n";

		for (int i = 0; i < ethnics.length; i++) {
			rScript += ("result <- singlesnpMeta(" + allDataByEthnic[i] + " SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE);\n"
						+ "write.table(result, \"" + outputDirAndRoot + ethnics[i] + "_SingleSNP.csv\", sep=\",\", row.names = F);\n"
						+ "result <- burdenMeta(" + allDataByEthnic[i] + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1);\n"
						+ "write.table(result, \"" + outputDirAndRoot + ethnics[i] + "_T5Count.csv\", sep=\",\", row.names = F);\n"
						+ "result <- skatMeta(" + allDataByEthnic[i] + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)});\n"
						+ "write.table(result, \"" + outputDirAndRoot + ethnics[i] + "_T5.csv\", sep=\",\", row.names = F);\n\n"
						);
		}

		allEthnics = "";
		allData = "";
		for (int i = 0; i < ethnics.length; i++) {
			allData += allDataByEthnic[i];
			allEthnics += ethnics[i];
		}
		rScript += ("result <- singlesnpMeta(" + allData + " SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE);\n"
					+ "write.table(result, \"" + outputDirAndRoot + allEthnics + "_SingleSNP.csv\", sep=\",\", row.names = F);\n"
					+ "result <- burdenMeta(" + allData + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1);\n"
					+ "write.table(result, \"" + outputDirAndRoot + allEthnics + "_T5Count.csv\", sep=\",\", row.names = F);\n"
					+ "result <- skatMeta(" + allData + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)});\n"
					+ "write.table(result, \"" + outputDirAndRoot + allEthnics + "_T5.csv\", sep=\",\", row.names = F);\n"
					);

		return rScript;
	}

	public static String getRScriptForMetaAnalysis1 (String snpInfoFile, String[] ethnics, String[][] cohortResultFilesnames, String conditionFile, String outputDirFilenameTemplate, Logger log) {
		String rScript;
		String[] allDataByEthnic;
		String allEthnics, allData;

		if (log == null) {
			log = new Logger();
		}
		if (ethnics.length != cohortResultFilesnames.length) {
			log.reportError("Error - array 'ethnics' size (" + ethnics.length + ") is different from array 'filenames' size (" + cohortResultFilesnames.length + ")");
			System.exit(0);
		}

		allDataByEthnic = new String[ethnics.length];
		rScript = "library(seqMeta)\nlibrary(\"methods\")\n\ntmp <- load(\"" + snpInfoFile + "\")\nSNPInfo <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n\n";
		for (int i = 0; i < cohortResultFilesnames.length; i++) {
			allDataByEthnic[i] = "";
			for (int j = 0; j < cohortResultFilesnames[i].length; j++) {
				if (cohortResultFilesnames[i][j] != null) {
					allDataByEthnic[i] += (ethnics[i] + j + ",");
					rScript += "tmp <- load(\"" + cohortResultFilesnames[i][j] + "\");\n" + ethnics[i] + j + " <- get(tmp);\nrm(list=tmp);\nrm(tmp);\n\n";
				}
			}
		}
		rScript += "cond <- read.table(\"" + conditionFile + "\", header=TRUE);\ngenes <- SNPInfo$SKATgene %in% cond$SKATgene;\n\n";

		for (int i = 0; i < ethnics.length; i++) {
			if (! allDataByEthnic[i].equals("")) {
				rScript += ("result <- singlesnpMeta(" + allDataByEthnic[i] + " SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE);\n"
							+ "write.table(result, \"" + outputDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]).replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP") + "\", sep=\",\", row.names = F);\n"
							+ "result <- burdenMeta(" + allDataByEthnic[i] + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1);\n"
							+ "write.table(result, \"" + outputDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]).replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5Count") + "\", sep=\",\", row.names = F);\n"
							+ "result <- skatMeta(" + allDataByEthnic[i] + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)});\n"
							+ "write.table(result, \"" + outputDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]).replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5") + "\", sep=\",\", row.names = F);\n\n"
							);
			}
		}

		allEthnics = "";
		allData = "";
		for (int i = 0; i < ethnics.length; i++) {
			allData += allDataByEthnic[i];
			allEthnics += ethnics[i];
		}
		rScript += ("result <- singlesnpMeta(" + allData + " SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE);\n"
					+ "write.table(result, \"" + outputDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, allEthnics).replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP") + "\", sep=\",\", row.names = F);\n"
					+ "result <- burdenMeta(" + allData + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1);\n"
					+ "write.table(result, \"" + outputDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, allEthnics).replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5Count") + "\", sep=\",\", row.names = F);\n"
					+ "result <- skatMeta(" + allData + " SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)});\n"
					+ "write.table(result, \"" + outputDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, allEthnics).replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5") + "\", sep=\",\", row.names = F);\n"
					);

		return rScript;
	}

	public static void generateSkatMetaRScriptsForSubFolders(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
		String[] folders;
		
		if (!Files.exists(rScriptDir)) {
			new File(rScriptDir).mkdirs();
		}

		folders = Files.listDirectories(sourceRDataFilesDir, false);
		if (ext.indexOfStr("EAAA", folders) != -1) {
			folders = Sort.putInOrder(folders, Sort.quicksort(folders, Sort.DESCENDING));
		}
		if (folders == null || folders.length < 1) {
			generateSkatMetaRScriptConditional(sourceRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
		} else {
			for (int i = 0; i < folders.length; i++) {
				generateSkatMetaRScriptsForSubFolders(sourceRDataFilesDir + folders[i] + "/", snpInfoFile, condFileDir, rScriptDir, resultsDir);
			}
		}
	}

	public static void qcScript(String sourceRDataFilesDir) {
		String[] files;
		PrintWriter writer;

		files = listFilesInDirAndAllSubDirs(sourceRDataFilesDir);
		if (files != null && files.length > 0) {
			try {
				writer = new PrintWriter(new FileOutputStream(sourceRDataFilesDir  + "qc.R"));
				writer.println("result <- c(\"\", \"\")\n");
				for (int i = 0; i < files.length; i++) {
//					writer.println("file <- \"" + files[i] + "\"\ntemp <- load(file)\ntmp <- get(temp)\nrm(list=temp)\nrm(temp)\nif(length(ls(tmp))>7) result <- c(result, c(file, length(ls(tmp))))\n");
					writer.println("file <- \"" + files[i] + "\"\ntemp <- load(file)\ntmp <- get(temp)\nrm(list=temp)\nrm(temp)\nresult <- c(result, c(file, length(ls(tmp))))\n");
				}
				writer.println("result <- matrix(result, ncol=2, byrow = TRUE)\nwrite.table(result, \"" + sourceRDataFilesDir + "qc.csv\", sep=\",\")");
				writer.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	public static String getRScriptForConditionalAnalysis (String condFileDirAndNameTemplate, int currentCondition, String phenoDirAndNameTemplate, String genoDirAndNameTemplate, String snpInfoDirAndNameTemplate, String pheno, String[] ethnics, String rScriptDir, String resultsDir, String rcommand, String fullpathToOutputNextRoundCondition, double pThreshold, Logger log) {
			String condFile, phenoFile, genoFile, snpInfoFile = null, outputFile, rScript;
			String[] chrs, tmp;
			String[][] tmp2;
	
			tmp2 = new String[ethnics.length][1];
			rScript = "";
			for (int j = 0; j < ethnics.length; j++) {
				condFile = condFileDirAndNameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
				chrs = loadFile(condFile, null, new String[] {"CHR"}, null, null, null).keySet().toArray(new String[0]);
				if (chrs != null && chrs.length > 0) {
					tmp = new String[chrs.length];
					for (int k = 0; k < chrs.length; k++) {
						snpInfoFile = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chrs[k]);
						genoFile = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chrs[k]);
						genoFile = genoFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
	//					phenoFile = phenoDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, pheno);
						phenoFile = phenoDirAndNameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
						outputFile = resultsDir + "ARIC_" + pheno + "_" + ethnics[j] + "_cond" + currentCondition + "_chr" + chrs[k] + "_prepCondScores_seqMeta.RData";
						tmp[k] = outputFile;
	//					tmp[k] = rScriptDir + phenos[i] + "_" + ethnics[j] + "_cond" + startCondition + "_cohort_chr" + chrs[?] + ".R";
	//					Files.write(getRScriptForCohort(snpInfoFile, genoFile, phenoFile, condFile, outputFile), tmp[l]);
						rScript += (getRScriptForCohort(snpInfoFile, genoFile, phenoFile, condFile, outputFile) + "\n");
					}
					tmp2[j][0] = resultsDir + pheno + "_" + ethnics[j] + "_cond" + currentCondition + "_cohort_merged.RData";
	//				Files.write(getRScriptForMergingChromosomes(tmp, tmp2[k][0]), resultsDir + phenos[i] + "_" + ethnics[k] + "_cond" + startCondition + "_cohort_merge.R");
					rScript += (getRScriptForMergingChromosomes(tmp, tmp2[j][0]) + "\n");
				}
			}
	//		Files.write(getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFileDirAndNameTemplate, rScriptDir + phenos[i] + "_cond" + startCondition + "_", null), rScriptDir + phenos[i] + "_cond" + startCondition + "_meta.R");
			rScript += (getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFileDirAndNameTemplate, resultsDir + pheno + "_cond" + currentCondition, null));
			Files.write(rScript, rScriptDir + pheno + "_cond" + currentCondition + "_all.R");
			Files.qsub(rScriptDir + "[%0]_cond[%1]_all", rcommand + " " + rScriptDir + "[%0]_cond[%1]_all.R", new String[][] {{pheno, currentCondition + ""}}, -1, 6, -1);
			new File(rScriptDir + "master.qsub").renameTo(new File(rScriptDir + pheno + "_cond" + currentCondition + "_all.sh"));
	
			return rScript;
		}

	public static String getRScriptForConditionalAnalysis (String condFileFullPath, String snpInfoOfTheChromosome, String phenoFullPath, String genoFullPath, String[] ethnics, String outputDirAndRootTemplate, Logger log) {
			String phenoFile, genoFile, outputFile, rScript, tmp1, tmp2, tmp3;
	
			tmp2 = "";
			tmp3 = "";
			rScript = "library(seqMeta);\nlibrary(\"methods\");\n\ntmp <- load(\"" + snpInfoOfTheChromosome + "\");\nSNPInfo <- get(tmp);\nrm(list=tmp);\nrm(tmp)";
			for (int i = 0; i < ethnics.length; i++) {
				phenoFile = phenoFullPath.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]);
				genoFile = genoFullPath.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]);
				outputFile = outputDirAndRootTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]);
				tmp1 = ext.rootOf(outputFile).replaceAll("_" + FILENAME_ANALYSIS_SEGMENT, "");
				rScript += (genoFile.endsWith(".RData")?
							  "\n\ntmp <- load(\"" + genoFile + "\")\nZ <- get(tmp)\nrm(list=tmp)\nrm(tmp)\n"
							: "\n\nZ <- t(read.csv(\"" + genoFile + "\", header=T, as.is=T, row.names=1));\n")
						+ "names <- colnames(Z);\n"
						+ "for (i in 1:ncol(Z)) {\n"
						+ "	tmp <- names[i];\n"
						+ "	if (\"_\" == substr(tmp, start=nchar(tmp)-1, stop=nchar(tmp)-1)) {\n"
						+ "		names[i] = substr(tmp, start=1, stop=nchar(tmp)-2);\n"
						+ "	}\n"
						+ "}\n"
						+ "colnames(Z) <- names;\n"
						+ "\npheno <- na.omit(read.csv(\"" + phenoFile + "\", header=T, as.is=T, row.names=1));\n"
						+ "pheno_Z_merged <- merge(pheno, Z, by=\"row.names\")\n"
						+ "mPheno <- pheno_Z_merged[, (1 : ncol(pheno)) + 1]\n"
						+ "names <- colnames(pheno)\n"
						+ "if (length(names)>1) {\n"
						+ "	formu <- paste(names[1], \"~\", names[2])\n"
						+ "	for (i in 3:length(names)) {\n"
						+ "		formu <- paste(formu, \"+\", names[i])\n"
						+ "	}\n"
						+ "} else {\n"
						+ "	len <- length(mPheno)\n"
						+ "	mPheno <- c(mPheno, rep(1, len))\n"
						+ "	dim(mPheno) <- c(len, 2)\n"
						+ "	mPheno <- as.data.frame(mPheno)\n"
						+ "	colnames(mPheno) <- c(names[1], \"dummy\")\n"
						+ "	formu <- paste(names[1], \"~ 1\")\n"
						+ "}\n"
						+ "\nmGeno <- pheno_Z_merged[,(1 : ncol(Z)) + (1 + ncol(pheno))]"
						//TODO cond for a specific region or chr
						+ "\n\ncond <- read.table(\"" + condFileFullPath.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]) + "\", header=TRUE)"
						+ "\n\n" + tmp1 + " <- prepCondScores(Z=mGeno, formula(formu), SNPInfo=SNPInfo, adjustments=cond, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)"
						+ "\nsave(" + tmp1 + ", file=\"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "prepCondScores_seqMeta") + ".RData\", compress=\"bzip2\");\n"
						+ "\ngenes <- SNPInfo$SKATgene %in% cond$SKATgene"
						+ "\n\nresult <- singlesnpMeta(" + tmp1 + ", SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE);"
						+ "\nwrite.table(result, \"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP") + "\", sep=\",\", row.names = F);"
						+ "\nresult <- burdenMeta(" + tmp1 + ", SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1);"
						+ "\nwrite.table(result, \"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5Count") + "\", sep=\",\", row.names = F);"
						+ "\nresult <- skatMeta(" + tmp1 + ", SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)});"
						+ "\nwrite.table(result, \"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5") + "\", sep=\",\", row.names = F);"
						;
				tmp2 += (tmp1 + ", ");
				tmp3 += ethnics[i];
			}
			outputFile = outputDirAndRootTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, tmp3);
			if (! condFileFullPath.contains(FILENAME_ETHNIC_SEGMENT)) {
				rScript += ("\n\nresult <- singlesnpMeta(" + tmp2 + "SNPInfo=SNPInfo[genes,], snpNames = \"SNP\", aggregateBy=\"SKATgene\", studyBetas = TRUE);"
						+ "\nwrite.table(result, \"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP") + "\", sep=\",\", row.names = F);"
						+ "\nresult <- burdenMeta(" + tmp2 + "SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1);"
						+ "\nwrite.table(result, \"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5Count") + "\", sep=\",\", row.names = F);"
						+ "\nresult <- skatMeta(" + tmp2 + "SNPInfo=subset(SNPInfo[genes,], sc_nonsynSplice==TRUE), snpNames = \"SNP\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = function(maf) { dbeta(maf, 1, 25)*(maf <= 0.05)});"
						+ "\nwrite.table(result, \"" + outputFile.replaceAll(FILENAME_ANALYSIS_SEGMENT, "T5") + "\", sep=\",\", row.names = F);"
						);
			}
	
	//		Files.write(getRScriptForMetaAnalysis(snpInfoFile, ethnics, tmp2, condFileDirAndNameTemplate, rScriptDir + phenos[i] + "_cond" + startCondition + "_", null), rScriptDir + phenos[i] + "_cond" + startCondition + "_meta.R");
	//		Files.write(rScript, rScriptDir + pheno + "_cond" + currentCondition + "_all.R");
	//		Files.qsub(rScriptDir + "[%0]_cond[%1]_all", rcommand + " " + rScriptDir + "[%0]_cond[%1]_all.R", new String[][] {{pheno, currentCondition + ""}}, -1, 6, -1);
	
			return rScript;
		}

	public static String[] listFilesInDirAndAllSubDirs(String dir) {
		String[] result;
		String[] folders;

		folders = Files.listDirectories(dir, false);
		if (folders == null || folders.length < 1) {
				return Files.list(dir, null, ".RData", false, false, true);
		} else {
			result = new String[0];
			for (int i = 0; i < folders.length; i++) {
				result = combine(result, listFilesInDirAndAllSubDirs(dir + folders[i] + "/"));
			}
			return result;
		}
	}

	public static String[] combine(String[] a, String[] b) {
		String[] result;
		int j;

		result = new String[a.length + b.length];
		for (int i = 0; i < a.length; i++) {
			result[i] = a[i];
		}
		j = a.length;
		for (int i = 0; i < b.length; i++) {
			result[j] = b[i];
			j ++;
		}
		return result;
	}

	/* Working version of the method with the same name
	public static void summary(String resultsDir, int[] columnIndeciesOfPhenoConditionEthnicAnalysis, double pThreshold, String[] ethnicList, String fullpathToSnpInfo, String summaryDir, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvalSummary;
		String[] phenoList, analysesList, columnsAsTheKey, otherColumnsNeeded, otherColumnsNeededFromUnconditional;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll;
		String[] geneSnpList;
		Hashtable<String, String[]> snpInfo;
		Hashtable<String, Double> snpsWithSignificantPval;

		if (log == null) {
			log = new Logger();
		}

		snpInfo = loadSnpInfo(fullpathToSnpInfo);
		phenoGroups = groupFileNames(Files.list(resultsDir, null, ".csv", false, false), columnIndeciesOfPhenoConditionEthnicAnalysis, log);
		if (ethnicList == null || ethnicList.length < 1) {
			ethnicList = getEthnicList(phenoGroups, log);
		}
		if (phenoGroups != null && phenoGroups.size() > 0) {
//			analysesList = new String[] {"T5Count", "SKAT_T5"};
			analysesList = new String[] {"T5Count", "T5"};
			genePvalSummary = summarizeGenePvalues(phenoGroups, analysesList, resultsDir);
			phenoList = getPhenoList(genePvalSummary, log);
//			getListsOfPhenosConditionsEthnicsAnalysesGenes(genePvalSummary, phenoList, conditionListAllPhenos, ethnicList, analysesList, geneListAllPhenos, log);
			printSummarizedGenePvalues(genePvalSummary, phenoList, ethnicList, analysesList, summaryDir + "summary_genePvalues.xln", log);

			analysesList = new String[] {"SingleSNP"};
			columnsAsTheKey = new String[] {"Name", "gene"};	//TODO "Chr", "Position"
			otherColumnsNeeded = new String[] {"beta", "se", "p"};
			otherColumnsNeededFromUnconditional = new String[] {"ntotal", "nmiss", "maf"};
			for (String pheno : phenoList) {
				snpResultsAll = summarizeSnpPvalues(phenoGroups.get(pheno), snpInfo, resultsDir, analysesList, columnsAsTheKey, otherColumnsNeeded, otherColumnsNeededFromUnconditional, log);
				geneSnpList = getListOfGeneSnps(snpResultsAll, log);
				snpsWithSignificantPval = printSnpResults(snpResultsAll, getConditionList(genePvalSummary, pheno, log), ethnicList, analysesList, geneSnpList, otherColumnsNeeded, otherColumnsNeededFromUnconditional, new String[] {"p"}, pThreshold, summaryDir + "summary_" + pheno + "_snps.xln", log);
				printSummarizedSnpResults(snpResultsAll, snpsWithSignificantPval, getConditionList(genePvalSummary, pheno, log), ethnicList, analysesList, otherColumnsNeeded, summaryDir + "summary_" + pheno + "_snpPvalues.xln", log);
			}
		}
	}
	*/

	public static void summary(String resultsDir, int[] columnIndeciesOfPhenoConditionEthnicAnalysis, double pThreshold, String[] ethnics, String snpInfoDirNameTemplate, String condFileDirNameTemplate, String summaryDir, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvalSummary;
		String[] phenoList, analysesList, columnsAsTheKey, otherColumnsNeeded, otherColumnsNeededFromUnconditional;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll;
		String[] geneSnpList;
		Hashtable<String, Double> snpsWithSignificantPval;
		String tmp;

		if (log == null) {
			log = new Logger();
		}

		phenoGroups = groupFileNames(Files.list(resultsDir, null, ".csv", false, false), columnIndeciesOfPhenoConditionEthnicAnalysis, log);
		if (ethnics == null || ethnics.length < 1) {
			ethnics = getEthnicList(phenoGroups, log);
		}
		if (phenoGroups != null && phenoGroups.size() > 0) {
//			analysesList = new String[] {"T5Count", "SKAT_T5"};
			analysesList = new String[] {"T5Count", "T5"};
			genePvalSummary = summarizeGenePvalues(phenoGroups, analysesList, resultsDir, null, log);
			phenoList = getPhenoList(genePvalSummary, log);
//			getListsOfPhenosConditionsEthnicsAnalysesGenes(genePvalSummary, phenoList, conditionListAllPhenos, ethnicList, analysesList, geneListAllPhenos, log);
			printSummarizedGenePvalues(genePvalSummary, phenoList, ethnics, analysesList, summaryDir + "summary_" + Array.toStr(phenoList, "_") + "_genePvalues.xln", log);

			analysesList = new String[] {"SingleSNP"};
			columnsAsTheKey = new String[] {"Name", "gene"};	//TODO "Chr", "Position"
			otherColumnsNeeded = new String[] {"beta", "se", "p"};
			otherColumnsNeededFromUnconditional = new String[] {"ntotal", "nmiss", "maf"};
			for (String pheno : phenoList) {
				snpResultsAll = summarizeSnpPvalues(phenoGroups.get(pheno), null, resultsDir, analysesList, columnsAsTheKey, otherColumnsNeeded, otherColumnsNeededFromUnconditional, log);
				geneSnpList = getListOfGeneSnps(snpResultsAll, snpInfoDirNameTemplate, condFileDirNameTemplate, log);
				snpsWithSignificantPval = printSnpResults(snpResultsAll, getConditionList(genePvalSummary, pheno, log), ethnics, analysesList, geneSnpList, otherColumnsNeeded, otherColumnsNeededFromUnconditional, new String[] {"p"}, pThreshold, summaryDir + "summary_" + pheno + "_snps.xln", log);
				printSummarizedSnpResults(snpResultsAll, snpsWithSignificantPval, getConditionList(genePvalSummary, pheno, log), ethnics, analysesList, otherColumnsNeeded, summaryDir + "summary_" + pheno + "_snpPvalues.xln", log);
			}
		}
	}

	//TODO need to input condition list, but different pheno have different number of conditions and thus making the code complex.
	public static void summary(String resultsDirFilenameTemplate, double pThreshold, String[] phenos, String[] ethnics, String snpInfoDirNameTemplate, String conditionFileDirNameTemplate, String summaryDir, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvalSummary;
		String[] phenoList, analysesList, columnsAsTheKey, otherColumnsNeeded, otherColumnsNeededFromUnconditional;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll;
		String[] geneSnpList, ethnics1;
		Hashtable<String, Double> snpsWithSignificantPval;
		String resultDir, allEthnics;
//		int[] columnIndeciesOfPhenoConditionEthnicAnalysis;
		Hashtable <String, Hashtable <Integer, Vector<String>>> phenoToRegionToGenes;

		if (log == null) {
			log = new Logger();
		}

		allEthnics = "";
		for (int i = 0; i < ethnics.length; i++) {
			allEthnics += ethnics[i];
		}
		ethnics1 = Array.addStrToArray(allEthnics, ethnics);

		resultsDirFilenameTemplate = resultsDirFilenameTemplate.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, "");
//		columnIndeciesOfPhenoConditionEthnicAnalysis = getIndicesOfPhenoConditionEthnicAnalysis(resultsDirFilenameTemplate, "_");
		resultDir = ext.parseDirectoryOfFile(resultsDirFilenameTemplate);
//		phenoGroups = groupFileNames(Files.list(resultDir, null, ".csv", false, false), columnIndeciesOfPhenoConditionEthnicAnalysis, phenos, log);

		phenoGroups = groupFileNames(resultsDirFilenameTemplate, phenos, ethnics1, new String[] {"T5Count", "T5", "SingleSNP"}, log);

//		if (ethnics == null || ethnics.length < 1) {
//			ethnics = getEthnicList(phenoGroups, log);
//		}

		if (phenoGroups != null && phenoGroups.size() > 0) {
//			analysesList = new String[] {"T5Count", "SKAT_T5"};
			analysesList = new String[] {"T5Count", "T5"};
			phenoList = phenoGroups.keySet().toArray(new String[0]);
			Arrays.sort(phenoList);
			phenoToRegionToGenes = getConditionsGeneRegions(conditionFileDirNameTemplate, phenoList, log);
			genePvalSummary = summarizeGenePvalues(phenoGroups, analysesList, resultDir, phenoToRegionToGenes, log);
//			getListsOfPhenosConditionsEthnicsAnalysesGenes(genePvalSummary, phenoList, conditionListAllPhenos, ethnicList, analysesList, geneListAllPhenos, log);
			printSummarizedGenePvalues(genePvalSummary, phenoList, ethnics1, analysesList, summaryDir + "summary_" + Array.toStr(phenoList, "_") + "_genePvalues.xln", log);

			analysesList = new String[] {"SingleSNP"};
			columnsAsTheKey = new String[] {"Name", "gene"};	//TODO "Chr", "Position"
			otherColumnsNeeded = new String[] {"beta", "se", "p"};
			otherColumnsNeededFromUnconditional = new String[] {"ntotal", "nmiss", "maf"};
			for (String pheno : phenoList) {
				snpResultsAll = summarizeSnpPvalues(phenoGroups.get(pheno), null, resultDir, analysesList, columnsAsTheKey, otherColumnsNeeded, otherColumnsNeededFromUnconditional, log);
//				geneSnpList = getListOfGeneSnps(snpResultsAll, snpInfoDirNameTemplate, condFileDirNameTemplate, log);
				geneSnpList = getListOfGeneSnps(snpResultsAll, snpInfoDirNameTemplate, null, log);
				snpsWithSignificantPval = printSnpResults(snpResultsAll, getConditionList(genePvalSummary, pheno, log), ethnics1, analysesList, geneSnpList, otherColumnsNeeded, otherColumnsNeededFromUnconditional, new String[] {"p"}, pThreshold, summaryDir + "summary_" + pheno + "_snps.xln", log);
				printSummarizedSnpResults(snpResultsAll, snpsWithSignificantPval, getConditionList(genePvalSummary, pheno, log), ethnics1, analysesList, otherColumnsNeeded, summaryDir + "summary_" + pheno + "_snpPvalues.xln", log);
			}
		}
	}

	public static Hashtable <String, Hashtable <Integer, Vector<String>>> getConditionsGeneRegions (String conditionFileDirNameTemplate, String[] phenos, Logger log) {
		Hashtable <String, Hashtable <Integer, Vector<String>>> result;
		Object[] test;
		String condFile;

		result = new Hashtable <String, Hashtable <Integer, Vector<String>>>();
		for (int i = 0; i < phenos.length; i++) {
			condFile = conditionFileDirNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]).replaceAll(FILENAME_CONDITION_SEGMENT, "cond1");
			result.put(phenos[i], (Hashtable<Integer, Vector<String>>) organizeConditionsIntoRegions(condFile, log)[1]);
		}

		return result;
	}

	public static int[] getIndicesOfPhenoConditionEthnicAnalysis(String fileDirAndNameTemplate, String delimiter) {
		int[] indices;
		String[] line, fileNameSegments;

		fileNameSegments = new String[] {FILENAME_PHENO_SEGMENT, FILENAME_CONDITION_SEGMENT, FILENAME_ETHNIC_SEGMENT, FILENAME_ANALYSIS_SEGMENT};
		line = ext.rootOf(fileDirAndNameTemplate).replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, "").split(delimiter);
		indices = ext.indexFactors(fileNameSegments, line, true, false);

		return indices;
	}

	public static Hashtable<String, String> lookupChrPosInHashFormatFromSnpInfoFiles (String snpInfoDirFilenameTemplate, String[] snps, String[] chromosomes, Logger log) {
		String[] chrPos;
		Hashtable<String, String> result;

		chrPos = lookupChrPosFromSnpInfoFiles (snpInfoDirFilenameTemplate, snps, chromosomes, log);
		result = new Hashtable<String, String>();
		if (chrPos != null) {
			for (int i = 0; i < chrPos.length; i++) {
				result.put(snps[i], chrPos[i]);
			}
		}

		return result;
	}

	public static String[] lookupChrPosFromSnpInfoFiles (String snpInfoDirFilenameTemplate, String[] snps, String[] chromosomes, Logger log) {
		String[] result;
		BufferedReader reader;
		String delimiter, snpInfoFile, tmp;
		int[] indices;
		String[] line;
		Vector<Integer> notYetFound;
		int loop;

		if (log == null) {
			log = new Logger();
		}

		result = new String[snps.length];
		notYetFound = new Vector<Integer>(result.length);
		for (int i = 0; i < result.length; i++) {
			snps[i] = snps[i].trim();
			line = snps[i].split(":");
			if (line.length <= 1) {
				notYetFound.add(i);
			} else {
				try {
					Integer.parseInt(line[1]);
//					result[i] = snps[i].replaceFirst(":", "\t");
					result[i] = line[0] + "\t" + line[1];
				} catch (NumberFormatException e) {
					notYetFound.add(i);
				}
			}
		}

		if (snpInfoDirFilenameTemplate != null && (snpInfoDirFilenameTemplate.endsWith(".csv") || snpInfoDirFilenameTemplate.endsWith(".xln") || snpInfoDirFilenameTemplate.endsWith(".csv.gz") || snpInfoDirFilenameTemplate.endsWith(".xln.gz"))) {
			if (chromosomes == null) {
				if (snpInfoDirFilenameTemplate != null && snpInfoDirFilenameTemplate.contains(FILENAME_CHROMOSOME_SEGMENT)) {
					chromosomes = new String[] {"1", "2", "3", "4", "5", "6", "7", 	"8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"};
				} else {
					chromosomes = new String[] {"Null"};
				}
			}
	
			for (int i = 0; notYetFound.size() > 0 && i < chromosomes.length; i++) {
				snpInfoFile = snpInfoDirFilenameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, chromosomes[i]);
				try {
					delimiter = Files.determineDelimiter(snpInfoFile, null);
					reader = Files.getAppropriateReader(snpInfoFile);
					line = reader.readLine().split(delimiter);
					indices = ext.indexFactors(SNPINFO_COLUMNS, line, false, true);
					while (notYetFound.size() > 0 && reader.ready()) {
						line = reader.readLine().split(delimiter);
						loop = notYetFound.size();
						for (int j = 0; j < loop; j++) {
							if (line[indices[0]].equalsIgnoreCase(snps[notYetFound.elementAt(j)])) {
								result[notYetFound.elementAt(j)] = line[indices[1]] + "\t" + line[indices[2]];
								notYetFound.remove(j);
								break;
							}
						}
					}
					reader.close();
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		} else {
			log.report("Warning - Snp Info File is not available for the summarize result process: " + snpInfoDirFilenameTemplate);
		}

		if (notYetFound.size() > 0) {
			loop = Math.min(notYetFound.size(), 5);
			tmp = "";
			for (int i = 0; i < loop; i++) {
				tmp += ("\n" + snps[notYetFound.elementAt(i)]);
			}
			if (notYetFound.size() > loop) {
				tmp += "\n...";
			}
			log.report("Warning - chromosomes and locations of " + notYetFound.size() + " SNPs are not found in the SNP Info files and chromosome range sepecified: " + tmp);
		}

		if (notYetFound.size() == snps.length) {
			return null;
		} else {
			return result;
		}
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> summarizeGenePvalues_20150109 (Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups, String[] analysesTypesToSelect, String resultsDir, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> summary = null;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> summaryConditionGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> summaryEthnicGroup;
		Hashtable<String, Hashtable<String, String>> summaryAnalysesGroup;
		Hashtable<String, String> summaryGeneGroup;
		BufferedReader reader;
		String[] line;
		int[] indices;
		int tmp = 0;
		boolean isError;

		if (log == null) {
			log = new Logger();
		}

		if (phenoGroups != null && phenoGroups.size() > 0) {
			summary = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>>();
			for (String pheno : phenoGroups.keySet()) {
				conditionGroup = phenoGroups.get(pheno);
				if (! summary.containsKey(pheno)) {
					summaryConditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>();
					summary.put(pheno, summaryConditionGroup);
				} else {
					summaryConditionGroup = summary.get(pheno);
				}

				for (String condition : conditionGroup.keySet()) {
					ethnicGroup = conditionGroup.get(condition);
					if (! summaryConditionGroup.containsKey(condition)) {
						summaryEthnicGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String>>>();
						summaryConditionGroup.put(condition, summaryEthnicGroup);
					} else {
						summaryEthnicGroup = summaryConditionGroup.get(condition);
					}

					for (String ethnic : ethnicGroup.keySet()) {
						analysesGroup = ethnicGroup.get(ethnic);
						if (! summaryEthnicGroup.containsKey(ethnic)) {
							summaryAnalysesGroup = new Hashtable<String, Hashtable<String, String>>();
							summaryEthnicGroup.put(ethnic, summaryAnalysesGroup);
						} else {
							summaryAnalysesGroup = summaryEthnicGroup.get(ethnic);
						}

						for (String analysis : analysesTypesToSelect) {
							if (analysesGroup.containsKey(analysis)) {
								if (! summaryAnalysesGroup.containsKey(analysis)) {
									summaryGeneGroup = new Hashtable<String, String>();
									summaryAnalysesGroup.put(analysis, summaryGeneGroup);
								} else {
									summaryGeneGroup = summaryAnalysesGroup.get(analysis);
								}
								try {
//									reader = new BufferedReader(new FileReader(resultsDir + analysesGroup.get(analysis)));
									reader = Files.getAppropriateReader(resultsDir + analysesGroup.get(analysis));
									indices = ext.indexFactors(GENE_RESULT_COLUMNS, reader.readLine().replaceAll("\"", "").split(","), false, false);
									isError = false;
									for (int i = 0; i < indices.length; i++) {
										if (indices[i] < 0) {
											tmp = i;
											isError = true;
											break;
										}
									}
									if (isError) {
										log.reportError("Error - Cannot find the column name '" + GENE_RESULT_COLUMNS[tmp] + "' in the following file:\n" + resultsDir + analysesGroup.get(analysis));
									} else {
										while (reader.ready()) {
											line = reader.readLine().replaceAll("\"", "").split(",");
											summaryGeneGroup.put(line[indices[0]], line[indices[1]]);
										}
									}
									reader.close();
								} catch (FileNotFoundException e) {
									e.printStackTrace();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						}
					}
				}
			}
		}
		
		return summary;
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> summarizeGenePvalues (Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups, String[] analysesTypesToSelect, String resultsDir, Hashtable<String, Hashtable<Integer, Vector<String>>> phenoToRegionToGenes, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> summary = null;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> summaryConditionGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> summaryEthnicGroup;
		Hashtable<String, Hashtable<String, String>> summaryAnalysesGroup;
		Hashtable<String, String> summaryGeneGroup;
		Hashtable<Integer, Vector<String>> regionToGenes;
		BufferedReader reader;
		String[] line;
		int[] indices;
		int tmp = 0;
		boolean isError;

		if (log == null) {
			log = new Logger();
		}

		if (phenoGroups != null && phenoGroups.size() > 0) {
			summary = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>>();
			for (String pheno : phenoGroups.keySet()) {
				conditionGroup = phenoGroups.get(pheno);
				if (! summary.containsKey(pheno)) {
					summaryConditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>();
					summary.put(pheno, summaryConditionGroup);
				} else {
					summaryConditionGroup = summary.get(pheno);
				}
				regionToGenes = phenoToRegionToGenes.get(pheno);

				for (String condition : conditionGroup.keySet()) {
					ethnicGroup = conditionGroup.get(condition);
					if (! summaryConditionGroup.containsKey(condition)) {
						summaryEthnicGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String>>>();
						summaryConditionGroup.put(condition, summaryEthnicGroup);
					} else {
						summaryEthnicGroup = summaryConditionGroup.get(condition);
					}

					for (String ethnic : ethnicGroup.keySet()) {
						analysesGroup = ethnicGroup.get(ethnic);
						if (! summaryEthnicGroup.containsKey(ethnic)) {
							summaryAnalysesGroup = new Hashtable<String, Hashtable<String, String>>();
							summaryEthnicGroup.put(ethnic, summaryAnalysesGroup);
						} else {
							summaryAnalysesGroup = summaryEthnicGroup.get(ethnic);
						}

						for (String analysis : analysesTypesToSelect) {
							if (analysesGroup.containsKey(analysis)) {
								if (! summaryAnalysesGroup.containsKey(analysis)) {
									summaryGeneGroup = new Hashtable<String, String>();
									summaryAnalysesGroup.put(analysis, summaryGeneGroup);
								} else {
									summaryGeneGroup = summaryAnalysesGroup.get(analysis);
								}
								try {
//									reader = new BufferedReader(new FileReader(resultsDir + analysesGroup.get(analysis)));
									reader = Files.getAppropriateReader(resultsDir + analysesGroup.get(analysis));
									indices = ext.indexFactors(GENE_RESULT_COLUMNS, reader.readLine().replaceAll("\"", "").split(","), false, false);
									isError = false;
									for (int i = 0; i < indices.length; i++) {
										if (indices[i] < 0) {
											tmp = i;
											isError = true;
											break;
										}
									}
									if (isError) {
										log.reportError("Error - Cannot find the column name '" + GENE_RESULT_COLUMNS[tmp] + "' in the following file:\n" + resultsDir + analysesGroup.get(analysis));
									} else {
										while (reader.ready()) {
											line = reader.readLine().replaceAll("\"", "").split(",");
											summaryGeneGroup.put((regionToGenes == null? "" : (getKeyOfValue(regionToGenes, line[indices[0]], log) + "\t")) + line[indices[0]], line[indices[1]]);
										}
									}
									reader.close();
								} catch (FileNotFoundException e) {
									e.printStackTrace();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						}
					}
				}
			}
		}
		
		return summary;
	}

	public static String[] getGeneList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvaluesSummary, String pheno, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> conditionGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> ethnicGroup;
		Hashtable<String, Hashtable<String, String>> analysesGroup;
		Hashtable<String, String> geneGroup;
		HashSet<String> geneList;
		String[] geneListSorted = null;

		if (genePvaluesSummary != null && genePvaluesSummary.size() > 0) {
			geneList = new HashSet<String> ();
			conditionGroup = genePvaluesSummary.get(pheno);

			for (String condition : conditionGroup.keySet()) {
				ethnicGroup = conditionGroup.get(condition);

				for (String ethnic : ethnicGroup.keySet()) {
					analysesGroup = ethnicGroup.get(ethnic);

					for (String analysis : analysesGroup.keySet()) {
						geneGroup = analysesGroup.get(analysis);

						for (String gene : geneGroup.keySet()) {
							if (! geneList.contains(gene)) {
								geneList.add(gene);
							}
						}
					}
				}
			}

			geneListSorted = geneList.toArray(new String[0]);
			Arrays.sort(geneListSorted);
		}
		return geneListSorted;
	}

	public static String[] getPhenoList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvaluesSummary, Logger log) {
		String[] result;

		result = genePvaluesSummary.keySet().toArray(new String[0]);
		Arrays.sort(result);

		return result;
	}

	public static String[] getConditionList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvaluesSummary, String pheno, Logger log) {
		String[] conditionList = null;

		if (genePvaluesSummary != null && genePvaluesSummary.size() > 0) {
			conditionList = genePvaluesSummary.get(pheno).keySet().toArray(new String[0]);
			Arrays.sort(conditionList);
		}

		return conditionList;
	}

	public static String[] getEthnicsList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvalueSummary, String pheno, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> conditionGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> ethnicGroup;
		HashSet<String> ethnicList;
		String[] a = null;

		if (genePvalueSummary != null && genePvalueSummary.size() > 0) {
			ethnicList = new HashSet<String> ();
			conditionGroup = genePvalueSummary.get(pheno);
			for (String condition : conditionGroup.keySet()) {
				ethnicGroup = conditionGroup.get(condition);
				for (String ethnic : ethnicGroup.keySet()) {
					ethnicList.add(ethnic);
				}
			}

			a = ethnicList.toArray(new String[0]);
			Arrays.sort(a);
		}
		
		return a;
	}

	public static String[] getEthnicList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups, Logger log) {
		Vector<String> result;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
	
		if (log == null) {
			log = new Logger();
		}
	
		result = new Vector<String> ();
		if (phenoGroups != null && phenoGroups.size() > 0) {
			for (String pheno : phenoGroups.keySet()) {
				conditionGroup = phenoGroups.get(pheno);
				for (String condition : conditionGroup.keySet()) {
					ethnicGroup = conditionGroup.get(condition);
					for (String ethnic : ethnicGroup.keySet()) {
						if (! result.contains(ethnic)) {
							result.add(ethnic);
						}
					}
				}
			}
		}
	
		if (! result.isEmpty()) {
			return result.toArray(new String[0]);
		} else {
			return null;
		}
	}

	public static String[] getEthnicsList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> ethnicGroup;
		HashSet<String> ethnicList;
		String[] ethnicListSorted = null;

		if (snpResultsAll != null && snpResultsAll.size() > 0) {
			ethnicList = new HashSet<String> ();
			for (String condition : snpResultsAll.keySet()) {
				ethnicGroup = snpResultsAll.get(condition);
				for (String ethnic : ethnicGroup.keySet()) {
					ethnicList.add(ethnic);
				}
			}

			ethnicListSorted = ethnicList.toArray(new String[0]);
			Arrays.sort(ethnicListSorted);
		}
		
		return ethnicListSorted;
	}

	public static String[] getAnalysesList(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvaluesSummary, String pheno, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> conditionGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> ethnicGroup;
		Hashtable<String, Hashtable<String, String>> analysesGroup;
		HashSet<String> analysesList;
		String[] analysesListSorted = null;

		if (genePvaluesSummary != null && genePvaluesSummary.size() > 0) {
			analysesList = new HashSet<String> ();
			conditionGroup = genePvaluesSummary.get(pheno);
			for (String condition : conditionGroup.keySet()) {
				ethnicGroup = conditionGroup.get(condition);
				for (String ethnic : ethnicGroup.keySet()) {
					analysesGroup = ethnicGroup.get(ethnic);
					for (String analysis : analysesGroup.keySet()) {
						analysesList.add(analysis);
					}
				}
			}

			analysesListSorted = analysesList.toArray(new String[0]);
			Arrays.sort(analysesListSorted);
		}

		return analysesListSorted;
	}

	public static void printSummarizedGenePvalues(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvaluesSummary, String[] phenos, String[] ethnics, String[] analyses, String fullPathOutFilename, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> conditionGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> ethnicGroup;
		Hashtable<String, Hashtable<String, String>> analysesGroup;
		Hashtable<String, String> geneGroup;
		String[] conditionList, geneList;
		PrintWriter writer;
		String line;

		if (log == null) {
			log = new Logger();
		}

		if (genePvaluesSummary != null && genePvaluesSummary.size() > 0) {
			try {
				writer = new PrintWriter(new FileOutputStream(fullPathOutFilename));
//				line = "pheno\tethnic\tregionId\tgene";
				for (String pheno : phenos) {
					if (ethnics == null) {
						ethnics = getEthnicsList(genePvaluesSummary, pheno, log);
					}
					if (analyses == null) {
						analyses = getAnalysesList(genePvaluesSummary, pheno, log);
					}
					writer.println(pheno);
					line = "ethnic\tregionId\tgene";
					conditionList = getConditionList(genePvaluesSummary, pheno, log);
					for (String condition : conditionList) {
						for (String analysis : analyses) {
							line += ("\t" + condition + "_" + analysis);
						}
					}
					writer.println(line);

					geneList = getGeneList(genePvaluesSummary, pheno, log);
					conditionGroup = genePvaluesSummary.get(pheno);
					for (String ethnic : ethnics) {
						for (String gene : geneList) {
							line = (ethnic + "\t" + gene);
							for (String condition : conditionList) {
								if (conditionGroup.containsKey(condition)) {
									ethnicGroup = conditionGroup.get(condition);
									if (ethnicGroup.containsKey(ethnic)) {
										analysesGroup = ethnicGroup.get(ethnic);
										for (String analysis : analyses) {
											geneGroup = analysesGroup.get(analysis);
											line += ("\t" + (geneGroup.containsKey(gene)? geneGroup.get(gene) : ""));
										}
									} else {
										for (int i = 0; i < analyses.length; i++) {
											line += "\t";
										}
									}
								}
							}
							writer.println(line);
						}
					}
					writer.println();
				}

				writer.close();
				log.report("Summary for gene p-values is ready at " + fullPathOutFilename);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> summarizeSnpPvalues(Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroups, Hashtable<String, String[]> snpInfo, String dir, String[] analysesNeeded, String[] columnsAsTheKey, String[] otherColumnsNeeded, String[] otherColumnsNeededForUnconditional, Logger log) {
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll = null;
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> snpResultsConditionGroup;
		Hashtable<String, Hashtable<String, String[]>> snpResultsEthnicGroup;
		Hashtable<String, String[]> snpResultsEachFile;
//		HashSet<String> geneSnpList;

		if (conditionGroups != null && conditionGroups.size() > 0) {
			snpResultsAll = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>>();
//			geneSnpList = new HashSet<String>();
			for (String condition : conditionGroups.keySet()) {
				ethnicGroup = conditionGroups.get(condition);
				if (! snpResultsAll.containsKey(condition)) {
					snpResultsConditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>();
					snpResultsAll.put(condition, snpResultsConditionGroup);
				} else {
					snpResultsConditionGroup = snpResultsAll.get(condition);
				}
				for (String ethnic : ethnicGroup.keySet()) {
					analysesGroup = ethnicGroup.get(ethnic);
					if (! snpResultsConditionGroup.containsKey(ethnic)) {
						snpResultsEthnicGroup = new Hashtable<String, Hashtable<String, String[]>>();
						snpResultsConditionGroup.put(ethnic, snpResultsEthnicGroup);
					} else {
						snpResultsEthnicGroup = snpResultsConditionGroup.get(ethnic);
					}
					for (String analysis : analysesNeeded) {	//analysesGroup.keySet()
						if (condition.equals("cond1")) {
							snpResultsEachFile = loadFile(dir + analysesGroup.get(analysis), snpInfo, columnsAsTheKey, otherColumnsNeededForUnconditional, null, null);
							snpResultsEthnicGroup.put(analysis + "_totals", snpResultsEachFile);
						}
						snpResultsEachFile = loadFile(dir + analysesGroup.get(analysis), snpInfo, columnsAsTheKey, otherColumnsNeeded, null, null);
						snpResultsEthnicGroup.put(analysis, snpResultsEachFile);
//						for (String geneSnp : snpResultsEachFile.keySet()) {
//							if(! geneSnpList.contains(geneSnp)) {
//								geneSnpList.add(geneSnp);
//							}
//						}
					}
				}
			}
		}

		return snpResultsAll;
	}

	public static Hashtable<String, String[]> loadFile1(String filefullpath, Hashtable<String, String[]> snpList, String[] columnsAsTheKey, String[] otherColumnsNeeded, String[] criteriaColumns, String[] criteria) {
		return null;
	}

	/**
	 * 
	 * @param filefullpath
	 * @param snpList	String[] to be included in the output of this method
	 * @param namesOfColumnsToBeIncludedInTheKeyOfTheOutputHash
	 * @param namesOfColumnsToBeIncludedInTheValueOfTheOutputHash
	 * @param criteria	new String[] {"(ARIC_Whites_p_SingleSNP || ARIC_Blacks_p_SingleSNP) <= 0.00001", "(ARIC_Whites_maf_SingleSNP || ARIC_Blacks_maf_SingleSNP) > 0"}
	 * 					Note: "AND" relationship among rows, and "OR" relationship within a row.
	 * @param log
	 * @return
	 */
	public static Hashtable<String, String[]> loadFile(String filefullpath, Hashtable<String, String[]> snpList, String[] namesOfColumnsToBeIncludedInTheKeyOfTheOutputHash, String[] namesOfColumnsToBeIncludedInTheValueOfTheOutputHash, String[] criteria, Logger log) {
		boolean isToInclude;
		int[] indicesKey, indicesOther = null;
		String key, delimiter, header;
		String[] line, tmp, criteriaOperators = null, criteriaValues = null;
		Hashtable<String, String[]> result;
		BufferedReader reader;
		Object[] tmp2;
		int[][] criteriaColumnIndices = null;

		if (log == null) {
			log = new Logger();
		}

		result = new Hashtable<String, String[]>();
		try {
			reader = new BufferedReader(new FileReader(filefullpath));
			header = reader.readLine();
			delimiter = ext.determineDelimiter(header);
			line = header.replaceAll("\"", "").split(delimiter);
			indicesKey = ext.indexFactors(namesOfColumnsToBeIncludedInTheKeyOfTheOutputHash, line, false, true);
			if (namesOfColumnsToBeIncludedInTheValueOfTheOutputHash != null) {
				indicesOther = ext.indexFactors(namesOfColumnsToBeIncludedInTheValueOfTheOutputHash, line, false, false);
			}
			if (criteria != null) {
				tmp2 = parseCriteria(line, criteria, log);
				criteriaColumnIndices = (int[][]) tmp2[0];
				criteriaOperators = (String[]) tmp2[1];
				criteriaValues = (String[]) tmp2[2];
			}

			while (reader.ready()) {
				line = reader.readLine().replaceAll("\"", "").split(delimiter);
				
				if (criteria != null) {
					isToInclude = screenWithCriteria(line, criteriaColumnIndices, criteriaOperators, criteriaValues, log);
				} else {
					isToInclude = true;
				}

				if (isToInclude) {
					key = line[indicesKey[0]];
					if (snpList != null) {
						tmp = snpList.get(key);
					} else {
						tmp = null;
					}
					for (int i = 1; i < indicesKey.length; i++) {
						key += ("\t" + line[indicesKey[i]]);
					}
					if (snpList != null) {
						for (int i = 0; i < tmp.length; i++) {
							key += ("\t" + tmp[i]);
						}
					}
					if (indicesOther != null) {
						tmp = new String[indicesOther.length];
					} else {
						tmp = new String[0];
					}
					for (int i = 0; i < tmp.length; i++) {
						if (indicesOther[i] >= 0) {
							tmp[i] = line[indicesOther[i]];
						} else {
							tmp[i] = "";
						}
					}
					result.put(key, tmp);
				}
			}
			reader.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return result;
	}

	
	public static Vector<String> loadFile(String filefullpath, int numHeaderLinesToSkip) {
		BufferedReader reader;
		Vector<String> result;

		result = new Vector<String>();
		try {
			reader = new BufferedReader(new FileReader(filefullpath));
			for (int i = 0; i < numHeaderLinesToSkip; i++) {
				reader.readLine();
			}
			while (reader.ready()) {
				result.add(reader.readLine());
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return result;
	}

	private static Hashtable<String, String[]> loadSnpInfo(String fullpathToSnpInfo) {
		BufferedReader reader;
		Hashtable<String, String[]> result;
		String delimiter;
		int[] indices;
		String[] line;

		result = new Hashtable<String, String[]>();
		try {
			delimiter = Files.determineDelimiter(fullpathToSnpInfo, null);
			reader = Files.getAppropriateReader(fullpathToSnpInfo);
			line = reader.readLine().split(delimiter);
			indices = ext.indexFactors(SNPINFO_COLUMNS, line, false, true);
			while (reader.ready()) {
				line = reader.readLine().split(delimiter);
				result.put(line[indices[0]], new String[] {line[indices[1]], line[indices[2]]});
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return result;
	}

	/**
	 * Parse the criteria input by the users.
	 * @param headerFromTheFile: new String[] {"p_min", "p_EA", "maf_EA", "p_AA", "p_AA"}
	 * @param criteria: new String[] {"maf==.5", "Ethnic==AA", "p_min<=.00005", "(p_min||p_EA||p_AA)<=.00005", "( p_min || p_EA || p_AA )<=.00005"}
	 * @return
	 */
	public static Object[] parseCriteria (String[] headerFromTheFile, String[] criteria, Logger log) {
		int end, a;
		int[][] criteriaColumnIndices;
		String[] criteriaOperators, criteriaValues, tmp;

		if (log == null) {
			log = new Logger();
		}

		criteriaColumnIndices = new int[criteria.length][];
		criteriaOperators = new String[criteria.length];
		criteriaValues = new String[criteria.length];

		for (int i = 0; i < criteria.length; i++) {
			for (int j = 0; j < OPERATORS.length; j++) {
				tmp = criteria[i].split(OPERATORS[j]);
				if (tmp.length > 1) {
					criteriaOperators[i] = OPERATORS[j];
					criteriaValues[i] = tmp[1].trim();

					tmp = tmp[0].split("\\|\\|");
					tmp[0] = tmp[0].substring(tmp[0].indexOf("(") + 1);
					end = tmp.length - 1;
					a = tmp[end].indexOf(")");
					if (a >= 0) {
						tmp[end] = tmp[end].substring(0, a);
					}
					for (int k = 0; k < tmp.length; k++) {
						tmp[k] = tmp[k].trim();
					}
					criteriaColumnIndices[i] = ext.indexFactors(tmp, headerFromTheFile, false, false);
					break;
				}
			}

			if (criteriaOperators[i] == null || criteriaOperators[i].equals("")) {
				log.reportError("Error - operator specified in the criteria (" + criteria[i] + ") is not recognized by the current version of software, and is ignored.");
			}
		}

		return new Object[] {criteriaColumnIndices, criteriaOperators, criteriaValues};
	}

	/**
	 * To test whether a row of data fit the criteria or not. 
	 * @param strToBeScreened: one row from the data as the example below
	 * 				(col#1)	(col#2)	(col#3)	(col#4)
	 * 		(row#1)	.013	11.4589	12		A
	 * 		(row#2)	.035	4.93	13		A
	 * 		(row#3)	.913	7.395	12		B
	 * 				...
	 * 
	 * To apply the following criteria to the above data:
	 * 				(col#1)	(col#2)	(col#3)	(col#4)
	 * 				<=.05 Or <=5.0
	 * 			And					=13
	 * 			And							=A
	 *
	 * @param criteriaColumnIndices:	new int[][] {{1,2}, {3}, {4}}	Note: "OR" relationships are put in same row, while "AND" relationships are put in different rows.
	 * @param criteriaOperator:			new String[] {"<=", "==", "=="}
	 * @param criteriaValues:			new String[] {".05", "13", "A"};
	 * @return
	 */
	public static boolean screenWithCriteria (String[] strToBeScreened, int[][] criteriaColumnIndices, String[] criteriaOperators, String[] criteriaValues, Logger log) {
		boolean isMatchCriteria, isMatchCurrentCriterion;
		double value, criteriaValue;

		if (log == null) {
			log = new Logger();
		}

		isMatchCriteria = true;
		for (int i = 0; i < criteriaValues.length; i++) {
			isMatchCurrentCriterion = false;
			for (int j = 0; j < criteriaColumnIndices[i].length; j++) {
				if (criteriaColumnIndices[i][j] >= 0) {
					try {
						criteriaValue = Double.parseDouble(criteriaValues[i]);
						value = Double.parseDouble(strToBeScreened[criteriaColumnIndices[i][j]]);
						if ((criteriaOperators[i].equalsIgnoreCase("<=") && value <= criteriaValue)
								|| (criteriaOperators[i].equalsIgnoreCase("==") && value == criteriaValue)
								|| (criteriaOperators[i].equalsIgnoreCase(">=") && value >= criteriaValue)
								|| (criteriaOperators[i].equalsIgnoreCase("!=") && value != criteriaValue)
								|| (criteriaOperators[i].equalsIgnoreCase("<") && value < criteriaValue)
								|| (criteriaOperators[i].equalsIgnoreCase(">") && value > criteriaValue)) {
							isMatchCurrentCriterion = true;
							break;
						}
					} catch (NumberFormatException error) {
						if (criteriaOperators[i].equalsIgnoreCase("==") && strToBeScreened[criteriaColumnIndices[i][j]].equalsIgnoreCase(criteriaValues[i])) {
							isMatchCurrentCriterion = true;
							break;
						} else if (criteriaOperators[i].equalsIgnoreCase("!=") && ! strToBeScreened[criteriaColumnIndices[i][j]].equalsIgnoreCase(criteriaValues[i])) {
							isMatchCurrentCriterion = true;
							break;
						}
					}
	
				}
			}
			if(! isMatchCurrentCriterion) {
				isMatchCriteria = false;
				break;
			}
		}

		return isMatchCriteria;
	}

	public static String[] getListOfGeneSnps(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, String snpInfoDirNameTemplate, String condFileDirNameTemplate, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> snpResultsCurrentCondition;
		Hashtable<String, Hashtable<String, String[]>> snpResultsCurrentEthnic;
		Hashtable<String, String[]> snpResultsCurrentAnalysis;
		HashSet<String> geneSnpList = null;
		String[] snpGeneList, snpList, chrList;
		Hashtable<String, String> chrPositions;
		String snpGeneChrPos;

		snpGeneList = getListOfGeneSnps(snpResultsAll, log);

		snpList = new String[snpGeneList.length];
		for (int i = 0; i < snpList.length; i++) {
			snpList[i] = snpGeneList[i].split("\t")[0];
		}
		chrPositions = lookupChrPosInHashFormatFromSnpInfoFiles(snpInfoDirNameTemplate, snpList, null, log);
		if (chrPositions.size() > 0) {
			if (snpResultsAll != null && snpResultsAll.size() > 0) {
				geneSnpList = new HashSet<String>();
				for (String condition : snpResultsAll.keySet()) {
					snpResultsCurrentCondition = snpResultsAll.get(condition);
					for (String ethnic : snpResultsCurrentCondition.keySet()) {
						snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
						for (String analysis : snpResultsCurrentEthnic.keySet()) {	//analysesGroup.keySet()
							snpResultsCurrentAnalysis = snpResultsCurrentEthnic.get(analysis);
							snpList = snpResultsCurrentAnalysis.keySet().toArray(new String[0]);
							for (int i = 0; i < snpList.length; i++) {
								snpGeneChrPos = snpList[i] + "\t" + chrPositions.get(snpList[i].split("\t")[0]);
								snpResultsCurrentAnalysis.put(snpGeneChrPos, snpResultsCurrentAnalysis.get(snpList[i]));
								snpResultsCurrentAnalysis.remove(snpList[i]);
								if(! geneSnpList.contains(snpGeneChrPos)) {
									geneSnpList.add(snpGeneChrPos);
								}
							}
						}
					}
				}
				snpGeneList = geneSnpList.toArray(new String[0]);
				Arrays.sort(snpGeneList);
			}
		}

		return snpGeneList;
	}

	public static String[] getListOfGeneSnps(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> snpResultsCurrentCondition;
		Hashtable<String, Hashtable<String, String[]>> snpResultsCurrentEthnic;
		Hashtable<String, String[]> snpResultsCurrentAnalysis;
		HashSet<String> geneSnpList;
		String[] result = null;

		if (snpResultsAll != null && snpResultsAll.size() > 0) {
			geneSnpList = new HashSet<String>();
			for (String condition : snpResultsAll.keySet()) {
				snpResultsCurrentCondition = snpResultsAll.get(condition);
				for (String ethnic : snpResultsCurrentCondition.keySet()) {
					snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
					for (String analysis : snpResultsCurrentEthnic.keySet()) {	//analysesGroup.keySet()
						snpResultsCurrentAnalysis = snpResultsCurrentEthnic.get(analysis);
						for (String geneSnp : snpResultsCurrentAnalysis.keySet()) {
							if(! geneSnpList.contains(geneSnp)) {
								geneSnpList.add(geneSnp);
							}
						}
					}
				}
			}
			result = geneSnpList.toArray(new String[0]);
			Arrays.sort(result);
		}

		return result;
	}

	public static Hashtable<String, Double> printSnpResults (Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, String[] conditionList, String[] ethnicList, String[] analysesList, String[] geneSnpList, String[] columnNamesOfTheData, String[] columnNamesOfTheDataForUnconditional, String[] columnNamesToOutput, Double threshold, String fullPathOutFilename, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> snpResultsCurrentCondition;
		Hashtable<String, Hashtable<String, String[]>> snpResultsCurrentEthnic;
		Hashtable<String, String[]> snpResultsEachFile;
		PrintWriter writer;
		String[] tmp;
		String line;
		int loop, index;
		double minPvalue;
		Hashtable<String, Double> snpsWithSignificantPval = null;
		int[] indices;

		if (log == null) {
			log = new Logger();
		}

		if (snpResultsAll != null && snpResultsAll.size() > 0) {
			if (ethnicList == null) {
				ethnicList = getEthnicsList(snpResultsAll, log);
			}
			index = ext.indexOfStr("p", columnNamesOfTheData, false, true);
			if (index < 0) {
				log.reportError("No column has the name 'p'. Program halted.");
				System.exit(0);
			}
			indices = ext.indexFactors(columnNamesToOutput, columnNamesOfTheData, false, true);
			snpsWithSignificantPval = new Hashtable<String, Double> ();

			try {
				writer = new PrintWriter(new FileOutputStream(fullPathOutFilename));
				line = "snp\tgene\tchr\tposition\tmin_p";
				for (String condition : conditionList) {
					snpResultsCurrentCondition = snpResultsAll.get(condition);
					for (String item : columnNamesOfTheDataForUnconditional) {
						for (String ethnic : ethnicList) {
							snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
							line += ("\t" + item + "_" + ethnic);
						}
					}
					break;
				}

				for (String condition : conditionList) {
					snpResultsCurrentCondition = snpResultsAll.get(condition);
					for (String ethnic : ethnicList) {
						snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
						for (String analysis : analysesList) {
							for (int i = 0; i < columnNamesToOutput.length; i++) {
								line += "\t" + condition + "_" + ethnic + "_" + analysis + "_" + columnNamesToOutput[i];
							}
						}
					}
				}
				writer.println(line);

				for (String geneSnp : geneSnpList) {
					line = "";
					snpResultsCurrentCondition = snpResultsAll.get("cond1");
					for (int i = 0; i < columnNamesOfTheDataForUnconditional.length; i ++) {
						for (String ethnic : ethnicList) {
							if (snpResultsCurrentCondition.containsKey(ethnic)) {
								snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
								for (String analysis : analysesList) {
									tmp = snpResultsCurrentEthnic.get(analysis + "_totals").get(geneSnp);
									line += ("\t" + tmp[i]);
								}
							}
						}
					}

					minPvalue = Double.MAX_VALUE;
					for (String condition : conditionList) {
						snpResultsCurrentCondition = snpResultsAll.get(condition);
						for (String ethnic : ethnicList) {
							if (snpResultsCurrentCondition.containsKey(ethnic)) {
								snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
								for (String analysis : analysesList) {
									snpResultsEachFile = snpResultsCurrentEthnic.get(analysis);
									tmp = snpResultsEachFile.get(geneSnp);
									if (tmp != null) {
										for (int i = 0; i < indices.length; i++) {
											line += "\t" + tmp[indices[i]];
										}
//										if (tmp[index].matches("-?\\d+(\\.\\d+)?") && Double.parseDouble(tmp[index]) < minPvalue) {
										if (tmp[index] != null && !tmp[index].equals("NA") && Double.parseDouble(tmp[index]) < minPvalue) {
											minPvalue = Double.parseDouble(tmp[index]);
										}
									} else {
										for (int i = 0; i < columnNamesToOutput.length; i++) {
											line += "\t";
										}
									}
								}
							} else {
								loop = analysesList.length * columnNamesToOutput.length;
								for (int i = 0; i < loop; i++) {
									line += "\t";
								}
							}
						}
					}
					writer.println(geneSnp + "\t" + (minPvalue == Double.MAX_VALUE? "NA" : minPvalue) + line);
					if (minPvalue < threshold) {
						snpsWithSignificantPval.put(geneSnp, minPvalue);
					}
				}
				writer.close();
				log.report("snp results is ready at " + fullPathOutFilename);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		return snpsWithSignificantPval;
	}

	public static void printSummarizedSnpResults (Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, Hashtable<String, Double> snpsWithSignificantPval, String[] conditionList, String[] ethnicList, String[] analysesList, String[] columnNamesOfTheData, String fullPathOutFilename, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> snpResultsCurrentCondition;
		Hashtable<String, Hashtable<String, String[]>> snpResultsCurrentEthnic;
		Hashtable<String, String[]> snpResultsEachFile;
		PrintWriter writer;
		String[] tmp, geneSnpList;
		String line;
		int loop;

		if (log == null) {
			log = new Logger();
		}

		if (snpResultsAll != null && snpResultsAll.size() > 0) {
			if (ethnicList == null) {
				ethnicList = getEthnicsList(snpResultsAll, log);
			}

			try {
				writer = new PrintWriter(new FileOutputStream(fullPathOutFilename));
				line = "snp\tgene\tchr\tposition\tmin_p";
				for (String condition : conditionList) {
					snpResultsCurrentCondition = snpResultsAll.get(condition);
					for (String ethnic : ethnicList) {
						snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
						for (String analysis : analysesList) {
							for (int i = 0; i < columnNamesOfTheData.length; i++) {
								line += "\t" + condition + "_" + ethnic + "_" + analysis + "_" + columnNamesOfTheData[i];
							}
						}
					}
				}
				writer.println(line);
				geneSnpList = snpsWithSignificantPval.keySet().toArray(new String[0]);
				Arrays.sort(geneSnpList);
				for (String geneSnp : geneSnpList) {
					line = "";
					for (String condition : conditionList) {
						snpResultsCurrentCondition = snpResultsAll.get(condition);
						for (String ethnic : ethnicList) {
							if (snpResultsCurrentCondition.containsKey(ethnic)) {
								snpResultsCurrentEthnic = snpResultsCurrentCondition.get(ethnic);
								for (String analysis : analysesList) {
									snpResultsEachFile = snpResultsCurrentEthnic.get(analysis);
									tmp = snpResultsEachFile.get(geneSnp);
									if (tmp != null) {
										for (int i = 0; i < tmp.length; i++) {
											line += "\t" + tmp[i];
										}
									} else {
										for (int i = 0; i < columnNamesOfTheData.length; i++) {
											line += "\t";
										}
									}
								}
							} else {
								loop = analysesList.length * columnNamesOfTheData.length;
								for (int i = 0; i < loop; i++) {
									line += "\t";
								}
							}
						}
					}
					writer.println(geneSnp + "\t" + snpsWithSignificantPval.get(geneSnp) + line);
				}
				writer.close();
				log.report("Summary for snp p-values is ready at " + fullPathOutFilename);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> groupFileNames(String[] filenames, int[] columnIndeciesOfPhenoConditionEthnicAnalysis, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroup = null;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		String[] filenameRoot;

		if (log == null) {
			log = new Logger();
		}

		if (filenames != null && filenames.length > 0) {
			phenoGroup = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>();
			for (int i = 0; i < filenames.length; i++) {
				filenameRoot = ext.rootOf(filenames[i]).split("_");
				if (phenoGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[0]])) {
					conditionGroup = phenoGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[0]]);
				} else {
					conditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String>>>();
					phenoGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[0]], conditionGroup);
				}

				if (conditionGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[1]])) {
					ethnicGroup = conditionGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[1]]);
				} else {
					ethnicGroup = new Hashtable<String, Hashtable<String, String>>();
					conditionGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[1]], ethnicGroup);
				}

				if (ethnicGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[2]])) {
					analysesGroup = ethnicGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[2]]);
				} else {
					analysesGroup = new Hashtable<String, String>();
					ethnicGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[2]], analysesGroup);
				}

//				if (filenameRoot.length > 4) {
//					filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]] = filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]] + "_" + filenameRoot[4];
//				}

				if (analysesGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]])) {
					log.reportError("Error - " + filenames[i] + " get duplicated with " + analysesGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]]) + ". \nSystem halted.");
					System.exit(0);
				} else {
					analysesGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]], filenames[i]);
				}
			}
		}

		return phenoGroup;
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> groupFileNames(String[] filenames, int[] columnIndeciesOfPhenoConditionEthnicAnalysis, String[] phenos, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroup = null;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		String[] filenameRoot;
		String tmp;
		int index;
		boolean found;

		if (log == null) {
			log = new Logger();
		}

		if (filenames != null && filenames.length > 0) {
			for (int i = 1; i < columnIndeciesOfPhenoConditionEthnicAnalysis.length; i++) {
				if (columnIndeciesOfPhenoConditionEthnicAnalysis[i] > columnIndeciesOfPhenoConditionEthnicAnalysis[0]) {
					columnIndeciesOfPhenoConditionEthnicAnalysis[i] --;
				}
			}
			phenoGroup = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>();
			for (int i = 0; i < filenames.length; i++) {
				tmp = ext.rootOf(filenames[i]);
				index = -1;
				for (int j = 0; j < phenos.length; j++) {
					if (tmp.contains("_" + phenos[j] + "_") || tmp.endsWith("_" + phenos[j])) {
						tmp = tmp.replace("_" + phenos[j], "");
						index = j;
					} else if (tmp.startsWith(phenos[j] + "_")) {
						tmp = tmp.replace(phenos[j] + "_", "");
						index = j;
					}
					if (index >= 0) {
						break;
					}
				}
				if (index < 0) {
					log.reportError("Error - the following file name does not contain any pheno from the list: " + Array.toStr(phenos) + "\n" + filenames[i] + "\nSystem halted due to error.");
					System.exit(0);
				}
				
				filenameRoot = tmp.split("_");
				if (phenoGroup.containsKey(phenos[index])) {
					conditionGroup = phenoGroup.get(phenos[index]);
				} else {
					conditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String>>>();
					phenoGroup.put(phenos[index], conditionGroup);
				}

				if (conditionGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[1]])) {
					ethnicGroup = conditionGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[1]]);
				} else {
					ethnicGroup = new Hashtable<String, Hashtable<String, String>>();
					conditionGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[1]], ethnicGroup);
				}

				if (ethnicGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[2]])) {
					analysesGroup = ethnicGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[2]]);
				} else {
					analysesGroup = new Hashtable<String, String>();
					ethnicGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[2]], analysesGroup);
				}

				if (analysesGroup.containsKey(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]])) {
					log.reportError("Error - " + filenames[i] + " get duplicated with " + analysesGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]]) + ". \nSystem halted.");
					System.exit(0);
				} else {
					analysesGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]], filenames[i]);
				}
			}
		}

		return phenoGroup;
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> groupFileNames(String fileDirNameTemplate, String[] phenos, String[] ethnics, String[] analyses, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroup = null;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		String filename, filename1, filename2, filename3;
		String condition, errorMessage, tmp;
		int conditionId;
		boolean hasMoreCondition;

		if (log == null) {
			log = new Logger();
		}

		errorMessage = "";
		phenoGroup = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>();
		for (String pheno : phenos) {
			filename = fileDirNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, pheno);
			for (String ethnic : ethnics) {
				filename1 = filename.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnic);
				hasMoreCondition = true;
				conditionId = 0;
				while (hasMoreCondition) {
					conditionId ++;
					condition = "cond" + conditionId;
					filename2 = filename1.replaceAll(FILENAME_CONDITION_SEGMENT, condition);
					hasMoreCondition = false;
					tmp = "";
					for (String analysis : analyses) {
						filename3 = filename2.replaceAll(FILENAME_ANALYSIS_SEGMENT, analysis);
						if (new File(filename3).exists()) {
							hasMoreCondition = true;

							if (phenoGroup.containsKey(pheno)) {
								conditionGroup = phenoGroup.get(pheno);
							} else {
								conditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String>>>();
								phenoGroup.put(pheno, conditionGroup);
							}

							if (conditionGroup.containsKey(condition)) {
								ethnicGroup = conditionGroup.get(condition);
							} else {
								ethnicGroup = new Hashtable<String, Hashtable<String, String>>();
								conditionGroup.put(condition, ethnicGroup);
							}

							if (ethnicGroup.containsKey(ethnic)) {
								analysesGroup = ethnicGroup.get(ethnic);
							} else {
								analysesGroup = new Hashtable<String, String>();
								ethnicGroup.put(ethnic, analysesGroup);
							}

							if (analysesGroup.containsKey(analysis)) {
								log.reportError("Error - when summaring results, the following two files' names get duplidated pheno, condition number, ethinc, and analysis type:\n  " + filename3 + "\n  " + analysesGroup.get(analysis) + "\n.System halted due to error.");
								System.exit(1);
							} else {
								analysesGroup.put(analysis, filename3.substring(filename3.lastIndexOf("/")));
							}

						} else {
							tmp += ("\n  " + filename3);
						}
					}
					if (hasMoreCondition) {
						errorMessage += tmp;
					}
				}
			}
		}
		
		if (! errorMessage.equals("")) {
			log.report("Warning - when summerizing results, the following file(s) are not found:" + errorMessage);
		}

		return phenoGroup;
	}

	public static void getAllRoundsConditions(String fullpathToPreviousRoundResult, String fullpathToPreviousRoundCondition, String fullpathToOutputNextRoundCondition, double pThreshold, Logger log) {
		int a = 1;

		while (getNextRoundCondition("D:/Inflammation/outputs/ARIC_AAEA_LpPLA2_activity_prepCondScores_seqMeta_cond" + a + "_SingleSNP.csv", "D:/Inflammation/tests/activity_cond" + a + ".txt", "D:/Inflammation/tests/activity_cond" + (a + 1) + ".txt", pThreshold, log)) {
			a ++;
		}
	}

	public static void conditionalAnalysisWholeProcessMultiplePhenos (String previousResultFileDirFileameTemplate, String condFileDirFilenameTemplate, String phenoDirFilenameTemplate, String genoDirFileameTemplate, String snpInfoDirFilenameTemplate, String[] phenos, String[] ethnics, String rScriptDir, String resultSummariesDir, String rcommand, double pThresholdLower, double pThresholdHigher, int regionSearchDistance, Logger log) {
		int currentCondition, startCondition = 1;
		String condFileDirAndNameTemplateAfterPhenoFilled, phenoDirAndNameTemplateAfterPhenoFilled, snpInfoFile = null, resultFile, resultsDir, condFile, allEthnics;
		String[] fromPreConditionResults;
		boolean isInitialConditionReady;

		if (! previousResultFileDirFileameTemplate.contains(FILENAME_CONDITION_SEGMENT)) {
			log.reportError("Error - result file name template does not contain a condition segment: " + previousResultFileDirFileameTemplate + "\nNote: even the existing file(s) do not have condition segment, please still specify the segment in the file name template for next rounds' files.");
			return;
		}
		allEthnics = Array.toStr(ethnics, "");

		for (int i = 0; i < phenos.length; i++) {
			log.report("Processing pheno " + phenos[i]);
			isInitialConditionReady = false;
			phenoDirAndNameTemplateAfterPhenoFilled = phenoDirFilenameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
			condFileDirAndNameTemplateAfterPhenoFilled = condFileDirFilenameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
			resultsDir = previousResultFileDirFileameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);

			resultFile = resultsDir.replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleVariant");
			resultFile = resultFile.replaceAll("_" + FILENAME_CONDITION_SEGMENT, "").replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, ""); //TODO
			if (condFileDirFilenameTemplate.contains(FILENAME_ETHNIC_SEGMENT)) {
				for (int j = 0; j < ethnics.length; j++) {	//TODO
					resultFile = resultFile.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
					condFile = condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
					if (! new File(condFile.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1")).exists()) {
						condFile = condFile.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1");
						fromPreConditionResults = getFirstCondition(resultFile, new String[] {"SingleVariant", "Chr", "Position"}, new String[] {"SKATgene", "ARIC_Whites_p_SingleSNP", "ARIC_Blacks_p_SingleSNP"}, new String[] {"(ARIC_Whites_p_SingleSNP || ARIC_Blacks_p_SingleSNP) <= 0.0000001"}, new String[] {"(ARIC_Whites_p_SingleSNP || ARIC_Blacks_p_SingleSNP) <= 0.00001", "(ARIC_Whites_maf_SingleSNP || ARIC_Blacks_maf_SingleSNP) > 0"}, regionSearchDistance, log);
						if (fromPreConditionResults != null) {
							Files.writeList(fromPreConditionResults, condFile);
							isInitialConditionReady = true;
							log.report("No existing initial condition found. Has developed new initial condition from:\n  " + resultFile);
						}
					} else {
						log.report("Found existing initial condition and will skip developing new initial condition:\n  " + condFile.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1"));
						isInitialConditionReady = true;
					}
				}
			} else {
				resultFile = resultFile.replaceAll("_" + FILENAME_ETHNIC_SEGMENT, "");
				if (! new File(condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1")).exists()) {
					fromPreConditionResults = getFirstCondition(resultFile, new String[] {"SingleVariant", "Chr", "Position"}, new String[] {"SKATgene", "ARIC_Whites_p_SingleSNP", "ARIC_Blacks_p_SingleSNP"}, new String[] {"(ARIC_Whites_p_SingleSNP || ARIC_Blacks_p_SingleSNP) <= 0.0000001"}, new String[] {"(ARIC_Whites_p_SingleSNP || ARIC_Blacks_p_SingleSNP) <= 0.00001", "(ARIC_Whites_maf_SingleSNP || ARIC_Blacks_maf_SingleSNP) > 0"}, regionSearchDistance, log);
					if (fromPreConditionResults != null) {
						Files.writeList(fromPreConditionResults, condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1"));
						isInitialConditionReady = true;
						log.report("No existing initial condition found. Has developed new initial condition from:\n  " + resultFile);
					}
				} else {
					log.report("Found existing initial condition and will skip developing new initial condition:\n  " + condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1"));
					isInitialConditionReady = true;
				}
			}

			if (isInitialConditionReady) {
				startCondition = 1;
	
				if (condFileDirFilenameTemplate.contains(FILENAME_ETHNIC_SEGMENT)) {
					for (int j = 0; j < ethnics.length; j++) {
						condFileDirAndNameTemplateAfterPhenoFilled = condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
						currentCondition = startCondition;
						while (runConditionalAnalysisWholeProcessOfOnePheno(condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + currentCondition + ""), phenoDirFilenameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]), genoDirFileameTemplate, snpInfoDirFilenameTemplate, new String[] {ethnics[j]}, rScriptDir, resultsDir.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + currentCondition), rcommand, pThresholdHigher, condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + (currentCondition + 1)), log)) {
							log.report("Finished condition " + currentCondition);
							currentCondition ++;
						}
					}
				} else {
					currentCondition = startCondition;
					while (runConditionalAnalysisWholeProcessOfOnePheno(condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + currentCondition + ""), phenoDirFilenameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]), genoDirFileameTemplate, snpInfoDirFilenameTemplate, ethnics, rScriptDir, resultsDir.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + currentCondition), rcommand, pThresholdHigher, condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + (currentCondition + 1)), log)) {
						log.report("Finished condition " + currentCondition);
						currentCondition ++;
					}
				}
			} else {
				log.reportError("Error - No input is available to initialize conditional alaysis for pheno " + phenos[i] + ". Need either one of the following files:\n  " + resultFile.replaceAll("_" + FILENAME_ETHNIC_SEGMENT, "") + "\n  " + condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond1"));
			}
		}

		summary(previousResultFileDirFileameTemplate, pThresholdHigher, phenos, ethnics, snpInfoDirFilenameTemplate, condFileDirFilenameTemplate, resultSummariesDir, log);
	}

	/* This is a working copy before the last modification of the method with the same name
	public static void conditionalAnalysisWholeProcessMultiplePhenos (String condFileDirAndNameTemplate, String phenoDirAndNameTemplate, String genoDirAndNameTemplate, String snpInfoDirAndNameTemplate, String[] phenos, String[] ethnics, int startCondition, String rScriptDir, String resultsDir, String resultSummariesDir, String rcommand, double pThreshold, Logger log) {
		String condFileDirAndNameTemplateAfterPhenoFilled, phenoDirAndNameTemplateAfterPhenoFilled, snpInfoFile = null, outputFile;
		int currentCondition;

		for (int i = 0; i < phenos.length; i++) {
			phenoDirAndNameTemplateAfterPhenoFilled = phenoDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
			condFileDirAndNameTemplateAfterPhenoFilled = condFileDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]);
			if (condFileDirAndNameTemplate.contains(FILENAME_ETHNIC_SEGMENT)) {
				for (int j = 0; j < ethnics.length; j++) {
					condFileDirAndNameTemplateAfterPhenoFilled = condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]);
					currentCondition = startCondition;
					while (conditionalAnalysisWholeProcessOfOnePheno(condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + currentCondition), null, phenoDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]), genoDirAndNameTemplate, snpInfoDirAndNameTemplate, new String[] {ethnics[j]}, rScriptDir, resultsDir, rcommand, pThreshold, condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + (currentCondition + 1)), log)) {
						currentCondition ++;
					}
				}
			} else {
				currentCondition = startCondition;
				while (conditionalAnalysisWholeProcessOfOnePheno(condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + currentCondition), null, phenoDirAndNameTemplate.replaceAll(FILENAME_PHENO_SEGMENT, phenos[i]), genoDirAndNameTemplate, snpInfoDirAndNameTemplate, ethnics, rScriptDir, resultsDir, rcommand, pThreshold, condFileDirAndNameTemplateAfterPhenoFilled.replaceAll(FILENAME_CONDITION_SEGMENT, "cond" + (currentCondition + 1)), log)) {
					currentCondition ++;
				}
			}
			summary(resultsDir, new int[] {0, 1, 2, 3}, pThreshold, null, null, null, resultSummariesDir, log);
		}
	}
	*/

	/* This is working copy before the 2nd draft of the method with the same name
	public static boolean developConditions (String fullpathToPreviousRoundResult, String fullpathToPreviousRoundCondition, String fullpathToOutputNextRoundCondition, double pThreshold, Logger log) {
		Hashtable <String, String[]> previousRoundResult, previousCondition;
		Vector <String> geneList = null, snpList = null, uniqueGenes = null, uniqueSnps = null, output;
		Hashtable <Integer, Double> significantSnpOfEachRegion_p;
		Hashtable <Integer, String> significantSnpOfEachRegion_snp;
		String[] line;
		String tmp, chr;
		int regionId;
		int[] regionIdsOrdered;
		Hashtable <Integer, Vector<String>> regionToGene, regionToSnp;

		if (log == null) {
			log = new Logger();
		}
		previousCondition = loadFile(fullpathToPreviousRoundCondition, null, new String[] {"SKATgene", "SNP"}, new String[] {"CHROM"}, null, null);
		uniqueGenes = new Vector<String> ();
		uniqueSnps = new Vector<String> ();
		for (String geneSnp : previousCondition.keySet()) {
			line = geneSnp.split("\t");
			if(! uniqueGenes.contains(line[0])) {
				uniqueGenes.add(line[0]);
			}
			if(! uniqueSnps.contains(line[1])) {
				uniqueSnps.add(line[1]);
			}
		}
		regionToGene = new Hashtable <Integer, Vector<String>>();
		regionToSnp = new Hashtable <Integer, Vector<String>>();
		for (String snp : uniqueSnps) {
			for (String gene : uniqueGenes) {
				if (previousCondition.containsKey(gene + "\t" + snp)) {
					regionId = -1;
					regionId = getKeyOfValue(regionToSnp, snp, log);
					if (regionId == -1) {
						regionId = getKeyOfValue(regionToGene, gene, log);
					}
		
					if (regionId > -1) {
						snpList = regionToSnp.get(regionId);
						geneList = regionToGene.get(regionId);
					} else {
						snpList = new Vector<String>();
						geneList = new Vector<String>();
						regionId = regionToSnp.size();
						regionToSnp.put(regionId, snpList);
						regionToGene.put(regionId, geneList);
					}
		
					if (! snpList.contains(snp)) {
						snpList.add(snp);
					}
					if (! geneList.contains(gene)) {
						geneList.add(gene);
					}
				}
			}
		}
		
		previousRoundResult = loadFile(fullpathToPreviousRoundResult, null, new String[] {"gene", "Name"}, new String[] {"p"}, new String[] {"p"}, new String[] {"<=" + pThreshold});
		if (previousRoundResult.size() > 0) {
			significantSnpOfEachRegion_p = new Hashtable <Integer, Double> ();
			significantSnpOfEachRegion_snp = new Hashtable <Integer, String> ();
			for (String geneSnp: previousRoundResult.keySet()) {
				line = geneSnp.split("\t");
				regionId = getKeyOfValue(regionToGene, line[0], log);
				if (regionId == -1) {
					log.reportError("Error - the Snp (" + line[0] + ", " + line[1] + ") from the following file cannot match the regions in the condition\n" + fullpathToPreviousRoundResult);
					System.exit(0);
				}
				if (! significantSnpOfEachRegion_p.containsKey(regionId) || significantSnpOfEachRegion_p.get(regionId) > Double.parseDouble(previousRoundResult.get(geneSnp)[0])) {
					significantSnpOfEachRegion_p.put(regionId, Double.parseDouble(previousRoundResult.get(geneSnp)[0]));
					significantSnpOfEachRegion_snp.put(regionId, line[1]);
				}
			}
	
			regionIdsOrdered = new int[significantSnpOfEachRegion_snp.size()];
			regionId = 0;
			for (int region: significantSnpOfEachRegion_snp.keySet()) {
				regionIdsOrdered[regionId] = region;
				regionId ++;
			}
			Sort.putInOrder(regionIdsOrdered);
	
			output = new Vector<String>();
			output.add("SNP\tSKATgene\tCHROM");
			for (int region : regionIdsOrdered) {
				snpList = regionToSnp.get(region);
				geneList = regionToGene.get(region);
				for (String snp : snpList) {
					for (String gene : geneList) {
						tmp = gene + "\t" + snp;
						if (previousCondition.containsKey(tmp)) {
							output.add(snp + "\t" + gene + "\t" + previousCondition.get(tmp)[0]);
						}
					}
				}
	
				tmp = significantSnpOfEachRegion_snp.get(region);
				for (String gene : geneList) {
					chr = "";
					for (String geneSnp : previousCondition.keySet()) {
						if (geneSnp.split("\t")[0].equalsIgnoreCase(gene)) {
							chr = previousCondition.get(geneSnp)[0];
							break;
						}
					}
					output.add(tmp + "\t" + gene + "\t" + chr);
				}
			}
	
			Files.writeList(output.toArray(new String[0]), fullpathToOutputNextRoundCondition);
			return true;
		} else {
			log.report("No marker is found with significan p-value, so no condition file is generated at " + fullpathToOutputNextRoundCondition);
			return false;
		}
	}
	*/

	public static boolean conditionalAnalysisWholeProcessOfOnePheno_20150109 (String condFileFullPath, String phenoDirAndNameTemplate, String genoDirAndNameTemplate, String snpInfoDirAndNameTemplate, String[] ethnics, String rScriptDir, String resultsDirFilenameTemplate, String rcommand, double pThreshold, String nextCondFileFullPath, Logger log) {
		Hashtable <String, String[]> previousCondition;
		Hashtable <Integer, Vector<String>> regionToGenes, regionToSnps;
		Hashtable <Integer, String> significantSnpOfEachRegion;
//		Hashtable <Integer, String> regionToChr;
		String[] chrList, tmp;
		Object[] test;
		Vector<String> list, chrsWithoutFiles;
		String snpInfoFullPath, genoFullPath, rScriptFullPath, allEthnics, byChrResultDirFilenameTemplate, resultDirFilenameByChr, resultDirFileNameByChrWithEthicAndAnalysis;
//		String chr;
		boolean result, isFileAvailable, found;
		int a;
		int[] b;

		if (log == null) {
			log = new Logger();
		}

//		test = organizeConditionsIntoSegments(condFileFullPath, log);
		test = organizeConditionsIntoRegions(condFileFullPath, log);
		previousCondition = (Hashtable<String, String[]>) test[0];
		regionToGenes = (Hashtable<Integer, Vector<String>>) test[1];
		regionToSnps = (Hashtable<Integer, Vector<String>>) test[2];
//		regionToChr = (Hashtable<Integer, String>) test[3];
		chrList = ((Vector<String>) test[4]).toArray(new String[0]);

		significantSnpOfEachRegion = new Hashtable <Integer, String>();
		a = resultsDirFilenameTemplate.lastIndexOf("/");
		byChrResultDirFilenameTemplate = resultsDirFilenameTemplate.substring(0, a) + "/byChr/";
		if (! new File(byChrResultDirFilenameTemplate).exists()) {
			new File(byChrResultDirFilenameTemplate).mkdir();
		}
		byChrResultDirFilenameTemplate += resultsDirFilenameTemplate.substring(a + 1);

		b = new int[2];
		chrsWithoutFiles = new Vector<String>();
		allEthnics = Array.toStr(ethnics, "");
		for (String chr : chrList) {
			isFileAvailable = true;
			resultDirFilenameByChr = byChrResultDirFilenameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
			resultDirFileNameByChrWithEthicAndAnalysis = resultDirFilenameByChr.replaceAll(FILENAME_ETHNIC_SEGMENT, allEthnics).replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP");
			if (! new File(resultDirFileNameByChrWithEthicAndAnalysis).exists()) {
				rScriptFullPath = rScriptDir + ext.rootOf(condFileFullPath) + "_chr" + chr + ".R";
				b[0] = ext.indexOfStr(chr, SPECIAL_CHROMOSOMES[0]);
				b[1] = ext.indexOfStr(chr, SPECIAL_CHROMOSOMES[1]);
				if (b[0] >= 0) {
					snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
					if (! new File(snpInfoFullPath).exists()) {
						snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[1][b[0]]);
						if (! new File(snpInfoFullPath).exists()) {
							isFileAvailable = false;
							chrsWithoutFiles.add(chr);
							log.reportError("Error - cannot find SNP info file for chr " + chr);
						}
					}
					genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
					if (! new File(genoFullPath).exists()) {
						genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[1][b[0]]);
						if (! new File(genoFullPath).exists()) {
							isFileAvailable = false;
							chrsWithoutFiles.add(chr);
							log.reportError("Error - cannot find geno file for chr " + chr);
						}
					}
				} else if (b[1] >= 0) {
					snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
					if (! new File(snpInfoFullPath).exists()) {
						snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[0][b[1]]);
						if (! new File(snpInfoFullPath).exists()) {
							isFileAvailable = false;
							chrsWithoutFiles.add(chr);
							log.reportError("Error - cannot find SNP info file for chr " + chr);
						}
					}
					genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
					if (! new File(genoFullPath).exists()) {
						genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[0][b[1]]);
						if (! new File(genoFullPath).exists()) {
							isFileAvailable = false;
							chrsWithoutFiles.add(chr);
							log.reportError("Error - cannot find geno file for chr " + chr);
						}
					}
				} else {
					snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
					genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
				}

				if (isFileAvailable) {
					Files.write(getRScriptForConditionalAnalysis(condFileFullPath,
																 snpInfoFullPath,
																 phenoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr),
																 genoFullPath,
																 ethnics,
																 resultDirFilenameByChr,
																 log),
								rScriptFullPath);
					CmdLine.run(rcommand + " " + rScriptFullPath, rScriptDir);
				}
			}
			
			if (new File(resultDirFileNameByChrWithEthicAndAnalysis).exists()) {
				significantSnpOfEachRegion.putAll(getSignificantSnpOfEachRegion(resultDirFilenameByChr.replaceAll(FILENAME_ETHNIC_SEGMENT, allEthnics).replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP"), regionToGenes, pThreshold, log));
			} else if (isFileAvailable) {
				log.reportError("Error - Potential R failure. The following R output is not found:\n" + resultDirFileNameByChrWithEthicAndAnalysis);
			}
		}

		if (significantSnpOfEachRegion.size() > 0) {
			Files.writeList(getNextCondition(significantSnpOfEachRegion, previousCondition, regionToGenes, regionToSnps, log), nextCondFileFullPath);
			result = true;
		} else {
			result = false;
		}

		if (chrsWithoutFiles.size() > 0) {
			tmp = new String[chrList.length - chrsWithoutFiles.size()];
			a = 0;
			for (int i = 0; i < chrList.length; i++) {
				found = false;
				for (int j = 0; j < chrsWithoutFiles.size(); j++) {
					if (chrList[i].equals(chrsWithoutFiles.elementAt(j))) {
						chrsWithoutFiles.remove(j);
						found = true;
						break;
					}
				}
				if (! found) {
					tmp[a] = chrList[i];
					a ++;
				}
			}
			chrList = tmp;
		}

		mergeSkatMetaResultFiles(byChrResultDirFilenameTemplate, Array.addStrToArray(allEthnics, ethnics), chrList, new String[] {"SingleSNP", "T5Count", "T5"}, resultsDirFilenameTemplate.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, ""), log);

		return result;
	}

	public static boolean runConditionalAnalysisWholeProcessOfOnePheno (String condFileFullPath, String phenoDirAndNameTemplate, String genoDirAndNameTemplate, String snpInfoDirAndNameTemplate, String[] ethnics, String rScriptDir, String resultsDirFilenameTemplate, String rcommand, double pThreshold, String nextCondFileFullPath, Logger log) {
		Hashtable <String, String[]> previousCondition;
		Hashtable <Integer, Vector<String>> regionToGenes, regionToSnps;
		Hashtable <Integer, String> significantSnpOfEachRegion;
//		Hashtable <Integer, String> regionToChr;
		String[] chrs, filenames;
		Object[] test;
		Vector<String> chrsToRemove, ethnicsToRemove;
		String snpInfoFile = null, genoFile, rScriptFile = null, allEthnics, byChrResultDirFilenameTemplate, resultDirFilenameTemplateByChr = null, resultFile1, resultFileForNextCond, tmp2;
		String[][] resultFiles;
//		String chr;
		boolean result = false, isAllRDataFilesForTheChromosomeMissing;
		boolean[] isAllRDataFilesForTheEthnicMissing;
		boolean[][] areResultsByChrAvailable;
		int a;
		int[] b;

		if (log == null) {
			log = new Logger();
		}

		test = organizeConditionsIntoRegions(condFileFullPath, log);
		previousCondition = (Hashtable<String, String[]>) test[0];
		regionToGenes = (Hashtable<Integer, Vector<String>>) test[1];
		regionToSnps = (Hashtable<Integer, Vector<String>>) test[2];
//		regionToChr = test[3];
		chrs = ((Vector <String>) test[4]).toArray(new String[0]);

		significantSnpOfEachRegion = new Hashtable <Integer, String>();

		b = new int[2];
		allEthnics = Array.toStr(ethnics, "");
		resultFileForNextCond = resultsDirFilenameTemplate.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, "").replaceAll(FILENAME_ETHNIC_SEGMENT, allEthnics).replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP");
		if (! new File(resultFileForNextCond).exists()) {
			if (genoDirAndNameTemplate.contains(FILENAME_CHROMOSOME_SEGMENT) || snpInfoDirAndNameTemplate.contains(FILENAME_CHROMOSOME_SEGMENT)) {
				resultDirFilenameTemplateByChr = createSubDir(resultsDirFilenameTemplate, "/byChr");
				areResultsByChrAvailable = areFilesAvailable(resultDirFilenameTemplateByChr.replaceAll(FILENAME_ANALYSIS_SEGMENT, "SingleSNP"), chrs, ethnics, log);//TODO need all the analyses
				chrsToRemove = new Vector<String>();
				isAllRDataFilesForTheEthnicMissing = new boolean[ethnics.length];
				Arrays.fill(isAllRDataFilesForTheEthnicMissing, true);
				for (int i = 0; i < chrs.length; i++) {
					isAllRDataFilesForTheChromosomeMissing = true;
					resultFiles = new String[ethnics.length][1];
					for (int j = 0; j < ethnics.length; j++) {
						if (! areResultsByChrAvailable[i][j]) {
							filenames = getFileNames(chrs[i], snpInfoDirAndNameTemplate, genoDirAndNameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]), log);
							snpInfoFile = filenames[0];
							genoFile = filenames[1];
							if (snpInfoFile == null || genoFile == null) {
								chrsToRemove.add(chrs[i]);
							} else {
								resultFiles[j][0] = getFileNameRootFromTemplate(resultDirFilenameTemplateByChr, chrs[i], ethnics[j], log) + ".RData";
								if (! new File(resultFiles[j][0]).exists()) {
									rScriptFile = rScriptDir + ext.rootOf(resultFiles[j][0]) + ".R";
									Files.write(getRScriptForCohort1(snpInfoFile,
																	genoFile,
																	phenoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chrs[i]).replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[j]),
																	condFileFullPath,
																	resultFiles[j][0]),
												rScriptFile);
									CmdLine.run(rcommand + " " + rScriptFile, rScriptDir);
	
									if (! new File(resultFiles[j][0]).exists()) {
										log.report("Warning - Potential failure of R. The following R output is not found and will be ignored:\n  " + resultFiles[j][0] + "\nCheck the following R script for detail:\n  " + rScriptFile);
										resultFiles[j][0] = null;
									} else {
										isAllRDataFilesForTheChromosomeMissing = false;
										isAllRDataFilesForTheEthnicMissing[j] = false;
									}
								} else {
									isAllRDataFilesForTheChromosomeMissing = false;
									isAllRDataFilesForTheEthnicMissing[j] = false;
								}
							}
						}
					}
					
					if (! isAllRDataFilesForTheChromosomeMissing) {
						rScriptFile = rScriptDir + ext.rootOf(getFileNameRootFromTemplate(resultDirFilenameTemplateByChr, chrs[i], allEthnics + "_meta", log)) + ".R";
						Files.write(getRScriptForMetaAnalysis1(snpInfoFile,
															  ethnics,
															  resultFiles,
															  condFileFullPath,
															  resultDirFilenameTemplateByChr.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chrs[i]),
															  log),
									rScriptFile);
						CmdLine.run(rcommand + " " + rScriptFile, rScriptDir);
					} else {
						chrsToRemove.add(chrs[i]);
					}
				}

				chrs = removeElementsFromAarray(chrsToRemove, chrs);
				ethnicsToRemove = new Vector<String> (ethnics.length);
				for (int j = 0; j < ethnics.length; j++) {
					if (isAllRDataFilesForTheEthnicMissing[j]) {
						ethnicsToRemove.add(ethnics[j]);
					}
				}

				if (chrsToRemove.size() != chrs.length) {
					mergeSkatMetaResultFiles(resultDirFilenameTemplateByChr, Array.addStrToArray(allEthnics, removeElementsFromAarray(ethnicsToRemove, ethnics)), chrs, new String[] {"SingleSNP", "T5Count", "T5"}, resultsDirFilenameTemplate.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, ""), log);
				}
			} else {
				//TODO
//				getRScript(chrs[i], ethnics[j], ...);
//				Files.write(getRScriptForConditionalAnalysis(condFileFullPath, snpInfoFullPath, phenoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr), genoFullPath, ethnics, resultDirFilenameByChr, log), rScriptFullPath);
				CmdLine.run(rcommand + " " + rScriptFile, rScriptDir);
				for (int i = 0; i < ethnics.length; i++) {
					if (! new File(resultDirFilenameTemplateByChr.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, "").replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i])).exists()) {
						// exclude chr from chrList; ignore ethnic from merged results
						log.report("Warning - Potential failure of R. The following R output is not found and will be ignored:\n  " + resultDirFilenameTemplateByChr.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, "").replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i]));
					}
				}
			}

		} else {
			log.report("Note - Will skip R step and use the following existing result file to develop the next round condition:\n  " + resultFileForNextCond);
		}

		if (new File(resultFileForNextCond).exists()) {
			significantSnpOfEachRegion = getSignificantSnpOfEachRegion(resultFileForNextCond, regionToGenes, pThreshold, log);
			if (significantSnpOfEachRegion.size() > 0) {
				Files.writeList(getNextCondition(significantSnpOfEachRegion, previousCondition, regionToGenes, regionToSnps, log), nextCondFileFullPath);
				result = true;
			}
		} else {
			log.reportError("Warning - Potential failure of the merging of R outputs of different chromosomes. The following file is not found and will be ignored:\n  " + resultFileForNextCond);
		}

		return result;
	}

	/**
	 * Set chr=null, when you want the results for non-chromosome specific results.
	 * @param resultsDir
	 * @param currentCondition
	 * @param phenos
	 * @param chr
	 * @param ethnics
	 * @param log
	 * @return
	 */
	public static boolean[][][] areFilesAvailable (String resultsDir, String[] phenos, String[] chrs, String[] ethnics, Logger log) {
		boolean[][][] result;
		result = new boolean[phenos.length][][];

		for (int i = 0; i < phenos.length; i++) {
			result[i] = areFilesAvailable(resultsDir, chrs, ethnics, log);
		}

		return result;
	}

	/**
	 * Set chr=null, when you want the results for non-chromosome specific results.
	 * @param resultsDir
	 * @param currentCondition
	 * @param phenos
	 * @param chrs
	 * @param ethnics
	 * @param log
	 * @return
	 */
	public static boolean[][] areFilesAvailable (String resultsDir, String[] chrs, String[] ethnics, Logger log) {
		boolean[][] result;
		result = new boolean[chrs.length][];

		for (int i = 0; i < chrs.length; i++) {
			result[i] = areFilesAvailable(resultsDir, chrs[i], ethnics, log);
		}

		return result;
	}

	/**
	 * Set chr=null, when you want the results for non-chromosome specific files.
	 * @param filesDirAndNameTemplate should only have %chr% and/or %ethnic% segment(s).
	 * @param currentCondition
	 * @param chr
	 * @param ethnics
	 * @param log
	 * @return
	 */
	public static boolean[] areFilesAvailable (String filesDirAndNameTemplate, String chr, String[] ethnics, Logger log) {
		boolean[] result;
		String file;

		result = new boolean[ethnics.length];
		if (chr == null) {
			file = filesDirAndNameTemplate.replaceAll("_" + FILENAME_CHROMOSOME_SEGMENT, "").replaceAll(FILENAME_CHROMOSOME_SEGMENT, "");	//TODO: need to provide filename segement separator
		} else {
			file = filesDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
		}

		for (int i = 0; i < ethnics.length; i++) {
			if (new File(file.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnics[i])).exists()) {
				result[i] = true;
			}
		}

		return result;
	}

	public static String createSubDir(String fileDirAndNameTemplate, String subDir) {
		int a;
		String tmp;

		if (subDir.startsWith("/")) {
			subDir = subDir.substring(1);
		}
		if (subDir.endsWith("/")) {
			subDir = subDir.substring(0, subDir.length() - 1);
		}
		a = fileDirAndNameTemplate.lastIndexOf("/");
		tmp = fileDirAndNameTemplate.substring(0, a) + "/" + subDir + "/";
		if (! new File(tmp).exists()) {
			new File(tmp).mkdir();
		}
		tmp += fileDirAndNameTemplate.substring(a + 1);
		
		return tmp;
	}

	public static String[] getFileNames (String chr, String snpInfoDirAndNameTemplate, String genoDirAndNameTemplate, Logger log) {
		String snpInfoFullPath, genoFullPath;
		int[] indices;

		if (log == null) {
			log = new Logger();
		}

		indices = new int[2];
		indices[0] = ext.indexOfStr(chr, SPECIAL_CHROMOSOMES[0]);
		indices[1] = ext.indexOfStr(chr, SPECIAL_CHROMOSOMES[1]);

		if (indices[0] >= 0) {
			snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
			if (! new File(snpInfoFullPath).exists()) {
				snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[1][indices[0]]);
				if (! new File(snpInfoFullPath).exists()) {
					snpInfoFullPath = null;
					log.reportError("Error - cannot find SNP info file for chr " + chr);
				}
			}
			genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
			if (! new File(genoFullPath).exists()) {
				genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[1][indices[0]]);
				if (! new File(genoFullPath).exists()) {
					genoFullPath = null;
					log.reportError("Error - cannot find geno file for chr " + chr);
				}
			}
		} else if (indices[1] >= 0) {
			snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
			if (! new File(snpInfoFullPath).exists()) {
				snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[0][indices[1]]);
				if (! new File(snpInfoFullPath).exists()) {
					snpInfoFullPath = null;
					log.reportError("Error - cannot find SNP info file for chr " + chr);
				}
			}
			genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
			if (! new File(genoFullPath).exists()) {
				genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + SPECIAL_CHROMOSOMES[0][indices[1]]);
				if (! new File(genoFullPath).exists()) {
					genoFullPath = null;
					log.reportError("Error - cannot find geno file for chr " + chr);
				}
			}
		} else {
			snpInfoFullPath = snpInfoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
			genoFullPath = genoDirAndNameTemplate.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
		}

		return new String[] {snpInfoFullPath, genoFullPath};
	}

	public static String getFileNameRootFromTemplate (String resultDirFilenameTemplateByChr, String chr, String ethnic, Logger log) {
		String resultFile;
		int a;

		resultFile = resultDirFilenameTemplateByChr.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr).replaceAll(FILENAME_ETHNIC_SEGMENT, ethnic).replaceAll("_" + FILENAME_ANALYSIS_SEGMENT, "");
		a = resultFile.lastIndexOf(".");
		if (a >= 0 && !resultFile.endsWith("/")) {
			resultFile = resultFile.substring(0, a);
		}

		return resultFile;
	}

	public static String[] removeElementsFromAarray (Vector<String> elementsOfTheArrayToBeRemoved, String[] array) {
		String[] tmp;
		boolean found;
		int a;
		
		tmp = new String[array.length - elementsOfTheArrayToBeRemoved.size()];
		a = 0;
		for (int i = 0; i < array.length; i++) {
			found = false;
			for (int j = 0; j < elementsOfTheArrayToBeRemoved.size(); j++) {
				if (array[i].equals(elementsOfTheArrayToBeRemoved.elementAt(j))) {
					elementsOfTheArrayToBeRemoved.remove(j);
					found = true;
					break;
				}
			}
			if (! found) {
				tmp[a] = array[i];
				a ++;
			}
		}

		return tmp;
	}

	public static void mergeSkatMetaResultFiles (String byChrResultsDirFilenameTemplate, String[] ethnics, String[] chrList, String[] analyses, String mergedDirFilenameTemplate, Logger log) {
		Vector<String> list;
		String inFilename, inFilenameAfterEthnic, inFilenameAfterEthnicAnalysis, allEthnics, outFilename1, outFilename2;
		boolean toIncludeHeaderLine;

		for (String ethnic : ethnics) {
			inFilenameAfterEthnic = byChrResultsDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnic);
			outFilename1 = mergedDirFilenameTemplate.replaceAll(FILENAME_ETHNIC_SEGMENT, ethnic);
			for (String analysis : analyses) {
				inFilenameAfterEthnicAnalysis = inFilenameAfterEthnic.replaceAll(FILENAME_ANALYSIS_SEGMENT, analysis);
				outFilename2 = outFilename1.replaceAll(FILENAME_ANALYSIS_SEGMENT, analysis);
				list = new Vector<String>();
				toIncludeHeaderLine = true;
				for (String chr : chrList) {
					inFilename = inFilenameAfterEthnicAnalysis.replaceAll(FILENAME_CHROMOSOME_SEGMENT, "chr" + chr);
					if (! new File(inFilename).exists()) {
						//TODO This check might be duplicated with the parent code. Could add boolean[][] areFilesAvailable as an input, so as to avoid reporting two warning for the same error.
						log.reportError("Note - when merging SkatMeta result files of different chromosomes, the following file is not found and will be ignored:\n  " + inFilename);
					} else {
						list.addAll(loadFile(inFilename, (toIncludeHeaderLine? 0 : 1)));
						toIncludeHeaderLine = false;
					}
				}
				if (list.size() > 0) {
					Files.writeList(list.toArray(new String[0]), outFilename2);
				}
			}
		}
	}

	public static boolean getNextRoundCondition(String fullpathToPreviousRoundResult, String fullpathToPreviousRoundCondition, String fullpathToOutputNextRoundCondition, double pThreshold, Logger log) {
		Hashtable <String, String[]> previousCondition;
		Hashtable <Integer, Vector<String>> regionToGenes, regionToSnps;
		Object[] test;
		String[] newCondition;
		String newConditionDirAndName;

		if (log == null) {
			log = new Logger();
		}

		test = organizeConditionsIntoRegions(fullpathToPreviousRoundCondition, log);
		previousCondition = (Hashtable<String, String[]>) test[0];
		regionToGenes = (Hashtable<Integer, Vector<String>>) test[1];
		regionToSnps = (Hashtable<Integer, Vector<String>>) test[2];

		newCondition = getNextCondition(fullpathToPreviousRoundResult, previousCondition, regionToGenes, regionToSnps, pThreshold, log);
		if (newCondition == null) {
			return false;
		} else {
			Files.writeList(newCondition, fullpathToOutputNextRoundCondition);
			return true;
		}
	}

	public static Hashtable[] organizeConditionsIntoRegions_20150109 (String input_fullpathToPreviousRoundCondition, Logger log) {
		Vector <String> geneList = null, snpList = null, uniqueGenes = null, uniqueSnps = null;
		Hashtable <Integer, Vector<String>> regionToGenes, regionToSnps;
		Hashtable <Integer, String> regionToChr; 
		Hashtable <String, String> chrList; 
		Hashtable <String, String[]> condition;
		String[] line;
		String chr;
		int regionId;

		if (log == null) {
			log = new Logger();
		}

		condition = loadFile(input_fullpathToPreviousRoundCondition, null, new String[] {"SKATgene", "SNP"}, new String[] {"Chr"}, null, null);
		uniqueGenes = new Vector<String> ();
		uniqueSnps = new Vector<String> ();
		for (String geneSnp : condition.keySet()) {
			line = geneSnp.split("\t");
			if(! uniqueGenes.contains(line[0])) {
				uniqueGenes.add(line[0]);
			}
			if(! uniqueSnps.contains(line[1])) {
				uniqueSnps.add(line[1]);
			}
		}

		regionToGenes = new Hashtable <Integer, Vector<String>>();
		regionToSnps = new Hashtable <Integer, Vector<String>>();
		regionToChr = new Hashtable <Integer, String>();
		chrList = new Hashtable <String, String>();
		for (String snp : uniqueSnps) {
			for (String gene : uniqueGenes) {
				if (condition.containsKey(gene + "\t" + snp)) {
					chr = condition.get(gene + "\t" + snp)[0];
					regionId = -1;
					regionId = getKeyOfValue(regionToSnps, snp, log);
					if (regionId == -1) {
						regionId = getKeyOfValue(regionToGenes, gene, log);
					}
		
					if (regionId > -1) {
						snpList = regionToSnps.get(regionId);
						geneList = regionToGenes.get(regionId);
						if (! regionToChr.get(regionId).equalsIgnoreCase(chr)) {
							log.reportError("Condition file error: the same region that includes (" + gene + " " + snp + ") appears in two different chromosomes: " + chr + " and " + regionToChr.get(regionId) + ".");
							System.exit(1);
						}
					} else {
						regionId = regionToSnps.size();
						snpList = new Vector<String>();
						regionToSnps.put(regionId, snpList);
						geneList = new Vector<String>();
						regionToGenes.put(regionId, geneList);
						regionToChr.put(regionId, chr);
						chrList.put(chr, "");
					}
		
					if (! snpList.contains(snp)) {
						snpList.add(snp);
					}
					if (! geneList.contains(gene)) {
						geneList.add(gene);
					}
				}
			}
		}

		return new Hashtable[] {condition, regionToGenes, regionToSnps, regionToChr, chrList};
	}

	public static Object[] organizeConditionsIntoRegions (String input_fullpathToPreviousRoundCondition, Logger log) {
		Vector <String> geneList = null, snpList = null, uniqueGenes = null, uniqueSnps = null;
		Hashtable <Integer, Vector<String>> regionToGenes, regionToSnps;
		Hashtable <Integer, String> regionToChr; 
		Vector <String> chrList; 
		Hashtable <String, String[]> condition;
		String[] line;
		String chr;
		int regionId;

		if (log == null) {
			log = new Logger();
		}

		condition = loadFile(input_fullpathToPreviousRoundCondition, null, new String[] {"SKATgene", "SNP"}, new String[] {"Chr", "Pos1", "Pos2"}, null, null);
		uniqueGenes = new Vector<String> ();
		uniqueSnps = new Vector<String> ();
		for (String geneSnp : condition.keySet()) {
			line = geneSnp.split("\t");
			if(! uniqueGenes.contains(line[0])) {
				uniqueGenes.add(line[0]);
			}
			if(! uniqueSnps.contains(line[1])) {
				uniqueSnps.add(line[1]);
			}
		}

		regionToGenes = new Hashtable <Integer, Vector<String>>();
		regionToSnps = new Hashtable <Integer, Vector<String>>();
		regionToChr = new Hashtable <Integer, String>();
		chrList = new Vector <String>();
		for (String snp : uniqueSnps) {
			for (String gene : uniqueGenes) {
				if (condition.containsKey(gene + "\t" + snp)) {
					chr = condition.get(gene + "\t" + snp)[0];
					regionId = -1;
					regionId = getKeyOfValue(regionToSnps, snp, log);
					if (regionId == -1) {
						regionId = getKeyOfValue(regionToGenes, gene, log);
					}
		
					if (regionId > -1) {
						snpList = regionToSnps.get(regionId);
						geneList = regionToGenes.get(regionId);
						if (! regionToChr.get(regionId).equalsIgnoreCase(chr)) {
							log.reportError("Condition file error: the same region that includes (" + gene + " " + snp + ") appears in two different chromosomes: " + chr + " and " + regionToChr.get(regionId) + ".");
							System.exit(1);
						}
					} else {
						regionId = regionToSnps.size();
						snpList = new Vector<String>();
						regionToSnps.put(regionId, snpList);
						geneList = new Vector<String>();
						regionToGenes.put(regionId, geneList);
						regionToChr.put(regionId, chr);
						chrList.add(chr);
					}
		
					if (! snpList.contains(snp)) {
						snpList.add(snp);
					}
					if (! geneList.contains(gene)) {
						geneList.add(gene);
					}
				}
			}
		}

		return new Object[] {condition, regionToGenes, regionToSnps, regionToChr, chrList};
	}

	public static Hashtable[] organizeConditionsIntoSegments (String input_fullpathToPreviousRoundCondition, Logger log) {
		Vector <String> geneList = null;
		Hashtable <String, Vector<String>> snpToGene;//, regionToSnp; 
		Hashtable <String, String[]> condition;
		String[] line;

		if (log == null) {
			log = new Logger();
		}

		condition = loadFile(input_fullpathToPreviousRoundCondition, null, new String[] {"SNP", "SKATgene"}, new String[] {"Chr"}, null, null);
		snpToGene = new Hashtable <String, Vector<String>>();
		for (String geneSnp : condition.keySet()) {
			line = geneSnp.split("\t");
			if (snpToGene.containsKey(line[0])) {
				geneList = snpToGene.get(line[0]);
			} else {
				geneList = new Vector<String>();
				snpToGene.put(line[0], geneList);
			}
			geneList.add(line[1] + "\t" + condition.get(geneSnp)[0]);
		}

		return new Hashtable[] {condition, snpToGene};
	}

	// This is a fully functioning version of the method with the same name. Works when all chromsomes' results are combined into one single result file.
	public static String[] getNextCondition (String fullpathToPreviousRoundResult, Hashtable <String, String[]> previousCondition, Hashtable <Integer, Vector<String>> regionToGenes, Hashtable <Integer, Vector<String>> regionToSnps, double pThreshold, Logger log) {
		Hashtable <String, String[]> previousRoundResult;
		Vector <String> geneList = null, snpList = null, output;
		Hashtable <Integer, Double> significantSnpOfEachRegion_p;
		Hashtable <Integer, String> significantSnpOfEachRegion_snp;
		String[] line;
		String tmp, chr;
		int regionId;
		int[] regionIdsOrdered;
		

		if (log == null) {
			log = new Logger();
		}

		previousRoundResult = loadFile(fullpathToPreviousRoundResult, null, new String[] {"gene", "Name"}, new String[] {"p"}, new String[] {"p<=" + pThreshold}, log);
		if (previousRoundResult.size() > 0) {
			significantSnpOfEachRegion_p = new Hashtable <Integer, Double> ();
			significantSnpOfEachRegion_snp = new Hashtable <Integer, String> ();
			for (String geneSnp: previousRoundResult.keySet()) {
				line = geneSnp.split("\t");
				regionId = getKeyOfValue(regionToGenes, line[0], log);
				if (regionId == -1) {
					log.reportError("Error - the Snp (" + line[0] + ", " + line[1] + ") from the following file cannot match the regions in the condition\n" + fullpathToPreviousRoundResult);
					System.exit(0);
				}
				if (! significantSnpOfEachRegion_p.containsKey(regionId) || significantSnpOfEachRegion_p.get(regionId) > Double.parseDouble(previousRoundResult.get(geneSnp)[0])) {
					significantSnpOfEachRegion_p.put(regionId, Double.parseDouble(previousRoundResult.get(geneSnp)[0]));
					significantSnpOfEachRegion_snp.put(regionId, line[1]);
				}
			}
	
			regionIdsOrdered = new int[significantSnpOfEachRegion_snp.size()];
			regionId = 0;
			for (int region: significantSnpOfEachRegion_snp.keySet()) {
				regionIdsOrdered[regionId] = region;
				regionId ++;
			}
			Sort.putInOrder(regionIdsOrdered);
	
			output = new Vector<String>();
			output.add("SNP\tSKATgene\tCHROM");
			for (int region : regionIdsOrdered) {
				snpList = regionToSnps.get(region);
				geneList = regionToGenes.get(region);
				for (String snp : snpList) {
					for (String gene : geneList) {
						tmp = gene + "\t" + snp;
						if (previousCondition.containsKey(tmp)) {
							output.add(snp + "\t" + gene + "\t" + previousCondition.get(tmp)[0]);
						}
					}
				}
	
				tmp = significantSnpOfEachRegion_snp.get(region);
				for (String gene : geneList) {
					chr = "";
					for (String geneSnp : previousCondition.keySet()) {
						if (geneSnp.split("\t")[0].equalsIgnoreCase(gene)) {
							chr = previousCondition.get(geneSnp)[0];
							break;
						}
					}
					output.add(tmp + "\t" + gene + "\t" + chr);
				}
			}
	
			return output.toArray(new String[0]);
		} else {
			log.report("No marker is found with significan p-value, so no new condition file is generated after " + previousCondition);
			return null;
		}
	}

	public static Hashtable <Integer, String> getSignificantSnpOfEachRegion (String fullpathToPreviousRoundResult, Hashtable <Integer, Vector<String>> regionToGenes, double pThreshold, Logger log) {
		Hashtable <String, String[]> previousRoundResult;
		Vector <String> geneList = null, snpList = null, output;
		Hashtable <Integer, Double> significantSnpOfEachRegion_p;
		Hashtable <Integer, String> significantSnpOfEachRegion_snp;
		String[] line;
		String tmp, chr;
		int regionId;
		int[] regionIdsOrdered;

		if (log == null) {
			log = new Logger();
		}

		significantSnpOfEachRegion_snp = new Hashtable <Integer, String> ();
		previousRoundResult = loadFile(fullpathToPreviousRoundResult, null, new String[] {"gene", "Name"}, new String[] {"p"}, new String[] {"p<=" + pThreshold}, log);
		if (previousRoundResult.size() > 0) {
			significantSnpOfEachRegion_p = new Hashtable <Integer, Double> ();
			for (String geneSnp: previousRoundResult.keySet()) {
				line = geneSnp.split("\t");
				regionId = getKeyOfValue(regionToGenes, line[0], log);
				if (regionId == -1) {
					log.reportError("Error - the Snp (" + line[0] + ", " + line[1] + ") from the following file cannot match the regions in the condition\n" + fullpathToPreviousRoundResult);
					System.exit(0);
				}
				if (! significantSnpOfEachRegion_p.containsKey(regionId) || significantSnpOfEachRegion_p.get(regionId) > Double.parseDouble(previousRoundResult.get(geneSnp)[0])) {
					significantSnpOfEachRegion_p.put(regionId, Double.parseDouble(previousRoundResult.get(geneSnp)[0]));
					significantSnpOfEachRegion_snp.put(regionId, line[1]);
				}
			}
		}

		return significantSnpOfEachRegion_snp;
	}

	public static String[] getNextCondition (Hashtable <Integer, String> significantSnpOfEachRegion, Hashtable <String, String[]> previousCondition, Hashtable <Integer, Vector<String>> regionToGenes, Hashtable <Integer, Vector<String>> regionToSnps, Logger log) {
		Vector <String> geneList = null, snpList = null, output;
		String[] line;
		String tmp, chr;
		int regionId;
		int[] regionIdsOrdered;

		if (log == null) {
			log = new Logger();
		}

		regionIdsOrdered = new int[significantSnpOfEachRegion.size()];
		regionId = 0;
		for (int region: significantSnpOfEachRegion.keySet()) {
			regionIdsOrdered[regionId] = region;
			regionId ++;
		}
		Sort.putInOrder(regionIdsOrdered);

		output = new Vector<String>();
		output.add("SNP\tSKATgene\tChr");

//		for (String geneSnp : previousCondition.keySet()) {
//			line = geneSnp.split("\t");
//			regionId = getKeyOfValue(regionToGenes, line[0], log);
//			if (significantSnpOfEachRegion.containsKey(regionId)) {
//				tmp = line[0] + "\t" + line[1];
//				if (previousCondition.containsKey(tmp)) {
//					output.add(line[1] + "\t" + line[0] + "\t" + previousCondition.get(tmp)[0]);
//				}
//			}
//		}

		for (int region : regionIdsOrdered) {
			snpList = regionToSnps.get(region);
			geneList = regionToGenes.get(region);
			for (String snp : snpList) {
				for (String gene : geneList) {
					tmp = gene + "\t" + snp;
					if (previousCondition.containsKey(tmp)) {
						output.add(snp + "\t" + gene + "\t" + previousCondition.get(tmp)[0]);
					}
				}
			}

			tmp = significantSnpOfEachRegion.get(region);
			for (String gene : geneList) {
				chr = "";
				for (String geneSnp : previousCondition.keySet()) {
					if (geneSnp.split("\t")[0].equalsIgnoreCase(gene)) {
						chr = previousCondition.get(geneSnp)[0];
						break;
					}
				}
				output.add(tmp + "\t" + gene + "\t" + chr);
			}
		}

		return output.toArray(new String[0]);
	}

	/**
	 * 
	 * @param resultFullPath
	 * @param namesOfColumnsOfSnpChrPos new String[] {"SNP_RS_NUmber", "Chr", "Position"}
	 * @param namesOfColumnsOfGeneAndPvalues new String[] {"SKATgene", "minPval"}
	 * @param criteriaForPValTigherThreshold
	 * @param criteriaForPValLooserThreshold
	 * @param regionSearchDistance
	 * @param log
	 * @return
	 */
	public static String[] getFirstCondition (String resultFullPath, String[] namesOfColumnsOfSnpChrPos, String[] namesOfColumnsOfGeneAndPvalues, String[] criteriaForPValTigherThreshold, String[] criteriaForPValLooserThreshold, int regionSearchDistance, Logger log) {
		Hashtable<String, String[]> snpsFromResultFile;
		Hashtable <Integer, String> regionToSnpWithMinP, regionToChr;
		Hashtable <Integer, Vector<String>> regionToGenes;
		Hashtable <Integer, int[]> regionToBoundaries;
		Hashtable <Integer, Double> regionToMinP;
		Hashtable<String, Integer> ttt;
		Vector <String> geneList = null, geneList2, output;
		Vector <Integer> regionsFound;
		String[] line, snpChrPosition, regionIdsOrdered;
		String snp, chr, gene;
		int regionId, position, regionId2;
		int[] boundaries, boundaries2;
		double p;
		Object[] tmp2;
		int[][] criteriaColumnIndices;
		String[] criteriaOperators, criteriaValues;

		if (log == null) {
			log = new Logger();
		}
		if (! new File(resultFullPath).exists()) {
			log.reportError("Error - the following file does not exist:\n  " + resultFullPath);
			return null;
		}

		regionToSnpWithMinP = new Hashtable <Integer, String>();
		regionToGenes = new Hashtable <Integer, Vector<String>>();
		regionToBoundaries = new Hashtable <Integer, int[]>();
		regionToMinP = new Hashtable <Integer, Double>();
		regionToChr = new Hashtable <Integer, String>();
		snpsFromResultFile = loadFile(resultFullPath, null, namesOfColumnsOfSnpChrPos, namesOfColumnsOfGeneAndPvalues, criteriaForPValLooserThreshold, log);

		tmp2 = parseCriteria(namesOfColumnsOfGeneAndPvalues, criteriaForPValTigherThreshold, log);
		criteriaColumnIndices = (int[][]) tmp2[0];
		criteriaOperators = (String[]) tmp2[1];
		criteriaValues = (String[]) tmp2[2];

		regionId = 0;
		for (String snpChrPos: snpsFromResultFile.keySet()) {
			line = snpsFromResultFile.get(snpChrPos);
//			if (Double.parseDouble(line[3]) <= pThresholdTighter) {
			if (screenWithCriteria(line, criteriaColumnIndices, criteriaOperators, criteriaValues, log)) {
				snpChrPosition = snpChrPos.split("\t");
//				chr = snpChrPosition[1];
				position = Integer.parseInt(snpChrPosition[2]);
				gene = line[0];
				p = 1;
				for (int i = 1; i < line.length; i++) {
					try {
						if (Double.parseDouble(line[i]) < p) {
							p = Double.parseDouble(line[i]);
						}
					} catch (NumberFormatException err) {
					}
				}

				regionsFound = new Vector<Integer>();
				for (int region : regionToChr.keySet()) {
					if (regionToChr.get(region).equals(snpChrPosition[1])) {
						boundaries = regionToBoundaries.get(region);
						if (position >= Math.max((boundaries[0] - regionSearchDistance), 0) && position <= (boundaries[0] + regionSearchDistance)) {
							regionsFound.addElement(region);
						}
					}
				}

				if (regionsFound.size() == 0) {
					regionToBoundaries.put(regionId, new int[] {position, position});
					regionToChr.put(regionId, snpChrPosition[1]);
					geneList = new Vector<String>();
					geneList.add(gene);
					regionToGenes.put(regionId, geneList);
					regionToMinP.put(regionId, p);
					regionToSnpWithMinP.put(regionId, snpChrPosition[0]);
					regionId ++;
				} else {
					regionId2 = regionsFound.elementAt(0);
					boundaries = regionToBoundaries.get(regionId2);
					boundaries[0] = Math.min(boundaries[0], position);
					boundaries[1] = Math.max(boundaries[1], position);
					if (regionToMinP.get(regionId2) > p) {
						regionToMinP.put(regionId2, p);
						regionToSnpWithMinP.put(regionId2, snpChrPosition[0]);
					}
					geneList = regionToGenes.get(regionId2);
					if (! geneList.contains(gene)) {
						geneList.add(gene);
					}

					p = regionToMinP.get(regionId2);
					regionsFound.removeElementAt(0);
					for (int region : regionsFound) {
						boundaries2 = regionToBoundaries.get(region);
						boundaries[0] = Math.min(boundaries[0], boundaries2[0]);
						boundaries[1] = Math.max(boundaries[1], boundaries2[1]);
						regionToBoundaries.remove(region);
						if (regionToMinP.get(region) < p) {
							regionToMinP.put(regionId2, regionToMinP.get(region));
							regionToSnpWithMinP.put(regionId2, regionToSnpWithMinP.get(region));
						}
						regionToMinP.remove(region);
						regionToSnpWithMinP.remove(region);
						geneList2 = regionToGenes.get(region);
						for (String gene1 : geneList2) {
							if (! geneList.contains(gene1)) {
								geneList.add(gene1);
							}
						}
						regionToGenes.remove(region);
						regionToChr.remove(region);
					}
				}
			}
		}

		for (String snpChrPos: snpsFromResultFile.keySet()) {
			line = snpsFromResultFile.get(snpChrPos);
//			if (Double.parseDouble(line[3]) > pThresholdTighter) {
			if (! screenWithCriteria(line, criteriaColumnIndices, criteriaOperators, criteriaValues, log)) {
				snpChrPosition = snpChrPos.split("\t");
//				chr = snpChrPosition[1];
				position = Integer.parseInt(snpChrPosition[2]);
				gene = line[0];
				p = 1;
				for (int i = 1; i < line.length; i++) {
					try {
						if (Double.parseDouble(line[i]) < p) {
							p = Double.parseDouble(line[i]);
						}
					} catch (NumberFormatException err) {
					}
				}

				regionsFound = new Vector<Integer>();
				for (int region : regionToChr.keySet()) {
					if (regionToChr.get(region).equals(snpChrPosition[1])) {
						boundaries = regionToBoundaries.get(region);
						if (position >= Math.max((boundaries[0] - regionSearchDistance), 0) && position <= (boundaries[0] + regionSearchDistance)) {
							regionsFound.addElement(region);
						}
					}
				}

				if (regionsFound.size() > 0) {
					regionId2 = regionsFound.elementAt(0);
					boundaries = regionToBoundaries.get(regionId2);
					boundaries[0] = Math.min(boundaries[0], position);
					boundaries[1] = Math.max(boundaries[1], position);
					if (regionToMinP.get(regionId2) > p) {
						regionToMinP.put(regionId2, p);
						regionToSnpWithMinP.put(regionId2, snpChrPos);
					}
					geneList = regionToGenes.get(regionId2);
					if (! geneList.contains(gene)) {
						geneList.add(gene);
					}

					p = regionToMinP.get(regionId2);
					regionsFound.removeElementAt(0);
					for (int region : regionsFound) {
						boundaries2 = regionToBoundaries.get(region);
						boundaries[0] = Math.min(boundaries[0], boundaries2[0]);
						boundaries[1] = Math.max(boundaries[1], boundaries2[1]);
						regionToBoundaries.remove(region);
						if (regionToMinP.get(region) < p) {
							regionToMinP.put(regionId2, regionToMinP.get(region));
							regionToSnpWithMinP.put(regionId2, regionToSnpWithMinP.get(region));
						}
						regionToMinP.remove(region);
						regionToSnpWithMinP.remove(region);
						geneList2 = regionToGenes.get(region);
						for (String gene1 : geneList2) {
							if (! geneList.contains(gene1)) {
								geneList.add(gene1);
							}
						}
						regionToGenes.remove(region);
						regionToChr.remove(region);
					}
				}
			}
		}

		if (regionToChr.size() > 0) {
			ttt = new Hashtable<String, Integer>();
			for (int region : regionToChr.keySet()) {
				ttt.put(regionToChr.get(region) + "\t" + regionToSnpWithMinP.get(region), region);
			}
			line = ttt.keySet().toArray(new String[0]);
			Sort.putInOrder(line);
	
			output = new Vector<String>();
			output.add("SNP\tSKATgene\tChr\tPos1\tPos2");
			for (String key : line) {
				regionId = ttt.get(key);
				geneList = regionToGenes.get(regionId);
				snp = regionToSnpWithMinP.get(regionId);
				chr = regionToChr.get(regionId);
				boundaries = regionToBoundaries.get(regionId);
				for (String gene1 : geneList) {
					output.add(snp + "\t" + gene1 + "\t" + chr + "\t" + boundaries[0] + "\t" + boundaries[1]);
				}
			}

			return output.toArray(new String[0]);
		} else {
			return null;
		}
	}

	public static int getKeyOfValue (Hashtable <Integer, Vector<String>> hash, String value, Logger log) {
		int keyFound;
		Vector<String> vector;

		keyFound = -1;
		for (int key : hash.keySet()) {
			vector = hash.get(key);
			if (vector.contains(value)) {
				keyFound = key;
				break;
			}
		}

		return keyFound;
	}

	/**
	 * This is a program specifically for 
	 * @param args
	 */
	public static void main(String[] args) {
		int numArgs = args.length;
		boolean isRScripts, isRScriptsSubFolders, isQcScript, isSummary, summarize, isConditional;
		String cohortRDataFilesDir = null, snpInfoDirFilenameTemplate, genoDirFilenameTemplate, phenoDirFilenameTemplate, condFileDirFilenameTemplate, rScriptDir, condResultsDir, originalResultsDir, summaryDir = null, fullPathToSnpInfo = null, phenos, races, startConditionIdNumber, rCommandLine;
		String commandConditional, commandRscript, commandRscriptSubDir, commandSummary, commandQcScript, commandPhenoDirFilenameTemplate, commandPhenos, commandGenoDirFilenameTemplate, commandRaces, commandSummaryDir, commandRCommand, commandPThreshold, commandPThresholdHigher;
		String[] commands;
		int[] columnIndeciesOfPhenoConditionEthnicAnalysis = null;
		double pThreshold, pThresholdHigher;
		int regionSearchDistance;
		Logger log;

		isRScripts = false;
		isRScriptsSubFolders = false;
		isQcScript = false;
		isSummary = false;
		isConditional = false;
		summarize = false;

		pThreshold = 0.0001;
		pThresholdHigher = 0.001;
		regionSearchDistance = 1000000;

////		cohortRDataFilesDir = "N:/statgen/CHARGE-S_conditionals/cohortRDataFiles/";
////		snpInfoFile = "D:/CHARGE-S_conditionals/snpInfos/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.RData";
////		condFileDir = "N:/statgen/CHARGE-S_conditionals/conditions/";
//		snpInfoFile = "/home/pankrat2/shared/skatMeta/snpInfos/exome_chip_v5_snpInfo_chr" + FILENAME_CHROMOSOME_SEGMENT + ".RData";
//		snpInfoFile = "/home/pankrat2/shared/skatMeta/snpInfos/snpInfoMinSubSet_CHARGES_ESP_RS_092413_chr" + FILENAME_CHROMOSOME_SEGMENT + ".RData";
		snpInfoDirFilenameTemplate = "/home/pankrat2/shared/skatMeta/freeze4/snpInfos/snpInfo_chr" + FILENAME_CHROMOSOME_SEGMENT + ".RData";
//		genoFile = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/exome_chip/genotypes_EA/EA_ARIC_noJHS_chr" + chr + "t.csv";
//		phenoFile = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/exome_chip/inflammation/ARIC_EA_LpPLA2_" + subpheno + ".csv";
		genoDirFilenameTemplate = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/freeze4_genotypes/" + FILENAME_ETHNIC_SEGMENT + "_all/" + FILENAME_ETHNIC_SEGMENT + "_ARIC_ExFrz41_all_" + FILENAME_CHROMOSOME_SEGMENT + ".RData";
		phenoDirFilenameTemplate = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/results_hemostasis/ARIC_" + FILENAME_ETHNIC_SEGMENT + "_" + FILENAME_PHENO_SEGMENT + "_ABO.csv";
//		phenoDirFilenameTemplate = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/results_hemostasis/ARIC_" + FILENAME_ETHNIC_SEGMENT + "_" + FILENAME_PHENO_SEGMENT + ".csv";
		condFileDirFilenameTemplate = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/conditional/conditions/" + FILENAME_PHENO_SEGMENT + "_" + FILENAME_CONDITION_SEGMENT + ".txt";
//		condFileDirFilenameTemplate = "D:/charges_conditional_new/conditions/" + FILENAME_PHENO_SEGMENT + "_" + FILENAME_CONDITION_SEGMENT + ".txt";

		phenos = "F8,VWF";
//		phenos = "F7,Fibrinogen";
		races = "AA,EA";
		startConditionIdNumber = "1";

//		rScriptDir = "N:/statgen/CHARGE-S_conditionals/scripts/selectedSnpInfo_MoreCohorts/";
//		resultsDir = "N:/statgen/CHARGE-S_conditionals/results/newFromSmallerSNPInfo/";
//		summaryDir = "N:/statgen/CHARGE-S_conditionals/results/summary/automated_summaries/";
		rScriptDir = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/conditional/scripts_dec2014/";
		condResultsDir = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/conditional/results_dec2014/";
		rCommandLine = "/soft/R/3.1.1/bin/Rscript";
		summaryDir = "/home/pankrat2/shared/skatMeta/aric_primary_analyses/charges/conditional/summaries/";
		columnIndeciesOfPhenoConditionEthnicAnalysis = new int[] {0, 1, 2, 3};
//		rScriptDir = "D:/charges_conditional_new/";
//		resultsDirFilenameTemplate = "D:/charges_conditional_new/outputs/";
		
//		resultsDirFilenameTemplate = "D:/charges_conditional_new/test1/";
//		resultsDirFilenameTemplate = "D:/charges_conditional_new/preConditionalResults/" + FILENAME_PHENO_SEGMENT + "_" + FILENAME_CONDITION_SEGMENT + "_" + FILENAME_ANALYSIS_SEGMENT + ".csv";
//		summaryDir = "D:/charges_conditional_new/summaries/";

//		cohortRDataFilesDir = null;
//		snpInfoFile = null;
//		condFileDir = null;
//		rScriptDir = null;
//		resultsDir = "D:/Inflammation/outputs/";
//		columnIndeciesOfPhenoConditionEthnicAnalysis = new int[] {3, 6, 1, 7};	//Indices of the following info appeared in a file name Pheno, ConditionID, Ethnic, and AnalysisName, where a file name reads like this: ARIC_AA_LpPLA2_activity_prepCondScores_seqMeta_cond1_SingleSNP.csv
//		pThreshold = 0.0001;
//		fullPathToSnpInfo = "N:/statgen/inflammation/summary/SNPInfo_ExomeChipV5.csv";
//		summaryDir = "D:/Inflammation/summary/";

//		isRScriptsSubFolders = true;
//		summarize = true;
//		cohortRDataFilesDir = "D:/CHARGE-S_conditionals/npReplicate/cohortRDataFiles/";
//		snpInfoDirFilenameTemplate = "D:/CHARGE-S_conditionals/npReplicate/snpInfos/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.Rdata";
//		condFileDirFilenameTemplate = "D:/CHARGE-S_conditionals/npReplicate/conditions/" + FILENAME_PHENO_SEGMENT + "_" + FILENAME_CONDITION_SEGMENT + ".txt";
//		rScriptDir = "D:/CHARGE-S_conditionals/npReplicate/scripts/";
//		condResultsDir = "D:/CHARGE-S_conditionals/npReplicate/results/";
		originalResultsDir = "D:/CHARGE-S_conditionals/npReplicate/originalResults/";
//		fullPathToSnpInfo = "D:/CHARGE-S_conditionals/npReplicate/snpInfos/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.Rdata";
//		summaryDir = "D:/CHARGE-S_conditionals/npReplicate/summaries/";
		
//		phenos = "Fibrinogen,F7,F8,VWF";
//		races = "EAAA,EA,AA";

		commands = new String[] {"", "rdatadir=", "snpinfofile=", "conditionfile=", "resultdir=", "scriptdir="};
		commandRscript = "-rscript";
		commandRscriptSubDir = "-rscriptsubdir";
		commandQcScript = "-qcscript";
		commandSummary = "-summary";
		commandConditional = "-conditional";
		commandPhenoDirFilenameTemplate = "phenofile=";
		commandPhenos = "phenos=";
		commandGenoDirFilenameTemplate = "genofile=";
		commandRaces = "races=";
		commandSummaryDir = "summarydir=";
		commandRCommand = "rcommand=";
		commandPThreshold = "pthreshold=";
		commandPThresholdHigher = "pthresholdhigh=";

		String usage = "\nTo generate Skat Meta R scripts for all the .RData files in a single directory:"
					+ "\n   (1) command for generating R scripts (i.e. " + commands[0] + " (default))"
					+ "\n   (2) directory of the .RData files (i.e. " + commands[1] + cohortRDataFilesDir + " (default))"
					+ "\n   (3) full path of the SNPInfo file (i.e. " + commands[2] + snpInfoDirFilenameTemplate + " (default))"
					+ "\n   (4) directory of the condition files (i.e. " + commands[3] + condFileDirFilenameTemplate + " (default))"
					+ "\n   (5) directory of results the condition files (i.e. " + commands[4] + condResultsDir + " (default))"
					+ "\n   (6) directory to output the R scripts (i.e. " + commands[5] + rScriptDir + " (default))"
					+ "\n"
					+ "\nTo generate Skat Meta R scripts for all the .RData files in a single directory and all its subdirectories:"
					+ "\n   (1) command for generating R scripts for a directory and all its subdirectories (i.e. " + commands[6] + " (default))"
					+ "\n   (2) directory of the .RData files (i.e. " + commands[1] + cohortRDataFilesDir + " (default))"
					+ "\n   (3) full path of the SNPInfo file (i.e. " + commands[2] + snpInfoDirFilenameTemplate + " (default))"
					+ "\n   (4) directory of the condition files (i.e. " + commands[3] + condFileDirFilenameTemplate + " (default))"
					+ "\n   (5) directory of results the condition files (i.e. " + commands[4] + condResultsDir + " (default))"
					+ "\n   (6) directory to output the R scripts (i.e. " + commands[5] + rScriptDir + " (default))"
					+ "\n"
					+ "\nTo generate QC R scripts for all the .RData files in a single directory:"
					+ "\n   (1) command for generating R scripts for a directory and all its subdirectories (i.e. " + commands[7] + " (default))"
					+ "\n   (2) directory of the .RData files (i.e. " + commands[1] + cohortRDataFilesDir + " (default))"
					+ "\n"
					+ "\nTo summarize Skat Meta results:"
					+ "\n   (1) command for summarizing Skat Meta results (i.e. " + commands[8] + " (default))"
					+ "\n   (2) directory of Skat Meta results (i.e. " + commands[4] + condResultsDir + " (default))"
					+ "\n   (3) directory of original meta-analysis results (i.e. originalResultsDir=" + originalResultsDir + " (default))"
					+ "\n   (4) directory to output summary (i.e. " + commandSummaryDir + summaryDir + " (default))"
					+ "\n"
					+ "\nTo run conditional analysis from results of pre-conditional analysis:"
					+ "\n   (1) command for conditional analysis (i.e. " + commandConditional + " (default))"
					+ "\n   (2) template of dir and file name of the results of pre-conditional analyses, which the 1st condition is going to based on (i.e. " + commands[4] + condResultsDir + " (default))"
					+ "\n   (3) phenos (i.e. " + commandPhenos + phenos + " (default))"
					+ "\n   (4) template of pheno files' dir and name (i.e. " + commandPhenoDirFilenameTemplate + phenoDirFilenameTemplate + " (default))"
					+ "\n   (5) template of geno files' dir and name (i.e. " + commandGenoDirFilenameTemplate + genoDirFilenameTemplate + " (default))"
					+ "\n   (6) template of SNPInfo files' dir and name (i.e. " + commands[2] + snpInfoDirFilenameTemplate + " (default))"
					+ "\n   (7) races (i.e. " + commandRaces + races + " (default))"
					+ "\n   (8) folder to output the R scripts (i.e. " + commands[5] + rScriptDir + " (default))"
					+ "\n   (9) template of folder and name of  (i.e. " + commands[3] + condFileDirFilenameTemplate + " (default))"
					+ "\n   (10) folder to output summaries of Skat Meta results (i.e. " + commandSummaryDir + summaryDir + " (default))"
					+ "\n   (11) command to launch R (i.e. " + commandRCommand + rCommandLine + " (default))"
					+ "\n   (12) p-value threshold for next conditions (i.e. " + commandPThreshold + pThreshold  + " (default))"
					+ "\n"
					+ "\nTo run conditional analysis with existing 1st condition:"
					+ "\n   (1) command for conditional analysis (i.e. " + commandConditional + " (default))"
					+ "\n   (2) template of condition files' directory and name (i.e. " + commands[3] + condFileDirFilenameTemplate + " (default))"
					+ "\n   (3) phenos (i.e. " + commandPhenos + phenos + " (default))"
					+ "\n   (4) template of pheno files' dir and name (i.e. " + commandPhenoDirFilenameTemplate + phenoDirFilenameTemplate + " (default))"
					+ "\n   (5) template of geno files' dir and name (i.e. " + commandGenoDirFilenameTemplate + genoDirFilenameTemplate + " (default))"
					+ "\n   (6) template of SNPInfo files' dir and name (i.e. " + commands[2] + snpInfoDirFilenameTemplate + " (default))"
					+ "\n   (7) races (i.e. " + commandRaces + races + " (default))"
					+ "\n   (8) folder to output the R scripts (i.e. " + commands[5] + rScriptDir + " (default))"
					+ "\n   (9) folder to output Skat Meta results (i.e. " + commands[4] + condResultsDir + " (default))"
					+ "\n   (10) folder to output summaries of Skat Meta results (i.e. " + commandSummaryDir + summaryDir + " (default))"
					+ "\n   (11) command to launch R (i.e. " + commandRCommand + rCommandLine + " (default))"
					+ "\n   (12) p-value threshold for next conditions (i.e. " + commandPThreshold + pThreshold  + " (default))"
					+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(commandRscript)) {
				isRScripts = true;
				numArgs--;
			} else if (args[i].startsWith(commands[1])) {
				cohortRDataFilesDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[2])) {
				snpInfoDirFilenameTemplate = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[3])) {
				condFileDirFilenameTemplate = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[4])) {
				condResultsDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("originalResultsDir=")) {
				originalResultsDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[5])) {
				rScriptDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandRscriptSubDir)) {
				isRScriptsSubFolders = true;
				numArgs--;
			} else if (args[i].startsWith(commandQcScript)) {
				isQcScript = true;
				numArgs--;
			} else if (args[i].startsWith(commandSummary)) {
				isSummary = true;
				numArgs--;
			} else if (args[i].startsWith(commandSummaryDir)) {
				summaryDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandConditional)) {
				isConditional = true;
				numArgs--;
			} else if (args[i].startsWith(commandPhenoDirFilenameTemplate)) {
				phenoDirFilenameTemplate = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandGenoDirFilenameTemplate)) {
				genoDirFilenameTemplate = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandPhenos)) {
				phenos = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandRaces)) {
				races = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandRCommand)) {
				rCommandLine = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandPThreshold)) {
				pThreshold = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commandPThresholdHigher)) {
				pThresholdHigher = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-summarize")) {
				summarize = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		if (isRScripts) {
			log = new Logger();
//			generateSkatMetaRScript(phenoFile, genoFile, condFileDir, snpInfoFile, new String[] {"F7", "F8_ABO", "VWF_ABO", "Fibrinogen"}, new String[] {"AA", "EA"}, new String[] {"1"}, rScriptDir, resultsDir, rCommandLine);
			generateSkatMetaRScript(phenoDirFilenameTemplate, genoDirFilenameTemplate, condFileDirFilenameTemplate, snpInfoDirFilenameTemplate, phenos.split(","), races.split(","), startConditionIdNumber, rScriptDir, condResultsDir, rCommandLine);
//			generateSkatMetaRScriptConditional(cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
		} else if (isRScriptsSubFolders) {
			log = new Logger();
			generateSkatMetaRScriptsForSubFolders(cohortRDataFilesDir, snpInfoDirFilenameTemplate, condFileDirFilenameTemplate, rScriptDir, condResultsDir);
		} else if (isQcScript) {
			log = new Logger();
			qcScript(cohortRDataFilesDir);
		} else if (isSummary) {
			log = new Logger(condResultsDir + "SuperNovo_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
			summary(condResultsDir, columnIndeciesOfPhenoConditionEthnicAnalysis, pThreshold, null, fullPathToSnpInfo, null, summaryDir, log);
		} else if (summarize) {
			log = new Logger(summaryDir + "Summarize_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
			summaryTable(phenos.split(","), races.split(","), condResultsDir, originalResultsDir, summaryDir, snpInfoDirFilenameTemplate, log);
		} else if (isConditional) {
			log = new Logger();
//			conditionalAnalysisWholeProcess(fullpathToPreviousRoundResult, fullpathToPreviousRoundCondition, fullpathToOutputNextRoundCondition, pThreshold, log);
//			Files.writeList(getFirstCondition("D:/charges_conditional_new/preConditionalResults/vWF_SingleVariant.csv", null, null, .0000001, .00001, 1000000, null), "D:/charges_conditional_new/preConditionalResults/vWF_cond1.txt");
//			conditionalAnalysisWholeProcessMultiplePhenos(condFileDirFilenameTemplate, phenoDirFilenameTemplate, genoDirFilenameTemplate, snpInfoDirFilenameTemplate, phenos.split(","), races.split(","), Integer.parseInt(startConditionIdNumber), rScriptDir, resultsDir, summaryDir, rCommandLine, pThreshold, log);
			conditionalAnalysisWholeProcessMultiplePhenos(condResultsDir, condFileDirFilenameTemplate, phenoDirFilenameTemplate, genoDirFilenameTemplate, snpInfoDirFilenameTemplate, phenos.split(","), races.split(","), rScriptDir, summaryDir, rCommandLine, pThreshold, pThresholdHigher, regionSearchDistance, log);
		} else {
			log = new Logger();
			log.reportError("No command executed, due to none of the following is specified: " + commandRscript + ", " + commandRscriptSubDir + ", " + commandQcScript + ", " + commandSummary + ", or " + commandConditional + ".\n" + usage);
		}
		
		log.report("Program completed.");
	}
}
