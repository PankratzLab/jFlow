package one;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

import common.Files;
import common.ext;

public class SkatMeta {
	public static void generateRScript(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
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

	public static void subFolders(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
		String[] folders;
		
		folders = Files.listDirectories(sourceRDataFilesDir, false);
		if (folders == null || folders.length < 1) {
				generateRScript(sourceRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
		} else {
			for (int i = 0; i < folders.length; i++) {
				subFolders(sourceRDataFilesDir + folders[i] + "/", snpInfoFile, condFileDir, rScriptDir, resultsDir);
			}
		}
	}

	public static void qcScript(String sourceRDataFilesDir) {
		String[] files;
		PrintWriter writer;

		files = listFilesInDirAndAllSubDirs(sourceRDataFilesDir);
		if (files != null || files.length > 0) {
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

	public static void main(String[] args) {
		String cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir;

//		sourceRDataFilesDir = "N:/statgen/CHARGE-S_conditionals/F7/cond1/";
//		rScriptDir = "N:/statgen/CHARGE-S_conditionals/scripts/";
//		resultsDir = "N:/statgen/CHARGE-S_conditionals/results/";
		cohortRDataFilesDir = "N:/statgen/CHARGE-S_conditionals/cohortRDataFiles/";
		snpInfoFile = "D:/CHARGE-S_conditionals/snpInfos/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.RData";
		condFileDir = "N:/statgen/CHARGE-S_conditionals/conditions/";
		rScriptDir = "N:/statgen/CHARGE-S_conditionals/scripts/selectedSnpInfo/";
		resultsDir = "N:/statgen/CHARGE-S_conditionals/results/newFromSmallerSNPInfo/";

//		generateRScript(cohortRDataFilesDir, snpInfoFile, rScriptDir, resultsDir);
//		subFolders(cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
		qcScript(cohortRDataFilesDir);
	}

}
