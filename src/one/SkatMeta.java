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

import common.Files;
import common.Logger;
import common.ext;

public class SkatMeta {
	public static final String[] GENE_RESULT_COLUMNS = new String[] {"gene", "p"};

	public static void generateSkatMetaRScript(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
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

	public static void generateSkatMetaRScriptsForSubFolders(String sourceRDataFilesDir, String snpInfoFile, String condFileDir, String rScriptDir, String resultsDir) {
		String[] folders;
		
		folders = Files.listDirectories(sourceRDataFilesDir, false);
		if (folders == null || folders.length < 1) {
				generateSkatMetaRScript(sourceRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
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

	public static void summary(String resultsDir, String summaryDir, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvalSummary;
		String[] phenoList, analysesList, columnsAsTheKey, otherColumnsNeeded;
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll;
		String[] geneSnpList;
		Hashtable<String, Double> snpsWithSignificantPval;

		if (log == null) {
			log = new Logger();
		}

		phenoGroups = groupFileNames(Files.list(resultsDir, null, ".csv", false, false), log);
		if (phenoGroups != null && phenoGroups.size() > 0) {
			analysesList = new String[] {"T5Count", "SKAT_T5"};
			genePvalSummary = summarizeGenePvalues(phenoGroups, analysesList, resultsDir);
			phenoList = getPhenoList(genePvalSummary, log);
//			getListsOfPhenosConditionsEthnicsAnalysesGenes(genePvalSummary, phenoList, conditionListAllPhenos, ethnicList, analysesList, geneListAllPhenos, log);
			printSummarizedGenePvalues(genePvalSummary, phenoList, new String[] {"EA", "AA", "EAAA"}, analysesList, summaryDir + "summary_genePvalues.txt", log);

			analysesList = new String[] {"SingleSNP"};
			columnsAsTheKey = new String[] {"gene", "Name"};
			otherColumnsNeeded = new String[] {"maf", "ntotal", "nmiss", "beta", "se", "p"};
			for (String pheno : phenoList) {
				snpResultsAll = summarizeSnpPvalues(phenoGroups.get(pheno), resultsDir, analysesList, columnsAsTheKey, otherColumnsNeeded, log);
				geneSnpList = getListOfGeneSnps(snpResultsAll, log);
				snpsWithSignificantPval = printSnpResults(snpResultsAll, getConditionList(genePvalSummary, pheno, log), new String[] {"EAAA", "EA", "AA"}, analysesList, geneSnpList, otherColumnsNeeded, new String[] {"maf", "ntotal", "nmiss", "p"}, 0.00000185, summaryDir + "summary_" + pheno + "_snps.txt", log);
				printSummarizedSnpResults(snpResultsAll, snpsWithSignificantPval, getConditionList(genePvalSummary, pheno, log), new String[] {"EAAA", "EA", "AA"}, analysesList, otherColumnsNeeded, summaryDir + "summary_" + pheno + "_snpPvalues.txt", log);
			}
		}
	}

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> summarizeGenePvalues(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroups, String[] analysesTypesToSelect, String resultsDir) {
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
									indices = ext.indexFactors(GENE_RESULT_COLUMNS, reader.readLine().replaceAll("\"", "").split(","), false, true);
									while (reader.ready()) {
										line = reader.readLine().replaceAll("\"", "").split(",");
										summaryGeneGroup.put(line[indices[0]], line[indices[1]]);
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

	public static void printSummarizedGenePvalues(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>> genePvaluesSummary, String[] phenoList, String[] ethnicList, String[] analysesList, String fullPathOutFilename, Logger log) {
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
				for (String pheno : phenoList) {
					if (ethnicList == null) {
						ethnicList = getEthnicsList(genePvaluesSummary, pheno, log);
					}
					if (analysesList == null) {
						analysesList = getAnalysesList(genePvaluesSummary, pheno, log);
					}
					writer.println(pheno);
					line = "ethnic\tgene";
					conditionList = getConditionList(genePvaluesSummary, pheno, log);
					for (String condition : conditionList) {
						for (String analysis : analysesList) {
							line += ("\t" + condition + "_" + analysis);
						}
					}
					writer.println(line);

					geneList = getGeneList(genePvaluesSummary, pheno, log);
					conditionGroup = genePvaluesSummary.get(pheno);
					for (String ethnic : ethnicList) {
						for (String gene : geneList) {
							line = (ethnic + "\t" + gene);
							for (String condition : conditionList) {
								if (conditionGroup.containsKey(condition)) {
									ethnicGroup = conditionGroup.get(condition);
									if (ethnicGroup.containsKey(ethnic)) {
										analysesGroup = ethnicGroup.get(ethnic);
										for (String analysis : analysesList) {
											geneGroup = analysesGroup.get(analysis);
											line += ("\t" + (geneGroup.containsKey(gene)? geneGroup.get(gene) : ""));
										}
									} else {
										for (int i = 0; i < analysesList.length; i++) {
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

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> summarizeSnpPvalues(Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroups, String dir, String[] analysesNeeded, String[] columnsAsTheKey, String[] otherColumnsNeeded, Logger log) {
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
						snpResultsEachFile = loadFile(dir + analysesGroup.get(analysis), columnsAsTheKey, otherColumnsNeeded, null, null);
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

	public static Hashtable<String, String[]> loadFile(String filefullpath, String[] columnsAsTheKey, String[] otherColumnsNeeded, String[] criteriaColumns, double[] criteria) {
		BufferedReader reader;
		Hashtable<String, String[]> result;
		String key;
		String[] line, tmp;
		int[] indicesKey, indicesOther;

		result = new Hashtable<String, String[]>();
		try {
			reader = new BufferedReader(new FileReader(filefullpath));
			line = reader.readLine().replaceAll("\"", "").split(",");
			indicesKey = ext.indexFactors(columnsAsTheKey, line, false, true);
			indicesOther = ext.indexFactors(otherColumnsNeeded, line, false, true);
			while (reader.ready()) {
				line = reader.readLine().replaceAll("\"", "").split(",");
				key = line[indicesKey[0]];
				for (int i = 1; i < indicesKey.length; i++) {
					key += ("\t" + line[indicesKey[i]]);
				}
				tmp = new String[indicesOther.length];
				for (int i = 0; i < indicesOther.length; i++) {
					tmp[i] = line[indicesOther[i]];
				}
				result.put(key, tmp);
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return result;
	}

	public static String[] getListOfGeneSnps(Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> snpResultsCurrentCondition;
		Hashtable<String, Hashtable<String, String[]>> snpResultsCurrentEthnic;
		Hashtable<String, String[]> snpResultsCurrentAnalysis;
		HashSet<String> geneSnpList = null;
		String[] result;

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
		}

		result = geneSnpList.toArray(new String[0]);
		Arrays.sort(result);
		return result;
	}

	public static Hashtable<String, Double> printSnpResults (Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>> snpResultsAll, String[] conditionList, String[] ethnicList, String[] analysesList, String[] geneSnpList, String[] columnNamesOfTheData, String[] columnNamesToOutput, Double threshold, String fullPathOutFilename, Logger log) {
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
				line = "gene\tsnp\tmin_p";
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
				line = "gene\tsnp\tmin_p";
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

	public static Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>  groupFileNames(String[] filenames, Logger log) {
		Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>> phenoGroup = null;
		Hashtable<String, Hashtable<String, Hashtable<String, String>>> conditionGroup;
		Hashtable<String, Hashtable<String, String>> ethnicGroup;
		Hashtable<String, String> analysesGroup;
		String[] tmp;

		if (log == null) {
			log = new Logger();
		}

		if (filenames != null && filenames.length > 0) {
			phenoGroup = new Hashtable<String, Hashtable<String, Hashtable<String, Hashtable<String, String>>>>();
			for (int i = 0; i < filenames.length; i++) {
				tmp = ext.rootOf(filenames[i]).split("_");
				if (phenoGroup.containsKey(tmp[0])) {
					conditionGroup = phenoGroup.get(tmp[0]);
				} else {
					conditionGroup = new Hashtable<String, Hashtable<String, Hashtable<String, String>>>();
					phenoGroup.put(tmp[0], conditionGroup);
				}

				if (conditionGroup.containsKey(tmp[1])) {
					ethnicGroup = conditionGroup.get(tmp[1]);
				} else {
					ethnicGroup = new Hashtable<String, Hashtable<String, String>>();
					conditionGroup.put(tmp[1], ethnicGroup);
				}

				if (ethnicGroup.containsKey(tmp[2])) {
					analysesGroup = ethnicGroup.get(tmp[2]);
				} else {
					analysesGroup = new Hashtable<String, String>();
					ethnicGroup.put(tmp[2], analysesGroup);
				}

				if (tmp.length > 4) {
					tmp[3] = tmp[3] + "_" + tmp[4];
				}

				if (analysesGroup.containsKey(tmp[3])) {
					analysesGroup = ethnicGroup.get(tmp[3]);
					log.reportError("Error - " + filenames[i] + " get duplicated with " + analysesGroup.get(tmp[3]) + ". \nSystem halted.");
				} else {
					analysesGroup.put(tmp[3], filenames[i]);
				}
			}
		}

		return phenoGroup;
	}

	/**
	 * This is a program specifically for 
	 * @param args
	 */
	public static void main(String[] args) {
		int numArgs = args.length;
		boolean isRScripts, isRScriptsSubFolders, isQcScript, isSummary;
		String cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir, summaryDir;
		String[] commands;
		Logger log;

		isRScripts = false;
		isRScriptsSubFolders = false;
		isQcScript = false;
		isSummary = false;
//		sourceRDataFilesDir = "N:/statgen/CHARGE-S_conditionals/F7/cond1/";
//		rScriptDir = "N:/statgen/CHARGE-S_conditionals/scripts/";
//		resultsDir = "N:/statgen/CHARGE-S_conditionals/results/";
		cohortRDataFilesDir = "N:/statgen/CHARGE-S_conditionals/cohortRDataFiles/";
		snpInfoFile = "D:/CHARGE-S_conditionals/snpInfos/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.RData";
		condFileDir = "N:/statgen/CHARGE-S_conditionals/conditions/";
		rScriptDir = "N:/statgen/CHARGE-S_conditionals/scripts/selectedSnpInfo_MoreCohorts/";
		resultsDir = "N:/statgen/CHARGE-S_conditionals/results/newFromSmallerSNPInfo/";
		summaryDir = "N:/statgen/CHARGE-S_conditionals/results/summary/";

		commands = new String[] {"-rscript", "rdatadir=", "snpinfo=", "conditionsdir=", "resultdir=", "scriptdir=", "-rscriptsubdir", "-qcscript", "-summary", "summarydir="};
		String usage = "\nTo generate Skat Meta R scripts for all the .RData files in a single directory:"
					+ "\n   (1) command for generating R scripts (i.e. " + commands[0] + " (default))"
					+ "\n   (2) directory of the .RData files (i.e. " + commands[1] + cohortRDataFilesDir + " (default))"
					+ "\n   (3) full path of the SNPInfo file (i.e. " + commands[2] + snpInfoFile + " (default))"
					+ "\n   (4) directory of the condition files (i.e. " + commands[3] + condFileDir + " (default))"
					+ "\n   (5) directory of results the condition files (i.e. " + commands[4] + resultsDir + " (default))"
					+ "\n   (6) directory to output the R scripts (i.e. " + commands[5] + rScriptDir + " (default))"
					+ "\n"
					+ "\nTo generate Skat Meta R scripts for all the .RData files in a single directory and all its subdirectories:"
					+ "\n   (1) command for generating R scripts for a directory and all its subdirectories (i.e. " + commands[6] + " (default))"
					+ "\n   (2) directory of the .RData files (i.e. " + commands[1] + cohortRDataFilesDir + " (default))"
					+ "\n   (3) full path of the SNPInfo file (i.e. " + commands[2] + snpInfoFile + " (default))"
					+ "\n   (4) directory of the condition files (i.e. " + commands[3] + condFileDir + " (default))"
					+ "\n   (5) directory of results the condition files (i.e. " + commands[4] + resultsDir + " (default))"
					+ "\n   (6) directory to output the R scripts (i.e. " + commands[5] + rScriptDir + " (default))"
					+ "\n"
					+ "\nTo generate QC R scripts for all the .RData files in a single directory:"
					+ "\n   (1) command for generating R scripts for a directory and all its subdirectories (i.e. " + commands[7] + " (default))"
					+ "\n   (2) directory of the .RData files (i.e. " + commands[1] + cohortRDataFilesDir + " (default))"
					+ "\n"
					+ "\nTo summerize Skat Meta results:"
					+ "\n   (1) command for summarizing Skat Meta results (i.e. " + commands[8] + " (default))"
					+ "\n   (2) directory of Skat Meta results (i.e. " + commands[4] + resultsDir + " (default))"
					+ "\n   (3) directory to output summary (i.e. " + commands[9] + resultsDir + " (default))"
					+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(commands[0])) {
				isRScripts = true;
				numArgs--;
			} else if (args[i].startsWith(commands[1])) {
				cohortRDataFilesDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[2])) {
				snpInfoFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[3])) {
				condFileDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[4])) {
				resultsDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[5])) {
				rScriptDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[6])) {
				isRScriptsSubFolders = true;
				numArgs--;
			} else if (args[i].startsWith(commands[7])) {
				isQcScript = true;
				numArgs--;
			} else if (args[i].startsWith(commands[8])) {
				isSummary = true;
				numArgs--;
			} else if (args[i].startsWith(commands[9])) {
				summaryDir = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

//		isRScriptsSubFolders = true;
//		isSummary = true;

		if (isRScripts) {
			log = new Logger();
			generateSkatMetaRScript(cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
		} else if (isRScriptsSubFolders) {
			log = new Logger();
			generateSkatMetaRScriptsForSubFolders(cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir);
		} else if (isQcScript) {
			log = new Logger();
			qcScript(cohortRDataFilesDir);
		} else if (isSummary) {
			log = new Logger(resultsDir + "SuperNovo_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
			summary(resultsDir, summaryDir, log);
		} else {
			log = new Logger();
			log.reportError("No command executed, due to none of the following is specified: " + commands[0] + ", " + commands[12] + ", " + commands[19] + ", or " + commands[20] + ".");
		}
		
		log.report("Program completed.");
	}

}
