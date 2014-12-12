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
import common.Sort;
import common.ext;

public class SkatMeta {
	public static final String[] GENE_RESULT_COLUMNS = new String[] {"gene", "p"};
	public static final String[] SNPINFO_COLUMNS = new String[] {"Name", "Chr", "MapInfo"};

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

	public static Hashtable<String, String[]> loadFile(String filefullpath, Hashtable<String, String[]> snpList, String[] columnsAsTheKey, String[] otherColumnsNeeded, String[] criteriaColumns, String[] criteria) {
		BufferedReader reader;
		Hashtable<String, String[]> result;
		String key, delimiter, header;
		String[] line, tmp;
		int[] indicesKey, indicesOther, indicesCriteria = null, indicesCriteriaGroup_LessThanOrEqualTo = null, indicesCriteriaGroup_EqualsString = null;
		Vector <Integer> criteriaGroup_LessThanOrEqualTo_tmp, criteriaGroup_EqualsString_tmp;
		String[] criteriaColumns_LessThanOrEqualTo, criteriaColumns_EqualsString;
		double[] criteriaGroup_LessThanOrEqualTo = null;
		String[] criteriaGroup_EqualsString = null;
		boolean isToInclude;

		result = new Hashtable<String, String[]>();
		try {
			reader = new BufferedReader(new FileReader(filefullpath));
			header = reader.readLine();
			delimiter = ext.determineDelimiter(header);
			line = header.replaceAll("\"", "").split(delimiter);
			indicesKey = ext.indexFactors(columnsAsTheKey, line, false, true);
			indicesOther = ext.indexFactors(otherColumnsNeeded, line, false, true);
			if (criteriaColumns != null) {
				criteriaGroup_LessThanOrEqualTo_tmp = new Vector<Integer>();
				criteriaGroup_EqualsString_tmp = new Vector<Integer>();
				for (int i = 0; i < criteriaColumns.length; i++) {
					if (criteria[i].startsWith("<=")) {
						criteriaGroup_LessThanOrEqualTo_tmp.add(i);
					} else if ((criteria[i].charAt(0) >= 65 && criteria[i].charAt(0) <= 90) || (criteria[i].charAt(0) >= 97 && criteria[i].charAt(0) <= 122)) {
						criteriaGroup_EqualsString_tmp.add(i);
					} else {
						System.out.println("Error - current version does not support the criteria " + criteriaColumns[i] + " " + criteria[i]);
						System.exit(0);
					}
				}
				if (criteriaGroup_LessThanOrEqualTo_tmp.size() > 0) {
					criteriaColumns_LessThanOrEqualTo = new String[criteriaGroup_LessThanOrEqualTo_tmp.size()];
					criteriaGroup_LessThanOrEqualTo = new double[criteriaGroup_LessThanOrEqualTo_tmp.size()];
					for (int i = 0; i < criteriaColumns_LessThanOrEqualTo.length; i++) {
						criteriaColumns_LessThanOrEqualTo[i] = criteriaColumns[criteriaGroup_LessThanOrEqualTo_tmp.elementAt(i)];
						criteriaGroup_LessThanOrEqualTo[i] = Double.parseDouble(criteria[criteriaGroup_LessThanOrEqualTo_tmp.elementAt(i)].substring(2));
					}
					indicesCriteriaGroup_LessThanOrEqualTo = ext.indexFactors(criteriaColumns_LessThanOrEqualTo, line, false, true);
				}
				if (criteriaGroup_EqualsString_tmp.size() > 0) {
					criteriaColumns_EqualsString = new String[criteriaGroup_EqualsString_tmp.size()];
					criteriaGroup_EqualsString = new String[criteriaGroup_EqualsString_tmp.size()];
					for (int i = 0; i < criteriaColumns_EqualsString.length; i++) {
						criteriaColumns_EqualsString[i] = criteriaColumns[criteriaGroup_EqualsString_tmp.elementAt(i)];
						criteriaGroup_EqualsString[i] = criteria[criteriaGroup_EqualsString_tmp.elementAt(i)];
					}
					indicesCriteriaGroup_EqualsString = ext.indexFactors(criteriaColumns_EqualsString, line, false, true);
				}
//				indicesCriteria = ext.indexFactors(criteriaColumns, line, false, true);
			}
			while (reader.ready()) {
				line = reader.readLine().replaceAll("\"", "").split(delimiter);
				isToInclude = true;
//				for (int i = 0; indicesCriteria != null && i < indicesCriteria.length; i++) {
//					if (! line[indicesCriteria[i]].equalsIgnoreCase(criteria[i])) {
//						isToInclude = false;
//						break;
//					}
//				}
				for (int i = 0; indicesCriteriaGroup_LessThanOrEqualTo != null && i < indicesCriteriaGroup_LessThanOrEqualTo.length; i++) {
					if (line[1].equals("exm-rs4765623")) {
						System.out.println("");
					}
					if (line[indicesCriteriaGroup_LessThanOrEqualTo[i]].equals("NA") || Double.parseDouble(line[indicesCriteriaGroup_LessThanOrEqualTo[i]]) > criteriaGroup_LessThanOrEqualTo[i]) {
						isToInclude = false;
						break;
					}
				}
				for (int i = 0; indicesCriteriaGroup_EqualsString != null && i < indicesCriteriaGroup_EqualsString.length; i++) {
					if (! line[indicesCriteriaGroup_EqualsString[i]].equalsIgnoreCase(criteriaGroup_EqualsString[i])) {
						isToInclude = false;
						break;
					}
				}
				if (isToInclude) {
					key = line[indicesKey[0]];
					for (int i = 1; i < indicesKey.length; i++) {
						key += ("\t" + line[indicesKey[i]]);
					}
					if (snpList != null) {
						tmp = snpList.get(key);
						for (int i = 0; i < tmp.length; i++) {
							key += ("\t" + tmp[i]);
						}
					}
					tmp = new String[indicesOther.length];
					for (int i = 0; i < indicesOther.length; i++) {
						tmp[i] = line[indicesOther[i]];
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
					analysesGroup = ethnicGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]]);
					log.reportError("Error - " + filenames[i] + " get duplicated with " + analysesGroup.get(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]]) + ". \nSystem halted.");
				} else {
					analysesGroup.put(filenameRoot[columnIndeciesOfPhenoConditionEthnicAnalysis[3]], filenames[i]);
				}
			}
		}

		return phenoGroup;
	}

	public static void conditionalAnalysisWholeProcess (String fullpathToPreviousRoundResult, String fullpathToPreviousRoundCondition, String fullpathToOutputNextRoundCondition, double pThreshold, Logger log) {
		int a = 1;

		while (developConditions("D:/Inflammation/outputs/ARIC_AAEA_LpPLA2_activity_prepCondScores_seqMeta_cond" + a + "_SingleSNP.csv", "D:/Inflammation/tests/activity_cond" + a + ".txt", "D:/Inflammation/tests/activity_cond" + (a + 1) + ".txt", pThreshold, log)) {
			//TODO run R script;
			a ++;
		}
	}

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
		boolean isRScripts, isRScriptsSubFolders, isQcScript, isSummary, isConditional;
		String cohortRDataFilesDir, snpInfoFile, condFileDir, rScriptDir, resultsDir, summaryDir, fullPathToSnpInfo;
		String command, commandConditional;
		String[] commands;
		int[] columnIndeciesOfPhenoConditionEthnicAnalysis;
		double pThreshold;
		Logger log;

		isRScripts = false;
		isRScriptsSubFolders = false;
		isQcScript = false;
		isSummary = false;
		isConditional = false;
		cohortRDataFilesDir = "N:/statgen/CHARGE-S_conditionals/cohortRDataFiles/";
		snpInfoFile = "D:/CHARGE-S_conditionals/snpInfos/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.RData";
		condFileDir = "N:/statgen/CHARGE-S_conditionals/conditions/";
		rScriptDir = "N:/statgen/CHARGE-S_conditionals/scripts/selectedSnpInfo_MoreCohorts/";
		resultsDir = "N:/statgen/CHARGE-S_conditionals/results/newFromSmallerSNPInfo/";
		summaryDir = "N:/statgen/CHARGE-S_conditionals/results/summary/automated_summaries/";

		commands = new String[] {"-rscript", "rdatadir=", "snpinfo=", "conditionsdir=", "resultdir=", "scriptdir=", "-rscriptsubdir", "-qcscript", "-summary", "summarydir="};
		commandConditional = "-conditional";
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

		cohortRDataFilesDir = null;
		snpInfoFile = null;
		condFileDir = null;
		rScriptDir = null;
		resultsDir = "D:/Inflammation/outputs/";
		columnIndeciesOfPhenoConditionEthnicAnalysis = new int[] {3, 6, 1, 7};	//Indices of the following info appeared in a file name Pheno, ConditionID, Ethnic, and AnalysisName, where a file name reads like this: ARIC_AA_LpPLA2_activity_prepCondScores_seqMeta_cond1_SingleSNP.csv
		pThreshold = 0.0001;
		fullPathToSnpInfo = "N:/statgen/inflammation/summary/SNPInfo_ExomeChipV5.csv";
		summaryDir = "D:/Inflammation/summary/";

//		isRScriptsSubFolders = true;
//		isSummary = true;
		isConditional = false;

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
			summary(resultsDir, columnIndeciesOfPhenoConditionEthnicAnalysis, pThreshold, null, fullPathToSnpInfo, summaryDir, log);
		} else if (isConditional) {
			log = new Logger();
//			conditionalAnalysisWholeProcess(fullpathToPreviousRoundResult, fullpathToPreviousRoundCondition, fullpathToOutputNextRoundCondition, pThreshold, log);
			conditionalAnalysisWholeProcess(null, null, null, pThreshold, log);
		} else {
			log = new Logger();
			log.reportError("No command executed, due to none of the following is specified: " + commands[0] + ", " + commands[12] + ", " + commands[19] + ", or " + commands[20] + ".");
		}
		
		log.report("Program completed.");
	}

}
