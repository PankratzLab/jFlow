package gwas;

import java.io.*;
import java.util.*;

import stats.ProbDist;
import stats.Rscript;

import common.*;
import filesys.*;

public class SeqMeta {
	public static final String[] ALGORITHMS = {
		"singlesnpMeta", 
		"burdenMeta", 
		"skatMeta", 
		"skatOMeta"
	};
	
	public static final String[][] UNIT_OF_ANALYSIS = {
		Aliases.MARKER_NAMES,
		Aliases.GENE_UNITS,
		Aliases.GENE_UNITS,
		Aliases.GENE_UNITS
	};

	public static final boolean[] SINGLE_VARIANTS = {
		true,
		false,
		false,
		false
	};
	
	public static final String[][] HEADER_TYPES = {
		{"gene", "Name", "p", "maf", "nmiss", "ntotal", "beta", "se"}, // Single SNP
		{"gene", "p", "beta", "se", "cmafTotal", "cmafUsed", "nsnpsTotal", "nsnpsUsed", "nmiss"}, // Burden Test
		{"gene", "p", "Qmeta", "cmaf", "nmiss", "nsnps"}, // SKAT test
		{"gene", "p", "Qmeta", "cmaf", "nmiss", "nsnps"} // SKAT-O test (not verified)
	};
	
	public static String getRscriptExecutable(MetaAnalysisParams maps, Logger log) {
		if (maps != null && maps.getRExec() != null) {
			return maps.getRExec();
		} else {
			return Rscript.getRscriptExecutable(log);
		}
	}
	
	private static void determineObjectNames(String dir, MetaAnalysisParams maps, Logger log) {
		String[] lines, files;
		String[][] iterations;
		String root, commands, filename;
		Vector<String> v, remaining;
		
		files = Files.list("./", null, ".rdata", false, false);
		log.report("There are "+files.length+" total .Rdata files");
		
		dir = new File(dir).getAbsolutePath()+"/";

		v = new Vector<String>();
		remaining = new Vector<String>();
		new File("batchChecks/").mkdir();
		for (int i = 0; i < files.length; i++) {
			root = ext.rootOf(files[i]);
			if (!Files.exists("batchChecks/"+root+".object")) {
				lines = new String[] {
						"load(\""+files[i]+"\")",
						"name <- ls()",
	//					"write.table(name, \"name.txt\", sep=\"\t\")",
						"fileConn<-file(\"batchChecks/"+root+".object\")",
						"writeLines(c(name), fileConn)",
						"close(fileConn)",
				};
				filename = "batchChecks/"+root+".R";
				Files.writeList(lines, filename);
				v.add(filename);
				remaining.add(files[i]);
			}
		}
		log.report("There are "+v.size()+" .Rdata files remaining to interrogate:\n"+Array.toStr(Array.toStringArray(v), "\n"));

		if (v.size() > 0) {
			commands = getRscriptExecutable(maps, log)+" --no-save [%0]";
			iterations = Matrix.toMatrix(Array.toStringArray(v));
			Files.qsub("batchChecks/checkObject", dir, -1, commands, iterations, 4000, 1);
			Files.batchIt("master.checkObjectAll", null, 1, commands, iterations);
		}
	}
	
	public static String getObjectName(String dir, String filename) {
		return HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(filename)+".object", false, new int[] {0}, false)[0];		
	}
	
	public static String[][][] identifySet(MetaAnalysisParams maps, String[] files, Logger log) {
		String[][][] finalSets;
		boolean[] picks, used;
		int numMatches, index;
		String[][] phenotypes, maskedPhenos;
		String[] studies;
		String[][] races;
		
		used = Array.booleanArray(files.length, false);
		
		index = ext.indexOfStr(maps.getSnpInfoFilename(), files);
		if (index >= 0) {
			used[index] = true;
		}
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		maskedPhenos = maps.getPhenotypesWithFilenameAliases(false);
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();

		finalSets = new String[phenotypes.length][][]; // [pheno][study][race] <- all files meeting criteria
		for (int i = 0; i < phenotypes.length; i++) {
			finalSets[i] = Matrix.stringMatrix(studies.length, races.length, "<missing>");
			log.report("For "+phenotypes[i][0]+" identified:", true, true);
			log.report("\tStudy\t"+Array.toStr(Matrix.extractColumn(races, 0)));
			for (int j = 0; j < studies.length; j++) {
				log.report("\t"+studies[j], false, true);
				for (int k = 0; k < races.length; k++) {
					picks = Array.booleanArray(files.length, false);
					for (int f = 0; f < files.length; f++) {
						if (files[f].contains(studies[j]) && ext.containsAny(files[f], maskedPhenos[i]) && ext.containsAny(files[f], races[k])) {
							picks[f] = true;
							if (finalSets[i][j][k].equals("<missing>")) {
								finalSets[i][j][k] = files[f];
							} else {
								finalSets[i][j][k] += ";"+files[f];
							}
							if (used[f]) {
								log.reportError("Error - file '"+files[f]+"' matches to "+studies[j]+"/"+phenotypes[i][0]+" but was already picked for another purpose");
							}
							used[f] = true;
						}
					}
					numMatches = Array.booleanArraySum(picks);
					if (numMatches == 0) {
//						log.reportError("Warning - could not find a match for "+studies[j]+"/"+phenotypes[i][0]+"/"+races[k][0]);
					} else if (numMatches > 1) {
//						log.reportError("Error - found multiple matched for "+studies[j]+"/"+phenotypes[i][0]);
//						log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
					}
//					log.report("   "+finalSets[i][j], true, false);
					log.report("\t"+numMatches, false, true);
				}
				log.report("");
			}
			log.report("");
		}

		numMatches = Array.booleanArraySum(used);
		if (numMatches != files.length) {
			log.reportError("Warning - did not find a match for the following file(s):");
			for (int i = 0; i < files.length; i++) {
				if (!used[i]) {
					log.reportError("  "+files[i]);
				}
			}
		}
		
		return finalSets;
	}
	
	public static int getMaxChr() {
		String chrom;
		int maxChr;
		
		maxChr = 0;
		for (int chr = 1; chr <= 24; chr++) {
			chrom = chr==23?"X":(chr==24?"Y":chr+"");
			if (Files.exists("snpInfos/snpInfo_chr"+chrom+".RData")) {
				maxChr = chr;
			}
		}

		return maxChr;
	}

	public static void splitAll(String dir, MetaAnalysisParams maps) {
		String[] files;
		String[][][] finalSets;
		Logger log;
		String localDir;
		String filename;
		Vector<String> toBeSplit, commands;
		String objectName, snpInfoName, chrom, subsetObject;
		String[][] phenotypes;
		String[] studies;
		String[][] races;
		String snpInfoFile;
		String chromName;
		String geneName;
		boolean problem;
		int maxChr;
		IntVector chrsToDo;
		IntVector jobSizes;
		Vector<String> jobNames; 
		
		problem = false;
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"splitAll.log");
		new File(dir+"batchSplits/").mkdir();
		files = Files.list(dir, null, ".Rdata", false, false);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		snpInfoFile = maps.getSnpInfoFilename();
		chromName = maps.getChromName();
		geneName = maps.getGeneName();
		
		if (ext.indexOfStr(snpInfoFile, files) == -1) {
			log.reportError("Error - could not find SNP Info file '"+snpInfoFile+"'; aborting");
			return;
		}
		
		if (Files.exists(dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object")) {
			snpInfoName = getObjectName(dir, snpInfoFile);
		} else {
			log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object"+"'");
			snpInfoName = "UNKNOWN_SNP_INFO_OBJECT_NAME";
			problem = true;
		}
		
		
		commands = new Vector<String>();
		commands.add("load(\""+dir+snpInfoFile+"\")");
		commands.add("ls()");
		
		commands.add("chroms <- unique("+snpInfoName+"$"+chromName+")");
		commands.add("write.table( chroms, \"chroms.csv\", sep=\",\")");

		commands.add("for (chr in chroms) {");
		commands.add("  snps_on_chr <- "+snpInfoName+"["+snpInfoName+"$CHROM == chr,]");
		commands.add("  filename <- paste(\"snpInfos/snpInfo_chr\", chr, \".RData\", sep='')");
		commands.add("  save(snps_on_chr, file=filename, compress=\"bzip2\")");
		commands.add("}");
		
		filename = dir+"batchSplits/splitChrs.R";
		Files.writeList(Array.toStringArray(commands), filename);

		new File(dir+"snpInfos/").mkdirs();
		Files.qsub(dir+"batchSplits/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save "+filename, 5000, 0.25, 1);
		
		toBeSplit = new Vector<String>();
		toBeSplit.add("# make sure to run \"qsub "+ext.rootOf(filename)+".qsub\" first!!!");
		toBeSplit.add("cd batchSplits/");

		jobNames = new Vector<String>();
		jobSizes = new IntVector();
		
		dir = ext.verifyDirFormat(dir);
		finalSets = identifySet(maps, files, log);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				for (int k = 0; k < races.length; k++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						localDir = dir+"objects/"+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+"/";
						new File(localDir).mkdirs();
						
						files = finalSets[i][j][k].split(";");						
						for (int f = 0; f < files.length; f++) {
							commands = new Vector<String>();
							commands.add("load(\""+dir+snpInfoFile+"\")");
							commands.add("load(\""+dir+files[f]+"\")");
							if (Files.exists(dir+"batchChecks/"+ext.rootOf(files[f])+".object")) {
								objectName = getObjectName(dir, files[f]);
							} else {
								log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(files[f])+".object"+"'");
								objectName = "UNKNOWN_PREP_SCORES_OBJECT_NAME";
								problem = true;
							}
							commands.add("ls()");
							commands.add("objectType <- class("+objectName+")");
							commands.add("objectType");

							chrsToDo = new IntVector();
							maxChr = getMaxChr();
							for (int chr = 1; chr <= maxChr; chr++) {
								chrom = chr==23?"X":(chr==24?"Y":chr+"");
								subsetObject = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr"+chrom;
								if (!Files.exists(localDir+subsetObject+"_f"+f+".RData") && !Files.exists(localDir+subsetObject+".RData")) {
									chrsToDo.add(chr); 
								}
							}
							if (chrsToDo.size() != 0 && chrsToDo.size() != maxChr) {
								log.reportError("Warning - for "+studies[j]+";"+races[k][0]+"/"+phenotypes[i][0]+", missing chr(s) "+ext.listWithCommas(Array.toStringArray(chrsToDo.toArray())));
								log.reportError("        - if batch job was killed in the middle, suggest deleting the last attempted chromosome, in case it was incomplete");
							}
							
							for (int c = 0; c < chrsToDo.size(); c++) {
								chrom = chrsToDo.elementAt(c)==23?"X":(chrsToDo.elementAt(c)==24?"Y":chrsToDo.elementAt(c)+"");
								subsetObject = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr"+chrom;
								
								// filter for the gene names present on the chromosome
								commands.add("genes <- unique("+snpInfoName+"["+snpInfoName+"$"+chromName+" == \""+chrom+"\", \""+geneName+"\"])");

								// take the intersect of those actually present in the prepScores object
								commands.add("idx <- intersect(genes, names("+objectName+"))");

								// create the prepScores subset
								commands.add(subsetObject+" <- "+objectName+"[idx]");

								// make sure the dataset has the appropriate parent skatCohort or seqMeta class
								commands.add("class("+subsetObject+") <- objectType");

								// save the new file
								commands.add("save("+subsetObject+", file=\""+localDir+subsetObject+"_f"+f+".RData\", compress=\"bzip2\")");
								
								// free up memory
								commands.add("rm("+subsetObject+")");
								commands.add("");
							}

							if (chrsToDo.size() > 0) {
								filename = dir+"batchSplits/"+studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_f"+f+".R";
								Files.writeList(Array.toStringArray(commands), filename);
	
								Files.qsub(dir+"batchSplits/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save "+filename, 10000, 0.5, 1);
								toBeSplit.add("qsub "+ext.rootOf(filename)+".qsub");
								jobNames.add(dir+"batchSplits/"+ext.rootOf(filename)+".qsub");
								jobSizes.add((int)(new File(dir+files[f]).length()));
//								jobSizes.add((int)(new File(dir+files[f]).length()+chrsToDo.size()*2/maxChr*new File(dir+files[f]).length()));
							}
						}
					}
				}
			}
		}
		
		if (problem) {
			log.reportError("   need to first run using the -determineObjectNames option");
			return;
		}
		
		toBeSplit.add("# make sure to run \"SeqMeta -consolidate\" after everything else is run!!!");
		Files.writeList(Array.toStringArray(toBeSplit), dir+"master.toBeSplit");
		Files.chmod(dir+"master.toBeSplit");
		
		log.report("");
		log.report("Make sure to run \"qsub splitChrs.qsub\" first!!!");
		
		Files.qsubMultiple(jobNames, jobSizes, "chunks/", "chunkSplit", 8, true, null, -1, 22000, 2);
	}
	
	public static void consolidate(String dir, MetaAnalysisParams maps) {
		String[] files;
		String[][][] finalSets;
		Logger log;
		String localDir;
		String filename;
		String chrom, subsetObject;
		String[][] phenotypes;
		String[] studies;
		String[][] races;
		boolean problem;
		long fileSize, largestFileSize;
		int[][][][] finalSelections;
		int count;
		int maxChr;
		
		problem = false;
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"consolidateAll.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		maxChr = getMaxChr();
		
		count = 0;
		dir = ext.verifyDirFormat(dir);
		finalSets = identifySet(maps, files, log);
		finalSelections = new int[finalSets.length][finalSets[0].length][finalSets[0][0].length][];
		for (int iter = 0; iter < 2; iter++) {
			for (int i = 0; i < phenotypes.length; i++) {
				System.out.println(phenotypes[i][0]);
				for (int j = 0; j < studies.length; j++) {
					System.out.println("  "+studies[j]);
					for (int k = 0; k < races.length; k++) {
						if (!finalSets[i][j][k].equals("<missing>")) {
							System.out.println("    "+races[k][0]);
							localDir = dir+"objects/"+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+"/";
							files = finalSets[i][j][k].split(";");
							if (iter == 0) {
								finalSelections[i][j][k] = Array.intArray(maxChr, -1);
							}
							
							for (int chr = 1; chr <= maxChr; chr++) {

								chrom = chr==23?"X":(chr==24?"Y":chr+"");
								subsetObject = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr"+chrom;
								
								largestFileSize = 0;
								for (int f = 0; f < files.length; f++) {
									filename = subsetObject+"_f"+f+".RData";
									
									if (iter == 0) {
										if (Files.exists(localDir+filename)) {
											fileSize = new File(localDir+filename).length();
											if (fileSize > largestFileSize) {
												largestFileSize = fileSize;
												finalSelections[i][j][k][chr-1] = f;
											}
											count++;
//										} else {
										} else if (!Files.exists(localDir+subsetObject+".RData")) {
											System.err.println("Error - could not find '"+subsetObject+"_f"+f+".RData' in "+localDir);
											problem = true;
										}
									} else {
										if (f == finalSelections[i][j][k][chr-1]) {
											new File(localDir+filename).renameTo(new File(localDir+subsetObject+".RData"));
										} else {
											new File(localDir+filename).delete();
										}
									}
								}
								

							}
						}
					}
				}
			}
			if (iter == 0) {
				if (problem) {
					if (count == 0) {
						log.reportError("\n   discrepancies found; possible explanations are that the R parser was not run or did not complete, or the original .RData files were either removed or added to; if the latter, then suggest rerunning the parsers for those study/pheno pairs; no consolidation will occur\n");
					} else {
						log.reportError("\n   did not find a single file with a _f# extension; either nothing has been run or everything has already been processed by this algorithm; check the objects/ directory\n");
					}
					return;
				} else {
					log.reportError("\nEverything seems to be in order, proceeding with consolidation\n");
				}
			}
		}
	}

	public static void runAll(String dir, MetaAnalysisParams maps, boolean forceMeta) {
		String[] files;
		String[][][] finalSets;
		Logger log;
		String localDir;
		String root, filename, objectFilename, outputFilename;
		Vector<String> toBeRunIndividually, toBeRunMetad, commands, objects;
		int count;
		String originalObjectName, objectName, snpInfoFile, snpInfoName, chrom;
		String[][] phenotypes, races;
		String[] studies;
		String snpName;
		String[][] methods;
		String functionFlagName, geneName;
		boolean runningByChr;
		IntVector jobSizes;
		Vector<String> jobNames; 
		int[] infoSizes;
		int maxChr;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"runAll.log");
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		snpName = maps.getVariantName();
		methods = maps.getMethods();
		functionFlagName = maps.getFunctionFlagName();
		geneName = maps.getGeneName();
		runningByChr = maps.runningByChr();
		snpInfoFile = maps.getSnpInfoFilename();
		
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(maps, files, log);
		
		maxChr = getMaxChr();
		jobSizes = new IntVector();
		jobNames = new Vector<String>();
		infoSizes = new int[maxChr+2];

		if (runningByChr) {
			for (int chr = 1; chr <= maxChr; chr++) {
				chrom = chr==23?"X":(chr==24?"Y":chr+"");
				filename = "snpInfos/snpInfo_chr"+chrom+".RData";
				if (!Files.exists(filename)) {
					log.reportError("Error - could not find SNP Info file '"+filename+"'; aborting");
					return;
				} else {
					infoSizes[chr] = (int)new File(filename).length();
				}
				
			}
			snpInfoName = "snps_on_chr";
		} else {
			if (!Files.exists(snpInfoFile)) {
				log.reportError("Error - could not find SNP Info file '"+snpInfoFile+"'; aborting");
				return;
			}
			if (Files.exists(dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object")) {
				snpInfoName = getObjectName(dir, snpInfoFile);
			} else {
				log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object"+"'");
				return;
			}
		}
		
		toBeRunIndividually = new Vector<String>();
		toBeRunIndividually.add("cd batchRuns/");
		toBeRunMetad = new Vector<String>();
		toBeRunMetad.add("cd batchRuns/");
		new File(dir+"batchRuns/").mkdir();
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < methods.length; j++) {
				localDir = dir+phenotypes[i][0]+"/"+methods[j][0]+"/";
				new File(localDir).mkdirs();
				for (int k = 0; k < races.length; k++) {
					localDir = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[j][0]+"/";
					new File(localDir).mkdirs();
				}
			}
		}
		
		// Primary analysis by study/race
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				for (int k = 0; k < races.length; k++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						localDir = dir+"objects/"+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+"/";
						
						for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
							chrom = chr==23?"X":(chr==24?"Y":chr+"");

							if (runningByChr) {
								objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr"+chrom;
								objectFilename = localDir+objectName+".RData";
								snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";
							} else {
								objectFilename = dir+finalSets[i][j][k];
								if (objectFilename.contains(";")) {
									log.reportError("Error - more than one file is mapped to "+studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+": "+objectFilename);
									return;
								} else if (Files.exists(dir+"batchChecks/"+ext.rootOf(objectFilename)+".object")) {
									objectName = getObjectName(dir, objectFilename);
								} else {
									log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(objectFilename)+".object"+"'");
									return;
								}
							}
							
							if (!Files.exists(objectFilename)) {
								log.reportError("Error - missing object file: '"+objectFilename+"'");
								if (Files.exists(ext.addToRoot(objectFilename, "_f0"))) {
									log.reportError("     - however did find '"+ext.removeDirectoryInfo(ext.addToRoot(objectFilename, "_f0"))+"'; so try running -consolidate");
									return;
								}
							}

							commands = new Vector<String>();
							commands.add("library(bdsmatrix)");
							commands.add("library(seqMeta)");
							commands.add("load(\""+dir+snpInfoFile+"\")");
							commands.add("load(\""+objectFilename+"\")");
							commands.add("ls()");
							if (!runningByChr) {
								originalObjectName = objectName;
								objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0];
								commands.add(objectName+" <- "+originalObjectName);
								commands.add("rm(\""+originalObjectName+"\")");
								commands.add("ls()");
							}
							
							count = 0;
							for (int m = 0; m < methods.length; m++) {
								root = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0];
								outputFilename = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/"+root+(runningByChr?"_chr"+chrom:"")+".csv";
								if (!Files.exists(outputFilename) || new File(outputFilename).length() == 0) {
									if (new File(objectFilename).length() > 1024) {
										commands.add("results <- "+methods[m][2]+"("+
												objectName+
												", SNPInfo="+(SINGLE_VARIANTS[ext.indexOfStr(methods[m][2], ALGORITHMS)]||functionFlagName==null
													?snpInfoName
													:"subset("+snpInfoName+", "+functionFlagName+"==TRUE)")+
												", snpNames = \""+snpName+"\""+
												", aggregateBy=\""+geneName+"\""+
												(methods[m].length > 3 && ext.isValidDouble(methods[m][3])
														?", mafRange = c(0,"+methods[m][3]+")"+(methods[m].length>4?", "+Array.toStr(Array.subArray(methods[m],  4), ", "):"")
														:(methods[m].length>3?", "+Array.toStr(Array.subArray(methods[m],  3), ", "):""))+")");
										commands.add("write.table( results, \""+outputFilename+"\", sep=\",\", row.names = F)");
										count++;
									} else {
										Files.write(Array.toStr(getHeaderForMethod(methods[m]), ","), outputFilename);
									}
								}
							}
							if (count > 0) {
								count = 0;
								do {
									filename = dir+"batchRuns/"+studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+(runningByChr?"_chr"+chrom:"")+(count==0?"":"_"+count)+".R";
									count++;
								} while (Files.exists(filename));
								Files.writeList(Array.toStringArray(commands), filename);
			
								Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save "+filename, 5000, 1, 1);
								toBeRunIndividually.add("qsub "+ext.rootOf(filename)+".qsub");
								jobNames.add(dir+"batchRuns/"+ext.rootOf(filename)+".qsub");
								jobSizes.add(infoSizes[chr]);
							}
						}
					}
				}
			}
		}

		Files.writeList(Array.toStringArray(toBeRunIndividually), dir+"master.toBeRunIndividually");
		Files.chmod(dir+"master.toBeRunIndividually");
		System.err.println("qsubing multiple individual runs");
		Files.qsubMultiple(jobNames, jobSizes, "chunks/", "chunkRun", 16, true, "sb", -1, 62000, 2);
		System.err.println("multiple individual runs done");

		
		jobNames = new Vector<String>();
		jobSizes = new IntVector();
		// Meta-analysis stratified by race
		for (int i = 0; i < phenotypes.length; i++) {
			for (int k = 0; k < races.length; k++) {
				for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
					chrom = chr==23?"X":(chr==24?"Y":chr+"");
					commands = new Vector<String>();
					commands.add("library(bdsmatrix)");
					commands.add("library(seqMeta)");
					if (runningByChr) {
						snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";
					}
					commands.add("load(\""+dir+snpInfoFile+"\")");
					
					objects = new Vector<String>();
					for (int j = 0; j < studies.length; j++) {
						if (!finalSets[i][j][k].equals("<missing>")) {
							if (runningByChr) {
								localDir = dir+"objects/"+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+"/";
								objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr"+chrom;
								if (new File(localDir+objectName+".RData").length() > 1024) {
									commands.add("load(\""+localDir+objectName+".RData"+"\")");
									objects.add(objectName);
								}
							} else {
								objectFilename = finalSets[i][j][k];
								originalObjectName = getObjectName(dir, objectFilename);
								objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0];
								if (new File(dir+objectFilename).length() > 1024) {
									commands.add("load(\""+dir+objectFilename+"\")");
									commands.add(objectName+" <- "+originalObjectName);
									commands.add("rm(\""+originalObjectName+"\")");
									objects.add(objectName);
								}
							}
						}
					}
					commands.add("ls()");
					commands.add("");
					count = 0;
					for (int m = 0; m < methods.length; m++) {
						root = races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0];
						outputFilename = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/"+root+(runningByChr?"_chr"+chrom:"")+".csv";
						if (forceMeta || !Files.exists(outputFilename) || new File(outputFilename).length() == 0) {
							if (objects.size() > 0) {
								commands.add("results <- "+methods[m][2]+"("+Array.toStr(
										Array.toStringArray(objects), ", ")+
										", SNPInfo="+(SINGLE_VARIANTS[ext.indexOfStr(methods[m][2], ALGORITHMS)]||functionFlagName==null
											?snpInfoName
											:"subset("+snpInfoName+", "+functionFlagName+"==TRUE)")+
										", snpNames = \""+snpName+"\""+
										", aggregateBy=\""+geneName+"\""+
										(methods[m].length > 3 && ext.isValidDouble(methods[m][3])
												?", mafRange = c(0,"+methods[m][3]+")"+(methods[m].length>4?", "+Array.toStr(Array.subArray(methods[m],  4), ", "):"")
														:(methods[m].length>3?", "+Array.toStr(Array.subArray(methods[m],  3), ", "):""))+")");
								commands.add("write.table( results, \""+outputFilename+"\", sep=\",\", row.names = F)");
								commands.add("");
								count++;
							} else {
								Files.write(Array.toStr(getHeaderForMethod(methods[m]), ","), outputFilename);
							}

						}
					}
					if (count > 0) {
						count = 0;
						do {
							filename = dir+"batchRuns/"+races[k][0]+"_"+phenotypes[i][0]+(runningByChr?"_chr"+chrom:"")+(count==0?"":"_"+count)+".R";
							count++;
						} while (Files.exists(filename));
						Files.writeList(Array.toStringArray(commands), filename);
			
						Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save "+filename, 25000, 2, 1);
						toBeRunMetad.add("qsub "+ext.rootOf(filename)+".qsub");
						jobNames.add(dir+"batchRuns/"+ext.rootOf(filename)+".qsub");
						jobSizes.add(infoSizes[chr]);
					}
				}
			}

			// Meta-analysis of all races
			for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
				chrom = chr==23?"X":(chr==24?"Y":chr+"");
				commands = new Vector<String>();
				commands.add("library(bdsmatrix)");
				commands.add("library(seqMeta)");
				if (runningByChr) {
					snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";
				}
				commands.add("load(\""+dir+snpInfoFile+"\")");
				
				objects = new Vector<String>();
				for (int j = 0; j < studies.length; j++) {
					for (int k = 0; k < races.length; k++) {
						if (!finalSets[i][j][k].equals("<missing>")) {
							if (runningByChr) {
								localDir = dir+"objects/"+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+"/";
								objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr"+chrom;
								if (new File(localDir+objectName+".RData").length() > 1024) {
									commands.add("load(\""+localDir+objectName+".RData"+"\")");
									objects.add(objectName);
								}
							} else {
								objectFilename = finalSets[i][j][k];
								originalObjectName = getObjectName(dir, objectFilename);
								objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0];
								if (new File(dir+objectFilename).length() > 1024) {
									commands.add("load(\""+dir+objectFilename+"\")");
									commands.add(objectName+" <- "+originalObjectName);
									commands.add("rm(\""+originalObjectName+"\")");
									objects.add(objectName);
								}
							}
						}
					}
				}
				commands.add("ls()");
				commands.add("");
				count = 0;
				for (int m = 0; m < methods.length; m++) {
					root = phenotypes[i][0]+"_"+methods[m][0];
					outputFilename = dir+phenotypes[i][0]+"/"+methods[m][0]+"/"+root+(runningByChr?"_chr"+chrom:"")+".csv";
					if (forceMeta || !Files.exists(outputFilename) || new File(outputFilename).length() == 0) {
						if (objects.size() > 0) {
							commands.add("results <- "+methods[m][2]+"("+Array.toStr(
									Array.toStringArray(objects), ", ")+
									", SNPInfo="+(SINGLE_VARIANTS[ext.indexOfStr(methods[m][2], ALGORITHMS)]||functionFlagName==null
										?snpInfoName
										:"subset("+snpInfoName+", "+functionFlagName+"==TRUE)")+
									", snpNames = \""+snpName+"\""+
									", aggregateBy=\""+geneName+"\""+
									(methods[m].length > 3 && ext.isValidDouble(methods[m][3])
											?", mafRange = c(0,"+methods[m][3]+")"+(methods[m].length>4?", "+Array.toStr(Array.subArray(methods[m],  4), ", "):"")
													:(methods[m].length>3?", "+Array.toStr(Array.subArray(methods[m],  3), ", "):""))+")");
							commands.add("write.table( results, \""+outputFilename+"\", sep=\",\", row.names = F)");
							commands.add("");
							count++;
						} else {
							Files.write(Array.toStr(getHeaderForMethod(methods[m]), ","), outputFilename);
						}
					}
				}
				if (count > 0) {
					count = 0;
					do {
						filename = dir+"batchRuns/"+phenotypes[i][0]+(runningByChr?"_chr"+chrom:"")+(count==0?"":"_"+count)+".R";
						count++;
					} while (Files.exists(filename));
					Files.writeList(Array.toStringArray(commands), filename);
		
					Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save "+filename, 30000, 2, 1);
					toBeRunMetad.add("qsub "+ext.rootOf(filename)+".qsub");
					jobNames.add(dir+"batchRuns/"+ext.rootOf(filename)+".qsub");
					jobSizes.add(infoSizes[chr]);
				}
			}
		}
		Files.writeList(Array.toStringArray(toBeRunMetad), dir+"master.toBeMetaAnalyzed");
		Files.chmod(dir+"master.toBeMetaAnalyzed");
		System.err.println("qsubing multiple meta runs");
		Files.qsubMultiple(jobNames, jobSizes, "chunks/", "chunkMeta", 16, true, "sb", -1, 62000, 2);
		System.err.println("multiple meta runs done");
	}

	public static String[] getHeaderForMethod(String[] method) {
		return HEADER_TYPES[ext.indexOfStr(method[2], ALGORITHMS)];
	}

	public static void parseAll(String dir, MetaAnalysisParams maps, boolean forceMeta) {
		String[] files;
		String[][][] finalSets;
		Logger log;
		String root;
		String[][] phenotypes, races, methods;
		String[] studies;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		
		log = new Logger(dir+"parseAll.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(maps, files, log);
		
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				for (int k = 0; k < races.length; k++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						for (int m = 0; m < methods.length; m++) {
							root = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0];
							if (!Files.exists(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/"+root+".csv") || new File(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/"+root+".csv").length() == 0) {
								stitch(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/", root+"_chr#.csv", root+".csv", log);
							}
						}
					}
				}
			}
			
			for (int k = 0; k < races.length; k++) {
				for (int m = 0; m < methods.length; m++) {
					root = races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0];
					if (forceMeta || !Files.exists(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/"+root+".csv") || new File(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/"+root+".csv").length() == 0) {
						stitch(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/", root+"_chr#.csv", root+".csv", log);
					}
				}
			}
			
			for (int m = 0; m < methods.length; m++) {
				root = phenotypes[i][0]+"_"+methods[m][0];
				if (forceMeta || !Files.exists(dir+phenotypes[i][0]+"/"+methods[m][0]+"/"+root+".csv") || new File(dir+phenotypes[i][0]+"/"+methods[m][0]+"/"+root+".csv").length() == 0) {
					stitch(dir+phenotypes[i][0]+"/"+methods[m][0]+"/", root+"_chr#.csv", root+".csv", log);
				}
			}
			
		}
	}

	public static void parseMetrics(String dir, MetaAnalysisParams maps) {
		PrintWriter writer;
		String[] line;
		String[] files;
		String[][][] finalSets;
		Logger log;
		String root, outputFilename;
		String[][] phenotypes, races;
		String[] studies;
		String[][] methods;
		String method;
		String[][][] sampleSizes;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"parseMetrics.log");
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(maps, files, log);
		sampleSizes = new String[races.length][studies.length][phenotypes.length];

		method = null;
		for (int j = 0; j < methods.length; j++) {
			if (methods[j][2].equals("singlesnpMeta")) {
				method = methods[j][0];
			}
		}
		if (method == null) {
			log.report("Error - cannot summarize if there is no method using singlesnpMeta");
			return;
		}

		try {
			writer = new PrintWriter(new FileWriter(dir+"summary_metrics.xln"));
			writer.println("Study\tRace\tPhenotype\tN_samples\tn_genotyped\tn_MAF1%\tLambda_MAF1%\tbetaMean_MAF1%\tbetaSD_MAF1%\tn_20count\tLambda_20count\tbetaMean_20count\tbetaSD_20count\tn_MAF5%\tLambda_MAF5%\tbetaMean_MAF5%\tbetaSD_MAF5%");
			for (int j = 0; j < studies.length; j++) {
				for (int k = 0; k < races.length; k++) {
					for (int i = 0; i < phenotypes.length; i++) {
						if (finalSets[i][j][k].equals("<missing>")) {
							sampleSizes[k][j][i] = "";
						} else {
							root = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+method;
							outputFilename = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+method+"/"+root+".csv";
							if (Files.exists(outputFilename)) {
								line = procFile(outputFilename, log);
								sampleSizes[k][j][i] = line[0];
								writer.println(studies[j]+"\t"+races[k][0]+"\t"+phenotypes[i][0]+"\t"+Array.toStr(line));
							} else {
								log.report("Error - missing expected file: "+outputFilename);
							}
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"summary_metrics.xln");
			e.printStackTrace();
		}
		
		
		Vector<int[]> ranges;
		int start, stop;
		
		start = 1;
		ranges = new Vector<int[]>();
		try {
			writer = new PrintWriter(new FileWriter(dir+"sampleSizes.xln"));
			for (int k = 0; k < races.length; k++) {
				writer.println(races[k][0]+"\t"+Array.toStr(Matrix.extractColumn(phenotypes, 0)));
				start++;
				stop = start;
				for (int j = 0; j < studies.length; j++) {
					writer.println(studies[j]+"\t"+Array.toStr(sampleSizes[k][j]));
					stop++;
				}
				writer.print("Subtotal");
				for (int i = 0; i < phenotypes.length; i++) {
					writer.print("\t=SUM("+ext.getExcelColumn(i+1)+start+":"+ext.getExcelColumn(i+1)+(stop-1)+")");
				}
				ranges.add(new int[] {start, stop-1});
				writer.println();
				writer.println();
				start = stop + 2;
			}

			writer.println();
			writer.print("Total");
			for (int i = 0; i < phenotypes.length; i++) {
				writer.print("\t=SUM(");
				for (int j = 0; j < ranges.size(); j++) {
					writer.print((j==0?"":",")+ext.getExcelColumn(i+1)+ranges.elementAt(j)[0]+":"+ext.getExcelColumn(i+1)+ranges.elementAt(j)[1]);
				}
				writer.print(")");
			}
			writer.println();
			
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"sampleSizes.xln");
			e.printStackTrace();
		}
	}

	public static String[] procFile(String outputFilename, Logger log) {
		String[] metrics;
		
		BufferedReader reader;
		String[] line;
		String temp;
		int count;
		String delimiter;
		DoubleVector[] dvs;
		int maxSamples;
		
		count = 0;
		maxSamples = 0;
		metrics = Array.stringArray(14, "");
		dvs = DoubleVector.newDoubleVectors(6);
		
		try {
			reader = new BufferedReader(new FileReader(outputFilename));
			temp = ext.replaceAllWith(reader.readLine(), "\"", "");
			delimiter = ext.determineDelimiter(temp);

			if (ext.checkHeader(temp.split(delimiter, -1), HEADER_TYPES[0], Array.intArray(HEADER_TYPES[0].length), false, log, false)) {
				while (reader.ready()) {
					if (delimiter.equals(",")) {
						line = ext.splitCommasIntelligently(reader.readLine().trim(), true, log);
					} else {
						line = reader.readLine().trim().split(delimiter, -1);
					}
					try {
						if (!line[5].equals("0")) {
							count++;
							if (Integer.parseInt(line[5]) > maxSamples) {
								maxSamples = Integer.parseInt(line[5]);
							}
							if (Double.parseDouble(line[3]) >= 0.01) {
								dvs[0].add(Double.parseDouble(line[2]));
								dvs[1].add(Math.abs(Double.parseDouble(line[6])));
							}
							if (Double.parseDouble(line[3])*Double.parseDouble(line[5]) >= 10) {
								dvs[2].add(Double.parseDouble(line[2]));
								dvs[3].add(Math.abs(Double.parseDouble(line[6])));
							}
							if (Double.parseDouble(line[3]) >= 0.05) {
								dvs[4].add(Double.parseDouble(line[2]));
								dvs[5].add(Math.abs(Double.parseDouble(line[6])));
							}
						}
					} catch (Exception e) {
						log.reportError("Error parsing line: "+Array.toStr(line, delimiter));
					}
				}
				metrics[0] = maxSamples+"";
				metrics[1] = count+"";
				metrics[2] = dvs[0].size()+"";
				metrics[3] = ext.formDeci(ProbDist.ChiDistReverse(Array.median(dvs[0].toArray()), 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
				metrics[4] = ext.formDeci(Array.mean(dvs[1].toArray()), 6);
				metrics[5] = ext.formDeci(Array.stdev(dvs[1].toArray()), 6);
				metrics[6] = dvs[2].size()+"";
				metrics[7] = ext.formDeci(ProbDist.ChiDistReverse(Array.median(dvs[2].toArray()), 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
				metrics[8] = ext.formDeci(Array.mean(dvs[3].toArray()), 6);
				metrics[9] = ext.formDeci(Array.stdev(dvs[3].toArray()), 6);
				metrics[10] = dvs[4].size()+"";
				metrics[11] = ext.formDeci(ProbDist.ChiDistReverse(Array.median(dvs[4].toArray()), 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
				metrics[12] = ext.formDeci(Array.mean(dvs[5].toArray()), 6);
				metrics[13] = ext.formDeci(Array.stdev(dvs[5].toArray()), 6);
			} else {
				log.reportError("Error - unexpected header for file '"+outputFilename+"' : "+temp);
				System.exit(1);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + outputFilename + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + outputFilename + "\"");
		}
		
		
		return metrics;
	}

	public static void doubleCheckNs(String dir, MetaAnalysisParams maps) {
		String[] files;
		String[][][] finalSets;
		Logger log;
		String localDir;
		Vector<String> commands;
		String objectName;
		String[][] phenotypes, races;
		String[] studies;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();

		log = new Logger(dir+"checkNs.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(maps, files, log);
		
		commands = new Vector<String>();
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				for (int k = 0; k < races.length; k++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						localDir = dir+"objects/"+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+"/";
						objectName = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_chr18";
						commands.add("load(\""+localDir+objectName+".RData"+"\")");
						commands.add("ls()");
						commands.add("head("+objectName+"$RBFA$n)");
						commands.add("head("+objectName+"$PTPRM$n)");
						commands.add("head("+objectName+"$SOGA2$n)");
						commands.add("head("+objectName+"$RGL3$n)");
						commands.add("head("+objectName+"$MAN2B1$n)");
						commands.add("remove("+objectName+")");
						commands.add("ls()");
					}
				}
			}
		}
		
		Files.writeList(Array.toStringArray(commands), dir+"checkNs.R");
		Files.qsub(dir+"checkNs.qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save checkNs.R", 5000, 1, 1);
	}
	
	public static void computeAllRelevantMACs(String dir, MetaAnalysisParams maps, Logger log) {
		String[][] methods;
		Vector<String> v;
		
		v = new Vector<String>();
		methods = maps.getMethods();
		
		for (int m = 0; m < methods.length; m++) {
			if (methods[m].length >= 3 && ext.isValidDouble(methods[m][3])) {
				HashVec.addIfAbsent(methods[m][3], v);
			}
		}
		
		log.report("Computing MACs for the following mafThresholds: "+ext.listWithCommas(Array.toStringArray(v)));
		for (int i = 0; i < v.size(); i++) {
			log.report("Starting MACs for maf<="+v.elementAt(i));
			computeMAC(dir, maps, v.elementAt(i), log);
			log.report("");
		}
	}

	public static void computeMAC(String dir, MetaAnalysisParams maps, String mafThreshold, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		long time;
		String[] files;
		String[][][] finalSets;
		String localDir;
		String filename;
		String[] header;
		int[] indices;
		Hashtable<String, String> snpGeneHash;
		Hashtable<String, String> snpGeneFunctionalHash;
		Hashtable<String, Vector<String>> geneLoci;
		String[] keys;
		Hashtable<String, int[]> macs, raceSpecificMacs;
		String gene;
		int[] counts;
		double mafThresholdDouble;
		String[][] phenotypes, races, methods;
		String[] studies;
		String snpInfoFile, functionFlagName;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}

		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		snpInfoFile = maps.getSnpInfoFilename();
		functionFlagName = maps.getFunctionFlagName();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		
		dir = ext.verifyDirFormat(dir);
		mafThresholdDouble = Double.parseDouble(mafThreshold);
		
		time = new Date().getTime();
		filename = dir+ext.rootOf(snpInfoFile)+".csv";
		if (!Files.exists(filename)) {
			Files.writeList(new String[] {"load(\""+dir+snpInfoFile+"\")",
					"write.table("+getObjectName(dir, snpInfoFile)+", \""+filename+"\", sep=\",\", row.names = F)",
					}, dir+"dump_"+ext.rootOf(snpInfoFile)+".R");
			Files.qsub("dump_"+ext.rootOf(snpInfoFile)+".qsub", "cd "+dir+"\n"+getRscriptExecutable(maps, log)+" --no-save dump_"+ext.rootOf(snpInfoFile)+".R", 5000, 0.5, 1);
			log.reportError("Error - need to dump snpInfoFile by running 'qsub dump_"+ext.rootOf(snpInfoFile)+".qsub' first");
			return;
		}
		if (Files.exists(filename+".mappings.ser") && Files.exists(filename+".maf"+mafThreshold+".functionalMappings.ser")) {
			snpGeneHash = SerialHash.loadSerializedStringHash(filename+".mappings.ser");
			snpGeneFunctionalHash = SerialHash.loadSerializedStringHash(filename+".maf"+mafThreshold+".functionalMappings.ser");
			log.report(ext.getTime()+"\tReloaded marker mappings in " + ext.getTimeElapsed(time));
		} else {
			snpGeneHash = new Hashtable<String, String>();
			snpGeneFunctionalHash = new Hashtable<String, String>();
			geneLoci = new Hashtable<String, Vector<String>>();

			try {
				reader = new BufferedReader(new FileReader(filename));
				header = ext.splitCommasIntelligently(reader.readLine(), true, log);
				indices = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.GENE_UNITS, new String[] {functionFlagName}, Aliases.CHRS, {"MAF"}}, header, false, true, true, log, true);
				while (reader.ready()) {
					line = ext.splitCommasIntelligently(reader.readLine(), true, log);
					snpGeneHash.put(line[indices[0]], line[indices[1]]);
					if (line[indices[2]].equals("TRUE") && !line[indices[4]].equals("NA") && Double.parseDouble(line[indices[4]]) <= mafThresholdDouble ) {
						snpGeneFunctionalHash.put(line[indices[0]], line[indices[1]]);
					}
					if (!geneLoci.containsKey(line[indices[1]])) {
						geneLoci.put(line[indices[1]], new Vector<String>());
					}
					HashVec.addIfAbsent(line[indices[3]], geneLoci.get(line[indices[1]]));
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + filename + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + filename + "\"");
				System.exit(2);
			}
			
			log.report(ext.getTime()+"\tLoaded data in " + ext.getTimeElapsed(time));
			
			keys = HashVec.getKeys(geneLoci);
			for (int i = 0; i < keys.length; i++) {
				if (geneLoci.get(keys[i]).size() > 1) {
					log.reportError("Gene '"+keys[i]+"' can be found on chromosomes "+ext.listWithCommas(Array.toStringArray(geneLoci.get(keys[i])), true));
				}
			}
			
			SerialHash.createSerializedStringHash(filename+".mappings.ser", snpGeneHash);
			SerialHash.createSerializedStringHash(filename+".maf"+mafThreshold+".functionalMappings.ser", snpGeneFunctionalHash);

			log.report(ext.getTime()+"\tFinished mapping markers to genes in " + ext.getTimeElapsed(time));
		}
		
		
		if (!methods[0][0].equals("SingleSNP")) {
			System.err.println("Error - this program erroneously assumed that the first model was SingleSNP and got confused (it's actually "+methods[0][0]+"); aborting");
			return;
		}

		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(maps, files, log);

		for (int i = 0; i < phenotypes.length; i++) {
			macs = new Hashtable<String, int[]>();
			for (int k = 0; k < races.length; k++) {
				raceSpecificMacs = new Hashtable<String, int[]>();
				localDir = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[0][0]+"/";
				for (int j = 0; j < studies.length; j++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						filename = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+methods[0][0]+".csv";
						log.report(ext.getTime()+"\tReading "+filename);
						
						try {
							reader = new BufferedReader(new FileReader(localDir+filename));
							header = ext.splitCommasIntelligently(reader.readLine(), true, log);
//							ext.checkHeader(header, HEADER_TYPES[Integer.parseInt(MODELS[0][4])], Array.intArray(expected.length), false, log, true);
							
							indices = ext.indexFactors(new String[] {"Name", "maf", "ntotal"}, header, false, true);
							
							while (reader.ready()) {
								line = ext.splitCommasIntelligently(reader.readLine(), true, log);
								if (!snpGeneHash.containsKey(line[indices[0]])) {
									log.reportError("Warning - variant '"+line[indices[0]]+"' was not found in the snpInfo file");
								}
								if (snpGeneFunctionalHash.containsKey(line[indices[0]])) {
									gene = snpGeneFunctionalHash.get(line[indices[0]]);
									
									// pan-ethnic
									if (macs.containsKey(gene)) {
										counts = macs.get(gene);
									} else {
										macs.put(gene, counts = new int[studies.length]);
									}
									if (!line[indices[1]].equals("NA")) {
										counts[j] += Math.round(Double.parseDouble(line[indices[1]]) * Double.parseDouble(line[indices[2]]) * 2);
									}

									// race-specific
									if (raceSpecificMacs.containsKey(gene)) {
										counts = raceSpecificMacs.get(gene);
									} else {
										raceSpecificMacs.put(gene, counts = new int[studies.length]);
									}
									if (!line[indices[1]].equals("NA")) {
										counts[j] += Math.round(Double.parseDouble(line[indices[1]]) * Double.parseDouble(line[indices[2]]) * 2);
									}
								}
							}
							reader.close();
						} catch (FileNotFoundException fnfe) {
							System.err.println("Error - could not find '"+localDir+filename+"'; aborting");
							return;
						} catch (IOException ioe) {
							System.err.println("Error reading file \"" + filename + "\"");
							return;
						}
					}
				}
				log.report("", true, false);

				try {
					writer = new PrintWriter(new FileWriter(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+"minorAlleleCounts.maf"+mafThreshold+".xln"));
					keys = HashVec.getKeys(raceSpecificMacs);
					writer.println("Gene\t"+Array.toStr(studies)+"\tTotal");
					for (int m = 0; m < keys.length; m++) {
						counts = raceSpecificMacs.get(keys[m]);
						writer.println(keys[m]+"\t"+Array.toStr(counts)+"\t"+Array.sum(counts));
					}
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + dir+phenotypes[i][0]+"/"+races[k][0]+"/"+"minorAlleleCounts.maf"+mafThreshold+".xln");
					e.printStackTrace();
				}
			}

			try {
				writer = new PrintWriter(new FileWriter(dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+mafThreshold+".xln"));
				keys = HashVec.getKeys(macs);
				writer.println("Gene\t"+Array.toStr(studies)+"\tTotal");
				for (int m = 0; m < keys.length; m++) {
					counts = macs.get(keys[m]);
					writer.println(keys[m]+"\t"+Array.toStr(counts)+"\t"+Array.sum(counts));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+mafThreshold+".xln");
				e.printStackTrace();
			}
		}
	}

	public static void assembleHits(String dir, MetaAnalysisParams maps, int macThresholdStudy, int macThresholdTotal) {
		String[] files;
		String[][][] finalSets;
		Logger log;
		String localDir, localRaceDir;
		Hashtable<String,Hits> groupHits;
		Hashtable<String,Vector<String>> groupParams;
		String[] groups, hits;
		String filename, pvalFile;
		String[] header;
		String filenames;
		Hashtable<String,Hashtable<String,Hashtable<String,String>>> macHashesHashByRace;
		Hashtable<String,String> macHash;
		String[] line;
		String temp;
		String[][] phenotypes, races;
		String[] studies;
		String[][] methods;
		String[][] groupAnnotationParams;
		int count;
		Vector<String> lineCounts;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		
		//String[][] phenotypes, String[] studies, int[][] defaultSampleSizes, String[][] groupAnnotationParams, String snpInfoFile, 

		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		groupAnnotationParams = maps.getGroupAnnotationParams();

		log = new Logger(dir+"assembleHits.log");

		// match all rdata files to pheno/study/race
		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(maps, files, log);

		lineCounts = new Vector<String>();
		for (int i = 0; i < phenotypes.length; i++) {
			groupHits = new Hashtable<String, Hits>();
			groupParams = new Hashtable<String, Vector<String>>();
			
//			MODELS = { // name, grouping, subroutine, arguments, header type, mafThreshold //, parameters for parsing
//			log.reportError("Error - a method must have at least 3 parameters: name, grouping, algorithm, (optional) MAF threshold, (optional) additional arguments such as weighting");
			
			macHashesHashByRace = new Hashtable<String, Hashtable<String,Hashtable<String,String>>>();
			for (int m = 0; m < methods.length; m++) {
				if (!groupHits.containsKey(methods[m][1])) {
					groupHits.put(methods[m][1], new Hits());
					groupParams.put(methods[m][1], new Vector<String>());
				}
			}
			for (int k = 0; k < races.length; k++) {
				macHashesHashByRace.put(races[k][0], new Hashtable<String, Hashtable<String,String>>());
			}
			macHashesHashByRace.put("PanEthnic", new Hashtable<String, Hashtable<String,String>>());
			for (int m = 0; m < methods.length; m++) {
				if (methods[m][1].equals("BurdenTests")) {
					// see if particular maf threshold has been introduced prior
					if (!macHashesHashByRace.get("PanEthnic").containsKey(methods[m][3])) {
						line = Array.addStrToArray("Total", studies);
						for (int k = 0; k < races.length; k++) {
							macHashesHashByRace.get(races[k][0]).put(methods[m][3], macHash = HashVec.loadFileToHashString(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+"minorAlleleCounts.maf"+methods[m][3]+".xln", "Gene", line, "\t"));
							macHash.put("studies", Array.toStr(line));
						}
						macHashesHashByRace.get("PanEthnic").put(methods[m][3], macHash = HashVec.loadFileToHashString(dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+methods[m][3]+".xln", "Gene", new String[] {"Total"}, "\t"));
						macHash.put("studies", "Total");
					}
				}
				
				filenames = "";
				localDir = dir+phenotypes[i][0]+"/"+methods[m][0]+"/";
				for (int k = 0; k < races.length; k++) {
					localRaceDir = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/";

					for (int j = 0; j < studies.length; j++) {
						if (!finalSets[i][j][k].equals("<missing>")) {
							filename = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0]+".csv";
							pvalFile = studies[j]+"_"+races[k][0]+"_pvals.dat";

							count = Files.parsePvals(localRaceDir+filename, localDir+pvalFile, studies[j], methods[m], macHashesHashByRace.get(races[k][0]), macThresholdStudy, log);
							if (count == -1) {
								return;
							}
							lineCounts.add(studies[j]+"\t"+races[k][0]+"\t"+phenotypes[i][0]+"\t"+methods[m][0]+"\t"+count);


							filenames += pvalFile+",1="+studies[j]+"_"+races[k][0]+";";
							
							header = Files.getHeaderOfFile(localRaceDir+filename, ",!", log);
							header[0] = "'"+header[0]+"'";
							for (int h = 1; h < header.length; h++) {
								header[h] = "'"+header[h]+"'="+studies[j]+"_"+races[k][0]+"_"+header[h]+"_"+methods[m][0];
							}
							if (methods[m][1].equals("BurdenTests")) {
								groupParams.get(methods[m][1]).add(localRaceDir+filename+" simplifyQuotes "+Array.toStr(header, " "));
							} else {
								groupParams.get(methods[m][1]).add(localRaceDir+filename+" simplifyQuotes "+Array.toStr(Array.subArray(header, 1), " "));
							}
						}
					}
					log.report("", true, false);

					filename = races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0]+".csv";
					pvalFile = "meta_"+races[k][0]+"_pvals.dat";

					count = Files.parsePvals(localRaceDir+filename, localDir+pvalFile, "Total", methods[m], macHashesHashByRace.get(races[k][0]), macThresholdTotal, log);
					if (count == -1) {
						return;
					}
					lineCounts.add("Meta\t"+races[k][0]+"\t"+phenotypes[i][0]+"\t"+methods[m][0]+"\t"+count);
					
					filenames += pvalFile+",1=Meta_"+races[k][0]+";";

					header = Files.getHeaderOfFile(localRaceDir+filename, ",!", log);
					header[0] = "'"+header[0]+"'";
					for (int h = 1; h < header.length; h++) {
						header[h] = "'"+header[h]+"'="+races[k][0]+"_"+header[h]+"_"+methods[m][0];
					}
					if (methods[m][1].equals("SingleVariant")) {
						temp = header[0];
						header[0] = header[1];
						header[1] = temp;
						header = Array.subArray(header, 0, getHeaderForMethod(methods[m]).length);
					}
					groupParams.get(methods[m][1]).add(k, localRaceDir+filename+" simplifyQuotes "+Array.toStr(header, " "));
					groupHits.get(methods[m][1]).incorporateFromFile(localDir+pvalFile, new int[] {0,1}, 0.001, log);
				}

				filename = phenotypes[i][0]+"_"+methods[m][0]+".csv";
				pvalFile = "meta_panEthnic_pvals.dat";
				
				count = Files.parsePvals(localDir+filename, localDir+pvalFile, "Total", methods[m], macHashesHashByRace.get("PanEthnic"), macThresholdTotal, log);
				if (count == -1) {
					return;
				}
				lineCounts.add("Meta\t"+"PanEthnic"+"\t"+phenotypes[i][0]+"\t"+methods[m][0]+"\t"+count);
				
				filenames += pvalFile+",1=Meta_PanEthnic;";
				
				header = Files.getHeaderOfFile(localDir+filename, ",!", log);
				header[0] = "'"+header[0]+"'";
				for (int h = 1; h < header.length; h++) {
					header[h] = "'"+header[h]+"'=PanEthnic_"+header[h]+"_"+methods[m][0];
				}
				if (methods[m][1].equals("SingleVariant")) {
					temp = header[0];
					header[0] = header[1];
					header[1] = temp;
					header = Array.subArray(header, 0, getHeaderForMethod(methods[m]).length);
				}
				groupParams.get(methods[m][1]).add(0, localDir+filename+" simplifyQuotes "+Array.toStr(header, " "));
				groupHits.get(methods[m][1]).incorporateFromFile(localDir+pvalFile, new int[] {0,1}, 0.001, log);
				
				

//				filename = PHENOTYPES[i][0]+"_"+METHODS[j][0]+".pval.metal";
//				if (!Files.exists(localDir+filename+"1.out") || new File(localDir+filename+"1.out").length() < 500) {
//					Metal.metaAnalyze(localDir, Array.toStringArray(locals), UNIT_OF_ANALYSIS[j], filename, Metal.PVAL_ANALYSIS, null, log);
//					running = true;
//				} else {
//					groupHits.get(GROUPS[j]).incorporateFromFile(localDir+filename+"1.out", 0.001, log);
//					groupParams.get(GROUPS[j]).add(2, localDir+filename+"1.out 'MarkerName' 'Allele1'=Meta_pA1_"+METHODS[j][0]+" 'Allele2'=Meta_pA2_"+METHODS[j][0]+" 'P-value'=Meta_Nweighted_pval_"+METHODS[j][0]);
//					if (SINGLE_VARIANTS[j]) {
//						groupParams.get(GROUPS[j]).add(0, localDir+filename+"1.out 'MarkerName' 'Freq1'=Meta_AF");
//						filenames += localDir+filename+"1.out,7=pMeta;";
//					} else {
//						filenames += localDir+filename+"1.out,5=pMeta;";
//					}
//				}
				if (filenames.length() == 0) {
					System.err.println("Error - why are there no files for "+phenotypes[i][0]+" "+methods[m][0]);
				}
				Files.write("java -cp /home/npankrat/vis.jar cnv.plots.QQPlot files=\""+filenames.substring(0, filenames.length()-1)+"\" maxToPlot=10", localDir+"plotQQs.bat");
			}
			
			groups = HashVec.getKeys(groupHits);
			for (int g = 0; g < groups.length; g++) {
				filename = dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+"_hits.dat";
				groupHits.get(groups[g]).writeHits(filename);
				hits = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
				groupParams.get(groups[g]).add(0, filename+" 0 1=minPval skip=0");
				if (groups[g].equals("BurdenTests")) {
					line = HashVec.getKeys(macHashesHashByRace.get("PanEthnic"));
					for (int k = 0; k < line.length; k++) {
						groupParams.get(groups[g]).add(1, dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+line[k]+".xln 0 'Total'=MAC<"+line[k]+"%");
					}
				}
				
//				if (groups[g].equals("SingleVariant")) {
//					filename = dir+"CHARGE_wGC/"+phenotypes[i][0]+"/SingleSNP/"+phenotypes[i][0]+"_SingleSNP.se.metal1.out";
//					groupParams.get(groups[g]).add(1, filename+" 0 5=CHARGE_pval");
//					filename = dir+"ESP_"+phenotypes[i][0]+"_SingleSNP.csv";
//					groupParams.get(groups[g]).add(2, filename+" 0 10=ESP_pval");
//				} else {
//					filename = dir+"CHARGE_wGC/"+phenotypes[i][0]+"/T5Count/"+phenotypes[i][0]+"_T5Count.se.metal1.out";
//					groupParams.get(groups[g]).add(1, filename+" 0 5=CHARGE_pval");
//					filename = dir+"ESP."+phenotypes[i][0]+".T5.csv";
//					groupParams.get(groups[g]).add(2, filename+" 0 'pval'=ESP_pval");
//				}
				count = 0;
				for (int l = 0; l < groupAnnotationParams.length; l++) {
					if (groupAnnotationParams[l][0].equals(groups[g])) {
						groupParams.get(groups[g]).add(count, dir+groupAnnotationParams[l][1]);
						count++;
					}
				}
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+"_parser.crf"));
					writer.println("hits");
					writer.println(phenotypes[i][0]+"_hitters.dat 0 out="+dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+".csv");
					writer.println(Array.toStr(Array.toStringArray(groupParams.get(groups[g])), "\n"));
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + "parser_"+groups[g]+"_"+phenotypes[i][0]+".crf");
					e.printStackTrace();
				}
				Files.writeList(hits, phenotypes[i][0]+"_hitters.dat");
				Files.combine(hits, Array.toStringArray(groupParams.get(groups[g])), null, groups[g], ".", dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+".csv", log, true, true, false);
				System.out.println(Array.toStr(Array.toStringArray(groupParams.get(groups[g])), "\n"));
			}
		}
		
		copyHits(dir, maps);
		delineateRegions(dir, maps);
		
		lineCounts.insertElementAt("Study\tRace\tPhenotype\tMethod\tCount", 0);
		Files.writeList(Array.toStringArray(lineCounts), dir+"lineCounts.xln");
		log.report("check lineCounts.xln for completeness");
	}

	public static void copyHits(String dir, MetaAnalysisParams maps) {
		String[] groups;
		String filename;
		String[][] phenotypes;
		String[][] methods;
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		methods = maps.getMethods();

		new File(dir+"hitsAssembled/").mkdirs();
		for (int i = 0; i < phenotypes.length; i++) {
			groups = new String[] {};
			
			for (int m = 0; m < methods.length; m++) {
				if (ext.indexOfStr(methods[m][1], groups) == -1) {
					groups = Array.addStrToArray(methods[m][1], groups);
				}
			}
			
			for (int j = 0; j < groups.length; j++) {
				filename = dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[j]+".csv";
				System.out.println("cp "+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[j]+".csv hitsAssembled/");
				Files.copyFile(filename, dir+"hitsAssembled/"+ext.removeDirectoryInfo(filename));
			}
		}
	}
	
	public static void parseGenePositions(String dir, MetaAnalysisParams maps) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		Hashtable<String, Segment> hash;
		Hashtable<String, String> errors;
		long time;
		
		String filename;
		String[] header;
		int[] indices;
		Logger log;
		String delimiter;
		Segment seg, trav;
		String[] geneNames;
		byte[] chrs;
		int[] positions;
		CountHash ch;
		
		time = new Date().getTime();
		System.out.println(ext.getTime(time)+" started parsing gene positions");
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir+"parseGenePositions.log");
		
		filename = maps.getSnpInfoFilename();
		filename = ext.rootOf(filename)+".csv";
		System.out.println("Attempting to parse "+dir+filename);

		ch = new CountHash();
		hash = new Hashtable<String, Segment>();
		errors = new Hashtable<String, String>();
		try {
			reader = Files.getAppropriateReader(filename);
			temp = reader.readLine();
			delimiter = ext.determineDelimiter(temp);
			if (delimiter.equals(",")) {
				header = ext.splitCommasIntelligently(temp, true, log);
			} else {
				header = temp.trim().split(delimiter);
			}
			indices = ext.indexFactors(new String[][] {Aliases.GENE_UNITS,  Aliases.CHRS,  Aliases.POSITIONS}, header, false, false, true, false, log, true);
			while (reader.ready()) {
				temp = reader.readLine();
				if (delimiter.equals(",")) {
					line = ext.splitCommasIntelligently(temp, true, log);
				} else {
					line = temp.trim().split(delimiter);
				}
				line = Array.subArray(line, indices);
				trav = new Segment(Positions.chromosomeNumber(line[1], log), Integer.parseInt(line[2]), Integer.parseInt(line[2]));
				if (hash.containsKey(line[0])) {
					seg = hash.get(line[0]);
					if (trav.getChr() != seg.getChr()) {
						if (!errors.containsKey(line[0]+"\t"+seg.getChr()+" and "+trav.getChr())) {
							System.err.println("Error - still have gene '"+line[0]+"' listed on two different chromosomes ("+seg.getChr()+" and "+trav.getChr()+")");
							errors.put(line[0]+"\t"+seg.getChr()+" and "+trav.getChr(), "");
						}
					} else if (trav.getStart() < seg.getStart()) {
						seg = new Segment(seg.getChr(), trav.getStart(), seg.getStop());
					} else if (trav.getStop() > seg.getStop()) {
						seg = new Segment(seg.getChr(), seg.getStart(), trav.getStop());
					}
				} else {
					seg = trav;
				}
				hash.put(line[0], seg);
				ch.add(line[0]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		System.out.println("Read in all data in " + ext.getTimeElapsed(time));
		
		geneNames = HashVec.getKeys(hash, false, false);
		chrs = new byte[geneNames.length];
		positions = new int[geneNames.length];
		for (int i = 0; i < geneNames.length; i++) {
			seg = hash.get(geneNames[i]);
			chrs[i] = seg.getChr();
			positions[i] = seg.getStart();
		}
		
		System.out.println("Coalated positions at " + ext.getTimeElapsed(time));

		indices = Sort.orderTwoLayers(chrs, positions);
		
		System.out.println("Finished sorting at " + ext.getTimeElapsed(time));
			
		try {
			writer = new PrintWriter(new FileWriter(dir+"gene_positions.xln"));
			writer.println("Gene\tChr\tPositionOfFirstMarkerInGene\tPositionOfLastMarkerInGene\tNumMarkersInGene");
			for (int i = 0; i < geneNames.length; i++) {
				seg = hash.get(geneNames[indices[i]]);
				writer.println(geneNames[indices[i]]+"\t"+seg.getChr()+"\t"+seg.getStart()+"\t"+seg.getStop()+"\t"+ch.getCount(geneNames[indices[i]]));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"gene_positions.xln");
			e.printStackTrace();
		}

		System.out.println("and done with a total time of " + ext.getTimeElapsed(time));
	}
	
	public static void delineateRegions(String dir, MetaAnalysisParams maps) {
		Vector<Vector<String>> filesToCat;
		int[] ns;
		float indexThreshold;
		Logger log, pvalThresholdsLog;
		int index;
		CountHash countHash;
		
		String[] groups, races;
		String filename;
		String[][] phenotypes, methods, results;
		Vector<String> additionalCols;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		
		log = new Logger(dir+"delineateRegions.log");
		
		if (!Files.exists(dir+"hitsAssembled/")) {
			copyHits(dir, maps);
		}
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		methods = maps.getMethods();
		races = Matrix.extractColumn(maps.getRacesWithFilenameAliases(), 0);

		groups = new String[] {};
		filesToCat = new Vector<Vector<String>>();
		countHash = new CountHash();
		for (int m = 0; m < methods.length; m++) {
			if (ext.indexOfStr(methods[m][1], groups) == -1) {
				groups = Array.addStrToArray(methods[m][1], groups);
				filesToCat.add(new Vector<String>());
			}
			countHash.add(methods[m][1]);
		}

		pvalThresholdsLog = new Logger(dir+"pval_thresholds.xln");
		pvalThresholdsLog.report("Trait\tN\tnTests\tBonferroni_pval");
		for (int i = 0; i < phenotypes.length; i++) {
			ns = Array.intArray(groups.length, -1);
			for (int m = 0; m < methods.length; m++) {
				filename = phenotypes[i][0]+"/"+methods[m][0]+"/meta_panEthnic_pvals.dat";
				index = ext.indexOfStr(methods[m][1], groups);
				if (Files.exists(dir+filename)) {
					ns[index] = Math.max(ns[index], Files.countLines(filename, true));
				}
			}
			for (int g = 0; g < groups.length; g++) {
				if (groups[g].equals("SingleVariant") && ns[g] == -1) {
					ns[g] = 1000000;
				} else if (groups[g].equals("BurdenTests") && ns[g] == -1) {
					ns[g] = 20000;
				}
				
				indexThreshold = (float)(0.05 / (double)ns[g] / (double)countHash.getCount(groups[g]));
//				indexThreshold = (float)(0.05 / (double)ns[g]);
				filename = phenotypes[i][0]+"_"+groups[g]+".csv";
				pvalThresholdsLog.report(ext.rootOf(filename)+"\t"+ns[g]+"\t"+countHash.getCount(groups[g])+"\t"+indexThreshold);
				additionalCols = new Vector<String>();
				if (groups[g].equals("SingleVariant")) {
					additionalCols.add("SKATgene");
					additionalCols.add("CHARGE_ALL_AF");
//					additionalCols.add("function");
					additionalCols.add("single_func_region");
				}
				if (groups[g].equals("BurdenTests")) {
					additionalCols.add("PanEthnic_nsnpsTotal_T5Count");
//					additionalCols.add("PanEthnic_Meta_nsnpsTotal_T5Count");
				}
				for (int j = 0; j < methods.length; j++) {
					if (methods[j][1].equals(groups[g])) {
						additionalCols.add("PanEthnic_p_"+methods[j][0]);
//						additionalCols.add("PanEthnic_Meta_p_"+methods[j][0]);
						for (int k = 0; k < races.length; k++) {
							additionalCols.add(races[k]+"_p_"+methods[j][0]);
//							additionalCols.add(races[k]+"_Meta_p_"+methods[j][0]);
						}
					}
				}
				if (Files.exists(dir+"hitsAssembled/"+filename)) {
					results = HitWindows.determine(dir+"hitsAssembled/"+filename, indexThreshold, 500000, indexThreshold*100, Array.toStringArray(additionalCols));
					Files.writeMatrix(results, ext.rootOf(dir+"hitsAssembled/"+filename, false)+"_regions.xln", "\t");
					String temp = phenotypes[i][0];
					Files.write(temp, dir+"hitsAssembled/"+temp+".tmp");
					filesToCat.elementAt(g).add(dir+"hitsAssembled/"+temp+".tmp");
					filesToCat.elementAt(g).add(ext.rootOf(dir+"hitsAssembled/"+filename, false)+"_regions.xln");
				} else {
					log.reportError("Error - could not find expected file: "+dir+"hitsAssembled/"+filename);
				}
			}
		}
		
		for (int g = 0; g < groups.length; g++) {
			Files.cat(Array.toStringArray(filesToCat.elementAt(g)), dir+"hitsAssembled/"+groups[g]+"_regions.xln", null, log);
		}
	}
	
	public static void makeQQplots(String dir, MetaAnalysisParams maps) {
//		String[] groups;
//		String filename;
		String[][] phenotypes;
		String[][] methods;
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		methods = maps.getMethods();

		new File(dir+"hitsAssembled/").mkdirs();
		for (int i = 0; i < phenotypes.length; i++) {
//			groups = new String[] {};
			
			for (int m = 0; m < methods.length; m++) {
//				cnv.plots.QQPlot
			}
		}
	}
	
	public static void stitch(String dir, String pattern, String fileout, Logger log) {
		String[] list;
		int[] skips;
		int maxChr;

		maxChr = getMaxChr();
		list = new String[maxChr];
		for (int chr = 1; chr <= maxChr; chr++) {
			list[chr-1] = dir+ext.replaceAllWith(pattern, "#", chr==23?"X":(chr==24?"Y":chr+""));
		}
		skips = Array.intArray(list.length, 1);
		skips[0] = 0;
		Files.cat(list, dir+fileout, skips, log);
	}
	
	public static void removeIndicesFromRdata(String filein, String fileout) {
		BufferedReader reader;
		PrintWriter writer;
		String temp;
		
		try {
			reader = new BufferedReader(new FileReader(filein));
			writer = new PrintWriter(new FileWriter(fileout));
			writer.println(reader.readLine());
			while (reader.ready()) {
				temp = reader.readLine();
				writer.println(temp.substring(temp.indexOf(",")+1));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filein + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filein + "\"");
			System.exit(2);
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String logfile = null;
		Logger log;
		String dir = "";
		boolean determineObjectNames = false;
		boolean splitAll = false;
		boolean runAll = false;
		boolean parseAll = false;
		boolean checkNs = false;
		boolean computeMACs = false;
		boolean hits = false;
		boolean regions = false;
		String mafThreshold = null;
		int macThresholdStudy = 5;
		int macThresholdTotal = 40;
		String mapsFile = "metaAnalysis.params";
		MetaAnalysisParams maps;
		boolean consolidate = false;
		boolean metrics = false;
		boolean forceMeta = false;
		boolean copy = false;
		boolean qq = false;
		boolean genePositions = false;

		String usage = "\n" + 
		"gwas.SeqMeta requires 0-1 arguments\n" + 
		"   (0) directory (i.e. dir=" + dir + " (default))\n" + 
		"   (1) filename of MetaAnalysisParameters (i.e. maps=" + mapsFile + " (default; create an empty file of this name to populate with examples))\n" + 
		" AND\n" + 
		"   (2) determine object names (i.e. -determineObjectNames (not the default))\n" + 
		" OR\n" + 
		"   (2) split all (i.e. -splitAll (not the default))\n" + 
		" OR\n" + 
		"   (2) consolidate split files (i.e. -consolidate (not the default))\n" + 
		" OR\n" + 
		"   (2) run all (i.e. -runAll (not the default))\n" +
		"   (4) force the meta-analysis to be redone even if meta-analysis output files exist (i.e. -forceMeta (not the default))\n" +
		" OR\n" + 
		"   (2) parse all runs (i.e. -parseAll (not the default))\n" + 
		"   (3) if the meta-analysis was forced, then this will restitch the files (i.e. -forceMeta (not the default))\n" +
		" OR\n" + 
		"   (2) get summary metrics for single SNP results (i.e. -metrics (not the default))\n" + 
		" OR\n" + 
		"   (2) check sample sizes for each chr18 data file (i.e. -checkNs (not the default))\n" + 
		" OR\n" + 
		"   (2) compute minor allele counts (i.e. -computeMACs (not the default))\n" + 
		"   (3) minor allele frequency threshold (i.e. mafThreshold=0.01 (not the default; if left out will compute for all relevant threhsolds))\n" + 
		" OR\n" + 
		"   (2) delineate start and stop positions for the genes (i.e. -genePositions (not the default))\n" + 
		" OR\n" + 
		"   (2) parse top hits (i.e. -hits (not the default))\n" + 
		"   (3) minor allele count threshold for a study (i.e. macThresholdStudy="+macThresholdStudy+" (default))\n" + 
		"   (4) minor allele count threshold for meta-analysis (i.e. macThresholdTotal="+macThresholdTotal+" (default))\n" + 
		" OR\n" + 
		"   (2) parse regions from top hits (i.e. -regions (not the default))\n" + 
		" OR\n" + 
		"   (2) copy top hit files (i.e. -copy (not the default))\n" + 
		" OR\n" + 
		"   (2) makeQ plots (i.e. -qq (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-determineObjectNames")) {
				determineObjectNames = true;
				numArgs--;
			} else if (args[i].startsWith("-splitAll")) {
				splitAll = true;
				numArgs--;
			} else if (args[i].startsWith("-consolidate")) {
				consolidate = true;
				numArgs--;
			} else if (args[i].startsWith("-runAll")) {
				runAll = true;
				numArgs--;
			} else if (args[i].startsWith("-forceMeta")) {
				forceMeta = true;
				numArgs--;
			} else if (args[i].startsWith("-parseAll")) {
				parseAll = true;
				numArgs--;
			} else if (args[i].startsWith("-metrics")) {
				metrics = true;
				numArgs--;
			} else if (args[i].startsWith("-checkNs")) {
				checkNs = true;
				numArgs--;
			} else if (args[i].startsWith("-computeMACs")) {
				computeMACs = true;
				numArgs--;
			} else if (args[i].startsWith("-genePositions")) {
				genePositions = true;
				numArgs--;
			} else if (args[i].startsWith("mafThreshold=")) {
				mafThreshold = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("macThresholdStudy=")) {
				macThresholdStudy = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("macThresholdTotal=")) {
				macThresholdTotal = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-hits")) {
				hits = true;
				numArgs--;
			} else if (args[i].startsWith("-regions")) {
				regions = true;
				numArgs--;
			} else if (args[i].startsWith("-copy")) {
				copy = true;
				numArgs--;
			} else if (args[i].startsWith("-qq")) {
				qq = true;
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

//		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/metaAnalysis/";
//		stitch(dir, "CHS_Fibrinogen_SingleSNP_chr#.csv");
//		removeIndicesFromRdata(dir+"CHS_Fibrinogen_SingleSNP.csv", dir+"CHS_Fibrinogen_SingleSNP_slim.csv");

//		stitch(dir, "CHS_Fibrinogen_SKAT_T5_chr#.csv");
//		removeIndicesFromRdata(dir+"CHS_Fibrinogen_SKAT_T5.csv", dir+"CHS_Fibrinogen_SKAT_T5_slim.csv");
		
//		System.exit(1);
		
//		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/metaAnalysis/";
//		splitAll = true;
//		runAll = true;
//		computeMAC(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, "0.01");
//		computeMAC(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, "0.05");
//		System.exit(1);
//		hits = true;
		
//		parseGenePositions("", new MetaAnalysisParams("metaAnalysis.params", new Logger()));
//		System.exit(1);

//		dir = "D:/ExomeChip/Hematology/00src/CHARGE-RBC/06_everythingIncludingYFS/";
//		delineateRegions(dir, new MetaAnalysisParams(dir+"metaAnalysis.params", new Logger()), 1E-7, 2E-6);
//		System.exit(1);
		
//		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/metaAnalysis2/";
//		delineateRegions(dir, new MetaAnalysisParams(dir+"metaAnalysis.params", new Logger()), 1E-7, 2E-6);
//		System.exit(1);
		
//		dir = "D:/ExomeChip/Hematology/00src/CHARGE-RBC/06_everythingIncludingYFS/";
//		delineateRegions(dir, new MetaAnalysisParams(dir+"metaAnalysis.params", new Logger()), 1E-7, 2E-6);
//		System.exit(1);

//		dir = "D:/ExomeChip/Hematology/results/07_outliersRemoved/allBurdens/";
//		delineateRegions(dir, new MetaAnalysisParams(dir+"metaAnalysis.params", new Logger()));
//		System.exit(1);
		
		try {
			log = new Logger(logfile);
			maps = new MetaAnalysisParams(mapsFile, log);
			if (determineObjectNames) {
				determineObjectNames(dir, maps, log);
			} else if (splitAll) {
				splitAll(dir, maps);
			} else if (consolidate) {
				consolidate(dir, maps);
			} else if (runAll) {
				runAll(dir, maps, forceMeta);
			} else if (parseAll) {
				parseAll(dir, maps, forceMeta);
			} else if (checkNs) {
				doubleCheckNs(dir, maps);
			} else if (metrics) {
				parseMetrics(dir, maps);
			} else if (computeMACs) {
				if (mafThreshold == null) {
					computeAllRelevantMACs(dir, maps, log);
				} else {
					computeMAC(dir, maps, mafThreshold, new Logger(dir+"computeMACs_forMAF_LTE_"+mafThreshold+".log"));
				}
			} else if (hits) {
				assembleHits(dir, maps, macThresholdStudy, macThresholdTotal);
//				metaAll(dir, PHENOTYPES, STUDIES, GROUPS, METHODS, UNIT_OF_ANALYSIS, DEFAULT_SAMPLE_SIZES, WEIGHTED, SINGLE_VARIANTS, GROUP_ANNOTATION_PARAMS);
			} else if (genePositions) {
				parseGenePositions(dir, maps);
			} else if (regions) {
				delineateRegions(dir, maps);
			} else if (copy) {
				copyHits(dir, maps);
			} else if (qq) {
				makeQQplots(dir, maps);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
