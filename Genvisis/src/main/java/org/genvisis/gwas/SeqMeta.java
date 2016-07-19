package org.genvisis.gwas;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.filesys.*;
import org.genvisis.mining.Transformations;
import org.genvisis.parse.GenParser;
import org.genvisis.parse.LookupTable;
import org.genvisis.stats.Correlation;
import org.genvisis.stats.ProbDist;
import org.genvisis.stats.Rscript;

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
	
	public static final String[] SUMMARY_INFO_HEADER = {"Study", "Ethnicity", "Units", "Trait", "meanTrait", "stdevTrait", "minTrait", "maxTrait", "numFemales", "numMales", "meanAge", "stdevAge", "minAge", "maxAge"};
	
	public static String getRscriptExecutable(MetaAnalysisParams maps, Logger log) {
		if (maps != null && maps.getRExec() != null) {
			return maps.getRExec();
		} else {
			return Rscript.getRscriptExecutable(log);
		}
	}
	
	public static void determineObjectNames(String dir, MetaAnalysisParams maps, Logger log) {
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
		log.report("There are "+v.size()+" .Rdata files remaining to interrogate:\n"+Array.toStr(Array.toStringArray(v), "\n")+"\n\n./master.checkObjectAll");

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
		boolean[] picks;
		int numMatches, index;
		String[][] phenotypes, maskedPhenos;
		String[] studies;
		String[][] races;
		String[] usedFor;
		boolean first;
		
		usedFor = Array.stringArray(files.length, null);
		
		index = ext.indexOfStr(maps.getSnpInfoFilename(), files);
		if (index >= 0) {
			usedFor[index] = "SnpInfo file";
		}
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		maskedPhenos = maps.getPhenotypesWithFilenameAliases(false);
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();

		finalSets = new String[phenotypes.length][][]; // [pheno][study][race] <- all files meeting criteria
		for (int i = 0; i < phenotypes.length; i++) {
			finalSets[i] = Matrix.stringMatrix(studies.length, races.length, "<missing>");
			if (log.getLevel() > 5) {
				log.report("For "+phenotypes[i][0]+" identified:", true, true);
				log.report("\tStudy\t"+Array.toStr(Matrix.extractColumn(races, 0)));
			}
			for (int j = 0; j < studies.length; j++) {
				log.report("\t"+studies[j], false, true, 5);
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
							if (usedFor[f] != null) {
								log.reportError("Error - file '"+files[f]+"' matches to "+studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0]+" but was already picked for "+usedFor[f]);
							}
							usedFor[f] = studies[j]+"/"+races[k][0]+"/"+phenotypes[i][0];
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
					log.report("\t"+numMatches, false, true, 5);
				}
				log.report("", true, true, 5);
			}
			log.report("", true, true, 5);
		}

		first = true;
		for (int i = 0; i < files.length; i++) {
			if (usedFor[i] == null) {
				if (first) {
					log.reportError("Warning - did not find a match for the following file(s):");
					first = false;
				}
				log.reportError("  "+files[i]);
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
		Vector<String> commands;
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
		
		jobNames = new Vector<String>();
		jobSizes = new IntVector();
		
		dir = ext.verifyDirFormat(dir);
		finalSets = identifySet(maps, files, log);
		for (int i = 0; i < phenotypes.length; i++) {
//			log.report("Phenotype: "+phenotypes[i][0]);
			for (int j = 0; j < studies.length; j++) {
//				log.report("  Study: "+studies[j]);
				for (int k = 0; k < races.length; k++) {
//					log.report("    Race: "+races[k][0]);
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
	
								jobNames.add(getRscriptExecutable(maps, log)+" --no-save "+filename);
								jobSizes.add((int)(new File(dir+files[f]).length()));
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
		
		new File("chunks/").mkdir();
		if (jobNames.size() > 0) { 
			Files.qsubExecutor(dir, jobNames, jobSizes, "chunks/chunkSplit", 24, 62000, 2);
		} else {
			log.report("\nLooks like everything has been split");
		}
		
		log.report("");
		log.report("Make sure to run \"qsub splitChrs.qsub\" first!!!");
		log.report("Then run \"qsub chunks/chunkSplit.pbs\"");
		log.report("Make sure to run \"seq -consolidate\" after everything else is run!!!");
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
		Vector<String> commands, objects;
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
							commands.add("print(.libPaths())");
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
										if (chr < 23) {
											log.report("Creating a dummy file for "+outputFilename+" because "+objectFilename+" has a filesize of "+new File(objectFilename).length());
										}
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
			
								jobNames.add(getRscriptExecutable(maps, log)+" --no-save "+filename);
								jobSizes.add(infoSizes[chr]);
							}
						}
					}
				}
			}
		}

		if (jobNames.size() == 0) {
			log.report("All cohorts have been run for all phenotypes");
		} else if (jobNames.size() == 1) {
			log.report("There is one remaining individual cohort analysis yet to be run using:   qsub chunkRun.pbs");
		} else {
			log.report("There are "+jobNames.size()+" individual cohort analyses yet to be run using:   qsub chunkRun.pbs");
		}
		new File("chunks/").mkdir();
		Files.qsubExecutor(dir, jobNames, jobSizes, "chunks/chunkRun", 24, 62000, 24);
		
		jobNames = new Vector<String>();
		jobSizes = new IntVector();
		for (int i = 0; i < phenotypes.length; i++) {
			// Meta-analysis stratified by race
			for (int k = 0; k < races.length; k++) {
				for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
					chrom = chr==23?"X":(chr==24?"Y":chr+"");
					commands = new Vector<String>();
					commands.add("print(.libPaths())");
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
			
						jobNames.add(getRscriptExecutable(maps, log)+" --no-save "+filename);
						jobSizes.add(infoSizes[chr]);
					}
				}
			}

			// Meta-analysis of all races
			for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
				chrom = chr==23?"X":(chr==24?"Y":chr+"");
				commands = new Vector<String>();
				commands.add("print(.libPaths())");
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
		
					jobNames.add(getRscriptExecutable(maps, log)+" --no-save "+filename);
					jobSizes.add(infoSizes[chr]);
				}
			}
		}
		
		if (jobNames.size() == 0) {
			log.report("All meta-analyses have been run for all phenotypes");
		} else if (jobNames.size() == 1) {
			log.report("There is one meta-analysis yet to be run using:   qsub chunkMeta.pbs");
		} else {
			log.report("There are "+jobNames.size()+" meta-analyses yet to be run using:   qsub chunkMeta.pbs");
		}
		Files.qsubExecutor(dir, jobNames, jobSizes, "chunks/chunkMeta", 24, 62000, 24);
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
								log.report(ext.getTime()+"\tStiching up "+root+".csv");
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
						log.report(ext.getTime()+"\tStiching up "+root+".csv");
						stitch(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/", root+"_chr#.csv", root+".csv", log);
					}
				}
			}
			
			for (int m = 0; m < methods.length; m++) {
				root = phenotypes[i][0]+"_"+methods[m][0];
				if (forceMeta || !Files.exists(dir+phenotypes[i][0]+"/"+methods[m][0]+"/"+root+".csv") || new File(dir+phenotypes[i][0]+"/"+methods[m][0]+"/"+root+".csv").length() == 0) {
					log.report(ext.getTime()+"\tStiching up "+root+".csv");
					stitch(dir+phenotypes[i][0]+"/"+methods[m][0]+"/", root+"_chr#.csv", root+".csv", log);
				}
			}
		}
	}

	public static void parseMetrics(String dir, MetaAnalysisParams maps) {
		Hashtable<String, String> hash;
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
		Vector<int[]> ranges;
		int start, stop;
		int count;
		String temp;
		
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

		hash = new Hashtable<String, String>();
		try {
			writer = new PrintWriter(new FileWriter(dir+"summary_metrics_list.xln"));
//			writer.println("Study\tRace\tPhenotype\tN_samples\tn_genotyped\tn_MAF1%\tLambda_MAF1%\tbetaMean_MAF1%\tbetaSD_MAF1%\tn_20count\tLambda_20count\tbetaMean_20count\tbetaSD_20count\tn_MAF5%\tLambda_MAF5%\tbetaMean_MAF5%\tbetaSD_MAF5%");
			writer.println("Study\tRace\tPhenotype\tN_samples\tn_genotyped\tn_MAF1%\tLambda_MAF1%\tbetaMedian_MAF1%\tbetaSD_MAF1%\tn_20count\tLambda_20count\tbetaMedian_20count\tbetaSD_20count\tn_MAF5%\tLambda_MAF5%\tbetaMedian_MAF5%\tbetaSD_MAF5%");
			for (int j = 0; j < studies.length; j++) {
				for (int k = 0; k < races.length; k++) {
					for (int i = 0; i < phenotypes.length; i++) {
						if (finalSets[i][j][k].equals("<missing>")) {
							sampleSizes[k][j][i] = "";
						} else {
							root = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+method;
							outputFilename = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+method+"/"+root+".csv";
							if (Files.exists(outputFilename)) {
								log.report("Processing "+outputFilename);
								line = procFile(outputFilename, log);
								sampleSizes[k][j][i] = line[0];
								writer.println(studies[j]+"\t"+races[k][0]+"\t"+phenotypes[i][0]+"\t"+Array.toStr(line));
								hash.put(studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0], line[7]);
							} else {
								log.report("Error - missing expected file: "+outputFilename);
							}
						}
					}
				}
			}

			for (int i = 0; i < phenotypes.length; i++) {
				for (int k = 0; k <= races.length; k++) {
					root = (k == races.length?"":races[k][0]+"_")+phenotypes[i][0]+"_"+method;
					outputFilename = dir+phenotypes[i][0]+"/"+(k == races.length?"":races[k][0]+"/")+method+"/"+root+".csv";
					if (Files.exists(outputFilename)) {
						log.report("Processing "+outputFilename);
						line = procFile(outputFilename, log);
						writer.println("Meta\t"+(k == races.length?"PanEthnic":races[k][0])+"\t"+phenotypes[i][0]+"\t"+Array.toStr(line));
						hash.put("Meta_"+(k == races.length?"PanEthnic":races[k][0])+"_"+phenotypes[i][0], line[7]);
					} else {
						log.report("Error - missing expected file: "+outputFilename);
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"summary_metrics_list.xln");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(dir+"summary_metrics_table.xln"));
			for (int k = 0; k <= races.length; k++) {
				writer.print(k==races.length?"PanEthnic":races[k][0]);
				for (int i = 0; i < phenotypes.length; i++) {
					writer.print("\t"+phenotypes[i][0]);
				}
				writer.println();
				for (int j = (k==races.length?studies.length:0); j <= studies.length; j++) {
					count = 0;
					temp = (j==studies.length?"Meta":studies[j]);
					for (int i = 0; i < phenotypes.length; i++) {
						if (k<races.length && j<studies.length && finalSets[i][j][k].equals("<missing>")) {
							temp += "\t";
						} else {
							root = (j==studies.length?"Meta":studies[j])+"_"+(k==races.length?"PanEthnic":races[k][0])+"_"+phenotypes[i][0];
							if (hash.containsKey(root)) {
								temp += "\t"+hash.get(root);
							} else {
								temp += "\t?!?!?!?!?!";
								log.report(root);
							}
							count++;
						}
					}
					if (count > 0) {
						writer.println(temp);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"summary_metrics_table.xln");
			e.printStackTrace();
		}
		
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
		BufferedReader reader;
		String[] header, line;
		String temp;
		int count;
		String delimiter;
		DoubleVector[] dvs;
		String[] metrics;
		int maxSamples;
		int pvalIndex, mafIndex, ntotalIndex, betaIndex;
		boolean[] keeps;
		boolean problem;
		
		problem = false;
		count = 0;
		maxSamples = 0;
		metrics = Array.stringArray(14, "");
		dvs = DoubleVector.newDoubleVectors(6);
		
		try {
			reader = new BufferedReader(new FileReader(outputFilename));
			temp = ext.replaceAllWith(reader.readLine(), "\"", "");
			delimiter = ext.determineDelimiter(temp);
			if (delimiter.equals(",")) {
				header = ext.splitCommasIntelligently(temp, true, log);
			} else {
				header = temp.split(delimiter, -1);
			}
			
			keeps = Array.booleanArray(header.length, true);
			if (header[4].equals("caf")) {
				keeps[4] = false;
			}

			if (ext.checkHeader(Array.subArray(header, keeps), HEADER_TYPES[0], Array.intArray(HEADER_TYPES[0].length), false, log, false)) {
				pvalIndex = ext.indexOfStr("p", header, false, true);
				mafIndex = ext.indexOfStr("maf", header, false, true);
				ntotalIndex = ext.indexOfStr("ntotal", header, false, true);
				betaIndex = ext.indexOfStr("beta", header, false, true);
				log.report(outputFilename);
				log.report(Array.toStr(header,"/"));
				log.report("p="+pvalIndex+"="+header[pvalIndex]+" maf="+mafIndex+"="+header[mafIndex]+" ntotal="+ntotalIndex+"="+header[ntotalIndex]+" beta="+betaIndex+"="+header[betaIndex]);
				while (reader.ready()) {
					temp = reader.readLine().trim();
					if (delimiter.equals(",")) {
						line = ext.splitCommasIntelligently(temp, true, log);
					} else {
						line = temp.split(delimiter, -1);
					}
					try {
						if (!line[ntotalIndex].equals("0") && !line[betaIndex].equals("NA")) {
							count++;
							if (Integer.parseInt(line[ntotalIndex]) > maxSamples) {
								maxSamples = Integer.parseInt(line[ntotalIndex]);
							}
							if (Double.parseDouble(line[mafIndex]) >= 0.01) {
								dvs[0].add(Double.parseDouble(line[pvalIndex]));
								dvs[1].add(Math.abs(Double.parseDouble(line[betaIndex])));
							}
							if (Double.parseDouble(line[mafIndex])*Double.parseDouble(line[ntotalIndex]) >= 10) {
								dvs[2].add(Double.parseDouble(line[pvalIndex]));
								dvs[3].add(Math.abs(Double.parseDouble(line[betaIndex])));
							}
							if (Double.parseDouble(line[mafIndex]) >= 0.05) {
								dvs[4].add(Double.parseDouble(line[pvalIndex]));
								dvs[5].add(Math.abs(Double.parseDouble(line[betaIndex])));
							}
						}
					} catch (Exception e) {
						if (!problem) {
							problem = true;
							log.reportError("Problem with file "+outputFilename);
						}
						log.reportError("Error parsing line: "+temp+" into "+Array.toStr(line, delimiter)+"   line[ntotalIndex/"+ntotalIndex+"]="+line[ntotalIndex]+"  line[betaIndex/"+betaIndex+"]="+line[betaIndex]);
					}
				}
				metrics[0] = maxSamples+"";
				metrics[1] = count+"";
				metrics[2] = dvs[0].size()+"";
				metrics[3] = ext.formDeci(ProbDist.ChiDistReverse(Array.median(dvs[0].toArray()), 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
				metrics[4] = ext.formDeci(Array.median(dvs[1].toArray()), 6);
				metrics[5] = ext.formDeci(Array.stdev(dvs[1].toArray()), 6);
				metrics[6] = dvs[2].size()+"";
				metrics[7] = ext.formDeci(ProbDist.ChiDistReverse(Array.median(dvs[2].toArray()), 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
				metrics[8] = ext.formDeci(Array.median(dvs[3].toArray()), 6);
				metrics[9] = ext.formDeci(Array.stdev(dvs[3].toArray()), 6);
				metrics[10] = dvs[4].size()+"";
				metrics[11] = ext.formDeci(ProbDist.ChiDistReverse(Array.median(dvs[4].toArray()), 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
				metrics[12] = ext.formDeci(Array.median(dvs[5].toArray()), 6);
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
			if (methods[m].length > 3 && ext.isValidDouble(methods[m][3])) {
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
		Hashtable<String, String> snpGeneFunctionalHash, snpGeneFunctionalHashPan, snpGeneFunctionalHashRaceSpecific;
		Hashtable<String, Vector<String>> geneLoci;
		String[] keys;
		Hashtable<String, int[]> macs, raceSpecificMacs;
		String gene;
		int[] counts;
		double mafThresholdDouble;
		String[][] phenotypes, races, methods;
		String[] studies;
		String snpInfoFile, functionFlagName;
		String[][] needs;
		
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
			log.report(ext.getTime()+"\tLoading "+filename+".mappings.ser");
			snpGeneHash = SerialHash.loadSerializedStringHash(filename+".mappings.ser");
			log.report(ext.getTime()+"\tReloaded marker mappings (n="+ext.addCommas(snpGeneHash.size())+" variants) in " + ext.getTimeElapsed(time));
			log.report(ext.getTime()+"\tLoading "+filename+".functionalMappings.ser");
			snpGeneFunctionalHash = SerialHash.loadSerializedStringHash(filename+".maf"+mafThreshold+".functionalMappings.ser");
			log.report(ext.getTime()+"\tReloaded marker mappings (n="+ext.addCommas(snpGeneFunctionalHash.size())+" functional variants) in " + ext.getTimeElapsed(time));
		} else {
			snpGeneHash = new Hashtable<String, String>();
			snpGeneFunctionalHash = new Hashtable<String, String>();
			geneLoci = new Hashtable<String, Vector<String>>();

			try {
				log.report(ext.getTime()+"\tReading in " + filename);
				reader = new BufferedReader(new FileReader(filename));
				header = ext.splitCommasIntelligently(reader.readLine(), true, log);
				needs = new String[][] {Aliases.MARKER_NAMES, Aliases.GENE_UNITS, new String[] {functionFlagName}, Aliases.CHRS};
				indices = ext.indexFactors(needs, header, false, true, true, log, false);
				if (Array.min(indices) == -1) {
					log.reportError("Improper header for file '"+filename+"', found: "+Array.toStr(header, "/")+"\nMissing one of these: "+Array.toStr(needs[ext.indexOfInt(-1, indices)], "/"));
					reader.close();
					return;
				}
				while (reader.ready()) {
					line = ext.splitCommasIntelligently(reader.readLine(), true, log);
					snpGeneHash.put(line[indices[0]], line[indices[1]]);
					if (line[indices[2]].equals("TRUE")) { // && !line[indices[4]].equals("NA") && Double.parseDouble(line[indices[4]]) <= mafThresholdDouble ) {
						snpGeneFunctionalHash.put(line[indices[0]], line[indices[1]]);
					}
					if (!geneLoci.containsKey(line[indices[1]])) {
						geneLoci.put(line[indices[1]], new Vector<String>());
					}
					HashVec.addIfAbsent(line[indices[3]], geneLoci.get(line[indices[1]]));
				}
				reader.close();
				log.report(ext.getTime()+"\tProcessed marker mappings (n="+ext.addCommas(snpGeneHash.size())+" variants; n="+ext.addCommas(snpGeneFunctionalHash.size())+" functional variants) in " + ext.getTimeElapsed(time));
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + filename + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + filename + "\"");
				System.exit(2);
			}
			
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
		
		
		if (!methods[0][0].startsWith("SingleSNP")) {
			System.err.println("Error - this program erroneously assumed that the first model was SingleSNP and got confused (it's actually "+methods[0][0]+"); aborting");
			return;
		}
		
		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(maps, files, log);

		for (int i = 0; i < phenotypes.length; i++) {
			log.report(ext.getTime()+"\tStarting calculations for "+phenotypes[i][0]);
			macs = new Hashtable<String, int[]>();
			snpGeneFunctionalHashPan = filterSnpGeneFunctionalHash(dir+phenotypes[i][0]+"/"+methods[0][0]+"/"+phenotypes[i][0]+"_"+methods[0][0]+".csv", snpGeneFunctionalHash, mafThresholdDouble, log);
			log.report(ext.getTime()+"\tThere are "+ext.addCommas(snpGeneFunctionalHashPan.size())+" functional variants remaining after MAF checks in the cross-ethnic meta-analysis");

			for (int k = 0; k < races.length; k++) {
				log.report(ext.getTime()+"\tStarting calculations for "+phenotypes[i][0]+" specifically for "+races[k][0]);
				raceSpecificMacs = new Hashtable<String, int[]>();
				snpGeneFunctionalHashRaceSpecific = filterSnpGeneFunctionalHash(dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[0][0]+"/"+races[k][0]+"_"+phenotypes[i][0]+"_"+methods[0][0]+".csv", snpGeneFunctionalHash, mafThresholdDouble, log);
				log.report(ext.getTime()+"\tThere are "+ext.addCommas(snpGeneFunctionalHashRaceSpecific.size())+" functional variants remaining after MAF checks in the "+races[k][0]+" meta-analysis");
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
								if (snpGeneFunctionalHashPan.containsKey(line[indices[0]])) {
									gene = snpGeneFunctionalHashPan.get(line[indices[0]]);
									
									// pan-ethnic
									if (macs.containsKey(gene)) {
										counts = macs.get(gene);
									} else {
										macs.put(gene, counts = new int[studies.length]);
									}
									if (!line[indices[1]].equals("NA")) {
										counts[j] += Math.round(Double.parseDouble(line[indices[1]]) * Double.parseDouble(line[indices[2]]) * 2);
									}
								}
								

								if (snpGeneFunctionalHashRaceSpecific.containsKey(line[indices[0]])) {
									gene = snpGeneFunctionalHashRaceSpecific.get(line[indices[0]]);

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
		
		log.report(ext.getTime()+"\tFinished verything in "+ext.getTimeElapsed(time));
	}
	
	public static Hashtable<String,String> filterSnpGeneFunctionalHash(String filename, Hashtable<String,String> snpGeneFunctionalHash, double mafThresholdDouble, Logger log) {
		Hashtable<String,String> snpGeneFunctionalHashFiltered;
		BufferedReader reader;
		String[] line, header;
		int[] indices;
		long time;
		
		time = new Date().getTime();
		log.report(ext.getTime()+"\tFiltering "+ext.removeDirectoryInfo(filename)+" based on a "+ext.formDeci(mafThresholdDouble*100, 10, false)+"% threshold");
		
		snpGeneFunctionalHashFiltered = new Hashtable<String, String>(snpGeneFunctionalHash.size()/2);
		try {
			reader = Files.getAppropriateReader(filename);
			header = ext.splitCommasIntelligently(reader.readLine(), true, log);
			indices = ext.indexFactors(new String[] {"Name", "maf"}, header, false, true);
			while (reader.ready()) {
				line = ext.splitCommasIntelligently(reader.readLine(), true, log);
				if (snpGeneFunctionalHash.containsKey(line[indices[0]]) && !line[indices[1]].equals("NA") && Double.parseDouble(line[indices[1]]) <= mafThresholdDouble) {
					snpGeneFunctionalHashFiltered.put(line[indices[0]], snpGeneFunctionalHash.get(line[indices[0]]));
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return null;
		}
		log.report("Finished reading in " + ext.getTimeElapsed(time));
		
		return snpGeneFunctionalHashFiltered;
	}

	public static void assembleHits(String dir, String hitsDirectory, MetaAnalysisParams maps, int macThresholdStudy, int macThresholdTotal, double mafThreshold) {
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
		boolean problem;
		
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
		log.report("Using MAC thresholds of "+macThresholdStudy+" for study / "+macThresholdTotal+" for meta-analysis");
		
		problem = false;
		if (groupAnnotationParams.length > 0) {
			log.report("Testing for presence of Group annotation files");
			for (int i = 0; i < groupAnnotationParams.length; i++) {
				temp = groupAnnotationParams[i][1].split("[\\s]+")[0];
				if (Files.exists(temp)) {
					log.report("Found ", false, true);
				} else {
					log.report("Could not find ", false, true);
					problem = true;
				}
				log.report(groupAnnotationParams[i][0]+" annotation file '"+temp+"' (needed for arguments '"+groupAnnotationParams[i][1].substring(temp.length()+1)+"')");
			}
		}

		if (problem) {
			System.exit(1);
		}


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
							pvalFile = studies[j]+"_"+races[k][0]+"_pvals_mac"+macThresholdTotal+".dat";

							count = parsePvals(localRaceDir+filename, localDir+pvalFile, studies[j], methods[m], macHashesHashByRace.get(races[k][0]), macThresholdStudy, mafThreshold, log);
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
								groupParams.get(methods[m][1]).add(localRaceDir+filename+" "+Array.toStr(header, " "));
							} else {
								groupParams.get(methods[m][1]).add(localRaceDir+filename+" "+Array.toStr(Array.subArray(header, 1), " "));
							}
						}
					}
					log.report("", true, false);

					filename = races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0]+".csv";
					pvalFile = "meta_"+races[k][0]+"_pvals_mac"+macThresholdTotal+".dat";

					count = parsePvals(localRaceDir+filename, localDir+pvalFile, "Total", methods[m], macHashesHashByRace.get(races[k][0]), macThresholdTotal, mafThreshold, log);
					log.report(count+" lines of pvalues for "+filename);
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
					groupParams.get(methods[m][1]).add(k, localRaceDir+filename+" "+Array.toStr(header, " "));
					groupHits.get(methods[m][1]).incorporateFromFile(localDir+pvalFile, new int[] {0,1}, 0.001, log);
				}

				filename = phenotypes[i][0]+"_"+methods[m][0]+".csv";
				pvalFile = "meta_panEthnic_pvals_mac"+macThresholdTotal+".dat";
				
				count = parsePvals(localDir+filename, localDir+pvalFile, "Total", methods[m], macHashesHashByRace.get("PanEthnic"), macThresholdTotal, mafThreshold, log);
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
				groupParams.get(methods[m][1]).add(0, localDir+filename+" "+Array.toStr(header, " "));
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
				Files.write("java -cp /home/npankrat/vis.jar cnv.plots.QQPlot files=\""+filenames.substring(0, filenames.length()-1)+"\" maxToPlot=10", localDir+"plotQQs_mac"+macThresholdTotal+".bat");
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
					writer.println("lookup");
					writer.println(phenotypes[i][0]+"_hitters.dat 0 out="+dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+".csv");
					writer.println(Array.toStr(Array.toStringArray(groupParams.get(groups[g])), "\n"));
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + "parser_"+groups[g]+"_"+phenotypes[i][0]+".crf");
					e.printStackTrace();
				}
				Files.writeList(hits, phenotypes[i][0]+"_hitters.dat");
				Files.combine(hits, Array.toStringArray(groupParams.get(groups[g])), null, groups[g], ".", dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+".csv", log, true, true, false);
				log.report("\ncrf language for the creation of "+phenotypes[i][0]+"_"+groups[g]+".csv");
				log.report(Array.toStr(Array.toStringArray(groupParams.get(groups[g])), "\n"), true, false);
				log.report("");
			}
		}
		
		hitsDirectory = ext.verifyDirFormat(hitsDirectory);
		copyHits(dir, hitsDirectory, maps);
		log.report("Copied to "+hitsDirectory+" and regions are being delineated");
		delineateRegions(dir, hitsDirectory, maps, macThresholdTotal);
		
		lineCounts.insertElementAt("Study\tRace\tPhenotype\tMethod\tCount", 0);
		Files.writeList(Array.toStringArray(lineCounts), dir+"lineCounts.xln");
		log.report("check lineCounts.xln for completeness");
	}
	
	public static int parsePvals(String filename, String pvalFile, String study, String[] method, Hashtable<String,Hashtable<String,String>> macHashes, int macThreshold, double mafThreshold, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] header, expected;
		int index, mafIndex, macIndex, ntotalIndex, geneIndex, variantIndex;
		Hashtable<String,String> macHash;
		String[] line;
//		double mafThreshold;
		int count;
		boolean[] keeps;
		
		count = -1;
		
		if (!Files.exists(filename)) {
			System.err.println("Error - could not find '"+filename+"'; aborting");
			return count;
		}
		
		header = Files.getHeaderOfFile(filename, ",!", log);
		expected = getHeaderForMethod(method);

		keeps = Array.booleanArray(header.length, true);
		if (header[4].equals("caf")) {
			keeps[4] = false;
		}

		if (!ext.checkHeader(Array.subArray(header, keeps), expected, Array.intArray(expected.length), false, log, false)) {
			log.reportError("Error - unexpected header for file "+filename);
			System.exit(1);
		}
		
		index = -1;
		for (int h = 1; h < header.length; h++) {
			if (ext.indexOfStr(header[h], Aliases.PVALUES, false, true) >= 0) {
				if (index == -1) {
					index = h;
				} else {
					System.err.println("Error - both "+header[index]+" and "+header[h]+" are considered column headers for the p-value; using the former");
				}
			}
		}
		// generate a file with p-values for Q-Q plots and for lambdas
		if (index >= 0) {
			try {
				reader = new BufferedReader(new FileReader(filename));
				reader.readLine();
				writer = new PrintWriter(new FileWriter(pvalFile));
				count = 0;
				if (method[1].equals("SingleVariant")) {
					mafIndex = ext.indexOfStr("maf", header, false, true);
					ntotalIndex = ext.indexOfStr("ntotal", header, false, true);
					geneIndex = ext.indexOfStr("gene", header, false, true);
					variantIndex = ext.indexOfStr("Name", header, false, true);
					if (mafIndex == -1) {
						log.reportError("Error - no maf listed in single gene test result: "+filename);
					} else if (ntotalIndex == -1) {
						log.reportError("Error - no ntotal listed in single gene test result: "+filename);
					} else {
						writer.println("Variant\tpval");
						while (reader.ready()) {
							line = ext.splitCommasIntelligently(reader.readLine(), true, log);
							if (!line[mafIndex].equals("NA") && Double.parseDouble(line[mafIndex])*Double.parseDouble(line[ntotalIndex])*2 >= macThreshold && Double.parseDouble(line[mafIndex]) >= mafThreshold) {
								writer.println(line[variantIndex]+"\t"+line[index]);
							}
							count++;
						}
					}
				} else if (method[1].equals("BurdenTests")) {
					geneIndex = ext.indexOfStr("gene", header, false, true);
					macHash = macHashes.get(method[3]);
					line = macHash.get("studies").split("\t");
					macIndex = ext.indexOfStr(study, line);
					if (macIndex == -1) {
						log.reportError("Error - no minor allele counts for "+ext.rootOf(filename));
					} else {
						writer.println("Gene\tpval");
						while (reader.ready()) {
							line = ext.splitCommasIntelligently(reader.readLine(), true, log);
							if (macHash.containsKey(line[geneIndex])) {
								if (Integer.parseInt(macHash.get(line[geneIndex]).split("\t")[macIndex]) >= macThreshold) {
									writer.println(line[geneIndex]+"\t"+line[index]);
//								} else {
//									log.report(Integer.parseInt(macHash.get(line[geneIndex]).split("\t")[macIndex]) +" < "+ macThreshold);
								}
							}
							count++;
						}
					}
				} else {
					log.reportError("Error - unknown grouping variable: "+method[1]);
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + filename + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + filename + "\"");
				System.exit(2);
			}
			
		} else {
			log.reportError("Error - could not find p-value column header for file "+filename);
		}
		
		return count;
	}
	
	public static void copyHits(String dir, String hitsDirectory, MetaAnalysisParams maps) {
		String[] groups;
		String filename;
		String[][] phenotypes;

		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		hitsDirectory = ext.verifyDirFormat(hitsDirectory);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		groups = maps.getGroups();
		
		new File(dir+hitsDirectory).mkdirs();
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < groups.length; j++) {
				filename = dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[j]+".csv";
				System.out.println("cp "+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[j]+".csv "+hitsDirectory);
				Files.copyFile(filename, dir+hitsDirectory+ext.removeDirectoryInfo(filename));
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

		indices = Sort.orderTwoLayers(chrs, positions, log);
		
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
	
	public static void delineateRegions(String dir, String hitsDirectory, MetaAnalysisParams maps, int macThresholdTotal) {
		Vector<Vector<String>> filesToCat;
		Vector<String> inputsToCat, batchesToCat;
		int[] ns;
		float indexThreshold;
		Logger log, pvalThresholdsLog;
		int index;
		CountHash countHash;
		
		String[] groups, races;
		String filename;
		String[][] phenotypes, methods, results, forestInputs;
		Vector<String> additionalCols;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		hitsDirectory = ext.verifyDirFormat(hitsDirectory);
		
		log = new Logger(dir+"delineateRegions.log");
		log.report("Computing Bonferroni p-value thresholds from the MAC>="+macThresholdTotal+" results.");
		
		if (!Files.exists(dir+ext.verifyDirFormat(hitsDirectory))) {
			copyHits(dir, hitsDirectory, maps);
		}
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		methods = maps.getMethods();
		races = Matrix.extractColumn(maps.getRacesWithFilenameAliases(), 0);

		groups = new String[] {};
		filesToCat = new Vector<Vector<String>>();
		inputsToCat = new Vector<String>();
		batchesToCat = new Vector<String>();
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
				filename = phenotypes[i][0]+"/"+methods[m][0]+"/meta_panEthnic_pvals_mac"+macThresholdTotal+".dat";
				index = ext.indexOfStr(methods[m][1], groups);
				if (Files.exists(dir+filename)) {
					ns[index] = Math.max(ns[index], Files.countLines(filename, 1));
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
					additionalCols.add("PanEthnic_beta_SingleSNP");
//					additionalCols.add("CHARGE_ALL_AF");
//					additionalCols.add("single_func_region");
					additionalCols.add("Function");
				}
				if (groups[g].equals("BurdenTests")) {
					additionalCols.add("PanEthnic_nsnpsTotal_T5Count");
				}
				for (int j = 0; j < methods.length; j++) {
					if (methods[j][1].equals(groups[g])) {
						additionalCols.add("PanEthnic_p_"+methods[j][0]);
						if (groups[g].equals("SingleVariant")) {
							additionalCols.add("PanEthnic_maf_"+methods[j][0]);
							additionalCols.add("PanEthnic_ntotal_"+methods[j][0]);
//							additionalCols.add("$PanEthnic_maf_"+methods[j][0]+"*PanEthnic_ntotal_"+methods[j][0]+"=MAC_PanEthnic");
						}
						for (int k = 0; k < races.length; k++) {
							additionalCols.insertElementAt(races[k]+"_p_"+methods[j][0], k+1+additionalCols.indexOf("PanEthnic_p_"+methods[j][0]));
							if (groups[g].equals("SingleVariant")) {
								additionalCols.insertElementAt(races[k]+"_maf_"+methods[j][0], k+1+additionalCols.indexOf("PanEthnic_maf_"+methods[j][0]));
								additionalCols.add(races[k]+"_ntotal_"+methods[j][0]);
							}
						}
					}
				}
				if (Files.exists(dir+hitsDirectory+filename)) {
					results = HitWindows.determine(dir+hitsDirectory+filename, indexThreshold, 500000, indexThreshold*100, Array.toStringArray(additionalCols), log);
					if (results == null) {
						log.reportError("HitWindows result from "+filename+" was null");
					} else {
						results[0] = Array.addStrToArray("Trait", results[0], 0);
						for (int j = 1; j < results.length; j++) {
							results[j] = Array.addStrToArray(phenotypes[i][0], results[j], 0);
						}
						
						// this won't be called if there is no SingleSNP method
						if (groups[g].equals("SingleVariant") && methods[0][0].equals("SingleSNP")) {
							forestInputs = new String[results.length-1][3];
							for (int j = 0; j < forestInputs.length; j++) {
								forestInputs[j][0] = results[j+1][2];
								forestInputs[j][1] = results[j+1][0]+"/"+methods[0][0]+"/"+results[j+1][0]+"_"+methods[0][0]+".csv";
								forestInputs[j][2] = "PanEthnic "+methods[0][0]+" for "+results[j+1][0]+" (p="+ext.prettyP(results[j+1][5])+")";
							}
							Files.writeMatrix(forestInputs, ext.rootOf(dir+hitsDirectory+filename, false)+"_forestPlot.input", "\t");
							inputsToCat.add(ext.rootOf(dir+hitsDirectory+filename, false)+"_forestPlot.input");
							batchesToCat.addElement("# jcp cnv.plots.ForestPlot markerList="+hitsDirectory+ext.rootOf(filename, false)+"_forestPlot.input");
						}
						Files.writeMatrix(results, ext.rootOf(dir+hitsDirectory+filename, false)+"_regions.xln", "\t");
	//					String temp = phenotypes[i][0];
	//					Files.write(temp, dir+hitsDirectory+temp+".tmp");
	//					filesToCat.elementAt(g).add(dir+hitsDirectory+temp+".tmp");
						filesToCat.elementAt(g).add(ext.rootOf(dir+hitsDirectory+filename, false)+"_regions.xln");
					}
				} else {
					log.reportError("Error - could not find expected file: "+dir+hitsDirectory+filename);
				}
			}
		}
		
		for (int g = 0; g < groups.length; g++) {
			Files.cat(Array.toStringArray(filesToCat.elementAt(g)), dir+hitsDirectory+groups[g]+"_regions.xln", new int[0], log);
		}
		if (inputsToCat.size() > 0) {
			Files.cat(Array.toStringArray(inputsToCat), dir+hitsDirectory+"allForestPlots.input", new int[0], log);
			batchesToCat.addElement("jcp cnv.plots.ForestPlot markerList="+hitsDirectory+"allForestPlots.input");
			Files.writeList(Array.toStringArray(batchesToCat), dir+"allForestPlots.bat");
		}
	}
	
	public static void makeTables(String dir, String hitsDirectory, MetaAnalysisParams maps, int sigfigs) {
		BufferedReader reader;
		PrintWriter writer, forestWriter;
		String[] line;
		Logger log;
		boolean[] keeps, pvals;
		String[] typesToSplit = new String[] {"_p_", "_maf_", "_ntotal_"};			
		boolean[][] typeMatrix;
		Vector<String> v;
		String[] phenotypes;
		boolean error;
		String filename, outputFile;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		phenotypes = Matrix.extractColumn(maps.getPhenotypesWithFilenameAliases(true), 0);
		hitsDirectory = ext.verifyDirFormat(hitsDirectory);

		log = new Logger(dir+"makeTables.log");

		error = false;
		if (!Files.exists(dir+hitsDirectory)) {
			log.reportError("Error - directory "+dir+hitsDirectory+" does not exist");
			error = true;
		}

		if (!Files.exists(dir+hitsDirectory+"SingleVariant_regions.xln")) {
			log.reportError("Error - SingleVariant_regions.xln does not exist");
			error = true;
		}
		
		if (!Files.exists(dir+hitsDirectory+"BurdenTests_regions.xln")) {
			log.reportError("Error - BurdenTests_regions.xln does not exist");
			error = true;
		}
		
		for (int i = 0; i < phenotypes.length; i++) {
			if (!Files.exists(dir+hitsDirectory+phenotypes[i]+"_SingleVariant.csv")) {
				log.reportError("Error - "+phenotypes[i]+"_SingleVariant.csv does not exist");
				error = true;
			}
			if (!Files.exists(dir+hitsDirectory+phenotypes[i]+"_BurdenTests.csv")) {
				log.reportError("Error - "+phenotypes[i]+"_BurdenTests.csv does not exist");
				error = true;
			}
		}

		if (error) {
			return;
		}
		
		try {
			reader = Files.getAppropriateReader(dir+hitsDirectory+"SingleVariant_regions.xln");
			writer = new PrintWriter(new FileWriter(dir+hitsDirectory+"SingleVariant_regions_processed.xln"));
			line = reader.readLine().trim().split("[\\s]+");
			keeps = Array.booleanArray(line.length, true);

//			Trait	Region	Gene	Chr	 Position 	p-value	Region+Window	RegionStart	RegionStop	NumSigMarkers	NumSuggestiveMarkers	NumTotalMarkers	SizeOfRegion	PanEthnic_nsnpsTotal_T5Count	PanEthnic_p_T1Count	Whites_p_T1Count	Blacks_p_T1Count	Hispanics_p_T1Count	PanEthnic_p_T5Count	Whites_p_T5Count	Blacks_p_T5Count	Hispanics_p_T5Count	PanEthnic_p_SKATwu5	Whites_p_SKATwu5	Blacks_p_SKATwu5	Hispanics_p_SKATwu5
//			Trait	Region	MarkerName	Chr	Position	p-value	Region+Window	RegionStart	RegionStop	NumSigMarkers	NumSuggestiveMarkers	NumTotalMarkers	SizeOfRegion	SKATgene	Function	PanEthnic_p_SingleSNP	Whites_p_SingleSNP	Blacks_p_SingleSNP	Hispanics_p_SingleSNP	PanEthnic_maf_SingleSNP	Whites_maf_SingleSNP	Blacks_maf_SingleSNP	Hispanics_maf_SingleSNP	PanEthnic_ntotal_SingleSNP	Whites_ntotal_SingleSNP	Blacks_ntotal_SingleSNP	Hispanics_ntotal_SingleSNP
			keeps[1] = false; // Region # by trait
			keeps[5] = false; // min pvalue
			keeps[6] = false; // region
			keeps[7] = false; // regionstart
			keeps[8] = false; // regionstop
			keeps[11] = false; // total markers/genes
			keeps[12] = false; // size
			line[13] = "Gene";

			line = Array.subArray(line, keeps);
			typeMatrix = Matrix.booleanMatrix(typesToSplit.length, line.length, false);			

			v = new Vector<String>();
			for (int i = 1; i < line.length; i++) {
				writer.print("\t");
				for (int j = 0; j < typesToSplit.length; j++) {
					if (line[i].contains(typesToSplit[j])) {
						typeMatrix[j][i] = true;
						writer.print(line[i].substring(line[i].indexOf(typesToSplit[j])+typesToSplit[j].length(), line[i].length()));
						line[i] = line[i].substring(0, line[i].indexOf(typesToSplit[j]));
						if (line[i].equalsIgnoreCase("Whites")) {
							line[i] = "EA";
							HashVec.addIfAbsent("EA", v);
						} else if (line[i].equalsIgnoreCase("Blacks")) {
							line[i] = "AA";
							HashVec.addIfAbsent("AA", v);
						} else if (line[i].equalsIgnoreCase("Hispanics")) {
							line[i] = "HA";
							HashVec.addIfAbsent("HA", v);
						} else if (line[i].equalsIgnoreCase("Asians")) {
							line[i] = "AS";
							HashVec.addIfAbsent("AS", v);
						}
					}
				}
			}
			for (int i = 0; i < line.length; i++) {
				if (line[i].equalsIgnoreCase("PanEthnic")) {
					line[i] = Array.toStr(Array.toStringArray(v), "+");
				}
			}
			writer.println();
			writer.println(Array.toStr(line));
			
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				line = Array.subArray(line, keeps);			
				for (int i = 0; i < line.length; i++) {
					for (int j = 0; j < typesToSplit.length; j++) {
						if (typeMatrix[j][i]) {
							if (typesToSplit[j].equals("_p_")) {
								line[i] = "=\""+ext.prettyP(line[i], 2, 4, sigfigs, true)+"\"";
							} else if (typesToSplit[j].equals("_maf_")) {
								line[i] = ext.isMissingValue(line[i])?line[i]:("=\""+ext.prettyP(Double.parseDouble(line[i])*100, 1, 6, 1, true)+"%\"");
							} 
						}
					}
				}
				writer.println(Array.toStr(line));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + "BurdenTests_regions.xln" + "\" not found in "+dir+hitsDirectory);
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + dir+hitsDirectory+"BurdenTests_regions.xln" + "\"");
			return;
		}
		
		
		try {
			reader = Files.getAppropriateReader(dir+hitsDirectory+"BurdenTests_regions.xln");
			writer = new PrintWriter(new FileWriter(dir+hitsDirectory+"BurdenTests_regions_processed.xln"));
			line = reader.readLine().trim().split("[\\s]+");
			keeps = Array.booleanArray(line.length, true);

//			Trait	Region	Gene	Chr	 Position 	p-value	Region+Window	RegionStart	RegionStop	NumSigMarkers	NumSuggestiveMarkers	NumTotalMarkers	SizeOfRegion	PanEthnic_nsnpsTotal_T5Count	PanEthnic_p_T1Count	Whites_p_T1Count	Blacks_p_T1Count	Hispanics_p_T1Count	PanEthnic_p_T5Count	Whites_p_T5Count	Blacks_p_T5Count	Hispanics_p_T5Count	PanEthnic_p_SKATwu5	Whites_p_SKATwu5	Blacks_p_SKATwu5	Hispanics_p_SKATwu5
			keeps[1] = false; // Region # by trait
			keeps[5] = false; // min pvalue
			keeps[6] = false; // region
			keeps[7] = false; // regionstart
			keeps[8] = false; // regionstop
			keeps[11] = false; // total markers/genes
			keeps[12] = false; // size
			line[13] = "#variants";

			line = Array.subArray(line, keeps);			
			pvals = Array.booleanArray(line.length, false);
			
			v = new Vector<String>();
			for (int i = 1; i < line.length; i++) {
				writer.print("\t");
				if (line[i].contains("_p_")) {
					pvals[i] = true;
					writer.print(line[i].substring(line[i].indexOf("_p_")+3, line[i].length()));
					line[i] = line[i].substring(0, line[i].indexOf("_p_"));
					if (line[i].equalsIgnoreCase("Whites")) {
						line[i] = "EA";
						HashVec.addIfAbsent("EA", v);
					} else if (line[i].equalsIgnoreCase("Blacks")) {
						line[i] = "AA";
						HashVec.addIfAbsent("AA", v);
					} else if (line[i].equalsIgnoreCase("Hispanics")) {
						line[i] = "HA";
						HashVec.addIfAbsent("HA", v);
					} else if (line[i].equalsIgnoreCase("Asians")) {
						line[i] = "AS";
						HashVec.addIfAbsent("AS", v);
					}
				}
			}
			for (int i = 0; i < line.length; i++) {
				if (line[i].equalsIgnoreCase("PanEthnic")) {
					line[i] = Array.toStr(Array.toStringArray(v), "+");
				}
			}
			writer.println();
			writer.println(Array.toStr(line));
			
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				line = Array.subArray(line, keeps);			
				for (int i = 0; i < pvals.length; i++) {
					if (pvals[i]) {
						line[i] = "\""+ext.prettyP(line[i], 2, 4, sigfigs, true)+"\"";
					}
				}
				writer.println(Array.toStr(line));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + "BurdenTests_regions.xln" + "\" not found in "+dir+hitsDirectory);
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + dir+hitsDirectory+"BurdenTests_regions.xln" + "\"");
			return;
		}
		
		outputFile = "SupplementalTable2.xln";
		try {
			writer = new PrintWriter(new FileWriter(dir+hitsDirectory+outputFile));
			forestWriter = new PrintWriter(new FileWriter(dir+hitsDirectory+"subthresholdForest.input"));

			for (int i = 0; i < phenotypes.length; i++) {
				filename = phenotypes[i]+"_SingleVariant.csv";
				try {
					reader = Files.getAppropriateReader(dir+hitsDirectory+filename);
					line = ext.splitCommasIntelligently(reader.readLine(), true, log);
					keeps = Array.booleanArray(line.length, true);
					
					for (int j = 0; j < line.length; j++) {
						if (line[j].equals("gene")) {
							keeps[j] = false;
						}
						if (line[j].equals("CHARGE_ALL_AF")) {
							keeps[j] = false;
						}
						
						if (line[j].endsWith("_SingleSNP")) {
							line[j] = line[j].substring(0, line[j].lastIndexOf("_SingleSNP"));
						}
						
						line[j] = ext.replaceAllWith(line[j], new String[][] {
								{"Whites_", "EA_"},	
								{"Blacks_", "AA_"},	
								{"Hispanics_", "HA_"},	
								{"PanEthnic_", "EA+AA_"},	
//								{"_p", "_pval"},	
						});
					}
					for (int j = 28; j < line.length; j++) {
						keeps[j] = false;
					}

					if (i == 0) {
						line = Array.subArray(line, keeps);
						writer.println("Trait\t"+Array.toStr(line));
					}
					
					while (reader.ready()) {
						line = ext.splitCommasIntelligently(reader.readLine(), true, log);
						line = Array.subArray(line, keeps);			
						writer.println(ext.replaceAllWith(phenotypes[i], new String[][] { {"F7", "Factor VII"}, {"F8", "Factor VIII"} }) +"\t"+Array.toStr(line));
						forestWriter.println(line[0]+"\t"+phenotypes[i]+"/SingleSNP/"+phenotypes[i]+"_SingleSNP.csv\tPanEthnic SingleSNP for "+phenotypes[i]+" (p="+ext.prettyP(line[7])+")");
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \"" + filename + "\" not found in "+dir+hitsDirectory);
					reader.close();
					writer.close();
					forestWriter.close();
					return;
				} catch (IOException ioe) {
					log.reportError("Error reading file \"" + dir+hitsDirectory+filename + "\"");
					reader.close();
					writer.close();
					forestWriter.close();
					return;
				}
			}
			writer.close();
			forestWriter.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: could not write to file \"" + outputFile + "\" in "+dir+hitsDirectory);
			log.reportException(fnfe);
			return;
		} catch (IOException ioe) {
			log.reportError("Error writing to file \"" + dir+hitsDirectory+outputFile + "\"");
			log.reportException(ioe);
			return;
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

	public static void stashMetaResults(String dir, MetaAnalysisParams maps) {
		String[][] phenotypes;
		String[][] methods;
		String[][] races;
		String[] groups;
		String pheno, race, method;
		Vector<String> v;
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		groups = maps.getGroups();

		v = new Vector<String>();
		for (int i = 0; i < phenotypes.length; i++) {
			pheno = phenotypes[i][0];
			for (int g = 0; g < groups.length; g++) {
				v.add(dir+pheno+"/"+pheno+"_"+groups[g]+".csv");
			}
			for (int j = 0; j < methods.length; j++) {
				method = methods[j][0];
				for (int k = 0; k < races.length; k++) {
					race = races[k][0];
					v.add(dir+pheno+"/"+race+"/"+method+"/"+race+"_"+pheno+"_"+method+".csv");
				}
				v.add(dir+pheno+"/"+method+"/"+pheno+"_"+method+".csv");
			}
		}
		Zip.zip(Array.toStringArray(v), Files.getNextAvailableFilename("metaResultsStash#.zip"), new Logger(), true);
	}

	private static void parseBetasFor(String dir, MetaAnalysisParams maps, String betasFor, boolean metasOnly, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String[] files;
		String[][][] finalSets;
		String localDir;
		String filename;
		HashSet<String> markersOfInterest;
		String[][] phenotypes, races, methods;
		String[] studies;
		String temp;
		long time;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}

		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		
		time = new Date().getTime();
		dir = ext.verifyDirFormat(dir);
		
		if (!methods[0][0].startsWith("SingleSNP")) {
			System.err.println("Error - this program erroneously assumed that the first model was SingleSNP and got confused (it's actually "+methods[0][0]+"); aborting");
			return;
		}

		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(maps, files, log);
		markersOfInterest = HashVec.loadFileToHashSet(betasFor, false);

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(betasFor)+"_summary.xln"));
			writer.println("Pheno\tRace\tStudy\t"+Array.toStr(HEADER_TYPES[ext.indexOfStr(methods[0][2], ALGORITHMS)]));
			for (int i = 0; i < phenotypes.length; i++) {
				for (int k = 0; k <= races.length; k++) {
					localDir = dir+phenotypes[i][0]+"/"+(k==races.length?"":races[k][0]+"/")+methods[0][0]+"/";
					for (int j = metasOnly?studies.length:0; j <= studies.length; j++) {
						if ((k==races.length && j==studies.length) || j==studies.length || (k<races.length && !finalSets[i][j][k].equals("<missing>"))) {
							filename = (j==studies.length?"":studies[j]+"_")+(k==races.length?"":races[k][0]+"_")+phenotypes[i][0]+"_"+methods[0][0]+".csv";
							log.report(ext.getTime()+"\tReading "+filename);
							
							try {
								reader = new BufferedReader(new FileReader(localDir+filename));
								while (reader.ready()) {
									temp = reader.readLine();
									line = ext.splitCommasIntelligently(temp, true, log);
									if (markersOfInterest.contains(line[1])) {
										writer.println(phenotypes[i][0]+"\t"+(k==races.length?"PanEthnic":races[k][0])+"\t"+(j==studies.length?"Meta":studies[j])+"\t"+Array.toStr(line));
									}
								}
								reader.close();
							} catch (FileNotFoundException fnfe) {
								System.err.println("Error - could not find '"+localDir+filename+"'; aborting");
								writer.close();
								return;
							} catch (IOException ioe) {
								System.err.println("Error reading file \"" + filename + "\"");
								writer.close();
								return;
							}
						}
					}
					log.report("", true, false);
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(betasFor)+"_summary.xln");
			e.printStackTrace();
		}
		
		log.report(ext.getTimeElapsed(time));
	}
	
	private static void makeSetOfMockRdataFiles(String dir, MetaAnalysisParams maps) {
		String[] files;
		String[][][] finalSets;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		
		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(maps, files, new Logger());
		new File(dir+"mockRdata/").mkdirs();
		for (int i = 0; i < finalSets.length; i++) {
			for (int j = 0; j < finalSets[i].length; j++) {
				for (int k = 0; k < finalSets[i][j].length; k++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						files = finalSets[i][j][k].split(";");
						for (int f = 0; f < files.length; f++) {
							Files.writeList(new String[0], dir+"mockRdata/"+files[f]);
						}
					}
				}

			}
		}
	}
	
	private static Hashtable<String,String[]> loadSummaryStats(String dir, String[] studies, String[][] races, String[][] phenotypes, Logger log) {
		BufferedReader reader;
		String[] header, line;
		String trav;
		Hashtable<String,String[]> statsForStudyRacePheno;
		Vector<String> warnings;
		int raceIndex, phenoIndex;
		
		statsForStudyRacePheno = new Hashtable<String, String[]>();
		for (int i = 0; i < studies.length; i++) {
			warnings = new Vector<String>();
			if (Files.exists(dir+"summary_stats/"+studies[i]+"_summary_stats.txt")) {
				try {
					reader = Files.getAppropriateReader(dir+"summary_stats/"+studies[i]+"_summary_stats.txt");
					header = reader.readLine().trim().split("[\\s]+");
					if (!ext.checkHeader(header, SUMMARY_INFO_HEADER, false)) {
						log.reportError("Warning - invalid header for the "+studies[i]+" summary stats file");
					}
					while (reader.ready()) {
						line = reader.readLine().trim().split("\t");
						if (line.length > 1) {
							raceIndex = ext.indexOfAnyStr(line[1], races, false);
							phenoIndex = ext.indexOfAnyStr(line[3], phenotypes, false);
							if (!line[0].equalsIgnoreCase(studies[i])) {
								HashVec.addIfAbsent("Error - file "+studies[i]+"_summary_stats.txt uses a different study name in column 1: "+line[0], warnings);
							} else if (raceIndex == -1) {
								HashVec.addIfAbsent("Error - file "+studies[i]+"_summary_stats.txt uses a different race/ethnicity listing ("+line[1]+") in column 2 than are defined in the metaAnalysis.parameters", warnings);
							} else if (phenoIndex == -1) {
								HashVec.addIfAbsent("Error - file "+studies[i]+"_summary_stats.txt uses a different name for a phenotype ("+line[3]+") in column 3 than are defined in the metaAnalysis.parameters", warnings);
							} else {
								trav = studies[i]+"\t"+races[raceIndex][0]+"\t"+phenotypes[phenoIndex][0];
								if (statsForStudyRacePheno.containsKey(trav)) {
									log.reportError("Error - duplicate summary stats entry for "+trav);
									log.reportError("Previous: "+Array.toStr(statsForStudyRacePheno.get(trav)));
									log.reportError("Current : "+Array.toStr(line));
								} else {
									statsForStudyRacePheno.put(trav, line);
								}
							}
						}
					}
					reader.close();
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"summary_stats/"+studies[i]+"_summary_stats.txt" + "\"");
					System.exit(2);
				}
			} else {
				warnings.add("Warning - there is no summary stats file available for "+studies[i]);
			}
			if (warnings.size() > 0) {
				log.reportError(Array.toStr(Array.toStringArray(warnings), "\n"));
			}
		}
		log.report("Loaded info for "+statsForStudyRacePheno.size()+" Study/Race/Pheno sets");
		
		return statsForStudyRacePheno;
	}
	
	private static void summarizePhenotypes(String dir, MetaAnalysisParams maps, boolean includeSD, boolean includeRange, int numSigFigs, Logger log) {
		String[] line, mainLine, unitsLine, countsLine;
		String[] files;
		String[][][] finalSets;
		String[][] phenotypes, races;
		String[] studies;
		boolean something;
		Hashtable<String,String[]> statsForStudyRacePheno;
		Vector<String> mainTable, unitsTable, countsTable;
		double[] ages, proportionFemale;
		int[] sampleSizes;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}

		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		
		dir = ext.verifyDirFormat(dir);
		
		files = Files.list(dir, ".Rdata", false);
		log.setLevel(4);
		finalSets = identifySet(maps, files, log);

		statsForStudyRacePheno = loadSummaryStats(dir, studies, races, phenotypes, log);

		mainTable = new Vector<String>();
		mainTable.add("Study\tRace\tAge\tminMean\tmaxMean\tN\t%Female\t"+Array.toStr(Matrix.extractColumn(phenotypes, 0)));
		unitsTable = new Vector<String>();
		unitsTable.add("Study\tRace\t"+Array.toStr(Matrix.extractColumn(phenotypes, 0)));
		countsTable = new Vector<String>();
		countsTable.add("Study\tRace\t"+Array.toStr(Matrix.extractColumn(phenotypes, 0)));
		for (int k = 0; k < races.length; k++) {
			for (int j = 0; j < studies.length; j++) {
				mainLine = Array.stringArray(phenotypes.length, "");
				unitsLine = Array.stringArray(phenotypes.length, "");
				countsLine = Array.stringArray(phenotypes.length, "");
				something = false;
				ages = new double[0];
				proportionFemale = new double[0];
				sampleSizes = new int[0];
				for (int i = 0; i < phenotypes.length; i++) {
					if (!finalSets[i][j][k].equals("<missing>")) {
						if (statsForStudyRacePheno.containsKey(studies[j]+"\t"+races[k][0]+"\t"+phenotypes[i][0])) {
							line = statsForStudyRacePheno.get(studies[j]+"\t"+races[k][0]+"\t"+phenotypes[i][0]);
							mainLine[i] = ext.formDeci(Double.parseDouble(line[4]), 3, true);
							if (includeSD || includeRange) {
								mainLine[i] += " ("+(includeSD?ext.formDeci(Double.parseDouble(line[5]), numSigFigs, true):"")+(includeSD&&includeRange?"; ":"")+(includeRange?ext.formDeci(Double.parseDouble(line[6]), numSigFigs, true)+" - "+ext.formDeci(Double.parseDouble(line[7]), numSigFigs, true):"")+")";
							}
							unitsLine[i] = line[2];
							sampleSizes = Array.addIntToArray(Integer.parseInt(line[8])+Integer.parseInt(line[9]), sampleSizes);
							proportionFemale = Array.addDoubleToArray(Double.parseDouble(line[8])/(Double.parseDouble(line[8])+Double.parseDouble(line[9])), proportionFemale);
							countsLine[i] = sampleSizes[sampleSizes.length-1]+"";
							ages = Array.addDoubleToArray(Double.parseDouble(line[10]), ages);
						} else {
							mainLine[i] = "!!MISSING!!";
							unitsLine[i] = "!!MISSING!!";
							countsLine[i] = "!!MISSING!!";
						}

						something = true;
					}
				}
				if (something) {
					mainTable.add(studies[j]+"\t"+races[k][0]+"\t"+ext.formDeci(Array.mean(ages), 2)+"\t"+ext.formDeci(Array.min(ages),2)+"\t"+ext.formDeci(Array.max(ages),2)+"\t"+Array.max(sampleSizes)+"\t"+Array.mean(proportionFemale)+"\t"+Array.toStr(mainLine));
					unitsTable.add(studies[j]+"\t"+races[k][0]+"\t"+Array.toStr(unitsLine));
					countsTable.add(studies[j]+"\t"+races[k][0]+"\t"+Array.toStr(countsLine));
				}
			}
		}

		Files.writeList(Array.toStringArray(mainTable), dir+"pheno_summary_means.xln");
		Files.writeList(Array.toStringArray(unitsTable), dir+"pheno_summary_units.xln");
		Files.writeList(Array.toStringArray(countsTable), dir+"pheno_summary_counts.xln");
	}

	public static void metalCohortSensitivity(String dir, String compareToMetal, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header, originalBits;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		double[][][] data;
		int count;
		long time, total;
		String root;
		boolean gcControlOn;
		String study;
		
		
		total = new Date().getTime();
		
		gcControlOn = false;
		root = ext.rootOf(compareToMetal);
		new File(dir+root).mkdirs();

		header = null;
		try {
			reader = Files.getAppropriateReader(dir+compareToMetal);
			header = ext.splitCommasIntelligently(reader.readLine(), true, log);
			ext.checkHeader(header, new String[] {"gene", "Name", "p", "maf", "nmiss", "ntotal", "beta", "se"}, new int[] {0,1,2,3,4,5,6,7}, false, log, true);
			if (!Files.exists(dir+root+"/"+root+".xln")) {
				writer = new PrintWriter(new FileWriter(dir+root+"/"+root+".xln"));
				header[4] = "A1";
				header[5] = "A2";
				writer.println(Array.toStr(header));
				while (reader.ready()) {
					line = ext.splitCommasIntelligently(reader.readLine(), true, log);
					line[4] = "A";
					line[5] = "C";
					writer.println(Array.toStr(line));
				}
				writer.close();
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir+compareToMetal + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir+compareToMetal + "\"");
			System.exit(2);
		}

		try {
			for (int i = 8; i < header.length; i+=2) {
				study = header[i].substring(5);
				if (!Files.exists(dir+root+"/"+root+"_wo_"+study+"1.out")) {
					time = new Date().getTime();
					writer = new PrintWriter(new FileWriter(dir+root+"/"+root+".metal_wo_"+study+"_script.in"));
					writer.println("GENOMICCONTROL "+(gcControlOn?"ON":"OFF"));
					writer.println("SCHEME STDERR");
					writer.println("MARKER Name");
					writer.println("ALLELE A1 A2");
					writer.println();
					for (int j = 8; j < header.length; j+=2) {
						if (i != j) {
							writer.println("EFFECT "+header[j]);
							writer.println("STDERR "+header[j+1]);
							writer.println("PROCESS "+root+".xln");
							writer.println();
						}
					}
					writer.println("OUTFILE "+root+"_wo_"+study+" .out");
					writer.println("ANALYZE");
					writer.println("");
					writer.println("QUIT");
					writer.close();
				
					log.report("Running metal leaving out "+study+"...", false, true);
					CmdLine.run("metal < "+root+".metal_wo_"+study+"_script.in", dir+root+"/", new PrintStream(new File(dir+root+"/"+root+".metal_wo_"+study+".log")), false);
					log.report("done in " + ext.getTimeElapsed(time));
				} else {
					log.report("Using existing "+root+"_wo_"+study+"1.out");
				}
			}
		} catch (Exception e) {
			log.reportException(e);
		}

		try {
			if (!Files.exists(dir+root+"_metal_comparison.xln") || new File(dir+root+"_metal_comparison.xln").length()==0) {
				writer = new PrintWriter(new FileWriter(dir+root+"_metal_comparison.xln"));
				writer.println("Study\tAnalysis\tCorrelationPvals\tBeta_ratio\tSE_ratio\tNet_lambda\tMedianBetaWith\tMedianBetaWithout\tMedianStderrWith\tMedianStderrWithout\tLambdaWith\tLambdaWithout");
				hash = HashVec.loadFileToHashString(dir+root+"/"+root+".xln", 1, new int[] {6,7,2}, "\t", false);
				for (int i = 8; i < header.length; i+=2) {
					study = header[i].substring(5);
					try {
						reader = Files.getAppropriateReader(dir+root+"/"+root+"_wo_"+study+"1.out");
						log.report("Processing "+root+"_wo_"+study+"1.out");
						reader.readLine();
						data = new double[2][3][Files.countLines(dir+root+"/"+root+"_wo_"+study+"1.out", 1)];
						count = 0;
						while (reader.ready()) {
							line = reader.readLine().trim().split("[\\s]+");
							for (int j = 0; j < 3; j++) {
								data[0][j][count] = Double.parseDouble(line[3+j]);
							}
							originalBits = hash.get(line[0]).split("\t");
							for (int j = 0; j < 3; j++) {
								if (ext.isMissingValue(originalBits[j])) {
									log.reportError("Error - value for "+line[0]+" was "+Array.toStr(Array.subArray(line, 3, 6), "/")+" without "+study+" and was "+Array.toStr(originalBits, "/")+" with it");
									data[1][j][count] = Double.NaN;
								} else {
									data[1][j][count] = Double.parseDouble(originalBits[j]);
								}
							}
							count++;
						}
						reader.close();
						
						writer.println(
								study+"\t"+
								root+"\t"+
								ext.formDeci(Correlation.Pearson(Transformations.negativeLog10Transform(data[0][2]), Transformations.negativeLog10Transform(data[1][2]))[0], 3)+"\t"+
								ext.formDeci(Array.median(data[0][0])/Array.median(data[1][0]), 2)+"\t"+
								ext.formDeci(Array.median(data[0][1])/Array.median(data[1][1]), 2)+"\t"+
								ext.formDeci(Array.lambda(data[0][2])-Array.lambda(data[1][2]), 3)+"\t"+
								Array.median(data[0][0])+"\t"+
								Array.median(data[1][0])+"\t"+
								Array.median(data[0][1])+"\t"+
								Array.median(data[1][1])+"\t"+
								Array.lambda(data[0][2])+"\t"+
								Array.lambda(data[1][2])
								);	

					} catch (FileNotFoundException fnfe) {
						log.reportError("Error: file \"" + dir+compareToMetal + "\" not found in current directory");
					} catch (IOException ioe) {
						log.reportError("Error reading file \"" + dir+compareToMetal + "\"");
					}
				}
				writer.close();
			} else {
				log.report(root+"_metal_comparison.xln already exists!");
			}
		} catch (Exception e) {
			log.reportException(e);
		}
		
		System.out.println("Finished everything in " + ext.getTimeElapsed(total));
	}
	
	public static void batchAllMetalSensitivityAnalyses(String dir, MetaAnalysisParams maps) {
		String[] files;
		Logger log;
		String root, outputFilename, processedFilename;
		String[][] phenotypes, races;
		String[][] methods;
		String method;
		boolean runningByChr;
		Vector<String> needToBeProcessed, readyToBeConcatenated;
		String localDir;
		int maxChr;
		String chrom;

		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"allMetalSensitivityAnalyses.log");
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		runningByChr = maps.runningByChr();
		maxChr = getMaxChr();

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

		needToBeProcessed = new Vector<String>();
		readyToBeConcatenated = new Vector<String>();
		for (int i = 0; i < phenotypes.length; i++) {
			// Meta-analysis stratified by race
			for (int k = 0; k < races.length; k++) {
				for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
					chrom = chr==23?"X":(chr==24?"Y":chr+"");
					root = races[k][0]+"_"+phenotypes[i][0]+"_"+method;
					outputFilename = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+method+"/"+root+(runningByChr?"_chr"+chrom:"")+".csv";
					processedFilename = ext.rootOf(outputFilename, false)+"_metal_comparison.xln";
					if (!Files.exists(processedFilename) || new File(processedFilename).length() == 0) {
						needToBeProcessed.add(outputFilename);
					} else {
						readyToBeConcatenated.add(processedFilename);
					}
				}
			}

			// Meta-analysis of all races
			for (int chr = 1; chr <= (runningByChr?maxChr:1); chr++) {
				chrom = chr==23?"X":(chr==24?"Y":chr+"");
				root = phenotypes[i][0]+"_"+method;
				outputFilename = dir+phenotypes[i][0]+"/"+method+"/"+root+(runningByChr?"_chr"+chrom:"")+".csv";
				processedFilename = ext.rootOf(outputFilename, false)+"_metal_comparison.xln";
				if (!Files.exists(processedFilename) || new File(processedFilename).length() == 0) {
					needToBeProcessed.add(outputFilename);
				} else {
					readyToBeConcatenated.add(processedFilename);
				}
			}
		}
		
		files = Array.toStringArray(needToBeProcessed);
		needToBeProcessed.clear();
		if (files.length > 0) {
			for (int i = 0; i < files.length; i++) {
				localDir = ext.parseDirectoryOfFile(files[i]);
				Files.qsub(dir+"batchRuns/"+ext.rootOf(files[i])+".qsub", "cd "+localDir+"\njava -cp ~/" + org.genvisis.common.PSF.Java.GENVISIS + " gwas.SeqMeta dir="+localDir+" metalSensitivity="+ext.removeDirectoryInfo(files[i]), 25000, 3, 1);
				needToBeProcessed.add("qsub batchRuns/"+ext.rootOf(files[i])+".qsub");
			}
			Files.writeList(Array.toStringArray(needToBeProcessed), dir+"master.toBeProcessed");
			Files.chmod(dir+"master.toBeProcessed");
			log.report("Next run:\n./master.toBeProcessed");
		} else {
			log.report("All files have been processed; allMetalSensitivityAnalyses.xln should now be complete");
			Files.cat(Array.toStringArray(readyToBeConcatenated), dir+"allMetalSensitivityAnalyses.xln", Array.intArrayStandarddSkips(readyToBeConcatenated.size()), log);
		}
	}

	
	public static void checkMAFs(String dir, MetaAnalysisParams maps) {
		Hashtable<String, Integer> indexHash;
		PrintWriter writer;
		String[] files;
		String[][][] finalSets;
		Logger log;
		String root, outputFilename, hashFilename;
		String[][] phenotypes, races;
		String[] studies;
		String[][] methods;
		String method;
		float[] mafs;
		String[] markerNames;
		Vector<String> filenames;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"checkMAFs.log");
		phenotypes = maps.getPhenotypesWithFilenameAliases();
		studies = maps.getStudies();
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(maps, files, log);
		
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

		markerNames = HashVec.loadFileToStringArray(maps.getSnpInfoFilename(), true, new int[] {0}, false);
		indexHash = HashVec.loadToHashIndices(markerNames, log);

		for (int i = 0; i < phenotypes.length; i++) {
			hashFilename = dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_mafs.ser";
			if (Files.exists(hashFilename)) {
				log.report("Reloading from "+hashFilename);
				mafs = SerialFloatArray.load(hashFilename, false).getArray();
			} else {
				mafs = null;
				for (int k = 0; k < races.length; k++) {
					filenames = new Vector<String>();
					for (int j = 0; j < studies.length; j++) {
						if (!finalSets[i][j][k].equals("<missing>")) {
							root = studies[j]+"_"+races[k][0]+"_"+phenotypes[i][0]+"_"+method;
							outputFilename = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+method+"/"+root+".csv";
							
							log.report("Loading maf info from "+outputFilename);
							if (Files.exists(outputFilename)) {
//								mafs = fillMAFs(outputFilename, indexHash, log);
							} else {
								log.report("Error - missing expected file: "+outputFilename);
							}
						}
					}
				}
				log.report("Saving to "+hashFilename);
				new SerialFloatArray(mafs).serialize(hashFilename);
			}
			
			try {
				writer = new PrintWriter(new FileWriter(dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_mafReport.xln"));

//				writer.println();
//				writer.println();
				
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_mafReport.xln");
				e.printStackTrace();
			}

		}

		
	}

	public static void prepMAFreport(String[] markerNames, String[] filenames, String outputFilename, Logger log) {
		BufferedReader[] readers;	
        PrintWriter writer;
        String[] line, colNames;
		String temp, header;
		String delimiter;
		float[] mafs;
		int count;
		
		try {
			delimiter = null;
			readers = new BufferedReader[filenames.length];
			for (int i = 0; i < filenames.length; i++) {
				readers[i] = Files.getAppropriateReader(filenames[i]);
				header = ext.replaceAllWith(readers[i].readLine(), "\"", "");
				temp = ext.determineDelimiter(header);
				if (delimiter == null) {
					delimiter = temp;
				} else if (!temp.equals(delimiter)) {
					log.reportError("Different delimiter found in "+filenames[i]+" when trying to create "+outputFilename);
				}
				if (!ext.checkHeader(header.split(delimiter, -1), HEADER_TYPES[0], Array.intArray(HEADER_TYPES[0].length), false, log, false)) {
					log.reportError("Error - unexpected header for file '"+outputFilename+"' : "+header);
					readers[i].close();
				}
			}
			count = 0;
			writer = new PrintWriter(new FileWriter(outputFilename));
			while (readers[0].ready()) {
				for (int i = 0; i < filenames.length; i++) {
					try {
						temp = readers[i].readLine();
					} catch (Exception e) {
						log.reportError("Error reading line "+count+" from file "+filenames[i]+"; aborting");
						writer.println("trunctated at marker "+count+" ("+markerNames[count]+") due to file "+filenames[i]);
						writer.close();
						return;
					}
					if (delimiter.equals(",")) {
						line = ext.splitCommasIntelligently(temp, true, log);
					} else {
						line = temp.trim().split(delimiter, -1);
					}
					if (!markerNames[count].equals(line[1])) {
						log.reportError("Error - unsynchronized files at marker "+count+" of file "+filenames[i]+" (found "+line[1]+"; expecting "+markerNames[count]+")");
						writer.println("trunctated at marker "+count+" ("+markerNames[count]+") due to file "+filenames[i]);
						writer.close();
						return;
					}
					
					// make sure to determine mafIndex, etc. instead of relying on numbers, as old seqMeta results files won't have caf
//					hash.get(line[1]).put(studyIndex+"\t"+raceIndex, new String[] {line[3], line[5]});
				}
			}
			Files.closeAll(readers);
			writer.close();
		} catch (IOException ioe) {
			log.reportError("Error while creating file \"" + outputFilename + "\"");
			log.reportException(ioe);
		}
	}

	public static void pleiotropy(String dir, MetaAnalysisParams maps, int macThresholdTotal, double maxPvalAllowedForCountColumn) {
		Logger log;
		String localDir, localRaceDir;
		String filename;
		Vector<String> filenames;
		String[][] phenotypes, races;
		String[][] methods;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		races = maps.getRacesWithFilenameAliases();
		methods = maps.getMethods();

		log = new Logger(dir+"assemblePleiotropyTable.log");

		new File("pleiotropy/").mkdirs();
		for (int m = 0; m < methods.length; m++) {
			for (int k = 0; k < races.length; k++) {
				filenames = new Vector<String>();
				for (int i = 0; i < phenotypes.length; i++) {
					localRaceDir = dir+phenotypes[i][0]+"/"+races[k][0]+"/"+methods[m][0]+"/";
					filename = races[k][0]+"_"+phenotypes[i][0]+"_"+methods[m][0]+".csv";

					if (Files.exists(localRaceDir+filename) && Files.countIfMoreThanNLines(localRaceDir+filename, 1)) {
						filenames.add(localRaceDir+filename);
					} else {
						System.out.println("Excluding "+localRaceDir+filename);
					}
				}
				lookForPleiotropy(Array.toStringArray(filenames), macThresholdTotal, maxPvalAllowedForCountColumn, "pleiotropy/"+methods[m][0]+"_"+races[k][0]+"_pleiotropy.xln", log);
			}

			filenames = new Vector<String>();
			for (int i = 0; i < phenotypes.length; i++) {
				localDir = dir+phenotypes[i][0]+"/"+methods[m][0]+"/";
				filename = phenotypes[i][0]+"_"+methods[m][0]+".csv";

				if (Files.exists(localDir+filename) && Files.countIfMoreThanNLines(localDir+filename, 1)) {
					filenames.add(localDir+filename);
				} else {
					System.out.println("Excluding "+localDir+filename);
				}
			}
			lookForPleiotropy(Array.toStringArray(filenames), macThresholdTotal, maxPvalAllowedForCountColumn, "pleiotropy/"+methods[m][0]+"_pleiotropy.xln", log);
		}
		
//		delineateRegionsOfPleiotropy(dir, maps);
	}

	public static void lookForPleiotropy(String[] filenames, int macThresholdTotal, double maxPvalAllowedForCountColumn, String outputFile, Logger log) {
		PrintWriter writer;
		String[] line, header, params;
		String trav, temp;
		int count;
		long time;
		GenParser[] parsers;
		int[] indices;
		DoubleVector dv;
		double[] pvals;
		
		time = new Date().getTime();
		log.report("Creating "+outputFile+" using a MAC threshold of "+macThresholdTotal+" and a maximum p-value threhsold of "+ext.prettyP(maxPvalAllowedForCountColumn));
		parsers = new GenParser[filenames.length];
		for (int i = 0; i < filenames.length; i++) {
			header = Files.getHeaderOfFile(filenames[i], ",!", log);
			indices = ext.indexFactors(new String[][] {{"Name", "SNP", "gene"}, {"p"}, {"cmafUsed", "cmaf", "maf"}, {"ntotal"}}, header, true, false, true, false, log, false);
			if (Array.min(indices) == -1) {
				log.reportError("Error - header for file "+filenames[i]+"; aborting pleiotropy for "+outputFile);
				return;
			}
			
			params = Array.toStringArray(indices);
//			params = Array.addStrToArray("simplifyQuotes", params);
			params = Array.insertStringAt(filenames[i], params, 0);
			parsers[i] = new GenParser(params, log);
		}
		
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.print("Name");
			for (int i = 0; i < filenames.length; i++) {
				writer.print("\t"+ext.rootOf(filenames[i]));
			}
			writer.println("\tIndexValue\tSumOfOtherValues\tSummaryCount");
			
			while (parsers[0].ready()) {
				count = 0;
				trav = null;
				dv = new DoubleVector();
				for (int i = 0; i < parsers.length; i++) {
					if (!parsers[i].ready()) {
						log.reportError("Error - file "+filenames[i]+" is truncated compared to the first ("+filenames[0]+"); aborting the creation of "+outputFile);
						writer.close();
						GenParser.closeAllParsers(parsers);
						return;
					}
					line = parsers[i].nextLine();
					if (i==0) {
						trav = line[0];
					} else if (!line[0].equals(trav)) {
						log.reportError("Error - mismatch in unit name in file "+filenames[i]+" ("+line[1]+") versus the rest ("+trav+"); aborting the creation of "+outputFile);
						writer.close();
						GenParser.closeAllParsers(parsers);
						return;
					}
					if (ext.isMissingValue(line[1])) {
						dv.add(11.0);
					} else if (Double.parseDouble(line[2])*Double.parseDouble(line[3])*2 < macThresholdTotal || Double.parseDouble(line[1]) > maxPvalAllowedForCountColumn) {
						dv.add(2.0);
					} else {
						dv.add(Double.parseDouble(line[1]));
					}
				}
				temp = trav;
				pvals = dv.toArray();
				for (int i = 0; i < pvals.length; i++) {
					if (pvals[i] > 10) {
						temp += "\tNA";
						pvals[i] = Double.NaN;
					} else if (pvals[i] > 1) {
						temp += "\t.";
						pvals[i] = Double.NaN;
					} else {
//						temp += "\t"+ext.prettyP(pvals[i]);
						temp += "\t"+pvals[i];
					}
				}
				count = 0;
				for (int i = 0; i < pvals.length; i++) {
					if (Double.isNaN(pvals[i])) {
//						temp += "\t.";
					} else {
//						temp += "\t"+ext.formDeci(-1*Math.log10(pvals[i]), 2);
						if (pvals[i] < maxPvalAllowedForCountColumn) {
							count++;
						}
						pvals[i] = -1*Math.log10(pvals[i]);
					}
				}
				pvals = Array.removeNaN(pvals);
				pvals = Sort.putInOrder(pvals);
				if (count > 1) {
					writer.println(temp+"\t"+pvals[pvals.length-1]+"\t"+Array.sum(Array.subArray(pvals, 0, pvals.length-1))+"\t"+(count-1));
				}
			}
			writer.close();
			GenParser.closeAllParsers(parsers);
		} catch (Exception e) {
			System.err.println("Error writing to " + outputFile);
			e.printStackTrace();
		}
		
		log.report("Finished in " + ext.getTimeElapsed(time));
	}

	public static void extractFromDumps(String dir, String genoPattern, String markerList, boolean plinkFormat) {
		BufferedReader reader;
		String[] line;
		String temp, trav;
		Vector<String> v = new Vector<String>();
		int count;
		long time;
		String[][] replacements;
		
		PrintWriter writer;
		Hashtable<String, Vector<String>> hash;
		String[] markerNames, header;
		int[] indices;
		Logger log;
		String params;
		boolean hitup;
		
		hitup = false;
		log = new Logger();
		hash = HashVec.loadFileToHashVec(dir+markerList, 1, new int[] {0}, null, false, true);
		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(markerList)+"_extract.crf"));
			writer.println("lookup");
			for (int chr = 1; chr <= 24; chr++) {
				if (hash.containsKey(chr+"")) {
					markerNames = Array.toStringArray(hash.get(chr+""));
					header = Files.getHeaderOfFile(dir+ext.insertNumbers(genoPattern, chr, 1), ",!", log);
					indices = ext.indexFactors(markerNames, header, false, false);
					params = dir+ext.insertNumbers(genoPattern, chr, 1)+" 0";
					if (!hitup) {
						writer.println(dir+ext.insertNumbers(genoPattern, chr, 1)+" 0 , header out="+ext.rootOf(markerList)+"_extracted.xln lessMemoryButSlower");
						hitup = true;
//						params += plinkFormat?" 0=IID $#0=FA $#0=MO $#1=Sex $#1=Aff":"";
						params += plinkFormat?" 0=IID $@0=FA $@0=MO $@1=Sex $@1=Aff":"";
					}
					for (int i = 0; i < indices.length; i++) {
						if (indices[i] != -1) {
							params += " "+(indices[i]+1)+"="+markerNames[i];
						}
					}
//					params += plinkFormat?" 0=>AA 1=>AC 2=>CC NA=>00":"";
					writer.println(params);
				}

			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + dir+ext.rootOf(markerList)+"_extract.crf");
			log.reportException(e);
		}
		
		LookupTable.fromParameters(dir+ext.rootOf(markerList)+"_extract.crf", log);
		
		if (plinkFormat) {
			try {
				reader = Files.getAppropriateReader(ext.rootOf(markerList)+"_extracted.xln");
				header = reader.readLine().trim().split("[\\s]+");
				writer = new PrintWriter(new FileWriter(ext.rootOf(markerList)+".map"));
				for (int i = 6; i < header.length; i++) {
					line = header[i].split(":");
					if (line[0].startsWith("chr")) {
						line[0] = line[0].substring(3);
					}
					writer.println(line[0]+"\t"+header[i]+"\t0\t"+line[1]);
				}			
				writer.close();
				
				replacements = new String[][] {
						{"2", "C\tC"},
						{"1", "A\tC"},
						{"0", "A\tA"},
				};
				
				writer = new PrintWriter(new FileWriter(ext.rootOf(markerList)+".ped"));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					for (int i = 6; i < line.length; i++) {
						line[i] = ext.replaceAllWith(line[i], replacements);
						line[i] = ext.replaceAllWith(line[i], "NA", "0\t0");
					}
					writer.println(Array.toStr(line));
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + ext.rootOf(markerList)+"_extracted.xln" + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + ext.rootOf(markerList)+"_extracted.xln" + "\"");
				return;
			}
		}
	}

	public static void conatenateHits(String dir, MetaAnalysisParams maps, String hitsDirectory, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] groups;
		String filename;
		String[][] phenotypes, methods;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		hitsDirectory = ext.verifyDirFormat(hitsDirectory);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		methods = maps.getMethods();

		groups = new String[] {};
		for (int m = 0; m < methods.length; m++) {
			if (ext.indexOfStr(methods[m][1], groups) == -1) {
				groups = Array.addStrToArray(methods[m][1], groups);
			}
		}

		for (int g = 0; g < groups.length; g++) {
			try {
				writer = new PrintWriter(new FileWriter(dir+hitsDirectory+"AllPhenotypes_"+groups[g]+".csv"));
				for (int i = 0; i < phenotypes.length; i++) {
					filename = phenotypes[i][0]+"_"+groups[g]+".csv";
					
					if (Files.exists(dir+hitsDirectory+filename)) {
						try {
							reader = Files.getAppropriateReader(dir+hitsDirectory+filename);
							if (i==0) {
								writer.println("Trait,FileLocation,"+reader.readLine());
							} else {
								reader.readLine();
							}
							while (reader.ready()) {
								writer.println(phenotypes[i][0]+","+phenotypes[i][0]+"/SingleSNP/"+phenotypes[i][0]+"_SingleSNP.csv,"+reader.readLine());
							}
							reader.close();
						} catch (IOException ioe) {
							log.reportError("Error reading file \"" + dir+hitsDirectory+filename + "\"");
							return;
						}
					} else {
						log.reportError("Error - could not find expected file: "+dir+hitsDirectory+filename);
					}
				}
				writer.close();
			} catch (IOException ioe) {
				log.reportException(ioe);
			}
		}
	}

	public static void prepForestFile(String dir, MetaAnalysisParams maps, String forestMarkerList, String hitsDirectory) {
		BufferedReader reader;
		PrintWriter writer;
		String[] groups, line;
		String filename;
		String[][] phenotypes, methods;
		Logger log;
		Hashtable<String,String> hash;
		
		log = new Logger();
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		dir = ext.verifyDirFormat(dir);
		hitsDirectory = ext.verifyDirFormat(hitsDirectory);
		
		phenotypes = maps.getPhenotypesWithFilenameAliases(true);
		methods = maps.getMethods();

		groups = new String[] {};
		for (int m = 0; m < methods.length; m++) {
			if (ext.indexOfStr(methods[m][1], groups) == -1) {
				groups = Array.addStrToArray(methods[m][1], groups);
			}
		}
		
		hash = HashVec.loadFileToHashString(forestMarkerList, new int[] {0}, null, false, null, false, false, false);

		for (int g = 0; g < groups.length; g++) {
			try {
				writer = new PrintWriter(new FileWriter(ext.rootOf(forestMarkerList)+"_"+groups[g]+"_forest.input"));
				for (int i = 0; i < phenotypes.length; i++) {
					filename = phenotypes[i][0]+"_"+groups[g]+".csv";
					
					if (Files.exists(dir+hitsDirectory+filename)) {
						try {
							reader = Files.getAppropriateReader(dir+hitsDirectory+filename);
							while (reader.ready()) {
								line = ext.splitCommasIntelligently(reader.readLine(), true, log);
								if (hash.containsKey(line[0])) {
									writer.println(line[0]+"\t"+phenotypes[i][0]+"/SingleSNP/"+phenotypes[i][0]+"_SingleSNP.csv\t"+phenotypes[i][0]+" in "+line[3]+" (All p="+ext.prettyP(line[9])+"; EA p="+ext.prettyP(line[16])+"; AA p="+ext.prettyP(line[23])+"; HA p="+ext.prettyP(line[30])+")");
									hash.put(line[0], "present");
								}
							}
							reader.close();
						} catch (IOException ioe) {
							log.reportError("Error reading file \"" + dir+hitsDirectory+filename + "\"");
							return;
						}
					} else {
						log.reportError("Error - could not find expected file: "+dir+hitsDirectory+filename);
					}
				}
				writer.close();
			} catch (IOException ioe) {
				log.reportException(ioe);
			}
		}
		
		line = HashVec.getKeys(hash);
		for (int i = 0; i < line.length; i++) {
			if (!hash.get(line[i]).equals("present")) {
				log.report("Did not find "+line[i]+" in any of the hits files");
			}
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
		double mafThreshold = -1;
		int macThresholdStudy = 5;
		int macThresholdTotal = 40;
		String mapsFile = "metaAnalysis.params";
		MetaAnalysisParams maps;
		boolean consolidate = false;
		boolean metrics = false;
		boolean forceMeta = false;
		boolean copy = false;
		boolean genePositions = false;
		boolean stashMetaResults = false;
		String betasFor = null;
		boolean summarizePhenotypes = false;
		boolean includeSD = true;
		boolean includeRange = true;
		boolean mockery = false;
		String metalSensitivity = null;
		boolean metalSensitivityAll = false;
		boolean checkMAFs = false;
		boolean pleiotropy = false;
		double maxPval = 0.0001;
		int numSigFigs = 3;
		boolean tables = false;
		String hitsDirectory = "hitsAssembled/";
		boolean extractMarkers = false;
		String markerList = null;
		String genoPattern = "chr#.csv.gz";
		boolean plinkFormat = false;
		boolean metasOnly = false;
		String forestMarkerList = null;
		boolean concatenateHits = false;
		
//		metalCohortSensitivity("D:/ExomeChip/Hematology/results/DecemberPresentation/", "Whites_Hct_SingleSNP_withLRGP.csv", new Logger());
//		metalCohortSensitivity("D:/ExomeChip/Hematology/results/DecemberPresentation/", "Hct_SingleSNP.csv", new Logger());
//		metalCohortSensitivity("D:/ExomeChip/Hematology/results/08_withFHS_and_others/", "Whites_Hb_SingleSNP.csv", new Logger());
//		System.exit(1);
		
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
		"   (5) minor allele frequency threshold for meta-analysis (i.e. mafThreshold="+mafThreshold+" (default))\n" +
		"   (6) directory containing the assembled hits (i.e. hitsDir="+hitsDirectory+" (default))\n" +
		" OR\n" + 
		"   (2) copy top hit files (i.e. -copy (not the default))\n" + 
		"   (3) directory containing the assembled hits (i.e. hitsDir="+hitsDirectory+" (default))\n" +
		" OR\n" + 
		"   (2) parse regions from top hits (i.e. -regions (not the default))\n" + 
		"   (3) minor allele count threshold for meta-analysis to determine Bonferroni correction (i.e. macThresholdTotal="+macThresholdTotal+" (default))\n" + 
		"   (4) directory containing the assembled hits (i.e. hitsDir="+hitsDirectory+" (default))\n" +
		" OR\n" + 
		"   (2) make tables from the final list (i.e. -tables (not the default))\n" + 
		"   (3) number of significant digits for scientific notation (i.e. numSigFigs="+numSigFigs+" (default))\n" +
		"   (4) directory containing the assembled hits (i.e. hitsDir="+hitsDirectory+" (default))\n" +
		" OR\n" + 
		"   (2) making QQ plots is now part of assembleHits, so this the option -qq no longer works on its own\n" + 
		" OR\n" + 
		"   (2) stash meta-analysis resluts (before adding/removing a cohort) (i.e. -stashMetaResults (not the default))\n" + 
		" OR\n" + 
		"   (2) parse all results for a list of markers (i.e. betasFor=markerList.txt (not the default))\n" + 
		"   (3) parse results for meta-analyses only (i.e. metasOnly="+metasOnly+" (not the default))\n" + 
		" OR\n" + 
		"   (2) create a mock set of the same RData files for summarizing or testing purposes elsewhere (i.e. -mock (not the default))\n" + 
		" OR\n" + 
		"   (2) summarize phenotypes into a set of tables (i.e. -sumPhenos (not the default))\n" + 
		"   (3) include standard deviation (i.e. includeSD="+includeSD+" (default))\n" + 
		"   (4) include the range of values (i.e. includeRange="+includeRange+" (default))\n" +
		"   (5) number of significant digits for SD and range (i.e. numSigFigs="+numSigFigs+" (default))\n" +
		" \n" + 
		" OR\n" + 
		"   (1) filename of singleSNP results to run through a series of leave one study out METAL analyses (i.e. metalSensitivity=Whites_Hct_SingleSNP.csv (not the default))\n" +
		" OR\n" + 
		"   (2) batch all singleSNP results to be run through the METAL sensitivity analyses (i.e. -allMetalSensitivity (not the default))\n" +
		"       ( run again to concatenate all result files )\n" +
		" OR\n" + 
		"   (2) check minor allele frequencies for summed differences from other cohorts (i.e. -checkMAFs (not the default))\n" +
		" OR\n" + 
		"   (2) line up negative log p-values to check for pleiotropy (i.e. -pleiotropy (not the default))\n" +
		"   (3) minor allele count threshold to be included in the table (i.e. macThresholdTotal="+macThresholdTotal+" (default))\n" + 
		"   (4) maximum p-value allowed to be counted as a trait providing evidence of pleiotropy in count column (i.e. maxP="+maxPval+" (default))\n" + 
		" OR\n" +
		"   (2) concatenate hits files across all phenotypes (i.e. -concat (not the default))\n" +
		"   (3) directory containing the assembled hits (i.e. hitsDir="+hitsDirectory+" (default))\n" +
		" OR\n" +
		"   (2) ForestPlot - list of markers to search hits and prep all associations (i.e. forest=markerList.txt (not the default))\n" +
		"   (3) directory containing the assembled hits (i.e. hitsDir="+hitsDirectory+" (default))\n" +
		"  [future work would be to have p-value and maf filters to automatically geneate list]\n" +
		
		" \n" + 
		"   USE gwas.SeqMetaPrimary TO RUN AND gwas.SeqMeta TO META-ANALYZE/SUMMARIZE\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("maps=")) {
				mapsFile = args[i].split("=")[1];
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
			} else if (args[i].startsWith("macThresholdStudy=")) {
				macThresholdStudy = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("macThresholdTotal=")) {
				macThresholdTotal = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("mafThreshold=")) {
				mafThreshold = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-hits")) {
				hits = true;
				numArgs--;
			} else if (args[i].startsWith("hitsDir=")) {
				hitsDirectory = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-copy")) {
				copy = true;
				numArgs--;
			} else if (args[i].startsWith("-tables")) {
				tables = true;
				numArgs--;
			} else if (args[i].startsWith("-regions")) {
				regions = true;
				numArgs--;
			} else if (args[i].startsWith("-stashMetaResults")) {
				stashMetaResults = true;
				numArgs--;
			} else if (args[i].startsWith("betasFor=")) {
				betasFor = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("metasOnly=")) {
				metasOnly = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-sumPhenos")) {
				summarizePhenotypes = true;
				numArgs--;
			} else if (args[i].startsWith("includeSD=")) {
				includeSD = ext.parseBooleanArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("includeRange=")) {
				includeRange = ext.parseBooleanArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("numSigFigs=")) {
				numSigFigs = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-mock")) {
				mockery = true;
				numArgs--;
			} else if (args[i].startsWith("metalSensitivity=")) {
				metalSensitivity = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-allMetalSensitivity")) {
				metalSensitivityAll = true;
				numArgs--;
			} else if (args[i].startsWith("-checkMAFs")) {
				checkMAFs = true;
				numArgs--;
			} else if (args[i].startsWith("-pleiotropy")) {
				pleiotropy = true;
				numArgs--;
			} else if (args[i].startsWith("maxP=")) {
				maxPval = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-extract")) {
				extractMarkers = true;
				numArgs--;
			} else if (args[i].startsWith("markerList=")) {
				markerList = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("genoPattern=")) {
				genoPattern = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-plinkFormat")) {
				plinkFormat = true;
				numArgs--;
			} else if (args[i].startsWith("forest=")) {
				forestMarkerList = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-concat")) {
				concatenateHits = true;
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

//		dir = "D:/ExomeChip/Hematology/00src/CHARGE-RBC/";
//		summarizePhenotypes(dir, new MetaAnalysisParams(dir+"metaAnalysis.params", new Logger()), true, true, new Logger());
//		System.exit(1);

//		dir = "D:/ExomeChip/Hematology/results/11_accurateMACs/";
//		log = new Logger(logfile);
//		computeMACs = true;
		
//		regions = true;
		
//		dir = "D:/LITE/CHARGE-S/aric_wes_freeze4/results/";
//		hitsDirectory = "hitsAssembled_final_wAccurateMACs/";
//		tables = true;
		
//		extractMarkers = true;
//		dir = "D:/LITE/CHARGE-S/aric_wes_freeze4/EA_all/";
//		markerList = "FibrinogenConditionals.txt";
//		plinkFormat = true;
		
		try {
			log = new Logger(logfile);
			if (metalSensitivity != null) {
				metalCohortSensitivity(dir, metalSensitivity, log);
			} else if (extractMarkers) {
				extractFromDumps(dir, genoPattern, markerList, plinkFormat);
			} else {
				maps = new MetaAnalysisParams(dir+mapsFile, log);
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
//					if (mafThreshold == null) {
//						computeAllRelevantMACs(dir, maps, log);
//					} else {
//						computeMAC(dir, maps, mafThreshold, new Logger(dir+"computeMACs_forMAF_LTE_"+mafThreshold+".log"));
//					}
					computeAllRelevantMACs(dir, maps, log);
				} else if (hits) {
					assembleHits(dir, hitsDirectory, maps, macThresholdStudy, macThresholdTotal, mafThreshold);
				} else if (genePositions) {
					parseGenePositions(dir, maps);
				} else if (copy) {
					copyHits(dir, hitsDirectory, maps);
				} else if (regions) {
					delineateRegions(dir, hitsDirectory, maps, macThresholdTotal);
				} else if (tables) {
					makeTables(dir, hitsDirectory, maps, numSigFigs);
				} else if (stashMetaResults) {
					stashMetaResults(dir, maps);
				} else if (betasFor != null) {
					parseBetasFor(dir, maps, betasFor, metasOnly, log);
				} else if (mockery) {
					makeSetOfMockRdataFiles(dir, maps);
				} else if (summarizePhenotypes) {
					summarizePhenotypes(dir, maps, includeSD, includeRange, numSigFigs, log);
				} else if (metalSensitivityAll) {
					batchAllMetalSensitivityAnalyses(dir, maps);
				} else if (checkMAFs) {
					checkMAFs(dir, maps);
				} else if (pleiotropy) {
					pleiotropy(dir, maps, macThresholdTotal, maxPval);
				} else if (concatenateHits) {
					conatenateHits(dir, maps, hitsDirectory, log);
				} else if (forestMarkerList != null) {
					prepForestFile(dir, maps, forestMarkerList, hitsDirectory);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
