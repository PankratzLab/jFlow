package gwas;

import filesys.Hits;
import filesys.SerialHash;

import java.io.*;
import java.util.*;

import one.ChargeS;
import parse.GenParser;

import common.*;

public class SkatMeta {
	public static final String[][] MODELS = { // name, grouping, subroutine, arguments, header type, mafThreshold //, parameters for parsing 
		{"SingleSNP", "SingleVariant", "singlesnpMeta", ", cohortBetas = TRUE", "0", "0.01"},
		{"SKAT_T5", "BurdenTests", "skatMeta", ", wts = 1, mafRange = c(0,0.05)", "2", "0.05"},
		{"T5Count", "BurdenTests", "burdenMeta", ", wts = 1, mafRange = c(0,0.05)", "1", "0.05"},
//		{"T5MB", "BurdenTests", "burdenMeta", ", wts = function(maf){1/(maf*(1-maf))}, mafRange = c(0,0.05)", "1", "0.05"},
		{"T1Count", "BurdenTests", "burdenMeta", ", wts = 1, mafRange = c(0,0.01)", "1", "0.01"},
//		{"T1MB", "BurdenTests", "burdenMeta", ", wts = function(maf){1/(maf*(1-maf))}, mafRange = c(0,0.01)", "1", "0.01"}
	};
	
	public static final String[][] HEADER_TYPES = {
		{"gene", "Name", "p", "maf", "nmiss", "ntotal", "beta", "se"}, // Single SNP
		{"gene", "p", "beta", "se", "cmafTotal", "cmafUsed", "nsnpsTotal", "nsnpsUsed", "nmiss"}, // Burden Test
		{"gene", "p", "Qmeta", "cmaf", "nmiss", "nsnps"} // SKAT test
	};
	
	public static final String[] FUNCTIONAL = {"sc_nonsynSplice", "sc_nonsyn"};

	private static void determineObjectNames(String dir, Logger log) {
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
			commands = "/share/apps/R-2.15.1/bin/Rscript --no-save [%0]";
			iterations = Matrix.toMatrix(Array.toStringArray(v));
			Files.qsub("batchChecks/checkObject", dir, -1, commands, iterations, 4000);
		}
	}
	
	public static String[][] identifySet(String[][] phenotypes, String[] studies, String[] files, String[] filesToIgnore, Logger log) {
		String[][] finalSets;
		boolean[] picks, used;
		int numMatches, index;
		
		used = Array.booleanArray(files.length, false);
		
		for (int i = 0; i < filesToIgnore.length; i++) {
			index = ext.indexOfStr(filesToIgnore[i], files);
			if (index >= 0) {
				used[index] = true;
			}
		}

		finalSets = Matrix.stringMatrix(phenotypes.length, studies.length, "<missing>");
		for (int i = 0; i < phenotypes.length; i++) {
			log.report("For "+phenotypes[i][0]+" identified:", true, false);
			for (int j = 0; j < studies.length; j++) {
				picks = Array.booleanArray(files.length, false);
				for (int k = 0; k < files.length; k++) {
					if (files[k].contains(studies[j]) && ext.containsAny(files[k], phenotypes[i])) {
						picks[k] = true;
						finalSets[i][j] = files[k];
						if (used[k]) {
							log.reportError("Error - file '"+files[k]+"' matches to "+studies[j]+"/"+phenotypes[i][0]+" but was already picked for another purpose");
						}
						used[k] = true;
					}
				}
				numMatches = Array.booleanArraySum(picks);
				if (numMatches == 0) {
					log.reportError("Warning - could not find a match for "+studies[j]+"/"+phenotypes[i][0]);
				} else if (numMatches > 1) {
					log.reportError("Error - found multiple matched for "+studies[j]+"/"+phenotypes[i][0]);
					log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
				}
				log.report("   "+finalSets[i][j], true, false);
			}
			log.report("", true, false);
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

	public static void splitAll(String dir, String[][] phenotypes, String[] studies, String snpInfoFile, String snpName, String chromName, String geneName) {
		String[] files;
		String[][] finalSets;
		Logger log;
		String localDir;
		String filename;
		Vector<String> toBeSplit, commands;
		String objectName, snpInfoName, chrom, subsetObject;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"splitAll.log");
		new File(dir+"batchSplits/").mkdir();
		files = Files.list(dir, null, ".Rdata", false, false);
		
		if (ext.indexOfStr(snpInfoFile, files) == -1) {
			log.reportError("Error - could not find SNP Info file '"+snpInfoFile+"'; aborting");
			return;
		}
		
		if (Files.exists(dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object")) {
			snpInfoName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object", false, new int[] {0}, false)[0];
		} else {
			log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object"+"'");
			log.reportError("   need to first run using the -determineObjectNames option");
			return;
		}
		
		
		commands = new Vector<String>();
		commands.add("load(\""+dir+snpInfoFile+"\")");
		commands.add("ls()");
		
		commands.add("chroms <- unique("+snpInfoName+"$"+chromName+")");
		commands.add("write.table( chroms, \"chroms.csv\", sep=\",\")");

		commands.add("for (chr in chroms) {");
		commands.add("  snps_on_chr <- snpinfo.ChargeSFreeze3.ESP[snpinfo.ChargeSFreeze3.ESP$CHROM == chr,]");
		commands.add("  filename <- paste(\"snpInfos/snpInfo_chr\", chr, \".RData\", sep='')");
		commands.add("  save(snps_on_chr, file=filename, compress=\"bzip2\")");
		commands.add("}");
		
		filename = dir+"batchSplits/splitChrs.R";
		Files.writeList(Array.toStringArray(commands), filename);

		Files.qsub(dir+"batchSplits/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"mkdir snpInfos/\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 30000);
		
		toBeSplit = new Vector<String>();
		toBeSplit.add("# make sure to run \"qsub "+ext.rootOf(filename)+".qsub\" first!!!");

		dir = ext.verifyDirFormat(dir);
		finalSets = identifySet(phenotypes, studies, files, new String[] {snpInfoFile}, log);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				if (!finalSets[i][j].equals("<missing>")) {
					localDir = dir+phenotypes[i][0]+"/"+studies[j]+"/";
					new File(localDir).mkdirs();

					commands = new Vector<String>();
					commands.add("load(\""+dir+snpInfoFile+"\")");
					commands.add("load(\""+dir+finalSets[i][j]+"\")");
					if (Files.exists(dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object")) {
						objectName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object", false, new int[] {0}, false)[0];
					} else {
						log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object"+"'");
						log.reportError("   need to first run using the -determineObjectNames option");
						return;
					}
					commands.add("ls()");
					
					for (int chr = 1; chr <= 23; chr++) {
						chrom = chr==23?"X":chr+"";
						subsetObject = studies[j]+"_"+phenotypes[i][0]+"_chr"+chrom;
						
						// filter for the gene names present on the chromosome
						commands.add("genes <- unique("+snpInfoName+"["+snpInfoName+"$"+chromName+" == \""+chrom+"\", \""+geneName+"\"])");

						// take the intersect of those actually present in the skatCohort object
						commands.add("idx <- intersect(genes, names("+objectName+"))");

						// create the skatCohort subset
						commands.add(subsetObject+" <- "+objectName+"[idx]");

						// make sure the dataset has the skatCohort class
						commands.add("class("+subsetObject+") <- \"skatCohort\"");

						// save the new file
						commands.add("save("+subsetObject+", file=\""+localDir+subsetObject+".RData\", compress=\"bzip2\")");
						commands.add("");
					}
					
					filename = dir+"batchSplits/"+studies[j]+"_"+phenotypes[i][0]+".R";
					Files.writeList(Array.toStringArray(commands), filename);

					Files.qsub(dir+"batchSplits/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 30000);
					toBeSplit.add("qsub "+ext.rootOf(filename)+".qsub");
				}
			}
		}
		Files.writeList(Array.toStringArray(toBeSplit), dir+"master.toBeSplit");
		Files.chmod(dir+"master.toBeSplit");
		
		log.report("");
		log.report("Make sure to run \"qsub splitChrs.qsub\" first!!!");
	}

	public static void runAll(String dir, String[][] phenotypes, String[] studies, String snpName) {
		String[] files;
		String[][] finalSets;
		Logger log;
		String localDir;
		String root, filename;
		Vector<String> toBeRunIndividually, toBeRunMetad, commands, objects;
		int count;
		String objectName, snpInfoFile, chrom;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"runAll.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(phenotypes, studies, files, new String[0], log);
		
		for (int chr = 1; chr <= 23; chr++) {
			chrom = chr==23?"X":chr+"";
			filename = "snpInfos/snpInfo_chr"+chrom+".RData";
			if (!Files.exists(filename)) {
				log.reportError("Error - could not find SNP Info file '"+filename+"'; aborting");
				return;
			}
		}
		
		toBeRunIndividually = new Vector<String>();
		toBeRunMetad = new Vector<String>();
		new File(dir+"batchRuns/").mkdir();
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < MODELS.length; j++) {
				localDir = dir+phenotypes[i][0]+"/"+MODELS[j][0]+"/";
				new File(localDir).mkdirs();
			}

			for (int j = 0; j < studies.length; j++) {
				if (!finalSets[i][j].equals("<missing>")) {
					localDir = dir+phenotypes[i][0]+"/"+studies[j]+"/";
					
					for (int chr = 1; chr <= 23; chr++) {
						chrom = chr==23?"X":chr+"";

						objectName = studies[j]+"_"+phenotypes[i][0]+"_chr"+chrom;
						snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";

						commands = new Vector<String>();
						commands.add("library(skatMeta)");
						commands.add("load(\""+dir+snpInfoFile+"\")");
						commands.add("load(\""+localDir+objectName+".RData"+"\")");
						commands.add("ls()");
						count = 0;
						for (int k = 0; k < MODELS.length; k++) {
							root = studies[j]+"_"+phenotypes[i][0]+"_"+MODELS[k][0];
							if (!Files.exists(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv")) {
								commands.add("results <- "+MODELS[k][2]+"("+objectName+", "+(MODELS[k][4].equals("1")?"SNPInfo=subset(snps_on_chr, sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\"":"SNPInfo=snps_on_chr")+", snpNames = \""+snpName+"\""+MODELS[k][3]+")");
								commands.add("write.table( results, \""+dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv\", sep=\",\", row.names = F)");
								count++;
							}
						}
						if (count > 0) {
							count = 0;
							do {
								filename = dir+"batchRuns/"+studies[j]+"_"+phenotypes[i][0]+"_chr"+chrom+(count==0?"":"_"+count)+".R";
								count++;
							} while (Files.exists(filename));
							Files.writeList(Array.toStringArray(commands), filename);
		
							Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 5000);
							toBeRunIndividually.add("qsub "+ext.rootOf(filename)+".qsub");
						}
					}
				}
			}
			
			for (int chr = 1; chr <= 23; chr++) {
				chrom = chr==23?"X":chr+"";
				commands = new Vector<String>();
				commands.add("library(skatMeta)");
				snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";
				commands.add("load(\""+dir+snpInfoFile+"\")");
				
				objects = new Vector<String>();
				for (int j = 0; j < studies.length; j++) {
					if (!finalSets[i][j].equals("<missing>")) {
						localDir = dir+phenotypes[i][0]+"/"+studies[j]+"/";
						objectName = studies[j]+"_"+phenotypes[i][0]+"_chr"+chrom;
						commands.add("load(\""+localDir+objectName+".RData"+"\")");
						objects.add(objectName);
					}
				}
				commands.add("ls()");
				commands.add("");
				count = 0;
				for (int k = 0; k < MODELS.length; k++) {
					root = phenotypes[i][0]+"_"+MODELS[k][0];
					if (!Files.exists(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv")) {
						commands.add("results <- "+MODELS[k][2]+"("+Array.toStr(Array.toStringArray(objects), ", ")+", "+(MODELS[k][4].equals("1")?"SNPInfo=subset(snps_on_chr, sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\"":"SNPInfo=snps_on_chr")+", snpNames = \""+snpName+"\""+MODELS[k][3]+")");
						commands.add("write.table( results, \""+dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv\", sep=\",\", row.names = F)");
						commands.add("");
						count++;
					}
				}
				if (count > 0) {
					count = 0;
					do {
						filename = dir+"batchRuns/"+phenotypes[i][0]+(count==0?"":"_"+count)+"_chr"+chrom+".R";
						count++;
					} while (Files.exists(filename));
					Files.writeList(Array.toStringArray(commands), filename);
		
					Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 10000);
					toBeRunMetad.add("qsub "+ext.rootOf(filename)+".qsub");
				}
			}
		}
		Files.writeList(Array.toStringArray(toBeRunIndividually), dir+"master.toBeRunIndividually");
		Files.chmod(dir+"master.toBeRunIndividually");
		Files.writeList(Array.toStringArray(toBeRunMetad), dir+"master.toBeMetaAnalyzed");
		Files.chmod(dir+"master.toBeMetaAnalyzed");
	}

	public static void parseAll(String dir, String[][] phenotypes, String[] studies) {
		String[] files;
		String[][] finalSets;
		Logger log;
		String root;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"parseAll.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(phenotypes, studies, files, new String[0], log);
		
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				if (!finalSets[i][j].equals("<missing>")) {
					for (int k = 0; k < MODELS.length; k++) {
						root = studies[j]+"_"+phenotypes[i][0]+"_"+MODELS[k][0];
						if (!Files.exists(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+".csv") || new File(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+".csv").length() == 0) {
							stitch(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/", root+"_chr#.csv", root+".csv");
						}
					}
				}
			}
			
			for (int k = 0; k < MODELS.length; k++) {
				root = phenotypes[i][0]+"_"+MODELS[k][0];
				if (!Files.exists(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+".csv") || new File(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/"+root+".csv").length() == 0) {
					stitch(dir+phenotypes[i][0]+"/"+MODELS[k][0]+"/", root+"_chr#.csv", root+".csv");
				}
			}
		}
	}

	public static void checkNs(String dir, String[][] phenotypes, String[] studies) {
		String[] files;
		String[][] finalSets;
		Logger log;
		String localDir;
		Vector<String> commands;
		String objectName;

		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"checkNs.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(phenotypes, studies, files, new String[0], log);
		
		
		commands = new Vector<String>();
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < phenotypes.length; i++) {
			for (int j = 0; j < studies.length; j++) {
				if (!finalSets[i][j].equals("<missing>")) {
					localDir = dir+phenotypes[i][0]+"/"+studies[j]+"/";
					objectName = studies[j]+"_"+phenotypes[i][0]+"_chr18";
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
		
		Files.writeList(Array.toStringArray(commands), dir+"checkNs.R");
		Files.qsub(dir+"checkNs.qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save checkNs.R", 10000);
	}

	// no longer being updated, as it requires too much memory for whole exome sequencing
	public static void runAllTogether(String dir, String[][] PHENOTYPES, String[] STUDIES, String snpInfoFile, String snpName) {
		String[] files;
		String[][] finalSets;
		Logger log;
		String localDir;
		String root, filename;
		Vector<String> toBeRun, commands, objects;
		int count;
		String objectName, snpInfoName;
		
		if (dir == null || dir.equals("")) {
			dir = new File("").getAbsolutePath()+"/";
		}
		
		log = new Logger(dir+"runAllTogether.log");
		files = Files.list(dir, null, ".Rdata", false, false);
		finalSets = identifySet(PHENOTYPES, STUDIES, files, new String[] {snpInfoFile}, log);
		
		if (ext.indexOfStr(snpInfoFile, files) == -1) {
			log.reportError("Error - could not find SNP Info file '"+snpInfoFile+"'; aborting");
			return;
		}
		
		if (Files.exists(dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object")) {
			snpInfoName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object", false, new int[] {0}, false)[0];
		} else {
			log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(snpInfoFile)+".object"+"'");
			log.reportError("   need to first run using the -determineObjectNames option");
			return;
		}
		
		
		toBeRun = new Vector<String>();
		new File(dir+"batchRunsSingles/").mkdir();
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < PHENOTYPES.length; i++) {
			for (int j = 0; j < MODELS.length; j++) {
				localDir = dir+PHENOTYPES[i][0]+"/"+MODELS[j][0]+"/";
				new File(localDir).mkdirs();
			}
			
			for (int j = 0; j < STUDIES.length; j++) {
				if (!finalSets[i][j].equals("<missing>")) {
					commands = new Vector<String>();
					
					commands.add("library(skatMeta)");
					commands.add("load(\""+dir+snpInfoFile+"\")");
					commands.add("load(\""+dir+finalSets[i][j]+"\")");
					if (Files.exists(dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object")) {
						objectName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object", false, new int[] {0}, false)[0];
					} else {
						log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object"+"'");
						log.reportError("   need to first run using the -determineObjectNames option");
						return;
					}
					commands.add("ls()");
					commands.add("class("+objectName+") <- \"skatCohort\"");
					for (int k = 0; k < MODELS.length; k++) {
						root = STUDIES[j]+"_"+PHENOTYPES[i][0]+"_"+MODELS[k][0];
						if (!Files.exists(dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv")) {
							commands.add(root+"_results <- "+MODELS[k][2]+"("+objectName+", SNPInfo=subset("+snpInfoName+", sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\", snpNames = \""+snpName+"\""+MODELS[k][3]+")");
							commands.add("write.table( "+root+"_results, \""+dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv\", sep=\",\", row.names = F)");

//							commands.add("fileConn_"+MODELS[k][0]+"<-file(\""+dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv\")");
//							commands.add("writeLines(c("+root+"_results), fileConn_"+MODELS[k][0]+")");
//							commands.add("close(fileConn_"+MODELS[k][0]+")");
						}
					}
					count = 0;
					do {
						filename = dir+"batchRunsSingles/"+STUDIES[j]+"_"+PHENOTYPES[i][0]+(count==0?"":"_"+count)+".R";
						count++;
					} while (Files.exists(filename));
					Files.writeList(Array.toStringArray(commands), filename);

					Files.qsub(dir+"batchRunsSingles/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 30000);
					toBeRun.add("qsub "+ext.rootOf(filename)+".qsub");
				}
			}

			commands = new Vector<String>();
			commands.add("library(skatMeta)");
			commands.add("load(\""+dir+snpInfoFile+"\")");
			
			objects = new Vector<String>();
			for (int j = 0; j < STUDIES.length; j++) {
				if (!finalSets[i][j].equals("<missing>")) {
					commands.add("load(\""+dir+finalSets[i][j]+"\")");
					objectName = null;
					if (Files.exists(dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object")) {
						objectName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object", false, new int[] {0}, false)[0];
					} else {
						log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(finalSets[i][j])+".object"+"'");
						log.reportError("   need to first run using the -determineObjectNames option");
						return;
					}
					commands.add("class("+objectName+") <- \"skatCohort\"");
					objects.add(objectName);
				}
			}
			commands.add("ls()");
			for (int k = 0; k < MODELS.length; k++) {
				root = PHENOTYPES[i][0]+"_"+MODELS[k][0];
				if (!Files.exists(dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv")) {
					
					commands.add(root+"_results <- "+MODELS[k][2]+"("+Array.toStr(Array.toStringArray(objects), ", ")+", SNPInfo=subset("+snpInfoName+", sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\", snpNames = \""+snpName+"\""+MODELS[k][3]+")");
					commands.add("write.table( "+root+"_results, \""+dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv\", sep=\",\", row.names = F)");

//					commands.add("fileConn_"+MODELS[k][0]+"<-file(\""+dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv\")");
//					commands.add("writeLines(c("+root+"_results), fileConn_"+MODELS[k][0]+")");
//					commands.add("close(fileConn_"+MODELS[k][0]+")");
				}
			}
			count = 0;
			do {
				filename = dir+"batchRunsSingles/"+PHENOTYPES[i][0]+(count==0?"":"_"+count)+".R";
				count++;
			} while (Files.exists(filename));
			Files.writeList(Array.toStringArray(commands), filename);

			Files.qsub(dir+"batchRunsSingles/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 60000);
			toBeRun.add("qsub "+ext.rootOf(filename)+".qsub");
		}
		Files.writeList(Array.toStringArray(toBeRun), dir+"master.toBeRun");
		Files.chmod(dir+"master.toBeRun");
	}

	public static void computeMAC(String dir, String[][] phenotypes, String[] studies, String snpInfoFile, String mafThreshold) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		long time;
		String[] files;
		String[][] finalSets;
		String localDir;
		String filename;
		String[] header;
		int[] indices;
		Logger log;
		Hashtable<String, String> snpGeneHash;
		Hashtable<String, String> snpGeneFunctionalHash;
		Hashtable<String, Hashtable<String, IntVector>> geneLoci;
		IntVector loci;
		String[] keys;
		Hashtable<String, Integer> largeGenes;
		int diff;
		Hashtable<String, int[]> macs;
		String gene;
		int[] counts;
		double mafThresholdDouble;
		
		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir+"computeMACs.log");
		mafThresholdDouble = Double.parseDouble(mafThreshold);
		
		log.report(ext.getTime()+"\tBegan");
		time = new Date().getTime();
		filename = dir+ext.rootOf(snpInfoFile)+".csv";	
		if (Files.exists(filename+".mappings.ser") && Files.exists(filename+".maf"+mafThreshold+".functionalMappings.ser")) {
			snpGeneHash = SerialHash.loadSerializedStringHash(filename+".mappings.ser");
			snpGeneFunctionalHash = SerialHash.loadSerializedStringHash(filename+".maf"+mafThreshold+".functionalMappings.ser");
			log.report(ext.getTime()+"\tReloaded marker mappings in " + ext.getTimeElapsed(time));
		} else {
			snpGeneHash = new Hashtable<String, String>();
			snpGeneFunctionalHash = new Hashtable<String, String>();
			geneLoci = new Hashtable<String, Hashtable<String, IntVector>>();

			try {
				reader = new BufferedReader(new FileReader(filename));
				header = ext.splitCommasIntelligently(reader.readLine(), true, log);
				indices = ext.indexFactors(new String[][] {Metal.MARKER_NAMES, Metal.GENE_UNITS, FUNCTIONAL, Metal.CHRS, Metal.POSITIONS, {"MAF"}}, header, false, true, true, log, true);
				while (reader.ready()) {
					line = ext.splitCommasIntelligently(reader.readLine(), true, log);
					snpGeneHash.put(line[indices[0]], line[indices[1]]);
					if (line[indices[2]].equals("TRUE") && !line[indices[5]].equals("NA") && Double.parseDouble(line[indices[5]]) <= mafThresholdDouble ) {
						snpGeneFunctionalHash.put(line[indices[0]], line[indices[1]]);
					}
					if (!geneLoci.containsKey(line[indices[1]])) {
						geneLoci.put(line[indices[1]], new Hashtable<String, IntVector>());
					}
					if (!geneLoci.get(line[indices[1]]).containsKey(line[indices[3]])) {
						geneLoci.get(line[indices[1]]).put(line[indices[3]], new IntVector());
					}
					geneLoci.get(line[indices[1]]).get(line[indices[3]]).add(Integer.parseInt(line[indices[4]]));
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
			
			largeGenes = new Hashtable<String, Integer>();
			keys = HashVec.getKeys(geneLoci);
			for (int i = 0; i < keys.length; i++) {
				if (geneLoci.get(keys[i]).size() > 1) {
					log.reportError("Gene '"+keys[i]+"' can be found on chromosomes "+ext.listWithCommas(HashVec.getKeys(geneLoci.get(keys[i])), true));
				} else {
					loci = geneLoci.get(keys[i]).get(HashVec.getKeys(geneLoci.get(keys[i]))[0]);
					for (int j = 0; j < loci.size(); j++) {
						for (int k = 0; k < loci.size(); k++) {
							diff = Math.abs(loci.elementAt(j)-loci.elementAt(k));
							if (diff > 2000000) {
								if (!largeGenes.contains(keys[i]+"' (n="+loci.size()+")") || diff > largeGenes.get(keys[i]+"' (n="+loci.size()+")")) {
									largeGenes.put(keys[i]+"' (n="+loci.size()+")", diff);
								}
							}
						}
					}
				}
			}
			keys = HashVec.getKeys(largeGenes);
			for (int i = 0; i < keys.length; i++) {
				if (i==0) {
					log.reportError("Genes containing variants mapped more than 2Mb apart:");
				}
				log.reportError(keys[i]+" Max: "+largeGenes.get(keys[i]));
			}
			
			SerialHash.createSerializedStringHash(filename+".maf"+mafThreshold+".mappings.ser", snpGeneHash);
			SerialHash.createSerializedStringHash(filename+".maf"+mafThreshold+".functionalMappings.ser", snpGeneFunctionalHash);

			log.report(ext.getTime()+"\tFinished mapping markers to genes in " + ext.getTimeElapsed(time));
		}
		
		
		if (!MODELS[0][0].equals("SingleSNP")) {
			System.err.println("Error - the program erroneously assumed that the first model was SingleSNP and got confused (it's actually "+MODELS[0][0]+"); aborting");
			return;
		}

		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(phenotypes, studies, files, new String[] {snpInfoFile}, log);

		for (int i = 0; i < phenotypes.length; i++) {
			macs = new Hashtable<String, int[]>();
			localDir = dir+phenotypes[i][0]+"/"+MODELS[0][0]+"/";
			for (int k = 0; k < studies.length; k++) {
				if (!finalSets[i][k].equals("<missing>")) {
					filename = studies[k]+"_"+phenotypes[i][0]+"_"+MODELS[0][0]+".csv";
					log.report(ext.getTime()+"\tReading "+filename);
					
					try {
						reader = new BufferedReader(new FileReader(localDir+filename));
						header = ext.splitCommasIntelligently(reader.readLine(), true, log);
//						ext.checkHeader(header, HEADER_TYPES[Integer.parseInt(MODELS[0][4])], Array.intArray(expected.length), false, log, true);
						
						indices = ext.indexFactors(new String[] {"Name", "maf", "ntotal"}, header, false, true);
						
						while (reader.ready()) {
							line = ext.splitCommasIntelligently(reader.readLine(), true, log);
							if (!snpGeneHash.containsKey(line[indices[0]])) {
								log.reportError("Warning - variant '"+line[indices[0]]+"' was not found in the snpInfo file");
							}
							if (snpGeneFunctionalHash.containsKey(line[indices[0]])) {
								gene = snpGeneFunctionalHash.get(line[indices[0]]);
								if (macs.containsKey(gene)) {
									counts = macs.get(gene);
								} else {
									macs.put(gene, counts = new int[studies.length]);
								}
								if (!line[indices[1]].equals("NA")) {
									counts[k] += Math.round(Double.parseDouble(line[indices[1]]) * Double.parseDouble(line[indices[2]]) * 2);
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
				writer = new PrintWriter(new FileWriter(dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+mafThreshold+".xln"));
				keys = HashVec.getKeys(macs);
				writer.println("Gene\t"+Array.toStr(studies)+"\tTotal");
				for (int j = 0; j < keys.length; j++) {
					counts = macs.get(keys[j]);
					writer.println(keys[j]+"\t"+Array.toStr(counts)+"\t"+Array.sum(counts));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+mafThreshold+".xln");
				e.printStackTrace();
			}
		}
	}

	public static void assembleHits(String dir, String[][] phenotypes, String[] studies, int[][] defaultSampleSizes, String[][] groupAnnotationParams, String snpInfoFile, int macThresholdStudy, int macThresholdTotal) {
		BufferedReader reader;
		PrintWriter writer;
		String[] files;
		String[][] finalSets;
		Logger log;
		String localDir;
//		Vector<String> locals;
		Hashtable<String,Hits> groupHits;
		Hashtable<String,Vector<String>> groupParams;
		String[] groups, hits;
		String filename;
		String[] header, expected;
		String filenames;
//		String trav; not yet used
//		int[] indices;
		int index, mafIndex, macIndex;
		Hashtable<String,Hashtable<String,String>> macHashes;
		Hashtable<String,String> macHash;
		String[] line;
		double threshold;
		String temp;
		
		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir+"metaAll.log");

		files = Files.list(dir, ".Rdata", false);
		finalSets = identifySet(phenotypes, studies, files, new String[] {snpInfoFile}, log);

		for (int i = 0; i < phenotypes.length; i++) {
			groupHits = new Hashtable<String, Hits>();
			groupParams = new Hashtable<String, Vector<String>>();
			for (int j = 0; j < MODELS.length; j++) {
				if (!groupHits.containsKey(MODELS[j][1])) {
					groupHits.put(MODELS[j][1], new Hits());
					groupParams.put(MODELS[j][1], new Vector<String>());
				}
			}
			macHashes = new Hashtable<String, Hashtable<String,String>>();
			for (int j = 0; j < MODELS.length; j++) {
				localDir = dir+phenotypes[i][0]+"/"+MODELS[j][0]+"/";
				
				if (MODELS[j][1].equals("BurdenTests")) {
					if (!macHashes.containsKey(MODELS[j][5])) {
						line = Array.addStrToArray("Total", studies);
						macHash = HashVec.loadFileToHashString(dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+MODELS[j][5]+".xln", "Gene", line, "\t");
						macHash.put("studies", Array.toStr(line));
						macHashes.put(MODELS[j][5], macHash);
					}
				}

//				locals = new Vector<String>();
				filenames = "";
				for (int k = 0; k < studies.length; k++) {
					if (!finalSets[i][k].equals("<missing>")) {
						filename = studies[k]+"_"+phenotypes[i][0]+"_"+MODELS[j][0]+".csv";
//						locals.add(filename);	// used originally for Metal analyses, may want to resurrect for checks or p-value weighted
						if (!Files.exists(localDir+filename)) {
							System.err.println("Error - could not find '"+localDir+filename+"'; aborting");
							return;
						}
						
						// not currently used, just includes all values 
//						trav = methods[j][5];
//						trav = ext.replaceAllWith(trav, "[%0]", studies[k]);

						header = Files.getHeaderOfFile(localDir+filename, ",!", log);
						expected = HEADER_TYPES[Integer.parseInt(MODELS[j][4])];
						ext.checkHeader(header, expected, Array.intArray(expected.length), false, log, true);
						
						index = -1;
						for (int h = 1; h < header.length; h++) {
							if (ext.indexOfStr(header[h], Metal.PVALUES, false, true) >= 0) {
								if (index == -1) {
									index = h;
								} else {
									System.err.println("Error - both "+header[index]+" and "+header[h]+" are considered column headers for the p-value; using the former");
								}
							}
						}
						if (index >= 0) {
							try {
								reader = new BufferedReader(new FileReader(localDir+filename));
								reader.readLine();
								writer = new PrintWriter(new FileWriter(localDir+studies[k]+"_pvals.dat"));
								if (MODELS[j][1].equals("SingleVariant")) {
									mafIndex = ext.indexOfStr("maf", header, false, true);
									if (mafIndex == -1) {
										log.reportError("Error - no maf listed in single gene test result: "+filename);
									} else {
										writer.println("Variant\tpval");
										threshold = Double.parseDouble(MODELS[j][5]);
										while (reader.ready()) {
											line = ext.splitCommasIntelligently(reader.readLine(), true, log);
											if (!line[mafIndex].equals("NA") && Double.parseDouble(line[mafIndex]) >= threshold) {
												writer.println(line[0]+"\t"+line[index]);
											}
										}
									}
								} else if (MODELS[j][1].equals("BurdenTests")) {
									macHash = macHashes.get(MODELS[j][5]);
									line = macHash.get("studies").split("\t");
									macIndex = ext.indexOfStr(studies[k], line);
									if (macIndex == -1) {
										log.reportError("Error - no minor allele counts for "+studies[k]+" "+MODELS[j][0]);
									} else {
										writer.println("Gene\tpval");
										while (reader.ready()) {
											line = ext.splitCommasIntelligently(reader.readLine(), true, log);
											if (macHash.containsKey(line[0])) {
												if (Integer.parseInt(macHash.get(line[0]).split("\t")[macIndex]) >= macThresholdStudy) {
													writer.println(line[0]+"\t"+line[index]);
												}
											}
										}
									}									
								} else {
									log.reportError("Error - unknown grouping variable: "+MODELS[j][1]);
								}
								reader.close();
								writer.close();
							} catch (FileNotFoundException fnfe) {
								System.err.println("Error: file \"" + localDir+filename + "\" not found in current directory");
								System.exit(1);
							} catch (IOException ioe) {
								System.err.println("Error reading file \"" + localDir+filename + "\"");
								System.exit(2);
							}
							filenames += localDir+studies[k]+"_pvals.dat,1="+studies[k]+";";
						} else {
							log.reportError("Error - could not find p-value column header for file "+filename);
						}
						header[0] = "'"+header[0]+"'";
						for (int h = 1; h < header.length; h++) {
							header[h] = "'"+header[h]+"'="+header[h]+"_"+studies[k]+"_"+MODELS[j][0];
						}
						if (MODELS[j][1].equals("BurdenTests")) {
							groupParams.get(MODELS[j][1]).add(localDir+filename+" simplifyQuotes "+Array.toStr(header, " "));
						}
					}
				}
				log.report("", true, false);

				filename = phenotypes[i][0]+"_"+MODELS[j][0]+".csv";
				if (!Files.exists(localDir+filename)) {
					System.err.println("Error - could not find '"+localDir+filename+"'; aborting");
					return;
				}

				
				header = Files.getHeaderOfFile(localDir+filename, ",!", log);
				expected = HEADER_TYPES[Integer.parseInt(MODELS[j][4])];
				ext.checkHeader(header, expected, Array.intArray(expected.length), false, log, true);
				
				index = -1;
				for (int h = 1; h < header.length; h++) {
					if (ext.indexOfStr(header[h], Metal.PVALUES, false, true) >= 0) {
						if (index == -1) {
							index = h;
						} else {
							System.err.println("Error - both "+header[index]+" and "+header[h]+" are considered column headers for the p-value; using the former");
						}
					}
				}
				if (index >= 0) {
					try {
						reader = new BufferedReader(new FileReader(localDir+filename));
						reader.readLine();
						writer = new PrintWriter(new FileWriter(localDir+"meta_pvals.dat"));
						if (MODELS[j][1].equals("SingleVariant")) {
							mafIndex = ext.indexOfStr("maf", header, false, true);
							if (mafIndex == -1) {
								log.reportError("Error - no maf listed in single gene test result: "+filename);
							} else {
								writer.println("Variant\tpval");
								threshold = Double.parseDouble(MODELS[j][5]);
								while (reader.ready()) {
									line = ext.splitCommasIntelligently(reader.readLine(), true, log);
									if (!line[mafIndex].equals("NA") && Double.parseDouble(line[mafIndex]) >= threshold) {
										writer.println(line[1]+"\t"+line[index]);
									}
								}
							}
						} else if (MODELS[j][1].equals("BurdenTests")) {
							macHash = macHashes.get(MODELS[j][5]);
							line = macHash.get("studies").split("\t");
							macIndex = ext.indexOfStr("Total", line);
							if (macIndex == -1) {
								log.reportError("Error - no minor allele counts at all");
							} else {
								writer.println("Gene\tpval");
								while (reader.ready()) {
									line = ext.splitCommasIntelligently(reader.readLine(), true, log);
									if (macHash.containsKey(line[0])) {
										if (Integer.parseInt(macHash.get(line[0]).split("\t")[macIndex]) >= macThresholdTotal) {
											if (line[0].equals("CLEC9A")) {
												System.out.println("hola");
												String temp2 = macHash.get(line[0]);
												System.out.println(temp2);
												String temp3 = temp2.split("\t")[macIndex];
												System.out.println(temp3);
												writer.println(line[0]+"\t"+line[index]);
											}
											writer.println(line[0]+"\t"+line[index]);
										}
									}
								}
							}									
						} else {
							log.reportError("Error - unknown grouping variable: "+MODELS[j][1]);
						}
						reader.close();
						writer.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \"" + localDir+filename + "\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \"" + localDir+filename + "\"");
						System.exit(2);
					}
					filenames += localDir+"meta_pvals.dat,1=Meta;";
				} else {
					log.reportError("Error - could not find p-value column header for file "+filename);
				}

				header[0] = "'"+header[0]+"'";
				for (int h = 1; h < header.length; h++) {
					header[h] = "'"+header[h]+"'=Meta_"+header[h]+"_"+MODELS[j][0];
				}
				if (MODELS[j][1].equals("SingleVariant")) {
					temp = header[0];
					header[0] = header[1];
					header[1] = temp;
				}
				groupParams.get(MODELS[j][1]).add(0, localDir+filename+" simplifyQuotes "+Array.toStr(header, " "));
				groupHits.get(MODELS[j][1]).incorporateFromFile(localDir+"meta_pvals.dat", new int[] {0,1}, 0.001, log);
				

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
					System.err.println("Error - why are there no files for "+phenotypes[i][0]+" "+MODELS[j][0]);
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
					line = HashVec.getKeys(macHashes);
					for (int k = 0; k < line.length; k++) {
						groupParams.get(groups[g]).add(1, dir+phenotypes[i][0]+"/"+"minorAlleleCounts.maf"+line[k]+".xln 0 'Total'=MAC>"+line[k]+"%");
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
				for (int l = 0; l < groupAnnotationParams[g].length; l++) {
					groupParams.get(groups[g]).add(l, dir+groupAnnotationParams[g][l]);
				}				
				Files.combine(hits, Array.toStringArray(groupParams.get(groups[g])), null, groups[g], ".", dir+phenotypes[i][0]+"/"+phenotypes[i][0]+"_"+groups[g]+".csv", log, true, true, false);
			}
		}
	}

	public static void stitch(String dir, String pattern, String fileout) {
		String[] list;
		int[] skips;
		
		list = new String[23];
		skips = new int[23];
		for (int chr = 1; chr <= 23; chr++) {
			list[chr-1] = dir+ext.replaceAllWith(pattern, "#", chr==23?"X":chr+"");
			skips[chr-1] = chr>1?1:0;
		}
		Files.cat(list, dir+fileout, skips, new Logger());
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
		boolean hits = false;
		boolean computeMACs = false;
		String mafThreshold = "0.05";
		int macThresholdStudy = 5;
		int macThresholdTotal = 40;

		String usage = "\n" + 
		"gwas.SkatMeta requires 0-1 arguments\n" + 
		"   (0) directory (i.e. dir=" + dir + " (default))\n" + 
		" AND\n" + 
		"   (1) determine object names (i.e. -determineObjectNames (not the default))\n" + 
		" OR\n" + 
		"   (1) split all (i.e. -splitAll (not the default))\n" + 
		" OR\n" + 
		"   (1) run all (i.e. -runAll (not the default))\n" + 
		" OR\n" + 
		"   (1) parse all runs (i.e. -parseAll (not the default))\n" + 
		" OR\n" + 
		"   (1) check sample sizes for each chr18 data file (i.e. -checkNs (not the default))\n" + 
		" OR\n" + 
		"   (1) compute minor allele counts (i.e. -computeMACs (not the default))\n" + 
		"   (2) minor allele frequency threshold (i.e. mafThreshold="+mafThreshold+" (default))\n" + 
		" OR\n" + 
		"   (1) parse top hits (i.e. -hits (not the default))\n" + 
		"   (2) minor allele count threshold for a study (i.e. macThresholdStudy="+macThresholdStudy+" (default))\n" + 
		"   (3) minor allele count threshold for meta-analysis (i.e. macThresholdTotal="+macThresholdTotal+" (default))\n" + 
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
			} else if (args[i].startsWith("-parseAll")) {
				parseAll = true;
				numArgs--;
			} else if (args[i].startsWith("-checkNs")) {
				checkNs = true;
				numArgs--;
			} else if (args[i].startsWith("-runAll")) {
				runAll = true;
				numArgs--;
			} else if (args[i].startsWith("-computeMACs")) {
				computeMACs = true;
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
		
		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/metaAnalysis/";
//		splitAll = true;
//		runAll = true;
//		computeMAC(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, "0.01");
//		computeMAC(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, "0.05");
//		System.exit(1);
		hits = true;
		
		try {
			log = new Logger(logfile);
			if (determineObjectNames) {
				determineObjectNames(dir, log);
			} else if (splitAll) {
				splitAll(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, ChargeS.SNP_NAMES, ChargeS.CHROM_NAME, ChargeS.GENE_NAME);
			} else if (parseAll) {
				parseAll(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES);
			} else if (checkNs) {
				checkNs(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES);
			} else if (runAll) {
				runAll(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_NAMES);
//				runAll(dir, PHENOTYPES, STUDIES);
			} else if (computeMACs) {
				computeMAC(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, mafThreshold);
			} else if (hits) {
				assembleHits(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.FREEZE3_SAMPLE_SIZES, ChargeS.GROUP_ANNOTATION_PARAMS, ChargeS.SNP_INFO_FILE, macThresholdStudy, macThresholdTotal);
//				metaAll(dir, PHENOTYPES, STUDIES, GROUPS, METHODS, UNIT_OF_ANALYSIS, DEFAULT_SAMPLE_SIZES, WEIGHTED, SINGLE_VARIANTS, GROUP_ANNOTATION_PARAMS);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
