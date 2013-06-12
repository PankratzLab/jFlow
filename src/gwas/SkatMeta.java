package gwas;

import filesys.Hits;

import java.io.*;
import java.util.*;

import one.ChargeS;

import common.*;

public class SkatMeta {
	public static final String[][] MODELS = { // name, subroutine, filter on nonSynonSplic
		{"SingleSNP", "singlesnpMeta", ", cohortBetas = TRUE", "0"},
		{"SKAT_T5", "skatMeta", ", wts = 1, mafRange = c(0,0.05)", "1"},
		{"T5Count", "burdenMeta", ", wts = 1, mafRange = c(0,0.05)", "1"},
		{"T5MB", "burdenMeta", ", wts = function(maf){1/(maf*(1-maf))}, mafRange = c(0,0.05)", "1"},
		{"T1Count", "burdenMeta", ", wts = 1, mafRange = c(0,0.01)", "1"},
		{"T1MB", "burdenMeta", ", wts = function(maf){1/(maf*(1-maf))}, mafRange = c(0,0.01)", "1"}
	};

	private static void determineObjectNames(Logger log) {
		String[] lines, files;
		String[][] iterations;
		String dir, root, commands, filename;
		Vector<String> v, remaining;
		
		files = Files.list("./", null, ".rdata", false, false);
		log.report("There are "+files.length+" total .Rdata files");
		
		
		dir = new File("").getAbsolutePath()+"/";

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
	
	public static void splitAll(String dir, String[][] PHENOTYPES, String[] STUDIES, String snpInfoFile, String snpName, String chromName, String geneName) {
		String[] files, finalSet;
		boolean[] picks, used;
		int numMatches;
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
		used = Array.booleanArray(files.length, false);
		
		if (ext.indexOfStr(snpInfoFile, files) == -1) {
			log.reportError("Error - could not find SNP Info file '"+snpInfoFile+"'; aborting");
			return;
		} else {
			used[ext.indexOfStr(snpInfoFile, files)] = true;
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
		for (int i = 0; i < PHENOTYPES.length; i++) {
			finalSet = Array.stringArray(STUDIES.length, "<missing>");
			for (int j = 0; j < STUDIES.length; j++) {
				picks = Array.booleanArray(files.length, false);
				for (int k = 0; k < files.length; k++) {
					if (files[k].contains(STUDIES[j]) && ext.containsAny(files[k], PHENOTYPES[i])) {
						picks[k] = true;
						finalSet[j] = files[k];
						if (used[k]) {
							log.reportError("Error - file '"+files[k]+"' matches to "+STUDIES[j]+"/"+PHENOTYPES[i][0]+" but was already picked for another purpose");
						}
						used[k] = true;
					}
				}
				numMatches = Array.booleanArraySum(picks);
				if (numMatches == 0) {
					log.reportError("Warning - could not find a match for "+STUDIES[j]+"/"+PHENOTYPES[i][0]);
				} else if (numMatches > 1) {
					log.reportError("Error - found multiple matched for "+STUDIES[j]+"/"+PHENOTYPES[i][0]);
					log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
				}
			}
			
			log.report("For "+PHENOTYPES[i][0]+" identified:", true, false);
			for (int j = 0; j < STUDIES.length; j++) {
				log.report("   "+finalSet[j], true, false);
				if (!finalSet[j].equals("<missing>")) {
					localDir = dir+PHENOTYPES[i][0]+"/"+STUDIES[j]+"/";
					new File(localDir).mkdirs();

					commands = new Vector<String>();
					commands.add("load(\""+dir+snpInfoFile+"\")");
					commands.add("load(\""+dir+finalSet[j]+"\")");
					if (Files.exists(dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object")) {
						objectName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object", false, new int[] {0}, false)[0];
					} else {
						log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object"+"'");
						log.reportError("   need to first run using the -determineObjectNames option");
						return;
					}
					commands.add("ls()");
					
					for (int chr = 1; chr <= 23; chr++) {
						chrom = chr==23?"X":chr+"";
						subsetObject = STUDIES[j]+"_"+PHENOTYPES[i][0]+"_chr"+chrom;
						
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
					
					filename = dir+"batchSplits/"+STUDIES[j]+"_"+PHENOTYPES[i][0]+".R";
					Files.writeList(Array.toStringArray(commands), filename);

					Files.qsub(dir+"batchSplits/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 30000);
					toBeSplit.add("qsub "+ext.rootOf(filename)+".qsub");
				}
			}
			log.report("", true, false);			
		}
		Files.writeList(Array.toStringArray(toBeSplit), dir+"master.toBeSplit");
		Files.chmod(dir+"master.toBeSplit");
		
		numMatches = Array.booleanArraySum(used);
		if (numMatches != files.length) {
			log.reportError("Warning - did not find a match for the following file(s):");
			for (int i = 0; i < files.length; i++) {
				if (!used[i]) {
					log.reportError("  "+files[i]);
				}
			}
		}
		
		log.report("");
		log.report("Make sure to run \"qsub splitChrs.qsub\" first!!!");
	}

	public static void runAll(String dir, String[][] PHENOTYPES, String[] STUDIES, String snpName) {
		String[] files, finalSet;
		boolean[] picks, used;
		int numMatches;
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
		used = Array.booleanArray(files.length, false);
		
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
		for (int i = 0; i < PHENOTYPES.length; i++) {
			for (int j = 0; j < MODELS.length; j++) {
				localDir = dir+PHENOTYPES[i][0]+"/"+MODELS[j][0]+"/";
				new File(localDir).mkdirs();
			}
			
			finalSet = Array.stringArray(STUDIES.length, "<missing>");
			for (int j = 0; j < STUDIES.length; j++) {
				picks = Array.booleanArray(files.length, false);
				for (int k = 0; k < files.length; k++) {
					if (files[k].contains(STUDIES[j]) && ext.containsAny(files[k], PHENOTYPES[i])) {
						picks[k] = true;
						finalSet[j] = files[k];
						if (used[k]) {
							log.reportError("Error - file '"+files[k]+"' matches to "+STUDIES[j]+"/"+PHENOTYPES[i][0]+" but was already picked for another purpose");
						}
						used[k] = true;
					}
				}
				numMatches = Array.booleanArraySum(picks);
				if (numMatches == 0) {
					log.reportError("Warning - could not find a match for "+STUDIES[j]+"/"+PHENOTYPES[i][0]);
				} else if (numMatches > 1) {
					log.reportError("Error - found multiple matched for "+STUDIES[j]+"/"+PHENOTYPES[i][0]);
					log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
				}
			}
			
			log.report("For "+PHENOTYPES[i][0]+" identified:", true, false);
			for (int j = 0; j < STUDIES.length; j++) {
				log.report("   "+finalSet[j], true, false);
				if (!finalSet[j].equals("<missing>")) {
					localDir = dir+PHENOTYPES[i][0]+"/"+STUDIES[j]+"/";
					
					for (int chr = 1; chr <= 23; chr++) {
						chrom = chr==23?"X":chr+"";

						objectName = STUDIES[j]+"_"+PHENOTYPES[i][0]+"_chr"+chrom;
						snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";

						commands = new Vector<String>();
						commands.add("library(skatMeta)");
						commands.add("load(\""+dir+snpInfoFile+"\")");
						commands.add("load(\""+localDir+objectName+".RData"+"\")");
						commands.add("ls()");
						for (int k = 0; k < MODELS.length; k++) {
							root = STUDIES[j]+"_"+PHENOTYPES[i][0]+"_"+MODELS[k][0];
							if (!Files.exists(dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv")) {
								commands.add("results <- "+MODELS[k][1]+"("+objectName+", "+(MODELS[k][2].equals("1")?"SNPInfo=subset(snps_on_chr, sc_nonsynSplice==TRUE)":"SNPInfo=subsetsnps_on_chr")+", aggregateBy=\"SKATgene\", snpNames = \""+snpName+"\""+MODELS[k][2]+")");
								commands.add("write.table( results, \""+dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv\", sep=\",\", row.names = F)");
							}
						}
						count = 0;
						do {
							filename = dir+"batchRuns/"+STUDIES[j]+"_"+PHENOTYPES[i][0]+"_chr"+chrom+(count==0?"":"_"+count)+".R";
							count++;
						} while (Files.exists(filename));
						Files.writeList(Array.toStringArray(commands), filename);
	
						Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 5000);
						toBeRunIndividually.add("qsub "+ext.rootOf(filename)+".qsub");
					}
				}
			}
			log.report("", true, false);
			
			
			for (int chr = 1; chr <= 23; chr++) {
				chrom = chr==23?"X":chr+"";
				commands = new Vector<String>();
				commands.add("library(skatMeta)");
				snpInfoFile = "snpInfos/snpInfo_chr"+chrom+".RData";
				commands.add("load(\""+dir+snpInfoFile+"\")");
				
				objects = new Vector<String>();
				for (int j = 0; j < STUDIES.length; j++) {
					if (!finalSet[j].equals("<missing>")) {
						localDir = dir+PHENOTYPES[i][0]+"/"+STUDIES[j]+"/";
						objectName = STUDIES[j]+"_"+PHENOTYPES[i][0]+"_chr"+chrom;
						commands.add("load(\""+localDir+objectName+".RData"+"\")");
						objects.add(objectName);
					}
				}
				commands.add("ls()");
				commands.add("");
				for (int k = 0; k < MODELS.length; k++) {
					root = PHENOTYPES[i][0]+"_"+MODELS[k][0];
					if (!Files.exists(dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv")) {
						commands.add("results <- "+MODELS[k][1]+"("+Array.toStr(Array.toStringArray(objects), ", ")+", SNPInfo=subset(snps_on_chr, sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\", snpNames = \""+snpName+"\""+MODELS[k][2]+")");
						commands.add("write.table( results, \""+dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+"_chr"+chrom+".csv\", sep=\",\", row.names = F)");
						commands.add("");
					}
				}
				count = 0;
				do {
					filename = dir+"batchRuns/"+PHENOTYPES[i][0]+(count==0?"":"_"+count)+"_chr"+chrom+".R";
					count++;
				} while (Files.exists(filename));
				Files.writeList(Array.toStringArray(commands), filename);
	
				Files.qsub(dir+"batchRuns/"+ext.rootOf(filename)+".qsub", "cd "+dir+"\n"+"/share/apps/R-2.15.1/bin/Rscript --no-save "+filename, 10000);
				toBeRunMetad.add("qsub "+ext.rootOf(filename)+".qsub");
			}
			log.report("", true, false);			
		}
		Files.writeList(Array.toStringArray(toBeRunIndividually), dir+"master.toBeRunIndividually");
		Files.chmod(dir+"master.toBeRunIndividually");
		Files.writeList(Array.toStringArray(toBeRunMetad), dir+"master.toBeMetaAnalyzed");
		Files.chmod(dir+"master.toBeMetaAnalyzed");
		
		numMatches = Array.booleanArraySum(used);
		if (numMatches != files.length) {
			log.reportError("Warning - did not find a match for the following file(s):");
			for (int i = 0; i < files.length; i++) {
				if (!used[i]) {
					log.reportError("  "+files[i]);
				}
			}
		}
	}

	// no longer being updated, as it requires too much memory for whole exome sequencing
	public static void runAllTogether(String dir, String[][] PHENOTYPES, String[] STUDIES, String snpInfoFile, String snpName) {
		String[] files, finalSet;
		boolean[] picks, used;
		int numMatches;
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
		used = Array.booleanArray(files.length, false);
		
		if (ext.indexOfStr(snpInfoFile, files) == -1) {
			log.reportError("Error - could not find SNP Info file '"+snpInfoFile+"'; aborting");
			return;
		} else {
			used[ext.indexOfStr(snpInfoFile, files)] = true;
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
			
			finalSet = Array.stringArray(STUDIES.length, "<missing>");
			for (int j = 0; j < STUDIES.length; j++) {
				picks = Array.booleanArray(files.length, false);
				for (int k = 0; k < files.length; k++) {
					if (files[k].contains(STUDIES[j]) && ext.containsAny(files[k], PHENOTYPES[i])) {
						picks[k] = true;
						finalSet[j] = files[k];
						if (used[k]) {
							log.reportError("Error - file '"+files[k]+"' matches to "+STUDIES[j]+"/"+PHENOTYPES[i][0]+" but was already picked for another purpose");
						}
						used[k] = true;
					}
				}
				numMatches = Array.booleanArraySum(picks);
				if (numMatches == 0) {
					log.reportError("Warning - could not find a match for "+STUDIES[j]+"/"+PHENOTYPES[i][0]);
				} else if (numMatches > 1) {
					log.reportError("Error - found multiple matched for "+STUDIES[j]+"/"+PHENOTYPES[i][0]);
					log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
				}
			}
			
			log.report("For "+PHENOTYPES[i][0]+" identified:", true, false);
			for (int j = 0; j < STUDIES.length; j++) {
				log.report("   "+finalSet[j], true, false);
				if (!finalSet[j].equals("<missing>")) {
					commands = new Vector<String>();
					
					commands.add("library(skatMeta)");
					commands.add("load(\""+dir+snpInfoFile+"\")");
					commands.add("load(\""+dir+finalSet[j]+"\")");
					if (Files.exists(dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object")) {
						objectName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object", false, new int[] {0}, false)[0];
					} else {
						log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object"+"'");
						log.reportError("   need to first run using the -determineObjectNames option");
						return;
					}
					commands.add("ls()");
					commands.add("class("+objectName+") <- \"skatCohort\"");
					for (int k = 0; k < MODELS.length; k++) {
						root = STUDIES[j]+"_"+PHENOTYPES[i][0]+"_"+MODELS[k][0];
						if (!Files.exists(dir+PHENOTYPES[i][0]+"/"+MODELS[k][0]+"/"+root+".csv")) {
							commands.add(root+"_results <- "+MODELS[k][1]+"("+objectName+", SNPInfo=subset("+snpInfoName+", sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\", snpNames = \""+snpName+"\""+MODELS[k][2]+")");
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
			log.report("", true, false);
			
			
			commands = new Vector<String>();
			commands.add("library(skatMeta)");
			commands.add("load(\""+dir+snpInfoFile+"\")");
			
			objects = new Vector<String>();
			for (int j = 0; j < STUDIES.length; j++) {
				if (!finalSet[j].equals("<missing>")) {
					commands.add("load(\""+dir+finalSet[j]+"\")");
					objectName = null;
					if (Files.exists(dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object")) {
						objectName = HashVec.loadFileToStringArray(dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object", false, new int[] {0}, false)[0];
					} else {
						log.reportError("Error - could not find file '"+dir+"batchChecks/"+ext.rootOf(finalSet[j])+".object"+"'");
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
					
					commands.add(root+"_results <- "+MODELS[k][1]+"("+Array.toStr(Array.toStringArray(objects), ", ")+", SNPInfo=subset("+snpInfoName+", sc_nonsynSplice==TRUE), aggregateBy=\"SKATgene\", snpNames = \""+snpName+"\""+MODELS[k][2]+")");
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
			log.report("", true, false);			
		}
		Files.writeList(Array.toStringArray(toBeRun), dir+"master.toBeRun");
		Files.chmod(dir+"master.toBeRun");
		
		numMatches = Array.booleanArraySum(used);
		if (numMatches != files.length) {
			log.reportError("Warning - did not find a match for the following file(s):");
			for (int i = 0; i < files.length; i++) {
				if (!used[i]) {
					log.reportError("  "+files[i]);
				}
			}
		}
	}

	public static void metaAll(String dir, String[][] PHENOTYPES, String[] STUDIES, String[] GROUPS, String[][] METHODS, String[][] UNIT_OF_ANALYSIS, int[][] DEFAULT_SAMPLE_SIZES, boolean[] WEIGHTED, boolean[] SINGLE_VARIANTS, String[][] GROUP_ANNOTATION_PARAMS) {
		String[] files, finalSet;
		boolean[] picks, used;
		int numMatches;
		Logger log;
		String localDir;
		Vector<String> locals;
		Hashtable<String,Hits> groupHits;
		Hashtable<String,Vector<String>> groupParams;
		String[] groups, hits;
		String filename;
		String[] header;
		String filenames;
		boolean running;
		
		log = new Logger(dir+"metaAll.log");
		files = Files.list(dir, ".csv", false);
		used = Array.booleanArray(files.length, false);
		dir = ext.verifyDirFormat(dir);
		for (int i = 0; i < PHENOTYPES.length; i++) {
			running = false;
			groupHits = new Hashtable<String, Hits>();
			groupParams = new Hashtable<String, Vector<String>>();
			for (int j = 0; j < GROUPS.length; j++) {
				groupHits.put(GROUPS[j], new Hits());
				groupParams.put(GROUPS[j], new Vector<String>());
			}
			for (int j = 0; j < METHODS.length; j++) {
				localDir = dir+PHENOTYPES[i][0]+"/"+METHODS[j][0]+"/";
				new File(localDir).mkdirs();
				finalSet = Array.stringArray(STUDIES.length, "<missing>");
				for (int k = 0; k < STUDIES.length; k++) {
					picks = Array.booleanArray(files.length, false);
					for (int l = 0; l < files.length; l++) {
						if (files[l].contains(STUDIES[k]) && ext.containsAny(files[l], PHENOTYPES[i]) && ext.containsAny(files[l], METHODS[j])) {
							picks[l] = true;
							finalSet[k] = files[l];
							if (used[l]) {
								log.reportError("Error - file '"+files[l]+"' matches to "+STUDIES[k]+"/"+PHENOTYPES[i][0]+"/"+METHODS[j][0]+" but was already picked for another purpose");
							}
							used[l] = true;
						}
					}
					numMatches = Array.booleanArraySum(picks);
					if (numMatches == 0) {
						log.reportError("Warning - could not find a match for "+STUDIES[k]+"/"+PHENOTYPES[i][0]+"/"+METHODS[j][0]);
					} else if (numMatches > 1) {
						log.reportError("Error - found multiple matched for "+STUDIES[k]+"/"+PHENOTYPES[i][0]+"/"+METHODS[j][0]+":");
						log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
					}
				}
				log.report("For "+PHENOTYPES[i][0]+"/"+METHODS[j][0]+" identified:", true, false);
				locals = new Vector<String>();
				filenames = "";
				for (int k = 0; k < STUDIES.length; k++) {
					log.report("   "+finalSet[k], true, false);
					if (!finalSet[k].equals("<missing>")) {
						filename = STUDIES[k]+"_"+PHENOTYPES[i][0]+"_"+METHODS[j][0]+".dat";
						locals.add(filename);
						if (!Files.exists(localDir+filename)) {
							Metal.reformatResults(dir+finalSet[k], UNIT_OF_ANALYSIS[j], DEFAULT_SAMPLE_SIZES[i][k], localDir+filename, log);
						}
//						if (groupParams.get(GROUPS[j]).size() == 0) {
//							groupParams.get(GROUPS[j]).add(localDir+filename+" '"+Array.toStr(Matrix.extractColumn(Metal.MARKER_NAMES_AND_ALLELES, 0), "' '")+"' ");
//						}
						header = Files.getHeaderOfFile(localDir+filename, log);
						header[0] = "'"+header[0]+"'";
						for (int h = 1; h < header.length; h++) {
							if (header[h].contains("pval")) {
								filenames += localDir+filename+","+h+"="+STUDIES[k]+"_"+header[h]+";";
							}
							header[h] = "'"+header[h]+"'="+header[h]+"_"+STUDIES[k]+"_"+METHODS[j][0];
						}
						groupParams.get(GROUPS[j]).add(localDir+filename+" tab "+Array.toStr(header, " "));
					}
				}
				log.report("", true, false);
				
				filename = PHENOTYPES[i][0]+"_"+METHODS[j][0]+".se.metal";
				if (!Files.exists(localDir+filename+"1.out") || new File(localDir+filename+"1.out").length() < 500) {
					Metal.metaAnalyze(localDir, Array.toStringArray(locals), UNIT_OF_ANALYSIS[j], filename, Metal.SE_ANALYSIS, null, log);
					running = true;
				} else {
					groupHits.get(GROUPS[j]).incorporateFromFile(localDir+filename+"1.out", 0.001, log);
					groupParams.get(GROUPS[j]).add(0, localDir+filename+"1.out 'MarkerName' 'Allele1'=Meta_A1_"+METHODS[j][0]+" 'Allele2'=Meta_A2_"+METHODS[j][0]+" 'Effect'=Meta_beta_"+METHODS[j][0]+" 'StdErr'=Meta_se_"+METHODS[j][0]+" 'P-value'=Meta_pval_"+METHODS[j][0]+" 'Direction'=Direction_"+METHODS[j][0]);
					filenames += localDir+filename+"1.out,5=seMeta;";
				}
				if (WEIGHTED[j]) {
					filename = PHENOTYPES[i][0]+"_"+METHODS[j][0]+".wse.metal";
					if (!Files.exists(localDir+filename+"1.out") || new File(localDir+filename+"1.out").length() < 500) {
						Metal.metaAnalyze(localDir, Array.toStringArray(locals), UNIT_OF_ANALYSIS[j], filename, Metal.WEIGHTED_SE_ANALYSIS, null, log);
						running = true;
					} else {
						groupHits.get(GROUPS[j]).incorporateFromFile(localDir+filename+"1.out", 0.001, log);
						groupParams.get(GROUPS[j]).add(1, localDir+filename+"1.out 'MarkerName' 'Allele1'=Meta_wA1_"+METHODS[j][0]+" 'Allele2'=Meta_wA2_"+METHODS[j][0]+" 'Effect'=Meta_wbeta_"+METHODS[j][0]+" 'StdErr'=Meta_wse_"+METHODS[j][0]+" 'P-value'=Meta_pval_"+METHODS[j][0]+" 'Direction'=Direction_"+METHODS[j][0]);
						filenames += localDir+filename+"1.out,5=wseMeta;";
					}
				}
				filename = PHENOTYPES[i][0]+"_"+METHODS[j][0]+".pval.metal";
				if (!Files.exists(localDir+filename+"1.out") || new File(localDir+filename+"1.out").length() < 500) {
					Metal.metaAnalyze(localDir, Array.toStringArray(locals), UNIT_OF_ANALYSIS[j], filename, Metal.PVAL_ANALYSIS, null, log);
					running = true;
				} else {
					groupHits.get(GROUPS[j]).incorporateFromFile(localDir+filename+"1.out", 0.001, log);
					groupParams.get(GROUPS[j]).add(2, localDir+filename+"1.out 'MarkerName' 'Allele1'=Meta_pA1_"+METHODS[j][0]+" 'Allele2'=Meta_pA2_"+METHODS[j][0]+" 'P-value'=Meta_Nweighted_pval_"+METHODS[j][0]);
					if (SINGLE_VARIANTS[j]) {
						groupParams.get(GROUPS[j]).add(0, localDir+filename+"1.out 'MarkerName' 'Freq1'=Meta_AF");
						filenames += localDir+filename+"1.out,7=pMeta;";
					} else {
						filenames += localDir+filename+"1.out,5=pMeta;";
					}
				}
				Files.write("java -cp /home/npankrat/vis.jar cnv.plots.QQPlot files=\""+filenames.substring(0, filenames.length()-1)+"\" maxToPlot=10", localDir+"plotQQs.bat");
			}
			
			if (!running) {
				groups = HashVec.getKeys(groupHits);
				for (int g = 0; g < groups.length; g++) {
					filename = dir+PHENOTYPES[i][0]+"/"+PHENOTYPES[i][0]+"_"+groups[g]+"_hits.dat";
					groupHits.get(groups[g]).writeHits(filename);
					hits = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
					groupParams.get(groups[g]).add(0, filename+" 0 1=minPval skip=0");
					if (groups[g].equals("SingleVariant")) {
						filename = dir+"CHARGE_wGC/"+PHENOTYPES[i][0]+"/SingleSNP/"+PHENOTYPES[i][0]+"_SingleSNP.se.metal1.out";
						groupParams.get(groups[g]).add(1, filename+" 0 5=CHARGE_pval");
						filename = dir+"ESP_"+PHENOTYPES[i][0]+"_SingleSNP.csv";
						groupParams.get(groups[g]).add(2, filename+" 0 10=ESP_pval");
					} else {
						filename = dir+"CHARGE_wGC/"+PHENOTYPES[i][0]+"/T5Count/"+PHENOTYPES[i][0]+"_T5Count.se.metal1.out";
						groupParams.get(groups[g]).add(1, filename+" 0 5=CHARGE_pval");
						filename = dir+"ESP."+PHENOTYPES[i][0]+".T5.csv";
						groupParams.get(groups[g]).add(2, filename+" 0 'pval'=ESP_pval");
					}
					for (int l = 0; l < GROUP_ANNOTATION_PARAMS[g].length; l++) {
						groupParams.get(groups[g]).add(l, dir+GROUP_ANNOTATION_PARAMS[g][l]);
					}				
					Files.combine(hits, Array.toStringArray(groupParams.get(groups[g])), null, groups[g], ".", dir+PHENOTYPES[i][0]+"/"+PHENOTYPES[i][0]+"_"+groups[g]+".csv", log, true, true, false);
				}
			}
//			System.exit(1);
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
	}

	public static void stitch(String dir, String pattern) {
		String[] list;
		int[] skips;
		
		list = new String[23];
		skips = new int[23];
		for (int chr = 1; chr <= 23; chr++) {
			list[chr-1] = dir+ext.replaceAllWith(pattern, "#", chr==23?"X":chr+"");
			skips[chr-1] = chr>1?1:0;
		}
		Files.cat(list, dir+ext.replaceAllWith(pattern, "#", "All"), skips, new Logger());
		removeIndicesFromRdata(dir+ext.replaceAllWith(pattern, "#", "All"), dir+ext.addToRoot(ext.replaceAllWith(pattern, "#", "All"), "_slim"));
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
		String filename = "SkatMeta.dat";
		String logfile = null;
		Logger log;
		String dir = "";
		boolean determineObjectNames = false;
		boolean splitAll = false;
		boolean runAll = false;
		boolean metaAll = false;

		String usage = "\n" + 
		"gwas.SkatMeta requires 0-1 arguments\n" + 
		"   (1) determine object names (i.e. -determineObjectNames (not the default))\n" + 
		" OR\n" + 
		"   (1) split all (i.e. -splitAll  (not the default))\n" + 
		"   (2) directory (i.e. dir=" + dir + " (default))\n" + 
		" OR\n" + 
		"   (1) run all (i.e. -runAll  (not the default))\n" + 
		"   (2) directory (i.e. dir=" + dir + " (default))\n" + 
		" OR\n" + 
		"   (1) meta all (i.e. -metaAll  (not the default))\n" + 
		"   (2) directory (i.e. dir=" + dir + " (default))\n" + 
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
			} else if (args[i].startsWith("-runAll")) {
				runAll = true;
				numArgs--;
			} else if (args[i].startsWith("-metaAll")) {
				metaAll = true;
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
		
		try {
			log = new Logger(logfile);
			if (determineObjectNames) {
				determineObjectNames(log);
			} else if (splitAll) {
				splitAll(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_INFO_FILE, ChargeS.SNP_NAMES, ChargeS.CHROM_NAME, ChargeS.GENE_NAME);
			} else if (runAll) {
				runAll(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.SNP_NAMES);
//				runAll(dir, PHENOTYPES, STUDIES);
			} else if (metaAll) {
				metaAll(dir, ChargeS.PHENOTYPES, ChargeS.STUDIES, ChargeS.GROUPS, ChargeS.METHODS, ChargeS.UNIT_OF_ANALYSIS, ChargeS.DEFAULT_SAMPLE_SIZES, ChargeS.WEIGHTED, ChargeS.SINGLE_VARIANTS, ChargeS.GROUP_ANNOTATION_PARAMS);
//				metaAll(dir, PHENOTYPES, STUDIES, GROUPS, METHODS, UNIT_OF_ANALYSIS, DEFAULT_SAMPLE_SIZES, WEIGHTED, SINGLE_VARIANTS, GROUP_ANNOTATION_PARAMS);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
