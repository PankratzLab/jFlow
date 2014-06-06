// -Xmx4096M
package gwas;

import java.io.*;
import java.util.*;

import common.*;
import filesys.Hits;

public class Plink {
//	public static final String[] CLUSTER_HEADER = {"FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2", "PI_HAT", "PHE", "IBS0", "IBS1", "IBS2", "DST", "P", "HOMHOM", "HETHET", "RATIO"};
	public static final String[] CLUSTER_HEADER = {"FID1", "IID1", "FID2", "IID2", "RT", "EZ", "Z0", "Z1", "Z2", "PI_HAT", "PHE", "DST", "PPC", "RATIO"};
	public static final int[] RELEVANT_CLUSTER_INDICES = {0, 1, 2, 3, 6, 7, 8, 9};
	public static final String[] IMISS_HEADER = {"FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS"};
	public static final String[] MPERM_HEADER = {"CHR", "SNP", "EMP1", "EMP2"};
	public static final String[] LOGISTIC_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "SE", "L95", "U95", "STAT", "P"};
	public static final String[] LINEAR_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "SE", "L95", "U95", "STAT", "P"};
	public static final int[][] BUILD_36_PARS = {{}};
	public static final int[][] BUILD_37_PARS = {{}};
	
	public static final String[] FLAGS = new String[] {"duplicate", "parent-offspring", "sibling", "avuncular,gg", "first cousins", "second cousins"};
	public static final double[][] THRESHOLDS = new double[][] {{0, 0, 0, 0.95}, // duplicate
											{0, 0.90, 0, 0.49}, // parent-offspring
											{0, 0.30, 0.10, 0.35}, // sibling
											{0, 0.40, 0, 0}, // second degree, avuncular, gg
											{0, 0.20, 0, 0}, // third degree, cousins
											{0, 0.10, 0, 0}, // fourth degree 
	};


	public static void batchGenome(String root, int threads, double minPiHatToKeep) {
		PrintWriter writer;
		Hashtable<String,String> imiss;
		String[] inds, line;
		int step;
		Vector<String> extraCommands;
		Vector<String[]> v;
		String commands;
		Vector<String> filter;
        int count;

        extraCommands = new Vector<String>();
		if (!new File("plink.frq").exists()) {
			System.err.println("Error - 'plink.frq' required for batched genome");
			extraCommands.add("plink --bfile "+root+" --freq --noweb &");
		}

		if (!new File("plink.imiss").exists()) {
			System.err.println("Sorry, I can't trust you any more, I need to see a 'plink.imiss' file to proceed");
			extraCommands.add("plink --bfile "+root+" --missing --noweb &");
		}
		
		if (extraCommands.size() > 0) {
			System.out.println("Run the following first:");
			System.out.println(Array.toStr(Array.toStringArray(extraCommands), "\n"));
			return;
		}
		
		imiss = loadImissFile("plink.imiss");
		
		inds = HashVec.loadFileToStringArray(root+".fam", false, false, new int[] {0, 1, 5}, false);
		filter = new Vector<String>();
		count = 0;
		for (int i = 0; i<inds.length; i++) {
			line = inds[i].split("[\\s]+");
			if (!imiss.containsKey(line[0]+"\t"+line[1])) {
				System.err.println("Error - 'plink.imiss' does not contain individual"+line[0]+","+line[1]);
				return;
			} else if (imiss.get(line[0]+"\t"+line[1]).equals("1")) {
				count++;
			} else {
				if (Double.parseDouble(imiss.get(line[0]+"\t"+line[1])) > 0.50) {
					System.err.println("Warning - "+line[0]+","+line[1]+" is missing data for more than 50% of the markers");
				}
				filter.add(line[0]+"\t"+line[1]);
			}
        }
		if (count > 0) {
			System.err.println("Batch files will not include the "+count+" indiviudals that are not genotyped");
		}
		
		inds = Array.toStringArray(filter);
		step = (int)Math.ceil((double)inds.length/(double)threads);

		for (int i = 0; i<threads; i++) {
			try {
				writer = new PrintWriter(new FileWriter("tmp.list"+ext.formNum(i, 3)));
				for (int j = i*step; j<Math.min((i+1)*step, inds.length); j++) {
					writer.println(inds[j]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to file");
			}
		}

		v = new Vector<String[]>();
		for (int i = 0; i<threads; i++) {
			for (int j = i; j<threads; j++) {
				v.add(new String[] {ext.formNum(i, 3), ext.formNum(j, 3), i+"."+j});
			}
		}

		commands = "cd "+ext.pwd()+"\n~/bin/plink --noweb --bfile "+root+" --read-freq plink.frq --genome"+(minPiHatToKeep>0?" --min "+minPiHatToKeep:"")+" --genome-lists tmp.list[%0] tmp.list[%1] --out data.sub.[%2]";
		Files.qsub("genom", commands, Matrix.toStringArrays(v));
		try {
	        writer = new PrintWriter(new FileWriter("master.compile"));
	        writer.println("head -n1 data.sub.0.0.genome > header");
	        writer.println("cat data.sub*genome | fgrep -v FID1 | cat header - > cluster.genome");
	        writer.println("rm tmp.list*");
	        writer.println("rm data.sub.*");
	        writer.println("mkdir batches/");
	        writer.println("mv *.qsub* batches/");
	        writer.println();
	        writer.println("jcp gwas.Plink relate=cluster.genome");
	        writer.println("mkdir unrelateds");
	        writer.println("cd unrelateds/");
	        writer.println("plink --bfile ../../plink --keep ../cluster.genome_keep.dat --make-bed --noweb");
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+"master");
	        e.printStackTrace();
        }
	}

	public static Hashtable<String,String> loadImissFile(String filename) {
		BufferedReader reader;
		Hashtable<String,String> imiss;
		String[] line;

		imiss = new Hashtable<String,String>();
		try {
			reader = new BufferedReader(new FileReader(filename));
	        ext.checkHeader(reader.readLine().trim().split("[\\s]+"), IMISS_HEADER, true);
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	imiss.put(line[0]+"\t"+line[1], line[5]);
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			return null;
        } catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			return null;
        }
        
        return imiss;
	}

	public static void addToGenome(String root, int threads) {
		BufferedReader reader;
		PrintWriter writer;
		String[] inds, line;
		Vector<String[]> v;
		String commands;
		Vector<String> filter;
		Hashtable<String, String> hash; 
		int listCount, count, step, newCount;
		
		Hashtable<String,String> imiss;

		if (!new File("plink.frq").exists()) {
			System.err.println("Error - 'plink.frq' required for batched genome");
			return;
		}

		imiss = new Hashtable<String,String>();
		try {
	        reader = new BufferedReader(new FileReader("plink.imiss"));
	        ext.checkHeader(reader.readLine().trim().split("[\\s]+"), IMISS_HEADER, true);
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	imiss.put(line[0]+"\t"+line[1], line[5]);
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
			System.err.println("Sorry, I can't trust you any more, I need to see a 'plink.imiss' file to proceed");
			return;
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+"plink.imiss"+"\"");
	        return;
        }
		
		hash = new Hashtable<String, String>();
		listCount = 0;
		step = -1;
		while (new File("tmp.list"+ext.formNum(listCount, 3)).exists()) {
			System.out.println("tmp.list"+ext.formNum(listCount, 3)+" exists");
			try {
	            reader = new BufferedReader(new FileReader("tmp.list"+ext.formNum(listCount, 3)));
				count = 0;
	            while (reader.ready()) {
	            	line = reader.readLine().trim().split("[\\s]+");
	            	hash.put(line[0]+"\t"+line[1], "");
	            	count++;
	            }
	            step = Math.max(step, count);
	            reader.close();
            } catch (FileNotFoundException fnfe) {
	            System.err.println("Error: file \""+"tmp.list"+ext.formNum(listCount, 3)+"\" not found in current directory");
	            return;
            } catch (IOException ioe) {
	            System.err.println("Error reading file \""+"tmp.list"+ext.formNum(listCount, 3)+"\"");
	            return;
            }
            listCount++;
		}
		if (listCount == 0) {
			System.err.println("Error - cannot add to genome if there are no previous tmp.list### files to work from...");
			return;
		}
		
		inds = HashVec.loadFileToStringArray(root+".fam", false, false, new int[] {0, 1, 5}, false);
		filter = new Vector<String>();
		count = 0;
		for (int i = 0; i<inds.length; i++) {
			line = inds[i].split("[\\s]+");
			if (!imiss.containsKey(line[0]+"\t"+line[1])) {
				System.err.println("Error - 'plink.imiss' does not contain individual"+line[0]+","+line[1]);
				return;
			} else if (imiss.get(line[0]+"\t"+line[1]).equals("1")) {
				count++;
			} else if (!hash.containsKey(line[0]+"\t"+line[1])) {
				if (Double.parseDouble(imiss.get(line[0]+"\t"+line[1])) > 0.50) {
					System.err.println("Warning - "+line[0]+","+line[1]+" is missing data for more than 50% of the markers");
				}
				filter.add(line[0]+"\t"+line[1]);
			}
        }
		if (count > 0) {
			System.err.println("Batch files will not include the "+count+" indiviudals that are not genotyped");
		}
		inds = Array.toStringArray(filter);

		newCount = (int)Math.ceil((double)inds.length/(double)step);
		for (int i = 0; i<newCount; i++) {
			try {
				writer = new PrintWriter(new FileWriter("tmp.list"+ext.formNum(listCount+i, 3)));
				for (int j = i*step; j<Math.min((i+1)*step, inds.length); j++) {
					writer.println(inds[j]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to file");
			}
		}

		v = new Vector<String[]>();
		for (int i = 0; i<listCount+newCount; i++) {
			for (int j = 0; j<newCount; j++) {
				v.add(new String[] {ext.formNum(i, 3), ext.formNum(listCount+j, 3), i+"."+(listCount+j)});
			}
		}

		commands = "plink --noweb --bfile "+root+" --read-freq plink.frq --genome --genome-lists tmp.list[%0] tmp.list[%1] --out data.sub.[%2]";
		Files.batchIt("addGenome", "", threads, commands, Matrix.toStringArrays(v));
		try {
	        writer = new PrintWriter(new FileWriter("master", true));
	        writer.println();
	        writer.println("head -n1 data.sub.0.0.genome > header");
	        writer.println("cat data.sub*genome | fgrep -v FID1 | cat header - > cluster.genome");
	        writer.println("rm tmp.list*");
	        writer.println("rm data.sub.*");
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+"master");
	        e.printStackTrace();
        }
		System.out.println("Remember update the freq file!!");
		System.out.println("Reminder to either code HapMap samples as controls or use the -noFilter flag for unphenotyped (but still genotyped) samples");
		System.out.println("Also not a bad idea to run --missing to check for low-genotyped indiviudals that will crash the scripts");
	}
	
	public static void filterGenome(String genomeFile, String ids, boolean filterPairs) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Hashtable<String,String> hash = new Hashtable<String,String>();
        
        hash = HashVec.loadFileToHashString(ids, filterPairs?new int[] {0,1,2,3}:new int[] {0,1}, null, false, "\t", false, false, false);
        
        try {
	        reader = new BufferedReader(new FileReader(genomeFile));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(ids)+"_"+genomeFile));
	        writer.println(Array.toStr(reader.readLine().trim().split("[\\s]+")));
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if ((filterPairs && (hash.containsKey(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]) || hash.containsKey(line[2]+"\t"+line[3]+"\t"+line[0]+"\t"+line[1])))
	        			|| (!filterPairs && hash.containsKey(line[0]+"\t"+line[1]) && hash.containsKey(line[2]+"\t"+line[3]))) {
	        		writer.println(Array.toStr(line));
	        	}
	        }
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+genomeFile+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+genomeFile+"\"");
	        System.exit(2);
        }
	}
	
	public static void flagRelateds(String genomeFile, String famFile, String iMissFile, String lrrFile, String[] flags, double[][] thresholds, int removeOutTo) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, metrics, ids;
		Hashtable<String, String> callrates, lrr_sds, in, out, trav;
		Logger log;
		int numExtras;
		double[] probs;
		int[] order, ranks, counts;
		int sel;
		String temp;
		String genomeFileRoot;

		log = new Logger("flaggingRelateds.out", true);
		if (iMissFile == null) {
			log.report("Warning - no .imiss file specified, indiviudals will not be preferentially selected based on call rate");
			callrates = new Hashtable<String, String>();
		} else if (!Files.exists(iMissFile, false)) {
			log.reportError("Error - specified .imiss file missing , indiviudals will not be preferentially selected based on call rate");
			callrates = new Hashtable<String, String>();
		} else {
	        ext.checkHeader(Files.getHeaderOfFile(iMissFile, "[\\s]+", log), IMISS_HEADER, true);
			callrates = HashVec.loadFileToHashString(iMissFile, new int[] {0,1}, new int[] {5}, false, "\t", false, false, false);
			log.report("Successfully loaded '"+iMissFile+"'");
		}

		if (lrrFile == null) {
			log.report("Warning - no LRR_SD file specified, indiviudals will not be preferentially selected based on LRR_SD");
			lrr_sds = new Hashtable<String, String>();
		} else if (!Files.exists(lrrFile, false)) {
			log.reportError("Error - specified LRR_SD file ("+lrrFile+") is missing, indiviudals will not be preferentially selected based on LRR_SD");
			lrr_sds = new Hashtable<String, String>();
		} else {
	        ext.checkHeader(Files.getHeaderOfFile(lrrFile, "[\\s]+", log), new String[] {"FID", "IID", "LRR_SD"}, new int[] {0,1,2}, false, log, true);
			lrr_sds = HashVec.loadFileToHashString(lrrFile, new int[] {0,1}, new int[] {2}, false, "\t", false, false, false);
			log.report("Successfully loaded '"+lrrFile+"'");
		}
		
		if (famFile == null) {
			log.report("Warning - no .fam SD file specified, assuming .genome file is complete");
			in = new Hashtable<String, String>();
		} else if (!Files.exists(famFile, false)) {
			log.reportError("Error - specified .fam file ("+famFile+") is missing, assuming .genome file is complete");
			in = new Hashtable<String, String>();
			famFile = null;
		} else {
			in = HashVec.loadFileToHashString(famFile, new int[] {0,1}, new int[] {-7}, false, "\t", false, false, false);
			log.report("Successfully loaded '"+famFile+"'");
		}
		
		log.report("Started with "+in.size()+" known samples; will remove one of each pair of: "+Array.toStr(Array.subArray(flags, 0, removeOutTo), "/"));
		
		
		if (genomeFile.endsWith(".gz")) {
			genomeFileRoot = genomeFile.substring(0, genomeFile.length()-3);
		} else {
			genomeFileRoot = genomeFile;
		}
		
		
		counts = new int[removeOutTo];
		numExtras = 0;
		out = new Hashtable<String, String>();
		try {
			reader = Files.getAppropriateReader(genomeFile);
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CLUSTER_HEADER, true);
			log.report(ext.getTime()+"\tProcessing "+ext.removeDirectoryInfo(genomeFile)+"...");
			

			writer = new PrintWriter(new FileWriter(genomeFileRoot+"_relateds.xln"));
			writer.println("FID1\tIID1\tCALLRATE1\tLRR_SD1\tFID2\tIID2\tCALLRATE2\tLRR_SD2\tP(IBD=0)\tP(IBD=1)\tP(IBD=2)\tPI_HAT\tType\tFID_KEEP\tIID_KEEP\tFID_DROP\tIID_DROP");

			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				for (int i = 0; i < 2; i++) {
					if (!in.containsKey(line[i*2+0]+"\t"+line[i*2+1]) && !out.containsKey(line[i*2+0]+"\t"+line[i*2+1])) {
						if (famFile != null && numExtras < 20) {
							log.reportError("Warning - indiviudal '"+line[i*2+0]+"-"+line[i*2+1]+"' was in the .genome file, but not listed in the .fam file; adding...");
						}
						in.put(line[i*2+0]+"\t"+line[i*2+1], in.size()+"");
						numExtras++;
					}
				}
				probs = Array.toDoubleArray(Array.subArray(line, 6, 10));
				for (int i = 0; i < removeOutTo; i++) {
					if (probs != null && probs[0] >= thresholds[i][0] && probs[1] >= thresholds[i][1] && probs[2] >= thresholds[i][2] && probs[3] >= thresholds[i][3]) {
						line = Array.subArray(line, RELEVANT_CLUSTER_INDICES);
						counts[i]++;

						metrics = Array.stringArray(4, ".");
						for (int k = 0; k < 2; k++) {
							if (callrates.containsKey(line[k*2+0]+"\t"+line[k*2+1])) {
								metrics[k*2+0] = callrates.get(line[k*2+0]+"\t"+line[k*2+1]);
							}
							if (lrr_sds.containsKey(line[k*2+0]+"\t"+line[k*2+1])) {
								metrics[k*2+1] = lrr_sds.get(line[k*2+0]+"\t"+line[k*2+1]);
							}
						}
						writer.print(line[0]+"\t"+line[1]+"\t"+metrics[0]+"\t"+metrics[1]+"\t"+line[2]+"\t"+line[3]+"\t"+metrics[2]+"\t"+metrics[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+flags[i]);
						if (out.containsKey(line[0]+"\t"+line[1]) && out.containsKey(line[2]+"\t"+line[3])) {
							writer.println("\t.\t.\t.\t.");
						} else if (out.containsKey(line[0]+"\t"+line[1])) {
							writer.println("\t"+line[2]+"\t"+line[3]+"\t.\t.");
						} else if (out.containsKey(line[2]+"\t"+line[3])) {
							writer.println("\t"+line[0]+"\t"+line[1]+"\t.\t.");
						} else {
							if (metrics[0].equals(".") || metrics[2].equals(".") || (Double.parseDouble(metrics[0]) == Double.parseDouble(metrics[2]))) {
								if (metrics[1].equals(".") || metrics[3].equals(".") || (Double.parseDouble(metrics[1]) == Double.parseDouble(metrics[3]))) {
									sel = 0;
								} else if (Double.parseDouble(metrics[1]) < Double.parseDouble(metrics[3])){
									sel = 0;
								} else {
									sel = 1;
								}
							} else if (Double.parseDouble(metrics[0]) < Double.parseDouble(metrics[2])) {
								sel = 0;
							} else {
								sel = 1;
							}
							
							if (sel == 0) {
								writer.println("\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
								out.put(line[2]+"\t"+line[3], in.remove(line[2]+"\t"+line[3]));
							} else {
								writer.println("\t"+line[2]+"\t"+line[3]+"\t"+line[0]+"\t"+line[1]);
								out.put(line[0]+"\t"+line[1], in.remove(line[0]+"\t"+line[1]));
							}
						}

						probs = null;
					}
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + genomeFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + genomeFile + "\"");
			System.exit(2);
		}
		if (numExtras >= 20) {
			log.report("There were "+numExtras+" samples that were in the .genome file but not predefined in a .fam file");
		}
		
		temp = Files.getBakFilename(genomeFileRoot+"_relateds.xln", "", true);
		try {
			reader = new BufferedReader(new FileReader(temp));
			writer = new PrintWriter(new FileWriter(genomeFileRoot+"_relateds.xln"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (out.containsKey(line[13]+"\t"+line[14])) {
					line[13] = line[14] = ".";
				}
				writer.println(Array.toStr(line));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + temp + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + temp + "\"");
			System.exit(2);
		}
		new File(temp).renameTo(new File("trash"));
		
		log.report("Kept "+in.size()+" samples and dropped "+out.size()+" samples");
		log.report("Found the following number of each pairing:");
		for (int i = 0; i < counts.length; i++) {
			log.report(counts[i]+"\t"+flags[i]);
		}
		log.report("");
		for (int i = 0; i < 2; i++) {
			try {
				writer = new PrintWriter(new FileWriter(genomeFileRoot+"_"+(i==0?"keep":"drop")+".dat"));
				trav = i==0?in:out;
				ids = HashVec.getKeys(trav, false, false);
				ranks = new int[ids.length]; 
				for (int j = 0; j < ids.length; j++) {
					ranks[j] = Integer.parseInt(trav.get(ids[j]));
				}
				order = Sort.quicksort(ranks);
				for (int j = 0; j < ids.length; j++) {
					writer.println(ids[order[j]]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + genomeFileRoot+"_"+(i==0?"keep":"drop")+".dat");
				e.printStackTrace();
			}
		}
		
		if (Files.exists("ofInterest.txt", false)) {
			trav = HashVec.loadFileToHashString("ofInterest.txt", new int[] {0,1}, new int[] {-7}, false, "\t", false, false, false);
			log.report("Found file 'ofInterest.txt' with "+trav.size()+" records");
			ids = HashVec.getKeys(trav, false, false);
			counts = new int[3]; 
			for (int i = 0; i < ids.length; i++) {
				if (in.containsKey(ids[i])) {
					counts[0]++;
				} else if (out.containsKey(ids[i])) {
					counts[1]++;
				} else {
					counts[2]++;
				}
			}
			log.report("...of which "+counts[0]+" remain in "+counts[1]+" have been removed"+(counts[2]>0?"( AND "+counts[2]+" ARE UNACCOUNTED FOR!)":""));
		} else {
			log.report("Failed to locate file 'ofInterest.txt'");
		}
		log.report(ext.getTime()+"\tDone!");
	}
	
	// pair with a .log file if it exists to check numRepsEach and replace if need be
	// use hashtable instead so that mperms with only the most significant markers can get their extra reps (requires numRepsTotal to be an array)
	public static void collapseMperms(String patternAndCount) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;
		String pattern;
		int numReps, numRepsEach, numRepsTotal, rep;
		boolean done;
		String filename;
		int[][] counts;
		String[] markerNames;
		
		line = patternAndCount.split(",");
		pattern = line[0];
		numRepsEach = Integer.parseInt(line[1]);
		rep = 0;
		numRepsTotal = 0;
		done = false;
		markerNames = null;
		counts = null;
		do {
			rep++;
			filename = ext.insertNumbers(pattern, rep);
			if (new File(filename).exists()) {
				if (rep == 1) {
					markerNames = new String[Files.countLines(filename, true)];
					counts = new int[markerNames.length][2]; 
				}
				numReps = numRepsEach;
				// pair with a .log file if it exists to check numRepsEach and replace if need be
				try {
					reader = new BufferedReader(new FileReader(filename));
					line = reader.readLine().trim().split("[\\s]+");
					ext.checkHeader(line, MPERM_HEADER, true);
					count = 0;
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (rep==1) {
							markerNames[count] = line[1];
						} else if (!line[1].equals(markerNames[count])) {
							System.err.println("Error - marker mismatch in "+filename+"; expecting "+markerNames[count]+" at line "+(count+2)+", found "+line[1]);
						}
						counts[count][0] += (int)(Double.parseDouble(line[2])*(numReps+1))-1;
						counts[count][1] += (int)(Double.parseDouble(line[3])*(numReps+1))-1;
						
						count++;
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + filename + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + filename + "\"");
					System.exit(2);
				}

				numRepsTotal += numReps;
			} else {
				done = true;
			}
			
		} while (!done);
		System.out.println("Parsed the results of "+(rep-1)+" replicates, summing over a total of "+numRepsTotal+" permutations");
		if (rep > 0) {
			try {
				filename = ext.replaceAllWith(pattern, "#", "_ALL");
				writer = new PrintWriter(new FileWriter(filename));
				writer.println("SNP\tEMP1\tEMP2");
				for (int i = 0; i < markerNames.length; i++) {
					writer.println(markerNames[i]+"\t"+((double)counts[i][0]+1)/((double)numRepsTotal+1)+"\t"+((double)counts[i][1]+1)/((double)numRepsTotal+1));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + filename);
				e.printStackTrace();
			}
		}
	}
	
	public static void parseDiffMode6(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, int[]> indHash, markerHash;
		int count;
		int index;
		int[] counts;
		String[] keys;
		
		count = 0;
		indHash = new Hashtable<String, int[]>();
		markerHash = new Hashtable<String, int[]>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+".xln"));
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), new String[] {"SNP", "FID", "IID", "NEW", "OLD"}, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				count++;
				if (count % 1000000 == 0) {
					System.out.println(count);
				}
				if (!line[1].endsWith("-1") && !line[1].endsWith("-2")) {
					if (line[3].equals("0/0")) {
						index = 0;
					} else if (line[4].equals("0/0")) {
						index = 1;
					} else {
						index = 2;
					}
					if (indHash.containsKey(line[1]+"\t"+line[2])) {
						counts = indHash.get(line[1]+"\t"+line[2]);
					} else {
						indHash.put(line[1]+"\t"+line[2], counts = new int[3]);
					}
					counts[index]++;
					if (markerHash.containsKey(line[0])) {
						counts = markerHash.get(line[0]);
					} else {
						markerHash.put(line[0], counts = new int[3]);
					}
					counts[index]++;
				}
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
		
		try {
			writer = new PrintWriter(new FileWriter(filename+"_ind.xln"));
			keys = HashVec.getKeys(indHash);
			for (int i = 0; i < keys.length; i++) {
				writer.println(keys[i]+"\t"+Array.toStr(indHash.get(keys[i])));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename+"_ind.xln");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(filename+"_markers.xln"));
			keys = HashVec.getKeys(markerHash);
			for (int i = 0; i < keys.length; i++) {
				writer.println(keys[i]+"\t"+Array.toStr(markerHash.get(keys[i])));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename+"markers.xln");
			e.printStackTrace();
		}
	}
	
	public static void shiftPAR(String filename, int build) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int[][] pars;
		
		
		if (build == 36 || build == 18) {
			pars = BUILD_36_PARS;
		} else {
			if (build != 37 && build != 19) {
				System.err.println("Warning - '"+build+"' is an unknown build of the genome; using default build 37 (hg19)");
			}
			pars = BUILD_37_PARS;
		}

		
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.addToRoot(filename, "_shifterPAR")));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[0].equals("23")) {
					System.out.println(pars); // TODO unfinished business
				}
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
	}
	
	public static void batchLD(String root, double minR2ToKeep) {
		Files.qsub("runLD", 1, 22, "cd "+ext.pwd()+"\n~/bin/plink --noweb --bfile "+root+" --r2 --chr # --ld-window 99999 --ld-window-r2 "+minR2ToKeep+" --out chr#", 2000, 2);
	}
	
	public static void parseHitsForAllFilesInDirectory(String dir, double pvalThreshold, Logger log) {
		String[] models = null;
		String root;
		Hits hits;
		String[] hitList, params;

//		args = new String[] {null, null, "'MarkerName'", "'Chr'", "'Position'", "'P-value'", "!'P-value'<0.001", "tab", "replace=."};
//		GenParser.parse(args, log);
		new File(dir + "topHits.xln").delete();

		models = Files.list(dir, ".assoc.linear", false);
		
		hits = new Hits();
		params = new String[models.length];
		for (int i=0; i<models.length; i++) {
			root = models[i].substring(0, models[i].lastIndexOf(".assoc"));
			params[i] = dir+models[i]+" 'SNP' 'P'="+root+"_pval !'TEST'=ADD";
			log.report(ext.getTime()+"\t Incorporating data from "+root);
			hits.incorporateFromFile(dir+models[i], pvalThreshold, log);
		}
		hits.writeHits(dir+"hitList.dat");
		hitList = HashVec.loadFileToStringArray(dir+"hitList.dat", false, new int[] {0}, false);
		Files.combine(hitList, params, null, "MarkerName", ".", dir+"topHits.xln", log, true, true, false);
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		int genom = 0;
		int addGenom = -1;
		String root = "../plink";
		int count;
		String[] genomeID_files = null;
		String genomeFile = null;
		String famFile = "plink.fam";
		String iMissFile = "plink.imiss";
		String lrrFile = "lrr_sd_with_FID_IID.xln";
		int removeOutTo = 4;
		boolean filterPairs = false;
		String mperm = null;
		String diff = null;
		String shift = null;
		int build = 37;
		double minPiHatToKeep = -1;
		double minR2ToKeep = 0.20;
		boolean batchLD = false;
		
		String usage = "\n"+
		"gwas.Plink requires 0-1 arguments\n"+
		"   (1) collapse parallelized .mperm files with the given pattern and the specified number of reps per file (i.e. mperm=perm#.assoc.mperm,100000 (not the default))\n"+
		"  OR\n"+
		"   (1) number of parts to break into --genome on (i.e. genome=8 (not the default))\n"+
		"          type \"-remind\" for a reminder of how many jobs the number of parts leads to\n"+
		"   (2) name of PLINK root for --genome runs (i.e. root="+root+" (default))\n"+
		"   (3) minimum pi_hat to retain in --genome runs (i.e. minPiHatToKeep=0.1 (not the default))\n"+
		"   Note: requires both plink.frq and plink.imiss\n"+

		"  OR\n"+
		"   (1) add more indiviudals to a .genome using N threads (i.e. addGenome=4 (not the default))\n"+
		"          type \"-remind\" for a reminder of how many jobs the number of parts leads to\n"+
		"   (2) name of PLINK root for --genome runs (i.e. root="+root+" (default))\n"+
		"  OR\n"+
		"   (1) filter genome file by ids (i.e. filter=plink.genome,ids.txt (not the default))\n"+
		"   (2) filter pairs of ids (i.e. -filterPairs (not the default))\n"+
		"  OR\n"+
		"   (1) delineate related indiviudals and select one to keep (i.e. relate=plink.genome.gz (not the default))\n"+
		"   (2) (optional) if genome is incomplete (e.g. truncated at PI_HAT>0.10) (i.e. fam="+famFile+" (default))\n"+
		"   (3) (optional) select sample to keep based on call rate (i.e. imiss="+iMissFile+" (default))\n"+
		"   (4) (optional) select on LRR_SD if call rates are within 0.2% (i.e. lrr_sd="+lrrFile+" (default; needs to have format FID\tIID\tLRR_SD))\n"+
		"   (5) (optional) only purge out to certain level (i.e. level="+removeOutTo+" (default)), where:\n"+
		"  OR\n"+
		"   (1) summarize diff mode 6 results into markers and inds (i.e. diff=comp.diff (not the default))\n"+
		"  OR\n"+
		"   (1) shift markers in the PAR regions to chr25 (i.e. shift=plink.bim (not the default))\n"+
		"   (2) genome build (i.e. build="+build+" (default))\n"+
		"  OR\n"+
		"   (1) batch pairwise LD computation of r^2 (i.e. -batchLD (not the default))\n"+
		"   (2) name of PLINK root for --genome runs (i.e. root="+root+" (default))\n"+
		"   (3) minimum r^2 value to retain in the output to save on space (i.e. minR2="+minR2ToKeep+" (default; set to zero to keep all))\n"+
		"";
		for (int i = 0; i < FLAGS.length; i++) {
			usage += "           level "+(i+1)+": "+FLAGS[i]+"\tP(IBD=0)>="+THRESHOLDS[i][0]+"\tP(IBD=1)>="+THRESHOLDS[i][1]+"\tP(IBD=2)>="+THRESHOLDS[i][2]+"\tPI_HAT>="+THRESHOLDS[i][3]+"\n";
		}

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("genome=")) {
				genom = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("minPiHatToKeep=")) {
				minPiHatToKeep = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("addGenome=")) {
				addGenom = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("root=")) {
				root = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-remind")) {
				for (int threads = 1; threads<=20; threads++) {
					count = 0;
					for (int a = 0; a<threads; a++) {
						for (int b = a; b<threads; b++) {
							count++;
						}
					}
					System.out.println(threads+": "+count);
		        }
				return;
			} else if (args[i].startsWith("filter=")) {
				genomeID_files = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("-filterPairs")) {
				filterPairs = true;
				numArgs--;
			} else if (args[i].startsWith("relate=")) {
				genomeFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("fam=")) {
				famFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("imiss=")) {
				iMissFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("lrr_sd=")) {
				lrrFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("level=")) {
				removeOutTo = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("mperm=")) {
				mperm = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("diff=")) {
				diff = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("shift=")) {
				shift = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("build=")) {
				build = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-batchLD")) {
				batchLD = true;
				numArgs--;
			} else if (args[i].startsWith("minR2=")) {
				minR2ToKeep = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: "+args[i]);
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

//		genomeFile = "plink.genome";
//		iMissFile = "plink.imiss.xln";
//		removeOutTo = 5;
		
//		mperm = "perm#.assoc.mperm,100000";
		
//		genomeID_files = new String[] {"cluster.genome"};
		
		parseHitsForAllFilesInDirectory("D:/Visualization/ARIC_primary/gwas/QT_interval/", 0.0001, new Logger());
		System.exit(1);
		
		try {
			if (addGenom>0) {
				addToGenome(root, genom);
			} else if (genom>0) {
				batchGenome(root, genom, minPiHatToKeep);
			} else if (genomeID_files != null) {
				filterGenome(genomeID_files[0], genomeID_files[1], filterPairs);
			} else if (genomeFile != null) {
				flagRelateds(genomeFile, famFile, iMissFile, lrrFile, FLAGS, THRESHOLDS, removeOutTo);
			} else if (mperm != null) {
				collapseMperms(mperm);
			} else if (diff != null) {
				parseDiffMode6(diff);
			} else if (shift != null) {
				shiftPAR(shift, build);
			} else if (batchLD) {
				batchLD(root, minR2ToKeep);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
