package org.genvisis.link;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class Phenotype {
	public static final String[] OTHER_FILES = {"map##.dat", "chr##.dat", "chr##.map", "chr##.freq", "run#.qsub", "chr#_vc.qsub"};
	
	private static void update(String dir, String pheno, String pattern, int index, String missingValue, Logger log) {
		BufferedReader reader;
		PrintWriter[] writers;
		String[] line, phenoNames, phenos;
		String trav;
		Hashtable<String, String> hash;
		String filename;
		IntVector skips;
		
		dir = ext.verifyDirFormat(dir);
		phenoNames = Files.getHeaderOfFile(dir+pheno, log);
		ext.checkHeader(phenoNames, new String[] {"FID", "IID"}, new int[] {0,1}, false, log, false);
		System.out.println("Loading "+pheno);
		hash = HashVec.loadFileToHashString(dir+pheno, new int[] {0,1}, Array.subArray(Array.intArray(phenoNames.length), 2), pheno.endsWith(".csv"), "\t", true, false, false);
		phenoNames = Array.subArray(phenoNames, 2);
		writers = new PrintWriter[phenoNames.length];
		for (int i = 0; i < writers.length; i++) {
			new File(dir+phenoNames[i]+"/").mkdirs();
		}
		
		skips = new IntVector();
		for (int chr = 1; chr <= 23; chr++) {
			filename = ext.insertNumbers(pattern, chr);
			if (Files.exists(dir+filename)) {
				try {
					reader = new BufferedReader(new FileReader(dir+filename));
					for (int i = 0; i < writers.length; i++) {
						writers[i] = new PrintWriter(new FileWriter(dir+phenoNames[i]+"/"+filename));
					}
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						trav = hash.get(line[0]+"\t"+line[1]);
						if (trav == null) {
							phenos = Array.stringArray(phenoNames.length, missingValue);
						} else {
							phenos = trav.split("[\\s]+");
						}
						for (int i = 0; i < writers.length; i++) {
							line[index] = phenos[i];
							writers[i].println(Array.toStr(line));
						}
						
					}
					reader.close();
					for (int i = 0; i < writers.length; i++) {
						writers[i].close();
						for (int j = 0; j < OTHER_FILES.length; j++) {
							filename = ext.insertNumbers(OTHER_FILES[j], chr);
							if (Files.exists(dir+filename)) {
								Files.copyFile(dir+filename, dir+phenoNames[i]+"/"+filename);
							}
						}
					}
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+filename + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+filename + "\"");
					System.exit(2);
				}
				log.report("Updating chr"+chr);
				
			} else {
				skips.add(chr);
			}
		}
		
		if (skips.size() > 0) {
			log.report("skipped chromosomes: "+ext.listRanges(skips.toArray()));
		}
		
		
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "";
		String pheno = "pheno.dat";
		String pattern = "re_chrom##.pre";
		int index = 5;
		String missingValue = "x";
		String logfile = null;
		Logger log;

		String usage = "\n" + 
		"link.Phenotype requires 0-1 arguments\n" + 
		"   (1) directory (i.e. dir=" + dir + " (default))\n" + 
		"   (2) phenotype filename (i.e. pheno=" + pheno + " (default))\n" + 
		"   (3) pattern of files to update (i.e. pattern=" + pattern + " (default))\n" + 
		"   (4) index of column to swap out (i.e. index=" + index + " (default))\n" +
		"   (5) missing value (i.e. missingValue=" + missingValue + " (default))\n" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				pheno = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("pattern=")) {
				pattern = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("index=")) {
				index = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("missingValue=")) {
				missingValue = ext.parseStringArg(args[i], "x");
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
		
		dir = "D:/BOSS/Linkage/PCA_hits/";
		dir = "D:/BOSS/Linkage/PCA_all_files/";
		
		try {
			log = new Logger(logfile);
			update(dir, pheno, pattern, index, missingValue, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
