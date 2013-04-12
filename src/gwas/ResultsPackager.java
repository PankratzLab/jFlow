package gwas;

import java.io.*;
import java.util.*;

import bioinformatics.Alleles;
import common.*;

public class ResultsPackager {
	public static final String[] IBC_OUTPUT_FORMAT1 = {"CHR", "POS", "SNP", "STRAND (Illumina)", "STRAND (HapMap)", "N", "EFFECT_ALLELE1", "NON_EFFECT_ALLELE", "EA_FREQ", "BETA", "SE", "P_VAL"};
	public static final String[] TRADITIONAL_OUTPUT_FORMAT = {"Chr", "Position", "MarkerName", "Strand", "HapMapStrand", "N", "Effect_allele", "Reference_allele", "Freq1", "BETA", "SE", "P-value"};
	public static final String[] STANDARD_OUTPUT_FORMAT = {"MarkerName", "Chr", "Position", "Effect_allele", "Reference_allele", "Effect_allele_frequency", "N", "BETA", "SE", "P-value"};

	public static final String[] PLINK_REQS = {"SNP", "A1", "TEST", "NMISS", "OR", "BETA", "SE", "P"};
	
	public static void parseIBCFormatFromGWAF(String dir, String resultsFile, String mapFile, String originalFrqFile, String customFrqFile, String markersToReport, double filter, String outfile, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, String> markerHash, mapHash, originalFreqHash, customFreqHash;
		String delimiter;
		String[] alleles;
		String freq;
		
		if (outfile == null) {
			outfile = ext.rootOf(resultsFile, false)+"_out.csv";
		}
		
		if (markersToReport != null) {
			markerHash = HashVec.loadFileToHashNull(dir+markersToReport, false);
		} else {
			markerHash = null;
		}
		
		mapHash = HashVec.loadFileToHashString(dir+mapFile, new int[] {1}, new int[] {0,3}, false, "\t", false, false, false);
		originalFreqHash = HashVec.loadFileToHashString(dir+originalFrqFile, new int[] {1}, new int[] {2,3}, false, "\t", false, false, false); // add 4 if you want global frequency
		
		if (customFrqFile != null) {
			System.err.println("Warning - use of custom freq file has not been tested properly; if it works then remove this warning");
			customFreqHash = HashVec.loadFileToHashString(originalFrqFile, new int[] {1}, new int[] {2,3,4}, false, "\t", false, false, false); // add 4 if you want global frequency instead of custom Freq
		} else {
			customFreqHash = null;
		}
		
		try {
			reader = Files.getAppropriateReader(dir+resultsFile);
			writer = Files.getAppropriateWriter(dir+outfile);
			temp = reader.readLine().trim();
			delimiter = ext.determineDelimiter(temp);
			line = temp.split(delimiter);
			ext.checkHeader(line, GWAF.HEADER_SUMMARY, true);
//			writer.println(Array.toStr(IBC_OUTPUT_FORMAT));
			writer.println(Array.toStr(TRADITIONAL_OUTPUT_FORMAT));
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				trav = line[0];
				if ((markerHash == null || markerHash.containsKey(trav)) && !line[3].equals("") && (filter >= 1 || (!ext.isMissingValue(line[3]) && Double.parseDouble(line[3]) <= filter))) {
					if (mapHash.containsKey(trav)) {
						writer.print(mapHash.get(trav));
					} else if (mapHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
						writer.print(mapHash.get(ext.replaceAllWith(trav, ".", "-")));
					} else {
						log.reportError("Error - no map position for "+trav);
						writer.print(".\t.");
					}
					writer.print("\t"+trav+"\tNA\t+\t"+line[4]);
					
					if (originalFreqHash.containsKey(trav)) {
						alleles = originalFreqHash.get(trav).split("[\\s]+");
					} else if (originalFreqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
						alleles = originalFreqHash.get(ext.replaceAllWith(trav, ".", "-")).split("[\\s]+");
					} else {
						alleles = new String[] {".", "."};
						log.reportError("Error - no alleles from original .frq file for "+trav);
					}
					writer.print("\t"+Array.toStr(alleles));

					freq = line[5];
					if (customFreqHash != null) {
						if (customFreqHash.containsKey(trav)) {
							freq = Alleles.getAlleleFreqForA1(alleles[0], customFreqHash.get(trav).split("\t"))[2];
						} else if (customFreqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
							freq = Alleles.getAlleleFreqForA1(alleles[0], customFreqHash.get(ext.replaceAllWith(trav, ".", "-")).split("\t"))[2];
						} else {
							log.reportError("Error - no alleles from custom .frq file for "+trav);
							freq = ".";
						}
					}
					writer.println("\t"+freq+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + resultsFile + "\"");
			return;
		}
	}
	
	public static void parseStdFormatFromPlink(String dir, String resultsFile, String test, String mapFile, String freqFile, String markersToReport, double filter, String outfile, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, String> markerHash, mapHash, freqHash; // , customFreqHash;
		String delimiter;
		int[] indices;
		boolean logistic;
				
		if (outfile == null) {
			outfile = ext.rootOf(resultsFile, false)+".out";
		}
		
		if (markersToReport != null) {
			markerHash = HashVec.loadFileToHashNull(dir+markersToReport, false);
		} else {
			markerHash = null;
		}

		if (!Files.exists(dir+mapFile)) {
			log.reportError("Error: could not find map file '"+dir+mapFile+"'");
			return;
		}
		if (!Files.exists(dir+freqFile)) {
			log.reportError("Error: could not find freq file '"+dir+freqFile+"'");
			return;
		}
		if (!Files.exists(dir+resultsFile)) {
			log.reportError("Error: could not find results file '"+dir+mapFile+"'");
			return;
		}
		mapHash = HashVec.loadFileToHashString(dir+mapFile, new int[] {1}, new int[] {0,3}, false, "\t", false, false, false);
		freqHash = HashVec.loadFileToHashString(dir+freqFile, new int[] {1}, new int[] {2,3,4}, false, "\t", false, false, false); // 4 gets global frequency
		
		try {
			reader = Files.getAppropriateReader(dir+resultsFile);
			new File(ext.parseDirectoryOfFile(dir+outfile)).mkdirs();
			writer = Files.getAppropriateWriter(dir+outfile);
			temp = reader.readLine().trim();
			delimiter = ext.determineDelimiter(temp);
			line = temp.split(delimiter);
			indices = ext.indexFactors(PLINK_REQS, line, false, log, false, false);
			if (indices[4] == -1 && indices[5] == -1) {
				log.reportError("Error - results file did not contain a column for 'OR' (logistic) or 'BETA' (linear); aborting");
				return;
			} else if (indices[4] != -1 && indices[5] != -1) {
				log.reportError("Error - results file contain a column for both 'OR' (logistic) and 'BETA' (linear); aborting");
				return;
			} else if (indices[4] != -1) {
				logistic = true;
			} else {
				logistic = false;
			}

			if (indices[6] == -1) {
				log.reportError("Warning - results file did not contain a column for 'StdErr/SE'; values will be set to NA");
			}
			
			
			writer.println(Array.toStr(STANDARD_OUTPUT_FORMAT));
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				trav = line[indices[0]];
				if ((markerHash == null || markerHash.containsKey(trav)) && !line[3].equals("") && line[indices[2]].equalsIgnoreCase(test) && (filter >= 1 || (!ext.isMissingValue(line[indices[7]]) && Double.parseDouble(line[indices[7]]) <= filter))) {
					writer.print(trav); // MarkerName
					if (mapHash.containsKey(trav)) {
						writer.print("\t"+mapHash.get(trav)); // chr, pos
					} else if (mapHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
						writer.print("\t"+mapHash.get(ext.replaceAllWith(trav, ".", "-"))); // chr, pos
					} else {
						log.reportError("Error - no map position for "+trav);
						writer.print("\t.\t."); // null, null
					}
					if (line[indices[1]].equals("NA")) {
						writer.print("\t"+line[indices[1]]+"\t0\t0"); // missing data, null, null
					} else {
//try {
						if (freqHash.containsKey(trav)) {
							writer.print("\t"+Array.toStr(Alleles.getAlleleFreqForA1(line[indices[1]], freqHash.get(trav).split("\t")))); // a1, a2, a1_freq
						} else if (freqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
							writer.print("\t"+Array.toStr(Alleles.getAlleleFreqForA1(line[indices[1]], freqHash.get(ext.replaceAllWith(trav, ".", "-")).split("\t")))); // a1, a2, a1_freq
						} else {
							log.reportError("Error - no frequency for "+trav);
							writer.print("\t"+line[indices[1]]+"\t.\t."); // a1, null, null
						}
//} catch (Exception e) {
//	System.out.println("hola");
//}
					}
					writer.print("\t"+line[indices[3]]); // NMISS
					if (logistic) {
						if (line[indices[4]].equals("NA")) {
							writer.print("\t"+line[indices[4]]); // NA OR -> NA beta
						} else {
							writer.print("\t"+ext.formDeci(Math.log(Double.parseDouble(line[indices[4]])), 5, true)); // OR -> beta
						}
					} else {
						writer.print("\t"+line[indices[5]]); // beta
					}
//					writer.println("\t"+line[indices[6]]+"\t"+line[indices[7]]); // se, pval
					writer.println("\t"+(indices[6] == -1?"NA":line[indices[6]])+"\t"+line[indices[7]]); // se, pval
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + resultsFile + "\"");
			return;
		}
	}
	
	public static void parseBP() {
		String[] phenos = {"DBP", "HTN", "MAP", "PP", "SBP"};
		String[] types = {"M", "F"};
		String dir;
		
		dir = "C:/Users/npankrat/Documents/BloodPressure/";
		dir = "D:/BOSS/IBC_meta_analyses/BloodPressure/";
		for (int i = 0; i < types.length; i++) {
			for (int j = 0; j < phenos.length; j++) {
//				parseIBCFormatFromGWAF(dir, phenos[j]+"."+types[i]+".results.csv", "plink.map", "PLINK.FRQ", null, null, "BOSS-EHLS."+phenos[j]+"."+types[i]+".txt");
				parseIBCFormatFromGWAF(dir, phenos[j]+"."+types[i]+".results.csv", "plink.map", "plink.frq", null, "list.txt", 1, "BOSS-EHLS_"+phenos[j]+"_"+types[i]+"_072512.txt", new Logger());
				System.out.print(";BOSS-EHLS_"+phenos[j]+"_"+types[i]+"_072512.txt,11");
			}
		}
	}
	
	public static void createFromParameters(String filename, Logger log) {
        Vector<String> params;

		params = Files.parseControlFile(filename, "results", new String[] {"dir=", "results=plink.assoc.linear or plink.assoc.logistic", "type=plink or gwaf", "map=plink.map", "freq=plink.frq", "customFreq=", "list=specificSNPs.txt # leave blank for all", "filter=1.0 # set to 0.001 to report only those variants with a p<=0.001", "out=finalProduct.out"}, log);

		if (params != null) {
			params.add("log="+log.getFilename());
			main(Array.toStringArray(params));
		}
	}	

	public static void main(String[] args) {
		int numArgs = args.length;
		String resultsFile = null;
		String mapFile = "plink.bim";
		String freqFile = "plink.frq";
		String customFreqFile = null;
		String markersToReport = null;
		String outfile = null;
		String dir = "";
		String type = "plink";
		Logger log;
		String logfile = null;
		double filter = 1;

		String usage = "\n" +
		"gwas.ResultsPackager requires 0-1 arguments\n" + 
		"   (0) name of directory of all other files (i.e. dir="+dir+" (default))\n" + 
		"   (1) name of results file (i.e. results=all_results.csv (not the default))\n" + 
		"   (2) type of gwas file (i.e. type=plink (default) other options include =gwaf)\n" + 
		"   (3) name of map file (i.e. map="+mapFile+" (default))\n" + 
		"   (4) name of original freq file used to generate the gwaf files (i.e. freq="+mapFile+" (default; needs to be original freq file from when the gwaf files were made))\n" + 
		"   (5) (optional) name of freq file limited to those indiviudals used in anlayses (i.e. customFreq=femalesOnly.frq (not the default; only used in gwaf))\n" + 
		"   (6) (optional) list of markers to include (i.e. list=list.txt (not the default))\n" + 
		"   (7) (optional) name of output file (i.e. out=[results file]_out.csv (default; when I get around to coding it, it will be comma instead of tab delimited if ending in .csv))\n" + 
		"   (8) (optional) limit to those results with a p-value less than the specified filter (i.e. filter=0.001 (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("results=")) {
				resultsFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("type=")) {
				type = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("map=")) {
				mapFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("freq=")) {
				freqFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("customFreq=")) {
				customFreqFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("list=")) {
				markersToReport = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("filter=")) {
				filter = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
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
		try {
//			parseBP();
//			parseIBCFormatFromGWAF("D:/BOSS/GWAF/reclustered/PTA/", "pta_HmPCs_results.csv", "plink.map", "plink.frq", null, "list.txt", "PTA_results.txt");
//			parseIBCFormatFromGWAF("D:/CARe/CARe_geno_data_and_misc/IBC/FHS/iSELECT/gwaf/", "pta_fhs_parsedResults.csv", "grk100.map", "grk100.frq", null, "list.txt", "PTA_fhs_results.txt");

//			parseStdFormatFromPlink("D:/Myron/CALICO/T2DM/", "t2dm.assoc.logistic", "ADD", "plink.bim", "t2dm.frq", null, "cardia_page_t2dm_results.txt", new Logger());
//			parseStdFormatFromPlink("D:/Myron/CALICO/AgeMenarche/", "menarche1.assoc.linear", "ADD", "plink.bim", "menarche.frq", null, "cardia_page_menarche1_results.txt", new Logger());
//			parseStdFormatFromPlink("D:/Myron/CALICO/AgeMenarche/", "menarche2.assoc.linear", "ADD", "plink.bim", "menarche.frq", null, "cardia_page_menarche2_results.txt", new Logger());
//			System.exit(1);

			log = new Logger(logfile);
			if (type.equalsIgnoreCase("gwaf")) {
				parseIBCFormatFromGWAF(dir, resultsFile, mapFile, freqFile, customFreqFile, markersToReport, filter, outfile, log);
			} else if (type.equalsIgnoreCase("plink")) {
				parseStdFormatFromPlink(dir, resultsFile, "Add", mapFile, freqFile, markersToReport, filter, outfile, log);
			} else {
				System.err.println("Error - unknown results type: '"+type+"'");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
