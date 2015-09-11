package one;

import filesys.SnpMarkerSet;
import gwas.*;

import java.io.*;
import java.util.*;

import parse.GenParser;
import common.*;

public class Mega_Analyses {
//	public static String DIR = "C:\\mega\\";
	public static String DIR = "D:\\mega\\";
//	public static final String DIR = "";
	
	public static final String IMPUTATION_MAP = "/home/npankrat/NCBI/1000G/EUR.map";

	public static final String[][] STUDIES = {
		// abbrev, original filename, total sample size, lambda, lambda1000, coreStudy/needed for conditionals
		{"Ashkenazi", "AJ_PD_autosome_impute.txt", "446", "1.006", "1.028", "0", "1.028047962"},
		{"CHS", "CHS.PD.1000G.autosomes.11082011.txt", "3271", "1.009", "1.043", "0", "1.043478325"},
		{"deCODE", "DeCodeResults.LessDirty", "5520", "1.061", "1.057", "0", "1.05670092"},
//		{"deCODE", "deCODE_results.txt", "5520", "1.061", "1.057", "0"},
//		{"deCODE", "deCODE_results_wPositions.txt", "5520", "1.061", "1.057", "0", "1.05670092"},
		{"Dutch", "IPDGC.DutchResults.10312011.txt", "2763", "1.061", "1.056", "1", "1.056101112"},
		{"FHS", "FHS.PD.association.result.20111105.txt", "3949", "0.995", "0.958", "0", "0.957690495"},
		{"France", "PD_FR_RESULTS_MEGA_META.txt", "2969", "1.025", "1.032", "0", "0.88909397"},
		{"Germany", "IPDGC.GermanResults.10312011.txt", "1604", "1.025", "1.032", "1", "1.032081078"},
		{"HIHG", "hihg_hg19_autosomes.txt", "1193", "0.998", "0.997", "1", "0.996642331"},
		{"NGRC", "ngrc_hg19_autosomes.txt", "3938", "1.013", "1.007", "1", "1.006602624"},
		{"NIA", "IPDGC.NiaResults.10312011.txt", "2833", "1.035", "1.028", "1", "1.027906585"},
		{"PGPD", "progeni_genepd_hg19_autosomes.txt", "1680", "1.009", "1.011", "1", "1.010716473"},
		{"RSI", "pd05112011.RSI.dat", "5755", "0.984", "0.944", "0", "0.9437792"},
//		{"UK", "IPDGC.UKResults.txt", "6905", "1.034", "1.013", "1"},
		{"UK", "IPDGC.UKResults_duplicatesCollapsed.txt", "6905", "1.034", "1.013", "1", "1.013239905"},
//		{"23andMe", "23andMe.PD_meta.autosomal.20111026.txt", "66164", "1.212", "1.027", "1"},
		{"23andMe_v2", "23andMe.PD_meta.autosomal.v2.20120531.txt", "32760", "1.092", "1.016", "1", "1.015665477"},
		{"23andMe_v3", "23andMe.PD_meta.autosomal.v3.20120531.txt", "33404", "1.044", "1.026", "1", "1.026080290"}, // true
//		{"23andMe_v3", "23andMe.PD_meta.autosomal.v3.20120531.txt", "33404", "1.044", "1.026", "1", "1.025506832"}, // mike's first typo
//		{"23andMe_v3", "23andMe.PD_meta.autosomal.v3.20120531.txt", "33404", "1.044", "1.026", "1", "1.027393174"}, // mike's second typo
	};

//	public static final String[][] STUDIES = {
//		{"UK", "IPDGC.UKResults.txt", "6905", "0.987", "1"},
//	};
	
	public static void deCodedeCode() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		String dir, filename;
		
		dir = DIR+"deCODE_results\\";
		
		try {
			writer = new PrintWriter(new FileWriter(dir.substring(0, dir.length()-1)+".txt"));
			writer.println("MarkerName\tLocational_ID\tChr\tPosition\tEffect_allele\tReference_allele\tEAF\tGenotyped_or_imputed\tInformation\tBETA\tSE\tP\tNCASES\tNCONTROLS");
			for (int chr = 1; chr <= 22; chr++) {
				hash = new Hashtable<String, String>();
				for (int i = 0; i < 2; i++) {
					filename = dir + "chr"+chr+(i==0?"p":"q")+".assoc.log";
					if (new File(filename).exists()) {
						try {
							reader = new BufferedReader(new FileReader(filename));
							do {
								temp = reader.readLine();
							} while (!temp.startsWith("TRAIT"));
							ext.checkHeader(temp.trim().split("[\\s]+"), Mach.MACH2DAT_HEADER, true);
							while (reader.ready()) {
								line = reader.readLine().trim().split("[\\s]+");
								if (line[0].equals("Newton:")) {
									continue;
								}
								if (line[0].equals("")) {
									break;
								}
								if (i==0) {
									try {
										hash.put(line[1], line[4]);
									} catch (Exception e) {
										System.err.println("Error in '"+filename+"' : "+Array.toStr(line, "|"));
										e.printStackTrace();
									}
								} else if (hash.containsKey(line[1])) {
									if ((Double.parseDouble(scrub(hash.get(line[1])))) > Double.parseDouble(scrub(line[4]))) {
										hash.put(line[1], "first");
									} else {
										hash.put(line[1], "second");
									}
								} else {
									hash.put(line[1], "second");
								}							
							}
							reader.close();
						} catch (FileNotFoundException fnfe) {
							System.err.println("Error: file \"" + filename + "\" not found in current directory");
							System.exit(1);
						} catch (IOException ioe) {
							System.err.println("Error reading file \"" + filename + "\"");
							System.exit(2);
						}
					}
				}

				for (int i = 0; i < 2; i++) {
					filename = dir + "chr"+chr+(i==0?"p":"q")+".assoc.log";
					if (new File(filename).exists()) {
						try {
							reader = new BufferedReader(new FileReader(filename));
							do {
								temp = reader.readLine();
							} while (!temp.startsWith("TRAIT"));
							ext.checkHeader(temp.trim().split("[\\s]+"), Mach.MACH2DAT_HEADER, true);
							while (reader.ready()) {
								line = reader.readLine().trim().split("[\\s]+");
								if (line[0].equals("Newton:")) {
									continue;
								}
								if (line[0].equals("")) {
									break;
								}
								trav = hash.get(line[1]);
//								if (trav == null) {
//									System.err.println("Error in '"+filename+"' : how can this be the first time we see "+line[1]);
//								}
								if (i==1) {
									if (trav.equals("second")) {
										writer.println(translateDeCode(line));
									}
								} else if (trav.equals("first")) {
									writer.println(translateDeCode(line));
								} else if (!trav.equals("second")) {
									writer.println(translateDeCode(line));
								}
							}
							reader.close();
						} catch (FileNotFoundException fnfe) {
							System.err.println("Error: file \"" + filename + "\" not found in current directory");
							System.exit(1);
						} catch (IOException ioe) {
							System.err.println("Error reading file \"" + filename + "\"");
							System.exit(2);
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir.substring(0, dir.length()-1)+".txt");
			e.printStackTrace();
		}
	}
	
	public static void compare_deCode() {
		PrintWriter writer;
		String trav1, trav2;
		Hashtable<String, String> hash1, hash2;
		Vector<String> missingFrom2;
		String[] line, keys;
		
		System.out.println("Loading _deCODE_results.txt");
//		hash1 = HashVec.loadFileToHashString(DIR+"files/_deCODE_results.txt", new int[] {0}, new int[] {9, 10}, false, "\t", true, false, false);
		hash1 = HashVec.loadFileToHashString(DIR+"deCODE_results.txt", new int[] {0}, new int[] {9, 10, 8}, false, "\t", true, false, false);
		System.out.println("Loading DeCodeResults.LessDirty");
		hash2 = HashVec.loadFileToHashString(DIR+"00src/DeCodeResults.LessDirty", new int[] {0}, new int[] {5, 6, 4}, false, "\t", true, false, false); // 7, 8 two different p-values
		
		System.out.println("Comparing");
		missingFrom2 = new Vector<String>();
		keys = HashVec.getKeys(hash1, false, false);
		try {
			writer = new PrintWriter(new FileWriter("discordant.out"));
			writer.println("MikeBeta\tMikeSE\tNathanBeta\tNathanSE");
			for (int i = 0; i < keys.length; i++) {
				trav1 = hash1.get(keys[i]);
				line = trav1.trim().split("[\\s]+");
				for (int j = 0; j < line.length; j++) {
					while ((line[j].contains(".") && line[j].endsWith("0")) || line[j].endsWith(".")) {
						line[j] = line[j].substring(0, line[j].length()-1);
					}
					if (line[j].equals("-0")) {
						line[j] = "0";
					}
				}

				trav1 = Array.toStr(line);
				
				if (hash2.containsKey(keys[i])) {
					trav2 = hash2.remove(keys[i]);
					if (!trav1.substring(0, trav1.lastIndexOf("\t")).equals(trav2.substring(0, trav2.lastIndexOf("\t")))) {
						writer.println(keys[i]+"\t"+trav1+"\t"+trav2);
					}
				} else {
					missingFrom2.add(keys[i]);
					trav2 = null;
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "discordant.out");
			e.printStackTrace();
		}
		Files.writeList(Array.toStringArray(missingFrom2), "missingFrom2.out");
		Files.writeList(HashVec.getKeys(hash2, false, false), "missingFrom1.out");
		
		
	}
	
	
	public static void parseSuperSet() {
		String[] filenames;
		
		filenames = new String[STUDIES.length];
		for (int i = 0; i < STUDIES.length; i++) {
			filenames[i] = DIR+"files/"+STUDIES[i][1];
		}
		Unique.proc(filenames, null, null, DIR+"allVariants.txt", DIR+"variantCounts.out", true);
	}
	
	public static String scrub(String str) {
		if (str.equals(".-000")) {
			return ".000";
		}

		return str;
	}

	public static String translateDeCode(String[] line) {
		return line[1]+"\t.\t.\t.\t"+line[2].charAt(0)+"\t"+line[2].charAt(2)+"\t"+line[3]+"\t.\t"+line[4]+"\t"+line[5]+"\t"+line[7]+"\t"+line[9]+"\t"+line[12]+"\t"+line[13];
	}
	
	public static void parseAlleles() {
		PrintWriter writer;
		String filename;
		
		filename = DIR+"parseAlleles.crf";
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("lookup");
			writer.println("allAutosomalVariants.txt out=alleleFreqInput.xln lessMemoryButSlower keepIntermediateFiles");
			for (int i = 0; i < STUDIES.length; i++) {
				writer.println("files/"+STUDIES[i][1]+" 0 4="+STUDIES[i][0]+"_A1 5="+STUDIES[i][0]+"_A2 6="+STUDIES[i][0]+"_freq $#"+STUDIES[i][2]+"="+STUDIES[i][0]+"_N 8="+STUDIES[i][0]+"_Rsq");
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static void filterMarkers() {
		Hashtable<String,String> threeOrMoreAndRef, notInRef;
		PrintWriter writer;
		String[] allFiles, keys, metasoftParams;
		Logger log;
		String subDir;

		log = new Logger(DIR+"filtered/parse.log");
		new File(DIR+"filtered/").mkdirs();
		try {
			writer = new PrintWriter(new FileWriter(DIR+"filtered/input.txt"));
			writer.println("MARKER MarkerName");
			writer.println("ALLELE Allele1 Allele2");
			writer.println("EFFECT Effect");
			writer.println("STDERR StdErr");
			writer.println("SCHEME STDERR");
			writer.println("GENOMICCONTROL ON");
			writer.println("");
			new File(DIR+"filtered/").mkdirs();
			for (int i = 0; i < STUDIES.length; i++) {	
// prev			GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "4=Allele1", "5=Allele2", "9=Effect", "10=StdErr", "11=Pval", "!8>0.30", "!10!NA", "!10!.", "!10!0", "!6>-5", "!6<5", "!4!NA", "!6>0.001", "!6<0.999"}, log);
// ideal		GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "4=Allele1", "5=Allele2", "9=Effect", "10=StdErr", "11=Pval", "!8>=0.30", "!10!NA", "!10!.", "!10!0", "!9>=-5", "!9<=5", "!6>=0.001", "!6<=0.999", "6=MAF", "8=Rsq"}, log);
// nalls1		GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "4=Allele1", "5=Allele2", "9=Effect", "10=StdErr", "11=Pval", "!8>=0.30", "!9>=-5", "!9<=5", "!6>=0.001", "!6<=0.999", "6=MAF", "8=Rsq"}, log);
// finalIndices	GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "4=Allele1", "5=Allele2", "9=Effect", "10=StdErr", "11=Pval", "!8>=0.30", "!10!NA", "!10!.", "!10!0", "!9>=-5", "!9<=5", "!6>=0.001", "!6<=0.999", "6=MAF", "8=Rsq", "!2<23"}, log);
// ideal, chr<23
//				if (STUDIES[i][1].equals("deCODE")) {
//					GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "'MarkerName'=MarkerName", "'Effect_allele'=Allele1", "'Reference_allele'=Allele2", "'BETA'=Effect", "'SE'=StdErr", "'P'=Pval", "!'Information'>=0.30", "!'SE'!NA", "!'SE'!.", "!'SE'!0", "!'BETA'>=-5", "!'BETA'<=5", "!'EAF'>=0.001", "!'EAF'<=0.999", "'EAF'=MAF", "'Information'=Rsq"}, log);
//				} else {
//					GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "'MarkerName'=MarkerName", "'Effect_allele'=Allele1", "'Reference_allele'=Allele2", "'BETA'=Effect", "'SE'=StdErr", "'P'=Pval", "!'Information'>=0.30", "!'SE'!NA", "!'SE'!.", "!'SE'!0", "!'BETA'>=-5", "!'BETA'<=5", "!'EAF'>=0.001", "!'EAF'<=0.999", "'EAF'=MAF", "'Information'=Rsq", "!'Chr'<23"}, log);
//				}
// ideal, final, no chr<23
				GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "'MarkerName'=MarkerName", "'Effect_allele'=Allele1", "'Reference_allele'=Allele2", "'BETA'=Effect", "'SE'=StdErr", "'P'=Pval", "!'Information'>=0.30", "!'SE'!NA", "!'SE'!.", "!'SE'!0", "!'BETA'>=-5", "!'BETA'<=5", "!'EAF'>=0.001", "!'EAF'<=0.999", "'EAF'=MAF", "'Information'=Rsq"}, log);
//	pval only	GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "11="+STUDIES[i][0], "!8>0.30", "!10!NA", "!10!.", "!10!0", "!6>-5", "!6<5", "!4!NA"}, log);
				writer.println("PROCESS "+STUDIES[i][1]);
				log.report(STUDIES[i][0]+"\t"+Files.countLines(DIR+"filtered/"+STUDIES[i][1], 1)+"\t"+STUDIES[i][1]);
			}
			writer.println("");
			writer.println("");
			writer.println("OUTFILE META_ANALYSIS_beta_se .tbl");
			writer.println("ANALYZE");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + DIR+"filtered/input.txt");
			e.printStackTrace();
		}
	
		allFiles = new String[STUDIES.length];
		for (int i = 0; i < allFiles.length; i++) {
			allFiles[i] = DIR+"filtered/"+STUDIES[i][1];
		}
		Unique.proc(allFiles, Array.intArray(STUDIES.length, 1), null, DIR+"filtered/allSNPs.txt", DIR+"filtered/allSNPcounts.txt", true);
		
		GenParser.parse(new String[] {DIR+"filtered/allSNPcounts.txt", "skip=2", "1", "!2>=3", "out="+DIR+"filtered/SNPs_in3orMoreStudies.txt"}, log);
//		GenParser.parse(new String[] {DIR+"filtered/allSNPcounts.txt", "skip=2", "1", "!2>=5", "out="+DIR+"filtered/SNPs_in3orMoreStudies.txt"}, log);
		
//		threeOrMore = HashVec.loadFileToHashNull(DIR+"filtered/SNPs_in3orMoreStudies.txt", false);
//		threeOrMoreAndRef = HashVec.loadFileToHashNull(DIR+"filtered/SNPs_in3orMoreStudies.txt", false);
		threeOrMoreAndRef = HashVec.loadFileToHashString(DIR+"filtered/allSNPs.txt", false);
		notInRef = HashVec.loadFileToHashString(DIR+"non1000G_markers/markersToDelete.txt", false);
		
		keys = HashVec.getKeys(notInRef);
		for (int i = 0; i < keys.length; i++) {
			threeOrMoreAndRef.remove(keys[i]);
		}		
		
//		subDir = "InThreeOrMoreStudies/";
		subDir = "OnlyThoseInRef/";
		
		new File(DIR+"filtered/"+subDir).mkdirs();
		try {
			writer = new PrintWriter(new FileWriter(DIR+"filtered/"+subDir+"input.txt"));
			writer.println("MARKER MarkerName");
			writer.println("ALLELE Allele1 Allele2");
			writer.println("EFFECT Effect");
			writer.println("STDERR StdErr");
			writer.println("SCHEME STDERR");
			writer.println("GENOMICCONTROL ON");
			writer.println("");
			new File(DIR+"filtered/").mkdirs();
			for (int i = 0; i < STUDIES.length; i++) {
				FilterByLists.process(DIR+"filtered/"+STUDIES[i][1], new int[] {0}, DIR+"filtered/"+subDir+STUDIES[i][1], threeOrMoreAndRef, null, true, false, false, false, log);
				writer.println("PROCESS "+STUDIES[i][1]);
				log.report(STUDIES[i][0]+"\t"+Files.countLines(DIR+"filtered/"+subDir+STUDIES[i][1], 1)+"\t"+STUDIES[i][1]);
				Zip.gzip(DIR+"filtered/"+subDir+STUDIES[i][1]);
			}
			writer.println("");
			writer.println("");
			writer.println("OUTFILE META_ANALYSIS_beta_se .tbl");
			writer.println("ANALYZE HETEROGENEITY");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + DIR+"filtered/"+subDir+"input.txt");
			e.printStackTrace();
		}
		System.exit(1);
		
		allFiles = new String[STUDIES.length];
		metasoftParams = new String[STUDIES.length];
		for (int i = 0; i < allFiles.length; i++) {
			allFiles[i] = DIR+"filtered/"+subDir+STUDIES[i][1];
			metasoftParams[i] = DIR+"filtered/"+subDir+STUDIES[i][1]+" 'Effect' 'StdErr'";
		}
		log.report(ext.getTime()+"\tGenerating Metasoft input file");
//		Metasoft.generateInputFile(HashVec.getKeys(threeOrMoreAndRef), metasoftParams, DIR+"filtered/"+subDir+"metasoft.input", log);
		log.report(ext.getTime()+"\tCounting number with three or more");
		Unique.proc(allFiles, Array.intArray(STUDIES.length, 1), null, DIR+"filtered/"+subDir+"checkSNPs.txt", DIR+"filtered/"+subDir+"checkSNPcounts.txt", true);
		
		GenParser.parse(new String[] {DIR+"filtered/"+subDir+"checkSNPcounts.txt", "skip=2", "1", "!2>=3", "out="+DIR+"filtered/"+subDir+"doubleCheckedCounts.txt"}, log);
		
		
		Vector<String> v;
		v = new Vector<String>();
		for (int i = 0; i < STUDIES.length; i++) {
			if (STUDIES[i][5].equals("1")) {
				v.add(DIR+"filtered/"+subDir+""+STUDIES[i][1]);
			}
		}
		Unique.proc(Array.toStringArray(v), Array.intArray(v.size(), 1), null, DIR+"filtered/"+subDir+"coreStudySNPs.txt", DIR+"filtered/"+subDir+"coreStudyCounts.txt", true);
	}
	
	// note this search turned up empty
	public static void findNon1000Gmarkers() {
		Hashtable<String,String> in1000G;
		String[] allFiles, allFilteredFiles;
		Logger log;

		allFiles = new String[STUDIES.length];
		allFilteredFiles = new String[STUDIES.length];
		new File(DIR+"non1000G_markers/filtered/").mkdirs();
		log = new Logger(DIR+"non1000G_markers/parse.log");
		log.report("Loading reference markers from "+IMPUTATION_MAP);
		in1000G = HashVec.loadFileToHashString(IMPUTATION_MAP, 1, null, "", false);
		
		for (int i = 0; i < STUDIES.length; i++) {
//			System.out.println("Parsing files for "+STUDIES[i][0]);
			try {
				log.report(STUDIES[i][0], false, true);
				FilterByLists.process(DIR+"files/"+STUDIES[i][1], new int[] {0}, DIR+"non1000G_markers/"+STUDIES[i][0]+"_extras.txt", null, in1000G, true, false, false, false, log);
				log.report("\t"+Files.countLines(DIR+"non1000G_markers/"+STUDIES[i][0]+"_extras.txt", 1), false, true);
				FilterByLists.process(DIR+"filtered/"+STUDIES[i][1], new int[] {0}, DIR+"non1000G_markers/filtered/"+STUDIES[i][0]+"_extras.txt", null, in1000G, true, false, false, false, log);
				log.report("\t"+Files.countLines(DIR+"non1000G_markers/filtered/"+STUDIES[i][0]+"_extras.txt", 1), true, true);
				allFiles[i] = DIR+"non1000G_markers/"+STUDIES[i][0]+"_extras.txt";
				allFilteredFiles[i] = DIR+"non1000G_markers/filtered/"+STUDIES[i][0]+"_extras.txt";
			} catch (Exception e) {
				System.err.println("Error filtering files for to " + STUDIES[i][0]);
				e.printStackTrace();
			}
		}
		
		Unique.proc(allFiles, Array.intArray(STUDIES.length, 1), null, DIR+"non1000G_markers/extraSNPs.txt", DIR+"non1000G_markers/extraSNPcounts.txt", true);
		Unique.proc(allFilteredFiles, Array.intArray(STUDIES.length, 1), null, DIR+"non1000G_markers/filtered/extraSNPs.txt", DIR+"non1000G_markers/filtered/extraSNPcounts.txt", true);
		
		GenParser.parse(new String[] {DIR+"non1000G_markers/extraSNPcounts.txt", "skip=2", "1", "!2>=3", "out="+DIR+"non1000G_markers/relevantSNP.txt"}, log);
		GenParser.parse(new String[] {DIR+"non1000G_markers/filtered/extraSNPcounts.txt", "skip=2", "1", "!2>=3", "out="+DIR+"non1000G_markers/filtered/relevantSNP.txt"}, log);
		
	}
	
	public static void findDuplicateMarkers() {
		BufferedReader reader;
		PrintWriter writer;
		String[] allFiles, allFilteredFiles, allSNPs;
		Hashtable<String, Vector<String>> hash;
		Logger log;

		new File(DIR+"duplicatedMarkers/filtered/").mkdirs();
		log = new Logger(DIR+"duplicatedMarkers/parse.log");

		allFiles = new String[STUDIES.length];
		allFilteredFiles = new String[STUDIES.length];
		for (int i = 0; i < STUDIES.length; i++) {
//			System.out.println("Parsing files for "+STUDIES[i][0]);
			try {
//				Unique.proc(new String[] {DIR+"files/"+STUDIES[i][1]}, new int[] {1}, null, null, DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_SNPcounts.txt", true);
//				GenParser.parse(new String[] {DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_SNPcounts.txt", "skip=2", "1", "!2>1", "out="+DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedSNPs.txt"}, log);
//				FilterByLists.process(DIR+"files/"+STUDIES[i][1], DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedSNPs.txt", null, new int[] {0}, DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedPositions.txt", true, log);

//				Unique.proc(new String[] {DIR+"filtered/"+STUDIES[i][1]}, new int[] {1}, null, null, DIR+"duplicatedMarkers/filtered/"+STUDIES[i][0]+"_SNPcounts.txt", true);
//				GenParser.parse(new String[] {DIR+"duplicatedMarkers/filtered/"+STUDIES[i][0]+"_SNPcounts.txt", "skip=2", "1", "!2>1", "out="+DIR+"duplicatedMarkers/filtered/"+STUDIES[i][0]+"_duplicatedSNPs.txt"}, log);
//				
//				log.report(STUDIES[i][0]+"\t"+Files.countLines(DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedSNPs.txt", true)+"\t"+Files.countLines(DIR+"duplicatedMarkers/filtered/"+STUDIES[i][0]+"_duplicatedSNPs.txt", true));
				allFiles[i] = DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedSNPs.txt";
				allFilteredFiles[i] = DIR+"duplicatedMarkers/filtered/"+STUDIES[i][0]+"_duplicatedSNPs.txt";
			} catch (Exception e) {
				System.err.println("Error filtering files for to " + STUDIES[i][0]);
				e.printStackTrace();
			}
		}
//		
		Unique.proc(allFiles, Array.intArray(STUDIES.length, 1), null, DIR+"duplicatedMarkers/allDuplicatedSNPs.txt", DIR+"duplicatedMarkers/allDuplicatedSNPcounts.txt", true);
//		Unique.proc(allFilteredFiles, Array.intArray(STUDIES.length, 1), null, DIR+"duplicatedMarkers/filtered/allDuplicatedSNPs.txt", DIR+"duplicatedMarkers/filtered/allDuplicatedSNPcounts.txt", true);
//		
//		GenParser.parse(new String[] {DIR+"duplicatedMarkers/allDuplicatedSNPcounts.txt", "skip=2", "1", "!2>=3", "out="+DIR+"duplicatedMarkers/relevantSNPs.txt"}, log);
//		GenParser.parse(new String[] {DIR+"duplicatedMarkers/filtered/allDuplicatedSNPcounts.txt", "skip=2", "1", "!2>=3", "out="+DIR+"duplicatedMarkers/filtered/relevantSNPs.txt"}, log);

		String filename = DIR+"duplicatedMarkers/allDuplicatedSNPs.txt";
//		String filename = DIR+"duplicatedMarkers/filtered/allDuplicatedSNPcounts.txt";
		allSNPs = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
		hash = new Hashtable<String, Vector<String>>();
		for (int i = 0; i < allSNPs.length; i++) {
			hash.put(allSNPs[i], new Vector<String>());
		}
		
		String[] line;
		for (int i = 0; i < STUDIES.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedPositions.txt"));
				System.err.println("Parsing duplicates in "+STUDIES[i][0]);
				ext.checkHeader(reader.readLine().trim().split("[\\s]+"), new String[] {"MarkerName", "Chr", "Position", "Genotyped_or_imputed"}, new int[] {0,2,3,7}, false, log, false);
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (hash.containsKey(line[0])) {
						hash.get(line[0]).add(STUDIES[i][0]+"/"+line[2]+"/"+line[3]+"/"+line[7]);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedPositions.txt" + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + DIR+"duplicatedMarkers/"+STUDIES[i][0]+"_duplicatedPositions.txt" + "\"");
				System.exit(2);
			}
		}
		
		StringVector keys;
		ByteVector chrs;
		IntVector positions;
		Vector<String> v;
		byte chr;
		int initSize, pos;
		Hashtable<String, String> refChrHash;
		int[] order;
		
		keys = new StringVector(HashVec.getKeys(hash, false, false));
		chrs = new ByteVector();
		positions = new IntVector();
		initSize = keys.size();
		for (int i = 0; i < initSize; i++) {
			v = hash.get(keys.elementAt(i));
			chr = -1;
			pos = -1;
			for (int j = 0; j < v.size(); j++) {
				line = v.elementAt(j).split("/");
				if (chr == -1 && pos == -1) {
					chrs.add(chr = Byte.parseByte(line[1]));
					positions.add(pos = Integer.parseInt(line[2]));
//				} else if (Byte.parseByte(line[1]) != chr || Integer.parseInt(line[2]) != pos) {
//					keys.add(keys.elementAt(i)); 				// correct, adds at end
//					chrs.add(Byte.parseByte(line[1]));			// incorrect, adds from the beginning 
//					positions.add(Integer.parseInt(line[2]));	// incorrect, adds from the beginning
				}
			}
		}

		System.err.println("Loading reference map from: "+IMPUTATION_MAP);
		refChrHash = new SnpMarkerSet(IMPUTATION_MAP, SnpMarkerSet.PLINK_MAP_FORMAT_WITHOUT_CM, false, log).getChrHash();
		
		order = Sort.orderTwoLayers(chrs.toArray(), positions.toArray(), log);
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_positions.xln"));
			writer.println("MarkerName\tChr\tPosition\tAltLoc\trefChr\trefPos\trefAgree\tStudy/chr/pos");
			for (int i = 0; i < order.length; i++) {
				v = hash.get(keys.elementAt(order[i]));
				writer.println(keys.elementAt(order[i])+"\t"+chrs.elementAt(order[i])+"\t"+positions.elementAt(order[i])+"\t"+(order[i]<initSize?0:1)
						+"\t"+(refChrHash.containsKey(keys.elementAt(order[i]))?refChrHash.get(keys.elementAt(order[i]))+"\t"+(refChrHash.get(keys.elementAt(order[i])).split("[\\s]+")[1].equals(positions.elementAt(order[i])+"")?"1":"0"):".\t.\t0")
						+"\t"+Array.toStr(Array.toStringArray(v)));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename, false)+"_positions.xln");
			e.printStackTrace();
		}
	}

	public static void countLines(String dir) {
		String[] files;
		Logger log;
	
		log = new Logger(dir+"counts.log");
		files = Files.list(dir, ".txt", false);
		for (int i = 0; i < files.length; i++) {
			log.report(files[i]+"\t"+Files.countLines(dir+files[i], 0));
		}
	}
	
	public static void parseCRF() {
		PrintWriter writer;
		String filename;
		
//		filename = DIR+"parseHits.crf";
		filename = DIR+"parsePs.crf";
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("lookup");
//			writer.println("hits.txt out=hitsLookSee.xln");
			writer.println("allAutosomalVariants.txt out=allPvals.xln");
			for (int i = 0; i < STUDIES.length; i++) {
//				writer.println(STUDIES[i][1]+" 0 4=Allele1 5=Allele2 9=Beta 10=SE 11=Pvalue 6=Freq1 8=Rsq");
				writer.println(STUDIES[i][1]+" 0 11="+STUDIES[i][0]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
		
	}
	
	public static void keepMarkers() {
		Logger log;
		
		log = new Logger(DIR+"kept/");
		new File(DIR+"kept/").mkdirs();
		for (int i = 0; i < STUDIES.length; i++) {
			FilterByLists.process(DIR+"files/"+STUDIES[i][1], DIR+"duplicatedMarkers.txt", null, new int[] {0}, DIR+"kept/"+STUDIES[i][1], true, false, false, log);
		}
	}

	public static void listQQfiles() {
		for (int i = 0; i < STUDIES.length; i++) {
//			System.out.print((i==0?"":";")+"files/"+STUDIES[i][1]+",11");
			System.out.print((i==0?"":";")+"filtered/"+STUDIES[i][1]+",5");
		}
		System.out.println();
	}
	
	public static void purifyAllDatasets() {
		String[] header;
		int[] indices;
		Logger log;
		
		log = new Logger(DIR+"files/purification.log");		
		for (int i = 0; i < STUDIES.length; i++) {
			header = Files.getHeaderOfFile(DIR+"filtered/"+STUDIES[i][1], log);
			indices = ext.indexFactors(new String[][] {{"MarkerName"}, {"Information", "Rsq"}}, header, false, true, true, true);
//			purifyDuplicateMarkers(DIR+"filtered/"+STUDIES[i][1], indices[0], indices[1], false, log);
//			log.report("Parsing "+STUDIES[i][0]+" using "+header[indices[0]]+" and "+header[indices[1]]);
//			purifyDuplicateMarkers(DIR+"files/"+STUDIES[i][1], indices[0], indices[1], false, log);
			purifyDuplicateMarkers(DIR+"filtered/OnlyThoseInRef/"+STUDIES[i][1], indices[0], indices[1], false, log);
		}
	}
	
	private static void purifyDuplicateMarkers(String filename, int indexOfDuplicatedFactor, int indexOfFilter, boolean minimizeInsteadofMaximize, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		Hashtable<String, String> hash, duplicates;
		double mostExtremeTie;
		long time;
		
		time = new Date().getTime();
		mostExtremeTie = minimizeInsteadofMaximize?Double.POSITIVE_INFINITY:Double.NEGATIVE_INFINITY;
		hash = new Hashtable<String, String>();
		duplicates = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().trim().split("[\\s]+");
			log.report("Filtering "+filename+" using '"+line[indexOfDuplicatedFactor]+"' and '"+line[indexOfFilter]+"'");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!hash.containsKey(line[indexOfDuplicatedFactor])) {
					hash.put(line[indexOfDuplicatedFactor], line[indexOfFilter]);
				} else {
					if (line[indexOfFilter].equals(hash.get(line[indexOfDuplicatedFactor]))) {
						System.err.println(line[indexOfDuplicatedFactor]+"\t"+line[indexOfFilter]+"\tduplicates with identical filter factor");
						if (( minimizeInsteadofMaximize && Double.parseDouble(line[indexOfFilter])<mostExtremeTie )
								|| ( !minimizeInsteadofMaximize && Double.parseDouble(line[indexOfFilter])>mostExtremeTie )) {
							mostExtremeTie = Double.parseDouble(line[indexOfFilter]);
						}
					}
					if (( minimizeInsteadofMaximize && Double.parseDouble(line[indexOfFilter])<Double.parseDouble(hash.get(line[indexOfDuplicatedFactor])) )
							|| ( !minimizeInsteadofMaximize && Double.parseDouble(line[indexOfFilter])>Double.parseDouble(hash.get(line[indexOfDuplicatedFactor])) )) {
						hash.put(line[indexOfDuplicatedFactor], line[indexOfFilter]);
						duplicates.put(line[indexOfDuplicatedFactor], line[indexOfFilter]);
					} else {
						duplicates.put(line[indexOfDuplicatedFactor], hash.get(line[indexOfDuplicatedFactor]));
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		if (duplicates.size() == 0) {
			System.err.println("There were no duplicated factors in this file");
		} else if (duplicates.size() == 1) {
			System.err.println("There was only a single duplicated factor");
		} else {
			System.err.println("There were "+duplicates.size()+" duplicated factors");
		}
		if (!Double.isInfinite(mostExtremeTie)) {
			System.err.println("The value of the most extreme tie was: "+mostExtremeTie);
		}
		System.out.println("Finished scanning in " + ext.getTimeElapsed(time));
		
		if (duplicates.size() > 0) {
			try {
				reader = new BufferedReader(new FileReader(filename));
				new File(ext.parseDirectoryOfFile(filename)+"purified/").mkdir();
				writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename)+"purified/"+ext.removeDirectoryInfo(filename)));
				writer.println(reader.readLine());
				while (reader.ready()) {
					temp = reader.readLine();
					line = temp.trim().split("[\\s]+");
					if (duplicates.containsKey(line[indexOfDuplicatedFactor])) {
						if (line[indexOfFilter].equals(duplicates.get(line[indexOfDuplicatedFactor]))) {
							writer.println(temp);
							duplicates.put(line[indexOfDuplicatedFactor], "used");
						}
					} else {
						writer.println(temp);
					}
				}
				writer.close();
				reader.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + ext.parseDirectoryOfFile(filename)+"purified/"+ext.rootOf(filename, true));
				e.printStackTrace();
			}
			System.out.println("Start to finish took " + ext.getTimeElapsed(time));
		}
	}
	
//	public static void parseByStudy(String markerFile) {
//		BufferedReader reader;
//		PrintWriter writer;
//		String[] line;
//		String temp, trav;
//		Hashtable<String, String> hash = new Hashtable<String, String>();
//		Vector<String> v = new Vector<String>();
//		int count;
//		long time;
//		String[] markerNames;
//		
//		
//		
//		
//	}
//	
	public static void main(String[] args) {
		if (new File("/software/").exists()) {
			DIR = "";
		}
		
//		Hashtable<String,String> xChrome = HashVec.loadFileToHashNull("all_chrXmarkers_in1000G.txt", false);
//		FilterByLists.process("META_ANALYSIS_beta_se1.tbl", new int[] {0}, "META_ANALYSIS_beta_se_noX1.tbl", null, xChrome, true, false, false, new Logger());

//		Hashtable<String,String> threeOrMore = HashVec.loadFileToHashNull("SNPs_in3orMoreStudies.txt", false);
//		FilterByLists.process("META_ANALYSIS_beta_se_noX1.tbl", new int[] {0}, "META_ANALYSIS_beta_se_Final1.tbl", threeOrMore, null, true, false, false, new Logger());
//		System.exit(1);

//		Hashtable<String,String> threeOrMore = HashVec.loadFileToHashNull("SNPs_in3orMoreStudies.txt", false);
//		FilterByLists.process("META_ANALYSIS_beta_se1.tbl", new int[] {0}, "META_ANALYSIS_beta_se_Final1.tbl", threeOrMore, null, true, false, false, new Logger());
//		System.exit(1);

		Hashtable<String,String> threeOrMore = HashVec.loadFileToHashString("SNPs_in3orMoreStudies.txt", false);
		FilterByLists.process(DIR+"filtered/OnlyThoseInRef/MetasoftResults.txt", new int[] {0}, DIR+"filtered/OnlyThoseInRef/MetasoftResults_Final.out", threeOrMore, null, true, false, false, false, new Logger());
		System.exit(1);
		
		String subDir = "OnlyThoseInRef/";
		Vector<String> v;
		v = new Vector<String>();
		for (int i = 0; i < STUDIES.length; i++) {
			if (STUDIES[i][5].equals("1")) {
				v.add(DIR+"filtered/"+subDir+""+STUDIES[i][1]);
			}
		}
		Unique.proc(Array.toStringArray(v), Array.intArray(v.size(), 1), null, DIR+"filtered/"+subDir+"coreStudySNPs.txt", DIR+"filtered/"+subDir+"coreStudyCounts.txt", true);
		System.exit(1);
		
		String[] metasoftParams;
//		String dir; //, subDir;
		String[] line;
		double[] lambdas;


//		dir = "filtered/";
		subDir = "OnlyThoseInRef/";
		metasoftParams = new String[STUDIES.length];
		lambdas = new double[STUDIES.length];
		String subset;
		subset = null;
//		subset = "D:/mega/filtered_no23for23/OnlyThoseInRef/ToNathan.txt";
//		subset = "D:/mega/filtered_corruptDeCode/MegaJune5th/RE2/GcResultTopHits.txt";
		if (Math.random() < 2 || subset != null) {
			line = HashVec.loadFileToStringArray(subset, true, new int[] {0}, false);
		} else {
			line = HashVec.getKeys(HashVec.loadFileToHashSet(DIR+"filtered/allSNPs.txt", false), false, false);
		}
		
		for (int i = 0; i < metasoftParams.length; i++) {
			metasoftParams[i] = DIR+"filtered/"+subDir+STUDIES[i][1]+" 'MarkerName' 'Allele1' 'Allele2' 'Effect'="+STUDIES[i][1]+"_Beta;NA 'StdErr'="+STUDIES[i][1]+"_SE;NA fail";
			lambdas[i] = Math.sqrt(Double.parseDouble(STUDIES[i][6]));
//			System.out.print("\t"+STUDIES[i][0]+"_Beta\t"+STUDIES[i][0]+"_SE");
		}
 		Metasoft.generateRawFile(line, metasoftParams, DIR+"filtered/"+subDir+"metasoftFirstK.raw", new Logger());
 		Metasoft.processRawToInputFile(DIR+"filtered/"+subDir+"metasoftFirstK.raw", DIR+"filtered/"+subDir+"metasoftFirstK_flipped.raw", new Logger());
 		Metasoft.applyLambdas(DIR+"filtered/"+subDir+"metasoftFirstK_flipped.raw", DIR+"filtered/"+subDir+"metasoftFirstK_wGC.input", lambdas, new Logger());
		System.exit(1);

		
		deCodedeCode();
		compare_deCode();
		parseAlleles();
		parseCRF();
		parseSuperSet();
		filterMarkers();
		findNon1000Gmarkers();
		findDuplicateMarkers();
//		purifyDuplicateMarkers("/mega/duplicatedMarkers/UK_duplicatedPositions.txt", 0, 8, false); // test on a small set of data, seems to work
//		purifyDuplicateMarkers(DIR+"files/IPDGC.UKResults.txt", 0, 8, false); // apply to full dataset
//		purifyDuplicateMarkers(DIR+"files/23andMe.PD_meta.autosomal.v2.20120531.txt", 0, 8, false);
		purifyAllDatasets();
		countLines(DIR+"filtered_almostNalls_plusSEfilters/");
		countLines(DIR+"filtered_almostNalls_plusSEfilters/InThreeOrMoreStudies/");
		listQQfiles();
		keepMarkers();
//		parseByStudy(DIR+"filtered/"+subDir+"minHits.xln"); // used Metasoft file instead
	}
}
