package one;

import gwas.Mach;

import java.io.*;
import java.util.*;
import parse.GenParser;
import common.*;

public class Mega_Analyses {
	public static String DIR = "D:\\mega\\";
//	public static final String DIR = "";

	public static final String[][] STUDIES = {
		{"Ashkenazi", "AJ_PD_autosome_impute.txt", "446", "0.972"},
//		{"deCODE", "DeCodeResults.LessDirty.txt", "5520", "0.943"},
		{"deCODE", "deCODE_results.txt", "5520", "0.943"},
		{"Dutch", "IPDGC.DutchResults.10312011.txt", "2763", "0.944"},
		{"France", "PD_FR_RESULTS_MEGA_META.txt", "2969", "1.111"},
		{"Germany", "IPDGC.GermanResults.10312011.txt", "1604", "0.968"},
		{"HIHG", "hihg_hg19_autosomes.txt", "1193", "1.003"},
		{"NGRC", "ngrc_hg19_autosomes.txt", "3938", "0.993"},
		{"NIA", "IPDGC.NiaResults.10312011.txt", "2833", "0.972"},
		{"PGPD", "progeni_genepd_hg19_autosomes.txt", "1680", "0.989"},
		{"UK", "IPDGC.UKResults.txt", "6905", "0.987"},
		{"23andMe", "23andMe.PD_meta.autosomal.20111026.txt", "66164", "0.973"},
		{"CHS", "CHS.PD.1000G.autosomes.11082011.txt", "000", "0.000"},
		{"FHS", "FHS.PD.association.result.20111105.txt", "000", "0.000"},
		{"RSI", "pd05112011.RSI.dat", "000", "0.000"},
	};

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
	
	public static void parseSuperSet() {
		String[] filenames;
		
		filenames = new String[STUDIES.length];
		for (int i = 0; i < STUDIES.length; i++) {
			filenames[i] = DIR+"files/"+STUDIES[i][1];
		}
		Unique.proc(filenames, DIR+"allVariants.txt", DIR+"variantCounts.out", true);
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
			writer.println("hits");
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
		PrintWriter writer; 
		Logger log;

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
			log = new Logger(DIR+"filtered/parse.log");
			for (int i = 0; i < STUDIES.length; i++) {
				GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "4=Allele1", "5=Allele2", "9=Effect", "10=StdErr", "11=Pval", "!8>0.30", "!10!NA", "!10!.", "!10!0", "!6>-5", "!6<5", "!4!NA", "!6>0.001", "!6<0.999"}, log);
//				GenParser.parse(new String[] {DIR+"files/"+STUDIES[i][1], "out="+DIR+"filtered/"+STUDIES[i][1], ".-000=>.000", "0=MarkerName", "11="+STUDIES[i][0], "!8>0.30", "!10!NA", "!10!.", "!10!0", "!6>-5", "!6<5", "!4!NA"}, log);
				writer.println("PROCESS "+STUDIES[i][1]);
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
	}
	
	public static void countLines(String dir) {
		String[] files;
		
		files = Files.list(dir, ".txt", false);
		for (int i = 0; i < files.length; i++) {
			System.out.println(files[i]+"\t"+Files.countLines(dir+files[i], true));
		}
	}
	
	public static void parseCRF() {
		PrintWriter writer;
		String filename;
		
//		filename = DIR+"parseHits.crf";
		filename = DIR+"parsePs.crf";
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("hits");
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
			FilterByLists.process(DIR+"files/"+STUDIES[i][1], DIR+"duplicatedMarkers.txt", null, 0, DIR+"kept/"+STUDIES[i][1], true, log);
		}
	}

	public static void listQQfiles() {
		for (int i = 0; i < STUDIES.length; i++) {
//			System.out.print((i==0?"":";")+"files/"+STUDIES[i][1]+",11");
			System.out.print((i==0?"":";")+"filtered/"+STUDIES[i][1]+",5");
		}
		System.out.println();
	}
	
	public static void main(String[] args) {
		if (new File("/software/").exists()) {
			DIR = "";
		}
		
//		deCodedeCode();
//		parseAlleles();
//		parseCRF();
//		parseSuperSet();
		filterMarkers();
//		countLines(DIR+"filtered/");
//		listQQfiles();
//		keepMarkers();
	}
}
