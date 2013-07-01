package one;

import java.io.*;
import java.util.*;

import common.*;
import filesys.BurdenMatrix;
import filesys.GenotypeMatrix;
import filesys.Hits;
import filesys.SerialHash;
import gwas.Metal;

public class ChargeS {
//	public static final String[][] MAP_REQUIREMENTS = {{"SNP"}, {"Chr", "CHROM"}, {"Position", "POS"}, {"ALT"}, {"REF"}};
//	public static final String[][] ANNOTATION_REQUIREMENTS = {{"MAF"}, {"gene"}};

//	public static final String[][] BURDEN_REQUIREMENTS = {{"SNP"}, {"Chr", "CHROM"}, {"Position", "POS"}, {"ALT"}, {"REF"}, {"gene"}, {"AAF"}, {"use"}}; // , {"Function", "RFG"} // , {"MAF"} doesn't flip alleles
	public static final String[][] BURDEN_REQUIREMENTS = {{"SNP"}, {"Chr", "CHROM"}, {"Position", "POS"}, {"ALT"}, {"REF"}, {"gene"}, {"AAF"}, {"sc_nonsynSplice"}}; // , {"Function", "RFG"} // , {"MAF"} doesn't flip alleles

//	public static final String[] STUDIES = {"ARIC", "CHS", "FHS"}; // +ESP
//	public static final String[] STUDIES = {"ARIC", "CHS", "FHS", "ESP"};
	public static final String[] STUDIES = {"ARIC", "CHS", "FHS", "ESP6800"};
	public static final String[] STUDY_GROUPS = {"CHARGE", "CHARGE", "CHARGE", "ESP"}; // +ESP
	
	public static final String[][] PHENOTYPES = {
		{"Fibrinogen", "fibrinogen", ".lnFB."}, 
		{"F7", ".FVII.", "_FVII_"}, 
		{"F8", ".FVIII.", "_FVIII_"}, 
		{"vWF", "VWF"}
	};

	public static final int[][] DEFAULT_SAMPLE_SIZES = {
		{985, 628, 255, 2011}, // Fibr
		{965, 630, 248, 1221}, // F7
		{984, 626, 0, 1296}, // F8
		{985, 0, 249, 1184}, // vWF
	};
	public static final int[][] FREEZE3_SAMPLE_SIZES = {
		{3168, 741, 499, 2013}, // Fibr
		{3092, 744, 416, 1222}, // F7
		{3166, 736, 0, 1296}, // F8
		{3168, 0, 416, 1185}, // vWF
	};
	public static final String[][] METHODS = {{"SingleSNP", "singleSnpRes", ".LR."}, {"T5Count", ".T5."}, {"T5MB", "MBT5"}};
	public static final String[][] UNIT_OF_ANALYSIS = {Metal.MARKER_NAMES, Metal.GENE_UNITS, Metal.GENE_UNITS};
	public static final boolean[] SINGLE_VARIANTS = {true, false, false};
	public static final boolean[] WEIGHTED = {true, false, false};
	public static final String[] GROUPS = {"SingleVariant", "BurdenTests", "BurdenTests"};	
	public static final String[][] GROUP_ANNOTATION_PARAMS = {
		{},
//		{"../SNPInfo_ExomeFreeze2_120810_aafSlim.csv 'SNP' 'gene' 'AAF'=CHARGE_AF"},
		{},
	};

	public static final String SNP_INFO_FILE = "snpinfo_ChargeSFreeze3_ESP_05212013.RData";
	public static final String SNP_NAMES = "SNP";
	public static final String CHROM_NAME = "CHROM";
	public static final String GENE_NAME = "SKATgene";
	
	
//	wts =1, mafRange = c(0,0.01),
			
	public static void metaAll(String dir) {
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
	
	public static void runAll(String phenoFile, String genoFile, String annotationFile) {
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Logger log;
		GenotypeMatrix gens;
		Hashtable<String,String> annotationHash;
		int[] cols;
		String[] markerNames, keys;
		boolean chrPosCombo;
		
		if (true) { // two classes
			chrPosCombo = true;
		}		
		
		log = new Logger(ext.rootOf(phenoFile, false)+".log");
		if (!Files.exists(genoFile+".burdenInfo")) {
			cols = ext.indexFactors(BURDEN_REQUIREMENTS, ext.replaceAllWith(Files.getFirstNLinesOfFile(annotationFile, 1, log)[0], "\"", "").split(Files.determineDelimiter(annotationFile, log)), false, true, true, log, true);
			log.report("Loading map positions");
			annotationHash = HashVec.loadFileToHashString(annotationFile, new int[] {cols[0]}, cols, annotationFile.endsWith(".csv"), "\t", true, false, false);
			log.report("Getting keys");
			keys = HashVec.getKeys(annotationHash, false, false);
			log.report("Removing quotes");
			for (int i = 0; i < keys.length; i++) {
				if (keys[i].contains("\"")) {
					annotationHash.put(ext.removeQuotes(keys[i]), annotationHash.remove(keys[i]));
				}
			}
			log.report("Reading marker names required");
			if (chrPosCombo) {
				markerNames = HashVec.loadFileToStringArray(genoFile, false, true, new int[] {0,1}, false, false, Files.determineDelimiter(genoFile, log));
				for (int i = 0; i < markerNames.length; i++) {
					markerNames[i] = "chr"+ext.replaceAllWith(markerNames[i], "\t", ":");
				}
			} else {
				markerNames = HashVec.loadFileToStringArray(genoFile, false, true, new int[] {0}, false, false, Files.determineDelimiter(genoFile, log));
			}
			log.report("Writing map to file");
			try {
				writer = new PrintWriter(new FileWriter(genoFile+".burdenInfo"));
				writer.println("Marker\tChr\tPosition\tREF\tALT\tgene\tAAF\tFunction");
				for (int i = 0; i < markerNames.length; i++) {
					trav = ext.removeQuotes(markerNames[i]);
					temp = annotationHash.get(trav);
					if (temp == null) {
						log.reportError("Error - no information available for marker\t"+trav);
//						writer.println("0\t"+trav+"\t0\t0\t?\t?");
						writer.println(trav+"\t0\t0\t?\t?\t.\t.\t.");
					} else {
						line = ext.replaceAllWith(temp, "\"", "").split("\t", -1);
//						writer.println(line[1]+"\t"+trav+"\t0\t"+line[2]+"\t"+line[3]+"\t"+line[4]);
						writer.println(Array.toStr(line));
					}
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + genoFile+".burdenInfo");
				e.printStackTrace();
			}
		}
		if (Files.exists(genoFile+".ser")) {
			System.out.println("Loading serialized version: "+genoFile+".ser");
			gens = GenotypeMatrix.load(genoFile+".ser", false);
		} else {
			System.out.println("Loading: "+genoFile);
			gens = new GenotypeMatrix(genoFile, null, genoFile+".burdenInfo", log);
			System.out.println("Saving: "+genoFile+".ser");
			gens.serialize(genoFile+".ser");
		}
		System.out.println("Analyzing single variants with "+phenoFile);
		gens.analyze(phenoFile, "NA", null, true, log);
		System.out.println("Generating T5 burden dataset");
		BurdenMatrix burdenMatrix = new BurdenMatrix(gens, 0.05, BurdenMatrix.DEFAULT_ANNOTATIONS_TO_INCLUDE, BurdenMatrix.ALL_KNOWN_ANNOTATIONS, false, 0.01, new String[] {}, null, log);
		System.out.println("Running T1 burden test with "+phenoFile);
		burdenMatrix.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5.EA."+ext.getDate(new Date(),"")+".csv", true, log);
	}

	public static void convertESP(String plinkResults, String plinkHWE, String outfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String[]> hash;
		double freq;
		String[] subline;
		Logger log;
		String delimiter;
		double se;
		
		System.out.println(ext.getTime());
		log = new Logger(plinkHWE+".log");
		
		if (Files.exists(plinkHWE+".ser")) {
			hash = SerialHash.loadSerializedStringArrayHash(plinkHWE+".ser");
		} else {
			hash = new Hashtable<String, String[]>(2000000);
			try {
				reader = new BufferedReader(new FileReader(plinkHWE));
				line = reader.readLine().trim().split("[\\s]+");
				ext.checkHeader(line, new String[] {"CHR", "SNP", "TEST", "A1", "A2", "GENO", "O(HET)", "E(HET)", "P"}, true);
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					subline = line[5].trim().split("/");
					try {
						freq = (Double.parseDouble(subline[0])*2+Double.parseDouble(subline[1])*1)/(Double.parseDouble(subline[0])*2+Double.parseDouble(subline[1])*2+Double.parseDouble(subline[2])*2);
					} catch (Exception e) {
						System.err.println("Error - parsing "+line[5]);
						freq = -999;
					}
					if (line[2].equals("ALL")) {
						hash.put(ext.replaceAllWith(line[1], "_", ":"), new String[] {line[4], line[3], (Integer.parseInt(subline[0])+Integer.parseInt(subline[1])+Integer.parseInt(subline[2]))+"", freq+"", subline[0], subline[1], subline[2]});
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + plinkHWE + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + plinkHWE + "\"");
				System.exit(2);
			}
			SerialHash.createSerializedStringArrayHash(plinkHWE+".ser", hash);
		}
		
		System.out.println(ext.getTime());
		
		try {
			reader = new BufferedReader(new FileReader(plinkResults));
			writer = new PrintWriter(new FileWriter(outfile));
			delimiter = Files.suggestDelimiter(outfile, log);
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, new String[] {"CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"}, true);
			writer.println(Array.toStr(new String[] {"snp", "noncoded_all", "coded_all", "sampleSumN", "AF", "sampleAA", "sampleAR", "sampleRR", "beta", "se", "pval", "NMISS", "wbeta", "wse"}, delimiter));
			
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[4].equals("ADD")) {
					line[1] = ext.replaceAllWith(line[1], "_", ":");
	//				if (hash.containsKey(line[1])) {
						subline = hash.get(line[1]);
						if (!line[3].equals(subline[1])) { // A1 alleles
							log.reportError("Error - mismatched A1 alleles for plink "+line[3]+" and hwe "+subline[0]+" in "+plinkResults);
						}
	//				} else {
	//					System.err.println("Error - no hwe info for "+line[1]);
	//				}
					try {
						se = Double.parseDouble(line[6])/Double.parseDouble(line[7]);
					} catch (Exception e) {
						log.reportError("Could not divide: "+line[6]+"/"+line[7]);
						se = -999;
					}
					writer.println(line[1]+delimiter+Array.toStr(subline, delimiter)+delimiter+line[6]+delimiter+se+delimiter+line[8]+delimiter+line[5]+delimiter+line[6]+delimiter+se);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + plinkResults + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + plinkResults + "\"");
			System.exit(2);
		}

		System.out.println(ext.getTime());
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "ChargeS.dat";
		String phenoFile = null;
		String genoFile = "C:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/";
		String annotationFile = "C:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/SNPInfo_ExomeFreeze2_120810.csv";

		String usage = "\n" + "one.ChargeS requires 0-1 arguments\n"
				+ "   (1) filename (i.e. file=" + filename + " (default))\n"
				+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				phenoFile = args[i].split("=")[1];
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
			String dir = "D:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/";
			genoFile = dir+"aric_genotypes_frz2_final.csv";
//			annotationFile = "SNPInfo_ExomeFreeze2_120810.csv";

//			genoFile = "slim.csv";
//			annotationFile = "slimMarkers2.csv";

//			genoFile = "chr1pter.csv";
//			genoFile = "samd.csv";
//			annotationFile = "chr1pter_SNPInfo.csv";
//			annotationFile = "SNPInfo_ExomeFreeze2_120810_min.csv";

			annotationFile = dir+"SNPInfo_ExomeFreeze2_120810_aafSlim.csv";
			
//			phenoFile = dir+"lnFibrinogen.csv";
			
			
//			convertESP(dir+"results/ESP/HeartGO_WHISP_fibri_ea.assoc.linear", dir+"results/ESP/HWE_ea_fibri.hwe", dir+"results/ESP/ESP_Fibrinogen_SingleSNP.csv");
//			convertESP(dir+"results/ESP/HeartGO_WHISP_factor7_ea.assoc.linear", dir+"results/ESP/HWE_ea_factor7.hwe", dir+"results/ESP/ESP_F7_SingleSNP.csv");
//			convertESP(dir+"results/ESP/HeartGO_WHISP_factor8_ea.assoc.linear", dir+"results/ESP/HWE_ea_factor8.hwe", dir+"results/ESP/ESP_F8_SingleSNP.csv");
//			convertESP(dir+"results/ESP/HeartGO_WHISP_vwf_ea.assoc.linear", dir+"results/ESP/HWE_ea_vwf.hwe", dir+"results/ESP/ESP_vWF_SingleSNP.csv");
//			System.exit(1);

//			dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/test/";
//			genoFile = dir+"ARIC_EA_Freeze3_Chrom17_Genotypes.txt.gz";
//			annotationFile = dir+"test_SNPInfo.csv";
//			phenoFile = dir+"F7_all.csv";
//			phenoFile = dir+"F7_all_frz2.csv";
//			phenoFile = dir+"F7_all_replIn3.csv";
//			
//			dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/test/";
//			genoFile = dir+"ARIC_EA_Freeze3_Chrom19_Genotypes.txt.gz";
//			annotationFile = dir+"chr19_SNPInfo.csv";
//			phenoFile = dir+"F7_all_replIn3.csv";
//			phenoFile = dir+"F8_all_replIn3.csv";
//			phenoFile = dir+"vWF_all_replIn3.csv";
//			phenoFile = dir+"F8_all.csv";
//			phenoFile = dir+"vWF_all.csv";
//						
//			
//			
//
//			if (phenoFile != null) {
//				runAll(phenoFile, genoFile, annotationFile);
//			}
			
//			metaAll(dir+"results/");
			
			
//			runAll(dir+"lnFibrinogen.csv", dir+genoFile, dir+annotationFile);
//			runAll(dir+"FVII.csv", dir+genoFile, dir+annotationFile);
//			runAll(dir+"FVIII.csv", dir+genoFile, dir+annotationFile);
//			runAll(dir+"VWF.csv", dir+genoFile, dir+annotationFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
