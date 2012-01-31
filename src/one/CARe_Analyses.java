package one;

import java.io.*;
import java.util.*;

import mining.Transformations;

import parse.GenParser;
import common.*;
import gwas.*;
import filesys.*;

public class CARe_Analyses {
	public static String DRIVE_ROOT = "not_set";

	public static final String[][] DEMO_NEEDS = {{"FID"}, {"IID"}, {"CAReID"}, {"icam"}, {"ln_icam"}, {"Age"}, {"Male"}, {"BMI"}, {"Race"}, {"rs5491"}, {"rs5491carrier"}};
	public static final String[][] DEMO_PARAMETERS = {
		{"# whites with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=White", "-blank"},
		{"mean ICAM levels in whites", "mean", "icam", "icam=ValidDouble", "Race=White", "-blank"},
		{"# blacks with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Black", "-blank"},
		{"mean ICAM levels in blacks", "mean", "icam", "icam=ValidDouble", "Race=Black", "-blank"},
		{"# Asians with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Asian", "-blank"},
		{"mean ICAM levels in Asians", "mean", "icam", "icam=ValidDouble", "Race=Asian", "-blank"},
		{"# Hispanics with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank"},
		{"mean ICAM levels in Hispanics", "mean", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank"},
	  };

	public static final String[][] COUNT_PARAMETERS = {
		{"# whites", "count", "icam", "Race=White", "-blank"},
		{"# whites with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=White", "-blank"},
		{"# whites / rs5491 carrier", "count", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=1"},
		{"# whites / rs5491 non-carrier", "count", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=0"},
		{"# whites / carrier status unknown", "count", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=."},

		{"# blacks", "count", "icam", "Race=Black", "-blank"},
		{"# blacks with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Black", "-blank"},
		{"# blacks / rs5491 carrier", "count", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=1"},
		{"# blacks / rs5491 non-carrier", "count", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=0"},
		{"# blacks / carrier status unknown", "count", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=."},

		{"# Asians", "count", "icam", "Race=Asian", "-blank"},
		{"# Asians with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Asian", "-blank"},
		{"# Asians / rs5491 carrier", "count", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=1"},
		{"# Asians / rs5491 non-carrier", "count", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=0"},
		{"# Asians / carrier status unknown", "count", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=."},

		{"# Hispanics", "count", "icam", "Race=Hispanic", "-blank"},
		{"# Hispanics with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank"},
		{"# Hispanics / rs5491 carrier", "count", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=1"},
		{"# Hispanics / rs5491 non-carrier", "count", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=0"},
		{"# Hispanics / carrier status unknown", "count", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=."},
	};

	public static final String[][] LEVEL_PARAMETERS = {
		{"mean ICAM levels in whites / rs5491 carrier", "mean", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=1"},
		{"mean ICAM levels in whites / rs5491 non-carrier", "mean", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=0"},
		{"mean ICAM levels in whites / carrier status unknown", "mean", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=."},

		{"mean ICAM levels in blacks / rs5491 carrier", "mean", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=1"},
		{"mean ICAM levels in blacks / rs5491 non-carrier", "mean", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=0"},
		{"mean ICAM levels in blacks / carrier status unknown", "mean", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=."},

		{"mean ICAM levels in Asians / rs5491 carrier", "mean", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=1"},
		{"mean ICAM levels in Asians / rs5491 non-carrier", "mean", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=0"},
		{"mean ICAM levels in Asians / carrier status unknown", "mean", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=."},

		{"mean ICAM levels in Hispanics / rs5491 carrier", "mean", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=1"},
		{"mean ICAM levels in Hispanics / rs5491 non-carrier", "mean", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=0"},
		{"mean ICAM levels in Hispanics / carrier status unknown", "mean", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=."},
	};

	public static final String[][] FINAL_PARAMETERS = {
		{"# whites with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=White", "-blank", "rs5491=0"},
		{"mean ICAM levels in whites", "mean", "icam", "icam=ValidDouble", "Race=White", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"mean age in whites", "mean", "Age", "icam=ValidDouble", "Race=White", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"% male in whites", "mean", "Male", "icam=ValidDouble", "Race=White", "-blank", "sf=2", "rs5491=0", "-percent"},
		{"mean BMI in whites", "mean", "BMI", "icam=ValidDouble", "Race=White", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"# blacks with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Black", "-blank", "rs5491=0"},
		{"mean ICAM levels in blacks", "mean", "icam", "icam=ValidDouble", "Race=Black", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"mean age in blacks", "mean", "Age", "icam=ValidDouble", "Race=Black", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"% male in blacks", "mean", "Male", "icam=ValidDouble", "Race=Black", "-blank", "sf=2", "rs5491=0", "-percent"},
		{"mean BMI in blacks", "mean", "BMI", "icam=ValidDouble", "Race=Black", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"# Asians with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "rs5491=0"},
		{"mean ICAM levels in Asians", "mean", "icam", "icam=ValidDouble", "Race=Asian", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"mean age in Asians", "mean", "Age", "icam=ValidDouble", "Race=Asian", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"% male in Asians", "mean", "Male", "icam=ValidDouble", "Race=Asian", "-blank", "sf=2", "rs5491=0", "-percent"},
		{"mean BMI in Asians", "mean", "BMI", "icam=ValidDouble", "Race=Asian", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"# Hispanics with ICAM levels", "count", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "rs5491=0"},
		{"mean ICAM levels in Hispanics", "mean", "icam", "icam=ValidDouble", "Race=Hispanic", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"mean age in Hispanics", "mean", "Age", "icam=ValidDouble", "Race=Hispanic", "-blank", "sf=2", "rs5491=0", "-stdev"},
		{"% male in Hispanics", "mean", "Male", "icam=ValidDouble", "Race=Hispanic", "-blank", "sf=2", "rs5491=0", "-percent"},
		{"mean BMI in Hispanics", "mean", "BMI", "icam=ValidDouble", "Race=Hispanic", "-blank", "sf=2", "rs5491=0", "-stdev"},
	};
	
	public static void generate_ABO_ICAM_datasets() {
		Logger log;
		long date;
		String dir, root, finalDir;
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		boolean[] familyBased = new boolean[] {false, false, true, false, true, false};
//		String[] studies = new String[] {"CFS", "FHS"};
//		boolean[] familyBased = new boolean[] {true, true};

		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};

		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/abo_parsing.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/"+studies[i]+"/IBC/"+races[j][0]+"/";
				root = dir+"leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr9";
				if (new File(root+".mlinfo").exists()) {
					new File(dir+"candis/").mkdirs();
					log.report("Parsing abo for "+studies[i]+" "+races[j][0]);
					date = new Date().getTime();
//					DosageData.convert(root+".gen", root+".pfam", root+".mlinfo", DosageData.GEN_FORMAT, dir+"candis/abo"+(familyBased[i]?".gwaf":".gen"), null, DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/abo_snplist.txt", familyBased[i]?DosageData.GWAF_FORMAT:DosageData.MACH_MLDOSE_FORMAT, false, true, log);
					DosageData.convert(root+".gen", root+".pfam", root+".mlinfo", DosageData.GEN_FORMAT, dir+"candis/abo.dosage", null, DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/abo_snplist.txt", DosageData.PLINK_FORMAT, false, true, log);
					log.report("Time elapsed: "+ext.getTimeElapsed(date));
				}
			}
		}
		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/icam_parsing.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/"+studies[i]+"/IBC/"+races[j][0]+"/";
				root = dir+"leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr19";
				if (new File(root+".mlinfo").exists()) {
					log.report("Parsing icam1 for "+studies[i]+" "+races[j][0]);
					date = new Date().getTime();
//					DosageData.convert(root+".gen", root+".pfam", root+".mlinfo", DosageData.GEN_FORMAT, dir+"candis/icam"+(familyBased[i]?".gwaf":".gen"), null, DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/icam1_snplist.txt", familyBased[i]?DosageData.GWAF_FORMAT:DosageData.MACH_MLDOSE_FORMAT, false, true, log);
					DosageData.convert(root+".gen", root+".pfam", root+".mlinfo", DosageData.GEN_FORMAT, dir+"candis/icam.dosage", null, DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/icam1_snplist.txt", DosageData.PLINK_FORMAT, false, true, log);
					log.report("Time elapsed: "+ext.getTimeElapsed(date));
				}
			}
		}
		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/concatenating.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/"+studies[i]+"/IBC/"+races[j][0]+"/";
				finalDir = DRIVE_ROOT+"ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				new File(finalDir).mkdirs();
				root = dir+"leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr19";
				if (new File(root+".gen").exists()) {
					new File(dir+"candis/abo.gen").delete();
					new File(dir+"candis/icam.gen").delete();
					new File(dir+"candis/abo_icam.gen").delete();
					new File(dir+"candis/abo_icam.mldose").delete();
					log.report("Concatenating abo and icam1 for "+studies[i]+" "+races[j][0]);
					Files.cat(new String[] {dir+"candis/abo.dosage", dir+"candis/icam.dosage"}, finalDir+"abo_icam.dosage", new int[] {0, 1}, log);
					Files.cat(new String[] {dir+"candis/abo.pinfo", dir+"candis/icam.pinfo"}, finalDir+"abo_icam.pinfo", new int[] {0, 1}, log);
					Files.copyFile(dir+"candis/abo.ids.fam", finalDir+"abo_icam.ids.fam");
					if (familyBased[i]) {
						new DosageData(finalDir+"abo_icam.dosage", finalDir+"abo_icam.ids.fam", finalDir+"abo_icam.pinfo", true, log).writeToFile(finalDir+"abo_icam.fhsR", finalDir+"abo_icam.mlinfo", null, DosageData.GWAF_FORMAT, log);
					}
				}
			}
		}
		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/coverting.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"abo_icam.dosage").exists()) {
					new DosageData(dir+"abo_icam.dosage", dir+"abo_icam.ids.fam", dir+"abo_icam.pinfo", true, log).writeToFile(dir+"abo_icam.mldose", dir+"abo_icam.mlinfo", null, DosageData.MACH_MLDOSE_FORMAT, log);
				}
			}
		}

		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/coverting.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"abo_icam.dosage").exists()) {
//					new File(dir+"abo_icam.mldose_list.out").delete();
					Mach.listDosage(dir, "abo_icam");
				}
			}
		}
	}
	
	public static void generate_ABO_ICAM_covariates() {
		Logger log;
		long date;
		String dir, root, finalDir;
		String[] hitlist;
//		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
//		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
//		boolean[] familyBased = new boolean[] {false, false, true, false, true, false};

		String[] studies = new String[] {"MESA"};
		String[][] races = new String[][] {{"asians", "ASN", "Asian", "asian"},{"hispanics", "HIS", "Hispanic", "hispanic"}};
		boolean[] familyBased = new boolean[] {false};


		log = new Logger(DRIVE_ROOT+"Analyses/ICAM/abo_icam_covars.log");
		Files.writeList(new String[] {"rs5498", "rs1799969", "rs651007"}, DRIVE_ROOT+"Analyses/ICAM/covarSNPs.txt");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"abo_icam.dosage").exists()) {
					log.report("Parsing abo/icam for "+studies[i]+" "+races[j][0]);
//					new DosageData(dir+"abo_icam.dosage", dir+"abo_icam.ids.fam", dir+"abo_icam.pinfo", DosageData.PLINK_FORMAT, true, log).writeToFile(dir+"abo_icam.dosage.csv", dir+"abo_icam.dosage.pinfo", DRIVE_ROOT+"Analyses/ICAM/covarSNPs.txt", true, DosageData.PARAMETERS[DosageData.GWAF_FORMAT], log);
					hitlist = new String[] {"hits", "abo_icam.ids.fam 1 hideIndex out=plink_phenoWithConditionals.dat", "abo_icam.ids.fam 1 0=FID 1=IID skip=0", "pheno/db_clean_wPCs.txt 1 4;-9 5;. 6;. 7;. fail", "abo_icam.dosage.csv , 0 "+Array.toStr(Array.stringArraySequence(new SnpMarkerSet(dir+"abo_icam.dosage.pinfo").getMarkerNames().length, ""), " "), "pheno/db_clean_wPCs.txt 1 8;. 9;. 10;. 11;. 12;. 13;. 14;. 15;. 16;. 17;. fail"};
					if (studies[i].equals("MESA") && (races[j][0].equals("asians") || races[j][0].equals("hispanics"))) {
						hitlist = Array.addStrToArray(DRIVE_ROOT+"Analyses/ICAM/IBC/whites/"+studies[i]+"/pheno/"+races[j][3]+"QC/rs1799969.xln 1 4=rs1799969c", hitlist, 5);
					}
					Files.writeList(hitlist, dir+"generatePhenoForPlinkWithConditionals.crf");
					CmdLine.run("java -cp C:/home/npankrat/park.jar -Xmx1024M Launch -suppress generatePhenoForPlinkWithConditionals.crf", dir);
					if (familyBased[i]) {
						Files.writeList(new String[] {"hits", "abo_icam.ids.fam 1 hideIndex out=phenoWithConditionals.csv", "leslie_lange."+studies[i]+".IBC."+races[j][1]+".Rlinker 0 2=id skip=0", "pheno/db_clean_wPCs.txt 1 4; 5; 6; 7; fail", "abo_icam.dosage.csv , 0 "+Array.toStr(Array.stringArraySequence(new SnpMarkerSet(dir+"abo_icam.dosage.pinfo").getMarkerNames().length, "")), "pheno/db_clean_wPCs.txt 1 8; 9; 10; 11; 12; 13; 14; 15; 16; 17; fail"}, dir+"generatePhenoForGWAFWithConditionals.crf");
						CmdLine.run("java -cp C:/home/npankrat/park.jar -Xmx1024M Launch -suppress generatePhenoForGWAFWithConditionals.crf", dir);
						new File("D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/").mkdirs();
						Files.copyFile(dir+"phenoWithConditionals.csv", "D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/phenoWithConditionals.csv");
					} else {
						new File("D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/").mkdirs();
						Files.copyFile(dir+"plink_phenoWithConditionals.dat", "D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/plink_phenoWithConditionals.dat");
					}
				}
			}
		}
		
		System.exit(1);
		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/icam_parsing.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/"+studies[i]+"/IBC/"+races[j][0]+"/";
				root = dir+"leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr19";
				if (new File(root+".mlinfo").exists()) {
					log.report("Parsing icam1 for "+studies[i]+" "+races[j][0]);
					date = new Date().getTime();
//					DosageData.convert(root+".gen", root+".pfam", root+".mlinfo", DosageData.GEN_FORMAT, dir+"candis/icam"+(familyBased[i]?".gwaf":".gen"), null, DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/icam1_snplist.txt", familyBased[i]?DosageData.GWAF_FORMAT:DosageData.MACH_MLDOSE_FORMAT, false, true, log);
					DosageData.convert(root+".gen", root+".pfam", root+".mlinfo", DosageData.GEN_FORMAT, dir+"candis/icam.dosage", null, DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/icam1_snplist.txt", DosageData.PLINK_FORMAT, false, true, log);
					log.report("Time elapsed: "+ext.getTimeElapsed(date));
				}
			}
		}
		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/concatenating.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/"+studies[i]+"/IBC/"+races[j][0]+"/";
				finalDir = DRIVE_ROOT+"ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				new File(finalDir).mkdirs();
				root = dir+"leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr19";
				if (new File(root+".gen").exists()) {
					new File(dir+"candis/abo.gen").delete();
					new File(dir+"candis/icam.gen").delete();
					new File(dir+"candis/abo_icam.gen").delete();
					new File(dir+"candis/abo_icam.mldose").delete();
					log.report("Concatenating abo and icam1 for "+studies[i]+" "+races[j][0]);
					Files.cat(new String[] {dir+"candis/abo.dosage", dir+"candis/icam.dosage"}, finalDir+"abo_icam.dosage", new int[] {0, 1}, log);
					Files.cat(new String[] {dir+"candis/abo.pinfo", dir+"candis/icam.pinfo"}, finalDir+"abo_icam.pinfo", new int[] {0, 1}, log);
					Files.copyFile(dir+"candis/abo.ids.fam", finalDir+"abo_icam.ids.fam");
					if (familyBased[i]) {
						new DosageData(finalDir+"abo_icam.dosage", finalDir+"abo_icam.ids.fam", finalDir+"abo_icam.pinfo", true, log).writeToFile(finalDir+"abo_icam.fhsR", finalDir+"abo_icam.mlinfo", null, DosageData.GWAF_FORMAT, log);
					}
				}
			}
		}
		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/coverting.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"abo_icam.dosage").exists()) {
					new DosageData(dir+"abo_icam.dosage", dir+"abo_icam.ids.fam", dir+"abo_icam.pinfo", true, log).writeToFile(dir+"abo_icam.mldose", dir+"abo_icam.mlinfo", null, DosageData.MACH_MLDOSE_FORMAT, log);
				}
			}
		}

		log = new Logger(DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/coverting.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"abo_icam.dosage").exists()) {
//					new File(dir+"abo_icam.mldose_list.out").delete();
					Mach.listDosage(dir, "abo_icam");
				}
			}
		}
	}
	
	public static void validateLines() {
		BufferedReader reader, mapReader, infoReader;
		String dir, root;
		int count, fileNum;
		String[] header, line;
		String marker;
		int numInds;
		Logger log;
		
		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/FHS/IBC/whites/";
		root = "leslie_lange.FHS.IBC.CEU.chr";
//		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/whites/";
//		root = "leslie_lange.CFS.IBC.CEU.chr";
//		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/blacks/";
//		root = "leslie_lange.CFS.IBC.AAM.chr";

//		dir = "";

		log = new Logger(dir+"validatingFiles.out");
		fileNum = 0;
		try {
			numInds = HashVec.loadFileToStringArray(dir+root+"1.pfam", false, new int[] {0,1}, false).length;
			for (int chr = 1; chr <= 22; chr++) {
//			for (int chr = 15; chr <= 15; chr++) {
				log.report("Checking chr"+chr);
				try {
					reader = new BufferedReader(new FileReader(dir+root+chr+".gen"));
					mapReader = new BufferedReader(new FileReader(dir+root+chr+".pmap"));
					infoReader = new BufferedReader(new FileReader(dir+root+chr+".mlinfo"));
					header = infoReader.readLine().trim().split("[\\s]+");
					ext.checkHeader(header, new String[] {"SNP", "Al1", "Al2", "Freq1", "MAF", "Quality", "Rsq"}, true);
					count = 0;
					while (reader.ready()) {
						count++;
						line = mapReader.readLine().trim().split("[\\s]+");
						marker = line[1];
						if (line.length != 4) {
							log.reportError("Error at record marker '"+marker+"' on line "+count+" in file "+root+chr+".pmap; expecting "+4+" columns, but found "+line.length);
						}
						
						line = infoReader.readLine().trim().split("[\\s]+");
						if (!line[0].equals(marker)) {
							log.reportError("Error - mismatched marker on line "+count+" in file "+root+chr+".mlinfo; expecting '"+marker+"' but got '"+line[0]+"'");
						}
						if (line.length != 7) {
							log.reportError("Error at record marker '"+marker+"' on line "+count+" in file "+root+chr+".mlinfo; expecting "+7+" columns, but found "+line.length);
						}

						line = reader.readLine().trim().split("[\\s]+");
						if (!line[1].equals(marker)) {
							log.reportError("Error - mismatched marker on line "+count+" in file "+root+chr+".geb; expecting '"+marker+"' but got '"+line[1]+"'");
						}
						if (line.length != numInds*3+5) {
							log.reportError("Error at record marker '"+marker+"' on line "+count+" in file "+root+chr+".mlinfo; expecting "+(numInds*3+5)+" columns, but found "+line.length);
							log.reportError("The next unexpected token is "+line[numInds*3+5]+" preceeded by "+line[numInds*3+5-1]);
							return;
						}
					}
					reader.close();
					mapReader.close();
					infoReader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+root+chr+".gen" + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"leslie_lange.FHS.IBC.CEU.chr1.gen" + "\"");
					System.exit(2);
				}
			}
		} catch (Exception e) {
			System.err.println("Error writing to file #" + fileNum);
			e.printStackTrace();
		}
	}

	public static void standardizeAndConvert(int binSize, int resumeAt) {
		BufferedReader reader, mapReader, infoReader;
		PrintWriter writer, mapWriter, infoWriter;
		String header, dir, root;
		int count, fileNum;
//		Logger log;
//		
//		log = new Logger();

		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/FHS/IBC/whites/";
		root = "leslie_lange.FHS.IBC.CEU.chr";
		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/whites/";
		root = "leslie_lange.CFS.IBC.CEU.chr";
//		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/blacks/";
//		root = "leslie_lange.CFS.IBC.AAM.chr";

//		dir = "";

		count = 0;
		fileNum = 0;
		new File(dir+"std/").mkdirs();
		new File(dir+"gwaf/").mkdirs();
		// make sure you convert CARe IDs to numerical values during the linker
//		Files.copyFile(dir+root+"1.pfam", dir+"std/ids.pfam"); // e.g. this won't work
		Files.combine(HashVec.loadFileToStringArray(dir+root+"1.pfam", false, new int[] {1}, false), new String[] {dir+ext.rootOf(root)+".Rlinker 0 1 2 skip=0"}, null, "nohead", dir+"std/ids.pfam", new Logger(), false, false, true, false);
		try {
			writer = new PrintWriter(new FileWriter(dir+"std/file"+fileNum+".gen"));
			mapWriter = new PrintWriter(new FileWriter(dir+"std/file"+fileNum+".pmap"));
			infoWriter = new PrintWriter(new FileWriter(dir+"std/file"+fileNum+".mlinfo"));
			for (int chr = 1; chr <= 22; chr++) {
				try {
					reader = new BufferedReader(new FileReader(dir+root+chr+".gen"));
					mapReader = new BufferedReader(new FileReader(dir+root+chr+".pmap"));
					infoReader = new BufferedReader(new FileReader(dir+root+chr+".mlinfo"));
					header = infoReader.readLine();
					while (reader.ready()) {
						if (count == 0) {
							infoWriter.println(header);
						}
						if (fileNum >= resumeAt) {
							writer.println(reader.readLine());
							mapWriter.println(mapReader.readLine());
							infoWriter.println(infoReader.readLine());
						} else {
							reader.readLine();
							mapReader.readLine();
							infoReader.readLine();
						}
						count++;
						if (count == binSize) {
							writer.close();
							mapWriter.close();
							infoWriter.close();
							
//							new DosageData(dir+"std/file"+fileNum+".gen", dir+"std/ids.pfam", dir+"std/file"+fileNum+".mlinfo", true, log).writeToFile(dir+"gwaf/file"+fileNum+".fhsR", dir+"gwaf/file"+fileNum+".pmap", null, DosageData.GWAF_FORMAT, log);
//							Zip.gzip(dir+"gwaf/file"+fileNum+".fhsR");
//							System.out.println(fileNum);
							// delete intermediate .gen file
							// delete unzipped .fhsR file
							
							fileNum++;
							writer = new PrintWriter(new FileWriter(dir+"std/file"+fileNum+".gen"));
							mapWriter = new PrintWriter(new FileWriter(dir+"std/file"+fileNum+".pmap"));
							infoWriter = new PrintWriter(new FileWriter(dir+"std/file"+fileNum+".mlinfo"));
							count = 0;
						}
					}
					reader.close();
					mapReader.close();
					infoReader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+root+chr+".gen" + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"leslie_lange.FHS.IBC.CEU.chr1.gen" + "\"");
					System.exit(2);
				}
			}
			writer.close();
			mapWriter.close();
			infoWriter.close();
		} catch (Exception e) {
			System.err.println("Error writing to file #" + fileNum);
			e.printStackTrace();
		}
	}

	public static void justConvert(int startAt) {
		String dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/FHS/IBC/whites/";
		Logger log = new Logger(dir+"conversion.log");
		int fileNum;
		
//		for (int chr = 22; chr >= 1; chr--) {
//			new DosageData(dir+"leslie_lange.FHS.IBC.CEU.chr"+chr+".gen", dir+"leslie_lange.FHS.IBC.CEU.chr"+chr+".pfam", dir+"leslie_lange.FHS.IBC.CEU.chr"+chr+".mlinfo", true, log).writeToFile(dir+"gwaf/"+"leslie_lange.FHS.IBC.CEU.chr"+chr+".fhsR", dir+"gwaf/"+"leslie_lange.FHS.IBC.CEU.chr"+chr+".pinfo", null, DosageData.GWAF_FORMAT, log);
//		}

		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/FHS/IBC/whites/";
		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/whites/";
//		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/blacks/";
//		dir = "";
		
		new File(dir+"gwaf/").mkdirs();
		fileNum = startAt;
		while (Files.exists(dir+"std/file"+fileNum+".gen", false)) {
			new DosageData(dir+"std/file"+fileNum+".gen", dir+"std/ids.pfam", dir+"std/file"+fileNum+".mlinfo", true, log).writeToFile(dir+"gwaf/file"+fileNum+".fhsR", dir+"gwaf/file"+fileNum+".pmap", null, DosageData.GWAF_FORMAT, log);
			Zip.gzip(dir+"gwaf/file"+fileNum+".fhsR");
			fileNum++;
		}
	}
	
	public static void batch(int numBatches) {
//		BufferedReader reader;
		PrintWriter[] writers;
//		String[] line;
//		String temp, trav;
		int count;
		String dir, root;
		
		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/FHS/IBC/whites/";
//		dir = "";
		root = "leslie_lange.FHS.IBC.CEU.chr";
		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/whites/";
		root = "leslie_lange.CFS.IBC.CEU.chr";
//		dir = DRIVE_ROOT+"CARe_imputed_all_llange_24mar2010/CFS/IBC/blacks/";
//		root = "leslie_lange.CFS.IBC.AAM.chr";
		

		new File(dir+"gwaf/").mkdirs();
		Files.copyFile(dir+ext.rootOf(root)+".pedfile", dir+"gwaf/"+ext.rootOf(root)+".pedfile");
		count = 0;
		try {
			writers = new PrintWriter[numBatches];
			for (int i = 1; i <= numBatches; i++) {
				writers[i-1] = new PrintWriter(new FileWriter(dir+"gwaf/run"+i+".bat"));
				writers[i-1].println("R --no-save < gwaf"+i+".R");
				writers[i-1].close();
				writers[i-1] = new PrintWriter(new FileWriter(dir+"gwaf/gwaf"+i+".R"));
				writers[i-1].println("library(kinship)");
				writers[i-1].println("library(GWAF)");
			}
			while (new File(dir+"gwaf/file"+count+".fhsR").exists()) {
				writers[count%numBatches].println("lme.batch.imputed(\"pheno.csv\", \"file"+count+".fhsR.gz\", \""+ext.rootOf(root)+".pedfile\", \"pheno\", \"kmat.Rfile\", covars=NULL, \"results"+count+".csv\", col.names=T, sep.ped=\",\", sep.phe=\",\", sep.gen=\",\")");
				count++;
			}
			for (int i = 1; i <= numBatches; i++) {
				writers[i-1].close();
			}
		} catch (Exception e) {
			System.err.println("Error queuing up file " + count);
			e.printStackTrace();
		}
	}
	
	public static void generateForests() {
		Hashtable<String, String> hash;
		BufferedReader reader;
		PrintWriter writer;
		String[] line, dirs;
		String dir;
//		Logger log;
		
		hash = HashVec.loadFileToHashNull("targets.txt", false);
		dirs = new String[] {"whites/ARIC/", "whites/CARDIA/", "whites/CFS/", "whites/CHS/", "whites/FHS/", "whites/MESA/", "blacks/CARDIA/", "blacks/CHS/", "blacks/CFS/", "blacks/MESA/", "asians/MESA/", "hispanics/MESA/", "whites_i3_", "blacks_i3_", ""};

//		log = new Logger("forests.log");
		for (int i = 0; i < dirs.length; i++) {
			dir = dirs[i];
			if (new File(dir+"plink_pheno_iteration1.InvVar1.out").exists()) {
				try {
					reader = new BufferedReader(new FileReader(dir+"plink_pheno_iteration1.InvVar1.out"));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (hash.containsKey(line[0])) {
							try {
								writer = new PrintWriter(new FileWriter(line[0]+".txt", true));
								if (new File(line[0]+".txt").length() == 0) {
									writer.println("Cohort\tSNP\tA1\tA2\tBeta\tSE");
								}
								writer.println((dir.equals("")?"Joint_analysis":ext.replaceAllWith(dir, "/", "_").substring(0, dir.length()-1))+"\t"+line[0]+"\t"+line[1].toUpperCase()+"\t"+line[2].toUpperCase()+"\t"+line[3]+"\t"+line[4]);
								writer.close();
							} catch (Exception e) {
								System.err.println("Error writing to something");
								e.printStackTrace();
							}
						}
						
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+"plink_pheno_iteration1.InvVar1.out"
							+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"plink_pheno_iteration1.InvVar1.out"
							+ "\"");
					System.exit(2);
				}
			}
			if (new File(dir+"plink_pheno_iteration1.se.metal").exists()) {
				try {
					reader = new BufferedReader(new FileReader(dir+"plink_pheno_iteration1.se.metal"));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (hash.containsKey(line[0])) {
							try {
								writer = new PrintWriter(new FileWriter(line[0]+".txt", true));
								if (new File(line[0]+".txt").length() == 0) {
									writer.println("Cohort\tSNP\tA1\tA2\tBeta\tSE");
								}
								writer.println((dir.equals("")?"Joint_analysis":ext.replaceAllWith(dir, "/", "_").substring(0, dir.length()-1))+"\t"+line[0]+"\t"+line[1].toUpperCase()+"\t"+line[2].toUpperCase()+"\t"+line[6]+"\t"+line[7]);
								writer.close();
							} catch (Exception e) {
								System.err.println("Error writing to something");
								e.printStackTrace();
							}
						}
						
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+"plink_pheno_iteration1.InvVar1.out"
							+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"plink_pheno_iteration1.InvVar1.out"
							+ "\"");
					System.exit(2);
				}
			}
		}
	}
	
	public static void parseDemographics() {
		Logger log;
		long date;
//		String[] header;
		String dir;
//		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String[][] races = new String[][] {{"whites", "CEU"}};
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
//		boolean[] familyBased = new boolean[] {false, false, true, false, true, false};
		Vector<String> files, fileDescriptions;

		date = new Date().getTime();
		files = new Vector<String>();
		fileDescriptions = new Vector<String>();
		log = new Logger(DRIVE_ROOT+"Analyses/ICAM/parseDemographics.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				if (new File(dir).exists()) {
					if (new File(dir+"pheno/db.txt").exists()) {
//						log.reportError("Found pheno/db.txt in "+races[j][0]+"/"+studies[i]+"");
//						header = Files.getHeaderOfFile(dir+"pheno/db.txt", "\t", log);
//						indices = ext.indexFactors(DEMO_NEEDS, header, true, true, true, false);
						files.add(dir+"pheno/db.txt");
//						fileDescriptions.add(studies[i]+" "+races[j][0]);
						fileDescriptions.add(studies[i]);
					} else {
						log.reportError("Found directory for "+races[j][0]+"/"+studies[i]+", but no pheno/db.txt");
					}
				}
			}
		}
		Files.generateTables(DRIVE_ROOT+"Analyses/ICAM/raw_demo.xln", Array.toStringArray(files), Array.toStringArray(fileDescriptions), DEMO_PARAMETERS, log);
		Files.generateTables(DRIVE_ROOT+"Analyses/ICAM/counts_comparisons.xln", Array.toStringArray(files), Array.toStringArray(fileDescriptions), COUNT_PARAMETERS, log);
		Files.generateTables(DRIVE_ROOT+"Analyses/ICAM/level_comparisons.xln", Array.toStringArray(files), Array.toStringArray(fileDescriptions), LEVEL_PARAMETERS, log);
		Files.generateTables(DRIVE_ROOT+"Analyses/ICAM/final_demographics.xln", Array.toStringArray(files), Array.toStringArray(fileDescriptions), FINAL_PARAMETERS, log);
		
		System.out.println("Time elapsed: "+ext.getTimeElapsed(date));
		
		
	}

	public static void analyzeGenotypedVariant(String variantID) {
		Logger log;
		String dir;
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		Hashtable<String, Vector<String>> vs;
//		Vector<String> v;

		log = new Logger(DRIVE_ROOT+"Analyses/running_"+variantID+".log", true);

//		for (int i = 0; i < studies.length; i++) {
//			for (int j = 0; j < races.length; j++) {
//				dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
//				if (new File(dir+"plink_pheno.dat").exists()) {
//					log.report("Running in "+studies[i]+" "+races[j][0]);
//					CmdLine.run("plink --noweb --bfile "+DRIVE_ROOT+"Analyses/ICAM/IBC/whites/"+studies[i]+"/pheno/full --snps "+variantID+" --linear --ci 0.95 --pheno plink_pheno.dat --covar plink_pheno.dat --covar-name Age,Male,BMI,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10 --out "+DRIVE_ROOT+"Analyses/ICAM/"+races[j][0]+"_"+studies[i], dir);
////					CmdLine.run("plink --noweb --bfile "+DRIVE_ROOT+"Analyses/ICAM/IBC/whites/"+studies[i]+"/pheno/full --snps "+variantID+" --linear --pheno plink_pheno.dat --out "+DRIVE_ROOT+"Analyses/ICAM/"+races[j][0]+"_"+studies[i]+"noCovar", dir);
//				}
//			}
//		}

		vs = new Hashtable<String, Vector<String>>();
//		v = new Vector<String>();
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"plink_pheno.dat").exists()) {
					log.report("Parsing "+studies[i]+" "+races[j][0]);
//					CmdLine.run("plink --noweb --bfile "+DRIVE_ROOT+"Analyses/ICAM/IBC/whites/"+studies[i]+"/pheno/full --snps "+variantID+" --freq --out "+DRIVE_ROOT+"Analyses/ICAM/"+races[j][0]+"_"+studies[i], dir);
//					Metal.convert(DRIVE_ROOT+"Analyses/ICAM/", races[j][0]+"_"+studies[i]+".assoc.linear", "ADD", "linear", races[j][0]+"_"+studies[i]+".frq", true, true, races[j][0]+"_"+studies[i]+".se.metal", false);
					HashVec.addToHashVec(vs, races[j][0], races[j][0]+"_"+studies[i]+".se.metal", false);
					HashVec.addToHashVec(vs, "Joint", races[j][0]+"_"+studies[i]+".se.metal", false);
				}
			}
		}
		
		Metal.metaAnalyze(DRIVE_ROOT+"Analyses/ICAM/", Array.toStringArray(vs.get("Joint")), "Joint_SE", true, null, log);
		Metal.metaAnalyze(DRIVE_ROOT+"Analyses/ICAM/", Array.toStringArray(vs.get("whites")), "Whites_SE", true, null, log);
		Metal.metaAnalyze(DRIVE_ROOT+"Analyses/ICAM/", Array.toStringArray(vs.get("blacks")), "Blacks_SE", true, null, log);
		if (true)
		return;

		Hashtable<String, String> hash;
		BufferedReader reader;
		PrintWriter writer;
		String[] line, dirs;
		
		hash = HashVec.loadFileToHashNull("targets.txt", false);
		dirs = new String[] {"whites/ARIC/", "whites/CARDIA/", "whites/CFS/", "whites/CHS/", "whites/FHS/", "whites/MESA/", "blacks/CARDIA/", "blacks/CHS/", "blacks/CFS/", "blacks/MESA/", "asians/MESA/", "hispanics/MESA/", "whites_i3_", "blacks_i3_", ""};

//		log = new Logger("forests.log");
		for (int i = 0; i < dirs.length; i++) {
			dir = dirs[i];
			if (new File(dir+"plink_pheno_iteration1.InvVar1.out").exists()) {
				try {
					reader = new BufferedReader(new FileReader(dir+"plink_pheno_iteration1.InvVar1.out"));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (hash.containsKey(line[0])) {
							try {
								writer = new PrintWriter(new FileWriter(line[0]+".txt", true));
								if (new File(line[0]+".txt").length() == 0) {
									writer.println("Cohort\tSNP\tA1\tA2\tBeta\tSE");
								}
								writer.println((dir.equals("")?"Joint_analysis":ext.replaceAllWith(dir, "/", "_").substring(0, dir.length()-1))+"\t"+line[0]+"\t"+line[1].toUpperCase()+"\t"+line[2].toUpperCase()+"\t"+line[3]+"\t"+line[4]);
								writer.close();
							} catch (Exception e) {
								System.err.println("Error writing to something");
								e.printStackTrace();
							}
						}
						
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+"plink_pheno_iteration1.InvVar1.out"
							+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"plink_pheno_iteration1.InvVar1.out"
							+ "\"");
					System.exit(2);
				}
			}
			if (new File(dir+"plink_pheno_iteration1.se.metal").exists()) {
				try {
					reader = new BufferedReader(new FileReader(dir+"plink_pheno_iteration1.se.metal"));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (hash.containsKey(line[0])) {
							try {
								writer = new PrintWriter(new FileWriter(line[0]+".txt", true));
								if (new File(line[0]+".txt").length() == 0) {
									writer.println("Cohort\tSNP\tA1\tA2\tBeta\tSE");
								}
								writer.println((dir.equals("")?"Joint_analysis":ext.replaceAllWith(dir, "/", "_").substring(0, dir.length()-1))+"\t"+line[0]+"\t"+line[1].toUpperCase()+"\t"+line[2].toUpperCase()+"\t"+line[6]+"\t"+line[7]);
								writer.close();
							} catch (Exception e) {
								System.err.println("Error writing to something");
								e.printStackTrace();
							}
						}
						
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+"plink_pheno_iteration1.InvVar1.out" + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"plink_pheno_iteration1.InvVar1.out" + "\"");
					System.exit(2);
				}
			}
		}	
	}
	
	public static void upToDate(String phenofile, String pheno, boolean conditionals) {
		PrintWriter writer;
		Vector<String> v;
		String dir;
		int count;

		dir = "";

		count = 0;
		v = new Vector<String>();
		try {
			writer = new PrintWriter(new FileWriter("catchUp."+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")));
			while (new File(dir+"gwaf/file"+count+".fhsR.gz").exists()) {
				if (!new File(ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_results"+count+".csv").exists()) {
					writer.println("qsub batches/"+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_file"+count+".qsub");
					v.add(count+"");
				}
				count++;
			}
			writer.close();
			Files.chmod("catchUp."+ext.rootOf(phenofile)+(conditionals?"_withCondi":""));
			count--;
			System.out.println("Files remaining: "+(v.size()==0?"none":ext.listRanges(Array.toIntArray(Array.toStringArray(v)))));
		} catch (Exception e) {
			System.err.println("Error writing to " + "catchUp");
			e.printStackTrace();
		}
	}
	
	public static void plinkBatches(String phenofile, String pheno, boolean conditionals) {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
//		String[] studies = new String[] {"MESA"};
//		String[][] races = new String[][] {{"blacks", "AAM"},{"whites", "CEU"}};
		String dir, commands;
		PrintWriter writer;
		
		try {
			writer = new PrintWriter(new FileWriter("batchAllPlink."+pheno));
			for (int i = 0; i < studies.length; i++) {
				for (int j = 0; j < races.length; j++) {
//					dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
					dir = studies[i]+"_"+races[j][0]+"/";
					if (new File(dir+phenofile).exists()) {
						commands = "/home/npankrat/bin/plink --noweb --fam leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr#.pfam "+
								"--map leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr#.pmap "+
								"--dosage 00src/leslie_lange."+studies[i]+".IBC."+races[j][1]+".chr#.gen.gz Zin skip0=1 skip1=1 skip2=0 format=3 noheader "+
								"--pheno "+phenofile+" --pheno-name "+pheno+" "+
//								"--covar "+phenofile+" --covar-name Age,Male,BMI"+(conditionals?",rs5498,rs651007"+(studies[i].equals("MESA")&&(races[j][0].equals("asians")||races[j][0].equals("hispanics"))?"":",rs1799969"):"")+",EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10 "+
								"--covar "+phenofile+" --covar-name Age,Male,BMI"+(conditionals?",rs5498c,rs651007c,rs1799969c":"")+",EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10 "+
								"--out "+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_chr#";
						Files.writeList(Files.qsub(dir, "chr#_"+studies[i], 1, 22, commands, null, -1, null), dir+"master."+pheno);
						Files.chmod(dir+"master."+pheno);
						writer.println("cd "+dir);
						writer.println("./master."+pheno);
						writer.println("cd ..");
					}
				}
			}
			writer.close();
			Files.chmod("batchAllPlink."+pheno);
		} catch (Exception e) {
			System.err.println("Error writing to " + "batchAllPlink."+pheno);
			e.printStackTrace();
		}
	}
	
	public static void gwafBatches(String phenofile, String pheno, boolean conditionals) {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String dir, commands;
		PrintWriter writer;
		
		try {
			writer = new PrintWriter(new FileWriter("batchAllGWAF."+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")));
			for (int i = 0; i < studies.length; i++) {
				for (int j = 0; j < races.length; j++) {
					dir = studies[i]+"_"+races[j][0]+"/";
					if (new File(dir+phenofile).exists()) {
						commands = "jcp gwas.GWAF "+
						"phenoFile="+phenofile+" "+
						"pheno="+pheno+" "+
						"genoPrimer=gwaf/file#.fhsR.gz "+
						"outfileTemplate="+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_results#.csv "+
						"qsubRoot=batches/"+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_file#.qsub "+
						"startAt=0 "+
						"imputed=true "+
						"pedfile=pedfile.csv "+
						"nodesToUse=compute-0-0.local,compute-0-3.local,compute-0-4.local "+ //
						"covars=Age,Male,BMI"+(conditionals?",rs651007c,rs1799969c,rs5498c":"")+",EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10";
						writer.println("cd "+dir);
						writer.println(commands);
						writer.println("./master."+pheno);
						writer.println("cd ..");
					}
				}
			}
			writer.close();
			Files.chmod("batchAllGWAF."+ext.rootOf(phenofile)+(conditionals?"_withCondi":""));
		} catch (Exception e) {
			System.err.println("Error writing to " + "batchAllGWAF."+pheno);
			e.printStackTrace();
		}
	}
	
	public static void parsePlinkResults(String phenofile, String pheno, boolean conditionals) {
//		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
//		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String[] studies = new String[] {"MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"hispanics", "HIS"}};
		String dir;
		String[] files;
		int[] skips;
		
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = studies[i]+"_"+races[j][0]+"/";
				if (new File(dir+phenofile).exists()) {
					files = new String[22];
					skips = new int[22];
					for (int chr = 1; chr <= 22; chr++) {
						files[chr-1] = dir+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_chr"+chr+".assoc.dosage";
						skips[chr-1] = chr==1?0:1;
					}
					Files.cat(files, ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_"+dir.substring(0, dir.length()-1)+"_dosage_results.out", skips, new Logger());
				}
			}
		}
	}
	
	public static void parseGWAFResults(String phenofile, String pheno, boolean conditionals) {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String dir;
		String[] files;
		int[] skips;
		
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				dir = studies[i]+"_"+races[j][0]+"/";
				if (new File(dir+phenofile).exists()) {
					files = new String[22];
					skips = new int[22];
					for (int chr = 1; chr <= 22; chr++) {
						files[chr-1] = dir+ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_chr"+chr+".assoc.dosage";
						skips[chr-1] = chr==1?0:1;
					}
					Files.cat(files, ext.rootOf(phenofile)+(conditionals?"_withCondi":"")+"_"+dir.substring(0, dir.length()-1)+"_dosage_results.out", skips, new Logger());
				}
			}
		}
	}
	
	public static void parseMetaFiles(String label, boolean conditionals) {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String dir;
		String filename;
		Logger log;
		int count;
		
		if (label != null) {
			label = ext.rootOf(label)+(conditionals?"_withCondi":"");
		}
//		dir = "C:\\CARe_data\\conditionalMeta\\";
		dir = "";
		count = 0;
		log = new Logger(dir+"parseMetaInputs"+(label==null?"":"_"+label)+".log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
//				filename = studies[i]+"_"+races[j][0]+"_dosage_results.out";
				filename = "plink_"+(label==null?"":label+"_")+studies[i]+"_"+races[j][0]+"_dosage_results.out";
				if (new File(dir+filename).exists()) {
					GenParser.parse(new String[] {dir+filename, "out="+dir+(label==null?"":label+"_")+studies[i]+"_"+races[j][0]+".input", "1=MarkerName", "3=A1", "4=A2", "!6>0.30", "7=beta", "8=StdErr", "!7<5", "!7>-5", "!8!NA", "!8!0"}, log);
					count++;
				}

//				filename = studies[i]+"_"+races[j][0]+"_results.csv"; inverseNormalizedPheno_withCondi_CFS_whites_results.csv
				filename = (label==null?"":label+"_")+studies[i]+"_"+races[j][0]+"_results.csv";
				if (new File(dir+filename).exists()) {
//	        		Files.combine(HashVec.loadFileToStringArray(dir+studies[i]+"_"+races[j][0]+"_snps.txt", true, new int[] {0}, false), new String[] {dir+studies[i]+"_"+races[j][0]+".mlinfo 0 1=A1 2=A2", dir+studies[i]+"_"+races[j][0]+"_results.csv , 1 5=beta 6=StdErr"}, null, "MarkerName", dir+studies[i]+"_"+races[j][0]+".input", log, false, true, false, false);
	        		Files.combine(HashVec.loadFileToStringArray(dir+studies[i].toLowerCase()+"_"+races[j][0].toLowerCase()+"_snps.txt", true, new int[] {0}, false), new String[] {dir+studies[i].toLowerCase()+"_"+races[j][0].toLowerCase()+".mlinfo 0 1=A1 2=A2", dir+filename+" , 0 1=beta;. 2=StdErr;. fail"}, null, "MarkerName", dir+(label==null?"":label+"_")+studies[i]+"_"+races[j][0]+".input", log, false, true, false, false);
	        		count++;
				}
			}
		}
		log.report("Prepped "+count+" files");
	}
	
	public static void metaAnalyze(String label, boolean conditionals) {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String dir;
		String filename;
		Logger log;
		Vector<String> all, race;
		
		if (label != null) {
			label = ext.rootOf(label)+(conditionals?"_withCondi":"");
		}
//		dir = "C:\\CARe_data\\conditionalMeta\\";
		dir = "";
		log = new Logger(dir+"metaAnalyzing"+(label==null?"":"_"+label)+".log");
		all = new Vector<String>();
		for (int j = 0; j < races.length; j++) {
			race = new Vector<String>();
			for (int i = 0; i < studies.length; i++) {
				filename = (label==null?"":label+"_")+studies[i]+"_"+races[j][0]+".input";
				if (new File(dir+filename).exists()) {
					all.add(filename);
					race.add(filename);
				}
			}
			Metal.metaAnalyze(dir, Array.toStringArray(race), (label==null?"":label+"_")+races[j][0]+"_metal", true, null, log);
		}
		Metal.metaAnalyze(dir, Array.toStringArray(all), (label==null?"":label+"_")+"all_metal", true, null, log);
	}
	
	public static void allMarkers() {
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String dir;
		Vector<String> all;
		
		dir = "C:\\CARe_data\\conditionalMeta\\";
		all = new Vector<String>();
		for (int j = 0; j < races.length; j++) {
			all.add(dir+races[j][0]+"_metal1.out");
		}
		Unique.proc(Array.toStringArray(all), dir+"allSNPs.txt", dir+"allSNP_counts.txt", true);
	}

	public static void parseFreqFiles() {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN"},{"blacks", "AAM"},{"hispanics", "HIS"},{"whites", "CEU"}};
		String dir;
		String filename;
		Logger log;
		
		dir = "C:\\CARe_data\\conditionalMeta\\";
		log = new Logger(dir+"parseFinals.log");
		for (int i = 0; i < studies.length; i++) {
			for (int j = 0; j < races.length; j++) {
				filename = studies[i]+"_"+races[j][0]+"_dosage_results.out";
				if (new File(dir+filename).exists()) {
					GenParser.parse(new String[] {dir+filename, "out="+dir+studies[i]+"_"+races[j][0]+".freq", "1=MarkerName", "3=A1", "4=A2", "5=Freq", "6=Rsq", "!6>0.30"}, log);
				}

				filename = studies[i]+"_"+races[j][0]+"_results.csv";
				if (new File(dir+filename).exists()) {
	        		Files.combine(HashVec.loadFileToStringArray(dir+studies[i]+"_"+races[j][0]+"_snps.txt", true, new int[] {0}, false), new String[] {dir+studies[i]+"_"+races[j][0]+".mlinfo 0 1=A1 2=A2 3=Freq 6=Rsq"}, null, "MarkerName", dir+studies[i]+"_"+races[j][0]+".freq", log, false, true, false, false);
				}
			}
		}
	}
	
	public static void normalizePhenos() {
//		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
//		String[][] races = new String[][] {{"asians", "ASN", "Asian"},{"blacks", "AAM", "Black"},{"hispanics", "HIS", "Hispanic"},{"whites", "CEU", "White"}};
//		boolean[] familyBased = new boolean[] {false, false, true, false, true, false};
		
		String[] studies = new String[] {"MESA"};
		String[][] races = new String[][] {{"asians", "ASN", "Asian", "asian"},{"hispanics", "HIS", "Hispanic", "hispanic"}};
		boolean[] familyBased = new boolean[] {false};
		
		String dir;
		Logger log;
		String[] conditionals, hitlist;
		
		log = new Logger("C:/CARe_data/Analyses/ICAM/"+"normalizePhenos.log");
		for (int i = 0; i < studies.length; i++) {
//			dir = "C:/CARe_data/Analyses/ICAM/IBC/whites/"+studies[i]+"/";
			for (int j = 0; j < races.length; j++) {
				dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
				if (new File(dir+"pheno/db.txt").exists()) {
					GenParser.parse(new String[] {dir+"pheno/db.txt", "0", "1", "2", "4", "8", "!4!.", "!9=0", "!8="+races[j][2], "out="+dir+"pheno/db_"+races[j][0]+".txt"}, log);
//				if (Files.countLines(dir+"pheno/db_"+races[j][0]+".txt", true) == 0) {
//					new File(dir+"pheno/db_"+races[j][0]+".txt").delete();
//				} else {
					Transformations.transformFile(dir+"pheno/db_"+races[j][0]+".txt", dir+"pheno/db_"+races[j][0]+"_normalized.txt", true, 3, false, false, Transformations.NORMALIZE, log);
					Transformations.transformFile(dir+"pheno/db_"+races[j][0]+".txt", dir+"pheno/db_"+races[j][0]+"_inverseNormalized.txt", true, 3, false, false, Transformations.INVERSE_NORMALIZE, log);

					dir = DRIVE_ROOT+"Analyses/ICAM/IBC/"+races[j][0]+"/"+studies[i]+"/";
					if (new File(dir+"abo_icam.dosage").exists()) {
//						log.report("Parsing abo/icam for "+studies[i]+" "+races[j][0]);
//						new DosageData(dir+"abo_icam.dosage", dir+"abo_icam.ids.fam", dir+"abo_icam.pinfo", DosageData.PLINK_FORMAT, true, log).writeToFile(dir+"abo_icam.dosage.csv", dir+"abo_icam.dosage.pinfo", DRIVE_ROOT+"Analyses/ICAM/covarSNPs.txt", true, DosageData.PARAMETERS[DosageData.GWAF_FORMAT], log);
						
						conditionals = new SnpMarkerSet(dir+"abo_icam.dosage.pinfo").getMarkerNames();
						for (int k = 0; k < conditionals.length; k++) {
							conditionals[k] = (k+1)+"="+conditionals[k]+"c";
						}
						
						hitlist = new String[] {"hits", "abo_icam.ids.fam 1 hideIndex out=plink_normalizedPheno.dat", "abo_icam.ids.fam 1 0=FID 1=IID skip=0", "pheno/db_"+races[j][0]+"_normalized.txt 1 4;-9 fail", "pheno/db_clean_wPCs.txt 2 5;. 6;. 7;. fail", "abo_icam.dosage.csv , 0 "+Array.toStr(conditionals), "pheno/db_clean_wPCs.txt 2 8;. 9;. 10;. 11;. 12;. 13;. 14;. 15;. 16;. 17;. fail"};
						if (studies[i].equals("MESA") && (races[j][0].equals("asians") || races[j][0].equals("hispanics"))) {
							hitlist = Array.addStrToArray(DRIVE_ROOT+"Analyses/ICAM/IBC/whites/"+studies[i]+"/pheno/"+races[j][3]+"QC/rs1799969.xln 1 4=rs1799969c", hitlist, 6);
						}
						Files.writeList(hitlist, dir+"generateNormalizedPhenoForPlinkWithConditionals.crf");
						CmdLine.run("java -cp C:/home/npankrat/park.jar -Xmx1024M Launch -suppress generateNormalizedPhenoForPlinkWithConditionals.crf", dir);

						hitlist = new String[] {"hits", "abo_icam.ids.fam 1 hideIndex out=plink_inverseNormalizedPheno.dat", "abo_icam.ids.fam 1 0=FID 1=IID skip=0", "pheno/db_"+races[j][0]+"_inverseNormalized.txt 1 4;-9 fail", "pheno/db_clean_wPCs.txt 2 5;. 6;. 7;. fail", "abo_icam.dosage.csv , 0 "+Array.toStr(conditionals), "pheno/db_clean_wPCs.txt 2 8;. 9;. 10;. 11;. 12;. 13;. 14;. 15;. 16;. 17;. fail"};
						if (studies[i].equals("MESA") && (races[j][0].equals("asians") || races[j][0].equals("hispanics"))) {
							hitlist = Array.addStrToArray(DRIVE_ROOT+"Analyses/ICAM/IBC/whites/"+studies[i]+"/pheno/"+races[j][3]+"QC/rs1799969.xln 1 4=rs1799969c", hitlist, 6);
						}
						Files.writeList(hitlist, dir+"generateInverseNormalizedPhenoForPlinkWithConditionals.crf");
						CmdLine.run("java -cp C:/home/npankrat/park.jar -Xmx1024M Launch -suppress generateInverseNormalizedPhenoForPlinkWithConditionals.crf", dir);
						if (familyBased[i]) {
							Files.writeList(new String[] {"hits", "abo_icam.ids.fam 1 hideIndex out=normalizedPhenoWithConditionals.csv", "leslie_lange."+studies[i]+".IBC."+races[j][1]+".Rlinker 0 2=id skip=0", "pheno/db_"+races[j][0]+"_normalized.txt 1 4; fail", "pheno/db_clean_wPCs.txt 2 5; 6; 7; fail #VALUE!=> tab", "abo_icam.dosage.csv , 0 "+Array.toStr(conditionals), "pheno/db_clean_wPCs.txt 2 8; 9; 10; 11; 12; 13; 14; 15; 16; 17; fail"}, dir+"generateNormalizedPhenoForGWAFWithConditionals.crf");
							CmdLine.run("java -cp C:/home/npankrat/park.jar -Xmx1024M Launch -suppress generateNormalizedPhenoForGWAFWithConditionals.crf", dir);
							Files.writeList(new String[] {"hits", "abo_icam.ids.fam 1 hideIndex out=inverseNormalizedPhenoWithConditionals.csv", "leslie_lange."+studies[i]+".IBC."+races[j][1]+".Rlinker 0 2=id skip=0", "pheno/db_"+races[j][0]+"_inverseNormalized.txt 1 4; fail", "pheno/db_clean_wPCs.txt 2 5; 6; 7; fail #VALUE!=> tab", "abo_icam.dosage.csv , 0 "+Array.toStr(conditionals), "pheno/db_clean_wPCs.txt 2 8; 9; 10; 11; 12; 13; 14; 15; 16; 17; fail"}, dir+"generateInverseNormalizedPhenoForGWAFWithConditionals.crf");
							CmdLine.run("java -cp C:/home/npankrat/park.jar -Xmx1024M Launch -suppress generateInverseNormalizedPhenoForGWAFWithConditionals.crf", dir);
							new File("D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/").mkdirs();
							Files.copyFile(dir+"normalizedPhenoWithConditionals.csv", "D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/normalizedPhenoWithConditionals.csv");
							Files.copyFile(dir+"inverseNormalizedPhenoWithConditionals.csv", "D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/inverseNormalizedPhenoWithConditionals.csv");
						} else {
							new File("D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/").mkdirs();
							Files.copyFile(dir+"plink_normalizedPheno.dat", "D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/plink_normalizedPheno.dat");
							Files.copyFile(dir+"plink_inverseNormalizedPheno.dat", "D:/upload/phenos/"+studies[i]+"_"+races[j][0]+"/plink_inverseNormalizedPheno.dat");
						}
					}
				
				}
			}
		}
	}
	
	public static void cleanUp(String note) {
		String[] studies = new String[] {"ARIC", "CARDIA", "CFS", "CHS", "FHS", "MESA"};
		String[][] races = new String[][] {{"asians", "ASN", "Asian"},{"blacks", "AAM", "Black"},{"hispanics", "HIS", "Hispanic"},{"whites", "CEU", "White"}};
		boolean[] familyBased = new boolean[] {false, false, true, false, true, false};
		String dir;
		PrintWriter writer;
		
		
		try {
			writer = new PrintWriter(new FileWriter("cleanup."+note));
			for (int i = 0; i < studies.length; i++) {
				for (int j = 0; j < races.length; j++) {
					dir = studies[i]+"_"+races[j][0]+"/";
					if (new File(dir).exists()) {
						new File(dir+"results/"+note+"/").mkdirs();
//						CmdLine.run("mv "+(familyBased[i]?"resul*.csv":"*.dosage")+" *.log *.qsub* results/"+note+"/", dir);
						writer.println("cd "+dir);
						writer.println("mv "+(familyBased[i]?"resul*.csv":"*.dosage")+" *.log *.qsub* results/"+note+"/");
						writer.println("cd ..");
					}
				}
			}
			writer.close();
			Files.chmod("cleanup."+note);
		} catch (Exception e) {
			System.err.println("Error writing to " + "cleanup."+note);
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "CARe_Analyses.dat";

		String usage = "\n" + "one.CARe_Analyses requires 0-1 arguments\n"
				+ "   (1) filename (i.e. file=" + filename + " (default))\n"
				+ "";

		if (args.length > 0) {
			if (args[0].equals("cleanup")) {
				System.out.println("Cleaning up and moving things to results/"+args[1]+"/");
				cleanUp(args[1]);
				return;
			}
			
			if (args[0].equals("batchPlink")) {
				System.out.println("Batching plink: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				plinkBatches(args[1], args[2], ext.parseBooleanArg(args[3]));
				return;
			}		
	
			if (args[0].equals("batchGWAF")) {
				System.out.println("Batching GWAF: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				gwafBatches(args[1], args[2], ext.parseBooleanArg(args[3]));
				return;
			}
			
			if (args[0].equals("parsePlink")) {
				System.out.println("Parsing plink results: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				parsePlinkResults(args[1], args[2], ext.parseBooleanArg(args[3]));
				return;
			}		
			
			if (args[0].equals("parseGWAF")) {
				System.out.println("Parsing plink results: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				parsePlinkResults(args[1], args[2], ext.parseBooleanArg(args[3]));
				return;
			}

			if (args[0].equals("upToDate")) {
				System.out.println("Seeing which GWAF runs have yet to be completed: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				upToDate(args[1], args[2], ext.parseBooleanArg(args[3]));
				return;
			}
			
			if (args[0].equals("prepMeta")) {
				System.out.println("Creating METAL input files: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				parseMetaFiles(args[1], ext.parseBooleanArg(args[3]));
				return;
			}
			
			if (args[0].equals("runMeta")) {
				System.out.println("Running METAL by race: file="+args[1]+" pheno="+args[2]+" includeConditionals="+args[3]);
				metaAnalyze(args[1], ext.parseBooleanArg(args[3]));
				return;
			}
			
		}
		
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help")
					|| args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}


		if (new File("C:/CARe_data/").exists()) {
			DRIVE_ROOT = "C:/CARe_data/";	
		} else {		
			DRIVE_ROOT = "D:/CARe/";
		}		
		
//		generate_ABO_ICAM_datasets();
//		validateLines();
//		standardizeAndConvert(10000, -1);
//		justConvert(0);
//		batch(5);

		
//		GWAF.batch("phenoWithConditionals.csv", "pheno", "", "gwaf/file#.fhsR.gz", 0, "pedfile.csv", "results#.csv", "compute-0-0.local,compute-0-3.local,compute-0-4.local", -1);

		
		
//		generateForests();
//		parseDemographics();
//		analyzeGenotypedVariant("rs5030400");
//		parseMetaFiles();
//		metaAnalyze();
//		allMarkers();
//		parseFreqFiles();
//		generate_ABO_ICAM_covariates();
//		normalizePhenos();

		
//		cd ../CFS_blacks/
//		tar -zcvf ../cfs_blacks_results.tar.gz resul*.csv

	}
}
