package org.genvisis.park;

import java.io.*;

import org.genvisis.common.*;
import org.genvisis.db.crfDB;

public class UpdateCRFdb {
	public static final String CRF_DIR = tools.CRF_DIR;

	public static final String BACKUP_DIR = CRF_DIR+"backup/";

	public static final String PLATE_LIST_DIR = tools.MASTER_PLATELIST_DIR;

	public static final String[] TRAITS = {"Depressed", "MajorDepression", "MinorDepression", "Depression", "AOO", "parkin", "AnyLRRK2", "GBA_carrier", "GBA_del"};

	public static void update() throws IOException {
		System.out.println("Parsing ninfos...");
		ParseNinfos.parseNinfos(ParseNinfos.DEFAULT_DIR, "ninfo1.txt", "ninfo2.txt", "ninfo5.txt", "ninfo6.txt");
		Files.backupAndMove("ninfo1.csv", tools.NINFO_DIR, CRF_DIR, BACKUP_DIR);
		Files.backupAndMove("ninfo2.csv", tools.NINFO_DIR, CRF_DIR, BACKUP_DIR);
		Files.backupAndMove("ninfo5.csv", tools.NINFO_DIR, CRF_DIR, BACKUP_DIR);
		Files.backupAndMove("ninfo1_BirthDates.csv", tools.NINFO_DIR, CRF_DIR, BACKUP_DIR);
		Files.backupAndMove("ninfo5_BirthDates.csv", tools.NINFO_DIR, CRF_DIR, BACKUP_DIR);
		System.out.println("Parsing affected relatives...");
		Additionals.affRent();
		Files.backup("affRent.csv", CRF_DIR, BACKUP_DIR, false);
		System.out.println("Parsing probands...");
		updateProbands(PLATE_LIST_DIR+"PROGENI/", "probands.xln");
		Files.backupAndMove("probands.csv", PLATE_LIST_DIR+"PROGENI/", CRF_DIR, BACKUP_DIR);
		System.out.println("Building database...");
		try {
			new crfDB(CRF_DIR, "crf_db_key.crf");
		} catch (Elision e) {
			System.err.println("Error creating "+"crf_db.dat");
			System.exit(1);
		}
		for (int i = 0; i<TRAITS.length; i++) {
			System.out.println("Parsing family history of "+TRAITS[i]+"...");
			FamilyHistory.proc(TRAITS[i]);
			Files.backup(TRAITS[i]+"_sibHistory.csv", CRF_DIR, BACKUP_DIR, false);
		}
		System.out.println("Building database again...");
		try {
			new crfDB(CRF_DIR, "crf_db_key.crf");
		} catch (Elision e) {
			System.err.println("Error creating "+"crf_db.dat");
		}
		try {
			new crfDB(CRF_DIR, "CARES_crf_db_key.crf");
		} catch (Elision e) {
			System.err.println("Error creating "+"CARES_crf_db.dat");
		}
	}

	public static void updateProbands(String dir, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(dir+filename));
			writer = new PrintWriter(new FileWriter(dir+"probands.csv"));
			writer.println("FamID,IndID,proband,VPDproband");
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(Array.toStr(tools.getFamID(line[1]), ",")+",1,"+(line[6].equals("VPD")||line[6].equals("CONF_PD")?"1":"0"));
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+filename+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;

		String usage = "\n"+"park.UpdateCRFdb requires no arguments\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			update();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
