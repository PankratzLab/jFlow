package park;

import java.io.*;
import java.util.*;

import common.*;

public class ParseNinfos {
	public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\00ninfos\\";

	// public static final String[] RAW_NINFO1_HEADER = {"FamNo", "IndNo",
	// "AgeOfOnset", "IllnessStatusCode", "DNA.", "Sex", "IllnessCode",
	// "MonthofBirth", "DayofBirth", "YearofBirth", "Q28", "Q29", "prox_aod",
	// "prox_aoo", "RaceCode", "Autopsy?", "AutopsyStatus", "CenterCode"};
	// public static final String[] RAW_NINFO1_HEADER = {"FamNo", "IndNo",
	// "AgeOfOnset", "IllnessStatusCode", "DNA", "Sex", "IllnessCode",
	// "MonthofBirth", "DayofBirth", "YearofBirth", "Q28", "Q29", "prox_aod",
	// "prox_aoo", "RaceCode", "Autopsy?", "AutopsyStatus", "CenterCode"};
	public static final String[] RAW_NINFO1_HEADER = {"FamNo", "IndNo", "AgeOfOnset", "IllnessStatusCode", "DNA", "Sex", "IllnessCode", "MonthofBirth", "DayofBirth", "YearofBirth", "Q28", "Q29", "prox_aod", "prox_aoo", "RaceCode", "Autopsy?", "AutopsyStatus", "CenterCode", "SDate", "Status"};

	public static final String[] RAW_NINFO5_HEADER = {"FamNo", "IndNo", "AgeOfOnset", "IllnessStatusCode", "DNANum", "Sex", "IllnessCode", "MonthofBirth", "DayofBirth", "YearofBirth", "Q28", "Q29", "prox_aod", "prox_aoo", "RaceCode", "Autopsy?", "AutopsyStatus", "SDate", "Status"};

	public static final String LAST_PERM_COL = "Autopsy?";

	public static final int[] FINAL_NINFO1_INDICES = {0, 1, 2, 3, 4, 5, 6, 14};

	public static final int[] DOB_INDICES = {7, 8, 9};

	public static void parseNinfos(String dir, String ninfo1, String ninfo2, String ninfo5, String ninfo6) throws IOException {
		Vector<String> inds;

		System.out.println("Processing ninfo1...");
		inds = procNinfoType1(dir+ninfo1, RAW_NINFO1_HEADER);
		Files.backup("ninfo1.dat", dir, dir+"backup/", false);

		System.out.println("Processing ninfo2...");
		procNinfoType2(dir+ninfo2, tools.NINFO2_HEADER, inds);
		Files.backup("ninfo2.dat", dir, dir+"backup/", false);
		if (new File(dir+"ninfo3.dat").exists()) {
			Files.backup("ninfo3.dat", dir, dir+"backup/, false", false);
		}

		System.out.println("Processing ninfo5...");
		inds = procNinfoType1(dir+ninfo5, RAW_NINFO5_HEADER);
		Files.backup("ninfo5.dat", dir, dir+"backup/", false);

		System.out.println("Processing ninfo6...");
		procNinfoType2(dir+ninfo6, tools.NINFO2_HEADER, inds);
		Files.backup("ninfo6.dat", dir, dir+"backup/", false);
	}

	public static Vector<String> procNinfoType1(String filename, String[] header) {
		BufferedReader reader = null;
		PrintWriter writer_CSV, writer_DAT, writer_DOB, writer;
		String[] line, master;
		String temp, trav;
		Hashtable<String,Vector<String[]>> hash = new Hashtable<String,Vector<String[]>>();
		Vector<String[]> v;
		Vector<String> inds = new Vector<String>();
		Vector<String> phenos, dnas;
		int lastCol;
		Vector<String> missingPheno = new Vector<String>();
		Vector<String> missingDOB = new Vector<String>();

		try {
			try {
				reader = new BufferedReader(new FileReader(filename));
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+filename+"\" not found in current directory");
				System.exit(1);
			}
			trav = filename.substring(filename.lastIndexOf("/")+1, filename.lastIndexOf("."));
			writer_CSV = new PrintWriter(new FileWriter(trav+".csv"));
			writer_DAT = new PrintWriter(new FileWriter(trav+".dat"));
			writer_DOB = new PrintWriter(new FileWriter(trav+"_BirthDates.csv"));

			line = reader.readLine().split("[\\s]+");
			ext.checkHeader(line, header, true);
			lastCol = ext.indexOfStr(LAST_PERM_COL, header);

			for (int i = 0; i<line.length; i++) {
				writer_CSV.print((i==0?"":",")+line[i]);
			}
			writer_CSV.println();

			for (int i = 0; i<FINAL_NINFO1_INDICES.length; i++) {
				writer_DAT.print((i==0?"":"\t")+line[FINAL_NINFO1_INDICES[i]]);
			}
			writer_DAT.println();

			writer_DOB.println("FamID,IndID,DOB");

			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split("\t");
				if (line.length!=header.length) {
					if (line.length>=lastCol&&line.length<=header.length&&(line[lastCol].equals("TRUE")||line[lastCol].equals("FALSE"))) {} else {
						System.err.println("Error - Expecting "+header.length+" columns, finding only "+line.length);
						System.err.println("      - "+temp);
						System.exit(3);
					}
				}
				if (hash.containsKey(line[0]+"\t"+line[1])) {
					v = hash.get(line[0]+"\t"+line[1]);
				} else {
					hash.put(line[0]+"\t"+line[1], v = new Vector<String[]>());
					inds.add(line[0]+"\t"+line[1]);
				}
				v.add(line);

				writer_DOB.println(line[0]+","+line[1]+","+procDOB(line, DOB_INDICES, missingDOB));
			}
			reader.close();
			writer_DOB.close();

			for (int i = 0; i<inds.size(); i++) {
				v = hash.get(inds.elementAt(i));
				master = new String[header.length];
				phenos = new Vector<String>();
				dnas = new Vector<String>();
				for (int j = 0; j<v.size(); j++) {
					line = v.elementAt(j);
					if (!phenos.contains(line[3]+"-"+line[6])) {
						phenos.add(line[3]+"-"+line[6]);
					}
					if (!dnas.contains(line[4])) {
						dnas.add(line[4]);
					}
					for (int k = 0; k<line.length; k++) {
						if (master[k]==null||!master[k].equals(line[k])) {
							if (j==0) {
								master[k] = line[k];
							} else {
								if (k<3&&k>6) {
									System.err.println("Error - Unexpectedly inconsistent data for "+line[0]+"-"+line[1]);
								}
							}
						}
					}
				}
				if (phenos.contains("CONF-PD")||phenos.contains("CONF-DLB")||phenos.contains("CONF-LB")) {
					master[3] = "CONF_PD";
					master[6] = "PD";
				} else if (phenos.contains("VPD-PD")) {
					master[3] = "VPD";
					master[6] = "PD";
				} else if (phenos.contains("NVPD-PD")) {
					master[3] = "NVPD";
					master[6] = "PD";
				} else if (phenos.contains("RPD-PD")) {
					master[3] = "NVPD";
					master[6] = "PD";
				} else if (phenos.contains("NRPD-PD")) {
					master[3] = "NRPD";
					master[6] = "PD";
				} else if (phenos.contains("NOEV-PD")) {
					master[3] = "NOEV";
					master[6] = "PD";
				} else if (phenos.contains("-")||phenos.contains("-PD")) {
					missingPheno.add(master[0]+"\t"+master[1]+"\t"+Array.toStr(Array.toStringArray(dnas), "\t"));
					master[3] = ".";
					master[6] = ".";
				} else if (phenos.contains("VPD-PD")&&phenos.contains("NVPD-PD")) {
					System.err.println("Error - "+master[0]+"-"+master[1]+" is listed as both VPD and NVPD");
				} else if (phenos.size()>1) {
					System.err.println("Error - "+master[0]+"-"+master[1]+" is listed as several things, none of which are PD ("+ext.listWithCommas(Array.toStringArray(phenos))+")");
				}
				if (dnas.size()>1) {
					master[4] = "";
					for (int j = 0; j<dnas.size(); j++) {
						master[4] += (j==0?"":":")+dnas.elementAt(j);
					}
				}
				for (int j = lastCol; j<master.length; j++) {
					if (master[j]==null) {
						master[j] = "";
					}
				}

				if (master[2].equals("")) {
					master[2] = ".";
				}
				for (int j = 0; j<master.length; j++) {
					writer_CSV.print((j==0?"":",")+master[j]);
				}
				for (int j = 0; j<FINAL_NINFO1_INDICES.length; j++) {
					writer_DAT.print((j==0?"":"\t")+(master[FINAL_NINFO1_INDICES[j]].equals("")?".":master[FINAL_NINFO1_INDICES[j]])); // altered
					// to .
					// out
					// blanks
				}
				writer_CSV.println();
				writer_DAT.println();
			}
			writer_CSV.close();
			writer_DAT.close();

			if (missingPheno.size()>0) {
				System.err.println("There was no phenotypic information for "+missingPheno.size()+" individual"+(missingPheno.size()>1?"s":"")+"; see missingPhenos.xls for details.");
				writer = new PrintWriter(new FileWriter(trav+"_missingPhenos.xln"));
				writer.println("FamID\tIndID\tDNAs");
				for (int i = 0; i<missingPheno.size(); i++) {
					writer.println(missingPheno.elementAt(i));
				}
				writer.close();
			}

			if (missingPheno.size()>0) {
				System.err.println("There was no DOB information for "+missingDOB.size()+" individual"+(missingDOB.size()>1?"s":"")+"; see missingDOBs.xls for details.");
				writer = new PrintWriter(new FileWriter(trav+"_missingDOBs.xln"));
				writer.println("FamID\tIndID\tDNA");
				for (int i = 0; i<missingDOB.size(); i++) {
					writer.println(missingDOB.elementAt(i));
				}
				writer.close();
			}
		} catch (IOException ioe) {
			System.err.println("Error parsing "+filename);

		}

		return inds;
	}

	public static void procNinfoType2(String filename, String[] header, Vector<String> inds) {
		BufferedReader reader = null;
		PrintWriter writer_CSV, writer_DAT;
		String[] line;
		String trav;
		Hashtable<String,Vector<String[]>> hash = new Hashtable<String,Vector<String[]>>();
		Hashtable<String,String> hashString;

		try {
			try {
				reader = new BufferedReader(new FileReader(filename));
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+filename+"\" not found in current directory");
				System.exit(1);
			}
			trav = filename.substring(filename.lastIndexOf("/")+1, filename.lastIndexOf("."));
			writer_CSV = new PrintWriter(new FileWriter(trav+".csv"));
			writer_DAT = new PrintWriter(new FileWriter(trav+".dat"));

			line = reader.readLine().split("[\\s]+");
			ext.checkHeader(line, header, true);
			for (int i = 0; i<line.length; i++) {
				writer_CSV.print((i==0?"":",")+line[i]);
				writer_DAT.print((i==0?"":"\t")+line[i]);
			}
			writer_CSV.println();
			writer_DAT.println();

			hashString = new Hashtable<String,String>();
			while (reader.ready()) {
				line = reader.readLine().split("\t", -1);
				line[4] = line[4].equals("")?"0":line[4];
				line[5] = line[5].equals("")?"0":line[5];
				for (int i = 0; i<line.length; i++) {
					if (inds.contains(line[0]+"\t"+line[1])&&!hash.containsKey(line[0]+"\t"+line[1])) {
						writer_CSV.print((i==0?"":",")+(line[i].equals("")?".":line[i]));
					}
					writer_DAT.print((i==0?"":"\t")+(line[i].equals("")?".":line[i]));
				}
				for (int i = line.length; i<8; i++) {
					if (inds.contains(line[0]+"\t"+line[1])&&!hash.containsKey(line[0]+"\t"+line[1])) {
						writer_CSV.print(",.");
					}
					writer_DAT.print("\t.");
				}
				if (inds.contains(line[0]+"\t"+line[1])&&!hash.containsKey(line[0]+"\t"+line[1])) {
					writer_CSV.println();
					hashString.put(line[0]+"\t"+line[1], "");
				}
				writer_DAT.println();
			}
			reader.close();
			writer_CSV.close();
			writer_DAT.close();
		} catch (IOException ioe) {
			System.err.println("Error parsing "+filename);

		}
	}

	public static String procDOB(String[] line, int[] indices, Vector<String> missingDOB) {
		String str = "";

		for (int i = 0; i<2; i++) {
			if (line[indices[i]].equals("")) {
				str += "UU";
			} else if (Integer.parseInt(line[indices[i]])<1||Integer.parseInt(line[indices[i]])>(i==0?12:31)) {
				System.err.println("Error - '"+line[indices[i]]+"' is not a valid "+(i==0?"":"day of the ")+"month to be born into ("+line[0]+"-"+line[1]+")");
				str += "UU";
			} else {
				try {
					str += ext.formNum(Integer.parseInt(line[indices[i]]), 2);
				} catch (NumberFormatException nfe) {
					System.err.println("Error - '"+line[indices[i]]+"' is not a valid "+(i==0?"":"day of the ")+"month to be born into ("+line[0]+"-"+line[1]+")");
					str += "UU";
				}
			}
		}

		if (line[indices[2]].equals("")) {
			missingDOB.add(line[0]+"\t"+line[1]+"\t"+line[4]);
			str += "UUUU";
		} else if (line[indices[2]].length()==4&&Integer.parseInt(line[indices[2]])<1900||Integer.parseInt(line[indices[2]])>2000) {
			System.err.println("Error - '"+line[indices[2]]+"' is not a valid year to be born into ("+line[0]+"-"+line[1]+")");
			str += "UUUU";
		} else {
			try {
				str += Integer.parseInt(line[indices[2]]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - '"+line[indices[2]]+"' is not a valid year to be born into ("+line[0]+"-"+line[1]+")");
				str += "UUUU";
			}
		}

		return str;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = DEFAULT_DIR;
		String ninfo1 = "ninfo1.txt";
		String ninfo2 = "ninfo2.txt";
		String ninfo5 = "ninfo5.txt";
		String ninfo6 = "ninfo6.txt";

		String usage = "\n"+"park.ParseNinfos requires 0-2 arguments\n"+"   (1) the tab-delimited Nathan's-info-1 file (default: ninfo1="+ninfo1+")\n"+"   (2) the tab-delimited Nathan's-info-2 file (default: ninfo2="+ninfo2+")\n"+"   (2) the tab-delimited Cares-Nathan's Info-1 file (default: ninfo5="+ninfo5+")\n"+"   (2) the tab-delimited Cares-Nathan's-info-2 file (default: ninfo6="+ninfo6+")\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ninfo1=")) {
				ninfo1 = args[i].split("=")[1];
				if (!new File(ninfo1).exists()) {
					System.err.println("Error - file '"+ninfo1+"' does not exist");
					System.exit(2);
				}
				numArgs--;
			} else if (args[i].startsWith("ninfo2=")) {
				ninfo2 = args[i].split("=")[1];
				if (!new File(ninfo2).exists()) {
					System.err.println("Error - file '"+ninfo2+"' does not exist");
					System.exit(2);
				}
				numArgs--;
			} else if (args[i].startsWith("ninfo5=")) {
				ninfo5 = args[i].split("=")[1];
				if (!new File(ninfo5).exists()) {
					System.err.println("Error - file '"+ninfo5+"' does not exist");
					System.exit(2);
				}
				numArgs--;
			} else if (args[i].startsWith("ninfo6=")) {
				ninfo6 = args[i].split("=")[1];
				if (!new File(ninfo6).exists()) {
					System.err.println("Error - file '"+ninfo6+"' does not exist");
					System.exit(2);
				}
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.out.println("Warning - using defaults (processing data in "+ninfo1+", "+ninfo2+", "+ninfo5+" and "+ninfo6+")");
		}
		try {
			parseNinfos(dir, ninfo1, ninfo2, ninfo5, ninfo6);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
