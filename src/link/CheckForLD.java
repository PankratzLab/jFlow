package link;

import java.io.*;
import java.util.*;

import common.*;
import filesys.LDdatabase;
import bioinformatics.HapMapParser;
import bioinformatics.ParseSNPlocations;
import link.LinkageToPlink;

public class CheckForLD {
	public static final String DBSNP_SOURCE = ParseSNPlocations.DEFAULT_B36_SOURCE_FILENAME;
	public static final String DBSNP_LOCAL = "6K_b129.bcp";
	public static final String[] LD_HEADER = {"L1", "L2", "D'", "LOD", "r^2", "CIlow", "CIhi", "Dist", "T-int"};
	public static final String[] CHECK_HEADER = {"#", "Name", "Position", "ObsHET", "PredHET", "HWpval", "%Geno", "FamTrio", "MendErr", "MAF", "Alleles", "Rating"};
	public static final double DEFAULT_MAX_DPRIME = 0.7;
	public static final double DEFAULT_R2 = 0.3;

	public static void createLD(String dir, String checkDir, boolean lastNotFirst) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, keys, data;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		Hashtable<String,String[]> markersByChrome = new Hashtable<String,String[]>();
		Hashtable<String,String> markerPositions = new Hashtable<String,String>();
		int start = 1;
		int stop = 22;

		new File(dir+checkDir).mkdirs();
		for (int i = start; i<=stop; i++) {
			hash.clear();
			try {
				// reader = new BufferedReader(new
				// FileReader(root+"mrkr"+ext.chrome(i)+".dat"));
				reader = new BufferedReader(new FileReader(dir+"re_chrom"+ext.chrome(i)+".pre"));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					// data = Array.subArray(line, 2);
					data = Array.subArray(line, 6);
					if (Array.sum(Array.toIntArray(data))>0&&(lastNotFirst||!hash.containsKey(line[0]))) {
						hash.put(line[0], line[0]+"\t"+line[1]+"\t0\t0\t1\t2\t"+Array.toStr(data));
					}
				}
				reader.close();

				keys = HashVec.getKeys(hash, true, true);
				writer = new PrintWriter(new FileWriter(dir+checkDir+"check"+ext.chrome(i)+".pre"));
				for (int j = 0; j<keys.length; j++) {
					writer.println(hash.get(keys[j]));
				}
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+dir+"re_chrom"+ext.chrome(i)+".pre"+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+dir+"re_chrom"+ext.chrome(i)+".pre"+"\"");
				System.exit(2);
			}
		}

		for (int i = start; i<=stop; i++) {
			markersByChrome.put(i+"", data = new LinkageMap(dir+"map"+ext.chrome(i)+".dat").getMarkerNames());
			for (int j = 0; j<data.length; j++) {
				markerPositions.put(data[j], "-1");
			}
		}

		if (!new File(dir+DBSNP_LOCAL).exists()) {
			parseLocalDBSNP(DBSNP_SOURCE, markerPositions, dir+DBSNP_LOCAL);
		}

		try {
			reader = new BufferedReader(new FileReader(dir+DBSNP_LOCAL));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (markerPositions.containsKey(line[0])) {
					markerPositions.put(line[0], line[1]+"\t"+line[2]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+DBSNP_SOURCE+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+DBSNP_SOURCE+"\"");
			System.exit(2);
		}

		for (int i = start; i<=stop; i++) {
			data = markersByChrome.get(i+"");
			try {
				writer = new PrintWriter(new FileWriter(dir+checkDir+"check"+ext.chrome(i)+".info"));
				for (int j = 0; j<data.length; j++) {
					line = markerPositions.get(data[j]).split("[\\s]+");
					if (line[0].equals("-1")) {
						System.err.println("Error - '"+data[j]+"' is supposed to be on chromosome "+i+", but was not found in the dbSNP database");
						writer.println(data[j]+"\t0");
					} else if (!line[0].equals(i+"")) {
						System.err.println("Error - '"+data[j]+"' was supposed to be on chromosome "+i+", but the dbSNP database places it on chr "+line[0]);
						writer.println(data[j]+"\t0");
					} else {
						writer.println(data[j]+"\t"+line[1]);
					}
					writer.flush();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to "+dir+checkDir+"check"+ext.chrome(i)+".info");
				e.printStackTrace();
			}
		}
	}

	public static void createHapMap(String root, String checkDir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, data;
		Hashtable<String,String[]> markersByChrome = new Hashtable<String,String[]>();
		Hashtable<String,String> markerPositions = new Hashtable<String,String>();
		int start = 1;
		int stop = 22;
		String chrome;

		new File(root+checkDir).mkdirs();
		for (int i = start; i<=stop; i++) {
			markersByChrome.put(i+"", data = new LinkageMap(root+"map"+ext.chrome(i)+".dat").getMarkerNames());
			for (int j = 0; j<data.length; j++) {
				markerPositions.put(data[j], "-1");
			}
		}

		if (!new File(root+DBSNP_LOCAL).exists()) {
			parseLocalDBSNP(DBSNP_SOURCE, markerPositions, root+DBSNP_LOCAL);
		}

		try {
			reader = new BufferedReader(new FileReader(root+DBSNP_LOCAL));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (markerPositions.containsKey(line[0])) {
					markerPositions.put(line[0], line[1]+"\t"+line[2]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+DBSNP_SOURCE+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+DBSNP_SOURCE+"\"");
			System.exit(2);
		}

		for (int i = start; i<=stop; i++) {
			data = markersByChrome.get(i+"");
			chrome = ext.chrome(i);
			try {
				writer = new PrintWriter(new FileWriter(root+checkDir+"hapmap"+chrome+".info"));
				for (int j = 0; j<data.length; j++) {
					line = markerPositions.get(data[j]).split("[\\s]+");
					if (line[0].equals("-1")) {
						System.err.println("Error - '"+data[j]+"' is supposed to be on chromosome "+i+", but was not found in the dbSNP database");
						writer.println(data[j]+"\t0");
					} else if (!line[0].equals(i+"")) {
						System.err.println("Error - '"+data[j]+"' was supposed to be on chromosome "+i+", but the dbSNP database places it on chr "+line[0]);
						writer.println(data[j]+"\t0");
					} else {
						writer.println(data[j]+"\t"+line[1]);
					}
					writer.flush();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to "+root+checkDir+"hapmap"+chrome+".info");
				e.printStackTrace();
			}

			// LDdatabase.MASTER_HAPMAP_ROOT used to be hard coded as /home/npankrat/NCBI/HapMap/hapmap-ceu-chr"+i
			CmdLine.run("plink --bfile "+LDdatabase.MASTER_HAPMAP_ROOT+" --extract hapmap"+chrome+".info --missing-phenotype 0 --recode --out hapmap"+chrome, root+checkDir);
			HapMapParser.plinkMapToHaploviewInfo(root+checkDir+"hapmap"+chrome+".map", root+checkDir+"hapmap"+chrome+".info");
			new File(root+checkDir+"hapmap"+chrome+".pre").delete();
			new File(root+checkDir+"hapmap"+chrome+".ped").renameTo(new File(root+checkDir+"hapmap"+chrome+".pre"));
			new File(root+checkDir+"hapmap"+chrome+".map").delete();
			new File(root+checkDir+"hapmap"+chrome+".log").delete();
			new File(root+checkDir+".pversion").delete();

		}
	}

	public static void parseLocalDBSNP(String source, Hashtable<String,String> markerPositions, String fileout) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		long time;

		try {
			System.out.println("Creating a faster local copy of the SNP positions found in "+DBSNP_SOURCE);
			time = new Date().getTime();
			reader = Files.getAppropriateReader(DBSNP_SOURCE);
			writer = new PrintWriter(new FileWriter(fileout));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (markerPositions.containsKey("rs"+line[0])) {
					try {
						writer.println("rs"+line[0]+"\t"+line[1]+"\t"+(Integer.parseInt(line[2])+1));
					} catch (Exception e) {
						// System.err.println(Array.toStr(line));
						// e.printStackTrace();
					}
				}
			}
			reader.close();
			writer.close();
			System.out.println("Finished in "+ext.getTimeElapsed(time));
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+DBSNP_SOURCE+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+DBSNP_SOURCE+"\"");
			System.exit(2);
		}

	}

	public static void checkLD(String root, String checkDir, String prefix) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(root+checkDir+prefix+"_haplo.bat"));
			for (int chr = 1; chr<=22; chr++) {
				writer.println("java -jar Haploview.jar -nogui -pedfile "+prefix+ext.chrome(chr)+".pre -info "+prefix+ext.chrome(chr)+".info -check -dprime");
			}
			writer.close();

		} catch (Exception e) {
			System.err.println("Error writing batch for Haploview");
			e.printStackTrace();
		}
		System.out.println("Don't forget to copy over Haploview.jar");
	}

	public static void parseLD(String root, String checkDir, String prefix, String hapmapDir, String hapmapPrefix, double maxDprime, double maxr2) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, subline;
		Hashtable<String,String> maxObsLD, maxHapmapLD, hapCheck;
		Vector<String> v = new Vector<String>();

		try {
			writer = new PrintWriter(new FileWriter(root+checkDir+prefix+"_summary.xln"));
			writer.println("SNP\tChr\tPosition\tMAF\tObsHET\tPredHET\tHWpval\tMax LD marker\tD'\tr2\tHapMap MAF\tHapMap HW\tHapMap Max LD marker\tHapMap D'\tHapMap r2");

			for (int chr = 1; chr<=22; chr++) {
				maxObsLD = parseMaxLD(root+checkDir+prefix+ext.chrome(chr)+".pre.LD");
				maxHapmapLD = parseMaxLD(root+hapmapDir+hapmapPrefix+ext.chrome(chr)+".pre.LD");
				hapCheck = new Hashtable<String,String>();
				try {
					reader = new BufferedReader(new FileReader(root+hapmapDir+hapmapPrefix+ext.chrome(chr)+".pre.CHECK"));
					ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CHECK_HEADER, true);
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						hapCheck.put(line[1], line[9]+"\t"+line[5]);
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+root+hapmapDir+hapmapPrefix+ext.chrome(chr)+".pre.CHECK"+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+root+hapmapDir+hapmapPrefix+ext.chrome(chr)+".pre.CHECK"+"\"");
					System.exit(2);
				}
				try {
					reader = new BufferedReader(new FileReader(root+checkDir+prefix+ext.chrome(chr)+".pre.CHECK"));
					ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CHECK_HEADER, true);
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						writer.println(line[1]+"\t"+chr+"\t"+line[2]+"\t"+line[9]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+(maxObsLD.containsKey(line[1])?maxObsLD.get(line[1]):".\t.\t.")+"\t"+(hapCheck.containsKey(line[1])?hapCheck.get(line[1]):".\t.")+"\t"+(maxHapmapLD.containsKey(line[1])?maxHapmapLD.get(line[1]):".\t.\t.")+"\t");
						if (maxObsLD.containsKey(line[1])) {
							subline = maxObsLD.get(line[1]).split("[\\s]+");
							if (Double.parseDouble(subline[1])>=maxDprime||Double.parseDouble(subline[2])>=maxr2) {
								v.add(line[1]+"\t"+Array.toStr(subline));
							}
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+root+checkDir+prefix+ext.chrome(chr)+".pre.CHECK"+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+root+checkDir+prefix+ext.chrome(chr)+".pre.CHECK"+"\"");
					System.exit(2);
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to "+root+checkDir+prefix+"_summary.xln");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(root+"dropList.dat"));
			for (int i = 0; i<v.size(); i++) {
				writer.println(v.elementAt(i));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to "+root+checkDir+prefix+"_summary.xln");
			e.printStackTrace();
		}

	}

	public static Hashtable<String,String> parseMaxLD(String filename) {
		BufferedReader reader;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		try {
			reader = new BufferedReader(new FileReader(filename));
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), LD_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!hash.containsKey(line[0])||Double.parseDouble(line[4])>Double.parseDouble(hash.get(line[0]).split("[\\s]+")[2])) {
					hash.put(line[0], line[1]+"\t"+line[2]+"\t"+line[4]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hash;
	}

	public static void plinkMethod(String dir, boolean vif) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] markerNames;

		if (!new File(dir+"plink.bed").exists()) {
			for (int i = 1; i<=23; i++) {
				if (new File(dir+"map"+ext.chrome(i)+".dat").exists()) {
					markerNames = new LinkageMap(dir+"map"+ext.chrome(i)+".dat").getMarkerNames();
					for (int j = 0; j<markerNames.length; j++) {
						hash.put(markerNames[j], "-1");
					}
				} else {
					System.err.println("skipping chromosome "+i+" (map"+ext.chrome(i)+".dat not found)");
				}
			}

			if (!new File(dir+DBSNP_LOCAL).exists()) {
				parseLocalDBSNP(DBSNP_SOURCE, hash, dir+DBSNP_LOCAL);
			}

			LinkageToPlink.convert(dir, dir+DBSNP_LOCAL);
		}
		System.out.println("Running plink's LD based SNP pruning method:");
		if (vif) {
			System.out.println("plink --bfile plink --indep 50 5 2");
			CmdLine.run("plink --bfile plink --indep 50 5 2", dir);
		} else {
			System.out.println("plink --bfile plink --indep-pairwise 50 5 0.3");
			CmdLine.run("plink --bfile plink --indep-pairwise 50 5 0.3", dir);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "CheckForLD.dat";
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\";
		// String dir = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\replicate\\";

		boolean createLD = false;
		boolean checkLD = false;
		boolean parseLD = false;
		boolean plinkMethod = true;
		String checkDir = "check/";
		// String hapmap = "hapmap/";
		double maxDprime = DEFAULT_MAX_DPRIME;
		double maxr2 = DEFAULT_R2;
		boolean vif = false;

		String usage = "\\n"+"link.CheckForLD requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (!new File(dir).exists()) {
			System.err.println("Error - using current directory instead of missing directory ("+dir+")");
			dir = "";
		}
		try {

			parseLocalDBSNP(DBSNP_SOURCE, HashVec.loadFileToHashString(dir+"hisplusours.txt", false), dir+DBSNP_LOCAL);
			System.exit(1);

			if (createLD) {
				// createLD(dir, checkDir, false);
				createLD(dir, checkDir, true);
				createHapMap(dir, checkDir);
			}
			if (checkLD) {
				checkLD(dir, checkDir, "check");
				checkLD(dir, checkDir, "hapmap");
			}
			if (parseLD) {
				parseLD(dir, checkDir, "check", checkDir, "hapmap", maxDprime, maxr2);
			}
			if (plinkMethod) {
				plinkMethod(dir, vif);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
