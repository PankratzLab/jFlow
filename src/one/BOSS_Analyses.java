package one;

import gwas.CreateDatabaseFromPlink;

import java.io.*;
//import java.util.*;
//import common.*;

import common.CmdLine;

public class BOSS_Analyses {
	public static void generate(String dir) {
		new File(dir+"split/").mkdirs();
		for (int chr = 1; chr <= 23; chr++) {
//			CmdLine.run("plink --bfile ../plink --chr "+chr+" --recode --maf 0.001 --geno 0.2 --out chr"+chr, dir+"split/");
			CmdLine.run("plink --file chr"+chr+" --freq --out chr"+chr, dir+"split/");
		}

		new File(dir+"gwaf/").mkdirs();
		for (int chr = 1; chr <= 23; chr++) {
			CreateDatabaseFromPlink.toGWAF(dir+"split/chr"+chr+".ped", dir+"split/chr"+chr+".map", dir+"split/chr"+chr+".frq", dir+"gwaf/chr"+chr+".fhsR");
		}
	}

	public static void batch(String dir) {
		
	}
	
	public static void standardize(String dir) {
		BufferedReader reader, mapReader, infoReader;
		PrintWriter writer, mapWriter, infoWriter;
		String header, root;
		int count, fileNum;
	
		root = "gwaf/chr";
	
		count = 0;
		fileNum = 0;
		new File(dir+"std/").mkdirs();
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
						writer.println(reader.readLine());
						mapWriter.println(mapReader.readLine());
						infoWriter.println(infoReader.readLine());
						count++;
						if (count == 10000) {
							writer.close();
							mapWriter.close();
							infoWriter.close();
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
					System.err.println("Error: file \"" + dir+"leslie_lange.FHS.IBC.CEU.chr1.gen" + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"leslie_lange.FHS.IBC.CEU.chr1.gen" + "\"");
					System.exit(2);
				}
			}
			writer.close();
			mapWriter.close();
		} catch (Exception e) {
			System.err.println("Error writing to file #" + fileNum);
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "BOSS_Analyses.dat";
	
		String usage = "\n" + 
		"one.BOSS_Analyses requires 0-1 arguments\n" + 
		"   (1) directory (i.e. dir=" + dir + " (default))\n" + 
		"";
	
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		dir = "D:\\BOSS\\GWAF\\";
		
		try {
			generate(dir);
//			standardize(dir);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
