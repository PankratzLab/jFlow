package one;

import java.io.*;

//import java.util.*;
//import CommonTools.*;

public class Gaw16 {
	public static final String RAW_GENOTYPES = "narac.csv";

	public static final String RAW_MAP = "narac.map";

	public static void createPlink(String genotypes, String mapfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;

		try {
			reader = new BufferedReader(new FileReader(genotypes));
			writer = new PrintWriter(new FileWriter("plink.ped"));
			reader.readLine();
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split(",");
				writer.print(++count+"\t"+line[0]+"\t0\t0\t"+(line[2].equals("F")?"2":(line[2].equals("M")?"1":"?"))+"\t"+(line[1].equals("1")?"2":(line[1].equals("0")?"1":"?")));
				for (int i = 9; i<line.length; i++) {
					if (line[i].contains("?")) {
						writer.print("\t0\t0");
					} else {
						writer.print("\t"+line[i].charAt(0)+"\t"+line[i].charAt(2));
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+genotypes+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+genotypes+"\"");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(mapfile));
			writer = new PrintWriter(new FileWriter("plink.map"));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(line[1]+"\t"+line[0]+"\t0\t"+line[2]);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+mapfile+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+mapfile+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String genotypes = RAW_GENOTYPES;
		String mapfile = RAW_MAP;

		String usage = "\\n"+"park.oneoff.Gaw16 requires 0-1 arguments\n"+"   (1) genotypes file (i.e. geno="+genotypes+" (default))\n"+"   (2) map file (i.e. map="+mapfile+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("geno=")) {
				genotypes = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("map=")) {
				mapfile = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (!new File("plink.ped").exists()&&!new File("plink.bed").exists()) {
				createPlink(genotypes, mapfile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
