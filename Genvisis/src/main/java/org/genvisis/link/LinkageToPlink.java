package link;

import java.io.*;
import java.util.*;
import common.*;

public class LinkageToPlink {
	// public static final String DEFAULT_DIRECTORY = "";
	public static final String DEFAULT_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\";

	public static final String DEFAULT_DATABASE = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\local_b129.bcp";

	public static void convert(String dir, String snp_database) {
		PrintWriter writer;
		Vector<String[]> v = new Vector<String[]>();
		String chrome;

		System.out.println("Creating files...");
		for (int chr = 1; chr<=23; chr++) {
			chrome = ext.chrome(chr);
			if (new File(dir+"map"+chrome+".dat").exists()) {
				new LinkageMap(dir, chr).createPlinkMap(dir+"plink"+chrome+".map", snp_database);
				v.add(new String[] {"re_chrom"+chrome+".pre", "plink"+chrome+".map"});
			} else {
				System.err.println("skipping chromosome "+chr+" (map"+chrome+".dat not found)");
			}
		}

		if (v.size()>1) {
			try {
				writer = new PrintWriter(new FileWriter(dir+"plink-merge.txt"));
				for (int i = 1; i<v.size(); i++) {
					writer.println(Array.toStr(v.elementAt(i)));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error creating merge list for plink");
				e.printStackTrace();
			}

			System.out.println("Merging files...");
			System.out.println("plink --ped "+v.elementAt(0)[0]+" --map "+v.elementAt(0)[1]+" --merge-list plink-merge.txt --make-bed");
			CmdLine.run("plink --ped "+v.elementAt(0)[0]+" --map "+v.elementAt(0)[1]+" --merge-list plink-merge.txt --make-bed", dir);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = DEFAULT_DIRECTORY;
		String db = DEFAULT_DATABASE;

		String usage = "\\n"+"link.LinkageToPlink requires 0-1 arguments\n"+"   (1) directory (i.e. dir="+dir+" (default))\n"+"   (2) database of SNP chromosomes and positions (i.e. db=null or db=filename (default is db="+db+"))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("db=")) {
				db = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			convert(dir, db);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
