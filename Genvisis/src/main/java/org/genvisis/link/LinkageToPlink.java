package org.genvisis.link;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.ext;

public class LinkageToPlink {
	// public static final String DEFAULT_DIRECTORY = "";
	public static final String DEFAULT_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\";

	public static final String DEFAULT_DATABASE =
																							"C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\local_b129.bcp";

	public static void convert(String dir, String snp_database) {
		PrintWriter writer;
		Vector<String[]> v = new Vector<String[]>();
		String chrome;

		System.out.println("Creating files...");
		for (int chr = 1; chr <= 23; chr++) {
			chrome = ext.chrome(chr);
			if (new File(dir + "map" + chrome + ".dat").exists()) {
				new LinkageMap(dir, chr).createPlinkMap(dir + "plink" + chrome + ".map", snp_database);
				v.add(new String[] {"re_chrom" + chrome + ".pre", "plink" + chrome + ".map"});
			} else {
				System.err.println("skipping chromosome " + chr + " (map" + chrome + ".dat not found)");
			}
		}

		if (v.size() > 1) {
			try {
				writer = new PrintWriter(new FileWriter(dir + "plink-merge.txt"));
				for (int i = 1; i < v.size(); i++) {
					writer.println(ArrayUtils.toStr(v.elementAt(i)));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error creating merge list for plink");
				e.printStackTrace();
			}

			System.out.println("Merging files...");
			System.out.println("plink --ped "	+ v.elementAt(0)[0] + " --map " + v.elementAt(0)[1]
													+ " --merge-list plink-merge.txt --make-bed");
			CmdLine.run("plink --ped "	+ v.elementAt(0)[0] + " --map " + v.elementAt(0)[1]
									+ " --merge-list plink-merge.txt --make-bed", dir);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = DEFAULT_DIRECTORY;
		String db = DEFAULT_DATABASE;

		String usage = "\\n"	+ "link.LinkageToPlink requires 0-1 arguments\n"
										+ "   (1) directory (i.e. dir=" + dir + " (default))\n"
										+ "   (2) database of SNP chromosomes and positions (i.e. db=null or db=filename (default is db="
										+ db + "))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("db=")) {
				db = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
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
