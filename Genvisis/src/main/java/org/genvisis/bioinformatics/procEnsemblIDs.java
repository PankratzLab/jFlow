package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

public class procEnsemblIDs {
	public static final String DEFAULT_GENES = "genes.xls";

	// public static final String DEFAULT_Q = "gene_seq.q";
	public static final String DEFAULT_Q = "C:\\Download\\seq_gene.q.shortcut.prn";

	public static final int[] COLS = {0, 1, 2};

	public static void runProcEnsemblIDs(String filename, String geneIDs) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String temp;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Hashtable<String, String> genes = new Hashtable<String, String>();

		try {
			reader = new BufferedReader(new FileReader(geneIDs));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				genes.put(line[1], line[0]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + geneIDs + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + geneIDs + "\"");
			System.exit(2);
		}

		int c1 = 0, c2 = 0;
		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			writer = new PrintWriter(new FileWriter(filename + ".out"));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				if (!line[0].equals("")) {
					line[0] = line[0].substring(7);
					if (hash.containsKey(line[1])) {
						if (line[0].equals(hash.get(line[1]))) {
							c1++;
						} else {
							if (genes.containsKey(line[1])) {
								writer.println(line[1]	+ "\t" + line[0] + "\t" + hash.get(line[1]) + "\t"
																+ genes.get(line[1]) + "\t" + genes.get(line[1]));
								hash.put(line[1], genes.get(line[1]));
							} else {
								if (Integer.parseInt(line[0]) < Integer.parseInt(hash.get(line[1]))) {
									writer.println(line[1]	+ "\t" + line[0] + "\t" + hash.get(line[1]) + "\t" + "#N/A"
																	+ "\t" + line[0]);
									hash.put(line[1], line[0]);
								} else {
									writer.println(line[1]	+ "\t" + line[0] + "\t" + hash.get(line[1]) + "\t" + "#N/A"
																	+ "\t" + hash.get(line[1]));
								}
							}
							c2++;
						}
					} else {
						hash.put(line[1], line[0]);
					}
				}
			}
			reader.close();
			System.out.println("Agrees " + c1);
			System.out.println("Disagrees " + c2);
			writer.close();

			writer = new PrintWriter(new FileWriter("Ensembl_lookup.xls"));
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			genes.clear();
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split("\t", -1);
				if (line[2].startsWith("ENSG") && hash.containsKey(line[1])) {
					writer.println(line[2] + "\t" + hash.get(line[1]));
				}
			}
			writer.close();

		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = DEFAULT_Q;
		String geneIDs = DEFAULT_GENES;

		String usage = "\n"	+ "park.procEnsemblIDs requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n"
										+ "   (2) filename with geneID/geneName in first two columns (i.e. genes="
										+ geneIDs + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("genes=")) {
				geneIDs = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			runProcEnsemblIDs(filename, geneIDs);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
