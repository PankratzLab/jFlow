package org.genvisis.cnv.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

public class ReferenceMaps {
	public static final String MCCAROLL_DIR =
																					"C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\McCaroll map\\";

	public static final String MCCAROLL_MAP = "ng.238-S2.txt";

	public static final String[] MCCAROLL_HEADER = {"CNP_id", "chr", "start_hg17", "end_hg17",
																									"start_hg18", "end_hg18", "CN classes observed"};

	public static final String[] MCCAROLL_POPULATIONS = {"CEU", "YRI", "CHB+JPT"};

	public static final int MAX_COPIES = 6;

	public static void createMcCarollMap() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, keys;
		String temp;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Hashtable<String, String> inds;
		int count;
		int[][][] counts;
		boolean[] founders;
		int index;

		try {
			reader = new BufferedReader(new FileReader(MCCAROLL_DIR + MCCAROLL_MAP));
			do {
				temp = reader.readLine();
			} while (!temp.startsWith("CNP_id"));
			ext.checkHeader(temp.trim().split("\t", -1), MCCAROLL_HEADER, true);
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[0], count	+ "\t" + line[1] + ":" + ext.replaceAllWith(line[4], "\"", "") + "-"
													+ ext.replaceAllWith(line[5], "\"", ""));
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ MCCAROLL_DIR + MCCAROLL_MAP
													+ "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + MCCAROLL_DIR + MCCAROLL_MAP + "\"");
			System.exit(2);
		}

		keys = HashVec.getKeys(hash);
		counts = new int[keys.length][MCCAROLL_POPULATIONS.length][MAX_COPIES + 2];
		for (int i = 0; i < MCCAROLL_POPULATIONS.length; i++) {
			inds = new Hashtable<String, String>();
			try {
				reader = new BufferedReader(new FileReader(MCCAROLL_DIR	+ "pedinfo2sample_"
																										+ MCCAROLL_POPULATIONS[i] + ".txt"));
				while (reader.ready()) {
					line = reader.readLine().split("\t", -1);
					if (line[2].equals("0") && line[3].equals("0")) {
						inds.put(line[6].trim().split(":")[4], "");
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""	+ MCCAROLL_DIR + "pedinfo2sample_"
														+ MCCAROLL_POPULATIONS[i] + ".txt"
														+ "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""	+ MCCAROLL_DIR + "pedinfo2sample_"
														+ MCCAROLL_POPULATIONS[i] + ".txt" + "\"");
				System.exit(2);
			}

			try {
				reader = new BufferedReader(new FileReader(MCCAROLL_DIR	+ "ng.238-S3_"
																										+ MCCAROLL_POPULATIONS[i] + ".txt"));
				line = reader.readLine().split("[\\s]+");
				if (!line[0].equals("CNP_id")) {
					System.err.println("Error - not the header I was expecting");
					System.exit(1);
				}
				founders = new boolean[line.length - 1];
				for (int j = 0; j < founders.length; j++) {
					founders[j] = inds.containsKey(line[j + 1]);
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					index = Integer.parseInt(hash.get(line[0]).split("[\\s]+")[0]);
					for (int j = 0; j < founders.length; j++) {
						if (founders[j] && !line[j + 1].equals("NA")) {
							counts[index][i][Integer.parseInt(line[j + 1])]++;
							counts[index][i][MAX_COPIES + 1]++;
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""	+ MCCAROLL_DIR + "ng.238-S3_" + MCCAROLL_POPULATIONS[i]
														+ ".txt" + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""	+ MCCAROLL_DIR + "ng.238-S3_"
														+ MCCAROLL_POPULATIONS[i] + ".txt" + "\"");
				System.exit(2);
			}

		}

		try {
			writer = new PrintWriter(new FileWriter(MCCAROLL_DIR + "file.xln"));
			writer.print("CNP_id\tloc\tindex");
			for (String element : MCCAROLL_POPULATIONS) {
				for (int j = 0; j <= MAX_COPIES; j++) {
					writer.print("\t" + element + "_" + j);
				}
			}
			writer.println();
			for (int i = 0; i < keys.length; i++) {
				writer.print(keys[i] + "\t" + hash.get(keys[i]));
				for (int j = 0; j < MCCAROLL_POPULATIONS.length; j++) {
					for (int k = 0; k <= MAX_COPIES; k++) {
						writer.print("\t" + ((double) counts[i][j][k] / (double) counts[i][j][MAX_COPIES + 1]));
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to file");
			e.printStackTrace();
		}
	}

	public static void generateBed(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		CNVariant cnv;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + ".cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			while (reader.ready()) {
				cnv = new CNVariant(reader.readLine());
				writer.println(cnv.toPlinkFormat());
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename =
										"C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\McCaroll map\\all_CNPs.txt";

		String usage = "\n"	+ "cnv.ReferenceMaps requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			generateBed(filename);
			// createMcCarollMap();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
