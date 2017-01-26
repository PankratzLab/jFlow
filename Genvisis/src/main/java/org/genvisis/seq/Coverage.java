package org.genvisis.seq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Coverage {
	public static final String[] COVERAGE_HEADER = {"Chromosome", "Position", "Reference base",
																									"Consensus base",
																									"Percent of reads on the forward strand",
																									"Percent of reads on the reverse strand",
																									"Percent of 'A' nucleotides",
																									"Count of 'A' nucleotides",
																									"Percent of 'T' nucleotides",
																									"Count of 'T' nucleotides",
																									"Percent of 'C' nucleotides",
																									"Count of 'C' nucleotides",
																									"Percent of 'G' nucleotides",
																									"Count of 'G' nucleotides",
																									"Percent of 'N' nucleotides (no calls)",
																									"Count of 'N' nucleotides (no calls)",
																									"Count of matches to the reference",
																									"Count of mismatches to the reference",
																									"Percentage of matches to the reference",
																									"Consensus Quality", "SNP quality",
																									"Root mean square", "read depth",
																									"read base string", "base quality string"};

	public static void parse(String filename, String dir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		// String temp, trav;
		Hashtable<String, String> hash;
		Hashtable<String, Hashtable<String, String>> hashes;
		// Vector<String> v = new Vector<String>();
		// int count;
		// long time;
		String[] files, variants, inds;
		int fileIndex;

		files = Files.list(dir, ".txt", false);
		ext.checkHeader(Files.getHeaderOfFile(filename, "\t", new Logger()),
										new String[] {"Sample", "Chr", "Position"}, new int[] {0, 1, 2}, false,
										new Logger(), true);
		variants = HashVec.loadFileToStringArray(filename, true, new int[] {0, 1, 2}, false);

		hashes = new Hashtable<String, Hashtable<String, String>>();
		for (int i = 0; i < variants.length; i++) {
			line = variants[i].split("[\\s]+");
			HashVec.addToHashHash(hashes, line[0], "chr" + line[1] + "_" + line[2], i + "");
		}

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename) + "_list.xln"));
			writer.println("Sample\tIndex\t" + ArrayUtils.toStr(COVERAGE_HEADER));
			inds = HashVec.getKeys(hashes, false);
			for (int i = 0; i < inds.length; i++) {
				System.out.println("Analyzing indiviudal " + (i + 1) + " of " + inds.length);
				hash = hashes.get(inds[i]);
				fileIndex = -1;
				for (int j = 0; j < files.length; j++) {
					if (files[j].contains(inds[i])) {
						if (fileIndex == -1) {
							fileIndex = j;
						} else if (fileIndex == -2) {
							System.err.println("     " + files[j]);
						} else {
							System.err.println("Error - multiple files contain data for '" + inds[i] + "'");
							System.err.println("     " + files[fileIndex]);
							System.err.println("     " + files[j]);
							fileIndex = -2;
						}
					}
				}
				if (fileIndex == -1) {
					System.err.println("Error - no file match for sample '" + inds[i] + "'");
				} else {
					try {
						reader = new BufferedReader(new FileReader(dir + files[fileIndex]));
						while (reader.ready()) {
							line = reader.readLine().trim().split("[\\s]+");
							if (hash.containsKey(line[0] + "_" + line[1])) {
								writer.println(inds[i]	+ "\t" + hash.get(line[0] + "_" + line[1]) + "\t"
																+ ArrayUtils.toStr(line));
								writer.flush();
								hash.put(line[0] + "_" + line[1], "dup");
							}
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""	+ dir + files[fileIndex]
																+ "\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \"" + dir + files[fileIndex] + "\"");
						System.exit(2);
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename) + "_list.xln");
			e.printStackTrace();
		}
	}

	public static void filter(String filesToFilter) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Hashtable<String, String> hash;
		String varFile;
		String pileupFile;

		line = filesToFilter.split(",");
		pileupFile = line[0];
		varFile = line[1];

		hash = HashVec.loadFileToHashString(varFile, new int[] {0, 1}, null, false, "\t", false, false,
																				false);
		try {
			reader = new BufferedReader(new FileReader(pileupFile));
			writer = new PrintWriter(new FileWriter(pileupFile + "_filtered.out"));
			while (reader.ready()) {
				trav = reader.readLine();
				line = trav.trim().split("[\\s]+");
				if (hash.containsKey(line[0] + "\t" + line[1])) {
					writer.println(trav);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + pileupFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + pileupFile + "\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "variants.txt";
		String dir = "/data/parkinson/pd_sequence/coverage_files/";
		boolean parseCoverageFiles = false;
		String filter = "";

		String usage = "\n"	+ "seq.Coverage requires 0-1 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n" + "   (2) directory of coverage files (i.e. cov="
										+ dir + " (default))\n" + " OR\n"
										+ "   (1) files to filter (i.e. filter=fileToFilter.pileup,variants.txt (not the default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("filter=")) {
				filter = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (!filter.equals("")) {
				filter(filter);
			} else if (parseCoverageFiles) {
				parse(filename, dir);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
