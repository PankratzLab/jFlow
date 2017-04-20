package org.genvisis.db;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;

public class DummyDataset {
	public static void createFromParameters(String filename, Logger log) {
		PrintWriter writer;
		String trav;
		String[] line, header;
		Vector<String[]> v;
		IntVector countV;
		int num;
		String outfile;
		int count;
		String[][] params;

		params = Files.parseControlFile(filename, true, "dummy",
																		new String[] {"outputFile.xln",
																									"header\tAffected\tAllele1\tAllele2",
																									"30\t1\tE2\tE3", "40\t0\tE2\tE3",
																									"500\t1\tE3\tE3", "400\t0\tE3\tE3",
																									"200\t1\tE3\tE4", "100\t0\tE3\tE4",
																									"10\t1\tE4\tE4", "2\t0\tE4\tE4"},
																		log);
		if (params != null) {
			outfile = params[0][0];
			header = null;
			v = new Vector<String[]>();
			countV = new IntVector();
			for (int i = 1; i < params.length; i++) {
				if (params[i][0].equalsIgnoreCase("header")) {
					header = ArrayUtils.subArray(params[i], 1);
				} else if (!params[i][0].startsWith("#") && !params[i][0].equals("")) {
					try {
						num = Integer.parseInt(params[i][0]);
						countV.add(num);
					} catch (Exception e) {
						log.reportError("Error - '" + params[i][0]
														+ "' is an invalid integer; first column contains the number of time series is repreated");
						return;
					}
					v.add(ArrayUtils.subArray(params[i], 1));
				}
			}

			try {
				count = 0;
				writer = Files.openAppropriateWriter(outfile);
				if (header != null) {
					writer.println(ArrayUtils.toStr(header));
					for (int i = 0; i < v.size(); i++) {
						num = countV.elementAt(i);
						line = v.elementAt(i);
						trav = ArrayUtils.toStr(line);
						for (int j = 0; j < num; j++) {
							writer.println(trav);
							count++;
						}
					}
				}
				writer.close();
				log.report("Final file contains " + count + " records"
									 + (header == null ? "" : "; plus a header"));
			} catch (Exception e) {
				System.err.println("Error writing to " + outfile);
				e.printStackTrace();
			}
		}
	}

	public static void createReverseFromParameters(String filename, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, keys;
		Vector<String> v;
		String infile, outfile;
		boolean commaDelimited, tabDelimited;
		Hashtable<String, Hashtable<String, String>> hashes;
		Hashtable<String, String> hash;
		String[][] params;

		params = Files.parseControlFile(filename, false, "counts",
																		new String[] {"input.txt out=outputFile.xln tab"}, log);
		if (params != null) {
			infile = params[0][0];
			outfile = null;
			commaDelimited = false;
			tabDelimited = false;
			for (int i = 1; i < params[0].length; i++) {
				if (params[0][i].startsWith("out=")) {
					outfile = params[0][i].split("=")[1];
				} else if (params[0][i].equals(",")) {
					commaDelimited = true;
				} else if (params[0][i].equals("tab")) {
					tabDelimited = true;
				}
			}

			hashes = new Hashtable<String, Hashtable<String, String>>();
			v = new Vector<String>();
			try {
				reader = new BufferedReader(new FileReader(infile));
				while (reader.ready()) {
					line = reader.readLine().trim()
											 .split(commaDelimited ? ","
																						 : (tabDelimited ? "\t" : PSF.Regex.GREEDY_WHITESPACE));
					if (hashes.containsKey(line[0])) {
						hash = hashes.get(line[0]);
					} else {
						hashes.put(line[0], hash = new Hashtable<String, String>());
					}
					if (hash.containsKey(line[1])) {
						hash.put(line[1], (Integer.parseInt(hash.get(line[1])) + 1) + "");
					} else {
						hash.put(line[1], "1");
						HashVec.addIfAbsent(line[1], v);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + infile + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + infile + "\"");
				System.exit(2);
			}

			try {
				writer = Files.openAppropriateWriter(outfile);
				line = ArrayUtils.toStringArray(v);
				writer.println("Unit\t" + ArrayUtils.toStr(line));
				keys = HashVec.getKeys(hashes);
				for (String key : keys) {
					writer.print(key);
					hash = hashes.get(key);
					for (String element : line) {
						writer.print("\t" + (hash.containsKey(element) ? hash.get(element) : "0"));
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + outfile);
				e.printStackTrace();
			}
		}
	}
}
