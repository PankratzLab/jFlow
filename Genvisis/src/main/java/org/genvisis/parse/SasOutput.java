// for parseMeans add final -----------
package org.genvisis.parse;

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

public class SasOutput {
	public static final String[] HEADERS = {"Variable", "Label", "N", "Mean", "Std Dev", "Minimum",
																					"Maximum"};

	public static void parseMeans(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		Hashtable<String, String> hash = new Hashtable<String, String>();

		String type, stratum;
		boolean inMeans;
		String[] keys;

		type = "unknown type";
		stratum = "unknown stratum";
		inMeans = false;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + ".xln"));
			writer.println(Array.toStr(HEADERS) + "\tType\tStratum");
			while (reader.ready()) {
				temp = reader.readLine().trim();
				if (temp.startsWith("For all")) {
					line = temp.trim().split("-");
					line = line[1].split("[\\s]+");
					type = Array.toStr(Array.subArray(line, 0, line.length - 1), "_");
				}
				if (temp.startsWith("-----")) {
					line = temp.trim().split("[\\s]+");
					// System.err.println(temp);
					stratum = Array.toStr(Array.subArray(line, 1, line.length - 1), "_");

					keys = HashVec.getKeys(hash);
					for (String key : keys) {
						writer.println(hash.get(key));
					}

					hash = new Hashtable<String, String>();
				}
				if (temp.startsWith("******")) {
					inMeans = !inMeans;
				} else if (inMeans) {
					line = temp.trim().split("[\\s]+");
					if (hash.containsKey(line[0])) {
						temp = hash.get(line[0]) + "\t" + Array.toStr(line);
					}
					line = temp.trim().split("[\\s]+");
					hash.put(line[0], Array.toStr(line) + "\t" + type + "\t" + stratum);
				}

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
		String filename = "Gwar.dat";

		String usage = "\n"	+ "parse.SasOutput requires 0-1 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		filename = "D:\\means.txt";
		try {
			parseMeans(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
