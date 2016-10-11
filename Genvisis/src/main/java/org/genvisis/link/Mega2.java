package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class Mega2 {

	public static void createFiles(int start, int stop) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		Vector<String> v = new Vector<String>();
		int numMarkers;
		String chrome;
		double sumPos;

		for (int chr = start; chr <= stop; chr++) {
			chrome = ext.chrome(chr);
			try {
				reader = new BufferedReader(new FileReader("map" + chrome + ".dat"));
				writer = new PrintWriter(new FileWriter("map." + chrome));
				writer.println("CHROMOSOME\tKOSAMBI\tNAME");

				numMarkers = Integer.parseInt(reader.readLine().trim().split("[\\s]+")[0]) - 1;
				reader.readLine();
				reader.readLine();
				if (!reader.readLine().contains("#")) {
					System.err.println("Error - Line 4 (underlying disease allele) must by named (i.e. \"# TRAIT\") or else Mega2 will fail");

				}
				for (int i = 0; i < (chr == 23 ? 4 : 3); i++) {
					reader.readLine();
				}
				for (int i = 0; i < numMarkers; i++) {
					temp = reader.readLine();
					try {
						v.add(temp.trim().split("[\\s]+")[3]);
					} catch (Exception e) {
						System.err.println("Error parsing line: " + temp);
						e.printStackTrace();
					}
					reader.readLine();
				}
				reader.readLine();
				line = reader.readLine().trim().split("[\\s]+");
				sumPos = -1 * Double.parseDouble(line[0]);
				for (int i = 0; i < numMarkers; i++) {
					sumPos += Double.parseDouble(line[i]);
					writer.println(chr + "\t" + sumPos + "\t" + v.elementAt(i));
				}
				reader.close();
				writer.close();

				Files.copyFile("re_chrom" + chrome + ".pre", "pedin." + chrome);
				Files.copyFile("map" + chrome + ".dat", "datain." + chrome);
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""	+ "map" + chrome + ".dat"
														+ "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + "map" + chrome + ".dat" + "\"");
				System.exit(2);
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		int start = 7;
		int stop = 7;

		String usage = "\\n"	+ "park.Mega2 requires 0-1 arguments\n"
										+ "   (1) start chromosome (i.e. start=" + start + " (default))\n"
										+ "   (2) stop chromosome (i.e. stop=" + stop + " (default))\n" + " OR\n"
										+ "   (1) chromsome # to do (i.e. chr=" + start + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("start=")) {
				start = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("stop=")) {
				stop = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("chr=")) {
				start = stop = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			createFiles(start, stop);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
