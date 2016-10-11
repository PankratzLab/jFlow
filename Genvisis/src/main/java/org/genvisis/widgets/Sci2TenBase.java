package org.genvisis.widgets;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.ext;

public class Sci2TenBase {
	public Sci2TenBase(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + "-better.txt"));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				for (int i = 0; i < line.length; i++) {
					writer.print((i == 0 ? "" : "\t") + ext.prettyP(line[i]));
				}
				writer.println();
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

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "Sci2TenBase.txt";
		String filename = "makeMePretty.txt";

		String usage = "\n"	+ "park.Sci2TenBase requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default)\n" + "";

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
			new Sci2TenBase(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
