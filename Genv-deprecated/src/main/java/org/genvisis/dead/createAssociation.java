package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

public class createAssociation {
	public createAssociation(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		String[] line, markers;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		String blank;

		try {
			reader = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException ex) {
			throw new RuntimeException(ex.getMessage());
		}

		reader.readLine();
		line = reader.readLine().split("[\\s]+");
		markers = new String[line.length - 3];
		blank = "";
		for (int i = 0; i < markers.length; i++) {
			markers[i] = line[3 + i];
			blank += "\t0\t0";
		}
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			hash.put(line[1] + "\t" + line[2], line);
		}
		reader.close();

		try {
			reader = new BufferedReader(new FileReader("struct.dat"));
			writer = new PrintWriter(new FileWriter(filename.substring(0, filename.indexOf("."))
																							+ ".pre"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				for (int i = 0; i < 6; i++) {
					writer.print((i == 0 ? "" : "\t") + line[i]);
				}
				if (hash.containsKey(line[0] + "\t" + line[1])) {
					line = hash.get(line[0] + "\t" + line[1]);
					for (int i = 3; i < line.length; i++) {
						writer.print("\t" + line[i]);
					}
					writer.println();
				} else {
					writer.println(blank);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - could not find " + "struct.dat" + " in current directory");
			System.exit(2);
		} catch (IOException ioe) {
			System.err.println("Error parsing " + "struct.dat" + "");
			System.exit(3);
		}

		writer = new PrintWriter(new FileWriter(filename.substring(0, filename.indexOf("."))
																						+ "-map.dat"));
		writer.println((markers.length + 1)
										+ " 0 0 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM");
		writer.println("0 0.0 0.0 0  << MUT LOCUS, MUT RATE, HAPLOTYPE FREQUENCIES (IF 1)");
		for (int i = 1; i <= markers.length + 1; i++) {
			writer.print((i == 1 ? "" : " ") + i);
		}
		writer.println();
		writer.println("1  2  << AFFECTATION, NO. OF ALLELES");
		writer.println("0.950000 0.0500000  << GENE FREQUENCIES");
		writer.println("1  << NO. OF LIABILITY CLASSES");
		writer.println("0.020000 0.750000 0.750000");
		for (String marker : markers) {
			writer.println("3 2  # " + marker);
			writer.println("0.50000 0.50000");
		}
		writer.println("0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)");
		writer.print("10.0");
		for (int i = 1; i < markers.length; i++) {
			writer.print(" 3.0");
		}
		writer.println("  << RECOMB VALUES");
		writer.println("1 0.1 0.45  << REC VARIED, INCREMENT, FINISHING VALUE");

		writer.close();

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "chromosome2.dat";

		String usage = "\n"	+ "createLinkage requires 0-3 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n" + "";

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
			new createAssociation(filename);
		} catch (RuntimeException re) {
			System.err.println("Error: " + re.getMessage());
			if (re.getMessage().length() < 10) {
				re.printStackTrace();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
