package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.ext;

public class parametric {

	public parametric(int chromosome, boolean dominant, double freq, double penetrance,
										double phenocopy) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		String chrome;

		if (chromosome < 10) {
			chrome = "0" + chromosome;
		} else {
			chrome = "" + chromosome;
		}

		try {
			reader = new BufferedReader(new FileReader("map" + chrome + ".dat"));
		} catch (Exception e) {
			try {
				reader = new BufferedReader(new FileReader("/home/npankrat/park/00masters/map"	+ chrome
																										+ ".dat"));
				System.err.println("Could not find map" + chrome + ".dat in the current directory");
				System.err.println("  using the one in /home/npankrat/park/00masters/");
			} catch (Exception e2) {
				System.err.println("Could not find map"	+ chrome
														+ ".dat in /home/npankrat/park/00masters/ or in the current directory");
				System.exit(1);
			}

		}

		if (dominant) {
			writer = new PrintWriter(new FileWriter("map" + chrome + ".D.dat"));
		} else {
			writer = new PrintWriter(new FileWriter("map" + chrome + ".R.dat"));
		}

		for (int i = 0; i < 3; i++) {
			writer.println(reader.readLine());
		}
		for (int i = 0; i < 4; i++) {
			reader.readLine();
		}

		writer.println("1  2  << AFFECTATION, NO. OF ALLELES");
		writer.println(ext.formDeci(1 - freq, 8, true)	+ " " + ext.formDeci(freq, 8, true)
										+ "  << GENE FREQUENCIES");
		writer.println("1  << NO. OF LIABILITY CLASSES");
		writer.println(ext.formDeci(phenocopy, 4, true)	+ " "
										+ ext.formDeci((dominant ? penetrance : phenocopy), 4, true) + " "
										+ ext.formDeci(penetrance, 4, true));

		if (chromosome == 23) {
			reader.readLine();
			writer.println(ext.formDeci(phenocopy, 4, true) + " " + ext.formDeci(penetrance, 4, true));
		}

		while (reader.ready()) {
			writer.println(reader.readLine());
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		boolean dominant = true;
		double freq = 0.005, penetrance = 0.80, phenocopy = 0.03;

		String usage = "\n"	+ "park.bat.parametric requires 2-5 arguments:\n"
										+ "   (1) a chromosome number (i.e. 2)\n"
										+ "   (2) dominant/recessive (i.e. dominant)\n"
										+ "   (3) disease allele freq (i.e. freq=" + freq + " (default))\n"
										+ "   (4) penetrance (i.e. pen=" + penetrance + " (default))\n"
										+ "   (5) phenocopy rate (i.e. phenocopy=" + phenocopy + " (default))\n" + "";

		if (!args[1].equals("dominant")) {
			freq = 0.10;
		}

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h")	|| args[i].equals("-help") || args[i].equals("/h")
					|| args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("freq=")) {
				freq = Double.valueOf(args[i].split("=")[1]).doubleValue();
				numArgs--;
			} else if (args[i].startsWith("pen=")) {
				penetrance = Double.valueOf(args[i].split("=")[1]).doubleValue();
				numArgs--;
			} else if (args[i].startsWith("phenocopy=")) {
				phenocopy = Double.valueOf(args[i].split("=")[1]).doubleValue();
				numArgs--;
			} else if (args[i].startsWith("dominant")) {
				dominant = true;
			} else if (args[i].startsWith("recessive")) {
				dominant = false;
			} else {
				int chr = -1;
				if (i == 0) {
					try {
						chr = Integer.valueOf(args[0]).intValue();
					} catch (Exception e) {
						System.err.println("Error - invalid chromosome number: " + args[0]);
						System.exit(3);
					}
				}
				if (chr < 1 || chr > 23) {
					System.err.println("Error - don't know what to do with the argument: " + args[i]);
					System.err.println("usage");
					System.exit(2);
				}
			}
		}
		if (numArgs != 2) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new parametric(Integer.valueOf(args[0]).intValue(), dominant, freq, penetrance, phenocopy);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
