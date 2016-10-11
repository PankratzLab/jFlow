package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Cyrillic {
	public static void format(String prefile, String mapfile, String outprefile, String outmapfile,
														int chr) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int numMrkrs, numClasses;

		try {
			reader = new BufferedReader(new FileReader(prefile));
			writer = new PrintWriter(new FileWriter(outprefile));
			while (reader.ready()) {
				writer.println(reader.readLine() + "  << ");
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - could not find the pre file \"" + prefile + "\"");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error parsing the pre file \"" + prefile + "\"");
			System.exit(1);
		}

		try {
			reader = new BufferedReader(new FileReader(mapfile));
			writer = new PrintWriter(new FileWriter(outmapfile));

			line = reader.readLine().trim().split("[\\s]+");
			numMrkrs = Integer.valueOf(line[0]).intValue();
			writer.print(numMrkrs);
			for (int j = 1; j <= 3; j++) {
				writer.print("  " + line[j]);
			}
			writer.println("<< chromosome: \"" + chr + "\"");
			line = reader.readLine().trim().split("[\\s]+");
			for (int j = 0; j < 4; j++) {
				writer.print(line[j] + "  ");
			}
			writer.println();

			writer.println(reader.readLine());
			reader.readLine();
			writer.println("1  2  << Name: \"Disease Name\"");
			line = reader.readLine().trim().split("[\\s]+");
			for (int j = 0; j < 2; j++) {
				writer.print("  " + line[j]);
			}
			writer.println();

			line = reader.readLine().trim().split("[\\s]+");
			numClasses = Integer.parseInt(line[0]);
			writer.println(numClasses);
			for (int j = 0; j < numClasses; j++) {
				line = reader.readLine().trim().split("[\\s]+");
				for (String element : line) {
					writer.print("  " + element);
				}
				writer.println();
			}

			for (int j = 0; j < numMrkrs - 1; j++) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.print(" " + line[0] + "  " + line[1]);
				writer.println(" << Name: \"" + (line.length > 2 ? line[3] : "Marker_" + (j + 1)) + "\"");
				writer.println(reader.readLine());
			}
			line = reader.readLine().trim().split("[\\s]+");
			for (int j = 0; j < 2; j++) {
				writer.print(line[j] + " ");
			}
			writer.println();
			line = reader.readLine().trim().split("[\\s]+");
			for (int j = 0; j < numMrkrs - 1; j++) {
				writer.print(line[j] + " ");
			}
			writer.println();
			line = reader.readLine().trim().split("[\\s]+");
			for (int j = 0; j < 3; j++) {
				writer.print("  " + line[j]);
			}
			writer.println();

			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: could not find the map file \"" + mapfile + "\"");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error: could not parse the map file \"" + mapfile + "\"");
			System.exit(1);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\David
		// Ginsburg\\line3\\2008.08.07\\chr12markers\\";
		String dir = "";
		int chr = 12;
		String prefile = "";
		String mapfile = "";

		String usage = "\n"	+ "park.newClass requires 0-3 arguments\n" + "   (0) directory (i.e. dir="
										+ dir + " (default))\n" + "   (1) chromosome number (i.e. chr=" + chr
										+ " (default))\n"
										+ "   (2) name of .pre file (i.e. pre=re_chromXX.pre (default))\n"
										+ "   (3) name of MAP file (i.e. map=mapXX.dat (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("pre=")) {
				prefile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("map=")) {
				mapfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("chr=")) {
				chr = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		if (prefile.equals("")) {
			prefile = "re_chrom" + ((chr < 10) ? "0" + chr : "" + chr) + ".pre";
		}
		if (mapfile.equals("")) {
			mapfile = "map" + ((chr < 10) ? "0" + chr : "" + chr) + ".dat";
		}

		try {
			format(dir + prefile, dir + mapfile, dir + "cy_" + prefile, dir + "cy_" + mapfile, chr);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
