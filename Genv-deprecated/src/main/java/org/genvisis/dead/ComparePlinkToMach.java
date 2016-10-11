package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class ComparePlinkToMach {

	public static void comp(String dir, int chr, String markersToDo) {
		BufferedReader reader1, reader2;
		// PrintWriter writer;
		String[] line;
		// String temp, trav;
		// Hashtable<String, String> hash = new Hashtable<String, String>();
		// Vector<String> v = new Vector<String>();
		// int count;

		// Vector<String> machPeopleIndicies;
		Vector<String> plinkPeopleIndicies;
		Hashtable<String, String[]> plinkData = new Hashtable<String, String[]>();
		String[] markers;

		markers = Array.toStringArray(HashVec.loadFileToVec(markersToDo, false, true, true));

		plinkPeopleIndicies = new Vector<String>();
		try {
			reader1 = new BufferedReader(new FileReader(dir + "chr" + chr + ".fam"));
			while (reader1.ready()) {
				line = reader1.readLine().trim().split("[\\s]+");
				if (Integer.parseInt(line[5]) > 0) {
					plinkPeopleIndicies.add(line[1]);
				}
			}
			reader1.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ dir + "chr" + chr + ".fam"
													+ "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + "chr" + chr + ".fam" + "\"");
			System.exit(2);
		}

		try {
			reader1 = new BufferedReader(new FileReader(dir + "chr" + chr + ".proxy.impute.dosage"));
			reader2 = new BufferedReader(new FileReader(dir + "chr" + chr + ".trans.mldose"));
			line = reader2.readLine().trim().split("[\\s]+");
			// machPeopleIndicies = Array.newVector(Array.subArray(line, 1));
			for (String marker : markers) {
				while (reader1.ready()) {
					line = reader1.readLine().trim().split("[\\s]+");
					if (ext.containsAny(line[0], markers)) {
						plinkData.put(line[0], line);
					}
				}
			}
			reader1.close();
			reader2.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ dir + "chr" + chr + ".proxy.impute.dosage"
													+ "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""	+ dir + "chr" + chr + ".proxy.impute.dosage"
													+ "\"");
			System.exit(2);
		}

	}

	public static void prepMach(int chr) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\imputation\\MACH comparison\\Mach_chr"
									+ chr + "\\";

		Vector<String> machMarkerIndices;

		try {
			if (new File(dir + "mach_step2_chr" + chr + ".mldose-transposeHuge.xln").exists()) {
				System.out.println("transposed file already exists for chr" + chr);
			} else {
				Files.transposeHuge(dir + "mach_step2_chr" + chr + ".mldose", 1000, new Logger());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		machMarkerIndices = HashVec.loadFileToVec(dir	+ "mach_step2_chr" + chr + ".mlinfo", true, true,
																							false);

		try {
			reader = new BufferedReader(new FileReader(dir	+ "mach_step2_chr" + chr
																									+ ".mldose-transposeHuge.xln"));
			writer = new PrintWriter(new FileWriter(dir + "mach_step2_chr" + chr + ".trans.mldose"));
			line = reader.readLine().trim().split("[\\s]+");
			writer.print("Marker");
			for (String element : line) {
				writer.print("\t" + element.substring(element.indexOf(">") + 1));
			}
			writer.println();
			reader.readLine();
			count = 0;
			while (reader.ready()) {
				writer.println(machMarkerIndices.elementAt(count) + "\t" + reader.readLine());
				count++;
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ dir + "mach_step2_chr" + chr
													+ ".mldose-transposeHuge.xln" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""	+ dir + "mach_step2_chr" + chr
													+ ".mldose-transposeHuge.xln" + "\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\imputation\\MACH comparison\\Mach_chr21\\";
		int chr = 21;
		String markersToDo = "markers.txt";

		String usage = "\\n"	+ "park.gwa.ComparePlinkToMach requires 0-1 arguments\n"
										+ "   (1) directory (i.e. dir=" + dir + " (default))\n"
										+ "   (2) chromosome (i.e. chr=" + chr + " (default))\n"
										+ "   (3) file with list of markers to do (i.e. file=" + chr + " (default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
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
		try {
			// prepMach(21);
			comp(dir, chr, markersToDo);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
