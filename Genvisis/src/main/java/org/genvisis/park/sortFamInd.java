// later we can update this to specify columns and number of header rows
package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class sortFamInd {

	public sortFamInd(int numLines_to_skip) throws IOException {
		for (int i = 1; i <= 23; i++) {
			new sortFamInd("chromosome" + i + ".dat", numLines_to_skip);
		}
	}

	public sortFamInd(String filename, int numLines_to_skip) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		StringTokenizer st;
		String temp, trav, dna;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Vector<String> vex = new Vector<String>();

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find " + filename + " in current directory");
			System.exit(2);
		}
		String bakFilename = Files.getBakFilename(filename, super.getClass().getName());
		(new File(filename)).renameTo(new File(bakFilename));
		reader = new BufferedReader(new FileReader(bakFilename));
		writer = new PrintWriter(new FileWriter(filename));

		for (int i = 0; i < numLines_to_skip; i++) {
			writer.println(reader.readLine());
		}

		while (reader.ready()) {
			temp = reader.readLine();
			st = new StringTokenizer(temp);
			dna = st.nextToken();
			trav = st.nextToken() + "-" + ext.formNum(st.nextToken(), 3) + "-" + dna;
			vex.add(trav);
			hash.put(trav, temp);
		}

		int[] key = Sort.quicksort(Array.toStringArray(vex));

		for (int element : key) {
			writer.println(hash.get(vex.elementAt(element)));
		}
		reader.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		int toSkip = 2, numArgs = args.length;
		for (int i = 0; i < args.length; i++) {
			if (args[i].startsWith("skip=")) {
				try {
					toSkip = Integer.valueOf(args[i].split("=")[1]).intValue();
					numArgs--;
				} catch (Exception e) {
					System.err.println("Error processing argument " + (i + 1) + " (" + args[i] + ").");
					System.err.println("Try 'skip=1' if you want to skip 1 line of the file.");
				}
			}
		}

		if (numArgs > 1) {
			System.out.println("Expecting 0-2 arguments: 1) filename (default: all 23 chromosome#.dat files).");
			System.out.println("                         2) skip=0 (number of lines to skip at beginning of file (default: 2)");
		} else {
			try {
				if (numArgs == 1) {
					new sortFamInd(args[0], toSkip);
				} else {
					new sortFamInd(toSkip);
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.err.println("Error in processing.");
			}
		}
	}
}
