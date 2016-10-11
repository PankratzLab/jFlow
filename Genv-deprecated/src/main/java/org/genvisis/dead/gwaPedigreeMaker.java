package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.ext;
import org.genvisis.park.tools;

public class gwaPedigreeMaker {
	public static final String[] HEADER = {	"FamID", "IndID", "UniqueID", "Father", "Mother", "Family",
																					"Individ", "Father", "Mother", "Sex", "DNA", "IRB",
																					"Shipment_Site", "DNA_Source", "DNA_ExtMeth", "Local_ID"};

	public gwaPedigreeMaker(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, trav = null;
		String prev;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Vector<String[]> v = new Vector<String[]>();
		int count;
		boolean done = false;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + "-out.xls"));
			line = reader.readLine().split("\t", -1);
			writer.println(Array.toStr(line));
			ext.checkHeader(line, HEADER, true);
			prev = "";
			while (!done) {
				if (reader.ready()) {
					trav = reader.readLine().split("\t");
				} else {
					prev = "";
					done = true;
				}
				if (!trav[0].equals(prev)) {
					if (v.size() == 1) {
						writer.println(Array.toStr(v.elementAt(0)));
					} else {
						count = 0;
						hash.clear();
						for (int i = 0; i < v.size(); i++) {
							line = v.elementAt(i);
							hash.put(line[1], line[6]);
						}
						for (int i = 0; i < v.size(); i++) {
							line = v.elementAt(i);
							for (int j = 3; j <= 4; j++) {
								if (!hash.containsKey(line[j])) {
									count++;
									writer.println(line[0]	+ "\t" + line[j] + "\t"
																	+ tools.getUniqueID(line[0], line[j]) + "\t.\t.\t" + line[5]
																	+ "\t" + line[5] + "1" + count + "\t0\t0\t" + (j - 2) + "\t0");
									hash.put(line[j], line[5] + "1" + count);
								}
								line[4 + j] = hash.get(line[j]);
							}
						}
						for (int i = 0; i < v.size(); i++) {
							writer.println(Array.toStr(v.elementAt(i)));
						}
					}

					v.removeAllElements();
				}
				v.add(trav);
				prev = trav[0];
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
		// String filename = "GWA_Sample_Pedigree_starter.txt";
		String filename = "Boston_GWA_Sample-Pedigree_file-starter.txt";

		String usage = "\n"	+ "park.gwaPedigreeMaker requires 0-1 arguments\n"
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
			new gwaPedigreeMaker(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
