package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

public class showMe {
	public showMe(String filename, boolean separate, boolean trim,
								boolean pelican) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, last;
		String temp, prev, aff;
		Vector<String> v = new Vector<String>(), veeps = new Vector<String>();
		Vector<String[]> peeps = new Vector<String[]>();
		Hashtable<String, String> ages = new Hashtable<String, String>(), affection;

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find " + filename + " in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(filename));
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			v.add(line[0]);
		}
		reader.close();

		reader = tools.getNinfoReader(1, false);
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			if (v.contains(line[0])) {
				veeps.add(line[0] + "\t" + line[1]);
				line[2] = line[2].equals(".") ? "?" : line[2];
				ages.put(line[0] + "\t" + line[1], line[2]);
			}
		}
		reader.close();

		affection = tools.getBestPDdx();

		reader = tools.getNinfoReader(2, false);
		while (reader.ready()) {
			temp = reader.readLine();
			line = temp.split("[\\s]+");
			if (v.contains(line[0])) {
				last = peeps.isEmpty() ? new String[] {"", "", "", "", "", ""} : peeps.lastElement();
				if ((last[0] + "\t" + last[1]).equals(line[0] + "\t" + line[1])) {

				} else {
					aff = tools.isAffected(affection, line[0] + "\t" + line[1]) ? "2" : "1";
					if (pelican) {
						peeps.add(new String[] {line[0], line[1], (line[4].equals(".") ? "0" : line[4]),
																		(line[5].equals(".") ? "0" : line[5]),
																		(line[2].toUpperCase().equals("M")	? "1"
																																				: line[2]	.toUpperCase()
																																									.equals("F")	? "2"
																																																: "0"),
																		aff, (veeps.contains(line[0] + "\t" + line[1]) ? "1" : "0"),
																		(veeps.contains(line[0] + "\t" + line[1]) ? "1" : "0"),
																		(ages.containsKey(line[0] + "\t" + line[1])
																			&& aff.equals("2")	? ages.get(line[0] + "\t" + line[1])
																													: "-"),
																		(affection.containsKey(line[0]	+ "\t"
																														+ line[1])	? affection.get(line[0]
																																													+ "\t"
																																												+ line[1])
																																				: "-"),
																		"<< <PelicanData>"													+ (line[3].equals("TRUE")	? "1"
																																																					: "0")
																																								+ (line[1].equals("1")	? "1"
																																																				: "0")
																																								+ aff
																																								+ "</PelicanData>"});
					} else {
						peeps.add(new String[] {line[0], line[1], (line[4].equals(".") ? "0" : line[4]),
																		(line[5].equals(".") ? "0" : line[5]),
																		(line[2].toUpperCase().equals("M")	? "1"
																																				: line[2]	.toUpperCase()
																																									.equals("F")	? "2"
																																																: "0"),
																		aff, (veeps.contains(line[0] + "\t" + line[1]) ? "1" : "0"),
																		(veeps.contains(line[0] + "\t" + line[1]) ? "1" : "0")});

					}
				}
			}
		}
		reader.close();

		// if (trim) {
		// System.err.println("Error - consol.prune disabled, see if you can get TrimFam to work");
		// // consol.prune.procPeeps(peeps, consol.prune.findKeepers(peeps,
		// // veeps));
		// }

		if (separate) {
			while (!peeps.isEmpty()) {
				writer = new PrintWriter(new FileWriter(peeps.elementAt(0)[0] + ".pre"));
				do {
					line = peeps.remove(0);
					writer.println(ArrayUtils.toStr(line));
					prev = line[0];
				} while (!peeps.isEmpty() && (separate && peeps.elementAt(0)[0].equals(prev)));
				writer.close();
			}
		} else {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename) + ".pre"));
			for (int i = 0; i < peeps.size(); i++) {
				writer.println(ArrayUtils.toStr(peeps.elementAt(i)));
			}
			writer.close();
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "fams.dat";
		// String filename = "SFHfams.dat";
		boolean separate = false;
		boolean trimmed = true;
		boolean pelican = true;

		String usage = "\n"	+ "park.showMe requires 0-1 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n"
										+ "   (2) put families in separate files (i.e. '-separate' (optional))\n"
										+ "   (3) trim families down to those genotyped (i.e. '-all' (optional))\n"
										+ "   (4) add in pelican info (i.e. '-pelican' (optional))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.equals("separate=")) {
				separate = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			} else if (arg.equals("trim=")) {
				trimmed = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			} else if (arg.equals("pelican=")) {
				pelican = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new showMe(filename, separate, trimmed, pelican);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
