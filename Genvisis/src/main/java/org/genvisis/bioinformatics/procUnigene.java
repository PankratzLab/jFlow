package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

public class procUnigene {
	public static final String[] VARS = {	"ID", "GENE ", "TITLE", "GENE_ID", "CHROMOSOME", "CYTOBAND",
																				"EXPRESS"};

	public procUnigene(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String temp = "", trav;
		Set<String> v = new LinkedHashSet<String>();
		Map<String, Integer> vCounts = new HashMap<String, Integer>();
		int count = 0;
		String[] data = ArrayUtils.stringArray(VARS.length + 1, ".");

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + ".xls"));
			writer.println("NM_ID\t" + ArrayUtils.toStr(VARS));
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.equals("//")) {
					writer.println(ArrayUtils.toStr(data));
					line = data[ext.indexOfStr("EXPRESS", VARS) + 1].split("\\|");
					if (!line[0].equals(".") && !data[3 + 1].equals(".")) {
						count++;
						for (String element : line) {
							trav = element.trim();
							if (v.contains(trav)) {
								vCounts.put(trav, vCounts.get(trav) + 1);
							} else {
								v.add(trav);
								vCounts.put(trav, 1);
							}
						}
					}
					data = ArrayUtils.stringArray(VARS.length + 1, ".");
				}
				for (int i = 0; i < VARS.length; i++) {
					if (temp.startsWith(VARS[i])) {
						data[i + 1] = temp.substring(VARS[i].length()).trim();
					}
				}
				if (temp.indexOf("ACC=NM_") > 0) {
					temp = temp.indexOf(";") > 0	? temp.substring(temp.indexOf("ACC=NM_")	+ 4,
																													temp.indexOf(";"))
																				: temp.substring(temp.indexOf("ACC=NM_") + 4)
																							.split("[\\s]+")[0];
					if (data[0].equals(".")) {
						data[0] = temp;
					} else {
						data[0] += ";" + temp;
					}
				}

			}
			writer.close();
			reader.close();

			writer = new PrintWriter(new FileWriter("possibleTissues.out"));
			List<String> sortedVals = new ArrayList<String>(v);
			Collections.sort(sortedVals);
			for (int i = 0; i < sortedVals.size(); i++) {
				String s = sortedVals.get(i);
				writer.println(s	+ "\t" + vCounts.get(s) + "\t"
												+ ((double) vCounts.get(s) / (double) count));
			}
			writer.println();
			writer.println("Total\t" + count);
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
		String filename = "C:\\Download\\Hs.data";

		String usage = "\n"	+ "park.procUnigene requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

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
			new procUnigene(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
