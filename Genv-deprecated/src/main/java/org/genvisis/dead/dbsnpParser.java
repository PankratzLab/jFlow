// used for MarkerAnalysis for 6K screen, not used for anything else
package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;

public class dbsnpParser {
	public dbsnpParser(String filename) throws IOException {
		Hashtable<String, String[]> hash = getHash(filename);
		String[] info, keys = HashVec.getKeys(hash);
		PrintWriter writer = new PrintWriter(new FileWriter(filename + "-select.xls"));

		for (String key : keys) {
			info = hash.get(key);
			writer.print(key);
			for (String element : info) {
				writer.print("\t" + element);
			}
			writer.println();
		}
		writer.close();

	}

	public static Hashtable<String, String[]> getHash(String filename) {
		BufferedReader reader = null;
		String[] line, info = null;
		String temp, trav = "nada";
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		int offset;
		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			reader.readLine();
			while (reader.ready()) {
				temp = reader.readLine();
				if (reader.ready()) {
					if (temp.equals("")) {
						hash.put(trav = reader.readLine().split("[\\s]+")[0],
										 info = ArrayUtils.stringArray(9, "-1"));
					}
					line = temp.split("\\|");
					if (temp.startsWith("SNP")) {
						if (temp.startsWith("SNP | alleles='")) {
							info[0] = temp.substring(15, 16);
							info[1] = temp.substring(17, 18);
						} else {
							System.err.println("Strange SNP line (no alleles): '" + temp + "'");
						}
						if (temp.indexOf("het=") > 0) {
							line[2] = line[2].trim().split("=")[1];
							info[2] = line[2].equals("?") ? "-1" : line[2];
						} else {
							System.err.println("Strange SNP line (no het): '" + temp + "'");
						}
					}
					if (temp.indexOf("assembly=reference") > 0 || temp.indexOf("assembly=Celera") > 0) {
						offset = (temp.indexOf("assembly=Celera") > 0 ? 1 : 0);
						if (line[2].startsWith(" chr=")) {
							line[2] = line[2].trim().split("=")[1];
							if (line[2].equals("X") || line[2].equals("Y")) {
								line[2] = "23";
							}
							if (line[2].equals("Un")) {
								line[2] = "?";
							}
							if (info[3 + offset].equals("-1") && !line[2].equals("?")) {
								info[3 + offset] = line[2];
							} else if (!line[2].equals("?") && !info[3 + offset].equals(line[2])) {
								System.err.println("Error - chromosomes don't agree between Celera/reference for marker "
																	 + trav + " (" + info[3 + offset] + " and " + line[2] + ")");
							}
						} else {
							System.err.println("Strange assembly=" + (offset == 1 ? "Celera" : "reference")
																 + " line (no chr): '" + temp + "'");
						}
						offset = (temp.indexOf("assembly=Celera") > 0 ? 2 : 0)
										 + (temp.indexOf("chr=Y") > 0 ? 1 : 0);
						if (line[3].startsWith(" chr-pos=")) {
							line[3] = line[3].trim().split("=")[1];
							if (info[5 + offset].equals("-1") && !line[3].equals("?")) {
								info[5 + offset] = line[3];
							} else if (!line[3].equals("?") && !line[3].equals(info[5 + offset])) {
								System.err.println("Error - more than one " + (offset == 1 ? "Celera" : "reference")
																	 + " chr-position for marker " + trav + " (" + info[5 + offset]
																	 + " and " + line[3] + ")");
							}
						} else {
							System.err.println("Strange assembly=" + (offset == 1 ? "Celera" : "reference")
																 + " line (no chr-pos): '" + temp + "'");
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.print("Marker position data source \"" + filename
											 + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return hash;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "070504153158_FLT.txt";

		String usage = "\n" + "park.dbsnpParser requires 0-1 arguments\n"
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
			new dbsnpParser(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
