package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ext;

public class FamilyHistory {
	public static final String CRF_DIR = tools.CRF_DIR;
	public static final String DEFAULT_DB_FILE = "crf_db.dat";

	public static void proc(String trait) {
		proc(DEFAULT_DB_FILE, trait);
	}

	public static void proc(String db_file, String trait) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, data, trav;
		Hashtable<String, String[]> hashStringArray = new Hashtable<String, String[]>();
		Hashtable<String, Vector<String[]>> hashVecStringArray = new Hashtable<String, Vector<String[]>>();
		int index;
		Vector<String[]> v;
		Vector<String> variance = new Vector<String>();
		boolean quant = false;
		double[] avg_max;
		int[] counts;

		try {
			reader = tools.getNinfoReader(2, false);
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				hashStringArray.put(line[0] + "\t" + line[1], new String[] {line[4], line[5]});
			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error reading fino2 file");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(CRF_DIR + db_file));
			index = ext.indexFactors(	new String[] {trait}, reader.readLine().split("\t", -1), true,
																true)[0];
			while (reader.ready()) {
				line = reader.readLine().split("\t", -1);
				if (hashVecStringArray.containsKey(line[1])) {
					v = hashVecStringArray.get(line[1]);
				} else {
					hashVecStringArray.put(line[1], v = new Vector<String[]>());
				}
				data = hashStringArray.get(line[1] + "\t" + line[2]);
				v.add(new String[] {line[2], line[index], data[0], data[1]});
				if (variance.size() < 3 && !line[index].equals(".") && !variance.contains(line[index])) {
					variance.add(line[index]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + db_file + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + db_file + "\"");
			System.exit(2);
		}

		quant = variance.size() > 2;

		try {
			reader = new BufferedReader(new FileReader(CRF_DIR + db_file));
			writer = new PrintWriter(new FileWriter(CRF_DIR + trait + "_sibHistory.csv"));
			if (quant) {
				writer.println("FamID,IndID,AvgRelative,MaxRelative,AvgFirstDegree,MaxFirstDegree,AvgWithinSibship,NumberOfSiblingsWithNonMissingData");
			} else {
				writer.println("FamID,IndID,NumRelativesWithHistory,NumFirstDegreesWithHistory,FamHist,FirstDegreeHist");
			}
			index = ext.indexFactors(	new String[] {trait}, reader.readLine().split("\t", -1), true,
																true)[0];
			while (reader.ready()) {
				line = reader.readLine().split("\t", -1);
				v = hashVecStringArray.get(line[1]);
				counts = new int[3];
				avg_max = new double[5];
				trav = null;
				for (int i = 0; i < v.size(); i++) {
					trav = v.elementAt(i)[0].equals(line[2]) ? v.elementAt(i) : trav;
				}
				for (int i = 0; i < v.size(); i++) {
					data = v.elementAt(i);
					if (!data[1].equals(".")) {
						// relatives calculations
						if (!data[0].equals(line[2])) {
							if (quant) {
								counts[0]++;
								avg_max[0] += Double.parseDouble(data[1]);
								avg_max[1] = Double.parseDouble(data[1]) > avg_max[1]	? Double.parseDouble(data[1])
																																			: avg_max[1];
								// first degree relatives calculations
								if ((data[2].equals(trav[2]) && data[3].equals(trav[3]))	|| data[0].equals(trav[2])
										|| data[0].equals(trav[3]) || data[2].equals(trav[0])
										|| data[3].equals(trav[0])) {
									counts[1]++;
									avg_max[2] += Double.parseDouble(data[1]);
									avg_max[3] = Double.parseDouble(data[1]) > avg_max[3]	? Double.parseDouble(data[1])
																																				: avg_max[3];
								}
							} else if (data[1].equals("1")) {
								counts[0]++;
								// first degree relatives calculations
								if ((data[2].equals(trav[2]) && data[3].equals(trav[3]))	|| data[0].equals(trav[2])
										|| data[0].equals(trav[3]) || data[2].equals(trav[0])
										|| data[3].equals(trav[0])) {
									counts[1]++;
								}

							}
						}

						// within sibship calculations
						if (quant && data[2].equals(trav[2]) && data[3].equals(trav[3])) {
							avg_max[4] += Double.parseDouble(data[1]);
							counts[2]++;
						}

					}
				}
				if (avg_max[2] == 0 && !trav[1].equals(".")) {
					avg_max[2] = Double.parseDouble(trav[1]);
				}
				if (quant) {
					writer.println(line[1]	+ "," + line[2] + ","
													+ (counts[0] > 0 ? (avg_max[0] / counts[0]) + "," + avg_max[1] : ".,.")
													+ ","
													+ (counts[1] > 0 ? (avg_max[2] / counts[1]) + "," + avg_max[3] : ".,.")
													+ "," + (counts[2] > 0 ? (avg_max[4] / counts[2]) : ".") + ","
													+ counts[2]);
				} else {
					writer.println(line[1]	+ "," + line[2] + "," + counts[0] + "," + counts[1] + ","
													+ (counts[0] > 0 ? 1 : 0) + "," + (counts[1] > 0 ? 1 : 0));
				}
			}
			reader.close();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing family history of " + trait);
			e.printStackTrace();
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String trait = "Depression";
		String trait = "AOO";
		String dbfile = DEFAULT_DB_FILE;

		String usage = "\n"	+ "db.FamilyHistory requires 0-1 arguments\n" + "   (1) trait (i.e. trait="
										+ trait + " (default)\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("db=")) {
				dbfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("trait=")) {
				trait = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proc(dbfile, trait);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
