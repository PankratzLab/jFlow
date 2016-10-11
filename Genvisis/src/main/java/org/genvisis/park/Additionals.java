package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Vectors;
import org.genvisis.common.ext;

public class Additionals {
	public static void affRent() {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String fam, id = "", trav, prev;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		boolean affy;

		try {
			reader = tools.getNinfoReader(2, false);
			line = reader.readLine().split("[\\s]+");
			ext.checkHeader(line, tools.NINFO2_HEADER, true);

			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (line.length < 8) {
					System.err.println("Error: make sure every entry has 8 columns (zero out founders' parents, '.' out missing values)");
					System.exit(1);
				}
				id = line[0] + "\t" + line[1];
				if (line[6].equals("PD")) {
					hash.put(id, "1");
				}
			}
			reader.close();
		} catch (Exception e) {
			System.err.println("Error reading ninfo1 file. Check for missing data, and replace with periods.");
		}

		try {
			reader = tools.getNinfoReader(2, false);
			writer = new PrintWriter(new FileWriter(tools.CRF_DIR + "ninfo2_affRent.xln"));
			writer.println("FamID\tIndID\tUniqueID\tAffFather\tAffMother\tAffParent");

			prev = "";
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				fam = line[0];
				id = line[1];
				trav = fam + "\t" + id;

				if (!trav.equals(prev)) {
					writer.print(trav + "\t" + (Integer.valueOf(fam).intValue() * 1000
																			+ Integer.valueOf(id).intValue()));
					affy = false;
					if (hash.containsKey(fam + "\t" + line[4])) {
						writer.print("\t1");
						affy = true;
					} else {
						writer.print("\t0");
					}
					if (hash.containsKey(fam + "\t" + line[5])) {
						writer.print("\t1");
						affy = true;
					} else {
						writer.print("\t0");
					}
					if (affy) {
						writer.println("\t1");
					} else {
						writer.println("\t0");
					}
				}
				prev = trav;
			}

			reader.close();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error reading ninfo2 file. Check for missing data, and replace with periods.");
		}

		try {
			hash = new Hashtable<String, String>();
			reader = Files.getReader("ninfo1.csv", tools.CRF_DIR); // (csv file
			// has 1
			// indiviual
			// per line)
			while (reader.ready()) {
				line = reader.readLine().split(",");
				hash.put(line[0] + ":" + line[1], "");
			}
			reader.close();

			reader = new BufferedReader(new FileReader(tools.CRF_DIR + "ninfo2_affRent.xln"));
			writer = new PrintWriter(new FileWriter(tools.CRF_DIR + "affRent.csv"));
			line = reader.readLine().split("[\\s]+");
			writer.println(line[0] + "," + line[1] + "," + line[3] + "," + line[4] + "," + line[5]);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (hash.containsKey(line[0] + ":" + line[1])) {
					writer.println(line[0] + "," + line[1] + "," + line[3] + "," + line[4] + "," + line[5]);
				}
			}
			reader.close();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error reading ninfo1.csv (1 indiviual per line) file. Check for missing data, and replace with periods.");
		}

	}

	public static void anyExtras() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav, fam;
		Vector<String> typed = new Vector<String>();
		Vector<String> fams = new Vector<String>();
		Hashtable<String, String> affecteds = new Hashtable<String, String>();
		Hashtable<String, Vector<Vector<String>>> hash = new Hashtable<String, Vector<Vector<String>>>();
		Vector<Vector<String>> data;
		int numDads, numMoms;

		try {
			reader = tools.getNinfoReader(1, false);
			line = reader.readLine().split("[\\s]+");
			for (int i = 0; i < line.length; i++) {
				if (!line[i].equals(tools.NINFO1_HEADER[i])) {
					System.err.println("Error - header has been changed; expecting: "	+ tools.NINFO1_HEADER[i]
															+ ", got: " + line[i]);
					System.exit(1);
				}
			}

			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (!typed.contains(line[0] + "\t" + line[1])) {
					typed.add(line[0] + "\t" + line[1]);
				}
				if (!fams.contains(line[0])) {
					fams.add(line[0]);
				}
			}
			reader.close();
		} catch (Exception e) {
			System.err.println("Error processing ninfo1 file.");
			e.printStackTrace();

		}

		affecteds = tools.getBestPDdx();

		try {
			reader = tools.getNinfoReader(2, false);
			reader.readLine();

			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trav = line[0] + "\t" + line[1];
				if (hash.containsKey(line[0])) {
					data = hash.get(line[0]);
				} else {
					hash.put(line[0], data = HashVec.newVecVecString(7));
				}

				if (tools.isAffected(affecteds, trav)) {
					if (typed.contains(trav)) {
						HashVec.addIfAbsent(line[1], data.elementAt(0));
						if (affecteds.get(trav).equals("VPD")) {
							HashVec.addIfAbsent(line[1], data.elementAt(6));
						}

						if (tools.isAffected(affecteds, line[0] + "\t" + line[4])) { // father
																																					// affected
							data.elementAt(4).add(line[4]);
						}
						if (tools.isAffected(affecteds, line[0] + "\t" + line[5])) { // mother
																																					// affected
							data.elementAt(5).add(line[5]);
						}
					} else {
						if (line[3].equals("TRUE")) { // are you dead?
							HashVec.addIfAbsent(line[1], data.elementAt(3));
						} else {
							HashVec.addIfAbsent(line[1], data.elementAt(2));
						}
					}
				} else {
					if (typed.contains(trav)) {
						HashVec.addIfAbsent(line[1], data.elementAt(1));
					}
				}
			}
			reader.close();
		} catch (Exception e) {
			System.err.println("Error parsing ninfo2 file.");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(tools.CRF_DIR + "anyExtras.xln"));
			writer.println("FamID\ttotalAffs\thaveAffs\tVPDs\thaveUnaffs\tabsAffs\tdeadAffs\ttotalGenotyped\taffDads\taffMoms\tgetThese\tVPD_IDs");
			for (int i = 0; i < fams.size(); i++) {
				fam = fams.elementAt(i);
				data = hash.get(fam);

				numDads = numMoms = 0;
				while (data.elementAt(4).size() > 1) {
					trav = data.elementAt(4).remove(0);
					if (data.elementAt(4).contains(trav)) {
						numDads++;
						while (data.elementAt(4).contains(trav)) {
							data.elementAt(4).remove(trav);
						}
					}
				}
				while (data.elementAt(5).size() > 1) {
					trav = data.elementAt(5).remove(0);
					if (data.elementAt(5).contains(trav)) {
						numMoms++;
						while (data.elementAt(5).contains(trav)) {
							data.elementAt(5).remove(trav);
						}
					}
				}

				writer.println(fam	+ "\t"
												+ (data.elementAt(0).size()	+ data.elementAt(2).size()
														+ data.elementAt(3).size())
												+ "\t" + (data.elementAt(0).size()) + "\t" + (data.elementAt(6).size())
												+ "\t" + (data.elementAt(1).size()) + "\t" + (data.elementAt(2).size())
												+ "\t" + (data.elementAt(3).size()) + "\t"
												+ (data.elementAt(0).size() + data.elementAt(1).size()) + "\t" + (numDads)
												+ "\t" + (numMoms) + "\t"
												+ (Array.toStr(Array.toStringArray(data.elementAt(2)), " ")) + "\t"
												+ (Array.toStr(Array.toStringArray(data.elementAt(6)), " ")));

			}

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results.");
			e.printStackTrace();
		}
	}

	public static void countFirstDegreeRelatives() {
		BufferedReader reader;
		PrintWriter writer, log;
		String[] line;
		String trav;
		Hashtable<String, String> affectionStatus;
		Hashtable<String, IntVector[]> hash = new Hashtable<String, IntVector[]>();
		Hashtable<String, IntVector[]> siblings = new Hashtable<String, IntVector[]>();
		Hashtable<String, String> processed = new Hashtable<String, String>();
		IntVector[] counts, sibs, rents; // parentAff, parentVPD, siblingAff,
		// siblingVPD, childAff, childVPD
		String affStat, parentalAffStat, sibsRents;

		affectionStatus = tools.getBestPDdx();

		try {
			reader = tools.getNinfoReader(2, true);
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trav = line[0] + "\t" + line[1];
				affStat = affectionStatus.get(trav);
				if (!processed.containsKey(trav)) {
					counts = hash.containsKey(trav)	? hash.get(trav)
																					: Vectors.initializedArray(IntVector.class, 6);
					for (int i = 0; i < 2; i++) {
						parentalAffStat = affectionStatus.get(line[0] + "\t" + line[4 + i]);
						if (parentalAffStat != null) {
							if (ext.indexOfStr(parentalAffStat, tools.AFFECTED) >= 0) {
								counts[0].add(Integer.parseInt(line[4 + i]));
							}
							if (parentalAffStat.equals("VPD")) {
								counts[1].add(Integer.parseInt(line[4 + i]));
							}
							rents =
										hash.containsKey(line[0] + "\t" + line[4 + i])
																																			? hash.get(line[0]	+ "\t"
																																								+ line[4 + i])
																																		: Vectors.initializedArray(	IntVector.class,
																																																6);

							if (affStat != null && ext.indexOfStr(affStat, tools.AFFECTED) >= 0) {
								rents[4].add(Integer.parseInt(line[1]));
							}
							if (affStat != null && affStat.equals("VPD")) {
								rents[5].add(Integer.parseInt(line[1]));
							}
							hash.put(line[0] + "\t" + line[4 + i], rents);
						}
					}

					if (affStat != null && !line[4].equals("0") && !line[5].equals("0")) {
						sibsRents = line[0] + ":" + line[4] + ":" + line[5];
						sibs = siblings.containsKey(sibsRents)	? siblings.get(sibsRents)
																										: Vectors.initializedArray(IntVector.class, 2);

						if (ext.indexOfStr(affStat, tools.AFFECTED) >= 0) {
							sibs[0].add(Integer.parseInt(line[1]));
						}
						if (affStat.equals("VPD")) {
							sibs[1].add(Integer.parseInt(line[1]));
						}

						siblings.put(sibsRents, sibs);
						counts[2] = sibs[0];
						counts[3] = sibs[1];
					}

					hash.put(trav, counts);
					processed.put(trav, "");
				}

			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error reading ninfo2 file.");
			System.exit(2);
		}

		try {
			reader = tools.getNinfoReader(1, true);
			writer = new PrintWriter(new FileWriter(tools.CRF_DIR + "countFirstDegreeRelatives.xln"));
			log = new PrintWriter(new FileWriter(tools.CRF_DIR + "countFirstDegreeRelatives.log"));
			writer.println("FamID\tIndID\tUniqueID\tNumAddFirstDegreeAff\tNumAddFirstDegreeVPD\tNumParentAff\tNumParentVPD\tNumSiblingAff\tNumSiblingVPD\tNumChildAff\tNumChildVPD");
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trav = line[0] + "\t" + line[1];
				if (hash.containsKey(trav)) {
					counts = hash.get(trav);
					writer.println(line[0]	+ "\t" + line[1] + "\t" + tools.getUniqueID(line[0], line[1])
													+ "\t" + (counts[0].size() + counts[2].size() + counts[4].size()) + "\t"
													+ (counts[1].size() + counts[3].size() + counts[5].size()) + "\t"
													+ counts[0].size() + "\t" + counts[1].size() + "\t" + counts[2].size()
													+ "\t" + counts[3].size() + "\t" + counts[4].size() + "\t"
													+ counts[5].size());
				} else {
					log.println("Error - Why is there no data on " + trav + "?");

				}
			}
			reader.close();
			writer.close();
			log.close();
		} catch (IOException ioe) {
			System.err.println("Error reading ninfo1 file.");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Additionals.dat";

		String usage = "\n"	+ "park.Additionals requires 0-1 arguments\n"
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
			// affRent();
			// anyExtras();
			// countFirstDegreeRelatives();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
