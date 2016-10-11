// should make an option of creating structs for linkage versus association
// redo this whole thing? ugly code, only looks for consecutive family members
// what was with the only taking affecteds, has that screwed up every analysis performed to date?
package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ext;
import org.genvisis.link.TrimFam;

public class phenoStruct {
	public phenoStruct(boolean include_mutations) throws IOException {
		BufferedReader genos = null;
		BufferedReader info = null;
		BufferedReader fams = null;
		BufferedReader extra = null;
		PrintWriter unused = null, writer[] = new PrintWriter[4], finalDx;
		String[] line, subline;
		String temp, id = "", trav, prev;
		boolean done = false;
		Vector<String> members = new Vector<String>(), v;
		Vector<String> VPDmems = new Vector<String>();
		Hashtable<String, String> diagnosis, sex, aoo, ethnicity, dna; // add
		// covariate
		// like
		// this
		Hashtable<String, String> confirmedDx = new Hashtable<String, String>();
		Hashtable<String, Vector<String>> VIPlists = new Hashtable<String, Vector<String>>();
		Hashtable<String, Vector<String>> VPDlists = new Hashtable<String, Vector<String>>();
		TrimFam tf;
		int count = 0;

		Vector<String> affected = new Vector<String>();
		Vector<String> unaffected = new Vector<String>();
		Vector<String> halfsibs = new Vector<String>();
		Vector<String> deleteThese = new Vector<String>();
		Vector<String> halfling = new Vector<String>();
		Vector<String> vips;

		affected.add("VPD");
		affected.add("NVPD");
		affected.add("CONF");
		affected.add("CONF_PD");
		affected.add("RPD");

		unaffected.add("NRPD");
		unaffected.add("NXAM");
		unaffected.add("NOEV");
		unaffected.add(".");

		while (count <= 23 && !new File("chromosome" + (++count) + ".dat").exists()) {
			;
		}
		if (count > 23) {
			count = 0;
			while (count <= 23 && !new File("/home/npankrat/park/00masters/chromosome"	+ (++count)
																			+ ".dat").exists()) {
				;
			}
			if (count > 23) {
				System.err.println("Could not find a chromosome#.dat file in this directory or /home/npankrat/park/00masters/");
				System.err.println("Please rectify");
			} else {
				genos = new BufferedReader(new FileReader("/home/npankrat/park/00masters/chromosome"	+ count
																									+ ".dat"));
			}
		} else {
			genos = new BufferedReader(new FileReader("chromosome" + count + ".dat"));
		}

		info = tools.getNinfoReader(1, false);
		fams = tools.getNinfoReader(2, false);
		extra = tools.getNinfoReader(3, false);

		if (genos == null || info == null || fams == null) {
			System.exit(0);
		}

		while (extra != null && extra.ready()) {
			line = extra.readLine().split("[\\s]+");
			if (line.length == 1 && line[0].equals("")) {
				System.err.println("Warning - extra whitespace present in ninfo3 file");
				continue;
			}

			if (line[0].equals("CONF")) {
				confirmedDx.put(line[1] + "\t" + line[2], line[3]);
			} else if (line[0].equals("HS")) {
				halfsibs.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("MZ")) {
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("MZLINKAGE")) { // recently added, never
				// confirmed
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("DEL")) {
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("DELDNA")) { // recently added, never
				// confirmed
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("DELIND")) { // recently added, never
				// confirmed
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("DEL4LINKAGE")) { // recently added, never
				// confirmed
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("POORYIELD")) { // recently added, never
				// confirmed
				deleteThese.add(line[1] + "\t" + line[2]);
			} else if (line[0].equals("SWAPPED")) {
			} else if (line[0].equals("MUT")) {
				if (!include_mutations) {
					deleteThese.add(line[1] + "\t" + line[2]);
				}
			} else {
				System.err.println("Error parsing ninfo3 file: '" + line[0] + "'");
				System.err.println("    Proper syntax includes:");
				System.err.println("        CONF\t70001\t1\tVPD (DLB)");
				System.err.println("               or");
				System.err.println("        MZ\t70002\t1\t70002\t2");
				System.err.println("               or");
				System.err.println("        HS\t70003\t1\t\t(the half sib that is odd man out)");
				System.err.println("               or");
				System.err.println("        DEL\t70004\t1\t\t(to delete a record)");
			}
		}

		// collect all phenotypic information from the ninfo1 file
		diagnosis = new Hashtable<String, String>();
		sex = new Hashtable<String, String>(); // add covariate like this
		aoo = new Hashtable<String, String>();
		ethnicity = new Hashtable<String, String>();
		dna = new Hashtable<String, String>();
		try {
			ext.checkHeader(info.readLine().split("[\\s]+"), tools.NINFO1_HEADER, true);

			while (info.ready()) {
				line = info.readLine().split("[\\s]+");
				id = line[0] + "\t" + line[1];
				if (!aoo.containsKey(id)) {
					aoo.put(id, line[2]);
					temp = line[3];
					if (confirmedDx.containsKey(id)) {
						diagnosis.put(id, confirmedDx.get(id));
					} else {
						if (temp.equals("CONF")) {
							System.err.println("Error - Diagnosis for "	+ id
																	+ " was CONF, but autopsy results were not found in ninfo3 file");
							System.err.println("        if came back VPD, add line to ninfo3 file saying 'CONF\t"
																	+ id + "\tVPD (DLB)' or as appropriate");
							System.err.println("        because for now it's assumed to be NVPD");
						}
						diagnosis.put(id, temp);
					}
					sex.put(id, line[5].equals("M") ? "1" : (line[5].equals("F") ? "2" : line[5]));

					ethnicity.put(id, line[7]);
				} else {
					if (!line[2].equals(".") && (aoo.get(id)).equals(".")) {
						System.err.println(id + "  age of onset missing");
						aoo.put(id, line[2]);
					}
					if (!(diagnosis.get(id)).equals(line[3])) {
						if (confirmedDx.containsKey(id)) {
							diagnosis.put(id, confirmedDx.get(id));
						} else {
							System.err.println("Warning - "	+ id + " has two different diagnoses ("
																	+ (diagnosis.get(id)) + " and " + line[3]
																	+ ") and is not listed as an autopsy in ninfo3 file");
							System.err.println("        - add line to ninfo3 file saying 'CONF\t"	+ id
																	+ "\tVPD (DLB)' or as appropriate");
							System.exit(1);
						}
					}
					if (!line[7].equals(".") && (ethnicity.get(id)).equals(".")) {
						System.err.println(id + "  ethnicity");
						ethnicity.put(id, line[7]);
					}
				}
			}
		} catch (Exception e) {
			System.err.println("Error reading ninfo1 file. Check for missing data, and replace with periods.");
			e.printStackTrace();
		}
		info.close();

		// collect all genotyped individuals from the marker file
		prev = "";

		finalDx = new PrintWriter(new FileWriter("finalDiagnoses.dat"));
		id = "null";
		genos.readLine(); // placeholder
		genos.readLine(); // marker names
		while (!done) {
			if (genos.ready()) {
				line = genos.readLine().split("[\\s]+");
				trav = line[1];
				id = line[2];
				dna.put(trav + "\t" + line[2], line[0]);
			} else {
				trav = "";
				line = null;
				done = true;
			}

			if (!trav.equals(prev)) {
				VIPlists.put(prev, members);
				VPDlists.put(prev, VPDmems);
				prev = trav;
				members = new Vector<String>();
				VPDmems = new Vector<String>();
			}

			if (!done) {
				if (!diagnosis.containsKey(trav + "\t" + id)) {
					System.err.println("Error - diagnosis wasn't present for individual " + trav + "\t" + id);
					System.err.println("      - could be that chromosome"	+ count
															+ ".dat is empty or corrupt?");
					System.exit(21);
				}
				// if
				// (!unaffected.contains((diagnosis.get(trav+"\t"+id))))
				// {
				members.add(id);
				if ((diagnosis.get(trav + "\t" + id)).equals("VPD")) {
					VPDmems.add(id);
				}
				finalDx.println(trav + "\t" + id + "\t" + (diagnosis.get(trav + "\t" + id)));
				// }
			}
		}
		genos.close();
		finalDx.close();

		// collect family relationship information from the ninfo2 file, and
		// create phenostruct
		writer[0] = new PrintWriter(new FileWriter("struct111+.dat"));
		writer[1] = new PrintWriter(new FileWriter("struct111-.dat"));
		writer[2] = new PrintWriter(new FileWriter("struct100+.dat"));
		writer[3] = new PrintWriter(new FileWriter("struct100-.dat"));
		// writer = new PrintWriter(new FileWriter("phenostruct.dat"));

		line = fams.readLine().split("[\\s]+");
		ext.checkHeader(line, tools.NINFO2_HEADER, true);

		done = false;
		prev = "";
		members = new Vector<String>();
		while (!done) {
			if (fams.ready()) {
				line = fams.readLine().split("[\\s]+");
				trav = line[0];
				if (prev.equals("")) {
					prev = trav;
				}
			} else {
				trav = "";
				done = true;
			}

			if (done || !trav.equals(prev)) {
				if (VIPlists.containsKey(prev)) {
					for (int allvpd = 0; allvpd <= 1; allvpd++) {
						vips = (allvpd == 0 ? VIPlists.get(prev) : VPDlists.get(prev));
						if (allvpd == 0 && vips.size() < 2) {
							System.err.println("Warning - dropping family "	+ prev
																	+ " because it is uninformative for linkage");
						}
						for (int i = 0; i < deleteThese.size(); i++) {
							if ((deleteThese.elementAt(i)).startsWith(prev)) {
								vips.remove((deleteThese.elementAt(i)).split("[\\s]+")[1]);
							}
						}
						for (int i = 0; i < halfling.size(); i++) {
							for (int j = 0; j < members.size(); j++) {
								if ((members.elementAt(j)).startsWith(halfling.elementAt(i) + "\t")) {
									subline = (members.elementAt(j)).split("[\\s]+");
									members.removeElementAt(j);
									members.add(j, prev + "\t" + (9990 + i) + "\t0\t0\t1");
									temp = subline[0] + "\t" + subline[1] + "\t" + (9990 + i);
									temp += "\t" + subline[3] + "\t" + subline[4];
									members.add(j, temp);
								}
							}
						}

						System.err.println("ERROR! YOU HAVE YET TO CONVERT TO THE NEW INPUT FOR TrimFam!");
						// tf = new park.TrimFam(members, vips);
						tf = new TrimFam(members);

						for (int extnuc = 0; extnuc <= 1; extnuc++) {
							v = (extnuc == 0	? tf.getExtendedFamilyInformation()
																: tf.getNuclearFamilyInformation());

							if (checkForParentOffspringSinglets(v, vips)) {
								v.removeAllElements();
							}

							for (int i = 0; i < v.size(); i++) {
								subline = (v.elementAt(i)).split("[\\s]+");
								id = subline[1];
								writer[2 * allvpd + extnuc].print(v.elementAt(i));
								if (vips.contains(id + "")) {
									id = prev + "\t" + id;
									if (diagnosis.get(id) != null) {
										temp = diagnosis.get(id);
										if (affected.contains(temp)) {
											writer[2 * allvpd + extnuc].print("\t2");
										} else if (unaffected.contains(temp)) {
											writer[2 * allvpd + extnuc].print("\t0");
										} else {
											System.err.println("Warning - Unknown phenotype for individual '"	+ id + "': "
																					+ temp);
											writer[2 * allvpd + extnuc].print("\t0");
										}
										writer[2 * allvpd + extnuc].print("\t" + dna.get(id));
										writer[2 * allvpd + extnuc].print("\t" + aoo.get(id));
										writer[2 * allvpd + extnuc].print("\t" + ethnicity.get(id));
									} else {
										writer[2 * allvpd + extnuc].print("\t0\t" + dna.get(id) + "\t.\t.");
									}
								} else {
									writer[2 * allvpd + extnuc].print("\t0\tNoDNAprsnt\t.\t.");
								}
								writer[2 * allvpd + extnuc].println();
							}
						}
						if (tf.hasUnused()) {
							if (unused == null) {
								unused = new PrintWriter(new FileWriter("unused_individuals-"
																												+ (allvpd == 0 ? "111" : "100") + ".dat"));
							}
							v = tf.getUnused();
							for (int i = 0; i < v.size(); i++) {
								unused.println(prev + "\t" + v.elementAt(i));
							}
						}
					}
				}

				prev = trav;
				members.removeAllElements();
				halfling.removeAllElements();
			}

			id = trav + "\t" + line[1];
			temp = id	+ "\t" + line[4] + "\t" + line[5] + "\t"
							+ (line[2].equals("M") ? "1" : (line[2].equals("F") ? "2" : "0"));
			members.add(temp);
			if (halfsibs.contains(id)) {
				halfling.add(id);
			}
		}
		fams.close();

		writer[0].close();
		writer[1].close();
		writer[2].close();
		writer[3].close();
		if (unused != null) {
			unused.close();
		}
	}

	public boolean checkForParentOffspringSinglets(Vector<String> v, Vector<String> vips) {
		Vector<String> presentVIPs = new Vector<String>();

		for (int i = 0; i < v.size(); i++) {
			if (vips.contains(v.elementAt(i).split("[\\s]+")[1])) {
				presentVIPs.add(v.elementAt(i));
			}
		}

		if (presentVIPs.size() == 2) {
			if ((v.elementAt(0)).split("[\\s]+")[1].equals((v.elementAt(1)).split("[\\s]+")[2])
						|| (v.elementAt(0)).split("[\\s]+")[1].equals((v.elementAt(1)).split("[\\s]+")[3])
					|| (v.elementAt(1)).split("[\\s]+")[1].equals((v.elementAt(1)).split("[\\s]+")[2])
					|| (v.elementAt(1)).split("[\\s]+")[1].equals((v.elementAt(1)).split("[\\s]+")[3])) {
				return true;
			}
		}

		return false;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		boolean include_mutations = true;

		String usage = "\n"	+ "park.phenoStruct requires 0-1 arguments:\n"
										+ "   (1) include mutations? (i.e. muts=" + include_mutations + " (default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.indexOf("nomuts") != -1) {
				include_mutations = false;
				numArgs--;
			} else if (arg.startsWith("muts=")) {
				if (arg.split("=")[1].toLowerCase().equals("true")) {
					include_mutations = true;
					numArgs--;
				} else if (arg.split("=")[1].toLowerCase().equals("false")) {
					include_mutations = false;
					numArgs--;
				} else {
					System.err.println("Error in syntax - '"	+ arg.split("=")[1]
															+ "' is not a valid flag for nomuts (use true/false)");
				}
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length == 0) {
			System.err.println("Using defaults (mutations "
													+ (include_mutations ? "included" : "excluded") + ")");
		}

		try {
			new phenoStruct(include_mutations);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
