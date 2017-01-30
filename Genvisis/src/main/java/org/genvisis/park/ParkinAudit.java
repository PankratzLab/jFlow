package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class ParkinAudit {
	public static final String[] ALTS = {tools.PARKIN_DIR};
	public static final int AOO_CUTOFF = 50;
	public static final String mlpa = "MLPA results.xln";
	public static final String seanList = "Sean sequence list.xln";
	public static final String dianeList = "Diane sequence list.xln";
	public static final String seanResults = "Sean sequence results.xln";
	public static final String dianeResults = "Diane sequence results.xln";
	// public static final String mlpa = "MLPA polymorphism results.xln";
	// public static final String seanList = "Sean polymorphism list.xln";
	// public static final String dianeList = "Diane polymorphism list.xln";
	// public static final String seanResults = "Sean polymorphism results.xln";
	// public static final String dianeResults = "Diane polymorphism results.xln";
	public static final String plateInfo = "plate list -one row per individual.xln";
	public static final String cherylSum = "Cheryl summary.xln";
	public static final String cherylGuess = "Cheryl's guess.xln";
	public static final String inLab = "DNAs in lab.xln";
	public static final String parkinLODs = "lods06_161.dat.xln";
	// public static final String controlInfo = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\tWork\\Global PD files\\control.db.xln";
	public static final String controlInfo = "C:\\Users\\npankrat\\Documents\\1_CRFdb\\Global PD files\\control.db.xln";
	public static final String[] CTRL_HEADER = {"DNA#", "FamID", "IndID", "Index", "UniqueID",
																							"Source", "Affected", "Male", "RaceCode", "Caucasian",
																							"DateOfBirth", "DateOfExam", "AgeAtExam", "Use",
																							"Comment"};
	// public static final String[] KNOWN_POLYMORPHISMS = {"S167N", "Ser167Asn", "V380L", "D394N"};
	public static final String[] KNOWN_POLYMORPHISMS = {"S167N", "Ser167Asn", "V380L", "D394N",
																											"A91A", "C238C", "H279H", "L174L", "L228L",
																											"L261L", "L307L", "P37P"};
	public static final String[][] CATEGORY_CRITERIA = {{	"KnownPolymorphisms", "equals", "S167N",
																												"Ser167Asn", "V380L", "D394N"},
																											{	"SynonymousMissense", "equals", "A91A",
																												"C238C", "H279H", "L174L", "L228L", "L261L",
																												"L307L", "P37P"},
																											{	"BenignMissense", "equals", "R33Q", "A82E",
																												"T83A", "M192V", "E310D", "T240M", "R256C",
																												"P437L"},
																											{	"DeleteriousMissense", "equals", "T173M",
																												"K211N", "R275W", "G430D"},
																											{"Nonsense", "endswith", "x"},
																											{"Frameshift", "contains", "fs"},
																											{"Deletion", "contains", "deletion"},
																											{	"Duplication", "contains", "duplication",
																												"triplication"},
																											{"Indel", "contains", "dup"},
																											{"SpliceSite", "contains", "IVS"}};
	public static final String[] MUT_CLASSES = {"nada", "het", "homo", "comp het", "odd"};
	public static final String[] DIANE_LIST_HEADER = {"DNA # ", "Fam. and Ind. #",
																										"Sequencing Plate"};
	public static final String[] DIANE_RESULTS_HEADER = {	"Exon", "DNA #", "Family and Ind. #",
																												"Nucleotide Change", "Amino Acid change"};
	public static final String[] SEAN_RESULTS_HEADER = {"Family", "ID", "UniqueID", "nummuts",
																											"AR lod", "Age of Onset", "dx",
																											"unique mutation", "mutation", "Exon", "type",
																											"Nucleotide Change", "Amino Acid Change",
																											"Homo/Hetero"};
	public static final String[] MLPA_HEADER = {"DNA #", "Fam.& Ind. #", "Exon 2", "Exon 3", "Exon 4",
																							"Exon 5", "Exon 6", "Exon 7", "Exon 8", "Exon 9",
																							"Exon 10", "Exon 11", "Exon 12"};
	public static final String[] PLATE_HEADER = {	"FamID", "IndID", "UniqueID", "time 1", "time 2",
																								"time 3"};
	public static final String[] LODS_HEADER = {"Family", "linearLOD", "NPL", "AD LOD", "AR LOD", "",
																							"positive"};
	public static final String[] INLAB_HEADER = {"FamNo", "IndNo", "DNA.", "CountOfIndNo"};

	// Recorded: Sean, Diane, MLPA
	// Change to: Sean, Diane, MLPA
	public static final String[][][][] INTERPRETATIONS = {{	{	{	"Gln57fs (+328 amino acids)",
																															"Gln57fs (+328 amino acids)"},
																														{}, {"deletion 3", "deletion 2-4"}},
																													{	{}, {},
																														{"deletion 2-3", "deletion 3-4"}}},
																												{	{	{	"Q57fs (+328 amino acids)",
																															"Q57fs (+328 amino acids)"},
																														{}, {"deletion 3", "deletion 2-4"}},
																													{	{}, {},
																														{"deletion 2-3", "deletion 3-4"}}},
																												{	{	{	"Gln57fs (+35 amino acids)",
																															"Gln57fs (+35 amino acids)"},
																														{}, {"deletion 3-4", "deletion 3-4"}},
																													{	{}, {},
																														{"deletion 3-4", "deletion 3-4"}}},
																												{	{	{	"Q57fs (+35 amino acids)",
																															"Q57fs (+35 amino acids)"},
																														{}, {"deletion 3-4", "deletion 3-4"}},
																													{	{}, {},
																														{"deletion 3-4", "deletion 3-4"}}},
																												{	{	{}, {"G430D"},
																														{"duplication 2-4", "duplication 2-4"}},
																													{{}, {"G430D"}, {"triplication 2-4"}}},
																												{	{	{}, {},
																														{"duplication 2-4", "duplication 2-4"}},
																													{{}, {}, {"triplication 2-4"}}},
																												{	{	{}, {"T240M", "T240M"},
																														{"deletion 5-6"}},
																													{{}, {"T240M"}, {"deletion 5-6"}}},
																												{	{{	"P113fs (+51 amino acids)",
																														"P113fs (+51 amino acids)"},
																													{}, {"deletion 3-4"}},
																													{	{"P113fs (+51 amino acids)"}, {},
																														{"deletion 3-4"}}},
																												{	{{},
																													{	"P113fs (+51 amino acids)",
																														"P113fs (+51 amino acids)"},
																													{"deletion 2-3"}},
																													{	{}, {"P113fs (+51 amino acids)"},
																														{"deletion 2-3"}}}};

	public static void audit() throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, trav;
		String temp;
		Hashtable<String, Vector<String>> famIDs = new Hashtable<String, Vector<String>>();
		Hashtable<String, MutationCarrier> indData = new Hashtable<String, MutationCarrier>();
		Hashtable<String, String> diagnoses = new Hashtable<String, String>(),
				lods = new Hashtable<String, String>();;
		Vector<String> inds = new Vector<String>(), fams = new Vector<String>(), members;
		CheckIDsAgainstDNAs idCheck = new CheckIDsAgainstDNAs();
		MutationCarrier mc, mem;
		int[] parkinExons;
		// int numControlsSequenced = 0, numControlsWithMLPA = 0;
		String uniqueID;
		int size;

		diagnoses = tools.getBestPDdx();

		try {
			reader = Files.getReader(parkinLODs, ALTS);
			ext.checkHeader(reader.readLine().split("\t"), LODS_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				lods.put(line[0], line[2]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + parkinLODs + "\" not found in current directory");
			System.err.println("       skipping linkage test of inclusion");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + parkinLODs + "\"");
			System.exit(2);
		}

		try {
			reader = tools.getNinfoReader(1, false);
			ext.checkHeader(reader.readLine().split("\t"), tools.NINFO1_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				mc = new MutationCarrier(line[0], line[1]);
				mc.AOO = line[2].equals(".") ? -1 : Integer.parseInt(line[2]);
				if (mc.AOO > 0 & mc.AOO < AOO_CUTOFF) {
					mc.shouldBe = true;
				}
				if (lods != null && lods.containsKey(line[0])) {
					mc.lod = Double.parseDouble(lods.get(line[0]));
				}
				indData.put(mc.UniqueID, mc);
				inds.add(mc.UniqueID);
				if (famIDs.containsKey(mc.FamID)) {
					members = famIDs.get(mc.FamID);
				} else {
					famIDs.put(mc.FamID, members = new Vector<String>());
					fams.add(mc.FamID);
				}
				members.add(mc.UniqueID);
			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error reading ninfo1 file");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(controlInfo));
			ext.checkHeader(reader.readLine().split("\t"), CTRL_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				mc = new MutationCarrier(line[0]);
				mc.AOO = -1;
				indData.put(mc.UniqueID, mc);
				inds.add(mc.UniqueID);
				if (famIDs.containsKey(mc.FamID)) {
					members = famIDs.get(mc.FamID);
				} else {
					famIDs.put(mc.FamID, members = new Vector<String>());
					fams.add(mc.FamID);
				}
				members.add(mc.UniqueID);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: could not find file \"" + controlInfo + "\"");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error processing file \"" + controlInfo + "\"");
			ioe.printStackTrace();
			System.exit(2);
		}

		try {
			reader = Files.getReader(plateInfo, ALTS);
			ext.checkHeader(reader.readLine().split("\t"), PLATE_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				uniqueID = unique(line[0], line[1]);
				if (indData.containsKey(uniqueID)) {
					indData.get(uniqueID).firstPlate = Integer.parseInt(line[3]);
				} else {
					System.err.println("Error - a nonexistent person was plated: " + line[0] + "-" + line[1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + plateInfo + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + plateInfo + "\"");
			System.exit(2);
		}

		try {
			reader = Files.getReader(inLab, ALTS);
			ext.checkHeader(reader.readLine().split("\t"), INLAB_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				uniqueID = unique(line[0], line[1]);
				if (indData.containsKey(uniqueID)) {
					indData.get(uniqueID).inLab = true;
				} else {
					System.err.println("Error - a nonexistent person was sent to lab: "	+ line[0] + "-"
															+ line[1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + inLab + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + inLab + "\"");
			System.exit(2);
		}

		try {
			reader = Files.getReader(seanList, ALTS);
			ext.checkHeader(reader.readLine().split("\t"), new String[] {"FamID", "IndID"}, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				uniqueID = unique(line[0], line[1]);
				if (indData.containsKey(uniqueID)) {
					indData.get(uniqueID).seanSeqd = true;
				} else {
					System.err.println("Error - Sean sequenced a nonexistent person: "	+ line[0] + "-t"
															+ line[1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + seanList + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + seanList + "\"");
			System.exit(2);
		}

		try {
			reader = Files.getReader(dianeList, ALTS);
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.startsWith(DIANE_LIST_HEADER[0]));
			if (!reader.ready()) {
				System.err.println("Error - Altered header for Diane's list of who was sequenced");
			}
			ext.checkHeader(temp.split("\t"), DIANE_LIST_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\\t", -1);
				if (line[0].equals("Seans_choice")
						&& !idCheck.checkPair(line[1], line[0], false).equals("")) {
					System.err.println("  (from " + dianeList + ")");
					// System.err.print(dianeList+"\t");
				}
				// if (line[0].startsWith("ND")) {
				// numControlsSequenced++;
				// }
				if (line[1].equals("")) {
					uniqueID = line[0];
				} else {
					trav = tools.getFamID(line[1]);
					uniqueID = unique(trav[0], trav[1]);
				}
				if (indData.containsKey(uniqueID)) {
					mc = indData.get(uniqueID);
					mc.dianeSeqd = true;
					mc.dna = line[0];
				} else {
					System.err.println("Error - Diane sequenced a nonexistent person: " + line[1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dianeList + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dianeList + "\"");
			System.exit(2);
		}

		try {
			reader = Files.getReader(seanResults, ALTS);
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.startsWith(SEAN_RESULTS_HEADER[0]));
			if (!reader.ready()) {
				System.err.println("Error - Altered header for Sean's sequence results");
			}
			ext.checkHeader(temp.split("\t"), SEAN_RESULTS_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				uniqueID = unique(line[0], line[1]);
				if (indData.containsKey(uniqueID)) {
					mc = indData.get(uniqueID);
					if (line[8].equals("point")) {
						for (int i = 0; i < (line.length > 13 && line[13].startsWith("Homo") ? 2 : 1); i++) {
							if (line.length > 12 && ext.indexOfStr(line[12], KNOWN_POLYMORPHISMS) != -1) {
								mc.knownPolymorphisms = ArrayUtils.addStrToArray(line[12], mc.knownPolymorphisms);
							} else {
								mc.seanMutations = ArrayUtils.addStrToArray(line.length > 12	&& !line[12].equals("")
																																																	? line[12]
																																																: line[11],
																												mc.seanMutations);
							}
						}
					}
				} else {
					System.err.println("Error - Sean results for a nonexistent person: "	+ line[0] + "-"
															+ line[1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + seanResults + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + seanResults + "\"");
			System.exit(2);
		}

		try {
			reader = Files.getReader(dianeResults, ALTS);
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.startsWith(DIANE_RESULTS_HEADER[0]));
			if (!reader.ready()) {
				System.err.println("Error - Altered header for Diane's sequence results");
			}
			ext.checkHeader(temp.split("\t"), DIANE_RESULTS_HEADER, true);
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split("\t");
				if (!idCheck.checkPair(line[2], line[1], false).equals("")) {
					System.err.println("  (from " + dianeResults + ": " + line[2] + " and " + line[1] + ")");
					// System.err.print(dianeResults+"\t");
				}
				if (line[2].equals("")) {
					uniqueID = line[1];
				} else {
					trav = tools.getFamID(line[2]);
					if (!idCheck.checkPair(line[2], line[1], false).equals("")) {
						System.err.println("  (from "	+ dianeResults + ": " + line[2] + " and " + line[1]
																+ ")");
					}
					uniqueID = unique(trav[0], trav[1]);
				}
				if (indData.containsKey(uniqueID)) {
					mc = indData.get(uniqueID);
					for (int i = 0; i < (temp.indexOf("Homozygous") > 0 ? 2 : 1); i++) {
						if (line.length > 4 && ext.indexOfStr(line[4], KNOWN_POLYMORPHISMS) != -1) {
							mc.knownPolymorphisms = ArrayUtils.addStrToArray(line[4], mc.knownPolymorphisms);
						} else {
							mc.dianeMutations =
																ArrayUtils.addStrToArray(line.length > 4	&& !line[4].equals("")
																																															? line[4]
																																														: line[3],
																										mc.dianeMutations);
						}
					}
					if (mc.dna.equals("")) {
						mc.dna = line[1];
					} else if (!mc.dna.equals(line[1])) {
						System.err.println("Different DNAs used between methods for "	+ mc.FamID + "-"
																+ mc.IndID + " (" + mc.dna + " and " + line[1] + ")");
						mc.dna += "/" + line[1];
					}
				} else {
					System.err.println("Error - Diane results for a nonexistent person: "	+ line[1] + "/"
															+ line[2]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dianeResults + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dianeResults + "\"");
			System.exit(2);
		}

		try {
			reader = Files.getReader(mlpa, ALTS);
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.startsWith(MLPA_HEADER[0]));
			if (!reader.ready()) {
				System.err.println("Error - Altered header for MLPA results");
			}
			ext.checkHeader(temp.split("\t"), MLPA_HEADER, true);
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split("\t", -1);
				if (line.length > 1) {
					if (line[1].equals("")) {
						uniqueID = line[0];
					} else {
						trav = tools.getFamID(line[1]);
						if (!idCheck.checkPair(line[1], line[0], false).equals("")) {
							System.err.println("  (from " + mlpa + ": " + line[1] + " and " + line[0] + ")");
							// System.err.print(mlpa+"\t");
						}
						uniqueID = unique(trav[0], trav[1]);
					}
					if (indData.containsKey(uniqueID)) {
						mc = indData.get(uniqueID);
						mc.MLPAd = true;
						if (mc.dna.equals("")) {
							mc.dna = line[0];
						} else if (!mc.dna.equals(line[0]) && !line[0].equals("")) {
							System.err.println("Different DNAs used between methods for "	+ mc.FamID + "-"
																	+ mc.IndID + " (" + mc.dna + " and " + line[0] + ")");
							mc.dna += "/" + line[0];
						}
						if (temp.substring(line[0].length() + 1 + line[1].length()).trim().length() > 0) {
							line = temp.split("\t", -1);
							if (temp.indexOf("Deletion") > 0 && temp.indexOf("Duplication") > 0) {
								System.err.println("Warning - Both a deletion and a duplication found in individual "
																		+ mc.FamID + "-" + mc.IndID);
							}
							if (temp.indexOf("Deletion") == -1 && temp.indexOf("Duplication") == -1) {
								System.err.println("Warning - MLPA result that is neither a deletion nor a duplication for individual "
																		+ mc.FamID + "-" + mc.IndID);
							}
							parkinExons = new int[12];
							try {
								for (int i = 1; i < (line.length > 12 ? 12 : line.length - 1); i++) {
									if (line[1 + i].equals("")) {
										parkinExons[i] = 0;
									} else if (line[1 + i].equals("Duplication")) {
										parkinExons[i] = 1;
									} else if (line[1 + i].equals("Homo. Duplication")
															|| line[1 + i].equals("Hom. Duplication")) {
										parkinExons[i] = 2;
									} else if (line[1 + i].equals("Deletion")) {
										parkinExons[i] = -1;
									} else if (line[1 + i].equals("Homo. Deletion")
															|| line[1 + i].equals("Hom. Deletion")) {
										parkinExons[i] = -2;
									} else {
										System.err.println("Error - invalid mutation for "	+ mc.FamID + "-" + mc.IndID
																				+ " in MLPA file: " + line[1 + i]);
									}
								}
							} catch (Exception e) {
								System.err.println("Error - invalid size for "	+ mc.FamID + "-" + mc.IndID
																		+ " in MLPA file");
							}
							trav = parseMLPAMutations(parkinExons);
							if (trav.length == 0) {
								System.err.println("Error - parse zero mutations for individual "	+ mc.FamID + "-"
																		+ mc.IndID + " (the file format expected some)");
							}
							for (String element : trav) {
								mc.MLPAresults = ArrayUtils.addStrToArray(element, mc.MLPAresults);
							}
							// System.err.println(mc.FamID+"-"+mc.IndID+"\t"+Array.toStr(mc.MLPAresults,
							// "\t\t") +"\t\t"+Array.toStr(parkinExons, " "));

						}
					} else {
						System.err.println("Error - MLPA results for a nonexistent person: "	+ line[0] + "/"
																+ line[1]);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + mlpa + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + mlpa + "\"");
			System.exit(2);
		}

		for (int i = 0; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i));

			transformMe(mc);

			for (int j = 0; j < 2; j++) {
				if (j == 0) {
					trav = mc.seanMutations.length > mc.dianeMutations.length	? mc.seanMutations
																																		: mc.dianeMutations;
				} else {
					trav = mc.MLPAresults;
				}
				for (String element : trav) {
					size = Matrix.getSize(mc.categories);
					for (int l = 0; l < CATEGORY_CRITERIA.length; l++) {
						if (CATEGORY_CRITERIA[l][1].equals("equals")) {
							for (int m = 2; m < CATEGORY_CRITERIA[l].length; m++) {
								if (element.toLowerCase().equals(CATEGORY_CRITERIA[l][m].toLowerCase())) {
									mc.categories[l] = ArrayUtils.addStrToArray(element, mc.categories[l]);
									m = CATEGORY_CRITERIA[l].length - 1;
									l = CATEGORY_CRITERIA.length - 1;
								}
							}
						} else if (CATEGORY_CRITERIA[l][1].equals("contains")) {
							for (int m = 2; m < CATEGORY_CRITERIA[l].length; m++) {
								if (element.toLowerCase().contains(CATEGORY_CRITERIA[l][m].toLowerCase())) {
									mc.categories[l] = ArrayUtils.addStrToArray(element, mc.categories[l]);
									m = CATEGORY_CRITERIA[l].length - 1;
									l = CATEGORY_CRITERIA.length - 1;
								}
							}
						} else if (CATEGORY_CRITERIA[l][1].equals("endswith")) {
							for (int m = 2; m < CATEGORY_CRITERIA[l].length; m++) {
								if (element.toLowerCase().endsWith(CATEGORY_CRITERIA[l][m].toLowerCase())) {
									mc.categories[l] = ArrayUtils.addStrToArray(element, mc.categories[l]);
									m = CATEGORY_CRITERIA[l].length - 1;
									l = CATEGORY_CRITERIA.length - 1;
								}
							}
						} else {
							System.err.println("Error - Don't know how to process category with operator: "
																	+ CATEGORY_CRITERIA[l][1]);
							System.exit(1);
						}
					}
					if (Matrix.getSize(mc.categories) == size) {
						System.err.println("Error - Could not find a category for: " + element);
					}
				}
			}

			mc.fullyChecked = (mc.seanSeqd || mc.dianeSeqd) && mc.MLPAd;
			mc.partiallyChecked = mc.seanSeqd || mc.dianeSeqd || mc.MLPAd;
			mc.numMutations =
											(mc.seanMutations.length > mc.dianeMutations.length	? mc.seanMutations.length
																																					: mc.dianeMutations.length)
												+ mc.MLPAresults.length;
			mc.numSeq = mc.seanMutations.length > mc.dianeMutations.length	? mc.seanMutations.length
																																			: mc.dianeMutations.length;
			mc.numDosage = mc.MLPAresults.length;

			if (mc.numMutations == 0) {
				mc.mutClass = 0;
			} else if (mc.numMutations == 1) {
				mc.mutClass = 1;
			} else {
				for (int j = 0; j < 3; j++) {
					line = j == 0 ? mc.seanMutations : (j == 1 ? mc.dianeMutations : mc.MLPAresults);
					if (line.length == 1 && mc.mutClass == 0) {
						mc.mutClass = 1;
					}
					if (line.length == 1 && mc.mutClass == 1) {
						mc.mutClass = 3;
					}
					if (line.length == 2) {
						if (line[0].equals(line[1])) {
							mc.mutClass = mc.mutClass == 0 ? 2 : 4;
						} else {
							mc.mutClass = mc.mutClass == 0 ? 3 : 4;
						}
					}
				}
			}

			members = famIDs.get(mc.FamID);
			mc.numberFamMems = members.size() - 1;
			for (int j = 0; j < members.size(); j++) {
				mem = indData.get(members.elementAt(j));
				if (mem != mc && mem.fullyChecked) {
					mc.familyMemberFullyChecked = true;
				}
				if (mem != mc && mem.partiallyChecked) {
					mc.familyMemberPartiallyChecked = true;
				}
			}
			if (!mc.shouldBe && mc.lod > 0) {
				mc.shouldBe = true;
				for (int j = 0; j < members.size(); j++) {
					mem = indData.get(members.elementAt(j));
					if (mem != mc && mem.shouldBe) {
						mc.shouldBe = false;
					}
				}
			}

			if (!mc.shouldBe	&& !mc.partiallyChecked && mc.lod < 0
					&& diagnoses.containsKey(mc.FamID + "\t" + mc.IndID)
					&& (diagnoses.get(mc.FamID + "\t" + mc.IndID)).equals("VPD")) {
				mc.randomCandidate = true;
				for (int j = 0; j < members.size(); j++) {
					mem = indData.get(members.elementAt(j));
					if (mem != mc && (mem.shouldBe || mem.partiallyChecked || mem.randomCandidate)) {
						mc.randomCandidate = false;
					}
				}
			}

			mc.butWasnt = mc.shouldBe && !mc.fullyChecked;
		}

		// System.out.println("There were "+numControlsWithSequenceMutations+"
		// of "+numControlsSequenced+" sequenced NINDS samples (i.e. \"ND\"
		// samples) with variants");

		writer = new PrintWriter(new FileWriter("park_audit.xln"));
		for (int i = -1; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i == -1 ? 0 : i));
			writer.println((i == -1 ? "UniqueID" : mc.UniqueID)	+ "\t"
											+ (i == -1 ? "DNA#" : (mc.dna.equals("") ? "" : mc.dna)) + "\t"
											+ (i == -1 ? "FamID" : mc.FamID) + "\t" + (i == -1 ? "IndID" : mc.IndID)
											+ "\t" + (i == -1 ? "Sean seq" : (mc.seanSeqd ? "1" : "0")) + "\t"
											+ (i == -1 ? "Sean Found" : mc.seanMutations.length) + "\t"
											+ (i == -1 ? "Diane Seq" : (mc.dianeSeqd ? "1" : "0")) + "\t"
											+ (i == -1 ? "Diane Found" : mc.dianeMutations.length) + "\t"
											+ (i == -1 ? "MLPA'd" : (mc.MLPAd ? "1" : "0")) + "\t"
											+ (i == -1 ? "MLPA found" : mc.MLPAresults.length) + "\t"
											+ (i == -1 ? "Num Mutations" : mc.numMutations) + "\t"
											+ (i == -1 ? "DNA sent" : (mc.inLab ? "1" : "0")) + "\t"
											+ (i == -1 ? "# additional fam members" : mc.numberFamMems) + "\t"
											+ (i == -1 ? "Partially Checked" : (mc.partiallyChecked ? "1" : "0")) + "\t"
											+ (i == -1 ? "Fully Checked" : (mc.fullyChecked ? "1" : "0")) + "\t"
											+ (i == -1	? "Family member fully checked"
																	: (mc.familyMemberFullyChecked ? "1" : "0"))
											+ "\t"
											+ (i == -1 ? "First Plate" : (mc.firstPlate == -1 ? "." : mc.firstPlate))
											+ "\t" + (i == -1 ? "NPL score" : (mc.lod == -999 ? "." : mc.lod)) + "\t"
											+ (i == -1	? "Dx"
																	: (diagnoses.containsKey(mc.FamID + "\t" + mc.IndID)
																																													? diagnoses.get(mc.FamID
																																																					+ "\t"
																																																				+ mc.IndID)
																																												: "."))
											+ "\t" + (i == -1 ? "AOO" : (mc.AOO == -1 ? "." : mc.AOO)) + "\t"
											+ (i == -1 ? "Random Candidate" : (mc.randomCandidate ? "1" : "0")) + "\t"
											+ (i == -1 ? "Should Be" : (mc.shouldBe ? "1" : "0")) + "\t"
											+ (i == -1 ? "But Wasn't" : (mc.butWasnt ? "1\t*" : "0")) + "\t" + "");
		}
		writer.close();

		writer = new PrintWriter(new FileWriter("parkin mutations.xln"));
		writer.println("DNA used\tUniqueID\tFamID\tIndID\tmut_count\ttype\tDx\tAOO\tMutations");
		for (int i = 0; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i == -1 ? 0 : i));
			if (mc.numMutations > 0) {
				writer.print((mc.dna.equals("") ? "unknown" : mc.dna)	+ "\t" + mc.UniqueID + "\t" + mc.FamID
											+ "\t" + mc.IndID);
				writer.print("\t"	+ mc.numMutations + "\t" + MUT_CLASSES[mc.mutClass] + "\t"
											+ (diagnoses.containsKey(mc.FamID + "\t" + mc.IndID)
																																							? diagnoses.get(mc.FamID
																																															+ "\t"
																																														+ mc.IndID)
																																						: ".")
											+ "\t" + (mc.AOO == -1 ? "." : mc.AOO));
				for (String seanMutation : mc.seanMutations) {
					writer.print("\tSean: " + seanMutation);
				}
				for (String dianeMutation : mc.dianeMutations) {
					writer.print("\tDiane: " + dianeMutation);
				}
				for (String mlpAresult : mc.MLPAresults) {
					writer.print("\tMLPA: " + mlpAresult);
				}
				writer.println();
			}
		}
		writer.close();

		writer = new PrintWriter(new FileWriter("parkinCategories.csv"));
		writer.print("DNA used,UniqueID,FamID,IndID,FullyChecked,NumMutations");
		for (String[] element : CATEGORY_CRITERIA) {
			writer.print(",#" + element[0]);
		}
		for (String[] element : CATEGORY_CRITERIA) {
			writer.print(",Just1" + element[0]);
		}
		for (String[] element : CATEGORY_CRITERIA) {
			writer.print("," + element[0]);
		}
		writer.println();
		for (int i = 0; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i == -1 ? 0 : i));
			if (mc.fullyChecked || mc.numMutations > 0 || mc.knownPolymorphisms.length > 0) {
				writer.print((mc.dna.equals("") ? "unknown" : mc.dna)	+ "," + mc.UniqueID + "," + mc.FamID
											+ "," + mc.IndID + "," + (mc.fullyChecked ? 1 : 0) + "," + mc.numMutations);
				for (int j = 0; j < CATEGORY_CRITERIA.length; j++) {
					writer.print("," + mc.categories[j].length);
				}
				for (int j = 0; j < CATEGORY_CRITERIA.length; j++) {
					writer.print("," + (mc.numMutations == 1 ? mc.categories[j].length : 0));
				}
				for (int j = 0; j < CATEGORY_CRITERIA.length; j++) {
					writer.print("," + ArrayUtils.toStr(mc.categories[j], "/"));
				}
				writer.println();
			}
		}
		writer.close();

		writer = new PrintWriter(new FileWriter("parkin polymorphisms.xln"));
		writer.println("DNA used\tUniqueID\tFamID\tIndID\tmut_count\ttype\tDx\tAOO\tnum_poly\tPolymorphisms");
		for (int i = 0; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i == -1 ? 0 : i));
			if (mc.knownPolymorphisms.length > 0) {
				writer.print((mc.dna.equals("") ? "unknown" : mc.dna)	+ "\t" + mc.UniqueID + "\t" + mc.FamID
											+ "\t" + mc.IndID);
				writer.print("\t"	+ mc.numMutations + "\t" + MUT_CLASSES[mc.mutClass] + "\t"
											+ (diagnoses.containsKey(mc.FamID + "\t" + mc.IndID)
																																							? diagnoses.get(mc.FamID
																																															+ "\t"
																																														+ mc.IndID)
																																						: ".")
											+ "\t" + (mc.AOO == -1 ? "." : mc.AOO) + "\t" + mc.knownPolymorphisms.length);
				for (String knownPolymorphism : mc.knownPolymorphisms) {
					writer.print("\t" + knownPolymorphism);
				}
				writer.println();
			}
		}
		writer.close();

		writer = new PrintWriter(new FileWriter("parkin.csv"));
		writer.println("DNA#,FamID,IndID,UniqueID,num_muts,carrier,2muts,1mut,num_knownPoly,numSeq,numDosage,carrierSeq,carrierDosage,justSeqHet,justDosageHet,justTwoMuts,FullyChecked");
		for (int i = 0; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i == -1 ? 0 : i));
			// if (mc.fullyChecked || mc.numMutations>0) {
			if (mc.partiallyChecked || mc.numMutations > 0) {
				writer.println((mc.dna.equals("") ? "." : mc.dna)	+ "," + mc.FamID + "," + mc.IndID + ","
												+ mc.UniqueID + "," + mc.numMutations + "," + (mc.numMutations > 0 ? 1 : 0)
												+ "," + (mc.numMutations > 1 ? 1 : 0) + "," + (mc.numMutations == 1 ? 1 : 0)
												+ "," + mc.knownPolymorphisms.length + "," + mc.numSeq + "," + mc.numDosage
												+ "," + (mc.numSeq > 0 ? 1 : 0) + "," + (mc.numDosage > 0 ? 1 : 0) + ","
												+ (mc.numMutations >= 2 || mc.numDosage >= 1	? "."
																																			: (mc.numSeq == 1 ? 1 : 0))
												+ ","
												+ (mc.numMutations >= 2 || mc.numSeq >= 1	? "."
																																	: (mc.numDosage == 1 ? 1 : 0))
												+ "," + (mc.numMutations == 1 ? "." : (mc.numMutations >= 2 ? 1 : 0)) + ","
												+ (mc.fullyChecked ? 1 : 0));
			}
		}
		writer.close();

		writer = new PrintWriter(new FileWriter("ninfo3_parkin_data.txt"));
		for (int i = 0; i < inds.size(); i++) {
			mc = indData.get(inds.elementAt(i == -1 ? 0 : i));
			if (mc.numMutations >= 2) {
				writer.print("MUT\t" + mc.FamID + "\t" + mc.IndID);
				writer.print("\t" + MUT_CLASSES[mc.mutClass] + " parkin mutation");
				for (String seanMutation : mc.seanMutations) {
					writer.print(" (" + seanMutation + ")");
				}
				for (String dianeMutation : mc.dianeMutations) {
					writer.print(" (" + dianeMutation + ")");
				}
				for (String mlpAresult : mc.MLPAresults) {
					writer.print(" (" + mlpAresult + ")");
				}
				writer.println();
			}
		}
		writer.close();

	}

	public static String unique(String famID, String indID) throws IOException {
		if (famID.startsWith("ND")) {
			return famID;
		} else {
			try {
				return Integer.parseInt(famID) * 1000 + Integer.parseInt(indID) + "";
			} catch (NumberFormatException nfe) {
				System.err.println("Error - could not form a unique identifier with Family ID = '"	+ famID
														+ "' and Individual ID = '" + indID + "'");
				return "error";
			}
		}
	}

	public static String[] parseMLPAMutations(int[] exons) throws IOException {
		String[] mutations = new String[0];
		String trav1, trav2;
		int count = 0;

		trav1 = trav2 = "";
		while (count < exons.length) {
			if (!trav2.equals("") && Math.abs(exons[count]) != 2) {
				if (trav2.endsWith("-")) {
					trav2 += count;
				}
				mutations = ArrayUtils.addStrToArray(trav2, mutations);
				trav2 = "";
			}
			if (!trav1.equals("") && exons[count] == 0) {
				if (trav1.endsWith("-")) {
					trav1 += count;
				}
				mutations = ArrayUtils.addStrToArray(trav1, mutations);
				trav1 = "";
			}
			if (Math.abs(exons[count]) == 2) {
				if (trav2.equals("")) {
					trav2 = (exons[count] < 0 ? "deletion" : "duplication") + " " + (count + 1);
				} else if (!trav2.endsWith("-")) {
					trav2 += "-";
				}
			}
			if (Math.abs(exons[count]) >= 1) {
				if (trav1.equals("")) {
					trav1 = (exons[count] < 0 ? "deletion" : "duplication") + " " + (count + 1);
				} else if (!trav1.endsWith("-")) {
					trav1 += "-";
				}
			}
			count++;
		}

		if (mutations.length > 2) {
			System.err.println("Error - more than 2 mutations parsed for: " + ArrayUtils.toStr(exons, " "));
		}
		return mutations;
	}

	public static void transformMe(MutationCarrier mc) throws IOException {
		for (String[][][] element : INTERPRETATIONS) {
			if (eqArrays(mc.seanMutations, element[0][0])	&& eqArrays(mc.dianeMutations, element[0][1])
					&& eqArrays(mc.MLPAresults, element[0][2])) {

				mc.seanMutations = element[1][0].clone();
				mc.dianeMutations = element[1][1].clone();
				mc.MLPAresults = element[1][2].clone();
			}
		}
	}

	public static boolean eqArrays(Object[] arr1, Object[] arr2) {
		if (arr1.length != arr2.length) {
			return false;
		}
		for (int i = 0; i < arr1.length; i++) {
			if (!arr1[i].equals(arr2[i])) {
				return false;
			}
		}

		return true;
	}

	public static void main(String[] args) throws IOException {
		try {
			audit();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}


class MutationCarrier {
	public String UniqueID;
	public String FamID;
	public String IndID;
	public int AOO;
	public int firstPlate;
	public double lod;
	public boolean inLab;
	public int numberFamMems;
	public boolean seanSeqd;
	public boolean dianeSeqd;
	public boolean MLPAd;
	public String[] seanMutations;
	public String[] dianeMutations;
	public String[] knownPolymorphisms;
	public String[] MLPAresults;
	public String[][] categories;
	public int numMutations;
	public int numSeq;
	public int numDosage;
	public int mutClass;
	public boolean partiallyChecked;
	public boolean fullyChecked;
	public boolean familyMemberFullyChecked;
	public boolean familyMemberPartiallyChecked;
	public boolean randomCandidate;
	public boolean shouldBe;
	public boolean butWasnt;
	public String dna;
	public Vector<String> mutations;

	public MutationCarrier() {
		firstPlate = -1;
		lod = -999;
		seanMutations = new String[0];
		dianeMutations = new String[0];
		knownPolymorphisms = new String[0];
		MLPAresults = new String[0];
		categories = new String[ParkinAudit.CATEGORY_CRITERIA.length][0];
		dna = "";
	}

	public MutationCarrier(String newUniqueID) {
		this();
		UniqueID = newUniqueID;
		try {
			FamID = Integer.parseInt(newUniqueID.substring(0, 5)) + "";
			IndID = Integer.parseInt(newUniqueID.substring(5, 8)) + "";
		} catch (NumberFormatException nfe) {
			FamID = newUniqueID;
			IndID = newUniqueID;
		}
	}

	public MutationCarrier(String newFamID, String newIndID) {
		this();
		FamID = Integer.parseInt(newFamID) + "";
		IndID = Integer.parseInt(newIndID) + "";
		try {
			UniqueID = (Integer.parseInt(FamID) * 1000 + Integer.parseInt(IndID)) + "";
		} catch (NumberFormatException nfe) {
			UniqueID = FamID;
			if (!newFamID.equals(newIndID)) {
				System.err.println("Warning - new MutationCarrier has either a FamID ("	+ newFamID
														+ ") or an IndID (" + newIndID
														+ ") that is non-numeric, and yet not identical. Declaring UniqueID = IndID.");
			}
		}
	}
}
