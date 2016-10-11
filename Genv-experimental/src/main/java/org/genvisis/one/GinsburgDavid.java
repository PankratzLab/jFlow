package org.genvisis.one;

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
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.link.LinkageFormat;
import org.genvisis.link.LinkageMap;
import org.genvisis.link.TrimFam;
import org.genvisis.link.bat.Mendel;

public class GinsburgDavid {
	public static final String LOG_FILE = "logfile.txt";
	public static final String DEFAULT_ROOT =
																					"C:\\Documents and Settings\\npankrat\\My Documents\\murine\\";
	// public static final String DEFAULT_ROOT = "C:\\transfer\\forSeattle\\murine\\";

	// public static final String DEFAULT_DIR = "DBALine\\doublecheck\\";
	// public static final String DEFAULT_DIR = "line1\\2008.08.04\\";
	// public static final String DEFAULT_DIR = "line1\\2008.08.27\\";
	// public static final String DEFAULT_DIR = "line1\\2008.08.27\\exploreModels\\";

	// public static final String DEFAULT_DIR = "line2\\2008.03.27\\";
	// public static final String DEFAULT_DIR = "line2\\2008.03.27\\exploreModel\\";
	// public static final String DEFAULT_DIR = "line2\\2008.07.29\\";
	// public static final String DEFAULT_DIR = "line2\\2008.07.29\\exploreModels\\";
	// public static final String DEFAULT_DIR = "line2\\2008.08.27\\";
	// public static final String DEFAULT_DIR = "line2\\2008.08.27\\exploreModels\\";
	// public static final String DEFAULT_DIR = "line2\\2009.03.11\\";

	// public static final String DEFAULT_DIR = "line3\\2008.07.17\\";
	// public static final String DEFAULT_DIR = "line3\\2008.07.17\\fineTuningModel\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.07\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.07\\deleteMarkers\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.07\\exploreModels\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.07\\chr12markers\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.27\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.27\\subset\\";
	// public static final String DEFAULT_DIR = "line3\\2008.08.27\\exploreModels\\";
	// public static final String DEFAULT_DIR = "line3\\2008.10.06\\";

	// public static final String DEFAULT_DIR = "lines2_3\\";

	// public static final String DEFAULT_DIR = "line4\\2008.07.29\\";
	// public static final String DEFAULT_DIR = "line4\\2008.07.29\\exploreResults\\";
	// public static final String DEFAULT_DIR = "line4\\2008.08.27\\";

	// public static final String DEFAULT_DIR = "line5\\2008.08.04\\";
	// public static final String DEFAULT_DIR = "line5\\2008.08.27\\";
	// public static final String DEFAULT_DIR = "line5\\2008.08.27\\exploreModels\\";
	// public static final String DEFAULT_DIR = "line5\\2009.03.11\\";

	// public static final String DEFAULT_DIR = "line6E\\2009.06.23\\";

	// public static final String DEFAULT_DIR = "line9\\2008.08.04\\";

	// public static final String DEFAULT_DIR = "line10\\2008.08.04\\";
	// public static final String DEFAULT_DIR = "line10\\2008.08.27\\";

	// public static final String DEFAULT_DIR = "line10E\\2009.06.12\\";
	// public static final String DEFAULT_DIR = "line10E\\2009.06.16\\";
	// public static final String DEFAULT_DIR = "line10E\\2009.11.06\\";

	// public static final String DEFAULT_DIR = "yossi\\2009.03.11\\";
	// public static final String DEFAULT_DIR = "yossi\\2009.03.11reversed\\";
	// public static final String DEFAULT_DIR = "yossi\\2010.02.03\\";
	// public static final String DEFAULT_DIR = "yossi\\2010.02.03reversed\\";
	// public static final String DEFAULT_DIR = "yossi\\2010.02.03mortality\\";
	// public static final String DEFAULT_DIR = "yossi\\2010.02.03mortalityModifier\\";
	public static final String DEFAULT_DIR = "yossi\\2010.02.03\\quantitative\\";


	public static final String[] PUNCTUATION = {"/", "\\", "(", "-", " "};
	public static final String[] IGNORES = {"%", "Total"};
	public static final String TRAIT_STRAIN_PREFIX = "999";
	public static final String BACKCROSS_STRAIN_PREFIX = "777";
	public static final String[] GENOTYPIC_HEADER = {"marker", "chr", "pos"};

	public static final double DEFAULT_FREQ = 0.01;
	public static final double[] DEFAULT_MODEL1 = {0.03, 0.80, 0.80};
	public static final double[] DEFAULT_MODEL2 = {0.03, 0.03, 0.80};
	// public static final double[] DEFAULT_MODEL = {0.0, 1.0, 1.0};
	// public static final double[] DEFAULT_MODEL = {0.0, 0.98, 0.98};
	// public static final double[] DEFAULT_MODEL = {0.0, 0.0, 0.98};

	public static final int EXPLORE_LIMIT = 12;

	// public static final double[] EXPLORE_PEN = {1.00, 0.99, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50,
	// 0.40, 0.30, 0.20, 0.10};
	// public static final double[] EXPLORE_HET = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15,
	// 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};

	public static final String DEFAULT_EXPLORE_FILE = "explore.txt";
	public static final String DEFAULT_EXPLORE_DIR = "exploreModels/";
	public static final double[] EXPLORE_PEN = {1.00, 0.99, 0.95, 0.90, 0.80, 0.70, 0.60};
	public static final double[] EXPLORE_HET = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20,
																							0.25};
	// public static final double[] EXPLORE_PEN = {1.00, 0.80};
	// public static final double[] EXPLORE_HET = {0.00, 0.03};

	public static void makePedigree(String dir, String pedigree, String genotypes) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, subline;
		Hashtable<String, String> lookup = new Hashtable<String, String>();
		Vector<String> dnaList = new Vector<String>();
		String temp, unparseable;
		int pedNum, mother, father, sex;
		String sStrain, bStrain;
		int sStrainCount, bStrainCount, numGenotyped;

		dnaList = getDNAlist(dir, genotypes, true);
		for (int i = 0; i < dnaList.size(); i++) {
			temp = dnaList.remove(i);
			if (dnaList.contains(temp)) {
				logIt("Warning - '" + temp + "' is duplicated, please collapse", true);
				try {
					writer = new PrintWriter(new FileWriter(dir + "DUPLICATE DNAS!!!.txt", true));
					writer.println(temp);
					writer.close();
				} catch (IOException ioe) {
					logIt("Error writing to file \"DUPLICATE DNAS!!!.txt\"", true);
				}
			}
			dnaList.insertElementAt(temp, i);
		}
		sStrain = dnaList.elementAt(0);
		bStrain = dnaList.elementAt(1);
		logIt("Assuming we're taking the F1s of a "	+ sStrain + " x " + bStrain
					+ " cross and backcrossing into " + bStrain, false);

		try {
			reader = new BufferedReader(new FileReader(dir + pedigree));
			reader.readLine();
			sStrainCount = bStrainCount = 0;
			unparseable = "";
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t");
				if (lookup.containsKey(line[1])) {
					logIt("Error - You have more than two individuals with the same ID ('"	+ line[1] + "')",
								true);
				}
				if (line[6].contains(sStrain) && line[6].contains(bStrain)) {
					logIt("Error - '"	+ line[6] + "' contains identifiers for both strains ('" + sStrain
								+ "' and '" + bStrain + "')", true);
					logIt("failed at " + ext.getTime(), true);
					System.exit(1);
				} else if (line[6].contains(sStrain)) {
					lookup.put(line[1], TRAIT_STRAIN_PREFIX + ext.formNum(++sStrainCount, 2));
				} else if (line[6].contains(bStrain)) {
					lookup.put(line[1], BACKCROSS_STRAIN_PREFIX + ext.formNum(++bStrainCount, 2));
				} else {
					temp = ext.removeQuotes(line[6]);
					subline = temp.trim().split("/");
					numGenotyped = 0;
					for (int i = 0; i < subline.length; i++) {
						for (String element : PUNCTUATION) {
							if (subline[i].contains(element) && i < 1) {
								unparseable += line[6] + "; ";
								subline[i] = subline[i].substring(0, temp.indexOf(element));
							}
						}
						if (dnaList.contains(subline[i])) {
							numGenotyped++;
						}
					}
					if (numGenotyped > 1) {
						logIt("Warning - make sure you check the duplicate sample '" + line[6] + "'", false);
					}
					if (subline.length > 1) {
						temp = temp.substring(0, temp.indexOf("/"));
					}
					if (temp.startsWith(TRAIT_STRAIN_PREFIX)) {
						logIt("Error - '"	+ temp + "' (parsed from '" + line[6]
									+ "') uses the reserved susceptibility strain prefix", true);
						logIt("failed at " + ext.getTime(), true);
						System.exit(1);
					}
					if (temp.startsWith(BACKCROSS_STRAIN_PREFIX)) {
						logIt("Error - '"	+ temp + "' (parsed from '" + line[6]
									+ "') uses the reserved backcross strain prefix", true);
						logIt("failed at " + ext.getTime(), true);
						System.exit(1);
					}

					try {
						Integer.parseInt(temp);
					} catch (NumberFormatException nfe) {
						logIt("Error - failed to parse strain from '"	+ temp + "' (for indiviudal " + line[1]
									+ "); not either of the ancestral strains (" + sStrain + " or " + bStrain + ")",
									true);
						logIt("failed at " + ext.getTime(), true);
						System.exit(1);
					}

					lookup.put(line[1], temp);
				}
			}
			reader.close();
			if (unparseable.length() > 0) {
				logIt("Error - are any of the following either the susceptible strain or the backcross strain?",
							true);
				logIt(unparseable, true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(1);
			}
		} catch (FileNotFoundException fnfe) {
			logIt("Error: file \"" + dir + pedigree + "\" not found in current directory", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		} catch (IOException ioe) {
			logIt("Error reading file \"" + dir + pedigree + "\"", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(dir + pedigree));
			writer = new PrintWriter(new FileWriter(dir + "untrimmed.pre"));
			line = reader.readLine().trim().split("\t");
			pedNum = ext.indexOfStr("Ped Num", line);
			if (pedNum == -1) {
				logIt("Error - could not find 'Ped Num' column in " + dir + pedigree, true);
			}
			mother = ext.indexOfStr("Mother", line);
			if (mother == -1) {
				logIt("Error - could not find 'Mother' column in " + dir + pedigree, true);
			}
			father = ext.indexOfStr("Father", line);
			if (father == -1) {
				logIt("Error - could not find 'Father' column in " + dir + pedigree, true);
			}
			sex = ext.indexOfStr("Sex", line);
			if (sex == -1) {
				logIt("Error - could not find 'Sex' column in " + dir + pedigree, true);
			}
			while (reader.ready()) {
				line = reader.readLine().split("\t", -1);
				writer.println("1\t"	+ lookup.get(line[pedNum]) + "\t"
												+ (line[father].equals("") ? "0" : lookup.get(line[father].trim())) + "\t"
												+ (line[mother].equals("") ? "0" : lookup.get(line[mother].trim())) + "\t"
												+ (line[sex].equals("m") ? "1" : (line[sex].equals("f") ? "2" : "0")) + "\t"
												+ (lookup.get(line[pedNum]).startsWith(TRAIT_STRAIN_PREFIX)
														|| dnaList.contains(lookup.get(line[pedNum]))	? "2"
																																					: lookup.get(line[pedNum])
																																									.startsWith(TRAIT_STRAIN_PREFIX)	? "1"
																																																										: "0"));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			logIt("Error: file \"" + dir + pedigree + "\" not found in current directory", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		} catch (IOException ioe) {
			logIt("Error reading file \"" + dir + pedigree + "\"", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(2);
		}
	}

	public static Vector<String> getDNAlist(String dir, String genotypes, boolean report) {
		BufferedReader reader;
		String[] line;
		Vector<String> dnaList = new Vector<String>();
		String temp;

		try {
			reader = new BufferedReader(new FileReader(dir + genotypes));
			line = reader.readLine().trim().split("\t");
			for (int i = 3; i < line.length; i++) {
				temp = line[i];
				for (String element : PUNCTUATION) {
					if (line[i].contains(element)) {
						if (!ext.containsAny(temp, IGNORES)) {
							if (dnaList.size() > 2 && report) {
								logIt(line[i] + "; ", false);
							}
							temp = temp.substring(0, temp.indexOf(element));
						}
					}
				}
				if (!ext.containsAny(temp, IGNORES)) {
					dnaList.add(temp);
					if (dnaList.size() > 2) {
						try {
							Integer.parseInt(temp);
						} catch (Exception e) {
							if (report) {
								logIt("\nError - '"	+ line[i]
											+ "' could not be parsed into an integer (closest I got was '" + temp + "')",
											true);
							}
						}
					}
				}
			}
			if (report) {
				logIt("", true);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			logIt("Error: file \"" + dir + genotypes + "\" not found in current directory", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		} catch (IOException ioe) {
			logIt("Error reading file \"" + dir + genotypes + "\"", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(2);
		}

		if (dnaList.size() < 3) {
			logIt("Error - parsed less than three DNAs from the genotype file (n="	+ dnaList.size() + ")",
						true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		}

		return dnaList;
	}

	public static void trimPedigree(String dir, String genotypes, String traitStrainPrefix,
																	String bcStrainPrefix) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Vector<String> v = new Vector<String>();
		Vector<String> preData = new Vector<String>();
		Vector<String> dnaList = new Vector<String>();
		String temp;

		dnaList = getDNAlist(dir, genotypes, false);

		try {
			reader = new BufferedReader(new FileReader(dir + "untrimmed.pre"));
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				if (line[1].startsWith(traitStrainPrefix)) {
					dnaList.add(line[1]);
				}
				preData.add(temp	+ "\t"
										+ (dnaList.contains(line[1]) || line[1].startsWith(traitStrainPrefix)	? "1"
																																													: "0"));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			logIt("Error: file \"" + dir + "untrimmed.pre" + "\" not found in current directory", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		} catch (IOException ioe) {
			logIt("Error reading file \"" + dir + "untrimmed.pre" + "\"", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(2);
		}

		// System.err.println("ERROR! YOU HAVE YET TO CONVERT TO THE NEW INPUT
		// FOR TrimFam!");
		// TrimFam trimmer = new TrimFam(preData, dnaList, true, false);
		TrimFam trimmer = new TrimFam(preData, true, false, false, TrimFam.SCORE_99_NAMING_SCHEME, 0,
																	false, false, new Logger());
		v = trimmer.getExtendedFamilyInformation();
		dnaList.remove(1);
		dnaList.remove(0);
		try {
			writer = new PrintWriter(new FileWriter(dir + "trimmed.pre"));
			for (int i = 0; i < v.size(); i++) {
				line = v.elementAt(i).split("[\\s]+");
				writer.println(line[0]	+ "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4]
												+ "\t" + (line[1].startsWith(bcStrainPrefix) ? "1" : "2"));
				dnaList.remove(line[1]);
			}
			writer.close();
		} catch (Exception e) {
		}

		if (dnaList.size() > 0) {
			logIt("The following dnas were never used: "
						+ Array.toStr(Array.toStringArray(dnaList), ", "), true);

			try {
				writer = new PrintWriter(new FileWriter(dir + "DNAS NEVER USED!!!.txt"));
				for (int i = 0; i < dnaList.size(); i++) {
					writer.println(dnaList.elementAt(i));
				}
				writer.close();
			} catch (IOException ioe) {
				logIt("Error writing to file \"DNAS NEVER USED!!!.txt\" ", true);
			}
		}
	}

	public static void addGenotypes(String dir, String genotypes, double maf,
																	String traitStrainPrefix, String bcStrainPrefix) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, inds = null;
		String temp;
		Vector<String> markers = new Vector<String>();
		int count, index;
		short[][] genos = null;
		boolean done;

		try {
			reader = new BufferedReader(new FileReader(dir + genotypes));
			line = reader.readLine().trim().split("[\\s]+");
			for (int i = 0; i < GENOTYPIC_HEADER.length; i++) {
				if (!line[i].toLowerCase().startsWith(GENOTYPIC_HEADER[i])) {
					logIt("Error - expecting the first three columns in the genotypes file to be "
									+ Array.toStr(GENOTYPIC_HEADER, " ") + " (found " + line[0] + " " + line[1] + " "
								+ line[2] + ")", true);
					logIt("failed at " + ext.getTime(), true);
					System.exit(1);
				}
			}
			done = false;
			while (!done) {
				if (reader.ready()) {
					line = reader.readLine().split("[\\s]+");
					if (line.length > 1) {
						markers.add(line[1] + "\t" + line[0] + "\t0\t" + line[2]);
					}
				} else {
					done = true;
				}
			}
			reader.close();

			reader = new BufferedReader(new FileReader(dir + genotypes));
			inds = Array.subArray(reader.readLine().trim().split("\t", -1), 3);
			for (int i = 0; i < inds.length; i++) {
				logIt((i + 1) + "\t" + inds[i], false);
			}
			genos = new short[inds.length + 1][markers.size()];
			count = 0;
			done = false;
			while (!done) {
				if (reader.ready()) {
					temp = reader.readLine();
					if (temp.trim().split("[\\s]+").length == 1) {
						done = true;
					} else {
						line = temp.trim().split("\t", -1);
						if (line.length != inds.length + 3) {
							logIt("Error - problem with marker " + line[0], true);
						}
						for (int i = 0; i < inds.length; i++) {
							temp = line[3 + i];
							if (temp.equals("AA")) {
								genos[i][count] = 2;
							} else if (temp.equals("AB")) {
								genos[i][count] = 1;
							} else if (temp.equals("BB")) {
								genos[i][count] = 0;
							} else if (temp.equals("oo")	|| temp.equals("0") || temp.equals("")
													|| temp.equals("--") || temp.equals("NoCall") || temp.equals("NoCAAll")
													|| temp.equals("X") || temp.equals("-") || temp.equals("?")) {
								genos[i][count] = -1;
							} else {
								logIt("Error - unknown genotype: " + temp, true);
								logIt("failed at " + ext.getTime(), true);
								System.exit(1);
							}
						}
						count++;
					}
				} else {
					done = true;
				}
			}
			reader.close();
			for (int i = 0; i < markers.size(); i++) {
				genos[inds.length][i] = -1;
			}
			writer = new PrintWriter(new FileWriter(dir + "map.dat"));
			writer.println((markers.size() + 1) + " 0 0 5");
			writer.println("0 0.0 0.0 0");
			writer.println(Array.toStr(Array.stringArraySequence(markers.size() + 1, ""), " "));
			writer.println("1 2 # TRAIT");
			writer.println((1 - DEFAULT_FREQ) + " " + DEFAULT_FREQ);
			writer.println("1");
			writer.println(DEFAULT_MODEL1[0] + " " + DEFAULT_MODEL1[1] + " " + DEFAULT_MODEL1[2]);
			for (int i = 0; i < markers.size(); i++) {
				writer.println("3 2 # " + markers.elementAt(i).split("[\\s]+")[1]);
				writer.println((1 - maf) + " " + maf);
			}
			writer.println("0 0");
			writer.println("10.0 " + Array.toStr(Array.stringArray(markers.size() + 1, "0.0"), " "));
			writer.println("1 0.1 0.45");

			writer.close();

			writer = new PrintWriter(new FileWriter(dir + "markers.dat"));
			for (int i = 0; i < markers.size(); i++) {
				writer.println(markers.elementAt(i));
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			logIt("Error: file \"" + dir + genotypes + "\" not found in current directory", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		} catch (IOException ioe) {
			logIt("Error reading file \"" + dir + genotypes + "\"", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(dir + "trimmed.pre"));
			writer = new PrintWriter(new FileWriter(dir + "pedfile.pre"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				index = ext.indexOfStr(line[1], inds);
				if (line[1].startsWith(traitStrainPrefix)) {
					index = 0;
				}
				if (line[1].startsWith(bcStrainPrefix)) {
					index = 1;
				}
				if (index == -1) {
					index = inds.length;
				}
				// logIt(line[1]+"\t"+index, false);
				writer.println(Array.toStr(line) + "\t" + translateAlleles(genos[index]));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			logIt("Error: file \"" + dir + "trimmed.pre" + "\" not found in current directory", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		} catch (IOException ioe) {
			logIt("Error reading file \"" + dir + "trimmed.pre" + "\"", true);
			logIt("failed at " + ext.getTime(), true);
			System.exit(2);
		}

	}

	public static String translateAlleles(short[] alleles) {
		String str = "";

		for (int i = 0; i < alleles.length; i++) {
			str += (i == 0 ? "" : "\t");
			if (alleles[i] == 2) {
				str += "2\t2";
			} else if (alleles[i] == 1) {
				str += "1\t2";
			} else if (alleles[i] == 0) {
				str += "1\t1";
			} else if (alleles[i] == -1) {
				str += "0\t0";
			} else {
				logIt("Problem with allele calls: '" + alleles[i] + "' is unrecognized", true);
			}
		}

		return str;
	}

	public static void runPedcheck(String dir, String pedfile, String output) throws IOException {
		CmdLine.run("pedcheck -p " + pedfile + " -d map.dat -o " + output + " -4", dir);
		logIt("Pedcheck finished successfully", false);
		new File(dir + "names.tmp").delete();
	}

	public static void mendErrorsFromList(String dir, String fixFile) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, fixes, fix;
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		Vector<String> v;
		int index;
		String[] markerNames;

		try {
			logIt("Mending errors in pedfile.pre using " + fixFile + " and map.dat", false);
			fixes = HashVec.loadFileToStringArray(dir + fixFile, false, new int[] {0, 1, 2, 3, 4}, false);
			markerNames = new LinkageMap(dir + "map.dat").getMarkerNames();

			for (String fixe : fixes) {
				line = fixe.trim().split("[\\s]+");
				index = ext.indexOfStr(line[2], markerNames, false, true);

				if (!hash.containsKey(line[0] + ":" + line[1])) {
					hash.put(line[0] + ":" + line[1], new Vector<String>());
				}
				hash.get(line[0] + ":" + line[1]).add(index + "\t" + line[3] + "\t" + line[4]);
			}

			try {
				reader = new BufferedReader(new FileReader(dir + "pedfile.pre"));
				writer = new PrintWriter(new FileWriter(dir + "mended_pedfile.pre"));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (hash.containsKey(line[0] + ":" + line[1])) {
						v = hash.get(line[0] + ":" + line[1]);
						for (int i = 0; i < v.size(); i++) {
							fix = v.elementAt(i).split("[\\s]+");
							index = Integer.parseInt(fix[0]);
							line[6 + index * 2 + 0] = fix[1];
							line[6 + index * 2 + 1] = fix[2];
						}
					}
					writer.println(Array.toStr(line));
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				logIt("Error: file \"" + dir + "pedfile.pre" + "\" not found in current directory", true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(1);
			} catch (IOException ioe) {
				logIt("Error reading file \"" + dir + "pedfile.pre" + "\"", true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(2);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void mendErrors(String dir) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, keys;
		String temp, fam;
		Vector<String> listOfErrors = new Vector<String>();
		boolean done;
		int locus, numLoci;
		String markerName, father, mother;
		Vector<String[]> errors = new Vector<String[]>();
		Hashtable<String, String> markerMap;
		Hashtable<String, IntVector> hash;
		IntVector iv;
		int level;

		try {
			logIt("Mending errors in pedfile.pre using errors.out and map.dat", false);
			// logIt("pedcheck -p pedfile.pre -d map.dat -o errors.out "+" -4", false);
			temp = "";
			listOfErrors.removeAllElements();

			reader = new BufferedReader(new FileReader(dir + "errors.out"));

			// while (reader.ready() && !reader.readLine().startsWith("
			// *********** LEVEL 1 ERRORS ****************"));
			while (reader.ready() && !reader.readLine().startsWith(" *********** LEVEL ")) {
				;
			}

			if (!reader.ready()) {
				logIt("No errors found in the dataset!", true);
				Files.copyFile(dir + "pedfile.pre", dir + "mended_pedfile.pre");
				reader.close();
				return;
			}

			done = false;
			level = 1;
			while (!done) {
				do {
					temp = reader.readLine();
					if (temp.startsWith(" *********** LEVEL 2 ERRORS ****************")) {
						level = 2;
						// do{
						// reader.readLine();
						// } while (reader.ready());
					}
					if (temp.startsWith(" *********** LEVEL 3 ERRORS ****************")) {
						level = 3;
						do {
							reader.readLine();
						} while (reader.ready());
					}
				} while (reader.ready() && !temp.startsWith("##### GENOTYPE ERROR:"));

				if (!reader.ready()) {
					if (errors.size() == 0) {
						logIt("No errors detected", false);
					}
					done = true;
				} else {
					line = temp.trim().split("[\\s]+");
					fam = line[4];
					locus = Integer.parseInt(line[6]);
					markerName = line[8];

					if (level == 1) {
						while (!reader.readLine().startsWith(" ORIGINAL SCORING:")) {
							;
						}

						line = reader.readLine().trim().split("[\\s]+");
						if (!line[1].endsWith(":") || !line[4].endsWith(":")) {
							logIt("Error - could not parse either father ("	+ line[1] + ") or mother (" + line[4]
										+ ")", true);

						}
						father = line[1].substring(0, line[1].indexOf(":"));
						mother = line[4].substring(0, line[4].indexOf(":"));
					} else { // if (level == 2) {
						while (!reader.readLine().startsWith("ORDERED GENOTYPE LISTS:")) {
							;
						}

						line = reader.readLine().trim().split("[\\s]+");
						father = line[2].substring(0, line[2].indexOf(":"));
						line = reader.readLine().trim().split("[\\s]+");
						mother = line[2].substring(0, line[2].indexOf(":"));
					}
					errors.add(new String[] {fam, locus + "", markerName, father, mother});

				}
			}
			reader.close();

			writer = new PrintWriter(new FileWriter(dir + "logfile of errors.out"));
			writer.println("FamID\tLocus#\tMarker\tFather\tMother");
			markerMap = new Hashtable<String, String>();
			hash = new Hashtable<String, IntVector>();
			for (int i = 0; i < errors.size(); i++) {
				line = errors.elementAt(i);
				writer.println(Array.toStr(line));

				if (!markerMap.containsKey(line[1])) {
					markerMap.put(line[1], line[2]);
				}
				if (!line[2].equals(markerMap.get(line[1]))) {
					logIt("Error - two markers ("	+ line[2] + " and " + markerMap.get(line[1])
								+ ") were mapped to the same locus (" + line[1]
								+ "); this really should be impossible", true);
				}

				if (!hash.containsKey(line[0] + ":" + line[3])) {
					hash.put(line[0] + ":" + line[3], new IntVector());
				}
				hash.get(line[0] + ":" + line[3]).add(Integer.parseInt(line[1]));
				if (!hash.containsKey(line[0] + ":" + line[4])) {
					hash.put(line[0] + ":" + line[4], new IntVector());
				}
				hash.get(line[0] + ":" + line[4]).add(Integer.parseInt(line[1]));
			}
			writer.println();
			writer.println();
			keys = HashVec.getKeys(hash);
			for (String key : keys) {
				writer.print(key);
				for (int j = 0; j < hash.get(key).size(); j++) {
					writer.print("\t" + hash.get(key).elementAt(j));
				}
				writer.println();
			}
			writer.close();

			try {
				reader = new BufferedReader(new FileReader(dir + "map.dat"));
				numLoci = Integer.parseInt(reader.readLine().trim().split("[\\s]+")[0]);
				for (int i = 0; i < 6; i++) {
					reader.readLine();
				}
				for (int i = 1; i <= numLoci; i++) {
					line = reader.readLine().trim().split("[\\s]+");
					if (markerMap.containsKey(i + "") && !markerMap.get(i + "").equals(line[3])) {
						logIt("Error - mismatched marker declarations - genotype files (pedfile.pre and map.dat) probably don't match up with the pedcheck output (errors.out)",
									true);
						logIt("failed at " + ext.getTime(), true);
						System.exit(1);
					}
					reader.readLine();
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				logIt("Error: file \"" + dir + "map.dat" + "\" not found in current directory", true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(1);
			} catch (IOException ioe) {
				logIt("Error reading file \"" + dir + "map.dat" + "\"", true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(2);
			}

			try {
				reader = new BufferedReader(new FileReader(dir + "pedfile.pre"));
				writer = new PrintWriter(new FileWriter(dir + "mended_pedfile.pre"));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					iv = hash.get(line[0] + ":" + line[1]);
					if (iv != null) {
						for (int i = 0; i < iv.size(); i++) {
							line[6 + (iv.elementAt(i) - 1) * 2 + 0] = "0";
							line[6 + (iv.elementAt(i) - 1) * 2 + 1] = "0";
						}
					}
					writer.println(Array.toStr(line));
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				logIt("Error: file \"" + dir + "pedfile.pre" + "\" not found in current directory", true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(1);
			} catch (IOException ioe) {
				logIt("Error reading file \"" + dir + "pedfile.pre" + "\"", true);
				logIt("failed at " + ext.getTime(), true);
				System.exit(2);
			}
		} catch (Exception e) {
			e.printStackTrace();

			logIt("Got an error processing pedcheck results, probably wasn't run. Have the map and pre files been made?",
						true);
			logIt("If that's not it, try deleting the files \'pedcheck.err\' and \'temp\'", true);
			logIt("If that's not it, make sure the first line of the .dat file is either blank or has more than one word.",
						true);

			logIt("", true);
		}
	}

	public static void runMendel(String dir, String filename) {
		PrintWriter writer;
		double[][] models;
		String[][][] allResults;
		SnpMarkerSet markerSet;
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		Vector<String> xMarkers;

		if (!new File(dir + filename).exists()) {
			logIt("Running Mendel using default models; expecting file '"	+ filename
						+ "' for custom models", false);
			models = new double[][] {	{DEFAULT_FREQ, DEFAULT_MODEL1[0], DEFAULT_MODEL1[1],
																DEFAULT_MODEL1[2]},
																{	DEFAULT_FREQ, DEFAULT_MODEL2[0], DEFAULT_MODEL2[1],
																	DEFAULT_MODEL2[2]}};
		} else {
			logIt("Running Mendel using models listed in " + filename, false);
			models = LinkageMap.parseModels(dir + filename);
		}

		markerSet = new SnpMarkerSet(dir	+ "markers.dat", SnpMarkerSet.PLINK_MAP_FORMAT, false,
																	new Logger());
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		xMarkers = new Vector<String>();
		for (int i = 0; i < markerNames.length; i++) {
			if (chrs[i] == 23) {
				xMarkers.add(markerNames[i]);
			}
		}

		allResults = new String[models.length][][];
		for (int i = 0; i < models.length; i++) {
			logIt("Running model " + (i + 1), false);
			allResults[i] = Mendel.runModel(dir, "mended_pedfile.pre", "map.dat", markerSet, models[i][0],
																			models[i][1], models[i][2], models[i][3]);
		}

		try {
			writer = new PrintWriter(new FileWriter(dir + "results.xln"));
			writer.print("Marker\tChromosome\tPosition");
			for (double[] model : models) {
				writer.print("\tModel: dxFreq="	+ model[0] + " penetrance=" + model[1] + "," + model[2]
											+ "," + model[3] + "\tMax Theta");
			}
			writer.println();
			writer.flush();
			for (int j = 0; j < markerNames.length; j++) {
				writer.print(markerNames[j] + "\t" + chrs[j] + "\t" + positions[j]);
				for (int i = 0; i < models.length; i++) {
					writer.print("\t" + allResults[i][j][1] + "\t" + allResults[i][j][2]);
				}
				writer.println();
			}
			writer.close();
		} catch (IOException ioe) {
			logIt("Error writing to results.xln", true);

		}
	}

	public static void exploreModels(String dir, String exploreDir, String exploreFile) {
		PrintWriter writer;
		String[][][][][] allResults;
		String[] markerNames;
		LinkageMap map;
		double[] array;
		int[] keys, indices;
		Hashtable<String, String> chrHash;
		SnpMarkerSet markerSet;
		String root;


		root = ext.rootOf(exploreFile, true);

		markerSet = new SnpMarkerSet(dir	+ "markers.dat", SnpMarkerSet.PLINK_MAP_FORMAT, false,
																	new Logger());
		chrHash = markerSet.getChrHash();

		// markerNames = Matrix.extractColumn(Markers.readIn(dir+"markers.dat"), 0);
		// if (markerNames.length>EXPLORE_LIMIT) {
		// logIt("Found "+markerNames.length+" markers to explore, so we won't explore all possible
		// models; the current limit is set to "+EXPLORE_LIMIT+" if you want to do that...", true);
		// logIt("Successfully completed at "+ext.getTime(), true);
		// System.exit(1);
		// }

		if (!new File(exploreFile).exists()) {
			logIt("Did not find a '"	+ exploreFile
						+ "' file; create this list of markers if you want to explore different models", false);
			return;
		}

		new File(dir + exploreDir).mkdirs();
		LinkageFormat.filterMarkers(dir	+ "mended_pedfile.pre", dir + exploreDir + root + ".pre",
																dir + "map.dat", dir + exploreDir + root + ".dat", exploreFile,
																null);

		dir = dir + exploreDir;

		map = new LinkageMap(dir + root + ".dat");
		markerNames = map.getMarkerNames();

		logIt("Exploring "	+ markerNames.length + " marker" + (markerNames.length > 1 ? "s" : "")
					+ " using 2x" + EXPLORE_PEN.length + "x" + EXPLORE_HET.length + " models using Mendel",
					false);
		allResults = new String[2][EXPLORE_PEN.length][EXPLORE_HET.length][][];

		for (int model = 0; model < 2; model++) {
			for (int i = 0; i < EXPLORE_PEN.length; i++) {
				System.out.print(".");
				for (int j = 0; j < EXPLORE_HET.length; j++) {
					allResults[model][i][j] = Mendel.runModel(dir, root + ".pre", root + ".dat", markerSet,
																										DEFAULT_FREQ, EXPLORE_HET[j],
																										model == 0 ? EXPLORE_PEN[i] : EXPLORE_HET[j],
																										EXPLORE_PEN[i]);
				}
			}
			System.out.println();
		}

		for (int k = 0; k < markerNames.length; k++) {
			try {
				writer = new PrintWriter(new FileWriter(dir + markerNames[k] + ".xln"));
				for (int model = 0; model < 2; model++) {
					writer.println((model == 0 ? "Dominant" : "Recessive") + " model");
					writer.print("Penetrance");
					for (double element : EXPLORE_HET) {
						writer.print("\tHET=" + element + "\tTheta");
					}
					writer.println();
					for (int i = 0; i < EXPLORE_PEN.length; i++) {
						writer.print("PEN=" + EXPLORE_PEN[i]);
						for (int j = 0; j < EXPLORE_HET.length; j++) {
							writer.print("\t"	+ allResults[model][i][j][k][1] + "\t"
														+ allResults[model][i][j][k][2]);
						}
						writer.println();
					}
					writer.println();
				}
				writer.close();
			} catch (IOException ioe) {
				logIt("Error writing to " + dir + markerNames[k] + ".xln", true);
			}

		}

		try {
			writer = new PrintWriter(new FileWriter(dir + "bestModels.xln"));
			for (int model = 0; model < 2; model++) {
				writer.println((model == 0 ? "Dominant" : "Recessive") + " model");
				writer.println("Marker\tChr\tPosition\t1st LOD\t1st model\t2nd LOD\t2nd model\t3rd LOD\t3rd model");
				for (int k = 0; k < markerNames.length; k++) {
					array = new double[EXPLORE_PEN.length * EXPLORE_HET.length];
					for (int i = 0; i < EXPLORE_PEN.length; i++) {
						for (int j = 0; j < EXPLORE_HET.length; j++) {
							array[Matrix.indexInArray(i, j,
																				EXPLORE_HET.length)] = Double.parseDouble(allResults[model][i][j][k][1]);
						}
					}
					keys = Sort.quicksort(array, Sort.DESCENDING);
					writer.print(markerNames[k] + "\t" + chrHash.get(markerNames[k]));
					for (int i = 0; i < 3; i++) {
						indices = Matrix.indicesInMatrix(keys[i], EXPLORE_HET.length);
						writer.print("\t"	+ allResults[model][indices[0]][indices[1]][k][1] + "\tTheta="
													+ allResults[model][indices[0]][indices[1]][k][2] + ", PEN="
													+ EXPLORE_PEN[indices[0]] + ", HET=" + EXPLORE_HET[indices[1]]);
					}
					writer.println();
				}
				writer.println();
			}
			writer.close();
		} catch (IOException ioe) {
			logIt("Error writing to bestModels.xln", true);
		}

	}

	public static boolean finishedExploration(String dir, String exploreDir, String exploreList) {
		SnpMarkerSet markerSet;
		String[] markerNames;

		if (!new File(dir + exploreDir + "bestModels.xln").exists()) {
			return false;
		}

		// markerNames = Matrix.extractColumn(Markers.readIn(dir+exploreList), 0);
		markerSet = new SnpMarkerSet(exploreList, SnpMarkerSet.NAMES_ONLY, false, new Logger());
		markerNames = markerSet.getMarkerNames();
		for (int i = 0; i < markerNames.length; i++) {
			if (!new File(dir + exploreDir + markerNames[i] + ".xln").exists()) {
				return false;
			}
		}
		return true;
	}

	public static void exploreSets(String dir, String listFile) {
		String[] sets;

		sets = HashVec.loadFileToStringArray(dir + listFile, false, new int[] {0}, false);

		for (int i = 0; i < sets.length; i++) {
			if (!finishedExploration(dir, sets[i] + "/", dir + sets[i] + ".txt")) {
				logIt("Exploring " + sets[i], false);
				exploreModels(dir, sets[i] + "/", dir + sets[i] + ".txt");
			}
		}
	}

	public static void logIt(String str, boolean error) {
		PrintWriter writer;

		if (error) {
			System.err.println(str);
		} else {
			System.out.println(str);
		}
		try {
			writer = new PrintWriter(new FileWriter(LOG_FILE, true));
			writer.println(str);
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + LOG_FILE);
			e.printStackTrace();
		}
	}

	public static void yossiPedigreeMaker(String phenoFile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		int[] indices;
		String cagePrefix;

		String[] headerTargets = new String[] {"IID", "Cage", "SEX", "AffectionStatus"};
		String[][] genderCodes = new String[][] {{"M", "1"}, {"F", "2"}};

		try {
			reader = new BufferedReader(new FileReader(phenoFile));
			indices = ext.indexFactors(	headerTargets, reader.readLine().trim().split("[\\s]+"), false,
																	true);
			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(phenoFile) + "trimmed.pre"));
			writer.println("1\t99901\t0\t0\t1\t2");
			writer.println("1\t77701\t0\t0\t2\t1");
			writer.println("1\t88801\t99901\t77701\t1\t2");
			writer.println("1\t77702\t0\t0\t2\t1");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				cagePrefix = (66000 + Integer.parseInt(line[indices[1]])) + "";
				if (!hash.containsKey(cagePrefix)) {
					writer.println("1\t" + cagePrefix + "1\t88801\t77702\t1\t0");
					writer.println("1\t" + cagePrefix + "2\t88801\t77702\t2\t0");
					hash.put(cagePrefix, "");
				}
				writer.println("1\t"	+ line[indices[0]] + "\t" + cagePrefix + "1\t" + cagePrefix + "2\t"
												+ ext.replaceAllWith(line[indices[2]], genderCodes) + "\t"
												+ line[indices[3]]);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + phenoFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + phenoFile + "\"");
			System.exit(2);
		}

		Files.copyFile(ext.parseDirectoryOfFile(phenoFile)	+ "trimmed.pre",
										ext.parseDirectoryOfFile(phenoFile) + "untrimmed.pre");

	}

	public static void main(String[] args) throws IOException {
		String dir;
		int numArgs = args.length;
		boolean jar = true;

		String pedigree = "pedigree.txt";
		String genotypes = "genotypes.txt";
		// String output = "trimmed.pre";
		// String pedfile = "pedfile.pre";
		// String mapfile = "map.dat";
		// String markersfile = "markers.dat";
		String models = "models.dat";
		String exploreDir = DEFAULT_EXPLORE_DIR;
		String exploreFile = DEFAULT_EXPLORE_FILE;
		double maf = 0.5;
		int option = -1;

		// yossiPedigreeMaker("C:\\Documents and Settings\\npankrat\\My
		// Documents\\murine\\yossi\\2010.02.03\\phenotype.dat");
		// System.exit(1);

		String usage = "\\n"	+ "oneoff.GinsburgPedigree requires 0-3 arguments\n"
										+ "   (1) pedigree file (i.e. ped=" + pedigree + " (default))\n"
										+ "   (2) genotype file (i.e. geno=" + genotypes + " (default))\n"
										+ "   (3) minor allele frequency (i.e. maf=" + maf + " (default))\n"
										+ "   (4) (optional) list of models (i.e. models=" + models + ")\n" + "\n"
										+ " Note: you can pick one of the following steps to perform; otherwise the step to be performed at each iteration will be determined by which files are already present:\n"
										+ "   -makePedigree\n" + "   -trimPedigree\n" + "   -addGenotypes\n"
										+ "   -mendErrors\n" + "   -runMendel\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				logIt("failed at " + ext.getTime(), true);
				System.exit(1);
			} else if (arg.startsWith("ped=")) {
				pedigree = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("geno=")) {
				genotypes = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("models=")) {
				models = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("maf=")) {
				maf = Double.parseDouble(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-makePedigree")) {
				option = 0;
				numArgs--;
			} else if (arg.startsWith("-trimPedigree")) {
				option = 1;
				numArgs--;
			} else if (arg.startsWith("-addGenotypes")) {
				option = 2;
				numArgs--;
			} else if (arg.startsWith("-mendErrors")) {
				option = 3;
				numArgs--;
			} else if (arg.startsWith("-runMendel")) {
				option = 4;
				numArgs--;
			} else if (arg.startsWith("-notJar")) {
				jar = false;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			logIt("failed at " + ext.getTime(), true);
			System.exit(1);
		}
		try {
			dir = jar	? "./"
								: (new File(DEFAULT_ROOT + DEFAULT_DIR).exists() ? DEFAULT_ROOT + DEFAULT_DIR : "");
			if (option == 0) {
				makePedigree(dir, pedigree, genotypes);
			} else if (option == 1) {
				trimPedigree(dir, genotypes, TRAIT_STRAIN_PREFIX, BACKCROSS_STRAIN_PREFIX);
			} else if (option == 2) {
				addGenotypes(dir, genotypes, maf, TRAIT_STRAIN_PREFIX, BACKCROSS_STRAIN_PREFIX);
			} else if (option == 3) {
				mendErrors(dir);
			} else if (option == 4) {
				runMendel(dir, models);
			} else if (option == -1) {
				logIt("Started new run on " + ext.getDate() + " at " + ext.getTime(), false);

				logIt("Determination of how to proceed next is based on files present in "
							+ (dir.equals("") ? "current directory" : dir), false);
				if (new File(dir + "untrimmed.pre").exists()) {
					logIt("Found untrimmed.pre, skipping makePedigree", false);
				} else {
					logIt("Making pedigree...", false);
					makePedigree(dir, pedigree, genotypes);
					logIt("", false);
				}
				if (new File(dir + "trimmed.pre").exists()) {
					logIt("Found trimmed.pre, skipping trimPedigree", false);
				} else {
					logIt("Trimming pedigree...", false);
					trimPedigree(dir, genotypes, TRAIT_STRAIN_PREFIX, BACKCROSS_STRAIN_PREFIX);
					logIt("", false);
				}
				if (new File(dir + "pedfile.pre").exists()) {
					logIt("Found pedfile.pre, skipping addGenotypes", false);
				} else {
					logIt("Adding Genotypes...", false);
					addGenotypes(dir, genotypes, maf, TRAIT_STRAIN_PREFIX, BACKCROSS_STRAIN_PREFIX);
					logIt("", false);
				}
				if (new File(dir + "errors.out").exists()) {
					logIt("Found errors.out, assuming pedcheck ran successfully the first time", false);
				} else {
					logIt("Running pedcheck -p pedfile.pre -d map.dat -o errors.out -4", false);
					runPedcheck(dir, "pedfile.pre", "errors.out");
					logIt("", false);
				}
				if (!new File(dir + "mended_pedfile.pre").exists()) {
					if (new File(dir + "fixes.xln").exists()) {
						mendErrorsFromList(dir, "fixes.xln");
					} else {
						mendErrors(dir);
					}
					logIt("", false);
				}
				if (!new File(dir + "doublecheck.out").exists()) {
					logIt("Double checking that all errors were caught", false);
					runPedcheck(dir, "mended_pedfile.pre", "doublecheck.out");
					logIt("", false);
				}
				if (!new File(dir + "results.xln").exists()) {
					runMendel(dir, models);
					logIt("", false);
				}
				if (new File(dir + "sets.txt").exists()) {
					exploreSets(dir, "sets.txt");
				} else if (!finishedExploration(dir, exploreDir, dir + exploreFile)) {
					exploreModels(dir, exploreDir, exploreFile);
					logIt("", false);
				} else {
					logIt("Looks like you're at the end of the line...", false);
				}
				logIt("Finished successfully at " + ext.getTime(), true);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
