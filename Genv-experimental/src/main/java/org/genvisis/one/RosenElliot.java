package org.genvisis.one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class RosenElliot {
	public static final String DEFAULT_ROOT =
																					"C:\\Documents and Settings\\npankrat\\My Documents\\Rosen, Elliot\\";
	// public static final String DEFAULT_DIR = "Data\\";
	// public static final String DEFAULT_DIR = "Data\\exploreModels\\";
	// public static final String DEFAULT_DIR = "Data\\longTerm\\";
	// public static final String DEFAULT_DIR = "Data\\wo6thlineage\\";
	public static final String DEFAULT_DIR = "Data\\swap6thlineage\\";
	public static final String PED_DATA = "pedData.xln";
	public static final String GENOTYPES_FILE = "genotypes.txt";
	public static final String[] PED_HEADER = {	"Animal", "Sire", "DOB", "DOD", "DaysAtDeath",
																							"DaysLivedSoFar", "Affection"};
	public static final String TRAIT_STRAIN = "129";
	public static final String BC_STRAIN = "C57";
	public static final String TRAIT_STRAIN_PREFIX = "99";
	public static final String BC_STRAIN_PREFIX = "77";
	// public static final String[] X_MARKERS_TO_KEEP = {"rs13484113"};
	public static final String[] PREFIXES = {"AS"};
	public static final String[] SUFFIXES = {"f"};
	public static final String[][] DIP_TRIP = {{"AA", "BB"}, {"AB", "AB"}, {"BB", "AA"}, {"0", "0"}};
	public static final int EXP_DNA_LEN = 4;
	public static final int AFF_THRESHOLD = 30;

	public static void decipherPedigree(String dir, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, subline;
		String temp, trav, dna;
		String[] tree;
		String sire, dob, dod;
		int daysAtDeath, daysLivedSoFar, aff;
		Date latest = new Date();
		Hashtable<String, String> hash = new Hashtable<String, String>();

		try {
			reader = new BufferedReader(new FileReader(dir + filename));
			writer = new PrintWriter(new FileWriter(dir + PED_DATA));
			writer.println("Animal\tSire\tDOB\tDOD\tDaysAtDeath\tDaysLivedSoFar\tAffection");
			if (!reader.readLine().contains("List")) {
				System.err.println("Error parsing pedigree file: expecting first line to contain the word List");
				System.exit(1);
			}
			line = reader.readLine().trim().split("[\\s]+");
			if (line.length != 1) {
				System.err.println("Error parsing pedigree file: expecting 2nd line to have the latest date in it");
				System.exit(1);
			}
			try {
				latest = ext.parseDate(line[0]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error parsing date: " + line[0]);
			}
			if (!reader.readLine().contains("DOB")) {
				System.err.println("Error parsing pedigree file: expecting 3rd line to have DOB in it");
				System.exit(1);
			}
			if (!reader.readLine().contains("(DOD)")) {
				System.err.println("Error parsing pedigree file: expecting 4th line to have (DOD) in it");
				System.exit(1);
			}
			tree = new String[10];
			trav = sire = dob = dod = "err";
			daysAtDeath = daysLivedSoFar = aff = -1;
			while (reader.ready()) {
				temp = reader.readLine();
				if (!temp.trim().equals("")) {
					line = temp.split("\\t", -1);
					for (int i = 0; i < line.length / 2; i++) {
						if (!line[i * 2].trim().equals("")) {
							trav = line[i * 2];
							dna = cleanUpID(trav);
							if (dna.startsWith(TRAIT_STRAIN_PREFIX) || dna.startsWith(BC_STRAIN_PREFIX)) {
								System.err.println("Error - sample '"	+ trav + "' has a dna '" + dna
																		+ "' that starts with a reserved prefix");
							}
							if (hash.containsKey(trav)) {
								System.err.println("Error - two samples with the same ID ('" + trav + "')");
							} else if (hash.containsKey(dna)) {
								System.err.println("Error - two samples with the same parsed DNA number ('"	+ dna
																		+ "')");
							}
							hash.put(trav, "");
							hash.put(dna, "");
							dob = line[i * 2 + 1];
							reader.mark(1000);
							subline = reader.readLine().split("\\t", -1);
							if (subline[i * 2 + 1].startsWith("(")) {
								dod = subline[i * 2 + 1].substring(1, subline[i * 2 + 1].length() - 1);
								daysAtDeath = ext.calcDays(ext.parseDate(dob), ext.parseDate(dod));
								daysLivedSoFar = -1;
								if (daysAtDeath >= AFF_THRESHOLD) {
									aff = 2;
								} else {
									aff = 1;
								}
							} else {
								dod = ".";
								daysLivedSoFar = ext.calcDays(ext.parseDate(dob), latest);
								daysAtDeath = -1;
								if (daysLivedSoFar >= AFF_THRESHOLD) {
									aff = 2;
								} else {
									aff = 0;
								}
								reader.reset();
							}
							tree[i] = trav;
							for (int j = i + 1; j < tree.length; j++) {
								tree[j] = "err";
							}
							if (i == 0) {
								sire = "Inbred";
							} else {
								sire = tree[i - 1];
							}
						}
					}
					writer.println(trav	+ "\t" + sire + "\t" + dob + "\t" + dod + "\t"
													+ (daysAtDeath == -1 ? "." : daysAtDeath) + "\t"
													+ (daysLivedSoFar == -1 ? "." : daysLivedSoFar) + "\t" + aff);
					daysAtDeath = daysLivedSoFar = aff = -1;
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}
	}

	public static void reformatGenotypes(String dir, String filename, String refGenotypes) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Vector<String> dnas = new Vector<String>();
		int index;
		boolean reversi;

		try {
			reader = new BufferedReader(new FileReader(dir + filename));
			writer = new PrintWriter(new FileWriter(dir + GENOTYPES_FILE));

			line = reader.readLine().trim().split("\\t");
			index = ext.indexOfStr(refGenotypes, line);
			if (index == -1) {
				System.err.println("Error - the reference dna '"	+ refGenotypes
														+ "' was not found in the header");
				System.exit(1);
			}

			writer.print(line[0]	+ "\t" + line[1] + "\t" + line[2] + "\t" + TRAIT_STRAIN + "\t"
										+ BC_STRAIN);
			for (int i = 3; i < line.length; i++) {
				if (i != index) {
					trav = cleanUpDNA(line[i]);
					writer.print("\t" + trav);
					if (dnas.indexOf(trav) != -1) {
						System.err.println("Warning - '" + trav + "' is duplicated, please collapse");
					}
					dnas.add(trav);
				}
			}
			writer.println();

			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				reversi = line[index].equals("AA");

				writer.print(line[0]	+ "\t" + line[1] + "\t" + line[2] + "\t"
											+ (reversi	? line[index] + "\t" + flip(line[index])
																	: flip(line[index]) + "\t" + line[index]));
				for (int i = 3; i < line.length; i++) {
					if (line[i].equals("H")) {
						line[i] = "AB";
					}
					if (i != index) {
						writer.print("\t" + (reversi ? flip(line[i]) : line[i]));
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}

	}

	public static String cleanUpID(String str) {
		str = str.trim();
		for (String element : PREFIXES) {
			if (str.startsWith(element)) {
				str = str.substring(element.length());
			}
		}
		for (String element : SUFFIXES) {
			if (str.endsWith(element)) {
				str = str.substring(0, str.length() - element.length());
			}
		}
		if (str.length() != EXP_DNA_LEN) {
			System.err.println("Error - parsed DNA# '"	+ str + "' was not of expected length: "
													+ EXP_DNA_LEN);
		}

		return str;
	}

	public static String cleanUpDNA(String str) {
		str = str.trim().split("[\\s]+")[0];
		if (str.startsWith("*")) {
			str = str.substring(1);
		}
		if (str.length() != EXP_DNA_LEN) {
			System.err.println("Error - parsed DNA# '"	+ str + "' was not of expected length: "
													+ EXP_DNA_LEN);
		}

		return str;
	}

	public static String flip(String str) {
		for (String[] element : DIP_TRIP) {
			if (str.equals(element[0])) {
				return element[1];
			}
		}

		System.err.println("Error - could not find the flip of '" + str + "'");
		return null;
	}

	public static void makePedigree(String dir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String dna, sire, dam;
		int count;

		Hashtable<String, String> lookup = new Hashtable<String, String>();
		Vector<String> dnaList = new Vector<String>();

		try {
			dnaList = GinsburgDavid.getDNAlist(dir, GENOTYPES_FILE, true);
			System.out.println("Assuming we're taking the F1s of a "	+ dnaList.elementAt(0) + " x "
													+ dnaList.elementAt(1) + " cross and backcrossing into "
													+ dnaList.elementAt(1));

			reader = new BufferedReader(new FileReader(dir + PED_DATA));
			writer = new PrintWriter(new FileWriter(dir + "untrimmed.pre"));
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), PED_HEADER, true);

			writer.println("1\t" + TRAIT_STRAIN_PREFIX + "99\t0\t0\t1\t2");
			writer.println("1\t" + BC_STRAIN_PREFIX + "01\t0\t0\t2\t1");
			count = 1;

			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				dna = cleanUpID(line[0]);
				if (line[1].equals("Inbred")) {
					writer.println("1\t"	+ dna + "\t" + TRAIT_STRAIN_PREFIX + "99\t" + BC_STRAIN_PREFIX
													+ "01\t1\t" + (dnaList.contains(dna) ? "2" : "0"));
				} else {
					sire = cleanUpID(line[1]);
					if (lookup.containsKey(sire)) {
						dam = lookup.get(sire);
					} else {
						lookup.put(sire, dam = BC_STRAIN_PREFIX + ext.formNum(++count, 2));
						writer.println("1\t" + dam + "\t0\t0\t2\t1");
					}
					writer.println("1\t"	+ dna + "\t" + sire + "\t" + dam + "\t1\t"
													+ (dnaList.contains(dna) ? "2" : "0"));
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + PED_DATA + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + PED_DATA + "\"");
			System.exit(2);
		}
	}

	public static void updateTrait(String dir, String traitfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav, filename;
		Hashtable<String, String> hash = new Hashtable<String, String>();

		try {
			reader = new BufferedReader(new FileReader(dir + traitfile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				trav = cleanUpID(line[0]);
				if (hash.containsKey(trav)) {
					System.err.println("Error - multiple entries for '" + trav + "'");
				}
				if (!line[1].equals("2") && !line[1].equals("1") && !line[1].equals("0")) {
					System.err.println("Error - the value for '"	+ trav + "' ('" + line[1]
															+ "') is not valid; should be either 0, 1, or 2");
				}
				hash.put(trav, line[1]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + traitfile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + traitfile + "\"");
			System.exit(2);
		}

		filename = Files.backup("mended_pedfile.pre", dir, dir);
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(dir + "mended_pedfile.pre"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[1].startsWith(TRAIT_STRAIN_PREFIX) || line[1].startsWith(BC_STRAIN_PREFIX)) {
				} else if (hash.containsKey(line[1])) {
					line[5] = hash.get(line[1]);
				} else {
					line[5] = "0";
				}
				writer.println(ArrayUtils.toStr(line));
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		String dir = new File(DEFAULT_ROOT + DEFAULT_DIR).exists() ? DEFAULT_ROOT + DEFAULT_DIR : "";
		int numArgs = args.length;
		String pedigreeData = "pedigreeData.txt";
		int option = -1;
		// int option = 0;
		String genotypicData = "genotypicData.txt";
		String models = "models.dat";
		String refGenotypes = "B6_1533360010_R004_C002";
		double maf = 0.5;
		String traitfile = "trait.dat";

		String usage = "\\n"	+ "park.oneoff.RosenElliot requires 0-1 arguments\n"
										+ "   (1) pedigree data (i.e. ped=" + pedigreeData + " (default))\n"
										+ "   (2) genotypic data (i.e. geno=" + genotypicData + " (default))\n"
										+ "   (3) minor allele frequency (i.e. maf=" + maf + " (default))\n"
										+ "   (4) (optional) list of models (i.e. models=" + models + ")\n"
										+ "   (5) (optional) new trait file (i.e. trait=" + traitfile + ")\n" + "\n"
										+ " Note: you can pick one of the following steps to perform; otherwise the step to be performed at each iteration will be determined by which files are already present:\n"
										+ "   -parsePedigree\n" + "   -trimPedigree\n" + "   -addGenotypes\n"
										+ "   -mendErrors\n" + "   -runMendel\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("ped=")) {
				pedigreeData = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("geno=")) {
				genotypicData = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("maf=")) {
				maf = Double.parseDouble(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("models=")) {
				models = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("trait=")) {
				traitfile = arg.split("=")[1];
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
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (option == 0) {
				decipherPedigree(dir, pedigreeData);
				reformatGenotypes(dir, genotypicData, refGenotypes);
				makePedigree(dir);
			} else if (option == 1) {
				GinsburgDavid.trimPedigree(dir, GENOTYPES_FILE, TRAIT_STRAIN_PREFIX, BC_STRAIN_PREFIX);
			} else if (option == 2) {
				GinsburgDavid.addGenotypes(dir, GENOTYPES_FILE, maf, TRAIT_STRAIN_PREFIX, BC_STRAIN_PREFIX);
			} else if (option == 3) {
				GinsburgDavid.mendErrors(dir);
			} else if (option == 4) {
				GinsburgDavid.runMendel(dir, models);
			} else if (option == -1) {
				System.out.println("Determination of how to proceed next is based on files present in "
														+ (dir.equals("") ? "current directory" : dir));
				if (!new File(dir + "pedData.xln").exists()) {
					System.out.println("Deciphering pedigree data...");
					decipherPedigree(dir, pedigreeData);
					System.out.println();
				}
				if (!new File(dir + "genotypes.txt").exists()) {
					System.out.println("Reformating genotypes...");
					reformatGenotypes(dir, genotypicData, refGenotypes);
					System.out.println();
				}
				if (!new File(dir + "untrimmed.pre").exists()) {
					System.out.println("Making pedigree...");
					makePedigree(dir);
					System.out.println();
				}
				if (!new File(dir + "trimmed.pre").exists()) {
					System.out.println("Trimming pedigree...");
					GinsburgDavid.trimPedigree(dir, GENOTYPES_FILE, TRAIT_STRAIN_PREFIX, BC_STRAIN_PREFIX);
					System.out.println();
				}
				if (!new File(dir + "pedfile.pre").exists()) {
					System.out.println("Adding Genotypes...");
					GinsburgDavid.addGenotypes(	dir, GENOTYPES_FILE, maf, TRAIT_STRAIN_PREFIX,
																			BC_STRAIN_PREFIX);
					System.out.println();
				}
				if (!new File(dir + "errors.out").exists()) {
					System.out.println("Running pedcheck -p pedfile.pre -d map.dat -o errors.out -4");
					GinsburgDavid.runPedcheck(dir, "pedfile.pre", "errors.out");
					System.out.println();
				}
				if (!new File(dir + "mended_pedfile.pre").exists()) {
					GinsburgDavid.mendErrors(dir);
					System.out.println();
				}
				if (!new File(dir + "doublecheck.out").exists()) {
					System.out.println("Double checking that all errors were caught");
					GinsburgDavid.runPedcheck(dir, "mended_pedfile.pre", "doublecheck.out");
					System.out.println();
				}
				if (new File(dir + traitfile).exists()) {
					System.out.println("Updating affection status using the phenotype in '"	+ traitfile
															+ "'");
					updateTrait(dir, traitfile);
					System.out.println();
				}
				if (!new File(dir + "results.xln").exists()) {
					GinsburgDavid.runMendel(dir, models);
					System.out.println();
				}
				if (!GinsburgDavid.finishedExploration(	dir, GinsburgDavid.DEFAULT_EXPLORE_DIR,
																								GinsburgDavid.DEFAULT_EXPLORE_FILE)) {
					GinsburgDavid.exploreModels(dir, GinsburgDavid.DEFAULT_EXPLORE_DIR,
																			GinsburgDavid.DEFAULT_EXPLORE_FILE);
					System.out.println();
				} else {
					System.out.println("Looks like you're at the end of the line...");
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
