package org.genvisis.link;

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
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class AltPheno {
	public static final String[] MERLIN_BINARY_HEADER = {	"CHR", "POS", "LABEL", "ANALYSIS", "ZSCORE",
																												"DELTA", "LOD", "PVALUE"};

	public static void rerunAnalyses(String dir, String filename) {
		BufferedReader reader;
		PrintWriter[] writers;
		String[] line;
		String trav;
		Hashtable<String, Hashtable<String, String>> hash =
																											new Hashtable<String, Hashtable<String, String>>();
		Vector<String> v;
		String[] chrs, phenoNames = null, fams, phenos;
		String classpath;

		v = new Vector<String>();
		for (int i = 1; i <= 23; i++) {
			if (new File(dir + "re_chrom" + ext.chrome(i) + ".pre").exists()) {
				if (new File(dir + "map" + ext.chrome(i) + ".dat").exists()) {
					v.add(ext.chrome(i));
				} else {
					System.err.println("Error - found re_chrom"	+ ext.chrome(i) + ".pre but not map"
															+ ext.chrome(i) + ".dat");
				}
			}
		}
		chrs = Array.toStringArray(v);
		System.out.println("Found data for chromosomes " + ext.listRanges(Array.toIntArray(chrs)));

		try {
			reader = new BufferedReader(new FileReader(dir + filename));
			line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equals("FID") || !line[1].equals("IID")) {
				System.err.println("Warning - assuming first two columns of "	+ filename
														+ " are FID and IID (found " + line[0] + " and " + line[1] + ")");
			}
			phenoNames = new String[line.length - 2];
			for (int i = 0; i < phenoNames.length; i++) {
				phenoNames[i] = line[i + 2];
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				HashVec.addToHashHash(hash, line[0], line[1], Array.toStr(Array.subArray(line, 2)));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}

		fams = HashVec.getKeys(hash);
		for (int i = 0; i < fams.length + 1; i++) {
			if (i == fams.length) {
				trav = dir + "allFams/";
			} else {
				trav = dir + fams[i] + "/";
			}
			new File(trav).mkdirs();

			for (String phenoName : phenoNames) {
				new File(trav + phenoName + "/").mkdirs();
			}
			for (String chr2 : chrs) {
				try {
					reader = new BufferedReader(new FileReader(dir + "re_chrom" + chr2 + ".pre"));
					writers = new PrintWriter[phenoNames.length];
					for (int j = 0; j < phenoNames.length; j++) {
						writers[j] = new PrintWriter(new FileWriter(trav	+ phenoNames[j] + "/re_chrom" + chr2
																												+ ".pre"));
						Files.copyFile(dir	+ "map" + chr2 + ".dat",
														trav + phenoNames[j] + "/map" + chr2 + ".dat");
					}
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (i == fams.length || line[0].equals(fams[i])) {
							phenos = hash	.get(line[0])
														.containsKey(line[1])	? hash.get(line[0]).get(line[1]).split("[\\s]+")
																									: Array.stringArray(phenoNames.length, "0");
							for (int j = 0; j < phenoNames.length; j++) {
								line[5] = phenos[j].equals("1") ? "2" : "0";
								writers[j].println(Array.toStr(line));
							}
						}
					}
					reader.close();
					for (int j = 0; j < phenoNames.length; j++) {
						writers[j].close();
					}
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ dir + "re_chrom" + chr2 + ".pre"
															+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir + "re_chrom" + chr2 + ".pre" + "\"");
					System.exit(2);
				}
			}
		}

		classpath = "-classpath /home/npankrat/" + org.genvisis.common.PSF.Java.GENVISIS + "";
		if (new File(".").getAbsolutePath().contains("bc2/pankratz")) {
			classpath = "-classpath /home/bc2/pankratz/" + org.genvisis.common.PSF.Java.GENVISIS + "";
		}
		try {
			writers = new PrintWriter[phenoNames.length];
			for (int j = 0; j < phenoNames.length; j++) {
				writers[j] = new PrintWriter(new FileWriter(dir + phenoNames[j] + ".batch"));
			}
			for (int i = 0; i < fams.length + 1; i++) {
				ext.writeToAll("cd " + (i == fams.length ? "allFams/" : fams[i] + "/"), writers);
				for (int j = 0; j < phenoNames.length; j++) {
					writers[j].println("cd " + phenoNames[j] + "/");
				}
				for (String chr2 : chrs) {
					ext.writeToAll("java "	+ classpath + " link.Merlin chr=" + Integer.parseInt(chr2),
													writers);
					ext.writeToAll((chr2.equals("23") ? "minx" : "merlin")	+ " -d chr" + chr2
													+ ".dat -p re_chrom" + chr2 + ".pre -m chr" + chr2 + ".map -f chr" + chr2
													+ ".freq --npl --tabulate --step 5 --markerNames --information --prefix merlin-chr"
													+ chr2 + "", writers);
				}
				ext.writeToAll("cd ..", writers);
				ext.writeToAll("cd ..", writers);
				ext.writeToAll("", writers);
			}

			for (int j = 0; j < phenoNames.length; j++) {
				writers[j].close();
				Files.chmod(dir + phenoNames[j] + ".batch");
			}
		} catch (Exception e) {
			System.err.println("Error writing to batch files");
			e.printStackTrace();
		}
	}

	public static void parse(String dir, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, percentages;
		Hashtable<String, Hashtable<String, Vector<String>>> hashes =
																																new Hashtable<String, Hashtable<String, Vector<String>>>();
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		Vector<String> v;
		String[] phenoNames = null, fams;
		double trav, max, maxPossible;

		try {
			reader = new BufferedReader(new FileReader(dir + filename));
			line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equals("FID") || !line[1].equals("IID")) {
				System.err.println("Warning - assuming first two columns of "	+ filename
														+ " are FID and IID (found " + line[0] + " and " + line[1] + ")");
			}
			phenoNames = new String[line.length - 2];
			for (int i = 0; i < phenoNames.length; i++) {
				phenoNames[i] = line[i + 2];
			}
			hashes.put("allFams", hash = new Hashtable<String, Vector<String>>());
			for (String phenoName : phenoNames) {
				hash.put(phenoName, new Vector<String>());
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (hashes.containsKey(line[0])) {
					hash = hashes.get(line[0]);
				} else {
					hashes.put(line[0], hash = new Hashtable<String, Vector<String>>());
					for (String phenoName : phenoNames) {
						hash.put(phenoName, new Vector<String>());
					}
				}
				for (int i = 0; i < phenoNames.length; i++) {
					if (line[2 + i].equals("1")) {
						hash.get(phenoNames[i]).add(line[1]);
						hashes.get("allFams").get(phenoNames[i]).add(line[1]);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}

		try {
			writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename) + "_summary.xln"));
			writer.println("FID\tPhenotype\tN affected\tAffecteds\t"
											+ Array.toStr(Array.stringArraySequence(23, "")));
			fams = HashVec.getKeys(hashes);
			for (String fam : fams) {
				writer.print(fam);
				for (String phenoName : phenoNames) {
					v = hashes.get(fam).get(phenoName);
					writer.print("\t"	+ phenoName + "\t" + v.size() + "\t"
												+ (fam.equals("allFams")	? "ALL"
																									: Array.toStr(Array.toStringArray(v), ",")));
					percentages = new String[23];
					for (int chr = 1; chr <= 23; chr++) {
						try {
							reader = new BufferedReader(new FileReader(dir	+ fam + "/" + phenoName + "/merlin-chr"
																													+ ext.chrome(chr)
																													+ "-nonparametric.tbl"));
							ext.checkHeader(reader.readLine().trim().split("[\\s]+"), MERLIN_BINARY_HEADER, true);
							if (reader.ready()) {
								reader.readLine();
								maxPossible = Double.parseDouble(reader.readLine().trim().split("\\t", -1)[6]);
								max = Double.NEGATIVE_INFINITY;
								while (reader.ready()) {
									trav = Double.parseDouble(reader.readLine().trim().split("\\t", -1)[6]);
									if (trav > max) {
										max = trav;
									}
								}
								writer.print("\t" + max);
								percentages[chr - 1] = ext.formDeci(max / maxPossible * 100, 1) + "%";
							} else {
								writer.print("\tx");
								percentages[chr - 1] = "x";
							}

							reader.close();
						} catch (FileNotFoundException fnfe) {
							System.err.println("Error: file \"merlin-chr"	+ ext.chrome(chr)
																	+ "-nonparametric.tbl\" not found for " + fam + phenoName + "/");
							writer.print("\t.");
						} catch (IOException ioe) {
							System.err.println("Error reading file \"merlin-chr"	+ ext.chrome(chr)
																	+ "-nonparametric.tbl\" not found for " + fam + phenoName + "/");
							writer.print("\t???");
						} catch (NullPointerException npe) {
							System.err.println("Error haven't finished writing file \"merlin-chr"
																		+ ext.chrome(chr) + "-nonparametric.tbl\" not found for " + fam
																	+ phenoName + "/");
							writer.print("\t?!?");
						}
					}
					writer.println();
					writer.println("\t\t\t\t" + Array.toStr(percentages));
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir + ext.rootOf(filename) + "_summary.xln");
			e.printStackTrace();
		}
	}

	public static void zip(String dir, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, Hashtable<String, Vector<String>>> hashes =
																																new Hashtable<String, Hashtable<String, Vector<String>>>();
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		String[] phenoNames = null, fams;

		try {
			reader = new BufferedReader(new FileReader(dir + filename));
			line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equals("FID") || !line[1].equals("IID")) {
				System.err.println("Warning - assuming first two columns of "	+ filename
														+ " are FID and IID (found " + line[0] + " and " + line[1] + ")");
			}
			phenoNames = new String[line.length - 2];
			for (int i = 0; i < phenoNames.length; i++) {
				phenoNames[i] = line[i + 2];
			}
			hashes.put("allFams", hash = new Hashtable<String, Vector<String>>());
			for (String phenoName : phenoNames) {
				hash.put(phenoName, new Vector<String>());
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (hashes.containsKey(line[0])) {
					hash = hashes.get(line[0]);
				} else {
					hashes.put(line[0], hash = new Hashtable<String, Vector<String>>());
					for (String phenoName : phenoNames) {
						hash.put(phenoName, new Vector<String>());
					}
				}
				for (int i = 0; i < phenoNames.length; i++) {
					if (line[2 + i].equals("1")) {
						hash.get(phenoNames[i]).add(line[1]);
						hashes.get("allFams").get(phenoNames[i]).add(line[1]);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}
		try {
			writer = new PrintWriter(new FileWriter(dir + "batchZip"));
			fams = HashVec.getKeys(hashes);
			for (String fam : fams) {
				writer.println("tar cvf - " + fam + "/*/merlin-* | gzip > " + fam + ".tar.gz");
			}
			writer.close();
			Files.chmod(dir + "batchZip");
		} catch (Exception e) {
			System.err.println("Error writing to " + dir + "batchZip");
			e.printStackTrace();
		}
	}

	public static void condense(String dir, String trait) {
		String[] dirs, files;

		dirs = Files.listDirectories(dir, false);
		new File(dir + "_" + trait).mkdirs();
		for (int i = 0; i < dirs.length; i++) {
			if (!dirs[i].startsWith("_")) {
				if (new File(dir + dirs[i] + "/" + trait).exists()) {
					System.out.println("Copying files from " + dirs[i] + "/");
					files = Files.list(dir + dirs[i] + "/" + trait, ".tbl", false);
					new File(dir + "_" + trait + "/" + dirs[i]).mkdirs();
					for (String file : files) {
						Files.copyFile(dir	+ dirs[i] + "/" + trait + "/" + file,
														dir + "_" + trait + "/" + dirs[i] + "/" + file);
					}
				} else {
					System.err.println("Error - directory '"	+ dirs[i]
															+ "/' did not contain a subdirectory called '" + trait + "/'");
				}
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\6K
		// Screens\\07_audit\\testBatch\\";
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\6K Screens\\08_AltPheno\\";
		// String dir = "";
		String filename = "AltPheno_6K.dat";
		boolean parse = false;
		boolean zip = false;
		String condense = "PDish";

		String usage = "\n"	+ "link.AltPheno requires 0-1 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n" + "   (2) directory (i.e. dir=" + dir
										+ " (default))\n" + "   (3) parse results (i.e. -parse (not the default))\n"
										+ "   (4) zip results (i.e. -zip (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-parse")) {
				parse = true;
				numArgs--;
			} else if (arg.startsWith("-zip")) {
				zip = true;
				numArgs--;
			} else if (arg.startsWith("condense=")) {
				condense = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (!dir.equals("") && (!new File(dir).exists() || !new File(dir).isDirectory())) {
			System.err.println("Warning - using current directory instead of " + dir);
			dir = "";
		}

		try {
			if (parse) {
				parse(dir, filename);
			} else if (zip) {
				zip(dir, filename);
			} else if (!condense.equals("")) {
				condense(dir, condense);
			} else {
				rerunAnalyses(dir, filename);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
