package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.parse.GenParser;
import org.genvisis.stats.ProbDist;

import com.google.common.primitives.Ints;

public class Probabel {
	public static final int PALOGIST = 0;
	public static final int PALINEAR = 1;
	public static final int PACOXPH = 2;
	// public static final String[] EXECS = {"palogist", "palinear", "pacoxph"};
	public static final String[] EXECS = {"/home/npankrat/bin/palogist",
																				"/home/npankrat/bin/palinear",
																				"/home/npankrat/bin/pacoxph"};
	public static final String[] LOGIST_OUTPUT_HEADER =
																										{	"name", "A1", "A2", "Freq1", "MAF", "Quality",
																											"Rsq", "n", "Mean_predictor_allele", "chrom",
																											"beta_SNP_add", "sebeta_SNP_add", "chi2_SNP"};

	public static void createQsubs(String pheno, int type, boolean minimac, int maxGb) {
		String commands;
		Vector<String> v;
		String[] markerNames;
		long size;
		int numSplits, increment, node, antinode;

		if (type == 3) {
			type = ArrayUtils.determineType(pheno, 1, false);
		}

		node = 0;
		antinode = 6;
		v = new Vector<String>();
		System.out.println("Limit: "	+ ext.addCommas((long) (maxGb * Math.pow(2, 30)) / 1000000)
												+ " megabytes");
		for (int chr = 1; chr <= 22; chr++) {
			size = new File("chr" + chr + "/" + (minimac	? "chr" + chr + ".dose"
																										: "MACH_step2_chr" + chr + ".mldose")).length();
			System.out.println("chr"	+ chr + "/"
													+ (minimac ? "chr" + chr + ".dose" : "MACH_step2_chr" + chr + ".mldose")
													+ ": " + ext.addCommas(size) + " bytes");
			if (size > maxGb * Math.pow(2, 30)) {
				numSplits = (int) Math.ceil(size / (maxGb * Math.pow(2, 30)));
				System.out.println("Splitting chr" + chr + " into " + numSplits + " pieces");
				if (minimac) {
					markerNames = new SnpMarkerSet("chr"	+ chr + "/chr" + chr + ".info",
																					SnpMarkerSet.MINIMAC_INFO_FORMAT, false,
																					new Logger()).getMarkerNames();
				} else {
					markerNames = new SnpMarkerSet("chr"	+ chr + "/MACH_step2_chr" + chr + ".mlinfo", false,
																					new Logger()).getMarkerNames();
				}
				increment = (int) Math.ceil((double) markerNames.length / (double) numSplits);
				for (int i = 0; i < numSplits; i++) {
					commands = "cd chr#\n";
					commands += "cp chr#.info chr#.minfo\n";
					if (!new File("quan" + (i + 1) + ".txt").exists()
								|| !new File("quan" + (i + 1) + ".pinfo").exists()
							|| new File("quan" + (i + 1) + ".pinfo").length() == 0) {
						if (chr > 5) {
							commands += "sed -n '"	+ (i * increment + 1) + ","
													+ (i == numSplits - 1 ? markerNames.length + 1 : (i + 1) * increment)
													+ " p' < " + Minimac.REF_SNPS_ROOT + " > quan" + (i + 1) + ".txt\n";
						}
						// if you get lots of errors about non-existent markers at this stage, check to make
						// sure the right reference panel is hard coded
						if (chr > 5) {
							commands += Files.getRunString()
														+ " filesys.DosageData dosageFile=chr#.dose idFile=plink.fam mapFile=chr#.minfo out=chr#_quan"
													+ (i + 1) + ".dose extract=quan" + (i + 1) + ".txt -awk\n";
						}
						// commands += "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}'
						// chr#_quan"+(i+1)+".minfo > chr#_quan"+(i+1)+".pinfo\n";
						commands += Files.getRunString()	+ " filesys.SnpMarkerSet file=chr#.minfo extract=quan"
												+ (i + 1) + ".txt out=chr#_quan" + (i + 1) + ".pinfo\n";
					}
					commands += "/home/npankrat/bin/palogist -p ../"	+ pheno + " -d chr#_quan" + (i + 1)
											+ ".dose -i chr#_quan" + (i + 1) + ".pinfo -c 1 -o " + ext.rootOf(pheno)
											+ "_chr#_quan" + (i + 1) + "\n";
					commands += "cd ..\n";

					v.add(Files.qsub(	"", "q" + (i + 1) + "_chr#_" + ext.rootOf(pheno), chr, chr, commands,
														null, 10000, 12, "compute-0-" + node + ".local")[0]);
					// node++;
					if (node > 4) {
						node = 0;
					}
				}
			} else {
				commands = "cd chr#\n";
				if (minimac) {
					commands +=
										"awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' chr#.info > chr#.pinfo\n";
				}
				commands += EXECS[type]	+ " -p ../" + pheno + " -d "
										+ (minimac ? "chr#.dose" : "MACH_step2_chr#.mldose") + " -i "
										+ (minimac ? "chr#.pinfo" : "MACH_step2_chr#.mlinfo") + " -c # -o "
										+ ext.rootOf(pheno) + "_chr#\n";
				commands += "cd ..\n";
				v.add(Files.qsub(	"", "chr#_" + ext.rootOf(pheno), chr, chr, commands, null, 10000, 12,
													"compute-0-" + antinode + ".local")[0]);
				antinode++;
				if (antinode > 12) {
					antinode = 6;
				}
			}
		}

		Files.writeArray(ArrayUtils.toStringArray(v), "master." + ext.rootOf(pheno));
		Files.chmod("master." + ext.rootOf(pheno));
	}

	public static void createQsubFilesForBlade(String filename) {
		String[][] parameters = new String[][] {{"1", "7", "4"}, {"2", "7", "4"}, {"3", "7", "4"},
																						{"4", "7", "4"}, {"5", "7", "4"}, {"6", "7", "4"},
																						{"7", "6", "4"}, {"8", "6", "4"}, {"9", "6", "1"},
																						{"10", "6", "1"}, {"11", "6", "1"}, {"12", "5", "1"},
																						{"13", "4", "1"}, {"14", "4", "1"}, {"15", "4", "1"},
																						{"16", "4", "1"}, {"17", "4", "1"}, {"18", "4", "1"},
																						{"19", "4", "1"}, {"20", "4", "1"}, {"21", "2", "1"},
																						{"22", "2", "1"}};
		String commands = "";

		commands += "#!/bin/bash -l\n";
		commands += "#PBS -l nodes=1:ppn=1,mem=[%1]gb,walltime=[%2]:00:00\n";
		commands += "#PBS -m abe\n";
		commands += "module load ProbABEL\n";
		commands += "echo \"start at: \" `date`\n";
		commands += "cd /scratch2/pankratz\n";
		commands += "pacoxph -p "	+ filename
								+ " -d chr[%0].aric.f3.mldose -i chr[%0].aric.f3.mlinfo -c [%0] -o chr[%0].out\n";
		commands += "echo \"end at: \" `date`\n";
		Files.batchIt("chr", "", 22, commands, parameters);

		for (int chr = 1; chr <= 22; chr++) {
			new File("chr." + chr).renameTo(new File("chr" + chr + ".qsub"));
		}
	}

	public static void createBatchFiles(String pheno, int type, int numBatches) {
		String[][] parameters;
		String commands = "";
		Vector<String[]> v = new Vector<String[]>();


		for (int i = 1; i <= 23; i++) {
			if (new File("chr" + i).exists()) {
				v.add(new String[] {i + ""});
			}
		}
		parameters = Matrix.toStringArrays(v);

		commands += "cp " + pheno + " chr[%0]\n";
		commands += "cd chr[%0]\n";
		commands += EXECS[type]	+ " -p " + pheno
								+ " -d MACH_step2_chr[%0].mldose -i MACH_step2_chr[%0].mlinfo -c [%0] -o "
								+ ext.rootOf(pheno) + "_chr[%0]\n";
		// commands += Files.getRunString() + " gwas.Probabel parse=[%0] pheno="+pheno+"\n";
		commands += "cd ..\n\n";
		Files.batchIt("prob", "", numBatches, commands, parameters);
	}

	public static void parsePvalues(String pheno, boolean minimac) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, markerNames;
		IntVector var, yok;
		int count, countErr;
		String filename;
		int quan;

		try {
			System.out.println("Collecting results from phenotype file '" + pheno + "'");
			writer = new PrintWriter(new FileWriter(ext.rootOf(pheno) + "_add.xln"));
			var = new IntVector();
			yok = new IntVector();
			for (int chr = 1; chr <= 22; chr++) {
				if (minimac) {
					if (new File("chr" + chr + "/" + "chr" + chr + ".info").exists()) {
						markerNames = HashVec.loadFileToStringArray("chr"	+ chr + "/" + "chr" + chr + ".info",
																												true, new int[] {0}, false);
					} else {
						System.err.println("Warning - failed to find chr"	+ chr + "/" + "chr" + chr
																+ ".info to check for completeness");
						markerNames = null;
					}
				} else {
					if (new File("chr" + chr + "/" + "MACH_step2_chr" + chr + ".mlinfo").exists()) {
						markerNames = HashVec.loadFileToStringArray("chr"	+ chr + "/" + "MACH_step2_chr" + chr
																												+ ".mlinfo", true, new int[] {0}, false);
					} else {
						System.err.println("Warning - failed to find chr"	+ chr + "/" + "MACH_step2_chr" + chr
																+ ".mlinfo to check for completeness");
						markerNames = null;
					}
				}
				filename = null;
				try {
					filename = "chr" + chr + "/" + ext.rootOf(pheno) + "_chr" + chr + "_add.out.txt";
					if (Files.exists(filename, false)) {
						quan = 0;
						reader = new BufferedReader(new FileReader(filename));
					} else {
						quan = 1;
						filename = "chr"	+ chr + "/" + ext.rootOf(pheno) + "_chr" + chr + "_quan" + quan
												+ "_add.out.txt";
						if (Files.exists(filename, false)) {
							reader = new BufferedReader(new FileReader(filename));
						} else {
							throw new FileNotFoundException();
						}
					}
					System.out.print("Parsing '" + filename + "'");
					line = reader.readLine().trim().split("[\\s]+");
					if (var.size() == 0) {
						writer.println(ArrayUtils.toStr(line) + "\tZscore_pval\tchi2_pval");
					}
					count = countErr = 0;
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (markerNames != null && count >= markerNames.length) {
							System.err.println("Error - There are more markers in '"	+ filename
																	+ "' than there are in "
																	+ (minimac	? "chr" + chr + ".info"
																							: "MACH_step2_chr" + chr + ".mlinfo"));
						} else if (markerNames != null && !line[0].equals(markerNames[count])) {
							if (countErr < 10) {
								System.err.println("Error - expecting marker '"	+ markerNames[count] + "' at line "
																		+ (count + 1) + " of " + filename + " but found '" + line[0]
																		+ "' instead");
							} else if (countErr == 10) {
								System.err.println("Error - etc...");
							}
							countErr++;
						}

						// writer.println(Array.toStr(line)+"\t"+(line[10].equals("nan")?"nan":ProbDist.NormDist(Double.parseDouble(line[10])/Double.parseDouble(line[11])))+"\t"+(line[12].equals("nan")||line[12].equals("-inf")?"nan":ProbDist.ChiDist(Double.parseDouble(line[12]),
						// 1)));
						writer.println(ArrayUtils.toStr(line)	+ "\t"
														+ (line[10].equals("nan")	? "nan"
																											: ProbDist.NormDist(Double.parseDouble(line[10])
																																					/ Double.parseDouble(line[11])))
														+ "\t"
														+ (line[12].equals("nan")
																|| line[12].equals("-inf")	? "nan"
																														: ProbDist.ChiDist(	Double.parseDouble(line[12]),
																																								1)));
						count++;
						if (!reader.ready() && count < markerNames.length && quan > 0) {
							quan++;
							filename = "chr"	+ chr + "/" + ext.rootOf(pheno) + "_chr" + chr + "_quan" + quan
													+ "_add.out.txt";
							if (Files.exists(filename, false)) {
								reader = new BufferedReader(new FileReader(filename));
								reader.readLine();
								System.out.print("Parsing '" + filename + "'");
							} else {
								System.err.println("Error - found results for earlier quantiles but could not find chr"
																			+ chr + "/" + ext.rootOf(pheno) + "_chr" + chr + "_quan" + quan
																		+ "_add.out.txt");
							}
						}
					}
					if (count < markerNames.length) {
						System.err.println("Warning - file(s) for chromosome " + chr + " may be truncated");
					}
					reader.close();
					var.add(chr);
				} catch (FileNotFoundException fnfe) {
					yok.add(chr);
				} catch (Exception e) {
					System.err.println("Error reading file \"" + filename + "\"");
				}
			}
			writer.close();
			System.out.println("Data for chromosomes " + ext.listRanges(Ints.toArray(var)));
			if (yok.size() > 0) {
				System.out.println("Missing data for chromosomes " + ext.listRanges(Ints.toArray(yok)));
			}
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(pheno) + "_add.xln");
			e.printStackTrace();
		}
	}

	public static void parseMetalFiles(String pheno) { // incomplete
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Vector<String> var, yok;

		try {
			System.out.println("Collecting results from phenotype file '" + pheno + "'");
			writer = new PrintWriter(new FileWriter(ext.rootOf(pheno) + "_add.xln"));
			var = new Vector<String>();
			yok = new Vector<String>();
			for (int chr = 1; chr <= 23; chr++) {
				try {
					reader = new BufferedReader(new FileReader("chr"	+ chr + "/" + ext.rootOf(pheno) + "_chr"
																											+ chr + "_add.out.txt"));
					line = reader.readLine().trim().split("[\\s]+");
					if (var.size() == 0) {
						writer.println(ArrayUtils.toStr(line) + "\tZscore_pval\tchi2_pval");
					}
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						writer.println(ArrayUtils.toStr(line)	+ "\t"
														+ (line[10].equals("nan")	? "nan"
																											: ProbDist.NormDist(Double.parseDouble(line[10])
																																					/ Double.parseDouble(line[11])))
														+ "\t"
														+ (line[12].equals("nan")	? "nan"
																											: ProbDist.ChiDist(	Double.parseDouble(line[12]),
																																					1)));
					}
					reader.close();
					var.add(chr + "");
				} catch (FileNotFoundException fnfe) {
					yok.add(chr + "");
				} catch (IOException ioe) {
					System.err.println("Error reading file \""	+ "chr" + chr + "/" + ext.rootOf(pheno)
															+ "_chr" + chr + "_add.out.txt" + "\"");
				}
			}
			writer.close();
			System.out.println("Data for chromosomes " + ext.listWithCommas(ArrayUtils.toStringArray(var)));
			if (yok.size() > 0) {
				System.out.println("Missing data for chromosomes "
														+ ext.listWithCommas(ArrayUtils.toStringArray(yok)));
			}
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(pheno) + "_add.xln");
			e.printStackTrace();
		}
	}

	public static void stratifyByClass(	String pheno, String optionalLookupFile, String classFile,
																			String selectedClass) {
		Hashtable<String, String> hash;
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;
		String[] lookup;

		hash = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(classFile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[0] + "\t" + line[1], line[2]);
				hash.put(line[1], line[2]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + classFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + classFile + "\"");
			System.exit(2);
		}

		if (optionalLookupFile == null) {
			lookup = null;
		} else if (Files.exists(optionalLookupFile, false)) {
			lookup = HashVec.loadFileToStringArray(optionalLookupFile, false, new int[] {0, 1}, false);
		} else {
			System.err.println("Error - could not find lookup file, hope everything is unique!");
			lookup = null;
		}

		try {
			reader = new BufferedReader(new FileReader(pheno));
			writer = new PrintWriter(new FileWriter(ext.rootOf(pheno, false)	+ "_" + selectedClass
																							+ ".dat"));
			writer.println(reader.readLine());
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (lookup == null) {
					if (hash.containsKey(line[0])) {
						if (!hash.get(line[0]).equalsIgnoreCase(selectedClass)) {
							line[1] = "NA";
						}
					} else {
						System.err.println("Error - no class information for indiviudal " + line[0]);
						line[1] = "NA";
					}
				} else {
					if (!line[0].equals(lookup[count].split("[\\s]+")[1])) {
						System.err.println("Error - lookup does not match phenotype file (at line "
																	+ (count + 1) + " found " + line[0] + ", expecting '"
																+ lookup[count].split("[\\s]+")[1] + "')");
						System.exit(1);
					}
					if (hash.containsKey(lookup[count])) {
						if (!hash.get(lookup[count]).equalsIgnoreCase(selectedClass)) {
							line[1] = "NA";
						}
					} else {
						System.err.println("Error - no class information for indiviudal "	+ lookup[count]
																+ " in the lookup file");
						line[1] = "NA";
					}
				}
				writer.println(ArrayUtils.toStr(line));
				count++;
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + pheno + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + pheno + "\"");
			System.exit(2);
		}
	}

	public static void compareClasses(String pheno, String optionalLookupFile, String classFile,
																		String compClasses) {
		String affectedClass, unaffectedClass;
		Hashtable<String, String> hash;
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;
		String[] lookup;

		line = compClasses.trim().split("\\|");
		affectedClass = line[0];
		unaffectedClass = line[1];

		hash = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(classFile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[0] + "\t" + line[1], line[2]);
				hash.put(line[1], line[2]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + classFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + classFile + "\"");
			System.exit(2);
		}

		if (optionalLookupFile == null) {
			lookup = null;
		} else if (Files.exists(optionalLookupFile, false)) {
			lookup = HashVec.loadFileToStringArray(optionalLookupFile, false, new int[] {0, 1}, false);
		} else {
			System.err.println("Error - could not find lookup file, hope everything is unique!");
			lookup = null;
		}

		try {
			reader = new BufferedReader(new FileReader(pheno));
			writer = new PrintWriter(new FileWriter(ext.rootOf(pheno, false)	+ "_" + affectedClass
																							+ "_vs_" + unaffectedClass + ".dat"));
			line = reader.readLine().trim().split("[\\s]+");
			writer.println(line[0] + "\t" + line[1]);
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (lookup == null) {
					if (hash.containsKey(line[0])) {
						if (hash.get(line[0]).equalsIgnoreCase(affectedClass)) {
							line[1] = "1";
						} else if (hash.get(line[0]).equalsIgnoreCase(unaffectedClass)) {
							line[1] = "0";
						} else {
							line[1] = "NA";
						}
					} else {
						System.err.println("Error - no class information for indiviudal " + line[0]);
						line[1] = "NA";
					}
				} else {
					if (!line[0].equals(lookup[count].split("[\\s]+")[1])) {
						System.err.println("Error - lookup does not match phenotype file (at line "
																	+ (count + 1) + " found " + line[0] + ", expecting '"
																+ lookup[count].split("[\\s]+")[1] + "')");
						System.exit(1);
					}
					if (hash.containsKey(lookup[count])) {
						if (hash.get(lookup[count]).equalsIgnoreCase(affectedClass)) {
							line[1] = "1";
						} else if (hash.get(lookup[count]).equalsIgnoreCase(unaffectedClass)) {
							line[1] = "0";
						} else {
							line[1] = "NA";
						}
					} else {
						System.err.println("Error - no class information for indiviudal "	+ lookup[count]
																+ " in the lookup file");
						line[1] = "NA";
					}
				}
				writer.println(line[0] + "\t" + line[1]);
				count++;
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + pheno + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + pheno + "\"");
			System.exit(2);
		}
	}

	public static void exploreStrata() {
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\analysisOfImputation\\local\\";
		// String dir = "";
		String[] studies = {"CIDR", "Fung", "Sing550", "LEAPS", "Miami", "NGRC"};
		String[] classes = {"Ashk", "Italian", "EastEur", "Scand", "Germ", "British"};
		String[] datasets = {"top4K"};
		PrintWriter writer;
		String[] line;
		String trav;

		for (String studie : studies) {
			for (int j = 0; j < classes.length; j++) {
				stratifyByClass(dir	+ studie + "/" + studie + "_Aff_final.dat", dir + studie + "/list.txt",
												dir + "EuropeanClasses.txt", classes[j]);
				for (String dataset : datasets) {
					CmdLine.run(EXECS[0]	+ " -p " + studie + "_Aff_final_" + classes[j] + ".dat" + " -d "
											+ dataset + ".mldose -i " + dataset + ".mlinfo -c 1 -o " + studie
											+ "_Aff_final_" + classes[j] + "_" + dataset + "\n", dir + studie);
					line = new String[] {dir	+ studie + "/" + studie + "_Aff_final_" + classes[j] + "_"
																+ dataset + "_add.out.txt",
																"out=" + dir + studie + "/" + dataset + "_" + classes[j] + ".metal",
																"!6>0.30", "!10>-5", "!10<5", "!11!nan", "!11!0", "0", "1", "2",
																"10", "11"};
					GenParser.parse(line,
													new Logger(dir + studie + "/" + dataset + "_" + classes[j] + ".log"));
				}

				if (j < classes.length - 1) {
					compareClasses(dir	+ studie + "/" + studie + "_Aff_final.dat", dir + studie + "/list.txt",
													dir + "EuropeanClasses.txt",
													classes[j] + "|" + classes[classes.length - 1]);
					for (String dataset : datasets) {
						CmdLine.run(EXECS[0]	+ " -p " + studie + "_Aff_final_" + classes[j] + "_vs_"
												+ classes[classes.length - 1] + ".dat" + " -d " + dataset + ".mldose -i "
												+ dataset + ".mlinfo -c 1 -o " + studie + "_Aff_final_" + classes[j]
												+ "_vs_" + classes[classes.length - 1] + "_" + dataset + "\n",
												dir + studie);
						line = new String[] {dir	+ studie + "/" + studie + "_Aff_final_" + classes[j] + "_vs_"
																	+ classes[classes.length - 1] + "_" + dataset
																	+ "_add.out.txt",
																	"out="						+ dir + studie + "/" + dataset + "_" + classes[j] + "_vs_"
																										+ classes[classes.length - 1] + ".metal",
																	"!6>0.30", "!10>-5", "!10<5", "!11!nan", "!11!0", "0", "1", "2",
																	"10", "11"};
						GenParser.parse(line, new Logger(dir	+ studie + "/" + dataset + "_" + classes[j]
																							+ "_vs_" + classes[classes.length - 1] + ".log"));
					}
				}
			}
		}
		for (String dataset : datasets) {
			for (int step = 0; step <= 1; step++) {
				for (int j = 0; j < (step == 0 ? classes.length : classes.length - 1); j++) {
					if (step == 0) {
						trav = dataset + "_" + classes[j];
					} else {
						trav = dataset + "_" + classes[j] + "_vs_" + classes[classes.length - 1];
					}
					try {
						writer = new PrintWriter(new FileWriter(dir + trav + ".script"));
						writer.println("MARKER name"	+ "\n" + "ALLELE A1 A2" + "\n" + "EFFECT beta_SNP_add"
														+ "\n" + "STDERR sebeta_SNP_add" + "\n" + "SCHEME STDERR" + "\n"
														+ "GENOMICCONTROL ON" + "\n");
						for (String studie : studies) {
							writer.println("PROCESS " + studie + "/" + trav + ".metal");
						}
						writer.println("\n" + "\n" + "OUTFILE " + trav + " .tbl" + "\n" + "ANALYZE" + "\n");
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing to " + dir + classes[j] + ".script");
						e.printStackTrace();
					}
					CmdLine.run("metal < " + trav + ".script", dir);
				}
			}
			try {
				writer = new PrintWriter(new FileWriter(dir + "parse_" + dataset + ".crf"));
				writer.println("hits\n" + dataset + ".txt 0");

				for (int j = 0; j < classes.length - 1; j++) {
					writer.println(dataset	+ "_" + classes[j] + "_vs_" + classes[classes.length - 1]
													+ "1.tbl 0 5=" + classes[j] + "_vs_" + classes[classes.length - 1]
													+ "_pval");
				}
				for (int j = 0; j < classes.length - 1; j++) {
					writer.println(dataset	+ "_" + classes[j] + "_vs_" + classes[classes.length - 1]
													+ "1.tbl 0 3=" + classes[j] + "_vs_" + classes[classes.length - 1]
													+ "_beta");
				}
				for (String classe : classes) {
					writer.println(dataset + "_" + classe + "1.tbl 0 3=" + classe + "_beta");
				}
				for (String classe : classes) {
					writer.println(dataset + "_" + classe + "1.tbl 0 5=" + classe + "_pval");
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + "parse_" + dataset + ".crf");
				e.printStackTrace();
			}
			System.out.println("still need to run db parse_" + dataset + ".crf");
		}
	}

	private static String[] getFileRoots(String sourceDir) {
		File dir = new File(sourceDir);
		if (!dir.exists()) {
			System.err.println(ext.getTime()	+ "Error - source directory {" + sourceDir
													+ "} doesn't exist");
			return null;
		}

		String[] files = dir.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				String root = ext.rootOf(name, true);
				return (new File(dir + root + ".dose")).exists()
								&& (new File(dir + root + ".info")).exists();
			}
		});
		HashSet<String> fileRoots = new HashSet<String>();
		for (String f : files) {
			fileRoots.add(ext.rootOf(f, true));
		}
		return fileRoots.toArray(new String[fileRoots.size()]);
	}

	private static String[] getFileRoots(String sourceDoseDir, String sourceInfoDir) {
		File doseDir = new File(sourceDoseDir);
		File infoDir = new File(sourceInfoDir);

		if (!doseDir.exists()) {
			System.err.println(ext.getTime()	+ "Error - dosage source directory {" + sourceDoseDir
													+ "} doesn't exist");
			return null;
		}
		if (!infoDir.exists()) {
			System.err.println(ext.getTime()	+ "Error - info source directory {" + sourceDoseDir
													+ "} doesn't exist");
			return null;
		}

		HashSet<String> dosageRoots = new HashSet<String>();
		HashSet<String> infoRoots = new HashSet<String>();
		HashSet<String> roots = new HashSet<String>();

		for (String f : doseDir.list()) {
			if (f.endsWith(".dose")) {
				dosageRoots.add(ext.rootOf(f, true));
			}
		}
		for (String f : infoDir.list()) {
			if (f.endsWith(".info")) {
				infoRoots.add(ext.rootOf(f, true));
			}
		}

		roots.addAll(dosageRoots);
		roots.retainAll(infoRoots);

		return roots.toArray(new String[roots.size()]);
	}

	/**
	 * Given a source file directory containing .dose and .info files, split files into chunks of
	 * 500,000 SNP dosage values
	 *
	 * @param fileDir
	 * @param chunkDir
	 * @param overwrite
	 */
	public static void chunkFiles(String fileDir, String chunkDir, boolean overwrite) {
		String[] roots = getFileRoots(fileDir);
		splitInfoFiles(fileDir, chunkDir, roots, overwrite);
		splitDoseFiles(fileDir, chunkDir, roots, overwrite);
	}

	/**
	 * Given two source directories, one containing .dose and the other .info files, split files into
	 * chunks of 500,000 SNP dosage values
	 *
	 * @param dosageFileDir
	 * @param infoFileDir
	 * @param chunkDir
	 * @param overwrite
	 */
	public static void chunkFiles(String dosageFileDir, String infoFileDir, String chunkDir,
																boolean overwrite) {
		String[] roots = getFileRoots(dosageFileDir, infoFileDir);
		splitInfoFiles(infoFileDir, chunkDir, roots, overwrite);
		splitDoseFiles(dosageFileDir, chunkDir, roots, overwrite);
	}

	private static void splitDoseFiles(	String fileDir, String outDir, String[] fileRoots,
																			boolean overwrite) {
		for (String fileRoot : fileRoots) {
			System.out.println(ext.getTime() + "]\tChunking .dose file " + fileRoot);
			try {
				BufferedReader reader = Files.getAppropriateReader(fileDir + fileRoot + ".dose");
				ArrayList<PrintWriter> dosageWriters = new ArrayList<PrintWriter>();
				String line = null;
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("[\\s]+");

					int chunks = (parts.length - 2) / 500000 + 1;
					if (dosageWriters.size() == 0) {
						for (int i = 0; i < chunks; i++) {
							String file = outDir + fileRoot + "_" + (i * 500000) + ".dose";
							if (overwrite || !Files.exists(file)) {
								dosageWriters.add(Files.getAppropriateWriter(outDir	+ fileRoot + "_" + (i * 500000)
																															+ ".dose"));
							} else {
								System.out.println(ext.getTime() + "]\tSkipping file " + file);
								dosageWriters.add(null);
							}
						}
						System.out.println(ext.getTime() + "]\tWriting " + chunks + " chunks...");
					}

					String prepend = parts[0] + "\t" + parts[1]; // id and DOSE

					for (int i = 0; i < chunks; i++) {
						if (dosageWriters.get(i) == null) {
							continue;
						}
						StringBuilder sb = new StringBuilder();
						sb.append(prepend);
						for (int j = i * 500000; j < (i + 1) * 500000 && j < (parts.length - 2); j++) {
							sb.append("\t").append(parts[j + 2]);
						}
						dosageWriters.get(i).println(sb.toString());
					}

				}
				for (PrintWriter writer : dosageWriters) {
					if (writer == null) {
						continue;
					}
					writer.flush();
					writer.close();
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void splitInfoFiles(	String dir, String outDir, String[] fileRoots,
																			boolean overwrite) {
		for (String fileRoot : fileRoots) {
			System.out.println(ext.getTime() + "]\tProcessing .info file {" + fileRoot + "}");
			int lines = Files.countLines(dir + fileRoot + ".info", 1);
			try {
				BufferedReader fileReader = Files.getAppropriateReader(dir + fileRoot + ".info");
				String header = fileReader.readLine();
				for (int i = 0; i < lines; i += 500000) {
					String file = outDir + fileRoot + "_" + i + ".info";
					if (overwrite || !Files.exists(file)) {
						System.out.println(ext.getTime() + "]\tWriting lines " + i + " up to " + (i + 500000));
						PrintWriter writer = Files.getAppropriateWriter(file);
						writer.println(header);
						String line = null;
						for (int j = 0; j < 500000 && (line = fileReader.readLine()) != null; j++) {
							writer.println(line);
						}
						writer.flush();
						writer.close();
					} else {
						System.out.println(ext.getTime() + "]\tSkipping file " + fileRoot);
					}
				}
				fileReader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println(ext.getTime() + "]\t.info file processing complete.");
	}

	public static void appendPValsToResults(final String dir1, final String filePattern,
																					final boolean suffixTruePrefixFalse) {
		String dir = ext.verifyDirFormat(dir1);
		String[] resultFiles = new File(dir).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return (suffixTruePrefixFalse ? name.endsWith(filePattern) : name.startsWith(filePattern));
			}
		});

		for (String file : resultFiles) {
			System.out.println(ext.getTime() + "]\tProcessing file " + file);
			String root = ext.rootOf(dir + file, false);
			String out = root + "_processed.out";

			try {
				BufferedReader reader = Files.getAppropriateReader(file);
				PrintWriter writer = Files.getAppropriateWriter(out);

				String hdr = reader.readLine();
				String delim = ext.determineDelimiter(hdr);
				writer.println(hdr.replaceAll(delim, "\t") + "\tzScore_Pval\tchi2_pval");
				String line = null;
				while ((line = reader.readLine()) != null) {

					String[] parts = line.trim().split("[\\s]+");
					writer.println(ArrayUtils.toStr(parts)	+ "\t"
													+ (parts[9].equals("nan")	? "nan"
																										: ProbDist.NormDist(Double.parseDouble(parts[9])
																																				/ Double.parseDouble(parts[10])))
													+ "\t"
													+ (parts[11].equals("nan")
															|| parts[11].equals("-inf")	? "nan"
																													: ProbDist.ChiDist(	Double.parseDouble(parts[11]),
																																							1)));

				}
				writer.flush();
				writer.close();
				reader.close();
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String pheno = "ARIC_pheno.dat";
		String pheno = "LOAD_Aff_Sex_Age_APOE.dat";
		boolean blade = false;
		boolean qsub = false;
		boolean minimac = false;
		boolean parse = false;
		boolean process = false;
		boolean suffix = true;
		int type = 3;
		int numBatches = -1;
		String classFile = "EuropeanClasses.txt";
		String lookupFile = null;
		String stratClass = "";
		String compClasses = "";
		int max = 12;

		String usage = "\n"	+ "gwas.Probabel requires 0-1 arguments\n"
										+ "  Analysis types include: 0=logistic, 1=linear, 2=coxproportionalhazards, 3=autodetermine between 0 and 1\n"
										+ "   (1) phenotype filename (i.e. pheno=" + pheno + " (default))\n"
										+ "   (2) use minimac naming scheme (i.e. -minimac (not the default))\n"
										+ "   (3) batch for blade (i.e. -blade (not the default))\n"
										+ "   (4) batch local qsub (i.e. -qsub (not the default))\n"
										+ "   (5) analysis type (i.e. type=" + type + " (default))\n"
										+ "   (6) number of batches (i.e. batches=8 (not the default))\n"
										+ "   (7) max number of Gb of dosage file before split (i.e. max=" + max
										+ " (default))\n" + "   (8) parse p-values (i.e. -parse (not the default))\n"
										+ "  OR\n" + "   (1) class file with which to compare or stratify (i.e. class="
										+ classFile + " (not the default))\n"
										+ "   (2) stratify sample by class (i.e. strat=Ashk (not the default))\n"
										+ "   (3) compare based on classes (i.e. comp=Ashk|British (not the default))\n"
										+ "   (4) (optional) lookup file in same order as phenotype file with FID and IID incase IIDs are not unique (i.e. lookup="
										+ lookupFile + " (default))\n" + "  OR\n"
										+ "   (1) append p-values to result files (i.e. -process (not the default))\n"
										+ "   (2) Result file directory (i.e. pheno= )\n"
										+ "   (3) Result file pattern to match, either a prefix or a suffix (i.e. lookup= )\n"
										+ "   (4) If the file pattern is a prefix (defaults to suffix) (i.e. -prefix (not the default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("pheno=")) {
				pheno = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("type=")) {
				type = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("batches=")) {
				numBatches = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-blade")) {
				blade = true;
				numArgs--;
			} else if (arg.startsWith("-qsub")) {
				qsub = true;
				numArgs--;
			} else if (arg.startsWith("-minimac")) {
				minimac = true;
				numArgs--;
			} else if (arg.startsWith("-parse")) {
				parse = true;
				numArgs--;
			} else if (arg.startsWith("-process")) {
				process = true;
				numArgs--;
			} else if (arg.startsWith("-prefix")) {
				suffix = false;
				numArgs--;
			} else if (arg.startsWith("strat=")) {
				stratClass = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("comp=")) {
				compClasses = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("lookup=")) {
				lookupFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("max=")) {
				max = ext.parseIntArg(arg);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		// for (int chr = 1; chr <= 22; chr++) { // if dosage files have been compressed, you'll need
		// the info file to double check order, etc.
		// System.out.println("zcat chr"+chr+".tar.gz | tar xvf - chr"+chr+"/chr"+chr+".info");
		// }
		// System.exit(1);

		try {
			if (!stratClass.equals("")) {
				stratifyByClass(pheno, lookupFile, classFile, stratClass);
			} else if (!compClasses.equals("")) {
				compareClasses(pheno, lookupFile, classFile, compClasses);
			} else if (blade) {
				if (minimac) {
					System.err.println("Error - minimac is only coded for regular qsub at the moment");
				} else {
					createQsubFilesForBlade(pheno);
				}
			} else if (qsub) {
				createQsubs(pheno, type, minimac, max);
			} else if (numBatches > 0) {
				createBatchFiles(pheno, type, numBatches);
			} else if (parse) {
				parsePvalues(pheno, minimac);
			} else if (process) {
				appendPValsToResults(pheno, lookupFile, suffix);
			} else {
				System.out.println("Warning - not enough flags to inititate a subroutine");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
