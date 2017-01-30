package org.genvisis.cnv.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.analysis.FilterCalls;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

public class QuantiSNP {
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\00src\\";
	public static final String LINUX_DIRECTORY = "/work/parkinsons/gwas/Genotype and intensity data files/Final_Reports";
	// public static final String LINUX_DIRECTORY = "/work/parkinsons/gwas/Genotype and intensity data
	// files/Final_Reports/test";
	// public static final String ROOT_DIR = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\CNV\\Software\\QuantiSNP\\QuantiSNP_v1.1_windows\\";
	// public static final String ROOT_DIR = "Q:\\parkinsons\\cnvs\\quantisnp\\";
	// public static final String ROOT_DIR = "C:\\QuantiSNP\\";
	public static final String ROOT_DIR =
																			"C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\quantisnp\\noGenderProblems\\";
	public static final String OUTPUT_DIR = "output\\";
	public static final String CNV_DIRECTORY = "quanti_data/";
	// public static final String DEFAULT_OUTPUT = "output_e10/";
	public static final String DEFAULT_OUTPUT = "output/";
	public static final String MARKER_FILE = "cnvMarkers.dat";
	public static final String MARKER_POSITIONS = "markerPositions.dat";
	public static final String LOOKUP_PLUS_GENDER_FILE = "LookupPlusGender.xln";
	public static final double BAYES_FACTOR_CUTOFF = 10;
	// public static final int EM_ITERATIONS = 10;
	public static final int EM_ITERATIONS = 25;
	public static final String[] FIELDS = {	"SNP Name", "Sample ID", "Sample Name", "GC Score",
																					"Allele1 - Forward", "Allele2 - Forward", "Allele1 - AB",
																					"Allele2 - AB", "B Allele Freq", "Log R Ratio"};
	public static final String[] LOOKUP_HEADER = {"Sample_ID", "FID", "IID", "Gender"};

	public static void createFiles() {
		BufferedReader reader;
		PrintWriter writer, cnvWriter = null;
		String[] line, snpNames = null;
		Vector<String> vNames = new Vector<String>();
		int count, version, index;
		String id = "", trav, temp;
		boolean pdgwas;
		int[] indices;
		Hashtable<String, String> snpPoslar = new Hashtable<String, String>();

		if (new File(WINDOWS_DIRECTORY).exists()) {
			trav = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			trav = LINUX_DIRECTORY;
		} else {
			trav = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked)");
			System.err.println(WINDOWS_DIRECTORY);
			System.err.println(LINUX_DIRECTORY);
			System.exit(1);
		}

		File[] files = new File(trav).listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return filename.startsWith("Myers") && filename.endsWith(".csv");
				// return filename.endsWith(".csv");
			}
		});

		System.out.println(ext.getTime());
		System.out.println("Found " + files.length + " files to parse");
		try {
			reader = new BufferedReader(new FileReader(files[0]));
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.contains("SNP Name") && !temp.contains("Sample ID"));

			indices = ext.indexFactors(FIELDS, temp.trim().split(","), false, true);

			index = indices[ext.indexOfStr("SNP Name", FIELDS)];
			while (reader.ready()) {
				vNames.add(reader.readLine().trim().split(",")[index]);
			}
			snpNames = ArrayUtils.toStringArray(vNames);
			writer = new PrintWriter(new FileWriter(MARKER_FILE));
			for (String snpName : snpNames) {
				writer.println(snpName);
			}
			writer.close();

			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + files[0] + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + files[0] + "\"");
			System.exit(2);
		}

		System.out.println("There are " + snpNames.length + " markers being processed");
		new File(CNV_DIRECTORY).mkdirs();
		snpPoslar = HashVec.loadFileToHashString(MARKER_POSITIONS, 0, new int[] {1, 2}, "\t", false);

		try {
			for (File file : files) {
				pdgwas = file.getName().contains("Myers");
				try {
					reader = new BufferedReader(new FileReader(file));
					do {
						temp = reader.readLine();
					} while (reader.ready() && !temp.contains("SNP Name") && !temp.contains("Sample ID"));

					if (!reader.ready()) {
						System.err.println("Error - went through enitre file without finding a line containing both 'SNP Name' and 'Sample Name/Sample ID'");
						System.exit(1);
					}
					indices = ext.indexFactors(FIELDS, temp.trim().split(","), false, true);

					count = 0;
					while (reader.ready()) {
						line = reader.readLine().split(",");
						if (pdgwas) {
							trav = line[indices[ext.indexOfStr(	"Sample Name",
																									FIELDS)]].substring(0,
																																			line[indices[ext.indexOfStr("Sample Name",
																																																	FIELDS)]].indexOf("@"));
						} else {
							trav = line[indices[ext.indexOfStr("Sample ID", FIELDS)]];
						}
						if (count == 0) {
							id = trav;
							version = 0;
							while (new File(CNV_DIRECTORY	+ id + (version == 0 ? "" : "." + version)
															+ ".qs").exists()) {
								version++;
							}
							cnvWriter = new PrintWriter(new FileWriter(CNV_DIRECTORY	+ id
																													+ (version == 0 ? "" : "." + version)
																													+ ".qs"));
							cnvWriter.println("Name\tChr\tPosition\t"	+ id + ".Log R Ratio\t" + id
																+ ".B Allele Freq");
						} else if (!trav.equals(id)) {
							System.err.println("Found "	+ trav + " -- expecting " + id + " in file "
																	+ file.getName());
							System.exit(1);
						}
						if (!snpNames[count].equals(line[indices[ext.indexOfStr("SNP Name", FIELDS)]])) {
							System.err.println("Found "	+ line[indices[ext.indexOfStr("SNP Name", FIELDS)]]
																	+ " -- expecting " + snpNames[count] + " in file "
																	+ file.getName());
							System.exit(1);
						}

						temp = line[indices[ext.indexOfStr("SNP Name", FIELDS)]];
						if (!line[indices[ext.indexOfStr("Log R Ratio", FIELDS)]].equals("NaN")) {
							cnvWriter.print(temp + "\t" + snpPoslar.get(temp));
							cnvWriter.print("\t" + line[indices[ext.indexOfStr("Log R Ratio", FIELDS)]]);
							cnvWriter.print("\t" + line[indices[ext.indexOfStr("B Allele Freq", FIELDS)]]);
							cnvWriter.println();
						}

						count++;
					}

					reader.close();
					cnvWriter.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ file.getName()
															+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + file.getName() + "\"");
					System.exit(2);
				}
			}

			System.out.println(ext.getTime());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createBatch(	String rootDirectory, String inputDirectory,
																	String outputDirectory, int numBatches) {
		Vector<String[]> v = new Vector<String[]>();
		Hashtable<String, String> genders;
		String[] inputs, outputs;
		String commands, gender;

		// genders = HashVec.loadFileToHashString(, "DNA", new String[] {"CLASS=Gender"}, "");

		genders = HashVec.loadFileToHashString(rootDirectory	+ LOOKUP_PLUS_GENDER_FILE, 0,
																						new int[] {3}, "\t", true);

		inputs = new File(rootDirectory + inputDirectory).list(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return filename.endsWith(".qs");
			}
		});

		outputs = new File(rootDirectory + outputDirectory).list(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return filename.endsWith("_output.out");
			}
		});

		System.out.println("Found "	+ inputs.length + " samples, as well as results for "
												+ outputs.length + " that have been done (not necessarily the same ones)");

		for (String input : inputs) {
			if (ext.indexOfStr(ext.rootOf(input) + "_output.out", outputs) == -1) {
				gender = genders.get(ext.rootOf(input));
				if (gender.equals("M")) {
					gender = "male";
				} else if (gender.equals("F")) {
					gender = "female";
					// } else if (gender==null) {
					// System.err.println("Error - no gender found for subject '"+ext.rootOf(inputs[i])+"'");
				} else {
					System.err.println("Error - '" + gender + "' is not a valid gender (expecting M/F)");
				}

				v.add(new String[] {ext.rootOf(input), gender});
			}
		}

		System.out.println("Made "	+ numBatches + " batch files that will take care of the " + v.size()
												+ " files yet to parse");

		commands = "quantisnp.exe --config ../windows/config.dat --output "	+ OUTPUT_DIR
								+ "[%0].out --sampleid AA --gender [%1] --emiters " + EM_ITERATIONS
								+ " --Lsetting 2000000 --maxcopy 3 --printRS --doGCcorrect --printGenotypes --gcdir ../gc/b36/ --input-files ../source/[%0].qs 300\n\n";
		Files.batchIt("batch", null, numBatches, commands, Matrix.toStringArrays(v));
		for (int i = 0; i < numBatches; i++) {
			try {
				new File("batch." + (i + 1)).renameTo(new File("batch" + (i + 1) + ".bat"));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public static void parseResults(String rootDirectory, String outputDirectory, String included,
																	double bayesCutoff) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> lookupDNAtoSubject = new Hashtable<String, String>();
		Hashtable<String, String> lookupSubjectToDNA = new Hashtable<String, String>();
		String[] outputs, inds;

		try {
			reader = new BufferedReader(new FileReader(rootDirectory + LOOKUP_PLUS_GENDER_FILE));
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), LOOKUP_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				lookupDNAtoSubject.put(line[0], line[1] + "\t" + line[2]);
				lookupSubjectToDNA.put(line[1] + "\t" + line[2], line[0]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ rootDirectory + LOOKUP_PLUS_GENDER_FILE
													+ "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + rootDirectory + LOOKUP_PLUS_GENDER_FILE + "\"");
			System.exit(2);
		}

		outputs = new File(rootDirectory + outputDirectory).list(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return filename.endsWith("_output.out");
			}
		});

		if (included == null) {
			System.out.println("Using all indivduals with result files (n=" + outputs.length + ")");
			inds = new String[outputs.length];
			for (int i = 0; i < outputs.length; i++) {
				inds[i] = outputs[i].substring(0, outputs[i].lastIndexOf("_"));
			}
		} else {
			inds = ArrayUtils.toStringArray(HashVec.loadFileToVec(	rootDirectory, false, new int[] {0, 1}, true,
																												false));
			System.out.println("Using all indivduals listed in '"	+ included + "' (n=" + inds.length
													+ ")");
			for (String ind : inds) {
				if (ext.indexOfStr(ind + "_output.out", outputs) == -1) {
					System.err.println("Error - '"	+ ind
															+ "' was found in the list of indiviudals to parse, but no results were found for this sample");
				}
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(rootDirectory + "plink.cnv"));
			writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER));
			for (String ind : inds) {
				try {
					reader = new BufferedReader(new FileReader(rootDirectory	+ outputDirectory + ind
																											+ "_output.out"));
					reader.readLine();
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (line[9].equals("Inf") || Double.parseDouble(line[9]) > bayesCutoff) {
							writer.println(lookupDNAtoSubject.get(ind)	+ "\t" + line[1] + "\t" + line[2] + "\t"
															+ line[3] + "\t" + line[8] + "\t"
															+ (line[9].equals("Inf") ? "999" : line[9]) + "\t" + line[7]);
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ rootDirectory + outputDirectory + ind
															+ "_output.out"
															+ "\" not found in directory... won't be in final plink.cnv file");
				} catch (IOException ioe) {
					System.err.println("Error reading file \""	+ rootDirectory + outputDirectory + ind
															+ "_output.out" + "\"");
					System.exit(2);
				}
				// writer.println(inds[i]);
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error writing to file");
		}
	}

	public static void filters(String rootDirectory) {
		FilterCalls.stdFilters(rootDirectory, "plink.cnv", true, null, 36);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String root = ROOT_DIR;
		String inputs = CNV_DIRECTORY;
		String outputs = DEFAULT_OUTPUT;
		String include = null;
		int numBatches = 6;
		boolean createFiles = false;
		boolean createBatch = true;
		boolean parseResults = false;
		boolean filter = false;
		double bayesCutoff = BAYES_FACTOR_CUTOFF;

		String usage = "\\n"	+ "park.cnv.QuantiSNP requires 0-1 arguments\n"
										+ "   (1) root directory of CNV files (i.e. root=" + root + " (default))\n"
										+ "   (2) -createFiles (not the default)\n"
										+ "   (3) number of batches to create (i.e. batches=" + numBatches
										+ " (not the default))\n" + "   (4) -parseResults (not the default)\n"
										+ "   (5) file with list of inidiuvdals to parse (i.e. include=keeps.txt (default is to include all with results))\n"
										+ "   (6) bayes factor threshold for inclusion (i.e. bf=" + bayesCutoff
										+ " (default))\n" + "   (7) -filter (not the default)\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("root=")) {
				root = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("batches=")) {
				numBatches = Integer.parseInt(arg.split("=")[1]);
				createBatch = true;
				numArgs--;
			} else if (arg.startsWith("-createFiles")) {
				createFiles = true;
				numArgs--;
			} else if (arg.startsWith("-parseResults")) {
				parseResults = true;
				numArgs--;
			} else if (arg.startsWith("include=")) {
				include = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("bf=")) {
				bayesCutoff = Double.parseDouble(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-filter")) {
				filter = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (createFiles) {
				createFiles();
			}
			if (createBatch) {
				createBatch(root, inputs, outputs, numBatches);
			}
			if (parseResults) {
				parseResults(root, outputs, include, bayesCutoff);
			}
			if (filter) {
				filters(root);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
