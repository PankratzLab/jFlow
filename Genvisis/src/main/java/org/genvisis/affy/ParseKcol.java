package org.genvisis.affy;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.common.*;

public class ParseKcol implements Runnable {
	public static final String[][] SNP_HEADER_OPTIONS = { { "SNP", "rsID", "ProbeSetName", "Name" } };
	public static final int EXP_NUM_MARKERS = 909622;
	public static final String[] FIELDS = { "Sample ID", "Sample Name" };
	public static final String[] DELIMITERS = { " " };
	public static final String OVERWRITE_OPTION_FILE = ".overwrite_option";
	public static final String CANCEL_OPTION_FILE = ".cancel_option";
	public static final String HOLD_OPTION_FILE = ".hold_option";
	public static final String FILENAME_AS_ID_OPTION = "[FILENAME_ROOT]";
	public static final String[][] SNP_TABLE_FIELDS = { { "Name", "SNP Name" }, { "Chr", "Chromosome" }, { "Position" } };
	public static final String[][] DATA_FIELDS = { { "GC Score", "GCscore", "Confidence" }, { "X Raw" }, { "Y Raw" }, { "X", "Chr", "Xvalue", "Log Ratio", "intensity_1", "Signal A" }, { "Y", "Yvalue", "Chr", "Strength", "intensity_2", "Signal B" }, { "Theta" }, { "R" }, { "B Allele Freq" }, { "Log R Ratio" } };
	public static final String[][] GENOTYPE_FIELDS = { { "Allele1 - Forward", "Position", "Allele1", "genotype1", "Call" }, { "Allele2 - Forward", "Position", "Allele2", "genotype2", "Forced Call" }, { "Allele1 - AB" }, { "Allele2 - AB" }, { "Forward Strand Base Calls" }, { "Call Codes" } };
	// Log R Ratio
	// B Allele Freq
	private Project proj;
	private String[] files;
	private String[] markerNames;
	private int[] keysKeys;
	private long fingerprint;
	private char[][] abLookup;
	private Hashtable<String, String> fixes;
	private long timeBegan;
	private int threadId;

	public ParseKcol(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, long fingerprint, Hashtable<String, String> fixes, long timeBegan) {
		this(proj, files, markerNames, keysKeys, abLookup, fingerprint, fixes, timeBegan, -1);
	}

	public ParseKcol(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, long fingerprint, Hashtable<String, String> fixes, long timeBegan, int threadId) {
		this.proj = proj;
		this.files = files;
		this.markerNames = markerNames;
		this.keysKeys = keysKeys;
		this.abLookup = abLookup;
		this.fingerprint = fingerprint;
		this.fixes = fixes;
		this.timeBegan = timeBegan;
		this.threadId = threadId;

		// TODO This seems not to be necessarily, because there is a check of "samples/.sampRAF".
		// if (Files.list(proj.getDir(proj.MARKER_DATA_DIRECTORY, true), MarkerData.MARKER_DATA_FILE_EXTENSION, proj.getJarStatus()).length>0) {
		// System.err.println("Error - Refusing to create new SampleList until the plots directory is either deleted or emptied; altering the SampleList will invalidate those files");
		// System.exit(1);
		// }

		// TODO There is a more complex scenario: what happens if you just want to add some new samples?
		// TransposeData.backupOlderFiles(proj.getDir(proj.SAMPLE_DIRECTORY, true), new String [] {"outliers.ser", ".sampRAF"}, true);
		// TransposeData.deleteOlderRafs(proj.getDir(proj.SAMPLE_DIRECTORY, true), new String[] {"outliers"}, new String[] {".ser"}, true);
	}

	public void run() {
		BufferedReader reader;
		String[] line;
		int count, snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		Sample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		String idHeader;
		String delimiter;
		String filename;
		Hashtable<String, Float> allOutliers;
		// String sourceExtension;
		// byte[] genos;

		// try {
		// PrintWriter writer = new PrintWriter(new FileWriter(files[0]+"_list.xln"));
		// for (int j = 0; j<files.length; j++) {
		// writer.println(files[j]);
		// }
		// writer.close();
		// return;
		// } catch (Exception e) {
		// System.err.println("Error writing to "+files[0]+"_list.xln");
		// e.printStackTrace();
		// }

		idHeader = proj.getProperty(proj.ID_HEADER);
		// sourceExtension = proj.getProperty(proj.SOURCE_FILENAME_EXTENSION);
		delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
		allOutliers = new Hashtable<String, Float>();
		try {
			for (int i = 0; i < files.length; i++) {
				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + CANCEL_OPTION_FILE).exists()) {
					return;
				}
				try {
					System.out.println(ext.getTime() + "\t" + (i + 1) + " of " + files.length);
					// reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));

					dataIndices = ext.indexFactors(DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(GENOTYPE_FIELDS, line, false, true, false, false);
					if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
						sampIndex = -7;
					} else {
						sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];
					}
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];

					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain " + Array.toStr(DATA_FIELDS[3], "/") + " and " + Array.toStr(DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[4] == -1 && (genotypeIndices[0] == -1 || genotypeIndices[1] == -1)) {
						System.err.println("Error - File format not consistent! The files need to contain " + Array.toStr(GENOTYPE_FIELDS[0], "/") + " and " + Array.toStr(GENOTYPE_FIELDS[1], "/") + " or " + GENOTYPE_FIELDS[4] + " (for single token calls)");
						return;
					}
					if (genotypeIndices[5] == -1 && (genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					sampleName = "";
					data = new float[DATA_FIELDS.length][];
					for (int j = 0; j < data.length; j++) {
						if (dataIndices[j] != -1) {
							data[j] = new float[markerNames.length];
						} else if (j == 7) {
							dataIndices[j] = 4;
							data[j] = new float[markerNames.length];
						} else if (j == 8) {
							dataIndices[j] = 3;
							data[j] = new float[markerNames.length];
						}
					}
					genotypes = new byte[2][];
					genotypes[0] = Array.byteArray(markerNames.length, (byte) 0); // used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
					if (!ignoreAB) {
						genotypes[1] = Array.byteArray(markerNames.length, (byte) -1); // used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
					}

					count = 0;
//					parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
					parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
					int linecount = 0;
					while (reader.ready()) {
						String testline = reader.readLine();
						linecount++;
						if (testline == null) {
							System.err.println("Warning - " + files[i] + " has a blank line on line " + linecount);
							break;
						}
						line = testline.split(delimiter);
						if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
							trav = files[i].substring(0, files[i].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
							;
						} else {
							if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
								System.err.println("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample");
								parseAtAt = false;
							}
							trav = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
						}

						if (count == 0) {
							sampleName = trav;
						} else if (!trav.equals(sampleName)) {
							System.err.println("Found more than one ID in file " + files[i] + "(found " + trav + ", expecting " + sampleName + ")");
							return;
						} else if (count > markerNames.length) {
							System.err.println("Error - expecting only " + markerNames.length + " markers and found more than that in file " + files[i]);
							return;
						} else if (!markerNames[count].equals(line[snpIndex]) && !markerNames[count].startsWith("Blank")) {
							System.err.println("Found " + line[snpIndex] + " at marker #" + (count + 1) + " in file " + files[i] + "; expecting " + markerNames[count]);
							return;
						}
						key = keysKeys[count];
						// System.out.println(count+"\t"+markerNames[count]+"\t"+key);
						for (int j = 0; j < DATA_FIELDS.length; j++) {
							try {
								if (dataIndices[j] != -1) {

									if (line[dataIndices[j]].equals("")) {
										data[j][key] = Float.NaN;
									} else if (j == 0) {
										data[j][key] = 0;
									} else if (j == 3 || j == 4) {
										data[j][key] = 0;
									} else if (j == 7) {
										if (line[dataIndices[j]].equals("2")) {
											data[j][key] = 0;
										} else {
											data[j][key] = Float.parseFloat(line[dataIndices[j]]);
										}
										// System.out.println(data[j][key]);
									} else if (j == 8) {
										data[j][key] = Float.parseFloat(line[dataIndices[j]]);
										// System.out.println(data[j][key]);
									} else {
										data[j][key] = 0;
									}
								}
							} catch (NumberFormatException nfe) {
								System.err.println("Error - failed to parse '" + line[dataIndices[j]] + "' into a valid " + Array.toStr(DATA_FIELDS[j], "/"));
								return;
							}
						}
						if (!ignoreAB) {

							// System.out.println("This is a problem at start:" + genotypeIndices[0] + "   " + genotypeIndices[1]);
							// String sta = line[genotypeIndices[0]] + " " + line[genotypeIndices[1]];
							String A = Character.toString(abLookup[count][0]);
							// String B = Character.toString(abLookup[count][1]);
							line[genotypeIndices[0]] = A;
							line[genotypeIndices[1]] = A;

						}
						genotypes[0][key] = (byte) ext.indexOfStr(line[genotypeIndices[0]] + line[genotypeIndices[1]], Sample.ALLELE_PAIRS);
						if (genotypes[0][key] == -1) {
							if (ext.indexOfStr(line[genotypeIndices[0]] + line[genotypeIndices[1]], Sample.ALT_NULLS) == -1) {
								System.err.println("Error - failed to lookup " + line[genotypeIndices[0]] + line[genotypeIndices[1]] + " for marker " + markerNames[count] + " of sample " + files[i]);
							} else {
								genotypes[0][key] = 0;
							}
						}
						if (ignoreAB) {
							// do nothing, will need to use these files to determine AB lookup table
						} else if (abLookup == null) {
							genotypes[1][key] = (byte) ext.indexOfStr(line[genotypeIndices[2]] + line[genotypeIndices[3]], Sample.AB_PAIRS);
						} else {
							if (genotypes[0][key] == 0) {
								genotypes[1][key] = -1;
							} else {
								genotypes[1][key] = 0;
								for (int j = 0; j < 2; j++) {
									if (line[genotypeIndices[j]].charAt(0) == abLookup[count][1]) {
										genotypes[1][key]++;
									} else if (line[genotypeIndices[j]].charAt(0) != abLookup[count][0]) {
										System.err.println("Error - alleles for individual '" + line[sampIndex] + "' (" + line[genotypeIndices[0]] + "/" + line[genotypeIndices[1]] + ") do not match up with the defined AB lookup alleles (" + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker " + markerNames[count]);
									}
								}
							}
						}

						// genos = parseGenotypes(line, genotypeIndices, ignoreAB, abLookup, count, sampleName, markerNames[count], files[i]);
						// genotypes[0][key] = genos[0];
						// if (!ignoreAB) {
						// genotypes[1][key] = genos[1];
						// }

						count++;
					}
					reader.close();
					if (count != markerNames.length) {
						System.err.println("Error - expecting " + markerNames.length + " markers and only found " + count + " in file " + files[i]);
						return;
					}

					if (fixes.containsKey(sampleName)) {
						sampleName = fixes.get(sampleName);
					}

					filename = determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true), sampleName, timeBegan);
					if (filename == null) {
						return;
					}

					samp = new Sample(sampleName, fingerprint, data, genotypes, true);
					// samp.serialize(proj.getDir(proj.SAMPLE_DIRECTORY, true) + trav + Sample.SAMPLE_DATA_FILE_EXTENSION);
					// samp.saveToRandomAccessFile(filename);
					samp.saveToRandomAccessFile(filename, allOutliers, sampleName); // TODO sampleIndex
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + files[i] + "\"");
					return;
				}
			}

			if (allOutliers.size() > 0) {
				if (threadId >= 0) {
					if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser").exists()) {
						System.err.println("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser");
						System.exit(1);
					} else {
						SerializedFiles.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser");
					}
				} else {
					if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser").exists()) {
						System.err.println("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser");
						System.exit(1);
					} else {
						SerializedFiles.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser");
					}
				}
			}

			SampleList.generateSampleList(proj).writeToTextFile(proj.PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");

			System.out.println(ext.getTime() + "\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static byte[] parseGenotypes(String[] line, int[] genotypeIndices, boolean ignoreAB, char[][] abLookup, int count, String sampleName, String markerName, String filename) {
		String genotype;
		byte genoForward, genoAB;

		if (genotypeIndices[4] >= 0) {
			genotype = line[genotypeIndices[4]];
		} else {
			genotype = line[genotypeIndices[0]] + line[genotypeIndices[1]];
		}
		genoForward = (byte) ext.indexOfStr(genotype, Sample.ALLELE_PAIRS);
		if (genoForward == -1) {
			if (ext.indexOfStr(genotype, Sample.ALT_NULLS) == -1) {
				// System.err.println("Error - failed to lookup "+genotype+" for marker "+markerName+" of sample "+filename);
				genoForward = -9;
			} else {
				genoForward = 0;
			}
		}
		if (ignoreAB) {
			genoAB = -1;
			// do nothing, will need to use these files to determine AB lookup table
		} else if (abLookup == null) {
			if (genotypeIndices[5] >= 0) {
				genotype = line[genotypeIndices[5]];
			} else {
				genotype = line[genotypeIndices[2]] + line[genotypeIndices[3]];
			}
			genoAB = (byte) ext.indexOfStr(genotype, Sample.AB_PAIRS);
			if (genoAB == -1) {
				if (ext.indexOfStr(genotype, Sample.ALT_NULLS) == -1) {
					System.err.println("Error - failed to lookup " + genotype + " for marker " + markerName + " of sample " + filename);
				}
			}
		} else {
			if (genoForward == 0) {
				genoAB = -1;
			} else {
				genoAB = 0;
				if (genotypeIndices[4] >= 0) {
					genotype = line[genotypeIndices[4]];
				} else {
					genotype = line[genotypeIndices[0]] + line[genotypeIndices[1]];
				}
				for (int j = 0; j < 2; j++) {
					if (genotype.charAt(j) == abLookup[count][1]) {
						genoAB++;
					} else if (genotype.charAt(j) != abLookup[count][0]) {
						System.err.println("Error - alleles for individual '" + sampleName + "' (" + genotype + ") do not match up with the defined AB lookup alleles (" + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker " + markerName);
					}
				}
			}
		}

		return new byte[] { genoForward == -9 ? 1 : genoForward, genoAB };
	}

	private static String determineFilename(String dir, String sampleName, long timeBegan) {
		String filename, trav;
		int version, versionToOverwrite;
		String[] overwriteOptions;
		int response;

		version = 0;
		versionToOverwrite = -1;
		do {
			trav = sampleName + (version == 0 ? "" : "." + version);
			version++;
			filename = dir + trav + Sample.SAMPLE_FILE_EXTENSION;
			if (new File(filename).exists() && new File(filename).lastModified() < timeBegan && versionToOverwrite == -1) {
				versionToOverwrite = version - 1;
			}
		} while (new File(filename).exists());

		overwriteOptions = new String[] { "Rename new file " + trav + Sample.SAMPLE_FILE_EXTENSION, "Overwrite existing file " + sampleName + (versionToOverwrite == 0 ? "" : "." + versionToOverwrite) + Sample.SAMPLE_FILE_EXTENSION, "Overwrite this and all future files", "Cancel parser" };

		if (versionToOverwrite != -1) {
			while (new File(dir + HOLD_OPTION_FILE).exists()) {
				try {
					Thread.sleep(500);
				} catch (InterruptedException ie) {
				}
			}
			if (new File(dir + CANCEL_OPTION_FILE).exists()) {
				return null;
			}

			Files.write("", dir + HOLD_OPTION_FILE);
			if (new File(dir + OVERWRITE_OPTION_FILE).exists()) {
				response = 1;
			} else {
				do {
					response = JOptionPane.showOptionDialog(null, "Error - the same sample name '" + sampleName + "' is being parsed again and the previous file existed before the current command began.\n" + "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted.\n" + "If you would like to start from scratch, the safest thing would be to cancel now and delete all files in the sample directory.\n" + "What would you like to do?", "What to do?", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[0]);
				} while (response == -1);
			}
			new File(dir + HOLD_OPTION_FILE).delete();
			switch (response) {
			case 0:
				return dir + trav + Sample.SAMPLE_FILE_EXTENSION;
			case 2:
				Files.write("", dir + OVERWRITE_OPTION_FILE);
			case 1:
				return dir + sampleName + (versionToOverwrite == 0 ? "" : "." + versionToOverwrite) + Sample.SAMPLE_FILE_EXTENSION;
			case 3:
				Files.write("", dir + CANCEL_OPTION_FILE);
				return null;
			default:
				JOptionPane.showMessageDialog(null, "Should be impossible to obtain this message (" + response + ")", "Error", JOptionPane.ERROR_MESSAGE);
				break;
			}
		}

		return dir + trav + Sample.SAMPLE_FILE_EXTENSION;
	}

	@SuppressWarnings("unchecked")
	public static void createFiles(Project proj, int numThreads) {
		BufferedReader reader;
		String[] line, markerNames, files;
		ArrayList<String> alNames;
		int snpIndex;
		int[] keys, keysKeys, indices;
		long fingerprint;
		Thread[] threads;
		Vector<Vector<String>> fileCabinet;
		String trav;
		int count;
		// Hashtable<String, Integer> markerNameHash;
		boolean abLookupRequired;
		char[][] lookup;
		Hashtable<String, String> fixes;
		String idHeader, delimiter;
		String temp;
		int[][] delimiterCounts;

		long timeBegan;
		boolean parseAtAt;
		int sampIndex;
		String sampleName;
		String[] overwriteOptions;
		int response;
		String[] filesToDelete;
		boolean complete;
		Hashtable<String, Float> allOutliers;

		timeBegan = new Date().getTime();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + OVERWRITE_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + HOLD_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + CANCEL_OPTION_FILE).delete();

		if (!proj.SOURCE_DIRECTORY.getValue(false, true).equals("") && !new File(proj.SOURCE_DIRECTORY.getValue(false, true)).exists()) {
			System.err.println("Error - the Project source location is invalid: " + proj.SOURCE_DIRECTORY.getValue(false, true));
			return;
		}

		if (!new File(proj.MARKER_POSITION_FILENAME.getValue(false, false)).exists()) {
			System.err.println("Error - missing markerPositions: " + proj.MARKER_POSITION_FILENAME.getValue(false, false));
			return;
		}

		delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
		idHeader = proj.getProperty(proj.ID_HEADER);
		System.out.println(ext.getTime() + "\tSearching for " + proj.getProperty(proj.SOURCE_FILENAME_EXTENSION) + " files in: " + proj.SOURCE_DIRECTORY.getValue(false, true));
		files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), "gw6_split", "", false, false);
		// Files.list(kColDir+dirList[0], prefix ,suffix,false ,false);
		System.out.println("\t\tFound " + files.length + " file" + (files.length == 1 ? "" : "s") + " with a " + proj.getProperty(proj.SOURCE_FILENAME_EXTENSION) + " extension");
		for (int i = 0; i < files.length; i++) {
			if (files[i].equals("Sample_Map.csv") || files[i].equals("SNP_Map.csv")) {
				files = Array.removeFromArray(files, i);
				i--;
			}
		}

		if (files.length == 0) {
			System.err.println("Error - no files to parse");
			return;
		}

		abLookupRequired = false;
		System.out.println("\t\tFound " + files.length + " file" + (files.length == 1 ? "" : "s") + " to parse");
		fixes = new Hashtable<String, String>();
		if (new File(proj.PROJECT_DIRECTORY.getValue() + "fixes.dat").exists()) {
			System.out.println("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
			fixes = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue() + "fixes.dat", false);
		} else {
			System.out.println("Did not find a 'fixes.dat' file; assuming you don't want to rename any IDs");
		}

		try {
			// reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[0]));
			reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
			System.out.println("Found appropriate reader for: " + proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
			count = 0;
			do {
				line = reader.readLine().trim().split(delimiter, -1);
				count++;
				// if (count < 20) {
				// System.out.println(Array.toStr(line));
				// }
			} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));

			if (!reader.ready()) {
				System.err.println("Error - reached the end of the file without finding a line with the following tokens: " + Array.toStr(SNP_HEADER_OPTIONS[0]));
				System.err.println("      - perhaps the delimiter is set incorrectly? Determing most stable delimiter...");

				reader.close();
				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
				delimiterCounts = new int[DELIMITERS.length][count];
				for (int i = 0; i < count; i++) {
					temp = reader.readLine();
					for (int j = 0; j < DELIMITERS.length; j++) {
						delimiterCounts[j][i] = ext.countInstancesOf(temp, DELIMITERS[j]);
					}
				}
				delimiter = null;
				for (int j = 0; j < DELIMITERS.length; j++) {
					if (Array.quantWithExtremeForTie(delimiterCounts[j], 0.5) > 4 && Array.quantWithExtremeForTie(delimiterCounts[j], 0.9) - Array.quantWithExtremeForTie(delimiterCounts[j], 0.1) == 0) {
						if (delimiter == null) {
							delimiter = DELIMITERS[j];
						} else {
							JOptionPane.showMessageDialog(null, "Could not auto-detect the delimiter used in the Final Reports file: could be '" + delimiter + "' or '" + DELIMITERS[j] + "'", "Error", JOptionPane.ERROR_MESSAGE);
							return;
						}
					}
				}
				reader.close();

				if (delimiter == null) {
					JOptionPane.showMessageDialog(null, "Failed to auto-detect the delimiter used in the Final Reports file; exitting", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}

				System.err.println("      - determined delimiter to be '" + delimiter + "'");

				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
				System.out.println(1);
			}

			// check immediately to make sure these fields are valid
			indices = ext.indexFactors(DATA_FIELDS, line, false, true, true, false); // dataIndices
			if (indices[3] == -1 || indices[4] == -1) {
				System.err.println("Error - at the very least the files need to contain " + Array.toStr(DATA_FIELDS[3], "/") + " and " + Array.toStr(DATA_FIELDS[4], "/"));
				System.err.println("      - failed to see that in " + files[0]);
				System.err.println(Array.toStr(line));

				return;
			}
			indices = ext.indexFactors(GENOTYPE_FIELDS, line, false, true, true, false); // genotypeIndices
			if (indices[4] == -1 && (indices[0] == -1 || indices[1] == -1)) {
				System.err.println("Error - the files need to contain " + Array.toStr(GENOTYPE_FIELDS[0], "/") + " and " + Array.toStr(GENOTYPE_FIELDS[1], "/") + " or " + GENOTYPE_FIELDS[4] + " (for single token calls)");
				return;
			}
			if (indices[5] == -1 && (indices[2] == -1 || indices[3] == -1)) {
				abLookupRequired = true;
			}

			snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, true, true)[0];

//			parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
			parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
			idHeader = proj.getProperty(proj.ID_HEADER);
			if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
				sampleName = files[0].substring(0, files[0].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
			} else {
				sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];
				reader.mark(1000);
				line = reader.readLine().split(delimiter);
				if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
					System.err.println("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample");
					parseAtAt = false;
				}
				sampleName = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
				reader.reset();
			}

			if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + sampleName + Sample.SAMPLE_FILE_EXTENSION).exists()) {

				overwriteOptions = new String[] { "Delete All", "Customize", "Cancel parser" };

				response = JOptionPane.showOptionDialog(null, "These data (at least the first sample '" + sampleName + "') have already been parsed.\n" + "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted.\n" + "If you would like to start from scratch, select \"Delete All\" earlier files.\n" + "Otherwise, cancel or you can \"Customize\" and determine what to do on a sample-by-sample basis.\n" + "What would you like to do?", "What to do?", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[2]);

				switch (response) {
				case -1:
					break;
				case 0:
					filesToDelete = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true), Sample.SAMPLE_FILE_EXTENSION, false);
					for (int i = 0; i < filesToDelete.length; i++) {
						new File(proj.SAMPLE_DIRECTORY.getValue(false, true) + filesToDelete[i]).delete();
					}
					new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").delete();
					break;
				case 1:
					// keep "outlier.ser"
					break;
				case 2:
					return;
				default:
					JOptionPane.showMessageDialog(null, "Should be impossible to obtain this message (" + response + ")", "Error", JOptionPane.ERROR_MESSAGE);
					break;
				}
			}
			TransposeData.deleteOlderRafs(proj.SAMPLE_DIRECTORY.getValue(true, true), new String[] { "outliers" }, new String[] { ".ser" }, true, new String[] { "outliers.ser" });

//			if (Boolean.parseBoolean(proj.getProperty(proj.LONG_FORMAT))) {
			if (proj.getProperty(proj.LONG_FORMAT)) {
				reader.close();
				createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter, abLookupRequired, timeBegan);
				return;
			}

			// count = 1;
			alNames = new ArrayList<String>(500000);
			int linecount = 0;
			while (reader.ready()) {
				String testline = reader.readLine();
				linecount++;
				if (testline == null) {
					System.err.println("Warning - " + files[0] + " has a blank line on line " + linecount);
					break;
				}
				// line = testline.split(delimiter);
				trav = testline.trim().split(delimiter)[snpIndex];
				if (trav.equals("") || trav.equals("0")) {
					trav = "Blank" + count;
					count++;
				}
				// if (markerNameHash.containsKey(trav)) {
				// done = true;
				// } else {
				// markerNameHash.put(trav, new Integer(markerNameHash.size()));
				// }
				alNames.add(trav);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + proj.SOURCE_DIRECTORY.getValue(false, true) + files[0] + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + proj.SOURCE_DIRECTORY.getValue(false, true) + files[0] + "\"");
			return;
		}

		// new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)).mkdirs();
		// new File(proj.getDir(proj.IND_DIRECTORY)).mkdirs();
		// new File(proj.getDir(proj.DATA_DIRECTORY)).mkdirs();

		// markerNames = Array.toStringArray(markerNameHash);
		markerNames = Array.toStringArray(alNames);
//		keys = Markers.orderMarkers(markerNames, proj.getFilename(proj.MARKER_POSITION_FILENAME), proj.getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
		keys = Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		if (keys == null) {
			return;
		}
		keysKeys = Sort.quicksort(keys); // very important
		fingerprint = proj.getMarkerSet().getFingerprint();
		System.out.println("There are " + markerNames.length + " markers being processed (fingerprint: " + fingerprint + ")");

		lookup = getABLookup(abLookupRequired, markerNames, proj);

		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i < numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i < files.length; i++) {
			fileCabinet.elementAt(i % numThreads).add(files[i]);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i < numThreads; i++) {
			threads[i] = new Thread(new ParseKcol(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), markerNames, keysKeys, lookup, fingerprint, fixes, timeBegan, i));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {
			}
		}

		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i < numThreads; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {
				}
			}
		}

		allOutliers = new Hashtable<String, Float>();
		for (int i = 0; i < numThreads; i++) {
			if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").exists()) {
				allOutliers.putAll((Hashtable<String, Float>) SerializedFiles.readSerial(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser"));
				new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").delete();
			}
		}
		if (allOutliers.size() > 0) {
			SerializedFiles.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		}

	}

	public static char[][] getABLookup(boolean abLookupRequired, String[] markerNames, Project proj) {
		ABLookup abLookup;
		char[][] lookup;

//		if (abLookupRequired && Files.exists(proj.getFilename(proj.AB_LOOKUP_FILENAME))) {
		if (abLookupRequired && Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
//			abLookup = new ABLookup(markerNames, proj.getFilename(proj.AB_LOOKUP_FILENAME), true, false, proj.getLog());
			abLookup = new ABLookup(markerNames, proj.AB_LOOKUP_FILENAME.getValue(), true, false, proj.getLog());
			lookup = abLookup.getLookup();
			if (lookup == null) {
				System.err.println("Warning - filed to provide columns \"" + GENOTYPE_FIELDS[2][0] + "\" / \"" + GENOTYPE_FIELDS[3][0] + "\" and the specificed AB_lookup file '" + proj.getProperty(proj.AB_LOOKUP_FILENAME) + "' does not exist; you'll need reconstruct the B allele for analysis");
			} else {
				abLookup.writeToFile(proj.PROJECT_DIRECTORY.getValue() + "checkAB.xln", proj.getLog());
			}
		} else {
			lookup = null;
		}

		return lookup;
	}

	public static void createFilesFromLongFormat(Project proj, String[] files, String idHeader, Hashtable<String, String> fixes, String delimiter, boolean abLookupRequired, long timeBegan) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, markerNames, list;
		long fingerprint;
		MarkerSet markerSet;
		Hashtable<String, Integer> markerIndices;
		int count;
		char[][] abLookup;
		String filename;
		Hashtable<String, Float> allOutliers;

		System.out.println("Parsing files using the Long Format algorithm");

//		Markers.orderMarkers(null, proj.getFilename(proj.MARKER_POSITION_FILENAME), proj.getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
		Markers.orderMarkers(null, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		fingerprint = proj.getMarkerSet().getFingerprint();

		abLookup = getABLookup(abLookupRequired, markerNames, proj);

		markerIndices = new Hashtable<String, Integer>();
		for (int i = 0; i < markerNames.length; i++) {
			markerIndices.put(markerNames[i], new Integer(i));
		}

		System.out.println("There were " + markerNames.length + " markers present in '" + proj.MARKERSET_FILENAME.getValue(true, true) + "' that will be processed from the source files (fingerprint: " + fingerprint + ")");

		int snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		Sample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		CountHash countHash, dupHash;
		boolean done;
		byte[] genos;

		count = 0;
		countHash = new CountHash();
		dupHash = new CountHash();
		allOutliers = new Hashtable<String, Float>();
		try {
			for (int i = 0; i < files.length; i++) {
				try {
					System.out.println(ext.getTime() + "\t" + (i + 1) + " of " + files.length + " (" + files[i] + ")");
					// reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || ext.indexOfStr(idHeader, line) == -1));

					System.err.println("Searching: " + Array.toStr(line));
					dataIndices = ext.indexFactors(DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];

					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain " + Array.toStr(DATA_FIELDS[3], "/") + " and " + Array.toStr(DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[4] == -1 && (genotypeIndices[0] == -1 || genotypeIndices[1] == -1)) {
						System.err.println("Error - File format not consistent! The files need to contain " + Array.toStr(GENOTYPE_FIELDS[0], "/") + " and " + Array.toStr(GENOTYPE_FIELDS[1], "/") + " or " + GENOTYPE_FIELDS[4] + " (for single token calls)");
						return;
					}
					if (genotypeIndices[5] == -1 && (genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					done = false;
					sampleName = "just starting";
					data = null;
					genotypes = null;
//					parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
					parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
					while (!done) {
						if (reader.ready()) {
							line = reader.readLine().split(delimiter);
							if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
								System.err.println("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample");
								parseAtAt = false;
							}
							trav = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
						} else {
							done = true;
							trav = null;
						}

						if (done || !trav.equals(sampleName)) {
							if (!sampleName.equals("just starting")) {
								if (fixes.containsKey(sampleName)) {
									sampleName = fixes.get(sampleName);
								}

								filename = determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true), sampleName, timeBegan);
								if (filename == null) {
									return;
								}

								samp = new Sample(sampleName, fingerprint, data, genotypes, true);
								samp.saveToRandomAccessFile(filename, allOutliers, sampleName);
							}
							if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + trav + Sample.SAMPLE_FILE_EXTENSION).exists()) {
								samp = Sample.loadFromRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(true, true) + (fixes.containsKey(trav) ? fixes.get(trav) : trav) + Sample.SAMPLE_FILE_EXTENSION, proj.JAR_STATUS.getValue());
								data = samp.getAllData();
								genotypes = samp.getAllGenotypes();
							} else {
								data = new float[DATA_FIELDS.length][];
								for (int j = 0; j < data.length; j++) {
									if (dataIndices[j] != -1) {
										data[j] = Array.floatArray(markerNames.length, Float.POSITIVE_INFINITY);
									}
								}
								genotypes = new byte[2][];
								genotypes[0] = Array.byteArray(markerNames.length, (byte) 0); // used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
								if (!ignoreAB) {
									genotypes[1] = Array.byteArray(markerNames.length, (byte) -1); // used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
								}
							}
							sampleName = trav;
						}

						if (!done) {
							countHash.add(line[snpIndex]);
							if (markerIndices.containsKey(line[snpIndex])) {
								key = markerIndices.get(line[snpIndex]);

								for (int j = 0; j < DATA_FIELDS.length; j++) {
									try {
										if (dataIndices[j] != -1) {
											if (!(data[j][key] + "").equals("Infinity")) {
												dupHash.add(line[snpIndex]);
												System.err.println("Sample " + trav + " already has data for marker " + line[snpIndex] + " (Was the parsing restarted? Delete the old directories first)");
											}
											data[j][key] = Float.parseFloat(line[dataIndices[j]]);
										}
									} catch (NumberFormatException nfe) {
										System.err.println("Error - failed to parse" + line[dataIndices[j]] + " into a valid " + DATA_FIELDS[j]);
										return;
									}
								}

								genos = parseGenotypes(line, genotypeIndices, ignoreAB, abLookup, count, sampleName, markerNames[count], files[i]);
								genotypes[0][key] = genos[0];
								if (!ignoreAB) {
									genotypes[1][key] = genos[1];
								}
							}
						}
						count++;
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + files[i] + "\"");
					return;
				}
			}

			if (allOutliers.size() > 0) {
				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").exists()) {
					System.err.println("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
					System.exit(1);
				} else {
					SerializedFiles.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
				}
			}

			System.out.println(ext.getTime() + "\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}

		System.out.println("Parsed " + count + " sample(s)");
		SampleList.generateSampleList(proj).writeToTextFile(proj.PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");

		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + "ListOfMarkers.txt"));
			writer.println("Marker\tExpected\tTimesSeen\tTimesDuplicated");
			for (int j = 0; j < markerNames.length; j++) {
				writer.println(markerNames[j] + "\t1\t" + countHash.getCount(markerNames[j]) + "\t" + dupHash.getCount(markerNames[j]));
				countHash.remove(markerNames[j], false);
			}
			list = countHash.getValues();
			for (int j = 0; j < list.length; j++) {
				writer.println(list[j] + "\t0\t" + countHash.getCount(markerNames[j]) + "\t.");
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "ListOfMarkers.txt");
			e.printStackTrace();
		}

		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + OVERWRITE_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + HOLD_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + CANCEL_OPTION_FILE).delete();
	}

	public static void mapFilenamesToSamples(Project proj, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int sampIndex;
		String[] files;
		String idHeader, delimiter;

		delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
		idHeader = proj.getProperty(proj.ID_HEADER);
		System.out.println(ext.getTime() + "\tSearching for " + proj.getProperty(proj.SOURCE_FILENAME_EXTENSION) + " files in: " + proj.SOURCE_DIRECTORY.getValue(false, true));
		files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		System.out.println("\t\tFound " + files.length + " file" + (files.length == 1 ? "" : "s") + " to parse");

		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + filename));
			for (int i = 0; i < files.length; i++) {
				try {
					// reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[i]);
					do {
						line = reader.readLine().trim().split(delimiter);
					} while (reader.ready() && (line.length < 3 || ext.indexOfStr(idHeader, line) == -1));

					if (!reader.ready()) {
						System.err.println("Error - went through enitre file without finding a line containing the user-defined ID header: " + idHeader);
						return;
					}
					sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];

					line = reader.readLine().split(delimiter);
					writer.println(files[i] + "\t" + line[sampIndex] + "\t" + (line[sampIndex].indexOf("@") >= 0 ? line[sampIndex].split("@")[0] : line[sampIndex]));
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + files[i] + "\" not found in " + proj.SOURCE_DIRECTORY.getValue(false, true));
					writer.close();
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + files[i] + "\"");
					writer.close();
					return;
				}
			}
			System.out.println(ext.getTime());
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void combineChpFiles(Project proj) {
		BufferedReader reader;
		PrintWriter writer;
		String idHeader, delimiter;
		idHeader = proj.getProperty(proj.ID_HEADER);
		delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
		String[] line;
		// common folder output by apt-genotype, present in each subdirectory of Source
		String commonSubFolder = "/cc-chp";
		// check source directory

		if (!proj.SOURCE_DIRECTORY.getValue(false, true).equals("") && !new File(proj.SOURCE_DIRECTORY.getValue(false, true)).exists()) {
			System.err.println("Error - the Project source location is invalid: " + proj.SOURCE_DIRECTORY.getValue(false, true));
			return;
		}

		String[] dirList = Files.listDirectories(proj.SOURCE_DIRECTORY.getValue(false, true), false);
		String[] chunkFiles = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true) + dirList[0] + commonSubFolder, proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);

		int counts = 0;

		for (int j = 0; j < chunkFiles.length; j++) {
			writer = Files.getAppropriateWriter(proj.SOURCE_DIRECTORY.getValue(false, true) + chunkFiles[j]);
			System.out.println("merging files " + (j + 1) + " of " + chunkFiles.length);

			for (int i = 0; i < dirList.length; i++) {
				try {
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + dirList[i] + commonSubFolder + "/" + chunkFiles[j]);
					// filter comments
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
					// if its the first directory, print the header

					if (i == 0) {
						writer.println(Array.toStr(line));
					}

					while (reader.ready()) {
						line = reader.readLine().trim().split(delimiter, -1);
						writer.println(Array.toStr(line));
						counts++;
					}
					System.out.println(counts);
					counts = 0;
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + chunkFiles[j] + "\" not found in " + proj.SOURCE_DIRECTORY.getValue(false, true) + dirList[i]);
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + proj.SOURCE_DIRECTORY.getValue(false, true) + dirList[i] + chunkFiles[j] + "\"");
					return;
				}
			}
			writer.close();
		}
	}

	// //iterate over directories in SOURCE_Directory

	// =
	// if (chunkFiles.length == 0) {
	// System.err.println("Error - no files to parse in "+srcDir+dirList[i]);
	// return;
	// }
	// //iterate over files in directories
	// for (int j =0; j< chunkFiles.length; j++){
	// try{
	// reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+dirList[i]+commonSubFolder+"/"+chunkFiles[j]);
	// //filter comments
	// do {
	// line = reader.readLine().trim().split(delimiter, -1);
	// } while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line)==-1)));
	// //if its the first directory, print the header
	// writer = Files.getAppropriateWriter(proj.getDir(proj.SOURCE_DIRECTORY)+chunkFiles[j]);
	// if(i ==0){
	// //System.out.println(Array.toStr(line));
	// numFiles = chunkFiles.length;
	// //System.out.println(proj.getDir(proj.SOURCE_DIRECTORY)+dirList[i]+"/"+chunkFiles[j]);
	//
	// writer.write(Array.toStr(line));
	// writer.write("\n");
	//
	// }
	// if(numFiles != chunkFiles.length){
	// System.err.println("There should be " +numFiles+ " files in each sub directory, there were only " + chunkFiles.length +" in " +proj.getDir(proj.SOURCE_DIRECTORY)+dirList[i]+commonSubFolder);
	// return;
	// }
	// //print lines following header from all file chunks
	// while (reader.ready()){
	// String aline = reader.readLine().trim();
	// writer.write(aline);
	// writer.write("\n");
	// }
	// reader.close();
	// writer.close();
	// }
	// catch (FileNotFoundException fnfe) {
	// System.err.println("Error: file \""+chunkFiles[j]+"\" not found in " +srcDir+dirList[i]);
	// return;
	// }
	// catch (IOException ioe) {
	// System.err.println("Error reading file \""+srcDir+dirList[i]+chunkFiles[j]+"\"");
	// return;
	// }
	// }
	// }
	//
	// }
	public static void parseAlleleLookupFromFinalReports(Project proj) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, String[]> hash;
		String[] files;
		int[] indices;
		int snpIndex;
		String delimiter, idHeader;
		String[] alleles;
		int expIndex;

		files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		if (files.length == 0) {
			System.err.println("Error - no files to parse");
			return;
		}
		System.out.println("\t\tFound " + files.length + " file" + (files.length == 1 ? "" : "s") + " to parse");

		idHeader = proj.getProperty(proj.ID_HEADER);
		delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
		hash = new Hashtable<String, String[]>();
		for (int i = 0; i < files.length; i++) {
			if (new File("report").exists()) {
				writeToLookupFile(proj, hash, i + 1);
				new File("report").delete();
			}
			try {
				System.out.println(ext.getTime() + "\t" + (i + 1) + " of " + files.length);
				// reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]));
				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[i]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));

				snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
				indices = ext.indexFactors(Sample.ALL_STANDARD_GENOTYPE_FIELDS, line, false, proj.getLog(), false, false);

				while (reader.ready()) {
					line = reader.readLine().split(delimiter);
					if (hash.containsKey(line[snpIndex])) {
						alleles = hash.get(line[snpIndex]);
					} else {
						if (i != 0) {
							System.err.println("Error - snp '" + line[snpIndex] + "' first seen in file #" + i + " (" + files[i] + ") and not earlier");
						}
						hash.put(line[snpIndex], alleles = new String[Sample.ALL_STANDARD_GENOTYPE_FIELDS.length]);
					}
					for (int j = 0; j < 2; j++) {
						if (ext.indexOfStr(line[indices[j]], Sample.ALT_NULL) == -1) {
							if (alleles[0] == null) {
								alleles[0] = line[indices[j]];
								expIndex = -1;
							} else if (line[indices[j]].equals(alleles[0])) {
								expIndex = 0;
							} else if (alleles[1] == null) {
								alleles[1] = line[indices[j]];
								expIndex = -2;
							} else if (line[indices[j]].equals(alleles[1])) {
								expIndex = 1;
							} else {
								System.err.println("Error - snp '" + line[snpIndex] + "' has a new allele in file #" + i + " (" + files[i] + "): " + line[indices[j]] + " (previously " + alleles[0] + "/" + alleles[1] + ")");
								expIndex = -9;
							}

							for (int k = 1; k < alleles.length / 2; k++) {
								switch (expIndex) {
								case -1:
									if (alleles[k * 2 + 0] != null) {
										System.err.println("Error - how can this be? -1");
									}
									alleles[k * 2 + 0] = line[indices[k * 2 + j]];
									break;
								case 0:
									if (!line[indices[k * 2 + j]].equals(alleles[k * 2 + 0])) {
										System.err.println("Error - how can this be? 0");
									}
									break;
								case -2:
									if (alleles[k * 2 + 1] != null) {
										System.err.println("Error - how can this be? -2");
									}
									alleles[k * 2 + 1] = line[indices[k * 2 + j]];
									break;
								case 1:
									if (!line[indices[k * 2 + j]].equals(alleles[k * 2 + 1])) {
										System.err.println("Error - how can this be? 1");
									}
									break;
								case -9:
									break;
								}
							}
						}
					}
				}
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + files[i] + "\"");
				return;
			}
		}
		writeToLookupFile(proj, hash, 0);
	}

	public static void writeToLookupFile(Project proj, Hashtable<String, String[]> hash, int fileNumber) {
		PrintWriter writer;
		String[] keys;

		System.out.print("Writing to file...");
		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + "alleleLookup" + (fileNumber > 0 ? "_atFile" + fileNumber : "") + ".xln"));
			keys = HashVec.getKeys(hash, false, false);
			writer.println("SNP\t" + Array.toStr(Sample.ALL_STANDARD_GENOTYPE_FIELDS));
			for (int k = 0; k < keys.length; k++) {
				writer.println(keys[k] + "\t" + Array.toStr(hash.get(keys[k])));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + "alleleLookup.xln");
			e.printStackTrace();
		}
		System.out.println("done");

	}

	// public static Hashtable filenamesUsedByExistingSampleFiless(Project proj) {
	// Hashtable <String, String> result = null;
	// String[] sampNames;
	// File[] files;
	// String[] fileNames;
	//
	// files = new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)).listFiles(new FilenameFilter() {
	// public boolean accept(File file, String filename) {
	// return filename.endsWith(Sample.SAMPLE_DATA_FILE_EXTENSION);
	// }
	// });
	//
	// if (files != null) {
	// sampNames = proj.getSamples();
	// fileNames = new String[files.length];
	// result = new Hashtable<String, String>();
	// for (int i = 0; i < files.length; i++) {
	// fileNames[i] = files[i].getName();
	// }
	// for (int i = 0; i < sampNames.length; i++) {
	// for (int j = 0; j < fileNames.length; j++) {
	// if (sampNames[i].equals(fileNames[j])) {
	// result.put(sampNames[i], files[j].getName());
	// }
	// }
	// }
	// }
	//
	// return result;
	// }
	//
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = null;
		boolean map = false;
		int numThreads = 1;
		// boolean parseABlookup = false;
		boolean parseAlleleLookupFromFinalReports = false;
		boolean combineChpFiles = false;
		String mapOutput = "filenamesMappedToSamples.txt";

		String usage = "\n" + 
		"cnv.manage.ParseKcol requires 0-1 arguments\n" + 
		"   (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) number of threads to use (i.e. threads=" + numThreads + " (default))\n" + 
		" OPTIONAL:\n" + 
		"   (3) map filenames to sample IDs (i.e. -mapFiles (" + (map ? "" : "not the ") + "default))\n" + 
		"   (4) output file for mappings (i.e. out=" + mapOutput + " (default))\n" +
		// " OR:\n"+
		// "   (1) parse AB lookup (i.e. --parseAB (not the default))\n"+
		" OR:\n" + 
		"   (1) parse Forward/TOP/AB/etc lookup (i.e. --parseAlleleLookup (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-mapFiles")) {
				map = true;
				numArgs--;
				// } else if (args[i].startsWith("-parseAB")) {
				// parseABlookup = true;
				// numArgs--;
			} else if (args[i].startsWith("-parseAlleleLookup")) {
				parseAlleleLookupFromFinalReports = true;
				numArgs--;
			} else if (args[i].startsWith("-combine")) {
				combineChpFiles = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}

		// filename = "C:/workspace/Genvisis/projects/KcolTest.properties";
		// proj = null;
		proj = new Project(filename, false);
		//
		// if (!proj.getDir(proj.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(proj.SOURCE_DIRECTORY)).exists()) {
		// System.err.println("Error - the project source location is invalid: "+proj.getDir(proj.SOURCE_DIRECTORY));
		// return;
		// }

		// fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.4pc_genotypes.txt";
		// fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.4pc_full_apoe_genotypes.txt";
		// fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.full_apoe_genotypes.txt";
		// fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.nocovars_genotypes.txt";

		// lookupFile = "D:\\LOAD\\alleleLookup.txt";

		try {
			// if (parseABlookup) {
			// ABLookup.parseABlookup(proj);
			// } else
			if (map) {
				mapFilenamesToSamples(proj, mapOutput);
			} else if (parseAlleleLookupFromFinalReports) {
				parseAlleleLookupFromFinalReports(proj);
			} else if (combineChpFiles) {
				combineChpFiles(proj);
			} else {
				// combineChpFiles(proj);
				createFiles(proj, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
