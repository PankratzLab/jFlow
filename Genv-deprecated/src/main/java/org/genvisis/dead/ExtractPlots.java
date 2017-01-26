// -Xms1024M -Xmx1024M
package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class ExtractPlots {
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\ADNI\\genotypes\\";
	public static final String ORIGINAL_LINUX_DIRECTORY = "/net/collabs/genetics_anals/LONI_subjs";
	public static final String LINUX_DIRECTORY = "/net/collabs/genetics_anals/LONI_subjs_fixed/";
	public static final String SLIM = "slim/";
	public static final String LOOKUP_DIR = "lookup/";
	public static final String PLOT_HEADER = "Sample\tGC Score\tTheta\tR\tX\tY\tX Raw\tY Raw\tB Allele Freq\tLog R Ratio\tAlleleCalls\tAlleleCount";
	public static final String MARKERS_IN_INDEX_ORDER = "markersInIndexOrder.dat";
	public static final String[] PLOT_NEEDS = {	"GC Score", "Theta", "R", "X", "Y", "X Raw", "Y Raw",
																							"B Allele Freq", "Log R Ratio", "Allele1 - Forward",
																							"Allele2 - Forward", "Allele1 - AB", "Allele2 - AB"};

	public static void extractPlotsInMemoryAndIndexed(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String trav, temp;
		Vector<String> markers = new Vector<String>();
		Vector<String> loaded = new Vector<String>();
		int count, version, countAll;
		String[] ids;
		int numMarkers, numSamples;
		int[] markerNameIndices, markerLocationIndices;
		String[][] data;
		int[] indices;
		boolean pdgwas;
		int snpNameIndex, sampleNameIndex;

		markers = HashVec.loadFileToVec(filename, false, false, true);
		numMarkers = markers.size();

		if (new File(SLIM).exists()) {
			trav = SLIM;
		} else if (new File(WINDOWS_DIRECTORY).exists()) {
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
				return filename.endsWith(".csv");
			}
		});
		numSamples = files.length;
		ids = ArrayUtils.stringArray(numSamples);

		markerNameIndices = ArrayUtils.intArray(numMarkers, -1);
		markerLocationIndices = ArrayUtils.intArray(numMarkers, -1);
		try {
			reader = new BufferedReader(new FileReader(files[0]));
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.contains("SNP Name"));
			snpNameIndex = ext.indexOfStr("SNP Name", temp.trim().split(","));
			count = countAll = 0;
			while (reader.ready()) {
				line = reader.readLine().split(",");
				if (markers.contains(line[snpNameIndex])) {
					markerNameIndices[count] = markers.indexOf(line[snpNameIndex]);
					markerLocationIndices[count] = countAll;
					count++;
				}
				countAll++;
			}
			reader.close();

			if (ArrayUtils.min(markerNameIndices) == -1) {
				for (int markerNameIndice : markerNameIndices) {
					if (markerNameIndice != -1) {
						markers.removeElementAt(markerNameIndice);
						markers.insertElementAt("", markerNameIndice);
					}
				}
				for (int i = 0; i < markers.size(); i++) {
					if (!markers.elementAt(i).equals("")) {
						System.err.println("Error - marker '"	+ markers.elementAt(i)
																+ "' was not found in the index file: " + files[0]);
						writer = new PrintWriter(new FileWriter("../MAJOR_PROBLEMS_WITH_INDICES!!.out", true));
						writer.println("Error - marker '"	+ markers.elementAt(i)
														+ "' was not found in the index file: " + files[0]);
						writer.close();

					}
				}
				System.exit(1);
			}

		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + files[0] + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + files[0] + "\"");
			System.exit(2);
		}

		System.out.println(ext.getTime()	+ " " + ext.getDate() + " (all in memory, but indexed; "
												+ filename + " had " + markers.size() + " markers)");
		Date date = new Date();
		data = new String[numSamples][numMarkers];
		try {
			for (int i = 0; i < numSamples; i++) {
				try {
					// System.out.println(files[i].getName());
					reader = new BufferedReader(new FileReader(files[i]));
					pdgwas = files[i].getName().contains("Myers");
					do {
						temp = reader.readLine();
					} while (reader.ready() && !temp.contains("SNP Name"));
					line = temp.split(",");
					indices = ext.indexFactors(PLOT_NEEDS, line, false, true);
					snpNameIndex = ext.indexOfStr("SNP Name", line);
					sampleNameIndex = ext.indexOfStr("Sample Name", line);

					count = 0;
					for (int j = 0; j < numMarkers; j++) {
						while (count < markerLocationIndices[j]) {
							reader.readLine();
							count++;
						}

						line = reader.readLine().split(",");
						try {
							if (!line[snpNameIndex].equals(markers.elementAt(markerNameIndices[j]))) {
								System.err.println("Error - out of sync in file "	+ files[i].getName()
																		+ ": expecting " + markers.elementAt(markerNameIndices[j])
																		+ ", found " + line[snpNameIndex]);
								System.exit(1);
							}
						} catch (Exception e) {
							System.err.println("Error - ");
							e.printStackTrace();
						}
						count++;
						if (ids[i].equals("")) {
							ids[i] = pdgwas	? line[sampleNameIndex].substring(0,
																																line[sampleNameIndex].indexOf("@"))
															: line[sampleNameIndex];
							version = 0;
							if (loaded.contains(ids[i] + (version == 0 ? "" : "__" + version))) {
								version++;
							}
							ids[i] += (version == 0 ? "" : "__" + version);
						}
						data[i][j] = line[indices[0]]	+ "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t"
													+ line[indices[3]] + "\t" + line[indices[4]] + "\t" + line[indices[5]]
													+ "\t" + line[indices[6]] + "\t" + line[indices[7]] + "\t"
													+ line[indices[8]] + "\t" + line[indices[9]] + line[indices[10]] + "\t"
													+ (line[indices[11]].equals("-")	? "-1"
																														: (line[indices[12]].equals("B")	? (line[indices[11]].equals("B")	? "2"
																																																																: "1")
																																															: "0"));
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ files[i].getName()
															+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + files[i].getName() + "\"");
					System.exit(2);
				}
			}
			for (int i = 0; i < numMarkers; i++) {
				writer = new PrintWriter(new FileWriter(translateName(markers.elementAt(markerNameIndices[i]))
																								+ "_plots.xls"));
				writer.println(PLOT_HEADER);
				for (int j = 0; j < numSamples; j++) {
					writer.println(ids[j] + "\t" + data[j][i]);
				}
				writer.close();
			}
			System.out.println(ext.getTime()	+ " " + ext.getDate() + " (finished in "
													+ ext.getTimeElapsed(date.getTime()));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void slimDownFiles() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav, temp;
		int numSamples;
		int[] indices;
		int snpNameIndex, sampleNameIndex;

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
				return filename.endsWith(".csv");
			}
		});
		numSamples = files.length;

		System.out.println(ext.getTime() + " " + ext.getDate() + " Slimming down all files");
		Date date = new Date();
		new File(SLIM).mkdirs();
		try {
			for (int i = 0; i < numSamples; i++) {
				try {
					System.out.println(files[i].getName());
					reader = new BufferedReader(new FileReader(files[i]));
					do {
						temp = reader.readLine();
					} while (reader.ready() && !temp.contains("SNP Name"));
					line = temp.split(",");
					indices = ext.indexFactors(PLOT_NEEDS, line, false, true);
					snpNameIndex = ext.indexOfStr("SNP Name", line);
					sampleNameIndex = ext.indexOfStr("Sample Name", line);

					writer = new PrintWriter(new FileWriter(SLIM + files[i].getName()));
					writer.println(line[sampleNameIndex]	+ "," + line[snpNameIndex] + ","
													+ ArrayUtils.toStr(PLOT_NEEDS, ","));
					while (reader.ready()) {
						line = reader.readLine().split(",");
						writer.print(line[sampleNameIndex] + "," + line[snpNameIndex]);
						for (int indice : indices) {
							writer.print("," + line[indice]);
						}
						writer.println();
					}
					reader.close();
					writer.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ files[i].getName()
															+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + files[i].getName() + "\"");
					System.exit(2);
				}
			}
			System.out.println(ext.getTime()	+ " " + ext.getDate() + " (finished in "
													+ ext.getTimeElapsed(date.getTime()));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void fixADNIformat() {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, header;
		String trav, temp;
		String[] markers;
		int count, numSamples;
		long time;

		markers =
						ArrayUtils.toStringArray(HashVec.loadFileToVec(MARKERS_IN_INDEX_ORDER, false, false, false));

		if (new File(WINDOWS_DIRECTORY).exists()) {
			trav = WINDOWS_DIRECTORY;
		} else if (new File(ORIGINAL_LINUX_DIRECTORY).exists()) {
			trav = ORIGINAL_LINUX_DIRECTORY;
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
				return filename.endsWith(".csv");
			}
		});
		numSamples = files.length;

		System.out.println(ext.getTime()	+ " " + ext.getDate() + " " + MARKERS_IN_INDEX_ORDER + " had "
												+ markers.length + " markers)");
		time = new Date().getTime();
		new File("fixed/").mkdirs();
		try {
			for (int i = 0; i < numSamples; i++) {
				try {
					if (!new File("fixed/" + files[i].getName()).exists()) {
						System.out.println(files[i].getName());
						reader = new BufferedReader(new FileReader(files[i]));
						writer = new PrintWriter(new FileWriter("fixed/" + files[i].getName()));

						do {
							temp = reader.readLine();
							writer.println(temp);
						} while (reader.ready() && !temp.contains("SNP Name"));
						header = temp.trim().split(",", -1);

						count = 0;
						while (reader.ready()) {
							line = reader.readLine().split(",", -1);
							if (line.length == header.length) {
								writer.println(ArrayUtils.toStr(line, ","));
							} else if (line.length == header.length * 2 - 1) {
								writer.println(ArrayUtils.toStr(ArrayUtils.subArray(line, 0, header.length - 1), ","));
								count++;
								writer.println(ArrayUtils.toStr(ArrayUtils.subArray(line, 0, 3), ",")	+ ","
																+ ArrayUtils.toStr(ArrayUtils.subArray(line, header.length - 1 + 3), ","));
							} else {
								System.out.println(count	+ " Expecting either " + header.length + " or "
																		+ (header.length * 2 - 1) + "; found: " + line.length);
							}
							count++;
						}
						if (count != markers.length) {
							System.err.println("Error - only recovered " + count + " lines ");
						}
						reader.close();
						writer.close();
					}
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ files[i].getName()
															+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + files[i].getName() + "\"");
					System.exit(2);
				}
			}
			System.out.println(ext.getTime()	+ " " + ext.getDate() + " (finished in "
													+ ext.getTimeElapsed(time) + ")");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void batchExtraction(String prefix, int reps, int repStep) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		int count, rep;
		boolean alreadyDone;

		new File(LOOKUP_DIR).mkdirs();
		File[] files = new File(LOOKUP_DIR).listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return filename.endsWith(".txt");
			}
		});
		Vector<Vector<String>> onesPicked = HashVec.newVecVecString(reps);

		for (int i = 0; i < files.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(files[i]));
				while (reader.ready()) {
					line = reader.readLine().split("[\\s]+");
					hash.put(line[0], i + "");
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + files[i] + "\"");
				System.exit(2);
			}
		}

		try {
			reader = new BufferedReader(new FileReader("../" + MARKERS_IN_INDEX_ORDER));
			writer = new PrintWriter(new FileWriter(LOOKUP_DIR + "markersDone.out"));
			rep = count = 0;
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				alreadyDone = false;
				for (int i = 0; i < files.length && !alreadyDone; i++) {
					if (hash.containsKey(line[0])) {
						alreadyDone = true;
					}
				}
				writer.println(line[0] + "\t" + (alreadyDone ? "1" : ""));
				if (rep < reps && !alreadyDone) {
					onesPicked.elementAt(rep).add(line[0]);
					count++;
					if (count == repStep) {
						rep++;
						count = 0;
					}
				}

			}
			writer.close();
			reader.close();

			for (int i = 0; i < reps; i++) {
				if (onesPicked.elementAt(i).size() > 0) {
					writer = new PrintWriter(new FileWriter(LOOKUP_DIR + prefix + (char) (i + 97) + ".txt"));
					for (int j = 0; j < onesPicked.elementAt(i).size(); j++) {
						writer.println(onesPicked.elementAt(i).elementAt(j));
					}
					writer.close();
				}
			}

			writer = new PrintWriter(new FileWriter("batch" + prefix));
			writer.println("mkdir " + prefix);
			writer.println("cp " + LOOKUP_DIR + prefix + "*.txt " + prefix + "/");
			writer.println("cd " + prefix);
			for (int i = 0; i < reps; i++) {
				// writer.println("jcp gwas.ExtractPlots
				// file="+prefix+(char)(i+97)+".txt");
				writer.println("java -jar /home/genanal/"	+ org.genvisis.common.PSF.Java.GENVISIS
												+ " gwas.ExtractPlots file=" + prefix + (char) (i + 97) + ".txt");

			}
			writer.println("cd ..");
			writer.close();

		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ "../" + MARKERS_IN_INDEX_ORDER
													+ "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + "../" + MARKERS_IN_INDEX_ORDER + "\"");
			System.exit(2);
		}
	}

	public static String translateName(String str) {
		return ext.replaceAllWith(str, ":", ".");
	}

	public static void batchAll() {
		BufferedReader reader;
		PrintWriter writer;
		int count;

		count = 0;
		do {
			count++;
			batchExtraction(count + "", 20, 2000);
		} while (new File(LOOKUP_DIR + count + "t.txt").exists());

		try {
			writer = new PrintWriter(new FileWriter("batchAll"));
			for (int i = 1; i <= count; i++) {
				try {
					reader = new BufferedReader(new FileReader("batch" + i));
					while (reader.ready()) {
						writer.println(reader.readLine());
					}
					reader.close();
					new File("batch" + i).delete();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + "batch" + i + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + "batch" + i + "\"");
					System.exit(2);
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to batchAll");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "1a.txt";
		boolean fixadni = false;
		boolean slim = false;
		boolean batchAll = false;

		String usage = "\\n"	+ "park.gwa.ExtractPlots requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n"
										+ "   (2) fix ADNI files (i.e. -fixadni (not the default))\n"
										+ "   (3) slim down files (i.e. -slim (not the default))\n"
										+ "   (4) batch extraction of all (i.e. -batch (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-fixadni")) {
				fixadni = true;
				numArgs--;
			} else if (arg.startsWith("-slim")) {
				slim = true;
				numArgs--;
			} else if (arg.startsWith("-batch")) {
				batchAll = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (fixadni) {
				fixADNIformat();
			} else if (slim) {
				slimDownFiles();
			} else if (batchAll) {
				batchAll();
			} else {
				extractPlotsInMemoryAndIndexed(filename);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
