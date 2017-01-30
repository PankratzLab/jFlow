package org.genvisis.kaput;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class FindNearestGenes {
	public static final String DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\merged\\results\\";
	public static final String MAP = "gwas.map";
	public static final String GENES = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\CommonTools\\genes.xls";
	public static final String DEFAULT_HIT_LIST = "marks2.txt";
	public static final int WINDOW = 3;

	public static void findNearest(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, travLine;
		Hashtable<String, String> poslar;
		String[] markers, genes;
		int[][] genePositions;
		int[] chr_start_stop, trav;
		Vector<int[]> queuePoslar;
		Vector<String> queueNames;
		String nearestMarker;
		int nearestPosition;

		System.out.println(ext.getTime());
		System.out.println("Loading markers...");
		markers = ArrayUtils.toStringArray(HashVec.loadFileToVec(DIR + filename, false, true, false));

		System.out.println("Loading map...");
		poslar = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(DIR + MAP));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (ext.indexOfStr(line[1], markers) >= 0) {
					poslar.put(line[1], line[0] + "\t" + line[3]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + DIR + MAP + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + DIR + MAP + "\"");
			System.exit(2);
		}

		System.out.println("Searching through genes...");
		genePositions = new int[markers.length][2];
		genes = ArrayUtils.stringArray(markers.length, "");
		for (int j = 0; j < markers.length; j++) {
			if (poslar.containsKey(markers[j])) {
				line = (poslar.get(markers[j])).split("[\\s]+");
				for (int k = 0; k < 2; k++) {
					genePositions[j][k] = Integer.parseInt(line[k]);
				}
			} else {
				System.err.println("Error - "	+ markers[j] + " was listed in " + filename
														+ " but was not found in the map loaded (" + MAP + ")");
				System.exit(1);
			}
		}
		try {
			reader = new BufferedReader(new FileReader(GENES));
			reader.readLine();
			chr_start_stop = new int[3];
			queuePoslar = new Vector<int[]>();
			queueNames = new Vector<String>();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				chr_start_stop = parseChrInfo(line);

				if (queuePoslar.size() > 0
						&& queuePoslar.elementAt(queuePoslar.size() - 1)[0] != chr_start_stop[0]) {
					queuePoslar.clear();
					queueNames.clear();
				}

				for (int j = 0; j < markers.length; j++) {
					if (genePositions[j][0] == chr_start_stop[0]	&& genePositions[j][1] < chr_start_stop[1]
							&& genes[j].equals("")) {
						if (queuePoslar.size() == 0) {
							genes[j] = "xx";
							nearestMarker = "xx";
							nearestPosition = -99999999;
						} else {
							nearestMarker = queueNames.elementAt(queuePoslar.size() - 1);
							nearestPosition =
															(genePositions[j][1] < queuePoslar.elementAt(queuePoslar.size()
																																						- 1)[2]	? 0
																																										: -1
																																											* (genePositions[j][1]
																																													- queuePoslar.elementAt(queuePoslar.size()
																																																									- 1)[2]));
							for (int i = 0; i < queuePoslar.size(); i++) {
								genes[j] += (genes[j].equals("") ? "" : "|")	+ queueNames.elementAt(i) + "["
														+ (queuePoslar.elementAt(i)[2] - queuePoslar.elementAt(i)[1]) + "] ("
														+ (genePositions[j][1] < queuePoslar.elementAt(i)[2]	? "within"
																																									: "-"
																																										+ (genePositions[j][1]
																																												- queuePoslar.elementAt(i)[2]))
														+ ")";
								if (genePositions[j][1] - queuePoslar.elementAt(i)[2] < Math.abs(nearestPosition)) {
									nearestMarker = queueNames.elementAt(i);
									nearestPosition =
																	(genePositions[j][1] < queuePoslar.elementAt(i)[2]	? 0
																																											: -1
																																												* (genePositions[j][1]
																																														- queuePoslar.elementAt(i)[2]));
								}

							}
						}
						reader.mark(1000);
						for (int i = 0; i < WINDOW; i++) {
							travLine = i == 0 ? line : reader.readLine().trim().split("[\\s]+");
							trav = parseChrInfo(travLine);
							if (genePositions[j][0] == trav[0]) {
								genes[j] += "|"	+ travLine[1] + "[" + (trav[2] - trav[1]) + "] (+"
														+ (trav[1] - genePositions[j][1]) + ")";
							}
							if (i == 0 && trav[1] - genePositions[j][1] < Math.abs(nearestPosition)) {
								nearestMarker = line[1];
								nearestPosition = trav[1] - genePositions[j][1];
							}
						}
						reader.reset();

						genes[j] = nearestMarker	+ "("
												+ (nearestPosition == 0	? "within"
																								: (nearestPosition < 0 ? "-" : "+")
																									+ ext.prettyUpDistance(	Math.abs(nearestPosition),
																																					0))
												+ ")\t" + genes[j];
					}
				}

				queuePoslar.add(chr_start_stop);
				queueNames.add(line[1]);
				if (queuePoslar.size() > WINDOW) {
					queuePoslar.removeElementAt(0);
					queueNames.removeElementAt(0);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + GENES + "\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + GENES + "\"");
			System.exit(2);
		}

		System.out.println("Writing to files...");
		try {
			writer = new PrintWriter(new FileWriter(DIR	+ filename.substring(0, filename.indexOf("."))
																							+ "_described.out"));
			writer.println("SNP\tChr\tloc_36.1\tGene");
			for (int j = 0; j < markers.length; j++) {
				writer.println(markers[j] + "\t" + poslar.get(markers[j]) + "\t" + genes[j]);
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing to file \"" + DIR + filename + "_described.out" + "\"");
			System.exit(2);
		}
		System.out.println(ext.getTime());

	}

	public static int[] parseChrInfo(String[] line) {
		int[] chr_start_stop = new int[3];

		chr_start_stop[0] =
											line[2].equals("X")	? 23
																					: (line[2].equals("Y")	? 24
																																	: (line[2].equals("XY")	? 25
																																													: (line[2].equals("MT")	? 26
																																																									: (line[2].equals("Un")	? 27
																																																																					: Integer.parseInt(line[2])))));
		chr_start_stop[1] = Integer.parseInt(line[3]);
		chr_start_stop[2] = Integer.parseInt(line[4]);

		return chr_start_stop;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = DEFAULT_HIT_LIST;

		String usage = "\\n"	+ "gwas.FindNearestGenes requires 0-1 arguments\n"
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
			findNearest(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
