package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CountVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.PSF;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class LinkageMap {
	public static final double DEFAULT_DX_ALLELE_FREQ = 0.01;
	public static final double[][] DEFAULT_AUTOSOMAL_DOMINANT_PENTRANCE = {{0.03, 0.80, 0.80}};
	public static final double[][] DEFAULT_AUTOSOMAL_RECESSIVE_PENTRANCE = {{0.03, 0.03, 0.80}};
	public static final double[][] DEFAULT_X_LINKED_DOMINANT_PENTRANCE = {{0.03, 0.80, 0.80},
																																				{0.03, 0.80}};
	public static final double[][] DEFAULT_X_LINKED_RECESSIVE_PENTRANCE = {{0.03, 0.03, 0.80},
																																				 {0.03, 0.80}};
	public static final double MIN_DIST_STEP = 0.0001;

	private int chr;
	private double dxAlleleFreq;
	private double[][] penetranceFunctions;
	private String[] markerNames;
	private String[] addInfo;
	private double[][] alleleFreqs;
	private double[] distances;
	private double[] positions;
	private final String filename;
	private boolean currentlyInMorgans;

	public LinkageMap(int chr, String[] markerNames, int numDefaultAlleles, double[] location_data,
										boolean currentlyInMorgans, boolean positionsInsteadOfDistances) {
		this(chr, markerNames, createDummyAlleleFreqs(markerNames.length, numDefaultAlleles),
				 location_data, currentlyInMorgans, positionsInsteadOfDistances);
	}

	public LinkageMap(int chr, String[] markerNames, double[][] alleleFreqs, double[] location_data,
										boolean currentlyInMorgans, boolean positionsInsteadOfDistances) {
		this.chr = chr;
		this.markerNames = markerNames;
		addInfo = new String[markerNames.length];
		this.alleleFreqs = alleleFreqs;
		dxAlleleFreq = DEFAULT_DX_ALLELE_FREQ;
		if (chr == 23) {
			penetranceFunctions = DEFAULT_X_LINKED_DOMINANT_PENTRANCE;
		} else {
			penetranceFunctions = DEFAULT_AUTOSOMAL_DOMINANT_PENTRANCE;
		}
		filename = null;

		this.currentlyInMorgans = currentlyInMorgans;
		if (positionsInsteadOfDistances) {
			positions = location_data;
			distances = computeDistancesFromPositions(location_data);
		} else {
			distances = location_data;
			positions = computeCumulativePositions(location_data);
		}

		if (markerNames.length != positions.length) {
			System.err.println("Error - " + markerNames.length + " marker names and " + positions.length
												 + " distance measures");
		}
	}

	public LinkageMap(int chr) {
		this("map" + ext.chrome(chr) + ".dat");
	}

	public LinkageMap(String dir, int chr) {
		this(dir + "map" + ext.chrome(chr) + ".dat");
	}

	public LinkageMap(String filename) {
		BufferedReader reader;
		String[] line;
		String trav;
		boolean noNamesFlag = false;
		int expNumAllele;
		boolean problem = false;
		int numPen;

		this.filename = filename;
		if (ext.rootOf(filename).startsWith("map") && filename.endsWith(".dat")
				&& ext.rootOf(filename).length() == 5) {
			try {
				chr = Integer.parseInt(ext.rootOf(filename).substring(3, 5));
			} catch (NumberFormatException e) {
				System.err.println("Warning - could not parse chromosome number for " + filename);
				if (filename.contains("23")) {
					chr = 23;
					System.err.println("           assuming it's the X chromosome");
				} else {
					chr = 0;
					System.err.println("           assuming it's an autosome");
				}
			}
		}

		try {
			reader = new BufferedReader(new FileReader(filename));

			markerNames = new String[Integer.valueOf(reader.readLine().trim()
																										 .split(PSF.Regex.GREEDY_WHITESPACE)[0])
																			.intValue()
															 - 1];
			for (int i = 0; i < 3; i++) {
				reader.readLine();
			}
			try {
				dxAlleleFreq = Double.parseDouble(reader.readLine().trim()
																								.split(PSF.Regex.GREEDY_WHITESPACE)[1]);
			} catch (NumberFormatException e) {
				System.err.println("Error - failed to parse disease allele frequency in file: " + filename);
				System.exit(1);
			}

			numPen = Integer.valueOf(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE)[0])
											.intValue();
			penetranceFunctions = new double[numPen][];
			for (int i = 0; i < numPen; i++) {
				try {
					penetranceFunctions[i] = ArrayUtils.toDoubleArray(reader.readLine().trim()
																																	.split(PSF.Regex.GREEDY_WHITESPACE));
				} catch (NumberFormatException e) {
					System.err.println("Error - failed to parse the penetrance model for file: " + filename);
					System.exit(1);
				}
			}
			if (chr == 23 && numPen == 1) {
				System.err.println("Error - X chromosome only has one penetrance function");
				System.exit(1);
			}

			alleleFreqs = new double[markerNames.length][];
			addInfo = new String[markerNames.length];
			for (int i = 0; i < markerNames.length; i++) {
				trav = reader.readLine();
				line = trav.trim().split(PSF.Regex.GREEDY_WHITESPACE);
				if (line.length < 4) {
					if (!noNamesFlag) {
						System.err.println("Warning - no marker names for at least a subset of the markers in "
															 + filename);
						noNamesFlag = true;
					}
					markerNames[i] = "Marker" + (i + 1);
				} else {
					markerNames[i] = line[3];
					addInfo[i] = trav.substring(trav.indexOf(line[3]) + line[3].length());
				}
				expNumAllele = Integer.parseInt(line[1]);
				line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
				try {
					if (line.length == expNumAllele) {
						alleleFreqs[i] = ArrayUtils.toDoubleArray(line);
					} else {
						if (line.length < expNumAllele) {
							System.err.println("Error - the number of allele frequencies provided for "
																 + markerNames[i] + " (n=" + alleleFreqs.length
																 + ") does not match up with the declared number (" + expNumAllele
																 + ")");
							problem = true;
						} else if (line[Integer.parseInt(line[1])].startsWith("<")) {
							alleleFreqs[i] = ArrayUtils.toDoubleArray(ArrayUtils.subArray(line, 0, expNumAllele));
						} else {
							System.err.println("Error - the number of allele frequencies for " + markerNames[i]
																 + " (n=" + alleleFreqs.length
																 + ") does not match up with the declared number ("
																 + Integer.parseInt(line[1]) + ")");
							problem = true;
						}
					}
				} catch (NumberFormatException e) {
					System.err.println("Error - failed to parse the allele frequencies for " + markerNames[i]
														 + " (" + ArrayUtils.toStr(line) + ")");
					problem = true;
				}
			}
			if (problem) {
				System.exit(1);
			}

			reader.readLine();
			line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
			distances = new double[markerNames.length];
			for (int i = 0; i < markerNames.length; i++) {
				distances[i] = Double.parseDouble(line[i]);
			}
			currentlyInMorgans = distances[0] < 1;
			positions = computeCumulativePositions(distances);
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error parsing " + filename);
		}
	}

	public void setChr(int chr) {
		this.chr = chr;
	}

	public int getChr() {
		return chr;
	}

	public String[] getMarkerNames() {
		return markerNames;
	}

	public double[][] getAlleleFreqs() {
		return alleleFreqs;
	}

	public double[] getDistances(boolean inMorgans) {
		return getCumulative(distances, inMorgans);
	}

	public boolean isCurrentlyInMorgans() {
		return currentlyInMorgans;
	}

	public double[] getCumulativePositions(boolean inMorgans) {
		return getCumulative(positions, inMorgans);
	}

	private double[] getCumulative(double[] base, boolean inMorgans) {
		double[] cumulative = new double[base.length];
		for (int i = 0; i < cumulative.length; i++) {
			cumulative[i] = base[i];
			if (inMorgans && !isCurrentlyInMorgans()) {
				cumulative[i] /= 100;
			} else if (!inMorgans && isCurrentlyInMorgans()) {
				cumulative[i] *= 100;
			}
		}
		return cumulative;
	}

	public static double[] computeDistancesFromPositions(double[] positions) {
		double[] distances = new double[positions.length];

		distances[0] = 10.0;
		for (int i = 1; i < positions.length; i++) {
			distances[i] = positions[i] - positions[i - 1];
		}
		return distances;
	}

	public static double[] computeCumulativePositions(double[] distances) {
		double[] cumulative = new double[distances.length];
		double sum;

		sum = distances[0] * -1;
		for (int i = 0; i < distances.length; i++) {
			sum += distances[i];
			cumulative[i] = sum;
		}
		return cumulative;
	}

	public void updateFile() {
		if (filename == null) {
			System.err.println("Error - can't update map file; no filename was stored; try createFile() for standard filename");
			System.exit(1);
		}
		createFile(filename, currentlyInMorgans);
	}

	public void createFileInDir(String dir) {
		createFile(dir + "map" + ext.chrome(chr) + ".dat", currentlyInMorgans);
	}

	public void createFile() {
		createFile("map" + ext.chrome(chr) + ".dat", currentlyInMorgans);
	}

	public void createFile(String filename) {
		createFile(filename, currentlyInMorgans);
	}

	public void createFile(String filename, boolean writeInMorgans) {
		PrintWriter writer;

		try {
			writer = Files.openAppropriateWriter(filename);
			writer.println((markerNames.length + 1) + " 0 " + (chr == 23 ? 1 : 0)
										 + " 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM");
			writer.println("0 0.0 0.0 0  << MUT LOCUS, MUT RATE, HAPLOTYPE FREQUENCIES (IF 1)");
			writer.println(ArrayUtils.toStr(ArrayUtils.stringArraySequence(markerNames.length + 1, ""),
																			" "));
			writer.println("1  2  << AFFECTATION, NO. OF ALLELES");
			writer.println(ext.formDeci(1 - dxAlleleFreq, 2, 5) + " " + ext.formDeci(dxAlleleFreq, 2, 5)
										 + " << GENE FREQUENCIES");
			writer.println(penetranceFunctions.length + "  << NO. OF LIABILITY CLASSES");
			for (double[] penetranceFunction : penetranceFunctions) {
				writer.println(ArrayUtils.toStr(penetranceFunction, 2, 5, " "));
			}

			for (int i = 0; i < markerNames.length; i++) {
				writer.println("3 " + alleleFreqs[i].length + "  # " + markerNames[i]
											 + (addInfo[i] == null ? "" : addInfo[i]));
				writer.println(ArrayUtils.toStr(alleleFreqs[i], 6, 6, " "));
			}

			writer.println("0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)");
			writer.print("0.10");
			// System.out.println("Currently in
			// "+(currentlyInMorgans?"":"centi")+"Morgans");
			// System.out.println("Writing in
			// "+(writeInMorgans?"":"centi")+"Morgans");
			for (int i = 0; i < markerNames.length - 1; i++) {
				writer.print(" " + ext.formDeci(
																				writeInMorgans ? (currentlyInMorgans ? distances[i + 1]
																																						 : distances[i + 1]
																																							 / 100)
																											 : (currentlyInMorgans ? distances[i + 1]
																																							 * 100
																																						 : distances[i + 1]),
																				8, false));
			}
			writer.println("  << RECOMB VALUES");
			writer.println("1 0.1 0.45  << REC VARIED, INCREMENT, FINISHING VALUE");

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing " + filename);
			e.printStackTrace();
		}
	}

	public void createPlinkMap(String filename) {
		createPlinkMap(filename, null);
	}

	public void createPlinkMap(String filename, String SNP_DB) {
		PrintWriter writer;
		String[] positions;
		Hashtable<String, String> hash = null;

		if (SNP_DB == null) {

		} else if (!new File(SNP_DB).exists()) {
			System.err.println("Error - could not find SNP database: " + SNP_DB);
		} else {
			hash = HashVec.loadFileToHashString(SNP_DB, 0, new int[] {1, 2}, "\t", false);
		}

		try {
			writer = Files.openAppropriateWriter(filename);
			for (int i = 0; i < markerNames.length; i++) {
				if (hash == null) {
					positions = new String[] {chr + "", (i + 1) + ""};
				} else if (!hash.containsKey(markerNames[i])) {
					System.err.println("Error - marker '" + markerNames[i] + "' not found in database: "
														 + SNP_DB);
					positions = new String[] {"-1", "-1"};
				} else {
					positions = hash.get(markerNames[i]).split(PSF.Regex.GREEDY_WHITESPACE);
					if (Integer.parseInt(positions[0]) != chr) {
						System.err.println("Error - " + SNP_DB + " places " + markerNames[i]
															 + " on a different chromosome (" + positions[0] + " instead of "
															 + chr + ")");
					}
				}
				writer.println(positions[0] + "\t" + markerNames[i] + "\t0\t" + positions[1]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing " + filename);
			e.printStackTrace();
		}
	}

	public void filter(boolean[] keeps) {
		int count;
		String[] oldMarkerNames;
		double[][] oldAlleleFreqs;
		double[] old_Cumulative_Positions;

		if (markerNames.length != keeps.length) {
			System.err.println("I do not think you are keeping the markers you think you are keeping...");
			System.exit(1);
		}

		oldMarkerNames = markerNames;
		oldAlleleFreqs = alleleFreqs;
		old_Cumulative_Positions = positions;

		count = ArrayUtils.booleanArraySum(keeps);
		markerNames = new String[count];
		alleleFreqs = new double[count][];
		positions = new double[count];
		count = 0;
		for (int i = 0; i < keeps.length; i++) {
			if (keeps[i]) {
				markerNames[count] = oldMarkerNames[i];
				alleleFreqs[count] = oldAlleleFreqs[i];
				positions[count] = old_Cumulative_Positions[i];
				count++;
			}
		}
		distances = computeDistancesFromPositions(positions);

		// Files.backup(filename, "", "");
		// createFile();
	}

	public void alterPenetrance(String fileout, double disease, double pp, double pq, double qq,
															boolean suppressCheck) {
		if (disease < 0.50 && qq < pp && !suppressCheck) {
			System.err.println("Error - you appear to have your penetrance function backward");
			System.exit(1);
		}

		dxAlleleFreq = disease;
		if (chr == 23) {
			penetranceFunctions = new double[][] {{pp, pq, qq}, {pp, qq}};
		} else {
			penetranceFunctions = new double[][] {{pp, pq, qq}};
		}
		createFile(fileout);
	}

	public void createDominantMap() {
		if (chr == 23) {
			penetranceFunctions = DEFAULT_X_LINKED_DOMINANT_PENTRANCE;
		} else {
			penetranceFunctions = DEFAULT_AUTOSOMAL_DOMINANT_PENTRANCE;
		}
		if (!filename.endsWith(".dat")) {
			System.err.println("Error - unexpected filename extension: " + filename);
		}
		createFile(filename.substring(0, filename.length() - 4) + ".D"
							 + filename.substring(filename.length() - 4));
	}

	public void createRecessiveMap() {
		if (chr == 23) {
			penetranceFunctions = DEFAULT_X_LINKED_RECESSIVE_PENTRANCE;
		} else {
			penetranceFunctions = DEFAULT_AUTOSOMAL_RECESSIVE_PENTRANCE;
		}
		if (!filename.endsWith(".dat")) {
			System.err.println("Error - unexpected filename extension: " + filename);
		}
		createFile(filename.substring(0, filename.length() - 4) + ".R"
							 + filename.substring(filename.length() - 4));
	}

	public void updatePositions(String markerPositionDatabase, int markerIndex, int chrIndex,
															int locIndex) {
		BufferedReader reader;
		String[] line;
		String trav;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		double[] newPositions;
		boolean problem = false;
		boolean wasInMorgans = currentlyInMorgans;

		for (String markerName : markerNames) {
			hash.put(markerName, "");
		}

		try {
			reader = new BufferedReader(new FileReader(markerPositionDatabase));
			while (reader.ready()) {
				line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
				if (hash.containsKey(line[markerIndex])) {
					if (Integer.parseInt(line[chrIndex]) != chr) {
						System.err.println("Error - " + markerPositionDatabase + " says that "
															 + line[markerIndex] + " is on " + Integer.parseInt(line[chrIndex])
															 + " and not " + chr);
					} else {
						hash.put(line[markerIndex], line[locIndex]);
					}
				}

			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + markerPositionDatabase
												 + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + markerPositionDatabase + "\"");
			System.exit(2);
		}

		newPositions = new double[markerNames.length];
		for (int i = 0; i < markerNames.length; i++) {
			trav = hash.get(markerNames[i]);
			if (trav.equals("")) {
				System.err.println("Error - marker " + markerNames[i] + " was not found in "
													 + markerPositionDatabase);
				problem = true;
			} else {
				newPositions[i] = Double.parseDouble(trav);
				if (i > 0 && newPositions[i] <= newPositions[i - 1]) {
					if (newPositions[i] - newPositions[i - 1] > MIN_DIST_STEP * 10) {
						System.err.println("Error - make sure that " + markerNames[i - 1] + " and "
															 + markerNames[i] + " have not flipped order");
					}
					newPositions[i] = newPositions[i - 1] + MIN_DIST_STEP;
				}
			}
		}
		if (problem) {
			System.err.println("  Distances were NOT updated");
		} else {
			positions = newPositions;
			distances = computeDistancesFromPositions(positions);
			currentlyInMorgans = false;
			Files.backup(filename, "", "");
			createFile(filename, wasInMorgans);
		}
	}

	public String[][] recode(CountVector[] cvs) {
		String[] alleleSizes;
		int[] alleleCounts, keys;
		int sum;
		String[][] orderedAlleles;
		int numMissing;

		if (markerNames.length != cvs.length) {
			System.err.println("Error - marker number mismatch recoding chromosome " + chr);
			System.exit(1);
		}

		orderedAlleles = new String[markerNames.length][];
		addInfo = new String[markerNames.length];
		alleleFreqs = new double[markerNames.length][];
		for (int i = 0; i < markerNames.length; i++) {
			alleleSizes = cvs[i].getValues();
			alleleCounts = cvs[i].getCounts();
			keys = Sort.getSortedIndices(alleleSizes);
			sum = 0;
			numMissing = alleleSizes[keys[0]].equals("0") ? 1 : 0;
			orderedAlleles[i] = new String[alleleSizes.length - numMissing];
			for (int j = numMissing; j < alleleSizes.length; j++) {
				orderedAlleles[i][j - numMissing] = alleleSizes[keys[j]];
				sum += alleleCounts[keys[j]];
			}
			addInfo[i] = "  recode : " + ArrayUtils.toStr(orderedAlleles[i], " ");
			alleleFreqs[i] = new double[alleleSizes.length - numMissing];
			for (int j = numMissing; j < alleleSizes.length; j++) {
				alleleFreqs[i][j - numMissing] = (double) alleleCounts[keys[j]] / sum;
			}
		}

		return orderedAlleles;
	}

	public void writeCumulative() {
		PrintWriter writer;

		try {
			writer = Files.openAppropriateWriter(ext.rootOf(filename) + "_map.xln");
			writer.println("MarkerName\tPosition");
			for (int i = 0; i < markerNames.length; i++) {
				writer.println(markerNames[i] + "\t" + positions[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing sumulative summary");
			e.printStackTrace();
		}
	}

	public static double[][] createDummyAlleleFreqs(int numMarkers, int numDefaultAlleles) {
		double[][] alleleFreqs = new double[numMarkers][];

		for (int i = 0; i < numMarkers; i++) {
			alleleFreqs[i] = new double[numDefaultAlleles];
			for (int j = 0; j < numDefaultAlleles; j++) {
				alleleFreqs[i][j] = 1.0 / numDefaultAlleles;
			}
		}

		return alleleFreqs;
	}

	public static double[][] parseModels(String filename) {
		String[] line;
		Vector<String> v;
		double[][] models;

		v = HashVec.loadFileToVec(filename, false, false, false);
		models = new double[v.size()][];
		for (int i = 0; i < v.size(); i++) {
			line = v.elementAt(i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
			if (line.length != 4) {
				System.err.println("Error - requires 4 parameters: freq(q)\tP(pp)\tP(pq)\tP(qq)");
				System.err.println("        found: " + v.elementAt(i));
				System.exit(1);
			}
			models[i] = ArrayUtils.toDoubleArray(line);
		}

		return models;
	}

	public static void mapStats() {
		double sumDist;
		int totalMarkers;
		LinkageMap map;
		double[] dists;
		IntVector iv = new IntVector();
		boolean currentlyInMorgans = true;

		sumDist = 0;
		totalMarkers = 0;
		for (int i = 1; i <= 23; i++) {
			if (new File("map" + ext.chrome(i) + ".dat").exists()) {
				iv.add(i);
				map = new LinkageMap(i);

				if (iv.size() == 1) {
					currentlyInMorgans = map.currentlyInMorgans;
				} else if (currentlyInMorgans != map.currentlyInMorgans) {
					System.err.println("Error - some chromosomes are in Morgans, some are in centiMorgans");
				}

				totalMarkers += map.getMarkerNames().length;
				dists = map.getDistances(currentlyInMorgans);
				for (int j = 1; j < dists.length; j++) {
					sumDist += dists[j];
				}
			}
		}
		System.out.println("There are " + totalMarkers + " markers at an average distance of "
											 + ext.formDeci(sumDist / (totalMarkers - iv.size()), 2)
											 + (currentlyInMorgans ? "cM" : "Morgans"));
	}

	public static void updateMaps(String dir, String db, int markerIndex, int chrIndex,
																int locIndex) {
		LinkageMap map;
		IntVector iv = new IntVector();

		for (int i = 1; i <= 23; i++) {
			if (new File(dir + "map" + ext.chrome(i) + ".dat").exists()) {
				map = new LinkageMap(dir, i);
				map.updatePositions(dir + db, markerIndex, chrIndex, locIndex);
			} else {
				iv.add(i);
			}
		}

		if (iv.size() > 0) {
			System.err.println("Skipped the following chromosomes as they could not be found: "
												 + ArrayUtils.toStr(iv.toArray()));

		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		double dx = 0.005, pp = 0.80, pq = 0.80, qq = 0.03;
		String filename = "map02.dat";
		boolean writeDist = false;

		String usage = "\n" + "link.LinkageMap requires 2-5 arguments:\n"
									 + "   (1) filename (i.e. file=" + filename + " (default))\n"
									 + "   (2) disease allele frequency (a.k.a. q) (i.e. dx=" + dx + " (default))\n"
									 + "   (3) penetrance of qq (i.e. qq=" + qq + " (default))\n"
									 + "   (4) penetrance of pq (i.e. pq=" + pq + " (default))\n"
									 + "   (5) penetrance of pp (i.e. pp=" + pp + " (default))\n" + " OR\n"
									 + "   (2) write distances (i.e. -dist (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("dx=")) {
				dx = Double.valueOf(arg.split("=")[1]).doubleValue();
				numArgs--;
			} else if (arg.startsWith("pp=")) {
				pp = Double.valueOf(arg.split("=")[1]).doubleValue();
				numArgs--;
			} else if (arg.startsWith("pq=")) {
				pq = Double.valueOf(arg.split("=")[1]).doubleValue();
				numArgs--;
			} else if (arg.startsWith("qq=")) {
				qq = Double.valueOf(arg.split("=")[1]).doubleValue();
				numArgs--;
			} else if (arg.startsWith("-dist")) {
				writeDist = true;
				numArgs--;
			} else {
				System.err.println("Error - bad argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (writeDist) {
				new LinkageMap(filename).writeCumulative();
			} else {
				new LinkageMap(filename).alterPenetrance(filename, dx, pp, pq, qq, false);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
