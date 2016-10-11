// -Xmx1024M
// added in centiMoragns to the 4th column, make sure that's not causing any trouble for you
package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.genvisis.bioinformatics.MapSNPsAndGenes;
import org.genvisis.bioinformatics.ParseSNPlocations;
import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.common.Unique;
import org.genvisis.common.ext;
import org.genvisis.link.LinkageMap;

import com.google.common.primitives.Ints;

public class SnpMarkerSet implements Serializable, PlainTextExport {
	public static final long serialVersionUID = 1L;
	public static final int NAMES_ONLY = 0;
	public static final int GENERIC_FORMAT = 1;
	public static final int PLINK_BIM_FORMAT = 2;
	public static final int PLINK_MAP_FORMAT = 3;
	public static final int EIGENSOFT_SNP_FORMAT = 4;
	public static final int GENERIC_FORMAT_IGNORE_FIRST_LINE = 5;
	public static final int GENERIC_FORMAT_ANNOTATED = 6;
	public static final int GENERIC_FORMAT_ANNOTATED_IGNORE_FIRST_LINE = 7;
	public static final int MACH_MLINFO_FORMAT = 8;
	public static final int HAPLOVIEW_INFO_FORMAT = 9;
	public static final int MINIMAC_INFO_FORMAT = 10;
	public static final int PLINK_MAP_FORMAT_WITHOUT_CM = 11;
	public static final int INFO_FOR_BURDEN_TESTING = 12;
	public static final int IMPUTE2_INFO_FORMAT = 13;
	public static final int FREEZE5_FORMAT = 14;

	// public static final String[][] NATURAL_PAIRINGS = {{"A", "T"}, {"C", "G"}};
	public static final String[] NULL_ALLLES = {"0", "-"};
	public static final String[] INDEL_ALLLES = {"I", "D"};

	public static final String[] HAPMAP_CM_SRC_HEADER = {	"position", "COMBINED_rate(cM/Mb)",
																												"Genetic_Map(cM)"};
	public static final String[][] HEADERS = {null, null, null, null, null,
																						{"Marker", "Chr", "Position"}, null,
																						{"Marker", "Chr", "Position", "Annotation"},
																						{"SNP", "Al1", "Al2", "Freq1", "MAF", "Quality", "Rsq"},
																						null,
																						{	"SNP", "Al1", "Al2", "Freq1", "MAF", "AvgCall", "Rsq",
																							"Genotyped", "LooRsq", "EmpR", "EmpRsq", "Dose1",
																							"Dose2"},
																						null,
																						{	"Marker", "Chr", "Position", "REF", "ALT", "gene",
																							"AAF", "Function"},
																						{	"snp_id", "rs_id", "position", "exp_freq_a1", "info",
																							"certainty", "type", "info_type0", "concord_type0",
																							"r2_type0"},
																						null,};

	public static final String[][] HEADER_ELEMEMTS = {Aliases.MARKER_NAMES, Aliases.CHRS,
																										Aliases.POSITIONS, Aliases.CENTIMORGANS,
																										Aliases.ALLELES[0], Aliases.ALLELES[1]};

	public static final int CHR_INFO_IN_FILENAME = -2;

	/** 0 1 2 3 4 5 6+ */
	/** Marker name, Chr, Position, centiMorgans, A1, A2, annotation */
	public static final int[][] INDICES = {	{0, -1, -1, -1, -1, -1}, {0, 1, 2, -1, -1, -1},
																					{1, 0, 3, 2, 4, 5}, {1, 0, 3, 2, -1, -1},
																					{0, 1, 3, 2, 4, 5}, {0, 1, 2, -1, -1, -1},
																					{0, 1, 2, -1, -1, -1, 3}, {0, 1, 2, -1, -1, -1, 3},
																					{0, -1, -1, -1, 1, 2, 3, 4, 5, 6}, {0, -1, 1, -1, -1, -1},
																					{0, -1, -1, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
																					{1, 0, 2, -1, -1, -1}, {0, 1, 2, -1, 3, 4, 5, 6, 7},
																					{	1, CHR_INFO_IN_FILENAME, 2, -1, -1, -1, 3, 4, 5, 6, 7,
																						8,
																						9},
																					{0, 1, 2, -1, 3, 4},
			// make sure to add an entry into HEADERS as well
	};

	private long fingerprint;
	private int[] rsNumbers;
	private byte[] chrs;
	private int[] positions;
	private double[] centiMorgans;
	private char[][] alleles;
	private String[][] annotation;
	// if any more information is added, BE SURE TO ADD TO sortMarkers()
	private Vector<String> nonRSmarkerNames;
	private int build;

	public SnpMarkerSet(String[] markerNames) {
		convertMarkerNamesToRSnumbers(markerNames, true, new Logger());
	}

	public SnpMarkerSet(String[] markerNames, byte[] chrs, int[] positions) {
		this(markerNames, chrs, positions, null, null, false, true, new Logger());
	}

	public SnpMarkerSet(String[] markerNames, byte[] chrs, int[] rawPositions, char[][] alleles,
											String[][] annotation, boolean sort, boolean verbose) {
		this(markerNames, chrs, rawPositions, alleles, annotation, sort, verbose, new Logger());
	}

	public SnpMarkerSet(String[] markerNames, byte[] chrs, int[] rawPositions, char[][] alleles,
											String[][] annotation, boolean sort, boolean verbose, Logger log) {
		if (markerNames.length != chrs.length	|| markerNames.length != rawPositions.length
				|| (annotation != null && annotation.length != markerNames.length)) {
			log.reportError("Error - mismatched number of markers and positions/annotations");
			System.exit(1);
		}

		nonRSmarkerNames = null;
		this.chrs = chrs;
		positions = rawPositions;
		this.alleles = alleles;
		this.annotation = annotation;
		build = ParseSNPlocations.DEFAULT_BUILD;

		convertMarkerNamesToRSnumbers(markerNames, verbose, log);
		fingerprint = fingerprint(rsNumbers);

		if (sort) { // used to be before conversion, but failed due to null rsNumber list
			sortMarkers();
		}
	}

	public SnpMarkerSet(String filename) {
		this(filename, true, new Logger());
	}

	public SnpMarkerSet(String filename, boolean verbose, Logger log) {
		this(filename, determineType(filename), verbose, log);
	}

	public SnpMarkerSet(String filename, int type, boolean verbose, Logger log) {
		this(	filename, type == -1 ? determineIndices(filename, log) : INDICES[type],
					type == -1 ? true : HEADERS[type] != null, verbose, log);
	}

	public static int[] determineIndices(String filename, Logger log) {
		log.report("Assuming there is a header");
		return ext.indexFactors(HEADER_ELEMEMTS, Files.getHeaderOfFile(filename, null, log), true, true,
														true, log, false);
	}

	public SnpMarkerSet(String filename, int[] indices, boolean header, boolean verbose, Logger log) {
		this(filename, Files.determineDelimiter(filename, log), indices, header, null, verbose, log);
	}

	public SnpMarkerSet(String filename, String delimiter, int[] indices, boolean header,
											String[][] stringReplacements, boolean verbose, Logger log) {
		String[] markerNames;


		if (indices == null) {
			String temp = Files.getFirstNLinesOfFile(filename, 1, log)[0];
			temp = ext.replaceAllWith(temp, stringReplacements);
			indices = ext.indexFactors(	SnpMarkerSet.HEADER_ELEMEMTS, temp.split(delimiter), true, true,
																	false, log, false);
			header = true;
		}

		String[] data = HashVec.loadFileToStringArray(filename, header, null, false);

		if (data == null) {
			return;
		}

		int count = data.length;

		markerNames = new String[count];
		if (indices[1] == CHR_INFO_IN_FILENAME) {
			Matcher m = Pattern.compile(DosageData.CHR_REGEX).matcher(filename);
			byte chr = -1;
			if (m.matches()) {
				chr = (byte) Integer.parseInt(m.group(1));
				log.report("Warning - the format given expects chromosome number to be part of the file name.  This was determined to be chr{"
										+ chr + "}.");
				chrs = Array.byteArray(count, chr);
			} else {
				log.reportError("Error - the format given expects chromosome number to be part of the file name, but no chromosome number was found.  Chromosome information will not be included.");
				chrs = new byte[count];
			}
		} else {
			chrs = new byte[count];
		}
		positions = new int[count];

		if (indices[3] == -1) {
			centiMorgans = null;
		} else {
			centiMorgans = new double[count];
		}
		if (indices[4] == -1 || indices[5] == -1) {
			alleles = null;
		} else {
			alleles = new char[count][2];
		}
		if (indices.length == 6) {
			annotation = null;
		} else {
			annotation = new String[count][indices.length - 6];
		}
		for (int i = 0; i < count; i++) {
			String temp = data[i];
			if (stringReplacements != null) {
				temp = ext.replaceAllWith(temp, stringReplacements);
			}
			String[] line = temp.split(delimiter, -1);
			markerNames[i] = line[indices[0]];
			if (indices[1] >= 0) {
				chrs[i] = Positions.chromosomeNumber(line[indices[1]]);
			}
			if (indices[2] != -1) {
				positions[i] = Integer.parseInt(line[indices[2]]);
			}
			if (indices[3] != -1) {
				centiMorgans[i] = Double.parseDouble(line[indices[3]]);
			}
			if (indices[4] != -1) {
				alleles[i][0] = line[indices[4]].charAt(0);
			}
			if (indices[5] != -1) {
				alleles[i][1] = line[indices[5]].charAt(0);
			}
			if (indices.length > 6) {
				for (int j = 6; j < indices.length; j++) {
					try {
						annotation[i][j - 6] = line[indices[j]];
					} catch (Exception e) {
						log.reportError("Unexpected exception in SnpMarker constructor");
					}
				}
			}
		}
		convertMarkerNamesToRSnumbers(markerNames, verbose, log);
	}

	public void convertMarkerNamesToRSnumbers(String[] markerNames, boolean verbose, Logger log) {
		int countNon;

		rsNumbers = new int[markerNames.length];
		nonRSmarkerNames = new Vector<String>();

		countNon = 0;
		for (int i = 0; i < markerNames.length; i++) {
			if (markerNames[i].toLowerCase().startsWith("rs")	&& markerNames[i].indexOf("_") == -1
					&& markerNames[i].indexOf(":") == -1) {
				try {
					rsNumbers[i] = Integer.parseInt(markerNames[i].substring(2));
				} catch (NumberFormatException nfe) {
					rsNumbers[i] = -1 * nonRSmarkerNames.size();
					nonRSmarkerNames.add(markerNames[i]);
					if (verbose && countNon < 20) {
						log.reportError("FYI, SNP with an addedum after the rs number: '"	+ markerNames[i]
														+ "'");
					}
					countNon++;
				}
			} else {
				rsNumbers[i] = -1 * nonRSmarkerNames.size();
				nonRSmarkerNames.add(markerNames[i]);
				if (verbose && countNon < 20) {
					log.reportError("FYI, SNP without an rs number: '" + markerNames[i] + "'");
				}
				countNon++;
			}
		}
		if (verbose && countNon >= 20) {
			log.reportError("FYI, total of " + countNon + " markers without an rs number");
		}
	}

	public int[] getRSnumbers() {
		return rsNumbers;
	}

	/**
	 * Sets the genome build for this marker set
	 *
	 * @param build build of the genome (36, 37, 38, etc)
	 */
	public void setBuild(int build) {
		this.build = build;
	}

	/**
	 * Returns the genome build for this marker set
	 *
	 * @return build of the genome (36, 37, 38, etc)
	 */
	public int getBuild() {
		return build;
	}

	public String[] getMarkerNames() {
		String[] markerNames = new String[rsNumbers.length];

		for (int i = 0; i < rsNumbers.length; i++) {
			if (rsNumbers[i] > 0) {
				markerNames[i] = "rs" + rsNumbers[i];
			} else {
				markerNames[i] = nonRSmarkerNames.elementAt(-1 * rsNumbers[i]);
			}
		}
		return markerNames;
	}

	public byte[] getChrs() {
		return chrs;
	}

	public int[] getPositions() {
		return positions;
	}

	/**
	 * Caution - if markers were not in order, positions will not be in order (meaning binarySearch
	 * will be inaccurate)
	 *
	 * @return
	 */
	public int[][] getPositionsByChr() {
		int[][] posByChr;

		HashMap<Integer, ArrayList<Integer>> posMap = new HashMap<Integer, ArrayList<Integer>>();
		for (int i = 0; i < 27; i++) {
			posMap.put(i, new ArrayList<Integer>());
		}
		for (int i = 0; i < chrs.length; i++) {
			posMap.get((int) chrs[i]).add(positions[i]);
		}

		posByChr = new int[27][0];
		for (int i = 0; i < 27; i++) {
			posByChr[i] = Ints.toArray(posMap.get(i));
		}

		return posByChr;
	}

	// public int[][] getPositionsByChr() {
	// IntVector iv;
	// byte chr;
	// int[][] positionsByChr;
	// boolean done;
	//
	// positionsByChr = new int[27][0];
	//
	// chr = 0;
	// iv = new IntVector(20000);
	// done = false;
	// for (int i = 0; !done; i++) {
	// if (i==chrs.length || chrs[i] != chr) {
	// positionsByChr[chr] = iv.toArray();
	// chr = i==chrs.length?0:chrs[i];
	// iv = new IntVector(20000);
	// }
	// if (i==chrs.length) {
	// done = true;
	// } else {
	// iv.add(positions[i]);
	// }
	// }
	//
	// return positionsByChr;
	// }

	public int[][] getChrAndPositionsAsInts() {
		int[][] chrPositions = new int[positions.length][];

		for (int i = 0; i < positions.length; i++) {
			chrPositions[i] = new int[] {chrs[i], positions[i]};
		}

		return chrPositions;
	}

	public String[] getChrAndPositions() {
		String[] chrPositions = new String[positions.length];

		for (int i = 0; i < positions.length; i++) {
			chrPositions[i] = chrs[i] + "\t" + positions[i];
		}

		return chrPositions;
	}

	public double[] getCentiMorgans() {
		return centiMorgans;
	}

	public char[][] getAlleles() {
		return alleles;
	}

	public String[][] getAnnotation() {
		return annotation;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	@Override
	public void exportToText(String outputFile, Logger log) {
		writeToFile(outputFile, determineType(outputFile), log);
	}

	public void writeToFile(String filename, int format, Logger log) {
		writeToFile(filename, INDICES[format], HEADERS[format], log);
	}

	public void writeToFile(String filename, int[] indices, String[] header, Logger log) {
		PrintWriter writer;
		String[] line, markerNames;

		markerNames = getMarkerNames();
		line = new String[Array.max(indices) + 1];
		try {
			writer = new PrintWriter(new FileWriter(filename));
			if (header != null) {
				writer.println(Array.toStr(header));
			}
			for (int i = 0; i < markerNames.length; i++) {
				for (int j = 0; j < line.length; j++) {
					line[j] = "0";
				}
				if (indices[0] != -1) {
					line[indices[0]] = markerNames[i];
				}
				if (indices[1] != -1) {
					line[indices[1]] = chrs == null ? "." : chrs[i] + "";
				}
				if (indices[2] != -1) {
					line[indices[2]] = positions == null ? "." : positions[i] + "";
				}
				if (indices[3] != -1) {
					line[indices[3]] = centiMorgans == null ? "0" : ext.formDeci(centiMorgans[i], 5, false);
				}
				if (indices[4] != -1) {
					line[indices[4]] = alleles == null ? "." : alleles[i][0] + "";
				}
				if (indices[5] != -1) {
					line[indices[5]] = alleles == null ? "." : alleles[i][1] + "";
				}
				if (indices.length > 6) {
					for (int j = 6; j < indices.length; j++) {
						line[indices[j]] = annotation == null ? "." : annotation[i][j - 6];
					}
				}
				writer.println(Array.toStr(line));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing " + filename);
			e.printStackTrace();
		}
	}

	public void autofillCentiMorgans() {
		centiMorgans = new double[positions.length];
		for (int i = 0; i < positions.length; i++) {
			centiMorgans[i] = (double) positions[i] / 100000000;
		}
	}

	public void interpolateCentiMorgans(SnpMarkerSet sourceSet, Logger log) {
		byte[] srcChrs;
		int[] srcPositions;
		double[] srcCentiMorgans;
		int[] starts, stops;
		int pos, index;
		byte chr;

		srcChrs = sourceSet.getChrs();
		srcPositions = sourceSet.getPositions();
		srcCentiMorgans = sourceSet.getCentiMorgans();

		starts = new int[27];
		stops = new int[27];

		chr = 0;
		pos = -1;
		for (int i = 0; i < srcPositions.length; i++) {
			if (srcChrs[i] != chr) {
				stops[chr] = i;
				chr = srcChrs[i];
				starts[chr] = i;
				pos = -1;
			}
			if (srcPositions[i] < pos) {
				log.reportError("Error - centiMorgans source MarkerSet was not in order starting at marker "
													+ sourceSet.getMarkerNames()[i] + " (" + srcPositions[i] + " < " + pos
												+ ")");
				return;
			}
			pos = srcPositions[i];
		}

		centiMorgans = new double[positions.length];
		int[] counts = new int[10];
		for (int i = 0; i < positions.length; i++) {
			index =
						Array.binarySearch(srcPositions, positions[i], starts[chrs[i]], stops[chrs[i]], false);
			if (chrs[i] == srcChrs[index] && positions[i] == srcPositions[index]) {
				counts[0]++;
				centiMorgans[i] = srcCentiMorgans[index];
			} else if (index == starts[chrs[i]]) {
				counts[1]++;
				// log.reportError(0+"\t"+srcPositions[index]);
				centiMorgans[i] = -1;
				// centiMorgans[i] =
				// Math.max(srcCentiMorgans[index]-(srcPositions[index]-positions[i])/1000000, 0);
			} else if (index == stops[chrs[i]] + 1) {
				counts[2]++;
				// centiMorgans[i] = srcCentiMorgans[index]+(positions[i]-srcPositions[index])/1000000;
				centiMorgans[i] = -2;
			} else if (chrs[i] == srcChrs[index]	&& positions[i] > srcPositions[index - 1]
									&& positions[i] < srcPositions[index]) {
				counts[3]++;
				// log.reportError(srcPositions[index-1]+"\t"+srcPositions[index]);
				centiMorgans[i] = srcCentiMorgans[index - 1]
													+ (srcCentiMorgans[index] - srcCentiMorgans[index - 1])
															* (positions[i] - srcPositions[index - 1])
														/ (srcPositions[index] - srcPositions[index - 1]);
			} else {
				log.reportError("Error - finding nearby markers");
			}
		}

		log.report("Out of " + chrs.length + " markers:");
		log.report(counts[0] + " matched exactly");
		log.report(counts[1] + " came before the first marker on the chromosome of the reference map");
		log.report(counts[2] + " came after last marker on the chromosome of the reference map");
		log.report(counts[3] + " were interpolated");
		log.report("All " + Array.sum(counts) + " should be accounted for");


	}

	public LinkageMap createLinkageMap() {
		String[] markerNames = getMarkerNames();
		double[] cM_data = new double[positions.length];
		int chr = chrs[0];

		for (int i = 0; i < markerNames.length; i++) {
			cM_data[i] = (double) positions[i] / 100000000; // actually Morgans
			if (chrs[i] != chr) {
				System.err.println("Error - rs"	+ rsNumbers[i] + " is on a different chromosome (" + chrs[i]
														+ ") than the others (" + chr + ")");
			}
		}

		if (Array.max(centiMorgans) > 0) {
			cM_data = centiMorgans;
		}

		return new LinkageMap(chr, markerNames, 2, cM_data, true, true);
	}

	public void listUnambiguousMarkers(	String filename, String excludeMarkersFile,
																			boolean autosomesOnly) {
		PrintWriter writer;
		int countAs_and_Ts, countGs_and_Cs;
		String[] markerNames;

		if (alleles == null) {
			System.err.println("Error - can't define unambiguous alleles if no alleles are listed");
			return;
		}

		Set<String> excludeMarkers = excludeMarkersFile == null	? new HashSet<String>()
																														: HashVec.loadFileToHashSet(excludeMarkersFile,
																																												false);

		try {
			writer = Files.getAppropriateWriter(filename);
			markerNames = getMarkerNames();
			for (int i = 0; i < markerNames.length; i++) {
				if (!excludeMarkers.contains(markerNames[i])) {
					countAs_and_Ts = countGs_and_Cs = 0;
					for (int j = 0; j < alleles[i].length; j++) {
						if (alleles[i][j] == 'A' || alleles[i][j] == 'T') {
							countAs_and_Ts++;
						} else if (alleles[i][j] == 'G' || alleles[i][j] == 'C') {
							countGs_and_Cs++;
						} else if (ext.indexOfStr(alleles[i][j] + "", NULL_ALLLES) == -1
												&& ext.indexOfStr(alleles[i][j] + "", INDEL_ALLLES) == -1) {
							System.err.println("Error - invalid allele '"	+ alleles[i][j] + "' for marker "
																	+ markerNames[i]);
						}
					}
					if (countAs_and_Ts == 1	&& countGs_and_Cs == 1 && chrs[i] > 0
							&& (!autosomesOnly || chrs[i] < 23)) {
						writer.println(markerNames[i]);
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static SnpMarkerSet load(String filename, boolean jar) {
		return load(filename, jar, new Logger());
	}

	public static SnpMarkerSet load(String filename, boolean jar, Logger log) {
		return (SnpMarkerSet) SerializedFiles.readSerial(filename, jar, log, true);
	}

	public static long fingerprint(int[] rsNumbers) {
		long sum;

		sum = 0;
		for (int i = 0; i < rsNumbers.length; i++) {
			sum += rsNumbers[i] * (i + 1);
		}

		return sum;
	}

	public String[] findNearby(SnpMarkerSet nearbyThese, int withinXbp) {
		int[] localOrder, anchorOrder;
		byte[] anchorChrs;
		int[] anchorPositions;
		int localIndex, anchorIndex;
		boolean[] overlaps;
		String[] markerNames, finalOverlaps;
		int count;

		anchorChrs = nearbyThese.getChrs();
		anchorPositions = nearbyThese.getPositions();

		anchorOrder = Sort.orderTwoLayers(anchorChrs, anchorPositions, new Logger());
		localOrder = Sort.orderTwoLayers(chrs, positions, new Logger());

		anchorIndex = localIndex = 0;
		overlaps = Array.booleanArray(positions.length, false);
		while (anchorIndex < anchorOrder.length && localIndex < localOrder.length) {
			if (chrs[localOrder[localIndex]] < anchorChrs[anchorOrder[anchorIndex]]
					|| (chrs[localOrder[localIndex]] == anchorChrs[anchorOrder[anchorIndex]]
							&& positions[localOrder[localIndex]] < anchorPositions[anchorOrder[anchorIndex]]
																											- withinXbp)) {
				localIndex++;
			} else if (chrs[localOrder[localIndex]] == anchorChrs[anchorOrder[anchorIndex]]
										&& positions[localOrder[localIndex]] >= anchorPositions[anchorOrder[anchorIndex]]
																													- withinXbp
									&& positions[localOrder[localIndex]] <= anchorPositions[anchorOrder[anchorIndex]]
																													+ withinXbp) {
				overlaps[localIndex] = true;
				localIndex++;
			} else {
				anchorIndex++;
			}
		}

		finalOverlaps = new String[Array.booleanArraySum(overlaps)];
		markerNames = getMarkerNames();
		count = 0;
		for (int i = 0; i < overlaps.length; i++) {
			if (overlaps[i]) {
				finalOverlaps[count] = markerNames[i];
				count++;
			}
		}

		return finalOverlaps;
	}

	public void parseSNPlocations() {
		parseSNPlocations(new Logger());
	}

	public void parseSNPlocations(Logger log) {
		parseSNPlocations(MapSNPsAndGenes.getSNPDB(build, log), MapSNPsAndGenes.getMergeDB(log), log);
	}

	public void setPositions(SnpMarkerSet refSet) {
		Hashtable<String, String> hash;
		String[] markerNames, line;

		hash = refSet.getChrHash();
		markerNames = getMarkerNames();
		for (int i = 0; i < markerNames.length; i++) {
			if (chrs[i] == 0 && positions[i] == 0 && hash.containsKey(markerNames[i])) {
				line = hash.get(markerNames[i]).split("[\\s]+");
				chrs[i] = Byte.parseByte(line[0]);
				positions[i] = Integer.parseInt(line[1]);
			}
		}
	}

	public void setAlleles(char[][] alleles2) {
		alleles = alleles2;
	}

	public void parseSNPlocations(String db, String mergeDB, Logger log) {
		SnpMarkerSet dbMarkerSet;
		int[] dbRSnumbers, dbPositions;
		String[] markerNames;
		byte[] dbChrs;
		int index;
		Hashtable<Integer, Integer> mergeHash;
		Integer trav, next;
		boolean annotateMerges;

		log.report(ext.getTime());
		markerNames = getMarkerNames();

		log.report("Loading database...");

		dbMarkerSet = SnpMarkerSet.load(db, false);
		dbRSnumbers = dbMarkerSet.getRSnumbers();
		dbPositions = dbMarkerSet.getPositions();
		dbChrs = dbMarkerSet.getChrs();

		log.report("Searching database...");

		annotateMerges = false;
		if (annotation == null) {
			annotateMerges = true;
			annotation = new String[markerNames.length][1];
		}

		mergeHash = null;
		chrs = new byte[markerNames.length];
		positions = new int[markerNames.length];
		for (int i = 0; i < markerNames.length; i++) {
			if (markerNames[i].startsWith("rs")) {
				index =
							Array.binarySearch(dbRSnumbers, Integer.parseInt(markerNames[i].substring(2)), true);
				if (index == -1) {
					if (mergeDB != null) {
						if (mergeHash == null) {
							log.report("Could not find a marker, loading merge database...", false, true);
							if (new File(mergeDB).exists()) {
								mergeHash = SerialHash.loadSerializedIntHash(mergeDB);
								log.report("done");
							} else {
								log.report("failed; " + mergeDB + " not found");
								mergeHash = new Hashtable<Integer, Integer>();
							}
						}
						trav = Integer.parseInt(markerNames[i].substring(2));
						next = trav;
						while (next != null) {
							trav = next;
							next = mergeHash.get(trav);
						}
						if (markerNames[i].equals("rs" + trav)) {
							log.reportError("\n\n****ERROR**** failed to find "	+ markerNames[i]
															+ " in any NCBI database ****ERROR****\n\n");
							if (annotateMerges) {
								annotation[i][0] = ".";
							}
						} else if (trav != null) {
							log.reportError("FYI - " + markerNames[i] + " has merged with rs" + trav);
							index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
							chrs[i] = dbChrs[index];
							positions[i] = dbPositions[index] + ParseSNPlocations.OFFSET;
							index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
							if (annotateMerges) {
								annotation[i][0] = "rs" + trav;
							}
						} else {
							log.reportError("Error - tried and failed to find a merge record for "
															+ markerNames[i]);
						}
					}
					if (index == -1) {
						log.reportError("Error - could not find " + markerNames[i] + " in " + db);
						chrs[i] = (byte) 0;
						positions[i] = 0;
						if (annotateMerges) {
							annotation[i][0] = ".";
						}
					}
				} else {
					if (dbChrs[index] == ParseSNPlocations.UNMAPPED) {
						log.reportError("Error - marker " + markerNames[i] + " is not mapped to a chromosome");
						chrs[i] = 0;
						positions[i] = -1;
					} else if (dbPositions[index] == ParseSNPlocations.MULTIPLE_POSITIONS) {
						log.reportError("Warning - marker "	+ markerNames[i]
														+ " likely has multiple positions on chromosome " + dbChrs[index]);
						chrs[i] = dbChrs[index];
						positions[i] = -1;
					} else {
						chrs[i] = dbChrs[index];
						positions[i] = dbPositions[index] + ParseSNPlocations.OFFSET;
					}
					if (annotateMerges) {
						annotation[i][0] = markerNames[i];
					}
				}
			} else {
				log.reportError("Error - can't look up a SNP without an rs number ("	+ markerNames[i]
												+ ")");
				chrs[i] = (byte) 0;
				positions[i] = 0;
				if (annotateMerges) {
					annotation[i][0] = ".";
				}
			}
		}

		log.report("Done with the positions! " + ext.getTime());
	}

	public void sortMarkers() {
		int[] keys;

		keys = Sort.orderTwoLayers(chrs, positions, new Logger());
		rsNumbers = Sort.putInOrder(rsNumbers, keys);
		chrs = Sort.putInOrder(chrs, keys);
		positions = Sort.putInOrder(positions, keys);
		if (alleles != null) {
			alleles = Sort.putInOrder(alleles, keys);
		}
		if (centiMorgans != null) {
			centiMorgans = Sort.putInOrder(centiMorgans, keys);
		}
		if (annotation != null) {
			annotation = Sort.putInOrder(annotation, keys);
		}
	}

	public Hashtable<String, String> getChrHash() {
		Hashtable<String, String> chrHash;
		String[] markerNames;
		byte[] chrs;
		int[] positions;

		markerNames = getMarkerNames();
		chrs = getChrs();
		positions = getPositions();
		chrHash = new Hashtable<String, String>(markerNames.length);
		for (int i = 0; i < markerNames.length; i++) {
			chrHash.put(markerNames[i], chrs[i] + "\t" + positions[i]);
		}

		return chrHash;
	}

	public Hashtable<String, Segment> getSegments() {
		Hashtable<String, Segment> segHash;
		String[] markerNames;
		byte[] chrs;
		int[] positions;

		markerNames = getMarkerNames();
		chrs = getChrs();
		positions = getPositions();
		segHash = new Hashtable<String, Segment>(markerNames.length);
		for (int i = 0; i < markerNames.length; i++) {
			segHash.put(markerNames[i], new Segment(chrs[i], positions[i], positions[i] + 1));
		}

		return segHash;
	}

	public static Hashtable<String, String> loadSnpMarkerSetToChrHash(String filenameSource) {
		return loadSnpMarkerSetToChrHash(filenameSource, determineType(filenameSource));
	}

	public static Hashtable<String, String> loadSnpMarkerSetToChrHash(String filenameSource,
																																		int snpMarkerSetType) {
		Hashtable<String, String> chrHash;
		long time;
		SnpMarkerSet markerSet;

		if (!Files.exists(filenameSource + ".hash.ser", false)) {
			System.out.print("Loading " + filenameSource);
			time = new Date().getTime();
			markerSet = new SnpMarkerSet(filenameSource, snpMarkerSetType, true, new Logger());
			System.out.println("...finished in " + ext.getTimeElapsed(time));
			time = new Date().getTime();
			System.out.print("Mapping master list to chromosomes");
			chrHash = markerSet.getChrHash();
			System.out.println("...finished in " + ext.getTimeElapsed(time));

			System.out.print("Saving " + filenameSource + ".hash.ser");
			time = new Date().getTime();
			SerialHash.createSerializedStringHash(filenameSource + ".hash.ser", chrHash);
			System.out.println("...finished in " + ext.getTimeElapsed(time));
		} else {
			time = new Date().getTime();
			System.out.print("Loading " + filenameSource + ".hash.ser");
			chrHash = SerialHash.loadSerializedStringHash(filenameSource + ".hash.ser");
			System.out.println("...loaded " + chrHash.size() + " markers in " + ext.getTimeElapsed(time));
		}

		return chrHash;
	}

	public SnpMarkerSet trim(int[][] rangesToKeep, boolean verbose, Logger log) {
		String[] markerNames, newMarkerNames;
		byte[] chrs, newChrs;
		int[] positions, newPositions;
		char[][] newAlleles;
		String[][] newAnnotation;

		markerNames = getMarkerNames();
		chrs = getChrs();
		positions = getPositions();
		boolean[] keep = Array.booleanArray(markerNames.length, false);
		for (int i = 0; i < markerNames.length; i++) {
			for (int[] range : rangesToKeep) {
				if (chrs[i] != range[0]) {
					continue; // wrong chr
				}
				if (positions[i] < range[1]) {
					continue; // smaller than range start
				}
				if (positions[i] > range[2]) {
					continue; // greater than range end
				}
				keep[i] = true;
				break;
			}
		}
		System.out.println("Trimming SnpMarkerSet: "	+ chrs.length + " chrs @ start, keeping "
												+ Array.booleanArraySum(keep));

		newMarkerNames = Array.subArray(markerNames, keep);
		newChrs = Array.subArray(chrs, keep);
		newPositions = Array.subArray(positions, keep);
		newAlleles = getAlleles() == null ? null : Array.subArray(getAlleles(), keep);
		newAnnotation = getAnnotation() == null ? null : Array.subArray(getAnnotation(), keep);

		return new SnpMarkerSet(newMarkerNames, newChrs, newPositions, newAlleles, newAnnotation, false,
														verbose, new Logger());
	}

	public SnpMarkerSet trim(	String[] markersToKeep, boolean allowIncompleteList, boolean verbose,
														Logger log) {
		String[] markerNames, newMarkerNames;
		byte[] newChrs;
		int[] newPositions;
		char[][] newAlleles;
		String[][] newAnnotation;
		boolean error;
		HashSet<String> hash, dupeHash = null;
		int count, numMissing;

		if (markersToKeep.length != Array.unique(markersToKeep).length) {
			log.reportError("Error - list of marker names to keep is not entirely unique; could cause problems with counts if not mirrored in map file");
			return null;
		}
		count = 0;
		error = false;
		numMissing = 0;
		markerNames = getMarkerNames();
		hash = HashVec.loadToHashSet(markerNames);
		if (hash.size() < markerNames.length) {
			log.reportTimeWarning((markerNames.length - hash.size())
														+ " duplicate markers detected in SnpMarkerSet");
			// reload marker hash and identify duplicate markers
			hash = new HashSet<String>((int) (1 + markerNames.length / .75));
			dupeHash = new HashSet<String>((int) (1 + (markerNames.length - hash.size()) / .75));
			for (String mkr : markerNames) {
				if (hash.contains(mkr)) {
					dupeHash.add(mkr);
				} else {
					hash.add(mkr);
				}
			}
		}
		for (String element : markersToKeep) {
			if (hash.contains(element)) {
				count++;
				if (dupeHash != null && dupeHash.contains(element)) {
					count++;
				}
			} else {
				if (!allowIncompleteList) {
					if (!error) {
						log.reportError("Error - the following markers were present in list but not in SnpMarkerSet:");
					}
					System.err.println(element);
				}
				error = true;
				numMissing++;
			}
		}
		if (numMissing > 0) {
			log.reportError("Warning - there "
												+ (numMissing == 1 ? "was 1 marker" : "were " + numMissing + " markers")
											+ " present in list but not in SnpMarkerSet");
		}
		if (!allowIncompleteList && error) {
			return null;
		}

		newMarkerNames = new String[count];
		if (chrs == null) {
			newChrs = null;
		} else {
			newChrs = new byte[count];
		}
		if (positions == null) {
			newPositions = null;
		} else {
			newPositions = new int[count];
		}
		if (alleles == null) {
			newAlleles = null;
		} else {
			newAlleles = new char[count][];
		}
		if (annotation == null) {
			newAnnotation = null;
		} else {
			newAnnotation = new String[count][];
		}
		hash = HashVec.loadToHashSet(markersToKeep);
		count = 0;
		for (int i = 0; i < markerNames.length; i++) {
			if (hash.contains(markerNames[i])) {
				if (count >= newMarkerNames.length) {
					log.reportError("Error - there are duplicate markerNames within this SnpMarkerSet which has caused trim() to fail");
					return null;
				}
				newMarkerNames[count] = markerNames[i];
				if (chrs != null) {
					newChrs[count] = chrs[i];
				}
				if (positions != null) {
					newPositions[count] = positions[i];
				}
				if (alleles != null) {
					newAlleles[count] = alleles[i];
				}
				if (annotation != null) {
					newAnnotation[count] = annotation[i];
				}
				count++;
			}
		}

		return new SnpMarkerSet(newMarkerNames, newChrs, newPositions, newAlleles, newAnnotation, false,
														verbose, log);
	}

	public static int determineType(String filename) {
		if (filename.endsWith(".bim")) {
			return PLINK_BIM_FORMAT;
		} else if (filename.endsWith(".map") || filename.endsWith(".pmap")) {
			return PLINK_MAP_FORMAT;
		} else if (filename.endsWith(".snp")) {
			return EIGENSOFT_SNP_FORMAT;
		} else if (filename.endsWith(".snps")) {
			return NAMES_ONLY;
		} else if (filename.endsWith(".info")) {
			if (Files.exists(filename, false)
					&& Files.getHeaderOfFile(filename, "[\\s]+", new Logger())[0].equals("SNP")) {
				return MACH_MLINFO_FORMAT;
			} else {
				return HAPLOVIEW_INFO_FORMAT;
			}
		} else if (filename.endsWith(".mlinfo") || filename.endsWith(".pinfo")) {
			return MACH_MLINFO_FORMAT;
		} else if (filename.endsWith(".minfo")) {
			return MINIMAC_INFO_FORMAT;
		} else if (filename.endsWith(".xln")	|| filename.endsWith(".txt") || filename.endsWith(".dat")
								|| filename.endsWith(".csv")) {
			System.err.println("Warning - assuming the map file has three columns with headers: MarkerName, Chr, Position");
			return GENERIC_FORMAT_IGNORE_FIRST_LINE;
		} else if (filename.endsWith(".burdenInfo")) {
			return INFO_FOR_BURDEN_TESTING;
		} else if (filename.endsWith(".impute2_info") || filename.endsWith(".imputed_info")) {
			return IMPUTE2_INFO_FORMAT;
		} else if (filename.endsWith(".freeze_info")) {
			return FREEZE5_FORMAT;
		} else {
			System.err.println("Warning - format of file ('"	+ filename
													+ "') could not be deduced solely by the filename extension");
			return -1;
		}
	}

	public String[] pickRandomMarkers(int numberOfMarkersToPick) {
		String[] markerNames, picks;
		int[] keys;

		markerNames = getMarkerNames();
		keys = Array.random(markerNames.length);
		picks = new String[numberOfMarkersToPick];

		for (int i = 0; i < picks.length; i++) {
			picks[i] = markerNames[keys[i]];
		}

		return picks;
	}

	public static void parseHapMap(String dir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;

		try {
			writer = new PrintWriter(new FileWriter(dir + "master.map"));
			count = 1;
			for (int chr = 1; chr <= 22; chr++) {
				try {
					reader = new BufferedReader(new FileReader(dir + "genetic_map_chr" + chr + "_b36.txt"));
					ext.checkHeader(reader.readLine().trim().split("[\\s]+"), HAPMAP_CM_SRC_HEADER, true);
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						writer.println(chr + "\trs" + count + "\t" + line[2] + "\t" + line[0]);
						count++;
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ dir + "genetic_map_chr" + chr + "_b36.txt"
															+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""	+ dir + "genetic_map_chr" + chr + "_b36.txt"
															+ "\"");
					System.exit(2);
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir + "master.map");
			e.printStackTrace();
		}
	}

	public static SnpMarkerSet merge(SnpMarkerSet... sets) {
		Hashtable<String, Integer> indices;
		String[][] arraysOfMarkerNames, annotation;
		boolean[] hasPositions, hasAlleles;
		String[] markerNames;
		byte[] chrs, travChrs;
		int[] rawPositions, travPositions;
		char[][] alleles, travAlleles;
		int index;

		arraysOfMarkerNames = new String[sets.length][];
		hasPositions = new boolean[sets.length];
		hasAlleles = new boolean[sets.length];
		for (int i = 0; i < sets.length; i++) {
			arraysOfMarkerNames[i] = sets[i].getMarkerNames();
			hasPositions[i] = sets[i].getChrs() != null;
			hasAlleles[i] = sets[i].getAlleles() != null;
		}
		annotation = Unique.proc(arraysOfMarkerNames, false);
		markerNames = Matrix.extractColumn(annotation, 0);
		annotation = Array.toMatrix(Matrix.extractColumns(annotation,
																											Array.subArray(	Array.arrayOfIndices(sets.length
																																														+ 2),
																																			1),
																											"\t"));

		indices = HashVec.loadToHashIndices(markerNames, new Logger());
		if (Array.booleanArraySum(hasPositions) > 0) {
			chrs = Array.byteArray(markerNames.length, (byte) -9);
			rawPositions = new int[markerNames.length];

			for (int i = 0; i < sets.length; i++) {
				if (hasPositions[i]) {
					travChrs = sets[i].getChrs();
					travPositions = sets[i].getPositions();
					System.out.println("Set "	+ i + " -- " + arraysOfMarkerNames[i].length + " markers, "
															+ travChrs.length + " chr values, " + travPositions.length
															+ " pos values");
					for (int j = 0; j < arraysOfMarkerNames[i].length; j++) {
						index = indices.get(arraysOfMarkerNames[i][j]);
						if (chrs[index] == (byte) -9) {
							chrs[index] = travChrs[j];
							rawPositions[index] = travPositions[j];
						} else if (travChrs[j] != 0 && (chrs[index] != travChrs[j]
																						|| rawPositions[index] != travPositions[j])) {
							System.err.println("Error - mismatched chr:position for marker "
																		+ arraysOfMarkerNames[i][j] + " (" + travChrs[j] + ":"
																	+ travPositions[j] + " but previously " + chrs[index] + ":"
																	+ rawPositions[index] + ")");
						}
					}
				}
			}

		} else {
			chrs = null;
			rawPositions = null;
		}
		if (Array.booleanArraySum(hasAlleles) > 0) {
			alleles = new char[markerNames.length][];

			for (int i = 0; i < sets.length; i++) {
				if (hasPositions[i]) {
					travAlleles = sets[i].getAlleles();
					for (int j = 0; j < arraysOfMarkerNames[i].length; j++) {
						index = indices.get(arraysOfMarkerNames[i][j]);
						if (alleles[index] == null) {
							alleles[index] = travAlleles[j];
						} else if ((travAlleles[j][0] != alleles[index][0]
												&& travAlleles[j][0] != alleles[index][1])
												|| (travAlleles[j][1] != alleles[index][0]
														&& travAlleles[j][1] != alleles[index][1])) {
							System.err.println("Error - mismatched alleles for marker "
																		+ arraysOfMarkerNames[i][j] + " (" + travAlleles[j][0] + "/"
																	+ travAlleles[j][1] + " but previously " + alleles[index][0] + "/"
																	+ alleles[index][1] + ")");
						}
					}
				}
			}
		} else {
			alleles = null;
		}

		return new SnpMarkerSet(markerNames, chrs, rawPositions, alleles, annotation, true, false,
														new Logger());
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "plink.bim";
		boolean noX = false;
		String excludeMarkers = null;
		boolean verbose = true;
		String source = "";
		SnpMarkerSet markerSet, sourceSet;
		String hapmap = "";
		String extract = null;
		String outfile = null;
		Logger log;

		String usage = "\n"	+ "filesys.SnpMarkerSet requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n"
										+ "   (2) autosomes only (i.e. -noX (not the default))\n"
										+ "   (3) list of markers to exclude (i.e. excludeMarkers=" + excludeMarkers
										+ " (default))\n" + "   (4) verbosity (i.e. verbose=" + verbose
										+ " (default))\n" + " OR\n"
										+ "   (1) interpolate centiMorgans (i.e. sourceMap=plink.bim (not the default; expecting a PLINK bim or map file))\n"
										+ " OR\n"
										+ "   (1) parse HapMap centiMorgans (i.e. parseHapMap=/home/directory/ (not the default))\n"
										+ " OR\n"
										+ "   (1) list of markers to extract (i.e. extract=list.txt (not the default))\n"
										+ "   (2) name of output file (i.e. out=chr1.pinfo (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-noX")) {
				noX = true;
				numArgs--;
			} else if (arg.startsWith("excludeMarkers=")) {
				excludeMarkers = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("verbose=")) {
				verbose = Boolean.valueOf(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("source=")) {
				source = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("parseHapMap=")) {
				hapmap = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("extract=")) {
				extract = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("out=")) {
				outfile = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		// filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\plink.bim";
		// source = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\NGRC.bim";
		// filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\maps\\Rutgers
		// smooth_map_b36\\plink.bim";
		// source = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\maps\\Rutgers
		// smooth_map_b36\\smoothRutgers.map";
		// filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\maps\\rutgers_map_b36\\plink.bim";
		// source = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\maps\\rutgers_map_b36\\rutgers.map";
		// filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\maps\\Rutgers Enhanced Linkage
		// Maps\\plink.bim";
		// source = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\analyses\\CIDR\\segments\\maps\\Rutgers Enhanced Linkage
		// Maps\\European.map";

		// filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\MasterSnpList\\centromereMidpoints.txt";
		// source = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\MasterSnpList\\allSNPs.xln";

		// hapmap = "C:\\Documents and Settings\\npankrat\\My Documents\\Downloads\\GERMLINE\\maps\\";

		// filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\Downloads\\GERMLINE\\maps\\plink.map";
		// source = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\Downloads\\GERMLINE\\maps\\master.map";

		// filename = "D:/BOSS/LinkageMergedIBC/conversion/plink.bim";
		// filename = "D:/BOSS/LinkageJustIBC/plink.bim";
		// source = "D:/BOSS/LinkageMergedIBC/conversion/European.map";

		log = new Logger();
		try {
			if (!hapmap.equals("")) {
				parseHapMap(hapmap);
			} else if (extract != null) {
				if (outfile == null) {
					outfile = ext.addToRoot(filename, ext.rootOf(extract));
				}
				new SnpMarkerSet(filename, false, log).trim(
																										HashVec.loadFileToStringArray(extract, false,
																																									new int[] {0},
																																									false),
																										true, false, log)
																							.writeToFile(	outfile,
																														SnpMarkerSet.determineType(outfile),
																														log);
			} else if (!source.equals("")) {
				if (new File(filename + ".ser").exists()) {
					markerSet = SnpMarkerSet.load(filename + ".ser", false);
				} else {
					markerSet = new SnpMarkerSet(filename, verbose, new Logger());
					markerSet.serialize(filename + ".ser");
				}
				if (new File(source + ".ser").exists()) {
					sourceSet = SnpMarkerSet.load(source + ".ser", false);
				} else {
					sourceSet = new SnpMarkerSet(source, verbose, new Logger());
					sourceSet.sortMarkers();
					sourceSet.serialize(source + ".ser");
				}
				markerSet.interpolateCentiMorgans(sourceSet, log);
				markerSet.writeToFile(ext.rootOf(filename, false)	+ "_with_centiMorgans.bim",
															determineType(filename), log);
			} else {
				new SnpMarkerSet(	filename, verbose,
													new Logger()).listUnambiguousMarkers(filename	+ "_unambiguous.txt",
																																excludeMarkers, noX);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
