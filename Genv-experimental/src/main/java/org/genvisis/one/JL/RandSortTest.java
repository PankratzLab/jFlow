package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.common.Aliases;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class RandSortTest {
	private static final int MAX_ERRORS_TO_REPORT = 10;

	public static void testEm() {
		String file = "D:/data/LLFS_GWAS/testSort/markerPositionsHG19.txt";
		orderMarkers(null, file, file + ".oradf", new Logger());

	}

	public static void orderMarkers(String[] markerNames, String markerDatabase, String output,
																	Logger log) {
		Hashtable<String, String> snpPositions;
		byte[] chrs;
		int[] positions;
		String[] line;
		Vector<String> v;
		long time;
		int[] chrCounts;

		// logLevel = log.getLevel();
		log.setLevel(9);
		time = new Date().getTime();
		log.report(ext.getTime() + "\tLoading marker data from " + markerDatabase);
		snpPositions = loadFileToHashString(markerDatabase, log);
		if (snpPositions == null) {
			return;
		}
		if (markerNames == null) {
			for (String element : Aliases.MARKER_NAMES) {
				if (snpPositions.containsKey(element)) {
					snpPositions.remove(element);
				}
			}
			log.report(ext.getTime() + "\tdone safsaf data from " + markerDatabase);
			// ArrayList<String> maTmp = new ArrayList<String>();
			// for (String ma : snpPositions.keySet()) {
			//
			// }
			markerNames = HashVec.getKeys(snpPositions, false, false);
		}
		log.report(ext.getTime() + "\tdone marker data from " + markerDatabase);

		v = new Vector<String>();
		log.report(ext.getTime() + "\tSorting markers by chromosome and position");
		chrs = new byte[markerNames.length];
		chrCounts = new int[Positions.CHR_CODES.length];
		positions = new int[markerNames.length];
		for (int i = 0; i < markerNames.length; i++) {
			if (snpPositions.containsKey(markerNames[i])) {
				line = snpPositions.get(markerNames[i]).split("[\\s]+");
				chrs[i] = Positions.chromosomeNumber(line[0], log);
				chrCounts[chrs[i]]++;
				positions[i] = Integer.parseInt(line[1]);
			} else {
				v.add(markerNames[i]);
			}
		}
		System.out.println("DHFSD");
		// long timse = System.currentTimeMillis();

		Posit[] posits = new Posit[chrs.length];
		for (int i = 0; i < posits.length; i++) {
			posits[i] = new Posit(chrs[i], positions[i], i);
		}
		Arrays.sort(posits, new PositComp());
		int[] t = new int[posits.length];
		for (int i = 0; i < t.length; i++) {
			t[i] = posits[i].getIndex();
		}
		System.out.println(ext.getTimeElapsed(time) + " for new");
		new MarkerSet(markerNames, chrs, positions, t).serialize(output);

		// timse = System.currentTimeMillis();
		Sort.orderTwoLayers(chrs, positions, log);
		System.out.println(ext.getTimeElapsed(time) + " for original");

	}

	private static class PositComp implements Comparator<Posit> {

		@Override
		public int compare(Posit o1, Posit o2) {
			int value1 = Byte.compare(o1.getChr(), o2.getChr());
			if (value1 == 0) {
				value1 = Integer.compare(o1.getPos(), o2.getPos());
			}
			return value1 == 0 ? Integer.compare(o1.getIndex(), o2.getIndex()) : value1;
		}
	}

	private static class Posit {
		final byte chr;
		final int pos;
		final int index;

		public Posit(byte chr, int pos, int index) {
			super();
			this.chr = chr;
			this.pos = pos;
			this.index = index;
		}

		public byte getChr() {
			return chr;
		}

		public int getPos() {
			return pos;
		}

		public int getIndex() {
			return index;
		}

	}

	public static Hashtable<String, String> loadFileToHashString(String filename, Logger log) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		String markerName, chr, position, delimiter, temp;
		byte chrValue;
		int count, countBad, numBlankNames, numBlankChrs, numBlankPositions, numRepeatedNames,
				numInvalidChrs, numInvalidPositions, numIncompleteLines;

		delimiter = Files.determineDelimiter(filename, log);

		count = countBad = 0;
		numBlankNames = numBlankChrs =
																	numBlankPositions =
																										numRepeatedNames = numInvalidChrs =
																																											numInvalidPositions =
																																																					numIncompleteLines = 0;
		try {
			reader = Files.getAppropriateReader(filename);
			while (reader.ready()) {
				temp = reader.readLine();
				if (delimiter.equals(",")) {
					line = ext.splitCommasIntelligently(temp, true, new Logger());
				} else if (temp.indexOf("\t") == -1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				if (count == 0 && ext.indexOfStr(line[0], Aliases.MARKER_NAMES) >= 0) {

				} else if (line.length < 3) {
					if (countBad < MAX_ERRORS_TO_REPORT) {
						log.report("Error - incomplete line at row "	+ count + " for marker \"" + line[0]
												+ "\"; line will not be added");
					}
					numIncompleteLines++;
				} else {
					markerName = line[0];
					chr = line[1];
					position = line[2];

					if (markerName.equals("")) {
						if (countBad < MAX_ERRORS_TO_REPORT) {
							log.reportError("Error - blank marker name at line " + count + " of " + filename);
						}
						numBlankNames++;
						countBad++;
					} else if (chr.equals("")) {
						if (countBad < MAX_ERRORS_TO_REPORT) {
							log.reportError("Error - blank chr for marker '"	+ markerName + "' at line " + count
															+ " of " + filename);
						}
						numBlankChrs++;
						countBad++;
					} else if (position.equals("")) {
						if (countBad < MAX_ERRORS_TO_REPORT) {
							log.reportError("Error - blank position for marker '"	+ markerName + "' at line "
															+ count + " of " + filename);
						}
						numBlankPositions++;
						countBad++;
					} else {
						if (hash.containsKey(markerName)) {
							log.reportError("Error - marker '"	+ markerName
															+ "' is already listed in the markerPositions file and is seen again at line "
															+ count + " of " + filename);
							numRepeatedNames++;
							countBad++;
						}
						chrValue = Positions.chromosomeNumber(chr, log);
						if (chrValue < 0 || chrValue > 26) {
							numInvalidChrs++;
							countBad++;
						}
						try {
							Integer.parseInt(position);
						} catch (NumberFormatException nfe) {
							if (countBad < 10) {
								log.reportError("Error - invalid position ("	+ position + ") for marker '"
																+ markerName + "' at line " + count + " of " + filename);
							}
							numInvalidPositions++;
							countBad++;
						}
					}

					hash.put(markerName, chr + "\t" + position);
				}
				if (countBad == 10) {
					log.reportError("...");
					countBad++;
				}
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		log.report("\nRead in " + ext.addCommas(count) + " lines from the markerPositions file");
		if (countBad > 10) {
			countBad--;
		}

		if (countBad > 0) {
			log.report("...with a total of "	+ ext.addCommas(countBad) + " problem"
									+ (countBad == 1 ? "" : "s"));
		}
		if (numIncompleteLines > 0) {
			log.report("...including "	+ ext.addCommas(numIncompleteLines) + " incomplete line"
									+ (numIncompleteLines == 1 ? "" : "s"));
		}
		log.report("Number of final valid marker positions: " + ext.addCommas(hash.size()));
		if (numBlankNames > 0) {
			log.report("Number of blank marker names: " + ext.addCommas(numBlankNames));
		}
		if (numBlankChrs > 0) {
			log.report("Number of blank chromosomes: " + ext.addCommas(numBlankChrs));
		}
		if (numBlankPositions > 0) {
			log.report("Number of blank marker positions: " + ext.addCommas(numBlankPositions));
		}
		if (numRepeatedNames > 0) {
			log.report("Number of repeated marker names: " + ext.addCommas(numRepeatedNames));
		}
		if (numInvalidChrs > 0) {
			log.report("Number of invalid chromosomes: " + ext.addCommas(numInvalidChrs));
		}
		if (numInvalidPositions > 0) {
			log.report("Number of invalid positions: " + ext.addCommas(numInvalidPositions));
		}
		log.report("");

		return hash;
	}

	public static void main(String[] args) {
		testEm();
	}

}
