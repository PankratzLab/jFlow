package org.genvisis.one;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class VTE_Analyses {
	public static String DIR = "D:/_transfer/umn/Folson/VTE_meta_analysis/finalAnalysis/";
	// public static final String DIR = "";

	public static final String IMPUTATION_MAP = "/home/npankrat/NCBI/1000G/EUR.map";

	public static final String[][] STUDIES = {{"ARIC", "ARIC_autosomes_table_withUpdatedChr2.xln",
																						 "\t"},
																						{"CHS", "GH-VTE-results.csv", ","},
																						{"GH", "GH-VTE-results.csv", ","},
																						{"RS1", "RS1VTE_CHARGE.chargefmt.RS.txt",
																						 PSF.Regex.GREEDY_WHITESPACE},
																						{"RS2", "RS2VTE_CHARGE.chargefmt.RS.txt",
																						 PSF.Regex.GREEDY_WHITESPACE},
																						{"WGHS", "WGHS_incident_VTE_oct_20_2009.txt",
																						 PSF.Regex.GREEDY_WHITESPACE},};

	public static void generateSuperMap() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		Logger log;
		String[] trav, ref;

		log = new Logger(DIR + "allSNPsMap.log");
		for (String[] element : STUDIES) {
			line = Files.getHeaderOfFile(DIR + element[1], element[2], log);
			// System.out.println(STUDIES[i][0]+"\t"+Array.toStr(line));
			ext.checkHeader(line, new String[] {"SNPID", "chr", "position"}, new int[] {0, 1, 2}, true,
											log, true);
		}
		for (String[] element : STUDIES) {
			System.out.println(element[0]);
			try {
				reader = new BufferedReader(new FileReader(DIR + element[1]));
				reader.readLine();
				while (reader.ready()) {
					line = reader.readLine().trim().split(element[2]);
					trav = new String[] {line[0], line[1], line[2]};
					if (hash.containsKey(line[0])) {
						ref = hash.get(line[0]);
						if (!ref[2].equals("NA") && (!line[1].equals(ref[1]) || !line[2].equals(ref[2]))) {
							log.reportError("Error - mismatched positions for marker '" + line[0] + "' (was "
															+ ref[1] + ":" + ref[2] + ", now " + line[1] + ":" + line[2] + ")");
						}
					}

					hash.put(line[0], trav);
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + DIR + element[1]
													 + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + DIR + element[1] + "\"");
				System.exit(2);
			}
		}

		String[] keys;
		byte[] chrs;
		int[] positions, order;

		keys = HashVec.getKeys(hash);
		chrs = new byte[keys.length];
		positions = new int[keys.length];
		for (int i = 0; i < keys.length; i++) {
			trav = hash.get(keys[i]);
			chrs[i] = Positions.chromosomeNumber(trav[1]);
			positions[i] = ext.isMissingValue(trav[2]) ? -1 : Integer.parseInt(trav[2]);
		}
		log.report(ext.getTime() + "\tSorting positions...");
		order = Sort.getSort2DIndices(chrs, positions);
		log.report(ext.getTime() + "\tWriting to file...");
		try {
			writer = Files.openAppropriateWriter(DIR + "allSNPs.map");
			writer.println("MarkerName\tChr\tPosition");
			for (int i = 0; i < keys.length; i++) {
				writer.println(keys[order[i]] + "\t" + chrs[order[i]] + "\t" + positions[order[i]]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + DIR + "allSNPs.map");
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		try {
			generateSuperMap();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
