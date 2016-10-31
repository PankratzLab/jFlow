package org.genvisis.kaput;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CountVector;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.park.CheckIDsAgainstDNAs;

public class CleanupSNPresults {
	public static final String[][] REQS = {	{"plate"}, {"well"}, {"dna"}, {"FamInd"}, {"result"},
																					{"call"}};

	public static final String[] MISSING_CALL_CODES = {"undetermined", "lost"};

	public static final String[] MISSING_GENOTYPE_CODES = {"1"};

	public static final String DEFAULT_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\";

	// public static final String DEFAULT_FILENAME = "G2019S_fixup.txt";
	public static final String DEFAULT_FILENAME = "G2019S_fixup_clean.xln";

	public static void clean(String dir, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, checkLine, hashKeys;
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		CountVector vicfam = new CountVector();
		Vector<String> v;
		int best;
		int[] indices, keys;
		CheckIDsAgainstDNAs check = new CheckIDsAgainstDNAs();
		String call;

		try {
			reader = new BufferedReader(new FileReader(dir + filename));
			indices = ext.indexFactors(	REQS, reader.readLine().trim().split("\t", -1), false, false, true,
																	true);
			while (reader.ready()) {
				line = reader.readLine().toUpperCase().split("\t", -1);
				checkLine = check.checkPair(line[indices[3]], line[indices[2]], false).split("[\\s]+");
				if (checkLine[0].equals("yearbug")) {
					// line[indices[2]] = checkLine[1];
					System.err.println("  " + Array.toStr(checkLine));

				}
				if (line[indices[0]].trim().equals("")) {
					line[indices[0]] = "unlisted";
				}
				HashVec.addToHashVec(	hash,
															line[indices[3]],
															line[indices[2]]	+ "\t" + line[indices[4]] + "\t" + line[indices[5]]
																								+ "\t" + line[indices[0]] + "\t" + line[indices[1]],
															false);
				vicfam.add(line[indices[4]] + " (" + line[indices[5]] + ")");
			}
			reader.close();

			writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename) + "_clean.xln"));
			writer.println("FamInd\tDNA\tResult\tCall\tPlate\tWell\t2nd_DNA\t2nd_Result\t2nd_Call\t2nd_Plate\t2nd_Well\t...");
			hashKeys = HashVec.getKeys(hash);
			keys = Sort.getSortedIndices(Array.toIntArray(hashKeys));
			for (int i = 0; i < hashKeys.length; i++) {
				v = hash.get(hashKeys[keys[i]]);
				best = -1;
				call = "nada";
				for (int j = 0; j < v.size(); j++) {
					line = v.elementAt(j).split("[\\s]+");
					if (ext.indexOfStr(line[1].toLowerCase(), MISSING_CALL_CODES) >= 0) {

					} else {
						if (call.equals("nada")) {
							call = line[1] + "/" + line[2];
						} else if (!call.equals(line[1] + "/" + line[2])) {
							System.err.println("Error - discrepant calls for "	+ hashKeys[keys[i]] + " (" + call
																	+ " and " + line[1] + "/" + line[2] + ")");
						}
						checkLine = check.checkPair(hashKeys[keys[i]], line[0], false).split("[\\s]+");
						if (!checkLine[0].equals("yearbug")) {
							best = j;
						}
					}
				}
				if (best >= 0) {
					v.insertElementAt(v.remove(best), 0);
				}
				writer.println(hashKeys[keys[i]] + "\t" + Array.toStr(Array.toStringArray(v)));
			}
			writer.close();

			for (int i = 0; i < vicfam.getCounts().length; i++) {
				System.out.println(vicfam.getCounts()[i] + "\t" + vicfam.getValues()[i]);
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + filename + "\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = DEFAULT_DIRECTORY;
		String filename = DEFAULT_FILENAME;

		String usage = "\\n"	+ "kaput.CleanupSNPresults requires 0-1 arguments\n"
										+ "   (1) directory (i.e. dir=" + dir + " (default))\n"
										+ "   (2) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
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
			clean(dir, filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
