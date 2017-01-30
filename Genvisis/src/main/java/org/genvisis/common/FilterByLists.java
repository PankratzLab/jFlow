package org.genvisis.common;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

public class FilterByLists {
	public static void fromParameters(String controlFile, Logger log) {
		String[] line;
		String filename, outputFilename;
		int col;
		Vector<String> paramV;
		boolean header;
		String keeps, deletes;
		boolean commaDelimited;

		paramV = Files.parseControlFile(controlFile, "filterByLists",
																		new String[] {"fileToFilter.txt header 0 out=file.out",
																									"keeps.txt", "deletes.txt"},
																		log);
		if (paramV == null) {
			return;
		}

		line = paramV.remove(0).trim().split("[\\s]+");
		filename = line[0];
		header = false;
		commaDelimited = false;
		col = -1;
		outputFilename = ext.rootOf(filename) + "_filtered.out";
		for (int j = 1; j < line.length; j++) {
			if (line[j].equals("header")) {
				header = true;
			} else if (line[j].equals(",")) {
				commaDelimited = true;
			} else if (line[j].startsWith("out=")) {
				outputFilename = line[j].split("=")[1];
			} else if (col == -1) {
				col = Integer.parseInt(line[j]);
			} else {
				System.err.println("Error - what am I supposed to do with '" + line[j] + "'?");
			}
		}
		if (col == -1) {
			System.err.println("Warning - assuming unique id can be found in column 0");
			col = 0;
		}
		log.report("Loading data from '" + filename + "'");
		keeps = paramV.remove(0);
		if (keeps.equalsIgnoreCase("null")) {
			keeps = null;
		}
		deletes = paramV.remove(0);
		if (deletes.equalsIgnoreCase("null")) {
			deletes = null;
		}

		process(filename, keeps, deletes, new int[] {col}, outputFilename, header, true, commaDelimited,
						log);
	}

	private static Hashtable<String, String> loadFileToHash(String filename, String type) {
		Hashtable<String, String> hash;
		String listFile;
		int[] listCols;

		if (filename == null) {
			hash = null;
		} else {
			if (filename.contains(";")) {
				listFile = filename.substring(0, filename.indexOf(";"));
				listCols = ArrayUtils.toIntArray(filename.substring(filename.indexOf(";") + 1).split(","));
			} else {
				listFile = filename;
				listCols = new int[] {0};
			}
			if (!new File(listFile).exists()) {
				System.err.println("Since '"	+ filename
														+ "' is not a filename, assuming this is the key to " + type);
				hash = new Hashtable<String, String>();
				hash.put(filename, "1");
			} else {
				hash = HashVec.loadFileToHashString(listFile, listCols, new int[] {-7}, false, null, false,
																						false, false);
			}
		}

		return hash;
	}

	public static void process(	String filename, String keeps, String deletes, int[] cols,
															String outfile, boolean keepFirstLine, boolean reportMissingElements,
															boolean commaDelimited, Logger log) {
		process(filename, cols, outfile, loadFileToHash(keeps, "keep"),
						loadFileToHash(deletes, "delete"), keepFirstLine, true, reportMissingElements,
						commaDelimited, log);
	}

	public static void process(	String filename, int[] cols, String outfile,
															Hashtable<String, String> keepsHash,
															Hashtable<String, String> deletesHash, boolean keepFirstLine,
															boolean checkForOverlap, boolean reportMissingElements,
															boolean commaDelimited, Logger log) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String temp, key;
		boolean isFirstLine;
		String[] keepKeys;
		String[] delKeys;
		Vector<String> v;
		IntVector iv;
		int[] order;

		if (keepsHash != null && deletesHash != null) {
			System.err.println("Error - you have specified both what to keep and what to delete. This is redudant if they are mirror opposites, and if there is overlap between the two, there is no way to decide which takes precendence. Use only one or the other.");
			return;
		}

		if (reportMissingElements && keepsHash != null) {
			keepKeys = HashVec.getKeys(keepsHash);
		} else {
			keepKeys = null;
		}

		if (reportMissingElements && deletesHash != null) {
			delKeys = HashVec.getKeys(deletesHash);
		} else {
			delKeys = null;
		}

		isFirstLine = true;
		if (outfile == null) {
			new File(filename).renameTo(new File(filename + ".bak"));
			outfile = filename;
			filename = filename + ".bak";
		}

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(outfile));
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split(commaDelimited ? "," : "[\\s]+");
				key = ArrayUtils.toStr(ArrayUtils.subArray(line, cols));
				if (isFirstLine && keepFirstLine) {
					writer.println(temp);
					isFirstLine = false;
				} else if (deletesHash != null && deletesHash.containsKey(key)) {
					deletesHash.put(key, "");
				} else if (keepsHash == null || keepsHash.containsKey(key)) {
					writer.println(temp);
					if (keepsHash != null) {
						keepsHash.put(key, "");
					}
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		if (reportMissingElements) {
			if (keepsHash != null) {
				v = new Vector<String>();
				for (int i = 0; i < keepKeys.length; i++) {
					if (!keepsHash.get(keepKeys[i]).equals("")) {
						v.add(keepKeys[i]);
					}
				}
				if (v.size() > 0) {
					System.err.println("Warning - the following were found in the keeps list but not in the data file:");
					Collections.sort(v);
					for (int i = 0; i < v.size(); i++) {
						System.err.println(v.get(i));
					}
				}
			}

			if (deletesHash != null) {
				v = new Vector<String>();
				for (int i = 0; i < delKeys.length; i++) {
					if (!deletesHash.get(delKeys[i]).equals("")) {
						v.add(delKeys[i]);
					}
				}
				if (v.size() > 0) {
					System.err.println("Warning - the following were found in the deletes list but not in the data file:");
					Collections.sort(v);
					for (int i = 0; i < v.size(); i++) {
						System.err.println(v.get(i));
					}
				}
			}
		}
	}



	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String keepsFile = null;
		String deletesFile = null;
		String outfile = null;
		boolean keepFirstLine = false;
		boolean reportMissingElements = true;
		int col = 0;
		boolean commaDelimited = false;

		String usage = "\n"	+ "common.FilterByLists requires 0-1 arguments\n"
										+ "   (1) input filename (i.e. file=" + filename + " (default))\n"
										+ "   (2) filename containing list of tokens to keep (i.e. keeps=keeps.txt (not the default; if token is not a filename, assuming it is the key to keep))\n"
										+ "   (3) filename containing list of tokens to remove (i.e. deletes=deltetes.txt (not the default; if token is not a filename, assuming it is the key to delete))\n"
										+ "   (4) index of the column of the token to match on (i.e. col=" + col
										+ " (default))\n"
										+ "   (5) output filename (i.e. out=output.out (not the default; default is to backup and replace the input file))\n"
										+ "   (6) always keep first line (i.e. -keepFirst (not the default))\n"
										+ "   (7) comma delimited file (i.e. -commaDelimited (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("keeps=")) {
				keepsFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("deletes=")) {
				deletesFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("col=")) {
				col = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("out=")) {
				outfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-keepFirst")) {
				keepFirstLine = true;
				numArgs--;
			} else if (arg.startsWith("-commaDelimited")) {
				commaDelimited = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			process(filename, keepsFile, deletesFile, new int[] {col}, outfile, keepFirstLine,
							reportMissingElements, commaDelimited, new Logger());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
