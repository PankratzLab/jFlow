package org.genvisis.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.filesys.SerialHash;

public class ConvertToBuild37 {
	public static final String HASH_MAP_FILE = "mappings.hash.ser";

	public static void convert(String dir, String suffix) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, lines;
		String temp;
		Hashtable<String, String> hash, err, newLocs;
		// Vector<String> v = new Vector<String>();
		int count;
		// long time;
		String[] files, locs, conversions;
		byte[] chrs;
		int[] order, positions, pos;
		String loc;

		if (new File(dir + HASH_MAP_FILE).exists()) {
			hash = SerialHash.loadSerializedStringHash(dir + HASH_MAP_FILE);
		} else {
			hash = new Hashtable<String, String>();
		}

		err = new Hashtable<String, String>();
		files = Files.list(dir, ".err", false);
		for (String file : files) {
			lines = HashVec.loadFileToStringArray(dir	+ file, false, false, new int[] {0}, true, false,
																						",");
			for (int j = 0; j < lines.length; j++) {
				if (lines[j].equals("#Deleted in new")) {
					err.put(lines[j + 1], "");
					j++;
				} else {
					System.err.println("Error - unknown error code: " + lines[j]);
					System.exit(1);
				}
			}
		}

		if (new File(dir + "liftOver.in").exists() && new File(dir + "liftOver.bed").exists()) {
			locs = HashVec.loadFileToStringArray(dir + "liftOver.in", false, new int[] {0}, false);
			conversions =
									HashVec.loadFileToStringArray(dir + "liftOver.bed", false, new int[] {0}, false);
			if (locs.length == conversions.length) {
				for (int i = 0; i < locs.length; i++) {
					hash.put(locs[i], conversions[i]);
				}
				new File(dir + "liftOver.bed").renameTo(new File(dir + "liftOver.bed.proccessed"));
				SerialHash.createSerializedStringHash(dir + HASH_MAP_FILE, hash);
				System.out.println("Successfully processed liftOver.bed");
			} else {
				System.out.println("Mismatched list-conversion lengths, assuming there are new errors that need to be parsed; parsing now...");
			}
		}

		count = 0;
		files = Files.list(dir, suffix, false);
		newLocs = new Hashtable<String, String>();
		for (String file : files) {
			try {
				reader = new BufferedReader(new FileReader(dir + file));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (!line[0].startsWith("#")) {
						loc = Positions.getUCSCformat(new String[] {line[1], line[2], line[2]});
						if (err.containsKey(loc)) {
							count++;
						} else {
							if (!hash.containsKey(loc)) {
								newLocs.put(line[1] + "\t" + line[2] + "\t" + line[2], "");
							}
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + dir + file + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + dir + file + "\"");
				System.exit(2);
			}
		}
		System.out.println();

		if (newLocs.size() > 0) {
			locs = HashVec.getKeys(newLocs, true, false);
			chrs = new byte[locs.length];
			positions = new int[locs.length];
			for (int i = 0; i < positions.length; i++) {
				line = locs[i].split("[\\s]+");
				chrs[i] = Positions.chromosomeNumber(line[0]);
				positions[i] = Integer.parseInt(line[1]);
				locs[i] = Positions.getUCSCformat(line);
			}
			order = Sort.orderTwoLayers(chrs, positions, new Logger());
			locs = Sort.putInOrder(locs, order);
			if (new File(dir + "liftOver.in").exists()
					&& Array.equals(locs, HashVec.loadFileToStringArray(dir	+ "liftOver.in", false,
																															new int[] {0}, false),
													false)) {
				System.err.println("Error - getting the same list of positions to convert and it's not the same number as is in liftOver.bed");
			} else {
				Files.writeArray(locs, dir + "liftOver.in");
				System.out.println("Found "	+ locs.length + " new positions that have yet to be parsed:\n"
														+ "     1) upload file 'liftOver.in' to liftOver on UCSC\n"
														+ "     2) copy any .err files to the current directory\n"
														+ "     3) copy the bed file to the current directory and rename to liftOver.bed"
														+ "     4) re-run once more (you'll still get this message if there were new errors being filtered out, simply re-run");
			}
			return;
		}

		System.out.println("Found no new positions that have yet to be parsed; generating new files for SeattleSeq");
		count = 0;
		new File(dir + "inputs/").mkdirs();
		for (String file : files) {
			try {
				reader = new BufferedReader(new FileReader(dir + file));
				writer = new PrintWriter(new FileWriter(dir	+ "inputs/"
																								+ file.substring(0, file.lastIndexOf("_"))
																								+ "_input.txt"));
				while (reader.ready()) {
					temp = reader.readLine();
					line = temp.trim().split("[\\s]+");
					if (!line[0].startsWith("#")) {
						if (line.length < 5) {
							System.err.println("Error - invalid line in '" + file + "': " + temp);
						}
						loc = Positions.getUCSCformat(new String[] {line[1], line[2], line[2]});
						pos = Positions.parseUCSClocation(hash.get(loc));
						if (err.containsKey(loc)) {
							count++;
						} else {
							if (pos.length < 3) {
								System.err.println("Error - invalid positions: " + Array.toStr(pos));
							}
							writer.println("chr"	+ pos[0] + "\t" + pos[1] + "\t" + pos[2] + "\t" + line[3] + "\t"
															+ line[4]);
						}
					}
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + dir + file + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + dir + file + "\"");
				System.exit(2);
			}
		}
		System.out.println("Ignored " + count + " positions not found on build 37.1");



	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "D:\\tWork\\2qSequencing\\OnTarget\\";
		String suffix = ".txt";

		String usage = "\n"	+ "seq.ConvertToBuild37 requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + dir + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				dir = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			convert(dir, suffix);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
