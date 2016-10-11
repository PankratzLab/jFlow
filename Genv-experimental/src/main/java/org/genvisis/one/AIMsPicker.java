package org.genvisis.one;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class AIMsPicker {
	private static void pick(String filename, int minDistance, double minScore, Logger log) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, Vector<String>> clusters;
		Hashtable<String, String[]> allInfo;
		Vector<String> order, cluster;
		String[] picks;
		byte[] chrs;
		byte chr;
		int[] positions;
		int position;
		boolean candidate;
		int count, round, randomPick, randomFail, pick;

		order = new Vector<String>();
		clusters = new Hashtable<String, Vector<String>>();
		allInfo = new Hashtable<String, String[]>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, new String[] {"Cluster", "Variant", "region", "gene", "pval", "Chr",
																					"Position", "Score"},
											true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!allInfo.containsKey(line[1]) && Double.parseDouble(line[7]) > minScore) {
					HashVec.addIfAbsent(line[0], order);
					HashVec.addToHashVec(clusters, line[0], line[1], true);
					allInfo.put(line[1], line);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		round = count = 0;
		picks = new String[order.size()];
		chrs = new byte[picks.length];
		positions = new int[picks.length];
		while (count != picks.length && round < 100) {
			for (int i = 0; i < picks.length; i++) {
				if (picks[i] == null) {
					cluster = clusters.get(order.elementAt(i));
					candidate = true;
					pick = -1;
					randomFail = -1;
					randomPick = (int) Math.floor(Math.random() * cluster.size());
					for (int j = 0; j < cluster.size(); j++) {
						line = allInfo.get(cluster.elementAt(j));
						chr = Byte.parseByte(line[5]);
						position = Integer.parseInt(line[6]);
						candidate = true;
						for (int k = 0; candidate && k < picks.length; k++) {
							if (chr == chrs[k] && Math.abs(position - positions[k]) < minDistance) {
								candidate = false;
								if (j == randomPick) {
									randomFail = k;
								}
							}
						}
						if (candidate) {
							pick = j;
							j = cluster.size();
						}
					}
					if (pick == -1) {
						picks[randomFail] = null;
						chrs[randomFail] = 0;
						positions[randomFail] = 0;
						count--;
						pick = randomPick;
					}
					line = allInfo.get(cluster.elementAt(pick));
					chr = Byte.parseByte(line[5]);
					position = Integer.parseInt(line[6]);
					picks[i] = line[1]	+ "\t" + chr + "\t" + position + "\t" + line[7] + "\t" + line[0] + "\t"
											+ (pick + 1) + " of " + cluster.size();
					chrs[i] = chr;
					positions[i] = position;
					count++;

					// if (!candidate) {
					// picks[i] = "none-valid\t-1\t-1\t-1\t"+order.elementAt(i)+"\ttried all "+cluster.size();
					// chrs[i] = -1;
					// positions[i] = -1;
					// }
				}
			}
			System.out.print(count + " of " + picks.length);
			for (int i = 0; i < picks.length; i++) {
				if (picks[i] == null) {
					System.out.print(" " + i);
				}
			}
			System.out.println();
			round++;
		}

		Files.writeArray(picks, ext.rootOf(filename, false) + "_picks.xln");
		System.out.println("Filled " + count + " of " + picks.length + " slots");
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "AIMsPicker.dat";
		String logfile = null;
		Logger log;

		String usage = "\n"	+ "one.AIMsPicker requires 0-1 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
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
			log = new Logger(logfile);
			filename = "D:/Myron/Indian_Diabetes/SequencingPilot/SingaporeSelections_Designed/aims/ALL_15000_AIMs_pos.txt";
			pick(filename, 500000, 1.0, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
