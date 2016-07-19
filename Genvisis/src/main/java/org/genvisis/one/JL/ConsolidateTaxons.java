package one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import common.Array;
import common.Files;
import common.Logger;
import common.ext;

public class ConsolidateTaxons {
	public static void consolidate(String rootDir) {
		String[] dirs = Files.listDirectories(rootDir, false);
		Logger log = new Logger();
		ArrayList<Hashtable<String, String[]>> taxa = new ArrayList<Hashtable<String, String[]>>();
		HashSet<String> allTaxa = new HashSet<String>();
		String output = rootDir + "taxaSummary.txt";
		log.reportTimeInfo("NUMDIRS = " + dirs.length);
		for (int i = 0; i < dirs.length; i++) {
			String[] contams = Files.listFullPaths(rootDir + dirs[i] + "/", "", false);
			log.reportTimeInfo("Current directory " + rootDir + dirs[i] + " Number of files " + contams.length);

			for (int j = 0; j < contams.length; j++) {

				if (contams[j].endsWith("contam")) {
					log.reportTimeInfo("Consolidating " + contams[j]);
					Hashtable<String, String[]> current = new Hashtable<String, String[]>();
					try {
						BufferedReader reader = Files.getAppropriateReader(contams[j]);
						while (reader.ready()) {
							String[] line = reader.readLine().trim().split("\t");
							allTaxa.add(line[0]);
							current.put(line[0], Array.subArray(line, 1, line.length));

							// System.out.println(Array.toStr(current.get(line[0])));
						}
						reader.close();
						taxa.add(current);
					} catch (FileNotFoundException fnfe) {
						log.reportError("Error: file \"" + contams[j] + "\" not found in current directory");
						return;
					} catch (IOException ioe) {
						log.reportError("Error reading file \"" + contams[j] + "\"");
						return;
					}
				}
			}
		}

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("Taxa");
			for (int i = 0; i < taxa.size(); i++) {
				String[] files = taxa.get(i).get("Taxa");
				for (int j = 0; j < files.length; j++) {
					writer.print("\t" + ext.rootOf(files[j]));
				}
			}
			writer.println();
			for (String ataxa : allTaxa) {
				if (!ataxa.equals("Taxa")) {
					ArrayList<Integer> counts = new ArrayList<Integer>();
					for (int i = 0; i < taxa.size(); i++) {
						int[] blankCounts = new int[taxa.get(i).get("Taxa").length];
						Arrays.fill(blankCounts, 0);
						String[] blanks = new String[taxa.get(i).get("Taxa").length];
						Arrays.fill(blanks, "0");
						if (taxa.get(i).containsKey(ataxa)) {
							String[] tmp = taxa.get(i).get(ataxa);
							for (int j = 0; j < tmp.length; j++) {
								counts.add( Integer.parseInt(tmp[j]));
							}
						} else {
							for (int j = 0; j < blankCounts.length; j++) {
								counts.add( blankCounts[j]);

							}
						}
					}
					int[] allCounts =Array.toIntArray(counts);
					if(Array.max(allCounts)>100){
						writer.println(ataxa+"\t"+Array.toStr(allCounts));
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	public static void main(String[] args) {
		//int numArgs = args.length;
		String rootDir = "/home/tsaim/shared/Project_Tsai_Project_021/contamination/";
		consolidate(rootDir);
	}
}
