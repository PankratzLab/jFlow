package org.genvisis.link.init;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

public class checkForTwins {
	public checkForTwins() throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String temp, chrome, trav, prev, id, famInfo, matchup;
		StringTokenizer st, famSt;
		int numInFam, A1, A2, B1, B2;
		int[] counts, geno;
		Vector<String> members, matchups = new Vector<String>();
		Hashtable<String, int[]> hash = new Hashtable<String, int[]>();
		int[][] genos;

		for (int chromosome = 1; chromosome <= 23; chromosome++) {
			chrome = (chromosome < 10) ? "0" + chromosome : "" + chromosome;
			try {
				reader = new BufferedReader(new FileReader("chrom" + chrome + ".pre"));
				prev = famInfo = "";
				numInFam = 0;
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					trav = st.nextToken();
					id = st.nextToken();
					if (!trav.equals(prev) || !reader.ready()) {
						if (!reader.ready() && Integer.valueOf(id).intValue() < 100) {
							famInfo += temp + "#";
							numInFam++;
						}

						genos = new int[numInFam][];
						famSt = new StringTokenizer(famInfo, "#");
						members = new Vector<String>();
						for (int i = 0; i < numInFam; i++) {
							st = new StringTokenizer(famSt.nextToken());
							members.add(st.nextToken() + "#" + st.nextToken());
							st.nextElement();
							st.nextElement();
							st.nextElement();
							st.nextElement();
							geno = new int[st.countTokens()];
							for (int j = 0; j < geno.length; j++) {
								geno[j] = Integer.valueOf(st.nextToken()).intValue();
							}
							genos[i] = geno;
						}
						for (int i = 0; i < numInFam; i++) {
							for (int j = i + 1; j < numInFam; j++) {
								matchup = members.elementAt(i) + "&" + members.elementAt(j);
								if (hash.containsKey(matchup)) {
									counts = hash.get(matchup);
								} else {
									counts = new int[2];
									hash.put(matchup, counts);
									matchups.add(matchup);
								}
								for (int k = 0; k < genos[1].length; k = k + 2) {
									A1 = genos[i][k];
									A2 = genos[i][k + 1];
									B1 = genos[j][k];
									B2 = genos[j][k + 1];

									if (A1 != 0 && A2 != 0 && B1 != 0 && B2 != 0) {
										if (A1 == B1) {
											counts[0]++;
										} else {
											if (A1 == B2) {
												counts[0]++;
											} else {
												counts[1]++;
											}
										}

										if (A2 == B1) {
											counts[0]++;
										} else {
											if (A2 == B2) {
												counts[0]++;
											} else {
												counts[1]++;
											}
										}

									}
								}
							}
						}

						if (Integer.valueOf(id).intValue() < 100) {
							famInfo = temp;
							numInFam = 1;
						} else {
							famInfo = "";
							numInFam = 0;
						}
						prev = trav;
					}
					if (Integer.valueOf(id).intValue() < 100) {
						famInfo += temp + "#";
						numInFam++;
					}
				}
				reader.close();
			} catch (IOException ioe) {
				System.err.println("Error processing file: chrom" + chrome + ".pre.");
			}

		}

		writer = new PrintWriter(new FileWriter("twinCheck.xls"));
		writer.println("FamID\tInd1\tInd2\t#Same\t#Diff\tRatio");
		for (int i = 0; i < matchups.size(); i++) {
			matchup = matchups.elementAt(i);
			counts = hash.get(matchup);
			famSt = new StringTokenizer(matchup, "&");
			st = new StringTokenizer(famSt.nextToken(), "#");
			writer.print(st.nextToken() + "\t" + st.nextToken());
			st = new StringTokenizer(famSt.nextToken(), "#");
			st.nextToken();
			writer.println("\t"	+ st.nextToken() + "\t" + counts[0] + "\t" + counts[1] + "\t"
											+ ((double) counts[0] / (double) counts[1]));
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length > 0) {
			System.out.println("Expecting no arguments.");
		}
		try {
			new checkForTwins();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
