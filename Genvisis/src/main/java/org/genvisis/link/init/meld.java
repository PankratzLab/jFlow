// this contained a Hashtable[] to Vector<Hashtable> conversion, beware
package org.genvisis.link.init;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class meld {
	public meld() throws IOException {
		BufferedReader[] readers = new BufferedReader[2];
		PrintWriter writer = null, dupes, mismatchedMarkers, mismatchGenotypes;
		String temp, id;
		Vector<String> individuals = new Vector<String>();
		int[] store;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Hashtable<String, String> differs = new Hashtable<String, String>();
		Vector<Hashtable<String, String>> offs;
		String[] line1, line2, markerNames;
		int diffs = 0, sames = 0;
		int off, offCount, totalOffCount, maxOffCount, maxOff;

		dupes = Files.openAppropriateWriter("melded_dupes.out");
		mismatchedMarkers = Files.openAppropriateWriter("mismatched_markers.out");
		mismatchGenotypes = Files.openAppropriateWriter("mismatched_genotypes.out");

		for (int chr = 1; chr <= 23; chr++) {
			writer = Files.openAppropriateWriter("chromosome" + chr + ".dat");
			try {
				individuals.removeAllElements();
				hash.clear();
				try {
					readers[0] = new BufferedReader(new FileReader("1/chromosome" + chr + ".dat"));
					readers[1] = new BufferedReader(new FileReader("2/chromosome" + chr + ".dat"));
				} catch (FileNotFoundException fnfe) {
					continue;
				}
				writer.println(readers[0].readLine());
				readers[1].readLine();
				temp = readers[0].readLine();
				line1 = temp.split("[\t]+");
				line2 = readers[1].readLine().split("[\t]+");
				if (line1.length != line2.length) {
					if (line1.length > line2.length) {
						System.err.println("Error - 1/chromosome" + chr + ".dat has "
															 + (line1.length - line2.length)
															 + " more columns (markers?) than 2/chromosome" + chr + ".dat");
					} else {
						System.err.println("Error - 2/chromosome" + chr + ".dat has "
															 + (line2.length - line1.length)
															 + " more columns (markers?) than 1/chromosome" + chr + ".dat");
					}
					System.exit(1);
				}
				store = new int[line1.length - 3];
				offs = new Vector<Hashtable<String, String>>();
				for (int i = 0; i < line1.length - 3; i++) {
					offs.add(new Hashtable<String, String>());
				}
				markerNames = line1;

				for (int j = 1; j < line1.length; j++) {
					if (!line1[j].equals(line2[j])) {
						System.err.println("Headers were not the same for chromosome " + chr);
						System.exit(2);
					}
				}

				System.out.println("Merging chromosome " + chr);
				writer.println(temp);

				for (int j = 0; j < 2; j++) {
					while (readers[j].ready()) {
						temp = readers[j].readLine();
						line1 = temp.split(PSF.Regex.GREEDY_WHITESPACE);
						id = ext.formNum(line1[1], line1[2], 5, 3);
						if (hash.containsKey(id)) {
							if (chr == 1) {
								dupes.println(id);
							}
							line1 = temp.split(PSF.Regex.GREEDY_WHITESPACE);
							line2 = (hash.get(id)).split(PSF.Regex.GREEDY_WHITESPACE);
							temp = line1[0] + "\t" + line1[1] + "\t" + line1[2];

							for (int k = 3; k < line1.length; k++) {
								if (line1[k].equals(line2[k])) {
									temp += "\t" + line1[k];
									sames++;
								} else if (line1[k].equals("0")) {
									temp += "\t" + line2[k];
								} else if (line2[k].equals("0")) {
									temp += "\t" + line1[k];
								} else {
									off = Integer.valueOf(line1[k]).intValue() - Integer.valueOf(line2[k]).intValue();
									store[(k - 3) / 2]++;
									if (offs.elementAt((k - 3) / 2).containsKey(off + "")) {
										offCount = Integer.valueOf(offs.elementAt((k - 3) / 2).get(off + "")).intValue()
															 + 1;
									} else {
										offCount = 1;
									}
									offs.elementAt((k - 3) / 2).put(off + "", offCount + "");

									diffs++;
									if (differs.containsKey(id)) {
										differs.put(id, (Integer.valueOf(differs.get(id)).intValue() + 1) + "");
									} else {
										differs.put(id, "1");
									}
									temp += "\t0";
								}
							}

							hash.put(id, temp);
						} else {
							hash.put(id, temp);
							individuals.add(id);
						}
					}
				}

				Collections.sort(individuals);

				for (String i : individuals) {
					writer.println(hash.get(i));
				}

				mismatchedMarkers.println("Chromosome " + chr);
				for (int j = 0; j < store.length; j++) {
					if (store[j] < 10) {
						temp = "-";
					} else {
						Enumeration<String> enumer = offs.elementAt(j).keys();
						totalOffCount = maxOffCount = maxOff = 0;
						while (enumer.hasMoreElements()) {
							off = Integer.valueOf(enumer.nextElement()).intValue();
							offCount = Integer.valueOf(offs.elementAt(j).get(off + "")).intValue();
							if (offCount > maxOffCount) {
								maxOffCount = offCount;
								maxOff = off;
							}
							totalOffCount += offCount;
						}
						temp = (int) (100 * (double) maxOffCount / totalOffCount) + "% shifted "
									 + (maxOff > 0 ? "+" + maxOff : "" + maxOff);
					}
					mismatchedMarkers.println(markerNames[j + 3] + "\t" + store[j] + "\t" + temp);
				}
				mismatchedMarkers.println();

				readers[0].close();
				readers[1].close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			writer.close();
		}

		Collections.sort(individuals);

		for (String i : individuals) {
			if (differs.containsKey(i)) {
				mismatchGenotypes.println(i + "\t" + differs.get(i));
			}
		}

		System.err.println("Overall error rate was " + diffs + "/" + (sames + diffs) + " ("
											 + ((double) diffs / (double) (sames + diffs) * 100) + "%)");

		mismatchedMarkers.close();
		mismatchGenotypes.close();
		dupes.close();

	}

	public static void main(String[] args) throws IOException {
		if (args.length != 0) {
			System.err.println("There's no need to argument.");
		}
		try {
			new meld();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
