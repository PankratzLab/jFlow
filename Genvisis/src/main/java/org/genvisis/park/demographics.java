// split hash of hashs into two different types of hashes

package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ext;

public class demographics {
	public int[] bb_bs_ss;

	public demographics(String filename) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		Hashtable<String, Hashtable<String, String[]>> hash =
																												new Hashtable<String, Hashtable<String, String[]>>();
		Hashtable<String, String> unaffs = new Hashtable<String, String>();
		Hashtable<String, String[]> h;
		Vector<String> fams = new Vector<String>();
		String[] line, info;
		String temp;
		int count = 0, total, max;
		int[] affPairs = new int[20];
		int[] affInds = new int[20];
		String[] listFamsPairs = new String[20];
		String[] listFamsAff = new String[20];
		Enumeration<String> enumer;
		int males = 0, females = 0, unaffected = 0, brothers, sisters, malesInFam, femalesInFam,
				unaffInFam;
		Vector<String> aoos = new Vector<String>();
		double meanAOO, stdev, d;

		bb_bs_ss = new int[3];

		for (int i = 0; i < 20; i++) {
			listFamsPairs[i] = "";
			listFamsAff[i] = "";
		}

		reader = new BufferedReader(new FileReader(filename));
		while (reader.ready()) {
			temp = reader.readLine();
			line = temp.split("[\\s]+");
			if (line.length < 8) {
				System.err.println("Error - requires at least 8 columns for every row: FamID IndID Father Mother Gender Affection DNA(not_used) AgeOfOnset");
				System.err.println("  got - " + temp);
			}
			if (hash.containsKey(line[0])) {
				h = hash.get(line[0]);
			} else {
				h = new Hashtable<String, String[]>();
				hash.put(line[0], h);
				fams.add(line[0]);
			}
			if (line[7].equals("-99")) {
				line[7] = ".";
			} else if (!line[7].equals(".")) {
				try {
					Double.valueOf(line[7]).doubleValue();
					aoos.add(line[7]);
				} catch (NumberFormatException nfe) {
					System.err.println("Error - Could not parse age of onset for "	+ line[0] + "-" + line[1]
															+ ": " + line[7]);
					System.err.println("        Looking for a number or a missing value character of '.' or '-99'");
					System.exit(3);
				}
			}
			if (line[5].equals("1")) {
				if (unaffs.containsKey(line[0])) {
					unaffs.put(line[0], (Integer.valueOf(unaffs.get(line[0])).intValue() + 1) + "");
				} else {
					unaffs.put(line[0], "1");
				}
			} else if (line[5].equals("2")) {
				if (h.containsKey(line[2] + "/" + line[3])) {
					info = h.get(line[2] + "/" + line[3]);
					info[Integer.valueOf(line[4]).intValue()
								- 1] = (Integer.valueOf(info[Integer.valueOf(line[4]).intValue() - 1]).intValue()
												+ 1) + "";
					total = Integer.valueOf(info[0]).intValue() + Integer.valueOf(info[1]).intValue();
				} else {
					info = new String[2];
					info[0] = info[1] = "0";
					info[Integer.valueOf(line[4]).intValue() - 1] = "1";
					h.put(line[2] + "/" + line[3], info);
				}
			}
		}
		for (int i = 0; i < fams.size(); i++) {
			h = hash.get(fams.elementAt(i));
			enumer = h.keys();
			brothers = sisters = total = count = max = malesInFam = femalesInFam = unaffInFam = 0;
			if (unaffs.containsKey(fams.elementAt(i))) {
				unaffInFam = Integer.valueOf(unaffs.get(fams.elementAt(i))).intValue();
			}
			while (enumer.hasMoreElements()) {
				info = h.get(enumer.nextElement());
				malesInFam += Integer.valueOf(info[0]).intValue();
				femalesInFam += Integer.valueOf(info[1]).intValue();
				count = Integer.valueOf(info[0]).intValue() + Integer.valueOf(info[1]).intValue();
				total += count;
				if (count > max) {
					max = count;
					brothers = Integer.valueOf(info[0]).intValue();
					sisters = Integer.valueOf(info[1]).intValue();
				}
			}
			if (max >= 20 || count >= 20) {
				System.err.println("Error - wasn't expecting a family with more than 20 affected, recompile with a higher limit");
			}

			if (total > 1) {
				males += malesInFam;
				females += femalesInFam;
				unaffected += unaffInFam;
				affPairs[max]++;
				affInds[total]++;
				listFamsPairs[max] += fams.elementAt(i) + "\n";
				listFamsAff[total] += fams.elementAt(i) + "\n";
				tally(brothers, sisters);
			}
		}

		writer = new PrintWriter(new FileWriter(filename + "-demographics.out"));
		writer.println((males + females) + " affected individuals");
		writer.println(unaffected + " unaffected individuals");
		writer.println((males + females + unaffected) + " genotipyped individuals");
		count = 0;
		for (int i = 0; i < 20; i++) {
			count += affPairs[i];
		}
		writer.println(count + " families");

		Collections.sort(aoos);
		meanAOO = 0;
		for (int i = 0; i < aoos.size(); i++) {
			meanAOO += Double.valueOf(aoos.elementAt(i)).doubleValue();
		}
		meanAOO /= aoos.size();
		stdev = 0;
		for (int i = 0; i < aoos.size(); i++) {
			d = Double.valueOf(aoos.elementAt(i)).doubleValue();
			stdev += (d - meanAOO) * (d - meanAOO);
		}
		stdev = Math.sqrt(stdev / (aoos.size() - 1));
		if (aoos.size() == 0) {
			writer.println("Age of onset information was not present in this stuct file");
		} else {
			writer.println("Age of onset information estimated from " + aoos.size() + " individuals:");
			writer.println(ext.formDeci(meanAOO, 1, true)	+ " \u00B1 " + ext.formDeci(stdev, 1, true)
											+ " (" + aoos.elementAt(0) + "-"
											+ aoos.elementAt(aoos.size() - 1) + ")");
		}

		total = bb_bs_ss[0] + bb_bs_ss[1] + bb_bs_ss[2];

		writer.println();
		writer.println(males	+ " males (" + ext.formDeci((double) males / (males + females), 2, true)
										+ ")");
		writer.println(females	+ " females ("
										+ ext.formDeci((double) females / (males + females), 2, true) + ")");
		writer.println();
		writer.println(bb_bs_ss[0]	+ " brother-brother pairs ("
										+ ext.formDeci((double) bb_bs_ss[0] / (double) total, 2, true) + ")");
		writer.println(bb_bs_ss[1]	+ " brother-sister pairs ("
										+ ext.formDeci((double) bb_bs_ss[1] / (double) total, 2, true) + ")");
		writer.println(bb_bs_ss[2]	+ " sister-sister pairs ("
										+ ext.formDeci((double) bb_bs_ss[2] / (double) total, 2, true) + ")");
		writer.println(total + " total pairs");

		writer.println();
		for (int i = 0; i < 20; i++) {
			if (affPairs[i] > 0) {
				writer.println("families with "
													+ (i == 1 ? " singlets/halfsibs = " : i + " affected sibpairs = ")
												+ affPairs[i]);
				if (i != 2) {
					writer.println(listFamsPairs[i]);
				}
			}
		}

		writer.println();
		for (int i = 0; i < 20; i++) {
			if (affInds[i] > 0) {
				writer.println("families with " + i + " affected individuals = " + affInds[i]);
				if (i != 2) {
					writer.println(listFamsAff[i]);
				}
			}
		}

		reader.close();
		writer.close();

	}

	public void tally(int Bs, int Ss) throws IOException {
		int[] sibs = new int[Bs + Ss];

		for (int i = Bs; i < Bs + Ss; i++) {
			sibs[i] = 1;
		}
		for (int i = 0; i < sibs.length - 1; i++) {
			for (int j = i + 1; j < sibs.length; j++) {
				bb_bs_ss[sibs[i] + sibs[j]]++;
			}
		}
	}

	public void tallyOld(int Bs, int Ss) throws IOException {
		if (Bs == 0 && Ss == 0) {
			return;
		} else if (Bs == 0) {
			Ss--;
			bb_bs_ss[2] += Ss;
			tally(Bs, Ss);
		} else {
			Bs--;
			bb_bs_ss[0] += Bs;
			bb_bs_ss[1] += Ss;
			tally(Bs, Ss);
		}
		return;
	}

	public static int combinations(int n, int r) {
		int p, c, rf;

		p = 1;
		for (int w = 1; w <= r - 1; w++) {
			p = (n - w) * p;
		}

		p = n * p;

		rf = 1;
		for (int w = 1; w <= r; w++) {
			rf = w * rf;
		}
		c = p / rf;

		return c;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "struct.dat";

		String usage = "\n"	+ "park.demographics requires 0-1 arguments\n"
										+ "   (1) a struct file (default: file=" + filename + ")\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				if (!new File(filename).exists()) {
					System.err.println("Error - file '" + filename + "' does not exist");
					System.exit(2);
				}
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new demographics(filename);
			// new demographics("struct111+.dat");
			// new demographics("progeni-broad.fam");
			// new demographics("gspd-broad.fam");
			// new demographics("merged-Broad.fam");
			// new demographics("progeni-narrow.fam");
			// new demographics("gspd-narrow.fam");
			// new demographics("merged-Narrow.fam");
			System.out.println(filename + "-demographics.out was created");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
