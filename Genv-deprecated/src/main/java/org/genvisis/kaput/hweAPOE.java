package org.genvisis.kaput;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

public class hweAPOE {
	public static String[] GENOTYPES = {"44", "34", "24", "33", "23", "22"};

	public hweAPOE(String filename, int numReps) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, trav, gender, alleles, affStat, mmse, aoo, demented, duration, smoked, data,
				ejukashun;
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		Vector<String> fams = new Vector<String>();
		Vector<String> v;
		int group;
		int[] totals = new int[GENOTYPES.length * 2], counts;
		double[] totalMeans = new double[GENOTYPES.length * 2], means;
		double affTotal, unaffTotal;

		try {
			reader = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - File '" + filename + "' not found in runtime directory");
			System.exit(1);
		}

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			trav = st.nextToken();
			temp = st.nextToken();
			st.nextToken();
			st.nextToken();

			gender = st.nextToken();
			if (gender.equals("2")) {
				gender = "0";
			}
			affStat = st.nextToken();
			alleles = st.nextToken();
			st.nextToken();
			aoo = st.nextToken();
			st.nextToken();
			mmse = st.nextToken();
			ejukashun = st.nextToken();
			duration = st.nextToken();
			if (duration.equals(".") && !aoo.equals(".") && !alleles.equals("0")) {
				System.err.println("Error - "	+ trav + "-" + temp
														+ " has an age of onset but not a duration; why?");
				System.err.println("  removing " + trav + "-" + temp + " from analysis");
				aoo = ".";
			}
			smoked = st.nextToken();

			group = -1;
			for (int i = 0; i < GENOTYPES.length; i++) {
				if (GENOTYPES[i].equals(alleles)) {
					group = i;
				}
			}

			if (!mmse.equals(".") && (group != -1) && affStat.equals("2") && !aoo.equals(".")) {
				if (smoked.equals(".")) {
					System.err.println("Warning: "	+ trav + "-" + temp
															+ " does not know if he/she ever smoked; assuming not");
					smoked = "0";
				}

				demented = "1";
				if ((ejukashun.equals(".") || Integer.valueOf(ejukashun).intValue() <= 8)
						&& Integer.valueOf(mmse).intValue() < 21) {
					demented = "0";
				} else if (Integer.valueOf(ejukashun).intValue() <= 12
										&& Integer.valueOf(mmse).intValue() < 23) {
					demented = "0";
				} else if (Integer.valueOf(mmse).intValue() < 24) {
					demented = "0";
				}
				data = group	+ "\t" + demented + "\t" + aoo + "\t" + gender + "\t" + duration + "\t"
								+ smoked + "\t" + ejukashun;

				if (addPerson(hash, trav, data)) {
					fams.add(trav);
				}

			}
		}
		reader.close();

		writer = new PrintWriter(new FileWriter("hweAPOE_MCMC.out"));
		writer.println("e4e4\t\t\t\te4e3\t\t\t\te4e2\t\t\t\te3e3\t\t\t\te3e2\t\t\t\te2e2\t\t\t");
		writer.println("aff\t\tunaff\t\taff\t\tunaff\t\taff\t\tunaff\t\taff\t\tunaff\t\taff\t\tunaff\t\taff\t\tunaff\t\t");

		for (int rep = 0; rep < numReps; rep++) {
			if (rep % 1000 == 0) {
				System.out.println("Rep " + rep);
			}

			counts = new int[GENOTYPES.length * 2];
			means = new double[GENOTYPES.length * 2];
			for (int j = 0; j < fams.size(); j++) {
				v = hash.get(fams.elementAt(j));
				temp = v.elementAt((int) (Math.random() * v.size()));
				st = new StringTokenizer(temp);

				group = Integer.valueOf(st.nextToken()).intValue();
				demented = st.nextToken();
				aoo = st.nextToken();

				counts[2 * group + Integer.valueOf(demented).intValue()]++;
				means[2 * group + Integer.valueOf(demented).intValue()] +=
																																Double.valueOf(aoo).doubleValue();
			}
			for (int j = 0; j < counts.length; j++) {
				means[j] /= counts[j];
				writer.print(counts[j] + "\t" + ext.formDeci(means[j], 2, true) + "\t");
				totals[j] += counts[j];
				totalMeans[j] += means[j];
			}
			writer.println();

		}
		writer.close();

		writer = new PrintWriter(new FileWriter("hweAPOE_MCMC-summary.out"));

		writer.println("Filename: " + filename);
		writer.println("Number of replicates: " + numReps);
		writer.println();
		writer.println();

		affTotal = unaffTotal = 0;
		means = new double[GENOTYPES.length * 2];
		for (int j = 0; j < GENOTYPES.length * 2; j++) {
			if (j % 2 == 0) {
				affTotal += (double) totals[j] / (double) numReps;
			} else {
				unaffTotal += (double) totals[j] / (double) numReps;
			}
		}
		writer.println("\tDemented\tNot Demented");
		for (int j = 0; j < GENOTYPES.length; j++) {
			means[2 * j] = (double) totals[2 * j] / (double) numReps;
			writer.print(GENOTYPES[j]	+ "\t" + ext.formDeci(means[2 * j], 1, true) + " ("
										+ ext.formDeci(100 * means[2 * j] / affTotal, 1, true) + "%)\t");
			means[2 * j + 1] = (double) totals[2 * j + 1] / (double) numReps;
			writer.println(ext.formDeci(means[2 * j + 1], 1, true)	+ " ("
											+ ext.formDeci(100 * means[2 * j + 1] / unaffTotal, 1, true) + "%)\t");
		}

		temp = apoeOR(means);
		System.out.println(temp);

		writer.println();
		writer.println(temp);
		writer.println();
		writer.close();
	}

	public String apoeOR(double[] multipleMeans) throws IOException {
		String line = "";
		double N21, N22, N11, N12;

		N11 = multipleMeans[6];
		N12 = multipleMeans[7];

		for (int i = 0; i < multipleMeans.length / 2; i++) {
			N21 = multipleMeans[i * 2 + 0];
			N22 = multipleMeans[i * 2 + 1];

			line += oddsRatio(N21, N22, N11, N12) + "\n";
		}

		return line;
	}

	public boolean addPerson(	Hashtable<String, Vector<String>> hash, String trav,
														String data) throws IOException {
		boolean newPerson = false;
		Vector<String> v;

		if (!hash.containsKey(trav)) {
			newPerson = true;
			hash.put(trav, v = new Vector<String>());
		} else {
			v = hash.get(trav);
		}
		v.add(data);

		return newPerson;
	}

	public String oddsRatio(double N21, double N22, double N11, double N12) throws IOException {
		String line;
		double OR, logOR, se_logOR;

		OR = (N21 / N22) / (N11 / N12);
		logOR = Math.log(OR);
		se_logOR = Math.sqrt(1 / N21 + 1 / N22 + 1 / N11 + 1 / N12);

		line = ext.formDeci(OR, 2, true)	+ " (95% CI: ("
						+ ext.formDeci(Math.exp(logOR - 1.96 * se_logOR), 2, true) + ", "
						+ ext.formDeci(Math.exp(logOR + 1.96 * se_logOR), 2, true) + "); p="
						+ ext.formDeci(ProbDist.NormDist(Math.abs(logOR / se_logOR)), 3, true) + ")";

		N21 += 0.5;
		N22 += 0.5;
		N11 += 0.5;
		N12 += 0.5;

		OR = (N21 / N22) / (N11 / N12);
		logOR = Math.log(OR);
		se_logOR = Math.sqrt(1 / N21 + 1 / N22 + 1 / N11 + 1 / N12);

		line += " corrected: "	+ ext.formDeci(OR, 2, true) + " (95% CI: ("
						+ ext.formDeci(Math.exp(logOR - 1.96 * se_logOR), 2, true) + ", "
						+ ext.formDeci(Math.exp(logOR + 1.96 * se_logOR), 2, true) + "); p="
						+ ext.formDeci(ProbDist.NormDist(Math.abs(logOR / se_logOR)), 3, true) + ")";

		return line;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "struct111+white(APOE-MMSE-education,duration,tobacky).pre";
		int numReps = 50000;

		String usage = "\n"	+ "park.hweAPOE requires 2 arguments:\n"
										+ "   (1) a pre-like file - FAMID INDID AFF/UNAFF group\n"
										+ "       (for an APOE analysis this would be 2, 3 or 4)\n"
										+ "       (i.e. file=" + filename + " (default))\n"
										+ "   (2) number of replicates to perform\n" + "       (i.e. reps=" + numReps
										+ " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("reps=")) {
				numReps = Integer.valueOf(arg.split("=")[1]).intValue();
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length == 0) {
			System.err.println("Warning: using defaults (file=" + filename + " reps=" + numReps + ")");
		}
		try {
			new hweAPOE(filename, numReps);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
