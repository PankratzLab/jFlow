package org.genvisis.kaput;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

public class bootstrapSNP {
	public static int NUM_BOOTSTRAP_REPS = 5000;

	public bootstrapSNP(String filename, boolean single) throws IOException {
		Hashtable<String, String> hash = new Hashtable<String, String>();
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		StringTokenizer st;
		double[] bootstrapped;

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			hash.put(hash.size() + "", st.nextToken());
		}
		reader.close();

		bootstrapped = Bootstrap(hash);

		System.out.println("bootstrapped mean (95% CI):\t"	+ ext.formDeci(bootstrapped[0], 4, true)
												+ " (" + ext.formDeci(bootstrapped[1], 4, true) + ","
												+ ext.formDeci(bootstrapped[2], 4, true) + ")");
		System.out.println("corresponding p-values:\t"
													+ ext.formDeci(ProbDist.ChiDist(bootstrapped[0], 1), 4, true) + " ("
												+ ext.formDeci(ProbDist.ChiDist(bootstrapped[1], 1), 4, true) + ","
												+ ext.formDeci(ProbDist.ChiDist(bootstrapped[2], 1), 4, true) + ")");
		System.out.println();

	}

	public bootstrapSNP(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st = null;
		String temp;
		double[] sums = new double[9];
		int count = 0;
		double[] bootstrapped;

		Hashtable<String, String> hash22 = new Hashtable<String, String>();
		Hashtable<String, String> hash12 = new Hashtable<String, String>();
		Hashtable<String, String> hash11 = new Hashtable<String, String>();

		reader = new BufferedReader(new FileReader(filename));
		writer = new PrintWriter(new FileWriter(filename.substring(0, filename.length() - 4)
																						+ "-summary.out"));

		reader.readLine();
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			count++;
			for (int i = 0; i < 6; i++) {
				sums[i] += Double.valueOf(st.nextToken()).doubleValue();
			}

			for (int i = 0; i < 4; i++) {
				st.nextToken();
			}

			temp = st.nextToken();
			sums[6] += Double.valueOf(temp).doubleValue();
			hash22.put(hash22.size() + "", temp);

			temp = st.nextToken();
			sums[7] += Double.valueOf(temp).doubleValue();
			hash12.put(hash12.size() + "", temp);

			temp = st.nextToken();
			sums[8] += Double.valueOf(temp).doubleValue();
			hash11.put(hash11.size() + "", temp);
		}
		reader.close();

		writer.println("Number of original replicates: " + count);
		writer.println("Number of bootstrapping iterations: " + NUM_BOOTSTRAP_REPS);
		writer.println();

		writer.println("average # cases with 22:\t" + ext.formDeci(sums[2] / count, 4, true));
		writer.println("average # cases with 12:\t" + ext.formDeci(sums[1] / count, 4, true));
		writer.println("average # cases with 11:\t" + ext.formDeci(sums[0] / count, 4, true));
		writer.println();
		writer.println("average # controls with 22:\t" + ext.formDeci(sums[5] / count, 4, true));
		writer.println("average # controls with 12:\t" + ext.formDeci(sums[4] / count, 4, true));
		writer.println("average # controls with 11:\t" + ext.formDeci(sums[3] / count, 4, true));
		writer.println();

		bootstrapped = Bootstrap(hash22);
		writer.println("average chi square for the 22 genotype:\t"
										+ ext.formDeci(sums[6] / count, 4, true));
		writer.println("\t\tbootstrapped mean (95% CI):\t"	+ ext.formDeci(bootstrapped[0], 4, true)
										+ " (" + ext.formDeci(bootstrapped[1], 4, true) + ","
										+ ext.formDeci(bootstrapped[2], 4, true) + ")");
		writer.println("\t\tcorresponding p-values:\t"
											+ ext.formDeci(ProbDist.ChiDist(bootstrapped[0], 1), 4, true) + " ("
										+ ext.formDeci(ProbDist.ChiDist(bootstrapped[1], 1), 4, true) + ","
										+ ext.formDeci(ProbDist.ChiDist(bootstrapped[2], 1), 4, true) + ")");
		writer.println();

		bootstrapped = Bootstrap(hash12);
		writer.println("average chi square for the 12 genotype:\t"
										+ ext.formDeci(sums[7] / count, 4, true));
		writer.println("\t\tbootstrapped mean (95% CI):\t"	+ ext.formDeci(bootstrapped[0], 4, true)
										+ " (" + ext.formDeci(bootstrapped[1], 4, true) + ","
										+ ext.formDeci(bootstrapped[2], 4, true) + ")");
		writer.println("\t\tcorresponding p-values:\t"
											+ ext.formDeci(ProbDist.ChiDist(bootstrapped[0], 1), 4, true) + " ("
										+ ext.formDeci(ProbDist.ChiDist(bootstrapped[1], 1), 4, true) + ","
										+ ext.formDeci(ProbDist.ChiDist(bootstrapped[2], 1), 4, true) + ")");
		writer.println();

		bootstrapped = Bootstrap(hash11);
		writer.println("average chi square for the 11 genotype:\t"
										+ ext.formDeci(sums[8] / count, 4, true));
		writer.println("\t\tbootstrapped mean (95% CI):\t"	+ ext.formDeci(bootstrapped[0], 4, true)
										+ " (" + ext.formDeci(bootstrapped[1], 4, true) + ","
										+ ext.formDeci(bootstrapped[2], 4, true) + ")");
		writer.println("\t\tcorresponding p-values:\t"
											+ ext.formDeci(ProbDist.ChiDist(bootstrapped[0], 1), 4, true) + " ("
										+ ext.formDeci(ProbDist.ChiDist(bootstrapped[1], 1), 4, true) + ","
										+ ext.formDeci(ProbDist.ChiDist(bootstrapped[2], 1), 4, true) + ")");
		writer.println();

		writer.close();
	}

	public double[] Bootstrap(Hashtable<String, String> hash) { // determines the median of the
																															// ordered distribution, 2.5%
																															// percentile, 97.5%
		double[] results = new double[3];
		double sumOfMeans;
		double[] replicates = new double[NUM_BOOTSTRAP_REPS];
		int[] keys;

		for (int i = 0; i < NUM_BOOTSTRAP_REPS; i++) {
			sumOfMeans = 0;
			for (int j = 0; j < hash.size(); j++) {
				sumOfMeans += Double.valueOf(hash.get("" + (int) (Math.random() * hash.size())))
														.doubleValue();
			}
			replicates[i] = sumOfMeans / hash.size();
		}
		keys = Sort.quicksort(replicates);

		results[0] = replicates[keys[(int) (NUM_BOOTSTRAP_REPS * 0.5)]];
		results[1] = replicates[keys[(int) (NUM_BOOTSTRAP_REPS * 0.025)]];
		results[2] = replicates[keys[(int) (NUM_BOOTSTRAP_REPS * 0.975)]];

		return results;
	}

	public static void main(String[] args) throws IOException {
		String usage = "Expecting 1-2 arguments: name of file created by randomSNP (default: \"randomSample.xls\")\n"
										+ "                       : or file with one column (add \"-single\" after filename).\n";
		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			}
		}
		if (args.length > 2) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (args.length == 0) {
				new bootstrapSNP("randomSample.xls");
			} else if (args.length == 1) {
				new bootstrapSNP(args[0]);
			}
			if (args.length == 2 && args[1].equals("-single")) {
				new bootstrapSNP(args[0], true);
			} else {
				System.err.println(usage);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
