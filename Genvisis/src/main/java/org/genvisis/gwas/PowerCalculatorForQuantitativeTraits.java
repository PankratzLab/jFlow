package org.genvisis.gwas;

// import java.io.*;
import java.util.Date;
import java.util.Hashtable;
import java.util.Map;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Internat;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.ProbDist;

public class PowerCalculatorForQuantitativeTraits {
	// public static final double[] MAFs = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
	// 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
	// public static final double[] MAFs = {0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25,
	// 0.30, 0.40, 0.50};
	public static final double[] MAFs = {	0.0001, 0.0002, 0.0004, 0.0005, 0.0007, 0.001, 0.002, 0.004,
																				0.005, 0.007, 0.01, 0.02, 0.04, 0.05, 0.07, 0.10, 0.20,
																				0.40, 0.50}; //
	public static final double[] ALPHAS = {0.05, 0.01, 0.0000025, 0.000000227, 0.00000005};
	public static final double[] SIGMAS = {0.10, 0.20, 0.40, 0.60, 0.80, 1.00};
	public static final double[] VARIANCES_EXPLAINED = {0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005,
																											0.01, 0.02, 0.05, 0.10};
	// public static final double[] MAFs = {0.20};
	public static final String[] FORMATTING_TO_REMOVE = {"<em><font color=\"navy\">", "</font></em>"};
	public static final double VARIANCE_INCREMENT = 0.001;

	public static int getSampleSize(int sampleSize, double varianceExplained, double maf,
																	double alpha) throws Exception {
		String[] results, line, cells;
		String trav;
		double[] powers;
		int[] sampleSizes;

		if (alpha < 1E-8) {
			System.err.println("Error - cannot set alpha to less than 1E-8; truncating");
			alpha = 1E-8;
		}

		Map<String, String> data = new Hashtable<String, String>();
		data.put("vq", ext.formDeci(varianceExplained, 4));
		data.put("da", "0");
		data.put("nodom", "TRUE");
		data.put("p", ext.prettyP(maf));
		data.put("m1", ext.prettyP(maf));
		data.put("dprime", "1.0");
		data.put("sibr", "1");
		data.put("n", sampleSize + "");
		data.put("s", "1"); // singletons, not siblings
		data.put("alpha", ext.prettyP(alpha, 2, 100, 2, true));
		data.put("power", "0.80");
		results = Internat.doSubmit("http://pngu.mgh.harvard.edu/~purcell/cgi-bin/qtlassoc.cgi", data,
																1000);

		Files.writeArray(results, "D:/test.html");

		if (results[0].contains("Overall ")) {
			trav = results[0].substring(results[0].indexOf("Overall "));
			for (String element : FORMATTING_TO_REMOVE) {
				trav = ext.replaceAllWith(trav, element, "");
			}
			line = trav.split("<tr>");

			powers = new double[5];
			sampleSizes = new int[5];
			for (int i = 0; i < 5; i++) {
				cells = line[i + 2].trim().split("<td>");
				for (int j = 0; j < 3; j++) {
					powers[i] = Double.parseDouble(cells[2].substring(0, cells[2].indexOf("<")).trim());
					// sampleSizes[i] = Integer.parseInt(cells[3].substring(0, cells[3].indexOf("<")).trim());
					trav = cells[3].substring(0, cells[3].indexOf("<")).trim();
					sampleSizes[i] = trav.equals("inf") ? Integer.MAX_VALUE : (int) Double.parseDouble(trav);
				}
				// System.out.println(powers[i]+"\t"+sampleSizes[i]);
			}
			// System.out.println("\t\t\t\t\t\t\t"+relativeRisk+"\t"+powers[4]+"\t"+sampleSizes[4]);
			return sampleSizes[4];
		}

		return -9;
	}

	public static double getVarianceExplainedAtEightyPercentPower(double maf, int sampleSize,
																																double alpha) throws Exception {
		boolean found;
		int index, prev;
		int[] array;
		double varianceExplained;


		found = false;
		index = 20;
		array = ArrayUtils.intArray(10000, -1);
		while (!found) {
			varianceExplained = 1 + index * VARIANCE_INCREMENT;
			if (array[index] == -1) {
				array[index] = getSampleSize(sampleSize, varianceExplained, maf, alpha);
			} else if (array[index] == -9) {
				return -9;
			} else if (array[index] <= sampleSize && array[index - 1] >= sampleSize) {
				return varianceExplained;
			} else if (array[index] > sampleSize) {
				prev = index;
				do {
					prev++;
				} while (prev - index < 20 && array[prev] == -1);
				index = (int) Math.ceil((prev - index) / 2.0) + index;
			} else {
				prev = index;
				do {
					prev--;
				} while (array[prev] == -1);
				index = (int) Math.floor((index - prev) / 2.0) + prev;
			}
		}

		return -7;
	}

	public static void rangeOfMaf(int sampleSize, int numTests) throws Exception {
		double alpha, varianceExplained;

		alpha = 0.05 / numTests;
		System.out.println("n = " + sampleSize);
		System.out.println("n tests = " + numTests + " (alpha=" + ext.prettyP(alpha) + ")");
		System.out.println();
		System.out.println("MinorAlleleFreq\tVarianceExplained @80% power");
		for (double maf : MAFs) {
			// System.out.println(MAFs[mafIndex]+"\t"+getSampleSize(prevalence, 1.6, MAFs[mafIndex],
			// numCases, numControls, alpha, false));
			varianceExplained = getVarianceExplainedAtEightyPercentPower(maf, sampleSize, alpha);
			System.out.println(maf	+ "\t"
													+ (varianceExplained == -9	? "failed"
																											: ext.formDeci(varianceExplained, 2)));
		}
	}

	public static void rangeOfVarianceExplained(int numTests) throws Exception {
		double alpha;
		int sampleSize;

		alpha = 0.05 / numTests;
		System.out.println("cells = the total sample size required");
		System.out.println("n tests = " + numTests + " (alpha=" + ext.prettyP(alpha) + ")");
		System.out.println();
		System.out.println("\tVariance Explained");
		System.out.print("MinorAlleleFreq");
		for (double element : VARIANCES_EXPLAINED) {
			System.out.print("\t" + ext.formDeci(element, 4, true));
		}
		System.out.println();
		// for (int mafIndex = 0; mafIndex < MAFs.length; mafIndex++) {
		// System.out.print(MAFs[mafIndex]);
		// for (int rrIndex = 0; rrIndex < VARIANCES_EXPLAINED.length; rrIndex++) {
		// System.out.print("\t"+getSampleSize(10000, VARIANCES_EXPLAINED[rrIndex], MAFs[mafIndex],
		// alpha));
		// }
		// System.out.println();
		// }
		for (double element : ALPHAS) {
			System.out.print("alpha=" + element);
			for (double element2 : VARIANCES_EXPLAINED) {
				sampleSize = getSampleSize(10000, element2, 0.10, element);
				System.out.print("\t" + (sampleSize == Integer.MAX_VALUE ? "Inf" : sampleSize));
			}
			System.out.println();
		}
	}

	public static double[] simulate(int n, double sigma, double maf) {
		// PrintWriter writer;
		double[] values, genotypes;

		values = new double[n];
		genotypes = new double[n];

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < 2; j++) {
				if (Math.random() < maf) {
					genotypes[i]++;
				}
			}
			values[i] = (Math.random() < 0.50 ? -1 : 1) * ProbDist.NormDistReverse(Math.random());
			values[i] += genotypes[i] * sigma;
		}

		LeastSquares reg;
		reg = new LeastSquares(values, Matrix.toMatrix(genotypes), true, false);
		// System.out.println(reg.getSummary());
		//
		// String str;
		// try {
		// writer = new PrintWriter(new FileWriter("dump.txt"));
		// writer.println("Value\tgenotype");
		// str = "Value\tgenotype\n";
		// for (int i = 0; i < n; i++) {
		// writer.println(values[i]+"\t"+genotypes[i]);
		// str += values[i]+"\t"+genotypes[i]+"\n";
		// }
		// writer.close();
		// ext.setClipboard(str);
		// } catch (Exception e) {
		// System.err.println("Error writing to " + "dump.xln");
		// e.printStackTrace();
		// }
		if (reg.analysisFailed()) {
			return new double[] {0, 1};
		}

		return new double[] {reg.getRsquare(), reg.getSigs()[1]};
	}

	public static double[] determineMeanRsqAndPower(int n, double sigma, double maf, int numReps,
																									double alpha) {
		double[] rsqs, power, values;

		rsqs = new double[numReps];
		power = new double[numReps];
		for (int i = 0; i < rsqs.length; i++) {
			// if ((i+1) % 10 == 0) {
			// System.out.print(".");
			// }
			values = simulate(n, sigma, maf);
			rsqs[i] = values[0];
			power[i] = values[1] < alpha ? 1 : 0;
		}
		// System.out.println();

		return new double[] {ArrayUtils.mean(rsqs), ArrayUtils.mean(power)};
	}

	public static void rangeOfSigmaShiftsAndMAFsViaSimulation(int n, int numReps,
																														int numTests) throws Exception {
		double[] meanRsqAndPower;
		double alpha;

		alpha = 0.05 / numTests;
		System.out.println("sample size (n) = " + ext.addCommas(n));
		System.out.println("num reps = " + numReps);
		System.out.println("n tests = " + numTests + " (alpha=" + ext.prettyP(alpha) + ")");

		System.out.println();
		System.out.println("\tSigma Shift Per Variant");
		System.out.print("MinorAlleleFreq");
		for (double element : SIGMAS) {
			System.out.print("\t" + ext.formDeci(element, 4, true));
		}
		System.out.println();
		for (double maf : MAFs) {
			System.out.print(maf);
			for (double element : SIGMAS) {
				meanRsqAndPower = determineMeanRsqAndPower(n, element, maf, numReps, alpha);
				System.out.print("\t"	+ ext.formDeci(meanRsqAndPower[0], 4, true) + ";"
													+ (int) (meanRsqAndPower[1] * 100) + "%");
			}
			System.out.println();
		}
	}

	public static void simulatedBetting(int n, int numReps, int numTests,
																			double proportionUpweighted) throws Exception {
		double[] meanRsqAndPower;
		double alpha;

		alpha = 0.05 / numTests;
		System.out.println("sample size (n) = " + ext.addCommas(n));
		System.out.println("num reps = " + numReps);
		System.out.println("n tests = " + numTests + " (alpha=" + ext.prettyP(alpha) + ")");

		System.out.println();
		System.out.println("\tSigma Shift Per Variant");
		System.out.print("MinorAlleleFreq");
		for (double element : SIGMAS) {
			System.out.print("\t" + ext.formDeci(element, 4, true));
		}
		System.out.println();
		for (double maf : MAFs) {
			System.out.print(maf);
			for (double element : SIGMAS) {
				meanRsqAndPower = determineMeanRsqAndPower(n, element, maf, numReps, alpha);
				System.out.print("\t"	+ ext.formDeci(meanRsqAndPower[0], 4, true) + ";"
													+ (int) (meanRsqAndPower[1] * 100) + "%");
			}
			System.out.println();
		}
	}


	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "PowerCalculator.dat";

		String usage = "\n"	+ "gwas.PowerCalculator requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
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
			// rangeOfMaf(0.15, 0.01, 200, 200, 6, false);
			// rangeOfVarianceExplained(20000);
			// getSampleSize();
			// System.out.println(Array.toStr(determineMeanRsqAndPower(10000,0.5,0.30, 100, 0.001)));
			long time;
			time = new Date().getTime();
			// rangeOfValues(10000, 100, 220000);
			// rangeOfSigmaShiftsAndMAFsViaSimulation(10000, 1000, 1000000); // inflammation discovery
			// rangeOfSigmaShiftsAndMAFsViaSimulation(5000, 1000, 100); // inflammation replication
			// rangeOfSigmaShiftsAndMAFsViaSimulation(200, 1000, 75); // hearing inflammation sequence
			// discovery
			// rangeOfSigmaShiftsAndMAFsViaSimulation(2500, 1000, 40); // hearing inflammation sequence
			// discovery
			// rangeOfSigmaShiftsAndMAFsViaSimulation(50000, 1000, 1000000); // Genvisis analyses
			rangeOfSigmaShiftsAndMAFsViaSimulation(5500, 1000, 7767); // Power for IISS analyses
			rangeOfSigmaShiftsAndMAFsViaSimulation(5500, 1000, 1000000); // Power for IISS analyses
			rangeOfSigmaShiftsAndMAFsViaSimulation(50, 100, 5); // Power for IISS analyses


			// simulatedBetting(10000, 1000, 1, 0.00);
			System.out.println("Finished in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
