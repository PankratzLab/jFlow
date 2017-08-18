package org.genvisis.stats;

import java.io.PrintWriter;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

/**
 * Class to compute quantiles of a data set, and facilitates quantile designation
 * http://en.wikipedia.org/wiki/Quantile
 */
public class Quantiles {
	public static double NAN_MEMBERSHIP = 0d;// set to 0 since this can't be a quantile, I think
	// private double[] data;
	private final double[] qs;
	private final double[] quantileMembership;
	private double[] quantiles;
	private String title;
	private final Logger log;

	public Quantiles(double[] data, int numQs, Logger log) {
		this(data, getQ_Quantiles(numQs, log), log);
	}

	public Quantiles(double[] data, double[] qs, Logger log) {
		super();
		this.qs = qs;
		quantiles = new double[qs.length];
		this.log = log;
		quantiles = initQuantiles(data, qs, log);
		quantileMembership = determineQuantileMembership(data, quantiles, qs, log);
		title = "Quantile";

	}

	private static double[] initQuantiles(double[] data, double[] qs, Logger log) {
		double[] quantiles = null;
		double[] tmp = ArrayUtils.removeNaN(data);
		if (tmp.length == 0) {
			log.reportError("Found all NaN values, setting all to quantiles 0");
			quantiles = new double[] {.5};
		} else {
			quantiles = ArrayUtils.quantsExclusive(tmp, qs);
		}
		return quantiles;

	}

	public Logger getLog() {
		return log;
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public double[] getQs() {
		return qs;
	}

	public int[] getQsAsRoundedInt() {
		return roundInt(qs, 10, log);
	}

	public double[] getQuantileMembership() {
		return quantileMembership;
	}

	public int[] getQuantileMembershipAsRoundedInt() {
		int mult = quantiles.length <= 10 ? 10
																			: quantiles.length <= 100 ? 100
																																: quantiles.length <= 1000 ? 1000
																																													 : 10000;
		return roundInt(quantileMembership, mult, log);
	}

	private static int[] roundInt(double[] quantileMembership, int mult, Logger log) {
		int[] qMembers = new int[quantileMembership.length];
		for (int i = 0; i < qMembers.length; i++) {
			if (mult > 10000) {
				log.reportTimeWarning("That is a lot of quintiles to be using this method");
			}
			qMembers[i] = (int) Math.round(10 * quantileMembership[i]);
		}
		return qMembers;
	}

	public double[] getQuantiles() {
		return quantiles;
	}

	private static double[] determineQuantileMembership(double[] data, double[] quantiles,
																											double[] qs, Logger log) {
		double[] quantileMembership = new double[data.length];
		int keysData[] = Sort.getSortedIndices(data);
		int keysQuant[] = Sort.getSortedIndices(quantiles);
		int currentQuantileIndex = 0;
		for (int i = 0; i < data.length; i++) {
			// System.out.println(data[keysData[i]]+"\t"+i);
			if (Double.isNaN(data[keysData[i]])) {
				quantileMembership[keysData[i]] = NAN_MEMBERSHIP;
			} else {
				if (currentQuantileIndex == quantiles.length - 1
						&& data[keysData[i]] > quantiles[keysQuant[currentQuantileIndex]]) {
					quantileMembership[keysData[i]] = 1;
				} else {
					if (data[keysData[i]] > quantiles[keysQuant[currentQuantileIndex]]) {
						while (data[keysData[i]] > quantiles[keysQuant[currentQuantileIndex]]) {
							// System.out.println(data[keysData[i]] + "\t" +
							// quantiles[keysQuant[currentQuantileIndex]] + "\t" + currentQuantileIndex);
							if (currentQuantileIndex == quantiles.length - 1) {
								break;
							} else {
								currentQuantileIndex++;
							}
						}
					}
					quantileMembership[keysData[i]] = qs[keysQuant[currentQuantileIndex]];
				}
			}
		}
		return quantileMembership;
	}

	private static double[] getQ_Quantiles(int numQ, Logger log) {
		if (numQ > 0) {
			double step = (double) 1 / numQ;
			double inc = step;
			double[] qs = new double[numQ - 1];
			for (int i = 0; i < numQ - 1; i++) {
				qs[i] = inc;
				inc += step;
			}
			return qs;
		} else {
			log.reportError("Number of Quantiles must be greater than 0");
			return null;
		}
	}

	public static Quantiles[] qetQuantilesFor(int numQ, double[][] variableDominantMatrix,
																						String[] variableTitles, Logger log) {
		if (variableTitles != null && variableTitles.length != variableDominantMatrix.length) {
			log.reportError("titles must be the same length as the data matrix, titles="
											+ variableTitles.length + " vs " + variableDominantMatrix.length);
			return null;
		}

		Quantiles[] quantiles = new Quantiles[variableDominantMatrix.length];

		for (int i = 0; i < quantiles.length; i++) {
			quantiles[i] = new Quantiles(variableDominantMatrix[i], numQ, log);
			quantiles[i].setTitle(variableTitles[i]);
		}
		return quantiles;
	}

	public static Quantiles[] qetQuantilesFor(int[] numQs, double[] variableDominant,
																						String variableTitle, Logger log) {
		Quantiles[] quantiles = new Quantiles[numQs.length];
		for (int i = 0; i < numQs.length; i++) {
			quantiles[i] = new Quantiles(variableDominant, numQs[i], log);
			if (variableDominant != null) {
				quantiles[i].setTitle(variableTitle);
			}
		}
		return quantiles;
	}

	public static void developQuantiles(String fileName, int[] toQuantileColumns, int numQ,
																			Logger log) {

		String[][] toQ = HashVec.loadFileToStringMatrix(fileName, false, null);// sample,data
		double[][] qData = new double[toQuantileColumns.length][toQ.length];
		for (int i = 0; i < toQ.length; i++) {// for sample
			for (int j = 0; j < toQuantileColumns.length; j++) {// for type
				qData[j][i] = Double.parseDouble(toQ[i][toQuantileColumns[j]]);
			}
		}
		// System.out.println(qData[0].length + "\t" + qData.length);
		int[][] memberships = new int[toQuantileColumns.length][];
		for (int i = 0; i < qData.length; i++) {
			double[] curData = qData[i];
			Quantiles quantiles = new Quantiles(curData, numQ, log);
			memberships[i] = quantiles.getQuantileMembershipAsRoundedInt();
			// System.out.println(Array.toStr(quantiles.getQuantiles()) + "\n" + Array.toStr(qData[i]));

		}
		String output = ext.addToRoot(fileName, ".quant");
		try {
			PrintWriter writer = Files.openAppropriateWriter(output);
			for (int i = 0; i < memberships[0].length; i++) {// for sample,
				writer.print(ArrayUtils.toStr(toQ[i]));
				for (int[] membership : memberships) {
					writer.print("\t" + membership[i]);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "D:/data/LLFS_GWAS/QPCR_MITO/QPCR_quantiles.txt";
		// String logfile = null;

		developQuantiles(filename, new int[] {1, 2}, 10, new Logger());
		String usage = "\n" + "stats.Quantiles requires 0-1 arguments\n" + "   (1) filename (i.e. file="
									 + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				// logfile = args[i].split("=")[1];
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
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
