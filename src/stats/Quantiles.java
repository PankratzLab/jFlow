package stats;

import common.Array;
import common.Logger;
import common.Sort;

/**
 * Class to compute quantiles of a data set, and facilitates quantile designation http://en.wikipedia.org/wiki/Quantile
 */
public class Quantiles {
	public static double NAN_MEMBERSHIP = 0d;// set to 0 since this can't be a quantile, I think
	// private double[] data;
	private double[] qs;
	private double[] quantileMembership;
	private double[] quantiles;
	private Logger log;

	public Quantiles(double[] data, int numQs, Logger log) {
		this(data, getQ_Quantiles(numQs, log), log);
	}

	public Quantiles(double[] data, double[] qs, Logger log) {
		super();
		// this.data = data;
		this.qs = qs;
		this.quantiles = new double[qs.length];
		this.log = log;
		this.quantiles = Array.quants(Array.removeNaN(data), qs);
		this.quantileMembership = determineQuantileMembership(data, quantiles, qs, log);
	}

	public Logger getLog() {
		return log;
	}

	public double[] getQs() {
		return qs;
	}

	public double[] getQuantileMembership() {
		return quantileMembership;
	}

	public int[] getQuantileMembershipAsRoundedInt() {
		int[] qMembers = new int[quantileMembership.length];
		for (int i = 0; i < qMembers.length; i++) {
			int mult = quantiles.length <= 10 ? 10 : quantiles.length <= 100 ? 100 : quantiles.length <= 1000 ? 1000 : 10000;
			if (mult > 10000) {
				log.reportTimeWarning("That is a lot of quintiles to be using this method");
			}
			qMembers[i] = (int) Math.round((float) 10 * quantileMembership[i]);
		}
		return qMembers;

	}

	public double[] getQuantiles() {
		return quantiles;
	}

	private static double[] determineQuantileMembership(double[] data, double[] quantiles, double[] qs, Logger log) {
		double[] quantileMembership = new double[data.length];
		int keysData[] = Sort.quicksort(data);
		int keysQuant[] = Sort.quicksort(quantiles);
		int currentQuantileIndex = 0;
		for (int i = 0; i < data.length; i++) {
			if (Double.isNaN(data[keysData[i]])) {
				quantileMembership[keysData[i]] = NAN_MEMBERSHIP;
			} else {
				if (currentQuantileIndex == quantiles.length - 1 && data[keysData[i]] > quantiles[keysQuant[currentQuantileIndex]]) {
					quantileMembership[keysData[i]] = 1;
				} else {
					if (data[keysData[i]] > quantiles[keysQuant[currentQuantileIndex]]) {
						while (data[keysData[i]] > quantiles[keysQuant[currentQuantileIndex]]) {
							currentQuantileIndex++;
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
			log.reportTimeError("Number of Quantiles must be greater than 0");
			return null;
		}
	}

}
