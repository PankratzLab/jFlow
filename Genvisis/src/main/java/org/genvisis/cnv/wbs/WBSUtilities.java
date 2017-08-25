/**
 * 
 */
package org.genvisis.cnv.wbs;

import org.apache.commons.math3.random.MersenneTwister;

/**
 * Helper functions for WBS
 */
class WBSUtilities {

	private WBSUtilities() {

	}

	/**
	 * @param n length of data
	 * @param m number of intervals
	 * @param seed to seed a {@link MersenneTwister} and obtain pseudo-random intervals
	 * @return
	 */
	static int[][] randomIntervals(int n, int m, int seed) {
		MersenneTwister mersenneTwister = new MersenneTwister(seed);
		int[][] intervals = new int[2][m];
		for (int i = 0; i < intervals[0].length; i++) {
			double first = Math.ceil(mersenneTwister.nextDouble() * (double) (n - 1));
			double second = first + Math.ceil(mersenneTwister.nextDouble() * (n - first));
			intervals[0][i] = (int) first;
			intervals[1][i] = (int) second;
		}
		return intervals;
	}



	// x[(1:n-lag)] - x[(lag:n)].
	// https://www.math.ucla.edu/~anderson/rw1001/library/base/html/diff.html
	static double[] computeLagOneDifference(double[] x) {
		double[] lag = new double[x.length - 1];
		for (int i = 0; i < x.length - 1; i++) {
			lag[i] = x[i + 1] - x[i];
		}
		return lag;
	}

	static double[] rootTwo(double[] x) {
		double[] rootTwo = new double[x.length];
		double div = Math.sqrt(2);
		for (int i = 0; i < rootTwo.length; i++) {
			rootTwo[i] = x[i] / div;
		}
		return rootTwo;

	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		randomIntervals(100, 5000, 42);

	}
}
