package org.genvisis.stats;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.common.ext;

/**
 * Utility class for computing
 * <a href="https://en.wikipedia.org/wiki/Fisher%27s_exact_test">Fisher's exact test</a> on a 2x2
 * matrix.
 */
public class FishersExact2by2Calculator {
	private FishersExact2by2Calculator() {
		// Prevent instantiation of static utility class
	}

	private static final List<Double> CACHED_LOG_FACTORIAL = new ArrayList<>();

	/**
	 * @param matrix A 2x2 matrix of counts
	 * @return The two-tailed p-value for the input 2x2 table
	 * @see #getPvalue(int[][], boolean)
	 */
	public static double getPvalue(int[][] matrix) {
		return getPvalue(matrix, false);
	}

	/**
	 * @return The two-tailed p-value for the four counts of a 2x2 table
	 * @see #getPvalue(int, int, int, int, boolean)
	 */
	public static double getPvalue(int a, int b, int c, int d) {
		return getPvalue(a, b, c, d, false);
	}

	/**
	 * @param matrix A 2x2 matrix of counts
	 * @param oneTailed If true, return value is one-tailed p-value. Otherwise, two-tailed.
	 * @return The one or two-tailed p-value for the input 2x2 table
	 * @see #getPvalue(int, int, int, int, boolean)
	 */
	public static double getPvalue(int[][] matrix, boolean oneTailed) {
		if (matrix.length != 2 || matrix[0].length != 2 || matrix[1].length != 2) {
			throw new IllegalArgumentException(FishersExact2by2Calculator.class.getName()
																				 + " requires a 2x2 input matrix to compute a p-value");
		}

		return getPvalue(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1], oneTailed);
	}

	/**
	 * @param oneTailed If true, return value is one-tailed p-value. Otherwise, two-tailed.
	 * @return The one or two-tailed p-value for the four counts of a 2x2 table
	 */
	public static double getPvalue(int a, int b, int c, int d, boolean oneTailed) {
		ensureCacheSize(a, b, c, d);

		int tmp;

		// The algorithm used requires the following assumptions to be true:
		// - the top left to bottom right diagonal contains the smaller values
		// - within both diagonals, the top row contains the smaller of the two values
		if (a * d > b * c) {
			tmp = a;
			a = b;
			b = tmp;

			tmp = c;
			c = d;
			d = tmp;
		}
		if (a > d) {
			tmp = a;
			a = d;
			d = tmp;
		}
		if (b > c) {
			tmp = b;
			b = c;
			c = tmp;
		}

		int original_a = a;
		double sumP = 0;

		// Here we iteratively compute the one-tailed p-value: the sum of all intermediate p-values
		// while iteratively shifting counts from the a/d diagonal to b/c.
		double p = computePval(a, b, c, d);
		double p_1 = p;

		while (a >= 0) {
			sumP += p;
			if (a == 0) {
				break;
			}
			a--;
			b++;
			c++;
			d--;
			p = computePval(a, b, c, d);
		}
		if (oneTailed) {
			return Math.min(sumP, 1);
		}

		// For the two-tailed statistic we swap the table diagonals and repeat the iterative process
		a = b;
		b = 0;
		c = c - a;
		d = d + a;
		p = computePval(a, b, c, d);

		while (p < p_1) {
			if (a == original_a) {
				break;
			}
			sumP += p;
			a--;
			b++;
			c++;
			d--;
			p = computePval(a, b, c, d);
		}
		return Math.min(sumP, 1);
	}

	/**
	 * @return A single iteration of the p-value computation for the given counts
	 */
	private static double computePval(int a, int b, int c, int d) {
		return Math.exp(
										(CACHED_LOG_FACTORIAL.get(a + b) +
										 CACHED_LOG_FACTORIAL.get(c + d) +
										 CACHED_LOG_FACTORIAL.get(a + c) +
										 CACHED_LOG_FACTORIAL.get(b + d))
										-
										(CACHED_LOG_FACTORIAL.get(a + b + c + d) +
										 CACHED_LOG_FACTORIAL.get(a) +
										 CACHED_LOG_FACTORIAL.get(b) +
										 CACHED_LOG_FACTORIAL.get(c) +
										 CACHED_LOG_FACTORIAL.get(d)));
	}

	/**
	 * Check if we have enough cached log factorial values for the given table counts. If not, expand
	 * the table cache.
	 */
	private static void ensureCacheSize(int n1, int n2, int n3, int n4) {
		int minSize = 1 + n1 + n2 + n3 + n4;
		if (CACHED_LOG_FACTORIAL.size() < minSize) {
			buildCachedEntries(minSize);
		}
	}

	/**
	 * @param minSize Number of log factorial entries to ensure are cached in our table
	 */
	private static void buildCachedEntries(int minSize) {
		synchronized (CACHED_LOG_FACTORIAL) {
			if (CACHED_LOG_FACTORIAL.isEmpty()) {
				// Seed the initial entry
				CACHED_LOG_FACTORIAL.add(0.0);
			}
			for (int i = CACHED_LOG_FACTORIAL.size(); i <= minSize; i++) {
				CACHED_LOG_FACTORIAL.add(CACHED_LOG_FACTORIAL.get(i - 1) + Math.log(i));
			}
		}
	}

	public static void main(String[] args) {
		double chiP, fishP;
		int[][][] testMatrices = new int[][][] {
																						{{0, 20}, {0, 25}},
																						{{2, 20}, {3, 25}},
																						{{0, 20}, {5, 20}},
																						{{0, 2100}, {0, 2500}},
																						{{0, 2100}, {5, 2500}},
																						{{0, 2100}, {10, 2500}},
																						{{210, 2100}, {250, 2500}},
		};

		System.out.println("chi^2 p\tFisher's Exact p");
		for (int i = 0; i < testMatrices.length; i++) {
			chiP = ProbDist.ChiDist(ContingencyTable.ChiSquare(testMatrices[i], false), 1);
			fishP = FishersExact2by2Calculator.getPvalue(testMatrices[i]);
			System.out.println(ext.prettyP(chiP) + "\t" + fishP);
		}
	}
}
