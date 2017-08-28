package org.genvisis.common;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import org.junit.Assert;

/**
 * Benchmark tests for String index methods in {@link ext} methods.
 */
public class TestIndex {

	private String[] superset;
	private String[][] targets;

	public void setup() {
		// NB: superSize is intended to be the square of targetSize, so that targets can be a NxN array
		final int superSize = 10000;
		final int targetSize = (int) Math.sqrt(superSize);
		superset = new String[superSize];
		targets = new String[targetSize][targetSize];
		// NB: since we are comparing expected hash results of a random array, we need to use a
		// consistent seed
		Random r = new Random(15l);

		// Build the targets and superset arrays
		for (int i = 0; i < targets.length; i++) {
			for (int j = 0; j < targets[i].length; j++) {
				String value = String.valueOf(r.nextInt(superSize));
				targets[i][j] = value;
				superset[(targets.length * i) + j] = value;
			}
		}
	}

	/**
	 * Demonstrate the "optimal" index lookup time using a pure map implementation
	 */
	public void naiveMap() {
		System.out.print("map creation and lookup finished in: ");
		long t = System.currentTimeMillis();
		Map<String, Integer> map = new HashMap<>();
		for (int i = 0; i < superset.length; i++) {
			map.put(superset[i].toLowerCase(), i);
		}
		int[][] result = new int[targets.length][];
		for (int i = 0; i < targets.length; i++) {
			result[i] = new int[targets[i].length];
			for (int j = 0; j < targets[i].length; j++) {
				result[i][j] = map.get(targets[i][j].toLowerCase());
			}
		}
		t = System.currentTimeMillis() - t;
		System.out.println(t + "ms");
	}

	/**
	 * Benchmark {@link ext#indexLargeFactors(String[], String[], boolean, Logger, boolean, boolean)}
	 */
	public void indexLargeFactors() {
		final int expectedHash = -828889133;
		long t = System.currentTimeMillis();
		int[] result = ext.indexLargeFactors(superset, superset, false, new Logger(), false,
																				 false);
		t = System.currentTimeMillis() - t;
		report("ext.indexLargeFactors", result, t, expectedHash);
	}

	/**
	 * Benchmark the {@link ext#indexFactors(String[], String[], boolean, Logger, boolean, boolean)}
	 * method (with a 1D input array)
	 */
	public void indexFactors1D() {
		final int expectedHash = -828889133;
		long t = System.currentTimeMillis();
		int[] result = ext.indexFactors(superset, superset, false, new Logger(), false,
																		false);
		t = System.currentTimeMillis() - t;
		report("ext.indexFactors - 1D", result, t, expectedHash);
	}

	/**
	 * Benchmark the
	 * {@link ext#indexFactors(String[][], String[], boolean, boolean, boolean, boolean, Logger, boolean)}
	 * method (with a 2D input array)
	 */
	public void indexFactors2D() {
		Logger log = new Logger();
		long t = System.currentTimeMillis();
		int[] result = ext.indexFactors(targets, superset, true, true, true, false, log, false);
		t = System.currentTimeMillis() - t;
		report("ext.indexFactors - 2D,target,exact", result, t, -1585103590);
		t = System.currentTimeMillis();
		result = ext.indexFactors(targets, superset, false, true, true, false, log, false);
		t = System.currentTimeMillis() - t;
		report("ext.indexFactors - 2D,super,exact", result, t, -404192925);
		t = System.currentTimeMillis();
		result = ext.indexFactors(targets, superset, true, true, false, false, log, false);
		t = System.currentTimeMillis() - t;
		report("ext.indexFactors - 2D,target,fuzzy", result, t, -842850808);
		t = System.currentTimeMillis();
		result = ext.indexFactors(targets, superset, false, true, false, false, log, false);
		t = System.currentTimeMillis() - t;
		report("ext.indexFactors - 2D,super,fuzzy", result, t, 254970431);
	}

	/**
	 * Helper method for reporting results of benchmark test
	 *
	 * @param method Method name that was tested
	 * @param result Resulting index array from method
	 * @param time How long it took to create the index array
	 * @param expectedHash Hash code of resultant array via original methods.
	 */
	private void report(String method, int[] result, long time, int expectedHash) {
		int hash = Arrays.hashCode(result);
		System.out.print(method + " - finished in: " + time + "ms. Expected hash: " + expectedHash
										 + " - got: " + hash);
		System.out.println();
		Assert.assertEquals(expectedHash, hash);
	}

	/**
	 * Entry point for manual testing
	 */
	public static void main(String... strings) {
		// NB: this class could easily be converted to JUnit tests by annotating:
		// - setup with @Before
		// - the testXXX methods with @Test
		// However, because part of this test involves benchmarking, it is not necessarily appropriate
		// for running with each build.
		TestIndex test = new TestIndex();
		test.setup();
		test.naiveMap();
		test.indexFactors1D();
		test.indexFactors2D();

		// NB: doesn't match output
		// test.indexLargeFactors();
	}
}
