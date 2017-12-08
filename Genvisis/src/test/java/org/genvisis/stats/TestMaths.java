package org.genvisis.stats;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

/**
 * Tests for the {@link Maths} class
 */
public class TestMaths {

	private static final List<Integer> POWERS_OF_10 = ImmutableList.of(1, 10, 100, 1000, 10000);
	private static final List<Integer> NOT_POWERS_OF_10 = ImmutableList.of(0, -1, -100, -2, 2, 15, 20,
																																				 110, 111, 101);

	/**
	 * Ensure that {@link Maths#isPowerOf10(int)} works as expected
	 * 
	 */
	@Test
	public void isPowerOf10Test() {
		for (int powerOf10 : POWERS_OF_10) {
			assertTrue(powerOf10 + " is a power of 10", Maths.isPowerOf10(powerOf10));
		}
		for (int notPowerOf10 : NOT_POWERS_OF_10) {
			assertFalse(notPowerOf10 + " is not a power of 10", Maths.isPowerOf10(notPowerOf10));
		}
	}

}
