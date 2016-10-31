package org.genvisis.common;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit tests for {@link Sort}.
 */
public class TestSort {

	/**
	 * Test {@link Sort#getReverseIndices} and {@link Sort#getSortedIndices}
	 */
	@Test
	public void indexTest() {
		int[] ints = {5, 7, 0, 3, 4, 1, 8, 2, 6, 9};

		// Reverse indices test
		int[] indices = Sort.getReverseIndices(ints);
		for (int i=0; i<indices.length; i++) {
			Assert.assertEquals(9 - i, ints[indices[i]]);
		}

		// Forward indices test
		indices = Sort.getSortedIndices(ints);
		for (int i=0; i<indices.length; i++) {
			Assert.assertEquals(i, ints[indices[i]]);
		}

		// Apply ordering
		int[] ordered = Sort.getOrdered(ints, indices);
		for (int i=0; i<ordered.length; i++) {
			Assert.assertEquals(i, ordered[i]);
		}
	}

	/**
	 * Test {@link Sort#getRanks}.
	 */
	@Test
	public void ranksTest() {
		double[] doubles = {40, 30, 50, 30, 30, 40, 10, 20, 10, 20};
		int[] expected = {1, 1, 2, 2, 3, 3, 3, 4, 4, 5};
		int[] ranks = Sort.getRanks(doubles, false);
		Assert.assertArrayEquals(expected, ranks);
	}

	/**
	 * Test {@link Sort#getSortedIndices}.
	 */
	@Test
	public void sort2dTest() {
		// ensure using secondary indices works for sorting
		int[] ints1 = {4, 2, 2, 1};
		int[] ints2 = {1, 2, 1, 4};
		int[] expected = {3, 2, 1, 0};
		int[] indices = Sort.getSort2DIndices(ints1, ints2);
		Assert.assertArrayEquals(expected, indices);
	}

	/**
		* Ensure strings can be sorted as numbers (with {@link SciStringComparator}
	 */
	@Test
	public void alphaTest() {
		String[] strings = {"1e3", "10000", "10"};
		int[] expected = {2, 0, 1};
		int[] indices = Sort.getSortedIndices(strings);
		Assert.assertArrayEquals(expected, indices);
	}
}
