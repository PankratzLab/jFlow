package org.genvisis.common;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.genvisis.cnv.util.Java6Helper;

public class Sort {

	/**
	 * As {@link #getSortedIndices(long[])} but returns in reverse order.
	 */
	public static int[] getReverseIndices(long[] arr) {
		return Array.reverseInPlace(getSortedIndices(arr));
	}

	/**
	 * Finds array of indices that would correspond to the sorted order of the given array.
	 *
	 * @param arr array of longs to sort
	 * @return array of ordered indices
	 */
	public static int[] getSortedIndices(long[] arr) {
		return getSortedIndices(new LongArrayWrapper(arr));
	}

	/**
	 * Finds array of indices that would correspond to the sorted order of the given array.
	 *
	 * @param arr array of longs to sort
	 * @return array of ordered indices
	 */
	public static int[] getSortedIndices(float[] arr) {
		return getSortedIndices(new FloatArrayWrapper(arr));
	}

	/**
	 * As {@link #getSortedIndices(int[])} but returns in reverse order.
	 */
	public static int[] getReverseIndices(int[] arr) {
		return Array.reverseInPlace(getSortedIndices(arr));
	}

	/**
	 * Finds array of indices that would correspond to the sorted order of the given array.
	 *
	 * @param arr array of integers to sort
	 * @return array of ordered indices
	 */
	public static int[] getSortedIndices(int[] arr) {
		return getSortedIndices(new IntArrayWrapper(arr));
	}

	/**
	 * As {@link #getSortedIndices(double[])} but returns in reverse order.
	 */
	public static int[] getReverseIndices(double[] arr) {
		return Array.reverseInPlace(getSortedIndices(arr));
	}

	/**
	 * Finds array of indices that would correspond to the sorted order of the given array.
	 *
	 * @param arr array of doubles to sort
	 * @return array of ordered indices
	 */
	public static int[] getSortedIndices(double[] arr) {
		return getSortedIndices(new DoubleArrayWrapper(arr));
	}

	public static int[] getSortedIndices(String[] arr) {
		return getSortedIndices(new StringArrayWrapper(arr));
	}

	public static <T extends Comparable<T>> int[] getReverseIndices(List<T> list) {
		return Array.reverseInPlace(getSortedIndices(list));
	}

	public static <T extends Comparable<T>> int[] getSortedIndices(List<T> list) {
		return getSortedIndices(new ComparableListWrapper<T>(list));
	}

	public static int[] getSort2DIndices(byte[] bytes, int[] ints) {
		IndexedByteInt[] ibi = new IndexedByteInt[bytes.length];
		for (int i = 0; i < bytes.length; i++) {
			ibi[i] = new IndexedByteInt(bytes[i], ints[i], i);
		}
		return getSort2DIndices(ibi);
	}

	public static int[] getSort2DIndices(List<int[]> matrix) {
		return getSort2DIndices(Matrix.extractColumn(matrix, 0), Matrix.extractColumn(matrix, 1));
	}

	public static int[] getSort2DIndices(int[] int1s, int[] int2s) {
		IndexedIntInt[] iii = new IndexedIntInt[int1s.length];
		for (int i = 0; i < int1s.length; i++) {
			iii[i] = new IndexedIntInt(int1s[i], int2s[i], i);
		}
		return getSort2DIndices(iii);
	}

	public static <C1 extends Comparable<C1>, C2 extends Comparable<C2>> int[] getSort2DIndices(List<C1> list1,
	                                                                                            List<C2> list2) {
		@SuppressWarnings("unchecked")
		IndexedCC<C1, C2>[] icc = new IndexedCC[list1.size()];
		for (int i = 0; i < icc.length; i++) {
			icc[i] = new IndexedCC<C1, C2>(list1.get(i), list2.get(i), i);
		}
		return getSort2DIndices(icc);
	}

	private static int[] getSortedIndices(SortableWrapper w) {
		int[] indices = new int[w.size()];
		if (w.size() == 0) {
		  return indices;
		}
		for (int i = 1; i < indices.length; i++) {
			indices[i] = i;
		}

		quicksort(indices, w, 0, indices.length-1);

		return indices;
	}

	/**
	 * from http://www.programcreek.com/2012/11/quicksort-array-in-java/
	 */
	private static void quicksort(int[] indices, SortableWrapper w, int lowIndex, int highIndex) {
		int i = lowIndex;
		int j = highIndex;
		// calculate pivot number
		w.markPivot(indices[(lowIndex + highIndex) / 2]);
		// Divide into two arrays
		while (i <= j) {
			while (w.ltPivot(indices[i])) {
				i++;
			}
			while (w.gtPivot(indices[j])) {
				j--;
			}
			if (i <= j) {
				exchangeNumbers(indices, i, j);
				// move index to next position on both sides
				i++;
				j--;
			}
		}
		// recurse to subarrays
		if (lowIndex < j) {
			quicksort(indices, w, lowIndex, j);
		}
		if (i < highIndex) {
			quicksort(indices, w, i, highIndex);
		}
	}

	private static int[] getSort2DIndices(IndexedWrapper[] wrappers) {
		Arrays.sort(wrappers);
		int[] indices = new int[wrappers.length];
		for (int i = 0; i < wrappers.length; i++) {
			indices[i] = wrappers[i].getIndex();
		}

		return indices;
	}

	/**
	 * Swaps the positions of data in an array
	 *
	 * @param array the array containg the data to swap
	 * @param index1 position 1 to swap
	 * @param index2 position 2 to swap
	 */
	private static void exchangeNumbers(int[] array, int index1, int index2) {
		int temp = array[index1];
		array[index1] = array[index2];
		array[index2] = temp;
	}

	/**
	 * Copies the given list and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static String[] getOrdered(List<String> unordered, int[] indices) {
		String[] ordered = new String[unordered.size()];

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered.get(indices[i]);
		}

		return ordered;
	}

	/**
	 * Copies the given array and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static String[] getOrdered(String[] unordered, int[] indices) {
		String[] ordered = new String[unordered.length];

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered[indices[i]];
		}

		return ordered;
	}

	/**
	 * Copies the given array and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static byte[] getOrdered(byte[] unordered, int[] indices) {
		byte[] ordered = new byte[unordered.length];

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered[indices[i]];
		}

		return ordered;
	}

	/**
	 * Copies the given array and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static int[] getOrdered(int[] unordered, int[] indices) {
		int[] ordered = new int[unordered.length];

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered[indices[i]];
		}

		return ordered;
	}

	/**
	 * Copies the given array and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static float[] getOrdered(float[] unordered, int[] indices) {
		float[] ordered = new float[unordered.length];

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered[indices[i]];
		}

		return ordered;
	}

	/**
	 * Copies the given array and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static double[] getOrdered(double[] unordered, int[] indices) {
		double[] ordered = new double[unordered.length];

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered[indices[i]];
		}

		return ordered;
	}

	/**
	 * Copies the given array and orders it using the specified indices
	 *
	 * @return the ordered copy
	 */
	public static <T> T[] getOrdered(T[] unordered, int[] indices) {
		T[] ordered = Arrays.copyOf(unordered, unordered.length);

		for (int i = 0; i < indices.length; i++) {
			ordered[i] = unordered[indices[i]];
		}

		return ordered;
	}

	public static void reverseSort(String[] array) {
		Arrays.sort(array);
		Array.reverseInPlace(array);
	}

	public static void reverseSort(long[] array) {
		Arrays.sort(array);
		Array.reverseInPlace(array);
	}

	public static void reverseSort(int[] array) {
		Arrays.sort(array);
		Array.reverseInPlace(array);
	}

	/**
	 * Ranks each element in the given array. Ranks will be returned in the sorted (or reverse sorted)
	 * order. For example: input of [4, 4, 6, 5] would produce ranks of [1, 1, 2, 3}]
	 *
	 * @param array
	 * @param reverse if true, ranks will be reversed
	 * @return array of ranks
	 */
	public static int[] getRanks(double[] array, boolean reverse) {
		int[] order = getSortedIndices(array);

		int[] ranks = new int[array.length];
		int curRank = 1;
		ranks[0] = curRank;

		for (int i = 1; i < array.length; i++) {
			if (Math.abs(array[order[i]] - array[order[i - 1]]) > 0.00001) {
				curRank++;
			}
			ranks[i] = curRank;
		}

		if (reverse) {
			Array.reverseInPlace(ranks);
		}

		return ranks;

	}

	/**
	 * Copied from http://stackoverflow.com/a/11648106/1027800
	 *
	 * @return A list of entries in the given map, sorted by their value.
	 */
	public static <K, V extends Comparable<? super V>> List<Entry<K, V>> entriesSortedByValues(Map<K, V> map) {

		List<Entry<K, V>> sortedEntries = new ArrayList<Entry<K, V>>(map.entrySet());

		Collections.sort(sortedEntries, new Comparator<Entry<K, V>>() {
			@Override
			public int compare(Entry<K, V> e1, Entry<K, V> e2) {
				return e2.getValue().compareTo(e1.getValue());
			}
		});

		return sortedEntries;
	}

	/**
	 * Store a byte, Integer pair with an external index, to facilitate sorting by byte first then
	 * int, and finally the original index
	 *
	 * @author lane0212
	 */
	private static class IndexedByteInt implements Comparable<IndexedByteInt>, IndexedWrapper {
		private final byte b;
		private final int i;
		private final int index;

		public IndexedByteInt(byte b, int i, int index) {
			super();
			this.b = b;
			this.i = i;
			this.index = index;
		}

		public int getIndex() {
			return index;
		}

		@Override
		public int compareTo(IndexedByteInt other) {
			int c = b - other.b;
			if (c == 0) {
				c = Java6Helper.compare(i, other.i);
			}
			if (c == 0) {
				c = Java6Helper.compare(index, other.index);
			}
			return c;
		}
	}

	/**
	 * Store two comparable objects with an external index, to facilitate sorting by the first object,
	 * then the second, and finally the original index
	 */
	private static class IndexedCC<C1 extends Comparable<C1>, C2 extends Comparable<C2>>
	                              implements Comparable<IndexedCC<C1, C2>>, IndexedWrapper {
		private final C1 c1;
		private final C2 c2;
		private final int index;

		public IndexedCC(C1 c1, C2 c2, int index) {
			super();
			this.c1 = c1;
			this.c2 = c2;
			this.index = index;
		}

		public int getIndex() {
			return index;
		}

		@Override
		public int compareTo(IndexedCC<C1, C2> other) {
			int c = c1.compareTo(other.c1);
			if (c == 0) {
				c = c2.compareTo(other.c2);
			}
			if (c == 0) {
				c = Java6Helper.compare(index, other.index);
			}
			return c;
		}
	}

	/**
	 * Store a Integer, Integer pair with an external index, to facilitate sorting by the first int
	 * then the second, and finally the original index
	 */
	private static class IndexedIntInt implements Comparable<IndexedIntInt>, IndexedWrapper {
		private final int i1;
		private final int i2;
		private final int index;

		public IndexedIntInt(int i1, int i2, int index) {
			super();
			this.i1 = i1;
			this.i2 = i2;
			this.index = index;
		}

		public int getIndex() {
			return index;
		}

		@Override
		public int compareTo(IndexedIntInt other) {
			int c = Java6Helper.compare(i1, other.i1);
			if (c == 0) {
				c = Java6Helper.compare(i2, other.i2);
			}
			if (c == 0) {
				c = Java6Helper.compare(index, other.index);
			}
			return c;
		}
	}

	/**
	 * Interface for wrapper classes that are used to sort two or more primitives
	 */
	private static interface IndexedWrapper {
		int getIndex();
	}

	private static class ComparableListWrapper<T extends Comparable<T>> implements SortableWrapper {
		private final List<T> list;
		private T a;

		public ComparableListWrapper(List<T> l) {
			list = l;
		}

		@Override
		public int size() {
			return list.size();
		}

		@Override
		public void markPivot(int index) {
			a = list.get(index);
		}

		@Override
		public boolean gtPivot(int index) {
			return list.get(index).compareTo(a) > 0;
		}

		@Override
		public boolean ltPivot(int index) {
			return list.get(index).compareTo(a) < 0;
		}
	}

	private static class IntArrayWrapper implements SortableWrapper {
		private final int[] array;
		private int a;

		public IntArrayWrapper(int[] arr) {
			array = arr;
		}

		@Override
		public int size() {
			return array.length;
		}

		@Override
		public void markPivot(int index) {
			a = array[index];
		}

		@Override
		public boolean gtPivot(int index) {
			return array[index] > a;
		}

		@Override
		public boolean ltPivot(int index) {
			return array[index] < a;
		}
	}

	private static class FloatArrayWrapper implements SortableWrapper {
		private final float[] array;
		private float a;

		public FloatArrayWrapper(float[] arr) {
			array = arr;
		}

		@Override
		public int size() {
			return array.length;
		}

		@Override
		public void markPivot(int index) {
			a = array[index];
		}

		@Override
		public boolean gtPivot(int index) {
			return array[index] > a;
		}

		@Override
		public boolean ltPivot(int index) {
			return array[index] < a;
		}
	}

	private static class LongArrayWrapper implements SortableWrapper {
		private final long[] array;
		private long a;

		public LongArrayWrapper(long[] arr) {
			array = arr;
		}

		@Override
		public int size() {
			return array.length;
		}

		@Override
		public void markPivot(int index) {
			a = array[index];
		}

		@Override
		public boolean gtPivot(int index) {
			return array[index] > a;
		}

		@Override
		public boolean ltPivot(int index) {
			return array[index] < a;
		}
	}

	private static class DoubleArrayWrapper implements SortableWrapper {
		private final double[] array;
		private double a;

		public DoubleArrayWrapper(double[] arr) {
			array = arr;
		}

		@Override
		public int size() {
			return array.length;
		}

		@Override
		public void markPivot(int index) {
			a = array[index];
		}

		@Override
		public boolean gtPivot(int index) {
			return array[index] > a;
		}

		@Override
		public boolean ltPivot(int index) {
			return array[index] < a;
		}
	}

	private static class StringArrayWrapper implements SortableWrapper {
		private final String[] array;
		private String a;
		private SciStringComparator ac = new SciStringComparator();

		public StringArrayWrapper(String[] arr) {
			array = arr;
		}

		@Override
		public int size() {
			return array.length;
		}

		@Override
		public void markPivot(int index) {
			a = array[index];
		}

		@Override
		public boolean gtPivot(int index) {
			return ac.compare(array[index], a) > 0;
		}

		@Override
		public boolean ltPivot(int index) {
			return ac.compare(array[index], a) < 0;
		}
	}

	/**
	 * Helper class to generalize sorting algorithms.
	 */
	private static interface SortableWrapper {

		/**
		 * @return Size of wrapped collection
		 */
		int size();

		/**
		 * Store the value at the specified index as the current pivot value
		 */
		void markPivot(int index);

		/**
		 * @return true if the value at the specified index is greater than the marked pivot
		 */
		boolean gtPivot(int index);

		/**
		 * @return true if the the value at the specified index is less than the marked pivot
		 */
		boolean ltPivot(int index);
	}

	public static void main(String... args) {
		System.out.println(Double.parseDouble("asdfas"));
	}
}
