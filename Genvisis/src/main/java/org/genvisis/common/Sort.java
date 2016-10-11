package org.genvisis.common;

import java.text.ParseException;
import java.text.RuleBasedCollator;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;

public class Sort {
	// Constants to define the order
	public static final int ASCENDING = 0;
	public static final int DESCENDING = 1;

	// Constants to define the comparison
	public static final int GREATER = 0;
	public static final int LESSTHAN = 1;
	public static final int GREATEROREQUAL = 2;
	public static final int LESSOREQUAL = 3;

	public static final String RULES = "<0<1<2<3<4<5<6<7<8<9"
																				+ "<a,A<b,B<c,C<d,D<e,E<f,F<g,G<h,H<i,I<j,J"
																			+ "<k,K<l,L<m,M<n,N<o,O<p,P<q,Q<r,R<s,S<t,T"
																			+ "<u,U<v,V<w,W<x,X<y,Y<z,Z";

	/**
	 * Reverses the order of an array key
	 *
	 * @param order array key
	 * @return reverse order of array key
	 */
	private static int[] reverse(int[] order) {
		int[] reverseOrder = new int[order.length];
		for (int i = 0; i < order.length; i++) {
			reverseOrder[i] = order[order.length - 1 - i];
		}
		return reverseOrder;
	}

	// /**
	// * Performs a quicksort or alphabetic sort of a Vector. The vector can be
	// * Strings, Integers, or Longs
	// *
	// * @param v
	// * the vector of data to sort
	// * @param direction
	// * either ASCENDING or DESCENDING
	// * @return array of ordered indices
	// */
	// public static int[] quicksort(Vector<String> v, int direction) {
	// return direction==ASCENDING?quicksort(v):reverse(quicksort(v));
	// }
	//
	// /**
	// * Performs a quicksort or alphabetic sort of a Vector. The vector can be
	// * Strings, Integers, or Longs
	// *
	// * @param v
	// * the vector of data to sort
	// * @return array of ordered indices
	// */
	// public static int[] quicksort(Vector<String> v) {
	// if (v==null||v.size()<=0) {
	// return new int[0];
	// }
	// if (v.elementAt(0) instanceof String) {
	// String[] s = new String[v.size()];
	// v.copyInto(s);
	// return quicksort(s);
	// }
	// if (v.elementAt(0) instanceof Integer) {
	// int[] num = new int[v.size()];
	// for (int i = 0; i<num.length; i++) {
	// num[i] = ((Integer)v.elementAt(i)).intValue();
	// }
	// return quicksort(num);
	// }
	// if (v.elementAt(0) instanceof Long) {
	// long[] num = new long[v.size()];
	// for (int i = 0; i<num.length; i++) {
	// num[i] = ((Long)v.elementAt(i)).longValue();
	// }
	// return quicksort(num);
	// }
	// return new int[0];
	// }

	/**
	 * Performs an alphabetic sort on the data
	 *
	 * @param s array of strings to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(String[] s, int direction) {
		return direction == ASCENDING ? quicksort(s) : reverse(quicksort(s));
	}

	/**
	 * Performs a quicksort on the data
	 *
	 * @param arr array of longs to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(float arr[], int direction) {
		return direction == ASCENDING ? quicksort(arr) : reverse(quicksort(arr));
	}

	/**
	 * Performs a quicksort on the data. Returns the data in ASCENDING order.
	 *
	 * @param arr array of longs to sort
	 * @return array of ordered indices
	 */
	public static int[] quicksort(float[] arr) {
		int n = arr.length;
		int[] indx = new int[arr.length];
		int i, indxt, ir = n - 1, j, k, l = 0;
		int jstack = 0, istack[];
		float a;

		istack = new int[arr.length + 1];
		for (j = 0; j < n; j++) {
			indx[j] = j;
		}
		for (;;) {
			if (ir - l < 7) {
				for (j = l + 1; j <= ir; j++) {
					indxt = indx[j];
					a = arr[indxt];
					for (i = j - 1; i >= 0; i--) {
						if (arr[indx[i]] <= a) {
							break;
						}
						indx[i + 1] = indx[i];
					}
					indx[i + 1] = indxt;
				}
				if (jstack == 0) {
					break;
				}
				ir = istack[jstack--];
				l = istack[jstack--];
			} else {
				k = (l + ir) >> 1;
				SWAP(indx, k, l + 1);
				if (arr[indx[l + 1]] > arr[indx[ir]]) {
					SWAP(indx, l + 1, ir);
				}
				if (arr[indx[l]] > arr[indx[ir]]) {
					SWAP(indx, l, ir);
				}
				if (arr[indx[l + 1]] > arr[indx[l]]) {
					SWAP(indx, l + 1, l);
				}
				i = l + 1;
				j = ir;
				indxt = indx[l];
				a = arr[indxt];
				for (;;) {
					do {
						i++;
					} while (arr[indx[i]] < a);
					do {
						j--;
					} while (arr[indx[j]] > a);
					if (j < i) {
						break;
					}
					SWAP(indx, i, j);
				}
				indx[l] = indx[j];
				indx[j] = indxt;
				jstack += 2;
				if (ir - i + 1 >= j - l) {
					istack[jstack] = ir;
					istack[jstack - 1] = i;
					ir = j - 1;
				} else {
					istack[jstack] = j - 1;
					istack[jstack - 1] = l;
					l = i;
				}
			}
		}
		return indx;
	}

	/**
	 * Performs a quicksort on the data
	 *
	 * @param arr array of longs to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(long arr[], int direction) {
		return direction == ASCENDING ? quicksort(arr) : reverse(quicksort(arr));
	}

	/**
	 * Performs a quicksort on the data. Returns the data in ASCENDING order.
	 *
	 * @param arr array of longs to sort
	 * @return array of ordered indices
	 */
	public static int[] quicksort(long[] arr) {
		int n = arr.length;
		int[] indx = new int[arr.length];
		int i, indxt, ir = n - 1, j, k, l = 0;
		int jstack = 0, istack[];
		long a;

		istack = new int[arr.length + 1];
		for (j = 0; j < n; j++) {
			indx[j] = j;
		}
		for (;;) {
			if (ir - l < 7) {
				for (j = l + 1; j <= ir; j++) {
					indxt = indx[j];
					a = arr[indxt];
					for (i = j - 1; i >= 0; i--) {
						if (arr[indx[i]] <= a) {
							break;
						}
						indx[i + 1] = indx[i];
					}
					indx[i + 1] = indxt;
				}
				if (jstack == 0) {
					break;
				}
				ir = istack[jstack--];
				l = istack[jstack--];
			} else {
				k = (l + ir) >> 1;
				SWAP(indx, k, l + 1);
				if (arr[indx[l + 1]] > arr[indx[ir]]) {
					SWAP(indx, l + 1, ir);
				}
				if (arr[indx[l]] > arr[indx[ir]]) {
					SWAP(indx, l, ir);
				}
				if (arr[indx[l + 1]] > arr[indx[l]]) {
					SWAP(indx, l + 1, l);
				}
				i = l + 1;
				j = ir;
				indxt = indx[l];
				a = arr[indxt];
				for (;;) {
					do {
						i++;
					} while (arr[indx[i]] < a);
					do {
						j--;
					} while (arr[indx[j]] > a);
					if (j < i) {
						break;
					}
					SWAP(indx, i, j);
				}
				indx[l] = indx[j];
				indx[j] = indxt;
				jstack += 2;
				if (ir - i + 1 >= j - l) {
					istack[jstack] = ir;
					istack[jstack - 1] = i;
					ir = j - 1;
				} else {
					istack[jstack] = j - 1;
					istack[jstack - 1] = l;
					l = i;
				}
			}
		}
		return indx;
	}

	/**
	 * Performs a quicksort on the data.
	 *
	 * @param arr array of integers to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(int arr[], int direction) {
		return direction == ASCENDING ? quicksort(arr) : reverse(quicksort(arr));
	}

	/**
	 * Performs a quicksort on the data. Returns the data in ASCENDING order.
	 *
	 * @param arr array of integers to sort
	 * @return array of ordered indices
	 */
	public static int[] quicksort(int[] arr) {
		int n = arr.length;
		int[] indx = new int[arr.length];
		int i, indxt, ir = n - 1, j, k, l = 0;
		int jstack = 0, istack[];
		int a;

		istack = new int[arr.length + 1];
		for (j = 0; j < n; j++) {
			indx[j] = j;
		}
		for (;;) {
			if (ir - l < 7) {
				for (j = l + 1; j <= ir; j++) {
					indxt = indx[j];
					a = arr[indxt];
					for (i = j - 1; i >= 0; i--) {
						if (arr[indx[i]] <= a) {
							break;
						}
						indx[i + 1] = indx[i];
					}
					indx[i + 1] = indxt;
				}
				if (jstack == 0) {
					break;
				}
				ir = istack[jstack--];
				l = istack[jstack--];
			} else {
				k = (l + ir) >> 1;
				SWAP(indx, k, l + 1);
				if (arr[indx[l + 1]] > arr[indx[ir]]) {
					SWAP(indx, l + 1, ir);
				}
				if (arr[indx[l]] > arr[indx[ir]]) {
					SWAP(indx, l, ir);
				}
				if (arr[indx[l + 1]] > arr[indx[l]]) {
					SWAP(indx, l + 1, l);
				}
				i = l + 1;
				j = ir;
				indxt = indx[l];
				a = arr[indxt];
				for (;;) {
					do {
						i++;
					} while (arr[indx[i]] < a);
					do {
						j--;
					} while (arr[indx[j]] > a);
					if (j < i) {
						break;
					}
					SWAP(indx, i, j);
				}
				indx[l] = indx[j];
				indx[j] = indxt;
				jstack += 2;
				if (ir - i + 1 >= j - l) {
					istack[jstack] = ir;
					istack[jstack - 1] = i;
					ir = j - 1;
				} else {
					istack[jstack] = j - 1;
					istack[jstack - 1] = l;
					l = i;
				}
			}
		}
		return indx;
	}

	/**
	 * Sort by doubles first, then by strings
	 *
	 */
	private static class ScoreDoubleStringIndex implements Comparable<ScoreDoubleStringIndex> {
		final double score;
		final String s;
		final int index;

		ScoreDoubleStringIndex(double score, String s, int index) {
			super();
			this.score = score;
			this.s = s;
			this.index = index;
		}

		public int getIndex() {
			return index;
		}

		@Override
		public int compareTo(ScoreDoubleStringIndex o) {
			int cmp = Double.compare(score, o.score);
			if (cmp == 0 && !s.equalsIgnoreCase(o.s)) {
				cmp = s.compareToIgnoreCase(o.s);
			}
			return cmp == 0 ? new Integer(index).compareTo(o.index) : cmp;
		}

	}

	private static class ScoreStringComp implements Comparator<ScoreDoubleStringIndex> {

		@Override
		public int compare(ScoreDoubleStringIndex o1, ScoreDoubleStringIndex o2) {
			int cmp = Double.compare(o1.score, o2.score);
			if (cmp == 0 && !o1.s.equalsIgnoreCase(o2.s)) {
				cmp = o1.s.compareToIgnoreCase(o2.s);
			}
			return cmp == 0 ? new Integer(o1.index).compareTo(o2.index) : cmp;
		}

	}

	/**
	 * uses {@link ScoreStringComp} to sort the double array first, then the string array within the
	 * double order, and get indices
	 *
	 * @param results
	 * @return
	 */
	public static int[] trickSort(double[] results, String[] s) {
		if (results.length != s.length) {
			throw new IllegalArgumentException("Double array and String array must be same size");
		}
		ScoreDoubleStringIndex[] sIndexs = new ScoreDoubleStringIndex[results.length];
		for (int i = 0; i < results.length; i++) {
			sIndexs[i] = new ScoreDoubleStringIndex(results[i], s[i], i);
		}
		Arrays.sort(sIndexs, new ScoreStringComp());
		int[] indexes = new int[results.length];
		for (int i = 0; i < sIndexs.length; i++) {
			indexes[i] = sIndexs[i].getIndex();
		}
		return indexes;
	}

	/**
	 * For sorting doubles and keeping the index, taken from
	 * http://stackoverflow.com/questions/14186529/java-array-of-sorted-indexes
	 */
	private static class ScoreDoubleIndex implements Comparable<ScoreDoubleIndex> {
		final double score;
		final int index;

		ScoreDoubleIndex(double score, int index) {
			this.score = score;
			this.index = index;
		}

		public double getScore() {
			return score;
		}

		public int getIndex() {
			return index;
		}

		@Override
		public int compareTo(ScoreDoubleIndex o) {
			int cmp = Double.compare(score, o.score);
			return cmp == 0 ? new Integer(index).compareTo(o.index) : cmp;
		}
	}

	private static class ScoreComp implements Comparator<ScoreDoubleIndex> {

		@Override
		public int compare(ScoreDoubleIndex o1, ScoreDoubleIndex o2) {
			int cmp = Double.compare(o1.getScore(), o2.getScore());
			return cmp == 0 ? new Integer(o1.getIndex()).compareTo(o2.getIndex()) : cmp;
		}

	}



	/**
	 * uses {@link ScoreDoubleIndex} to sort the double array and get indices
	 *
	 * @param results
	 * @return
	 */
	public static int[] trickSort(double[] results) {
		ScoreDoubleIndex[] sIndexs = new ScoreDoubleIndex[results.length];
		for (int i = 0; i < results.length; i++) {
			sIndexs[i] = new ScoreDoubleIndex(results[i], i);
		}
		Arrays.sort(sIndexs, new ScoreComp());
		int[] indexes = new int[results.length];
		for (int i = 0; i < sIndexs.length; i++) {
			indexes[i] = sIndexs[i].getIndex();
		}
		return indexes;
	}

	/**
	 * Performs a quicksort on the data.
	 *
	 * @param arr array of doubles to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(double arr[], int direction) {
		return direction == ASCENDING ? quicksort(arr) : reverse(quicksort(arr));
	}

	/**
	 * Performs a quicksort on the data. Returns the data in ASCENDING order.
	 *
	 * @param arr array of doubles to sort
	 * @return array of ordered indices
	 */
	public static int[] quicksort(double[] arr) {
		int n = arr.length;
		int[] indx = new int[arr.length];
		int i, indxt, ir = n - 1, j, k, l = 0;
		int jstack = 0, istack[];
		double a;

		istack = new int[arr.length + 1];
		for (j = 0; j < n; j++) {
			indx[j] = j;
		}
		for (;;) {
			if (ir - l < 7) {
				for (j = l + 1; j <= ir; j++) {
					indxt = indx[j];
					a = arr[indxt];
					for (i = j - 1; i >= 0; i--) {
						if (arr[indx[i]] <= a) {
							break;
						}
						indx[i + 1] = indx[i];
					}
					indx[i + 1] = indxt;
				}
				if (jstack == 0) {
					break;
				}
				ir = istack[jstack--];
				l = istack[jstack--];
			} else {
				k = (l + ir) >> 1;
				SWAP(indx, k, l + 1);
				if (arr[indx[l + 1]] > arr[indx[ir]]) {
					SWAP(indx, l + 1, ir);
				}
				if (arr[indx[l]] > arr[indx[ir]]) {
					SWAP(indx, l, ir);
				}
				if (arr[indx[l + 1]] > arr[indx[l]]) {
					SWAP(indx, l + 1, l);
				}
				i = l + 1;
				j = ir;
				indxt = indx[l];
				a = arr[indxt];
				for (;;) {
					do {
						i++;
					} while (arr[indx[i]] < a);
					do {
						j--;
					} while (arr[indx[j]] > a);
					if (j < i) {
						break;
					}
					SWAP(indx, i, j);
				}
				indx[l] = indx[j];
				indx[j] = indxt;
				jstack += 2;
				if (ir - i + 1 >= j - l) {
					istack[jstack] = ir;
					istack[jstack - 1] = i;
					ir = j - 1;
				} else {
					istack[jstack] = j - 1;
					istack[jstack - 1] = l;
					l = i;
				}
			}
		}
		return indx;
	}

	/**
	 * Performs a quicksort on the data.
	 *
	 * @param iv IntVector to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(IntVector iv, int direction) {
		return direction == ASCENDING ? quicksort(iv) : reverse(quicksort(iv));
	}

	/**
	 * Performs a quicksort on the data. Returns the data in ASCENDING order.
	 *
	 * @param iv IntVector to sort
	 * @return array of ordered indices
	 */
	public static int[] quicksort(IntVector iv) {
		int[] ia = new int[iv.size()];
		for (int i = 0; i < ia.length; i++) {
			ia[i] = iv.elementAt(i);
		}
		return quicksort(ia);
	}

	/**
	 * Performs a quicksort on the data.
	 *
	 * @param dv DoubleVector to sort
	 * @param direction either ASCENDING or DESCENDING
	 * @return array of ordered indices
	 */
	public static int[] quicksort(DoubleVector dv, int direction) {
		return direction == ASCENDING ? quicksort(dv) : reverse(quicksort(dv));
	}

	/**
	 * Performs a quicksort on the data. Returns the data in ASCENDING order.
	 *
	 * @param dv DoubleVector to sort
	 * @return array of ordered indices
	 */
	public static int[] quicksort(DoubleVector dv) {
		double[] da = new double[dv.size()];
		for (int i = 0; i < da.length; i++) {
			da[i] = dv.elementAt(i);
		}
		return quicksort(da);
	}

	/**
	 * Sort an array of Strings into ASCENDING alphabetical order
	 *
	 * @param arr the array of Strings to index
	 * @return an integer array of the sorted indices
	 */
	public static int[] quicksort(String[] arr) {
		int n = arr.length;
		int[] indx = new int[arr.length];
		int i, indxt, ir = n - 1, j, k, l = 0;
		int jstack = 0, istack[];
		String a;

		istack = new int[arr.length + 1];

		RuleBasedCollator theCollation;
		try {
			theCollation = new RuleBasedCollator(RULES);
			for (j = 0; j < n; j++) {
				indx[j] = j;
			}
			for (;;) {
				if (ir - l < 7) {
					for (j = l + 1; j <= ir; j++) {
						indxt = indx[j];
						a = arr[indxt];
						for (i = j - 1; i >= 0; i--) {
							if (compare(theCollation, arr[indx[i]], a, LESSOREQUAL)) {
								break;
							}
							indx[i + 1] = indx[i];
						}
						indx[i + 1] = indxt;
					}
					if (jstack == 0) {
						break;
					}
					ir = istack[jstack--];
					l = istack[jstack--];
				} else {
					k = (l + ir) >> 1;
					SWAP(indx, k, l + 1);
					if (compare(theCollation, arr[indx[l + 1]], arr[indx[ir]], GREATER)) {
						SWAP(indx, l + 1, ir);
					}
					if (compare(theCollation, arr[indx[l]], arr[indx[ir]], GREATER)) {
						SWAP(indx, l, ir);
					}
					if (compare(theCollation, arr[indx[l + 1]], arr[indx[l]], GREATER)) {
						SWAP(indx, l + 1, l);
					}
					i = l + 1;
					j = ir;
					indxt = indx[l];
					a = arr[indxt];
					for (;;) {
						do {
							i++;
						} while (compare(theCollation, arr[indx[i]], a, LESSTHAN));
						do {
							j--;
						} while (compare(theCollation, arr[indx[j]], a, GREATER));
						if (j < i) {
							break;
						}
						SWAP(indx, i, j);
					}
					indx[l] = indx[j];
					indx[j] = indxt;
					jstack += 2;
					if (ir - i + 1 >= j - l) {
						istack[jstack] = ir;
						istack[jstack - 1] = i;
						ir = j - 1;
					} else {
						istack[jstack] = j - 1;
						istack[jstack - 1] = l;
						l = i;
					}
				}
			}

		} catch (ParseException pe) {
			System.err.println("Error initializing RuleBasedCollator for quicksort(String[] arr)");
			pe.printStackTrace();
			System.exit(-1);
		}

		return indx;
	}

	/**
	 * Rule comparitor to compare 2 strings and return if it matches the comparison type
	 *
	 * @param collator the rule collator for the comparison
	 * @param arg1 string 1
	 * @param arg2 string 2
	 * @param type the type of comparison GREATER/LESSTHAN/GREATEROREQUAL/LESSOREQUAL
	 * @return the result of the comparison
	 */
	private static boolean compare(RuleBasedCollator collator, String arg1, String arg2, int type) {
		int compareResult = collator.compare(arg1, arg2);
		switch (type) {
			case (GREATER):
				if (compareResult > 0) {
					return true;
				}
				break;
			case (LESSTHAN):
				if (compareResult < 0) {
					return true;
				}
				break;
			case (GREATEROREQUAL):
				if (compareResult >= 0) {
					return true;
				}
				break;
			case (LESSOREQUAL):
				if (compareResult <= 0) {
					return true;
				}
				break;
		}
		return false;
	}

	/**
	 * Swaps the positions of data in an array
	 *
	 * @param array the array containg the data to swap
	 * @param index1 position 1 to swap
	 * @param index2 position 2 to swap
	 */
	private static void SWAP(int[] array, int index1, int index2) {
		int temp = array[index1];
		array[index1] = array[index2];
		array[index2] = temp;
	}

	/**
	 * Sorts an array of integers and returns the sorted array
	 *
	 * @param array the array to be sorted
	 * @return the sorted array
	 */
	public static int[] putInOrder(int[] array) {
		int[] order = quicksort(array);
		int[] newArray = new int[array.length];

		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	/**
	 * Sorts an array of numbers and returns the sorted array
	 *
	 * @param array the array to be sorted
	 * @return the sorted array
	 */
	public static double[] putInOrder(double[] array) {
		int[] order = quicksort(array);
		double[] newArray = new double[array.length];

		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	/**
	 * Sorts an array of strings and returns the sorted array
	 *
	 * @param array the array to be sorted
	 * @return the sorted array
	 */
	public static String[] putInOrder(String[] array) {
		return putInOrder(array, false);
	}

	/**
	 * Sorts an array of strings and returns the sorted array
	 *
	 * @param array the array to be sorted
	 * @return the sorted array
	 */
	public static String[] putInOrder(String[] array, boolean treatAsNumbers) {
		String[] newArray;
		int[] order;

		if (treatAsNumbers) {
			for (int i = 0; treatAsNumbers && i < array.length; i++) {
				if (!ext.isValidDouble(array[i])) {
					treatAsNumbers = false;
				}
			}
		}

		if (treatAsNumbers) {
			order = quicksort(Array.toDoubleArray(array));
		} else {
			order = quicksort(array);
		}

		newArray = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	/**
	 * Sorts by the first column first and then by the second; returns the order
	 *
	 * @param matrix matrix to be sorted
	 * @return array of sorted indices
	 */
	public static int[] orderTwoLayers(int[][] matrix, Logger log) {
		return orderTwoLayers(Array.toByteArray(Matrix.extractColumn(matrix, 0)),
													Matrix.extractColumn(matrix, 1), log);
	}

	public static int[] orderTwoLayers(byte[] first, int[] second, Logger log) {
		return orderTwoLayers(first, second, false, log);
	}

	/**
	 * Class copied from Java 1.8
	 */
	private static int byteCompare(byte x, byte y) {
		return x - y;
	}

	/**
	 * Class copied from Java 1.8
	 */
	public static int intCompare(int x, int y) {
		return (x < y) ? -1 : ((x == y) ? 0 : 1);
	}


	/**
	 * @author lane0212 {@link Comparator} for {@link BII}
	 *
	 */
	private static class BIIComp implements Comparator<BII> {

		@Override
		public int compare(BII o1, BII o2) {
			// int value1 = Byte.compare(o1.getB(), o2.getB()); // Not compatible with Java version 1.6
			int value1 = byteCompare(o1.getB(), o2.getB());
			if (value1 == 0) {
				// value1 = Integer.compare(o1.getI(), o2.getI()); // Not compatible with Java version 1.6
				value1 = intCompare(o1.getI(), o2.getI());
			}
			// return value1 == 0 ? Integer.compare(o1.getIndex(), o2.getIndex()) : value1; // Not
			// compatible with Java version 1.6
			return value1 == 0 ? intCompare(o1.getIndex(), o2.getIndex()) : value1;
		}
	}

	/**
	 * @author lane0212 Store a byte, Integer, Integer(Index)
	 */
	private static class BII {
		final byte b;
		final int i;
		final int index;

		public BII(byte b, int i, int index) {
			super();
			this.b = b;
			this.i = i;
			this.index = index;
		}

		public byte getB() {
			return b;
		}

		public int getI() {
			return i;
		}

		public int getIndex() {
			return index;
		}

	}

	private static int[] trickSortOrderTwoLayers(	byte[] first, int[] second, boolean verbose,
																								Logger log) {
		BII[] biisSorted = trickSortTwoLayers(first, second, verbose, log);
		if (biisSorted.length > 0) {
			int[] order = new int[biisSorted.length];
			byte tmpFirst = biisSorted[0].getB();
			int tmpSecond = biisSorted[0].getI();
			for (int i = 0; i < order.length; i++) {
				// TODO, these checks are because I am not positive of the Comparator implementation
				if (biisSorted[i].getI() < tmpSecond && biisSorted[i].getB() <= tmpFirst) {
					String error = "Invalid sorting of two layers";
					log.reportTimeError(error);
					System.out.println(tmpSecond + " - > " + biisSorted[i].getI());
					System.out.println(tmpFirst + " - > " + biisSorted[i].getB());

					throw new IllegalStateException(error);
				}
				if (biisSorted[i].getI() > tmpSecond || biisSorted[i].getB() >= tmpFirst) {
					tmpSecond = biisSorted[i].getI();
				}
				if (biisSorted[i].getB() < tmpFirst) {
					String error = "Invalid sorting of two layers";
					log.reportTimeError(error);
					throw new IllegalStateException(error);
				}
				if (biisSorted[i].getB() > tmpFirst) {
					tmpFirst = biisSorted[i].getB();
				}
				order[i] = biisSorted[i].getIndex();
			}
			return order;
		} else {
			return new int[] {};
		}
	}

	private static BII[] trickSortTwoLayers(byte[] first, int[] second, boolean verbose, Logger log) {
		if (first.length != second.length) {
			String error = "Error - Can't sort arrays if the number of entries do not match up";
			log.reportError(error);
			throw new IllegalArgumentException(error);
		} else {
			BII[] biis = new BII[first.length];
			for (int i = 0; i < biis.length; i++) {
				biis[i] = new BII(first[i], second[i], i);
			}
			Arrays.sort(biis, new BIIComp());
			return biis;
		}
	}

	/**
	 * Sorts by the first array first and then by the second; returns the order
	 *
	 * @param first first order array
	 * @param second second order array
	 *
	 * @return array of sorted indices
	 */
	public static int[] orderTwoLayers(byte[] first, int[] second, boolean verbose, Logger log) {
		return orderTwoLayers(first, second, verbose, true, log);

	}

	/**
	 * Sorts by the first array first and then by the second; returns the order
	 *
	 * @param first first order array
	 * @param second second order array
	 * @param trickSort use built in java sorting
	 * @return array of sorted indices
	 */
	public static int[] orderTwoLayers(	byte[] first, int[] second, boolean verbose, boolean trickSort,
																			Logger log) {
		if (trickSort) {
			return trickSortOrderTwoLayers(first, second, verbose, log);
		} else {
			String[] primaryKeys, secondaryKeys;
			int count;
			int[] values, finalIndices, primaryIndices, secondaryIndices;
			Hashtable<String, Hashtable<String, String>> mapFirstToSecond;
			Hashtable<String, String> hash, finalKeyHash;
			boolean inOrder;

			if (first.length != second.length) {
				log.reportError("Error - Can't sort markers if the number of chromosome numbers and positions don't match up");
				System.exit(1);
			}

			inOrder = true;
			for (int i = 1; i < first.length; i++) {
				if (first[i] < first[i - 1] || (second[i] < second[i - 1] && first[i] == first[i - 1])) {
					// log.report("chr"+first[i]+":"+second[i]+" < chr"+first[i-1]+":"+second[i-1]);
					inOrder = false;
				}
			}
			if (inOrder) {
				if (verbose) {
					log.report("Markers were already in order", true, true, 10);
				}
				return Array.arrayOfIndices(second.length);
			}

			mapFirstToSecond = new Hashtable<String, Hashtable<String, String>>();
			for (int i = 0; i < first.length; i++) {
				HashVec.addToHashHash(mapFirstToSecond, first[i] + "", i + "", second[i] + "");
			}

			// if (mapFirstToSecond.size() == 1) {
			// inOrder = true;
			// count = 1;
			// while (inOrder && count < second.length) {
			// if (second[count-1] > second[count]) {
			// inOrder = false;
			// }
			// count++;
			// }
			// if (inOrder) {
			// return Array.intArray(second.length);
			// }
			// }

			count = 0;
			finalKeyHash = new Hashtable<String, String>();
			values = new int[mapFirstToSecond.size()];
			primaryKeys = HashVec.getKeys(mapFirstToSecond);
			for (int i = 0; i < primaryKeys.length; i++) {
				values[i] = Integer.parseInt(primaryKeys[i]);
			}
			primaryIndices = quicksort(values);
			for (int i = 0; i < primaryKeys.length; i++) {
				hash = mapFirstToSecond.get(primaryKeys[primaryIndices[i]]);
				if (hash != null) {
					values = new int[hash.size()];
					secondaryKeys = HashVec.getKeys(hash);
					for (int j = 0; j < secondaryKeys.length; j++) {
						values[j] = Integer.parseInt(hash.get(secondaryKeys[j]));
					}
					secondaryIndices = quicksort(values);
					for (int secondaryIndice : secondaryIndices) {
						finalKeyHash.put(secondaryKeys[secondaryIndice], count + "");
						count++;
					}
				}
			}

			finalIndices = new int[first.length];
			for (int i = 0; i < first.length; i++) {
				finalIndices[Integer.parseInt(finalKeyHash.get(i + ""))] = i;
			}

			return finalIndices;
		}
	}

	/**
	 * Sorts by the first array first and then by the second; returns the order
	 *
	 * @param first first order array
	 * @param second second order array
	 * @return array of sorted indices
	 */
	public static int[] orderTwoLayers(	int[] first, int firstDirection, String[] second,
																			int secondDirection) {
		String[] primaryKeys, secondaryKeys;
		int count;
		int[] primaryValues, finalIndices, primaryIndices, secondaryIndices;
		String[] secondaryValues;
		Hashtable<String, Hashtable<String, String>> mapFirstToSecond;
		Hashtable<String, String> hash, finalKeyHash;

		if (first.length != second.length) {
			System.err.println("Error - Can't sort markers if the number of chromosome numbers and positions don't match up");
			System.exit(1);
		}

		mapFirstToSecond = new Hashtable<String, Hashtable<String, String>>();
		for (int i = 0; i < first.length; i++) {
			HashVec.addToHashHash(mapFirstToSecond, first[i] + "", i + "", second[i] + "");
		}

		count = 0;
		finalKeyHash = new Hashtable<String, String>();
		primaryValues = new int[mapFirstToSecond.size()];
		primaryKeys = HashVec.getKeys(mapFirstToSecond);
		for (int i = 0; i < primaryKeys.length; i++) {
			primaryValues[i] = Integer.parseInt(primaryKeys[i]);
		}
		primaryIndices = quicksort(primaryValues, firstDirection);
		for (int i = 0; i < primaryKeys.length; i++) {
			hash = mapFirstToSecond.get(primaryKeys[primaryIndices[i]]);
			if (hash != null) {
				secondaryValues = new String[hash.size()];
				secondaryKeys = HashVec.getKeys(hash);
				for (int j = 0; j < secondaryKeys.length; j++) {
					secondaryValues[j] = hash.get(secondaryKeys[j]);
				}
				secondaryIndices = quicksort(secondaryValues, secondDirection);
				for (int secondaryIndice : secondaryIndices) {
					finalKeyHash.put(secondaryKeys[secondaryIndice], count + "");
					count++;
				}
			}
		}

		finalIndices = new int[first.length];
		for (int i = 0; i < first.length; i++) {
			finalIndices[Integer.parseInt(finalKeyHash.get(i + ""))] = i;
		}

		return finalIndices;
	}

	public static int[] ranks(double[] array, int direction) {
		int[] ranks;
		int[] order;
		int start;

		order = quicksort(array, direction);
		ranks = new int[array.length];
		start = 0;
		for (int i = 1; i < array.length; i++) {
			if ((direction == ASCENDING && array[order[start]] < array[order[i]])
					|| (direction == DESCENDING && array[order[start]] > array[order[i]])) {
				for (int j = start; j < i; j++) {
					ranks[order[j]] = start + 1;
				}
				start = i;
			}
		}

		return ranks;

	}

	@SuppressWarnings("unchecked")
	public static <T> T[] putInOrder(T[] array, int[] order) {
		if (array == null || array.length == 0) {
			return array;
		}
		T[] newArray = null;

		for (T element : array) {
			if (element != null) {
				newArray = (T[]) java.lang.reflect.Array.newInstance(element.getClass(), array.length);
				break;
			}
		}

		if (newArray == null) {
			return array; // all elements of array are null
		}

		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static byte[] putInOrder(byte[] array, int[] order) {
		byte[] newArray;

		newArray = new byte[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static int[] putInOrder(int[] array, int[] order) {
		int[] newArray;

		newArray = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static double[] putInOrder(double[] array, int[] order) {
		double[] newArray;

		newArray = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static float[] putInOrder(float[] array, int[] order) {
		float[] newArray;

		newArray = new float[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static void main(String[] args) {
		String[][] test = {{"in"}, {"to"}, {"new"}, {"directory"}, {"be"}};
		int[] order = {1, 4, 0, 2, 3};

		putInOrder(test, order);
	}
}
