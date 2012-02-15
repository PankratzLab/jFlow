package common;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

import stats.ProbDist;

public class Array {

	/**
	 * Return the minimum in an array of integers
	 * 
	 * @param array
	 *            array of integers
	 * @return the minimum
	 */
	public static int min(int[] array) {
		int min;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return -999;
		}
		min = array[0];
		for (int i = 1; i<array.length; i++) {
			min = Math.min(array[i], min);
		}
		return min;
	}

	/**
	 * Return the minimum in an array of integers
	 * 
	 * @param array
	 *            array of integers
	 * @return the minimum
	 */
	public static double min(double[] array) {
		double min;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return Double.NaN;
		}
		min = array[0];
		for (int i = 1; i<array.length; i++) {
			min = Math.min(array[i], min);
		}
		return min;
	}

	/**
	 * Return the minimum in an array of floats
	 * 
	 * @param array
	 *            array of floats
	 * @return the minimum
	 */
	public static float min(float[] array) {
		float min;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return -999;
		}
		min = array[0];
		for (int i = 1; i<array.length; i++) {
			if (array[i]!=Float.NaN&&array[i]<min) {
				min = array[i];
			}
		}
		return min;
	}

	/**
	 * Return the index of the minimum in an array of floats
	 * 
	 * @param array
	 *            array of floats
	 * @return the index of the minimum
	 */
	public static int indexOfMin(float[] array) {
		float min;
		int index;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return -999;
		}
		min = array[0];
		index = 0;
		for (int i = 1; i<array.length; i++) {
			if (array[i]!=Float.NaN&&array[i]<min) {
				min = array[i];
				index = i;
			}
		}
		return index;
	}

	/**
	 * Return the index of the minimum in an array of integers
	 * 
	 * @param array
	 *            array of integers
	 * @return the minimum
	 */
	public static int minIndex(double[] array) {
		double min;
		int index;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return -1;
		}
		min = array[0];
		index = 0;
		for (int i = 1; i<array.length; i++) {
			if (array[i]<min) {
				min = array[i];
				index = i;
			}
		}
		return index;
	}

	/**
	 * Return the maximum in an array of integers
	 * 
	 * @param array
	 *            array of integers
	 * @return the maximum
	 */
	public static int max(int[] array) {
		int max;

		if (array.length==0) {
			System.err.println("Error - impossible to find the max of an empty array");
			return -999;
		}
		max = array[0];
		for (int i = 1; i<array.length; i++) {
			max = Math.max(array[i], max);
		}
		return max;
	}

	/**
	 * Return the maximum in an array of floats
	 * 
	 * @param array
	 *            array of floats
	 * @return the minimum
	 */
	public static float max(float[] array) {
		float max;

		if (array.length==0) {
			System.err.println("Error - impossible to find the max of an empty array");
			return -999;
		}
		max = array[0];
		for (int i = 1; i<array.length; i++) {
			if (array[i]!=Float.NaN&&array[i]>max) {
				max = array[i];
			}
		}

		return max;
	}

	/**
	 * Return the index of the maximum in an array of integers
	 * 
	 * @param array
	 *            array of integers
	 * @return the minimum
	 */
	public static int maxIndex(double[] array) {
		double max;
		int index;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return -1;
		}
		max = array[0];
		index = 0;
		for (int i = 1; i<array.length; i++) {
			if (array[i]>max) {
				max = array[i];
				index = i;
			}
		}
		return index;
	}

	/**
	 * Return the maximum in an array of numbers
	 * 
	 * @param array
	 *            array of numbers
	 * @return the maximum
	 */
	public static double max(double[] array) {
		double max;

		if (array.length==0) {
			System.err.println("Error - impossible to find the max of an empty array");
			return Double.NaN;
		}
		max = array[0];
		for (int i = 1; i<array.length; i++) {
			max = Math.max(array[i], max);
		}
		return max;
	}

	/**
	 * Creates an integer array of given size and initializes values to their
	 * respective indices
	 * 
	 * @param size
	 *            size of array
	 * @return array of integers with values initialized to their respective
	 *         indices
	 */
	public static int[] intArray(int size) {
		int[] arr = new int[size];
		for (int i = 0; i<size; i++) {
			arr[i] = i;
		}
		return arr;
	}

	/**
	 * Creates an integer array of given size and initializes each element with
	 * the given value
	 * 
	 * @param size
	 *            size of array
	 * @param initValue
	 *            initial value of each element
	 * @return array of integers initialized to the given value
	 */
	public static int[] intArray(int size, int initValue) {
		int[] arr = new int[size];
		for (int i = 0; i<size; i++) {
			arr[i] = initValue;
		}
		return arr;
	}

	/**
	 * Creates a long array of given size and initializes each element with
	 * the given value
	 * 
	 * @param size
	 *            size of array
	 * @param initValue
	 *            initial value of each element
	 * @return array of longs initialized to the given value
	 */
	public static long[] longArray(int size, long initValue) {
		long[] arr = new long[size];
		for (int i = 0; i<size; i++) {
			arr[i] = initValue;
		}
		return arr;
	}

	/**
	 * Creates an integer array from the contents of a string array
	 * 
	 * @param array
	 *            array of Strings to be converted
	 * @return array of the converted integers
	 */
	public static int[] toIntArray(String[] array) {
		int[] arr = new int[array.length];
		for (int i = 0; i<array.length; i++) {
			try {
				arr[i] = Integer.parseInt(array[i]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+array[i]+"' into an integer");
			}
		}
		return arr;
	}

	/**
	 * Creates an array of numbers from the contents of a string array
	 * 
	 * @param array
	 *            array of Strings to be converted
	 * @return array of the converted numbers
	 */
	public static double[] toDoubleArray(String[] array) {
		double[] arr = new double[array.length];
		for (int i = 0; i<array.length; i++) {
			try {
				arr[i] = Double.parseDouble(array[i]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+array[i]+"' into a number");
			}
		}
		return arr;
	}

	/**
	 * Creates an array of doubles from the contents of a float array
	 * 
	 * @param array
	 *            array of floats to be converted
	 * @return array of the converted doubles
	 */
	public static double[] toDoubleArray(float[] array) {
		double[] arr = new double[array.length];
		for (int i = 0; i<array.length; i++) {
			try {
				arr[i] = (double)array[i];
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+array[i]+"' into a double");
			}
		}
		return arr;
	}

	/**
	 * Creates an array of doubles from the contents of a byte array
	 * 
	 * @param array
	 *            array of floats to be converted
	 * @return array of the converted doubles
	 */
	public static double[] toDoubleArray(byte[] array) {
		double[] arr = new double[array.length];
		for (int i = 0; i<array.length; i++) {
			try {
				arr[i] = (double)array[i];
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+array[i]+"' into a double");
			}
		}
		return arr;
	}
	
	/**
	 * Creates an array of doubles and copies the contents of an int array into
	 * it
	 * 
	 * @param array
	 *            an array of integers
	 * @return an array of doubles copied from an array of integers
	 */
	public static double[] toDoubleArray(int[] array) {
		double[] arr = new double[array.length];
		for (int i = 0; i<array.length; i++) {
			arr[i] = array[i];
		}
		return arr;
	}

	/**
	 * Creates an array of doubles and copies the contents of a long array into it
	 * 
	 * @param array
	 *            an array of longs
	 * @return an array of doubles copied from an array of longs
	 */
	public static double[] toDoubleArray(long[] array) {
		double[] arr = new double[array.length];
		for (int i = 0; i<array.length; i++) {
			arr[i] = array[i];
		}
		return arr;
	}

	/**
	 * Creates a double array of given size and initializes each element with
	 * the given value
	 * 
	 * @param size
	 *            size of array
	 * @param initValue
	 *            initial value of each element
	 * @return array of numbers initialized to the given value
	 */
	public static double[] doubleArray(int size, double initValue) {
		double[] arr = new double[size];
		for (int i = 0; i<size; i++) {
			arr[i] = initValue;
		}
		return arr;
	}

	/**
	 * Creates a float array of given size and initializes each element with the
	 * given value
	 * 
	 * @param size
	 *            size of array
	 * @param initValue
	 *            initial value of each element
	 * @return array of floats initialized to the given value
	 */
	public static float[] floatArray(int size, float initValue) {
		float[] arr = new float[size];
		for (int i = 0; i<size; i++) {
			arr[i] = initValue;
		}
		return arr;
	}

	/**
	 * Creates a byte array of given size and initializes each element with the
	 * given value
	 * 
	 * @param size
	 *            size of array
	 * @param initValue
	 *            initial value of each element
	 * @return array of bytes initialized to the given value
	 */
	public static byte[] byteArray(int size, byte initValue) {
		byte[] arr = new byte[size];
		for (int i = 0; i<size; i++) {
			arr[i] = initValue;
		}
		return arr;
	}

	/**
	 * Creates a byte array from the contents of an integer array
	 * 
	 * @param array
	 *            array of ints to be converted
	 * @return array of the converted bytes
	 */
	public static byte[] toByteArray(int[] array) {
		byte[] arr = new byte[array.length];
		for (int i = 0; i<array.length; i++) {
			try {
				arr[i] = (byte)array[i];
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+array[i]+"' into a byte");
			}
		}
		return arr;
	}

	/**
	 * Creates a String array of given size and initializes each element with
	 * the given value
	 * 
	 * @param size
	 *            size of array
	 * @param initValue
	 *            initial value of each element
	 * @return array of Strings initialized to the given value
	 */
	public static String[] stringArray(int size, String initValue) {
		String[] array = new String[size];
		for (int i = 0; i<size; i++) {
			array[i] = initValue;
		}
		return array;
	}

	/**
	 * Creates a String array of given size and initializes each element with
	 * the "base value" +(index+1)
	 * 
	 * @param size
	 *            size of array
	 * @param base
	 *            base value of each element
	 * @return array of Strings initialized to the given indexed values
	 */
	public static String[] stringArraySequence(int size, String base) {
		return stringArraySequence(size, base, "");
	}

	/**
	 * Creates a String array of given size and initializes each element with
	 * the "prefix"+(index+1)+"suffix"
	 * 
	 * @param size
	 *            size of array
	 * @param prefix
	 *            prefix for each element
	 * @param suffix
	 *            suffix for each element
	 * @return array of Strings initialized to the given indexed values
	 */
	public static String[] stringArraySequence(int size, String prefix, String suffix) {
		String[] array = new String[size];
		for (int i = 0; i<size; i++) {
			array[i] = prefix+(i+1)+suffix;
		}
		return array;
	}

	/**
	 * Creates an array of blank Strings
	 * 
	 * @param size
	 *            size of array
	 * @return array of blank Strings
	 */
	public static String[] stringArray(int size) {
		return stringArray(size, "");
	}

	/**
	 * Creates an integer array of given size and initializes values by randomly
	 * sampling zero to size without replacement
	 * 
	 * @param size
	 *            size of array
	 * @return array of random indices for the given size
	 */
	public static int[] random(int size) {
		int[] array = new int[size];
		int index, num;

		for(int i=0; i<size; i++) {
			array[i] = i;
		}
		
		for(int i=size-1; i>=0; i--) {
			index = (int)(Math.random()*(i+1));
			num = array[i];
			array[i] = array[index];
			array[index] = num;
		}

		return array;
	}

	/**
	 * Creates an integer array of given size and initializes values by randomly
	 * sampling zero to size with replacement
	 * 
	 * @param size
	 *            size of array
	 * @return array of random indices for the given size
	 */
	public static int[] randomWithReplacement(int size, int numSelections) {
		int[] array = new int[numSelections];

		for(int i=0; i<numSelections; i++) {
			array[i] = (int)(Math.random()*(size+1));
		}

		return array;
	}
	
	/**
	 * Creates an integer array of given size and initializes values by randomly
	 * sampling zero to size without replacement
	 * 
	 * @param size
	 *            size of array
	 * @return array of random indices for the given size
	 */
	public static int[] random(int size, int numSelections) {
		int[] array = new int[size];
		int[] selections = new int[numSelections];
		int index;

		for(int i=0; i<size; i++) {
			array[i] = i;
		}
		
		for(int i=size-1; i>=0 && i>=size-numSelections; i--) {
			index = (int)(Math.random()*(i+1));
			selections[size-i-1] = array[index]; 
			array[index] = array[i];
		}

		return selections;
	}

	/**
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return sum of the array
	 */
	public static double sum(double[] array) {
		double sum = 0;

		for (int i = 0; i<array.length; i++) {
//			if ((array[i]+"").equals("NaN")) {
//				System.err.println("Are you sure you want to sum: "+array[i]);
//			}
			sum += array[i];
		}

		return sum;
	}

	/**
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return sum of the array
	 */
	public static double sumInOrder(double[] array, int[] order) {
		double sum = 0;

		for (int i = 0; i<array.length; i++) {
			sum += array[order[i]];
		}

		return sum;
	}

	/**
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return sum of the array
	 */
	public static double sumExactUsingBigDecimal(double[] array) {
		BigDecimal sum = BigDecimal.ZERO;

		for (int i = 0; i<array.length; i++) {
			sum = sum.add(BigDecimal.valueOf(array[i]));
		}

		return sum.doubleValue();
	}

	/**
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return sum of the array
	 */
	public static float sum(float[] array) {
		float sum = 0;

		for (int i = 0; i<array.length; i++) {
			sum += array[i];
		}

		return sum;
	}

	/**
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of integers
	 * @return sum of the array
	 */
	public static int sum(int[] array) {
		int sum = 0;

		for (int i = 0; i<array.length; i++) {
			sum += array[i];
		}

		return sum;
	}

	/**
	 * Calculates the mean of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static double mean(double[] array) {
		return sum(array)/array.length;
	}

	/**
	 * Calculates the mean of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static double mean(double[] array, boolean ignoreNaN) {
		double sum;
		int count;
		
		sum = 0;
		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (!Double.isNaN(array[i]) || !ignoreNaN) {
				sum += array[i];
				count++;
			}
		}
		
		if (count == 0) {
			return Double.NaN;
		}
		
		return sum/(double)count;
	}

	/**
	 * Calculates the mean of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static double mean(int[] array) {
		return (double)sum(array)/array.length;
	}

	/**
	 * Calculates the mean of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static double meanIf(double[] array, double[] filter, double filterValue, boolean ignoreNaN) {
		double sum;
		int count;
		
		if (array.length != filter.length) {
			System.err.println("Error - filter size does not match array size");
			return Double.NEGATIVE_INFINITY;
		}
		
		sum = 0;
		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (filter[i] == filterValue && (!ignoreNaN || !(array[i]+"").equals("NaN"))) {
				sum += array[i];
				count++;
			}
        }
		
		return sum/count;
	}

	/**
	 * Calculates the mean of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static float mean(float[] array) {
		return sum(array)/array.length;
	}

	/**
	 * Recreates array with only the first n values
	 * 
	 * @param array
	 *            an array of numbers
	 * @param n
	 *            the number of values to include
	 * @return an array of the first n values of the original
	 */
	public static double[] firstN(double[] array, int n) {
		double[] arr = new double[n];

		if (array.length<n) {
			System.err.println("Error - said to use the first "+n+" numbers, but the array is only size "+array.length);
			System.exit(1);
		}

		for (int i = 0; i<n; i++) {
			arr[i] = array[i];
		}

		return arr;
	}

	/**
	 * Calculates the variance of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return variance of the array
	 */
	public static double variance(double[] array) {
		double avg = mean(array);
		double sum = 0;

		for (int i = 0; i<array.length; i++) {
			sum += Math.pow(avg-array[i], 2);
		}

		return sum/(array.length-1);
	}

	/**
	 * Calculates the variance of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return variance of the array
	 */
	public static float variance(float[] array) {
		double sum, avg; // allows for larger arrays
		
		sum = 0;
		for (int i = 0; i<array.length; i++) {
			sum += array[i];
		}
		avg = sum/(double)array.length;

		for (int i = 0; i<array.length; i++) {
			sum += Math.pow(avg-array[i], 2);
		}

		return (float)(sum/(double)(array.length-1));
	}

	/**
	 * Calculates the standard deviation of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return standard deviation of the array
	 */
	public static double stdev(double[] array) {
		return Math.sqrt(variance(array));
	}

	/**
	 * Calculates the standard deviation of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @param removeNaN
	 *            remove any value that is not a number
	 * @return standard deviation of the [filtered] array
	 */
	public static float stdev(float[] array, boolean removeNaN) {
		double sum, avg;
		int count;
		
		sum = 0;
		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (!Float.isNaN(array[i]) || !removeNaN) {
				sum += array[i];
				count++;
			}
		}
		avg = (sum/(double)count);

		sum = 0;
		for (int i = 0; i<array.length; i++) {
			if (!Float.isNaN(array[i])) {
				sum += Math.pow(avg-array[i], 2);
			}
		}

		return (float)Math.sqrt(sum/(double)(count-1));
	}

	/**
	 * Calculates the standard deviation of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @param removeNaN
	 *            remove any value that is not a number
	 * @return standard deviation of the [filtered] array
	 */
	public static double stdev(double[] array, boolean removeNaN) {
		double sum, avg;
		int count;
		
		sum = 0;
		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (!Double.isNaN(array[i]) || !removeNaN) {
				sum += array[i];
				count++;
			}
		}
		avg = (sum/(double)count);

		sum = 0;
		for (int i = 0; i<array.length; i++) {
			if (!Double.isNaN(array[i])) {
				sum += Math.pow(avg-array[i], 2);
			}
		}

		return (float)Math.sqrt(sum/(double)(count-1));
	}

	/**
	 * Calculates the standard deviation of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return standard deviation of the array
	 */
	public static float stdevBigDecimal(float[] array, boolean removeNaN) {
		BigDecimal sum;
		int count;
		float avg;
		boolean first;
		
		if (true) {
			System.err.println("Error - the BigDecimal version of this method doesn't currently work");
			System.exit(1);
		}
		
		sum = BigDecimal.ZERO; // can't start with zero (see workaround below), or else BigDecimal won't add properly because it will assumes scale that the scale (i.e. number of decimal places) is zero
		count = 0;
		first = true;
		for (int i = 1; i<array.length; i++) {
			if (!Float.isNaN(array[i]) || !removeNaN) {
				if (first) {
					sum = BigDecimal.valueOf(array[i]);
				} else {
					sum.add(BigDecimal.valueOf(array[i]));
				}
				count++;
			}
		}
		avg = sum.divide(BigDecimal.valueOf(count), 8, RoundingMode.HALF_UP).floatValue();

		sum = BigDecimal.ZERO;
		for (int i = 0; i<array.length; i++) {
			if (!Float.isNaN(array[i])) {
				if (first) {
					sum = BigDecimal.valueOf(Math.pow(avg-array[i], 2));
				} else {
					sum.add(BigDecimal.valueOf(Math.pow(avg-array[i], 2)));
				}				
			}
		}

		return sum.divide(BigDecimal.valueOf(count-1), 8, RoundingMode.HALF_UP).floatValue();
	}
	
	/**
	 * Normalizes (calculates z-scores) for an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of z-scores
	 */
	public static double[] normalize(double[] array) {
		double[] newData = new double[array.length];
		double mean = Array.mean(array);
		double stdev = Array.stdev(array);

		for (int i = 0; i<newData.length; i++) {
			newData[i] = (array[i]-mean)/stdev;
		}

		return newData;
	}
	
	/**
	 * Inverse-normalizes an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of inverse-normalized values
	 */
	public static double[] inverseNormalize(double[] array) {
		double[] probits = new double[array.length];
		double rankDiv;
		int[] order;
		
		order = Sort.quicksort(array);

		for (int i = 0; i<probits.length; i++) {
			rankDiv = ((double)i+1)/((double)probits.length+1);
			if (rankDiv<0.5) {
				probits[order[i]] = ProbDist.NormDistReverse(rankDiv*2)*-1;
			} else {
				rankDiv = 1-rankDiv;
				probits[order[i]] = ProbDist.NormDistReverse(rankDiv*2)*1;
			}
		}

		return probits;
	}
	
	/**
	 * Returns the bootstrapped median of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @param numReps
	 *            number of replicates to perform
	 * @return array of the median, the 2.5 and the 97.5 bootstrap percentile
	 *         intervals
	 */
	public static double[] bootstrap(double[] array, int numReps, boolean verbose) {
		double[] results = new double[3];
		double[] replicates = new double[numReps];
		int[] keys;
		int progress = 0;
		double sum;

		if (verbose) {
			System.out.print("  Bootstrap rep 0");
		}
		for (int i = 0; i<numReps; i++) {
			if (verbose&&i+1==progress+numReps/10) {
				for (int j = 0; j<(progress+"").length(); j++) {
					System.out.print("\b");
				}
				progress += numReps/10;
				System.out.print(progress);
			}

			sum = 0;
			for (int j = 0; j<array.length; j++) {
				sum += array[(int)(Math.random()*array.length)];
			}
			replicates[i] = sum/array.length;
		}
		if (verbose) {
			System.out.println();
		}

		keys = Sort.quicksort(replicates);
		results[0] = replicates[keys[(int)((double)numReps*0.5)]];
		results[1] = replicates[keys[(int)((double)numReps*0.025)]];
		results[2] = replicates[keys[(int)((double)numReps*0.975)]];

		return results;
	}

	/**
	 * Determines the specified quantile of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @param q
	 *            quantile to be determined
	 * @return specified quantile of the array
	 */
	public static double quant(double[] array, double q) {
		int keys[] = Sort.quicksort(array);

		try {
			if (q>1||q<0) {
				return (0);
			} else {
				double index = (array.length+1)*q;
				if (index-(int)index==0) {
					return array[keys[(int)index-1]];
				} else {
					return q*array[keys[(int)Math.floor(index)-1]]+(1-q)*array[keys[(int)Math.ceil(index)-1]];
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return -1234567890;
		}
	}

	/**
	 * Determines the specified quantiles of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @param q
	 *            quantiles to be determined
	 * @return specified quantiles of the array
	 */
	public static float[] quants(float[] array, double[] qs) {
		int keys[] = Sort.quicksort(array);
		float[] quantiles;
		
		quantiles = new float[qs.length];
		for (int i = 0; i < quantiles.length; i++) {
			try {
				if (qs[i]>1||qs[i]<0) {
					quantiles[i] = -1;
				} else {
					double index = (array.length+1)*qs[i];
					if (index-(int)index==0) {
						quantiles[i] = array[keys[(int)index-1]];
					} else {
						quantiles[i] = (float)(qs[i]*array[keys[(int)Math.floor(index)-1]]+(1-qs[i])*array[keys[(int)Math.ceil(index)-1]]);
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
				quantiles[i] = -999;
			}
		}
		
		return quantiles;
	}

	/**
	 * Determines the median of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return median of the array
	 */
	public static double median(double[] array) {
		return (quant(array, 0.50));
	}

	/**
	 * Prints an array of objects separated by a tab
	 * 
	 * @param array
	 *            an array of objects
	 * @return String of printed objects
	 */
	public static String toStr(String[] array) {
		return toStr(array, null, "\t", null);
	}

	/**
	 * Prints an array of objects separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of objects
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed objects
	 */
	public static String toStr(String[] array, String delimiter) {
		return toStr(array, null, delimiter, null);
	}

	/**
	 * Prints an array of objects separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of objects
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed objects
	 */
	public static String toStr(String[] array, boolean[] display, String delimiter, String nullValue) {
		String str = "";
		int count;

		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (display == null || display[i]) {
				str += (count==0?"":delimiter)+(array[i]==null?nullValue:array[i]);
				count++;
			}
		}

		return str;
	}

	/**
	 * Prints the first element of each array in an array of arrays
	 * 
	 * @param array
	 *            an array of arrays
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed objects
	 */
	public static String toStr(String[][] array, String delimiter) {
		String str = "";

		for (int i = 0; i<array.length; i++) {
			str += (i==0?"":delimiter)+array[i][0];
		}

		return str;
	}

	/**
	 * Prints an array of integers separated by a tab
	 * 
	 * @param array
	 *            an array of integers
	 * @return String of printed integers
	 */
	public static String toStr(int[] array) {
		return toStr(array, "\t");
	}

	/**
	 * Prints an array of integers separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of integers
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed integers
	 */
	public static String toStr(int[] array, String delimiter) {
		String str = "";

		for (int i = 0; i<array.length; i++) {
			str += (i==0?"":delimiter)+array[i];
		}

		return str;
	}

	/**
	 * Prints an array of bytes separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of bytes
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed bytes
	 */
	public static String toStr(byte[] array, String delimiter) {
		String str = "";

		for (int i = 0; i<array.length; i++) {
			str += (i==0?"":delimiter)+array[i];
		}

		return str;
	}

	/**
	 * Prints an array of numbers with as many sigfigs as necessary, each
	 * separated by a tab
	 * 
	 * @param array
	 *            an array of numbers
	 * @return String of printed numbers
	 */
	public static String toStr(double[] array) {
		return toStr(array, -1, -1, "\t");
	}

	/**
	 * Prints an array of numbers separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of numbers
	 * @param sigfigs
	 *            number of significant digits
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed numbers
	 */
	public static String toStr(double[] array, int minSigFigs, int maxSigFigs, String delimiter) {
		String str = "";

		for (int i = 0; i<array.length; i++) {
			str += (i==0?"":delimiter)+(maxSigFigs==-1?ext.formDeci(array[i], 10):ext.formDeci(array[i], minSigFigs, maxSigFigs));
		}

		return str;
	}

	/**
	 * Prints an array of numbers with as many sigfigs as necessary, each
	 * separated by a tab
	 * 
	 * @param array
	 *            an array of numbers
	 * @return String of printed numbers
	 */
	public static String toStr(float[] array) {
		return toStr(array, -1, -1, "\t");
	}

	/**
	 * Prints an array of numbers separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of numbers
	 * @param sigfigs
	 *            number of significant digits
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed numbers
	 */
	public static String toStr(float[] array, int minSigFigs, int maxSigFigs, String delimiter) {
		String str = "";

		for (int i = 0; i<array.length; i++) {
			str += (i==0?"":delimiter)+(maxSigFigs==-1?ext.formDeci(array[i], 10):ext.formDeci(array[i], minSigFigs, maxSigFigs));
		}

		return str;
	}

	/**
	 * Returns an array splitting a number as equally as possible into different
	 * amounts
	 * 
	 * @param total
	 *            number to be split into groups
	 * @param numSplits
	 *            number of groups to split total into
	 * @return array of the numbers for each group
	 */
	public static int[] splitUp(int total, int numSplits) {
		int[] splits = new int[numSplits];

		for (int i = 0; i<numSplits-1; i++) {
			splits[i] = (int)Math.floor((double)total/(double)numSplits);
		}
		splits[numSplits-1] = total-(numSplits-1)*(int)Math.floor((double)total/(double)numSplits);

		return splits;
	}

	/**
	 * Creates an array of Strings and copies the contents of a Vector into it
	 * 
	 * @param v
	 *            vector of Strings
	 * @return an array of Strings from the Vector
	 */
	public static String[] toStringArray(Vector<String> v) {
		return v.toArray(new String[v.size()]);
	}

	/**
	 * Creates an array of Strings and copies the contents of a Vector into it in the specified order
	 * 
	 * @param v
	 *            vector of Strings
	 * @param order
	 *            order of elementes
	 * @return an array of ordered Strings from the Vector
	 */
	public static String[] toStringArray(Vector<String> v, int[] order) {
		String[] array;
		
		array = new String[v.size()];
		if (order.length != array.length) {
			System.err.println("Error - order does not have the same number of elements (n="+order.length+") as the Vector (n="+array.length+")");
			return null;
		}
		for (int i = 0; i<array.length; i++) {
			array[i] = v.elementAt(order[i]);
        }
		return array;
	}

	/**
	 * Creates an array of Strings and copies the contents of an ArrayList into
	 * it
	 * 
	 * @param v
	 *            vector of Strings
	 * @return an array of Strings from the ArrayList
	 */
	public static String[] toStringArray(ArrayList<String> al) {
		return al.toArray(new String[al.size()]);
	}

	/**
	 * Creates an array of Strings and copies the contents of an array of booleans
	 * 
	 * @param array
	 *            array of boolean
	 * @param onesAndZeros
	 *            array of boolean
	 * @return an array of the converted Strings
	 */
	public static String[] toStringArray(boolean[] array, boolean onesAndZeros) {
		String[] new_array;
		
		new_array = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			new_array[i] = onesAndZeros?(array[i]?"1":"0"):array[i]+"";
		}
		
		return new_array;
	}

	/**
	 * Creates an array of Strings and copies the contents of an array of longs
	 * 
	 * @param array
	 *            array of long
	 * @return an array of the converted Strings
	 */
	public static String[] toStringArray(long[] array) {
		String[] new_array;
		
		new_array = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			new_array[i] = array[i]+"";
		}
		
		return new_array;
	}

	/**
	 * Creates an array of Strings and copies the contents of an array of int
	 * 
	 * @param array
	 *            array of int
	 * @return an array of the converted Strings
	 */
	public static String[] toStringArray(int[] array) {
		String[] new_array;
		
		new_array = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			new_array[i] = array[i]+"";
		}
		
		return new_array;
	}

	/**
	 * Creates an array of Strings and copies the contents of an array of doubles
	 * 
	 * @param array
	 *            array of double
	 * @return an array of the converted Strings
	 */
	public static String[] toStringArray(double[] array) {
		String[] new_array;
		
		new_array = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			new_array[i] = array[i]+"";
		}
		
		return new_array;
	}
	
	/**
	 * Prints an array of Strings (culled from a Vector) separated by a tab
	 * 
	 * @param v
	 *            a Vector of Strings
	 * @return String of printed objects
	 */
	public static String toStr(Vector<String> v) {
		return toStr(toStringArray(v));
	}

	/**
	 * Creates a Vector and copies the conents of a String array into it
	 * 
	 * @param array
	 *            an array of Strings
	 * @return a vector of Strings
	 */
	public static Vector<String> toStringVector(String[] array) {
		Vector<String> v = new Vector<String>();
//		for (int i = 0; array != null && i<array.length; i++) {
		for (int i = 0; i<array.length; i++) {
			v.add(array[i]);
		}
		return v;
	}

	/**
	 * Returns an array of the unique Strings
	 * 
	 * @param array
	 *            an array of Strings
	 * @return array of the unique Strings
	 */
	public static String[] unique(String[] array) {
		Hashtable<String, String> hash = new Hashtable<String, String>();
		String[] newArray = new String[array.length];
		int count = 0;
		
		for (int i = 0; i<array.length; i++) {
			if (!hash.containsKey(array[i])) {
				newArray[count] = array[i];
				count++;
			}
		}

		return Array.subArray(newArray, 0, count);
	}

	/**
	 * Returns an array of the unique Strings
	 * 
	 * @param array
	 *            an array of Strings
	 * @return array of the unique Strings
	 */
	public static String[] uniqueOldSlowForLargeArrays(String[] array) {
		Vector<String> v = new Vector<String>(array.length);

		for (int i = 0; i<array.length; i++) {
			HashVec.addIfAbsent(array[i], v);
		}

		return Array.toStringArray(v);
	}

//	/**
//	 * Returns an array minus any entries with a null or value listed as invalid
//	 * 
//	 * @param array
//	 *            an array of Strings
//	 * @param thoseToBeRemoved
//	 *            an array of Strings representing invalid values
//	 * @return array of valid Strings
//	 */
//	public static String[] removeThese(String[] array, String[] thoseToBeRemoved) {
//		Vector<String> v = new Vector<String>();
//
//		for (int i = 0; i<array.length; i++) {
//			if (array[i]!=null&&ext.indexOfStr(array[i], thoseToBeRemoved)==-1) {
//				v.add(array[i]);
//			}
//		}
//
//		return Array.toStringArray(v);
//	}

	/**
	 * Inserts the specified String into an array the specified position
	 * 
	 * @param str
	 *            String to be inserted
	 * @param array
	 *            an array of Strings
	 * @param pos
	 *            position for str to be inserted
	 * @return altered array of Strings
	 */
	public static String[] insertStringAt(String str, String[] array, int pos) {
		Vector<String> v;

		if (pos<0||pos>array.length) {
			throw new ArrayIndexOutOfBoundsException(pos);
		}

		v = toStringVector(array);
		v.insertElementAt(str, pos);

		return Array.toStringArray(v);
	}

	/**
	 * Populate an array of doubles by converting the indicated members of a
	 * String array
	 * 
	 * @param line
	 *            array of Strings
	 * @param indices
	 *            indices of the Strings to be converted to doubles
	 * @return the resulting double array
	 */
	public static double[] populateDoubles(String[] line, int[] indices) {
		double[] array = new double[indices.length];

		for (int i = 0; i<indices.length; i++) {
			array[i] = line[indices[i]].equals(".")?-999:Double.parseDouble(line[indices[i]]);
		}

		return array;
	}

	/**
	 * Create a boolean array and make all states TRUE
	 * 
	 * @param int
	 *            size of array
	 * @return the boolean array
	 */
	public static boolean[] booleanArray(int size, boolean initValue) {
		boolean[] array = new boolean[size];

		for (int i = 0; i<size; i++) {
			array[i] = initValue;
		}

		return array;
	}

//	/**
//	 * Create an array of initialized IntVectors
//	 * 
//	 * @param int
//	 *            size of array
//	 * @return the array of IntVectors
//	 */
//	public static IntVector[] intVectorArray(int size) {
//		IntVector[] array = new IntVector[size];
//
//		for (int i = 0; i<size; i++) {
//			array[i] = new IntVector();
//		}
//
//		return array;
//	}
//
//	/**
//	 * Create an array of initialized DoubleVectors
//	 * 
//	 * @param int
//	 *            size of array
//	 * @return the array of Vectors
//	 */
//	public static DoubleVector[] doubleVectorArray(int size) {
//		DoubleVector[] array = new DoubleVector[size];
//
//		for (int i = 0; i<size; i++) {
//			array[i] = new DoubleVector();
//		}
//
//		return array;
//	}

	/**
	 * Find all instances of the values found in the first column of the matrix
	 * and replace them with the values from the second column
	 * 
	 * @param array
	 *            an array of integers
	 * @param swaps
	 *            a matrix of old and new values
	 * @return the same array but with the swapped values
	 */
	public static int[] findReplace(int[] array, int[][] swaps) {
		int[] newArray = array.clone();

		for (int i = 0; i<array.length; i++) {
			for (int j = 0; j<swaps.length; j++) {
				if (array[i]==swaps[j][0]) {
					newArray[i] = swaps[j][1];
				}
			}
		}

		return newArray;
	}

	/**
	 * Tries to finds the first instance of a given number within an array and
	 * returns either the index or -1 if not found
	 * 
	 * @param array
	 *            an array of integers
	 * @param target
	 *            the number to find
	 * @return the index or -1 if not found
	 */
	public static int indexOfInt(int[] array, int target) {
		for (int i = 0; i<array.length; i++) {
			if (array[i]==target) {
				return i;
			}
		}

		return -1;
	}

	/**
	 * Tries to finds the first instance of a given number within an array and
	 * returns either the index or -1 if not found
	 * 
	 * @param array
	 *            an array of bytes
	 * @param target
	 *            the number to find
	 * @return the index or -1 if not found
	 */
	public static int indexOfByte(byte[] array, byte target) {
		for (int i = 0; i<array.length; i++) {
			if (array[i]==target) {
				return i;
			}
		}

		return -1;
	}
	
	/**
	 * Creates a new array using only the indices between start and stop
	 * 
	 * @param array
	 *            an array of doubles
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static int[] subArray(int[] array, int start, int stopBefore) {
		int[] arr;

		if (start<0||stopBefore>array.length||stopBefore<=start) {
			System.err.println("Error - invalid start ("+start+") and stopBefore ("+stopBefore+") indicies for an array");
		}
		arr = new int[stopBefore-start];
		for (int i = start; i<stopBefore; i++) {
			arr[i-start] = array[i];
		}

		return arr;
	}

	/**
	 * Creates a new array using only the indices at and after start
	 * 
	 * @param array
	 *            an array of doubles
	 * @param start
	 *            first index to use
	 * @return the subset of the original array
	 */
	public static int[] subArray(int[] array, int start) {
		return subArray(array, start, array.length);
	}

	/**
	 * Creates a new array using only the indices between start and stop
	 * 
	 * @param array
	 *            an array of doubles
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static double[] subArray(double[] array, int start, int stopBefore) {
		double[] arr;

		if (start<0||stopBefore>array.length||stopBefore<=start) {
			System.err.println("Error - invalid start ("+start+") and stopBefore ("+stopBefore+") indicies for an array");
		}
		arr = new double[stopBefore-start];
		for (int i = start; i<stopBefore; i++) {
			arr[i-start] = array[i];
		}

		return arr;
	}

	/**
	 * Creates a new array using only the indices at and after start
	 * 
	 * @param array
	 *            an array of doubles
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static double[] subArray(double[] array, int start) {
		return subArray(array, start);
	}

	/**
	 * Creates a new array using only the indices between start and stop
	 * 
	 * @param array
	 *            an array of floats
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static float[] subArray(float[] array, int start, int stopBefore) {
		float[] arr;

		if (start<0||stopBefore>array.length||stopBefore<=start) {
			System.err.println("Error - invalid start ("+start+") and stopBefore ("+stopBefore+") indicies for an array");
		}
		arr = new float[stopBefore-start];
		for (int i = start; i<stopBefore; i++) {
			arr[i-start] = array[i];
		}

		return arr;
	}

	/**
	 * Creates a new array using only the indices between start and stop
	 * 
	 * @param array
	 *            an array of Strings
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static String[] subArray(String[] array, int start, int stopBefore) {
		String[] arr;

		if (start<0||stopBefore>array.length||stopBefore<start) {
			System.err.println("Error - invalid start ("+start+") and stopBefore ("+stopBefore+") indicies for an array of size "+array.length);
		}
		arr = new String[stopBefore-start];
		for (int i = start; i<stopBefore; i++) {
			arr[i-start] = array[i];
		}

		return arr;
	}

	/**
	 * Creates a new array using only the indices at and after start
	 * 
	 * @param array
	 *            an array of Strings
	 * @param start
	 *            first index to use
	 * @return the subset of the original array
	 */
	public static String[] subArray(String[] array, int start) {
		return subArray(array, start, array.length);
	}

	/**
	 * Creates a new array using only the indices provided
	 * 
	 * @param array
	 *            an array of Strings
	 * @param indices
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static String[] subArray(String[] array, int[] indices) {
		return subArray(array, replaceAllWith(indices, -1, -2), null);
	}

	/**
	 * Creates a new array using only the indices provided
	 * 
	 * @param array
	 *            an array of Strings
	 * @param indices
	 *            indices to use
	 * @param missingValue
	 *            if an index equals -1, then missingValue will be inserted
	 * @return the subset of the original array
	 */
	public static String[] subArray(String[] array, int[] indices, String missingValue) {
		String[] strs = new String[indices.length];

		if (indices.length == 0) {
			return strs;
		}
		if (min(indices)<-1) {
			System.err.println("Error missing an index; subarray will fail");
			System.exit(1);
		}
		for (int i = 0; i<strs.length; i++) {
			if (indices[i] == -1) {
				strs[i] = missingValue;
			} else {
				strs[i] = array[indices[i]];
			}
		}

		return strs;
	}

	/**
	 * Increments a boolean array
	 * 
	 * @param array
	 *            an array of boolean
	 */
	public static void incBoolean(boolean[] array) {
		for (int i = array.length-1; i>=0; i--) {
			if (!array[i]) {
				array[i] = true;
				for (int j = i+1; j<array.length; j++) {
					array[j] = false;
				}
				return;
			}
		}

	}

	/**
	 * Provides a frequency distribution of an array of Strings
	 * 
	 * @param array
	 *            an array of Strings
	 * @return the unique elements and their counts
	 */
	public static String[][] frequency(String[] array) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int count;
		String[] keys;
		String[][] summary;

		for (int i = 0; i<array.length; i++) {
			if (hash.containsKey(array[i])) {
				count = Integer.parseInt(hash.get(array[i]));
			} else {
				count = 0;
			}
			count++;
			hash.put(array[i], count+"");
		}

		keys = Sort.putInOrder(HashVec.getKeys(hash));
		summary = new String[keys.length][2];
		for (int i = 0; i<summary.length; i++) {
			summary[i][0] = keys[i];
			summary[i][1] = hash.get(keys[i]);
		}

		return summary;
	}

	public static int booleanArraySum(boolean[] array) {
		int count = 0;

		for (int i = 0; i<array.length; i++) {
			if (array[i]) {
				count++;
			}
		}

		return count;
	}

	public static int binarySearch(int[] array, int value, boolean exact) {
		return binarySearch(array, value, 0, array.length-1, exact);
	}

	public static int binarySearch(int[] array, int value, int low, int high, boolean exact) {
		int mid;

		while (low<=high) {
//			System.out.println(array[low]+"\t"+array[high]);
			mid = low+(high-low)/2;
			if (array[mid]>value) {
				high = mid-1;
			} else if (array[mid]<value) {
				low = mid+1;
			} else {
				return mid;
			}
		}
		if (exact) {
			return -1;
		} else {
			return low;
		}
	}

	public static String[] booleanArrayToStringArray(boolean[] array) {
		String[] strs = new String[array.length];

		for (int i = 0; i<array.length; i++) {
			strs[i] = array[i]?"1":"0";
		}

		return strs;
	}

	public static int countIf(String[] array, String target) {
		int count = 0;

		for (int i = 0; i<array.length; i++) {
			if (array[i].equals(target)) {
				count++;
			}
		}

		return count;
	}

	public static int countIf(int[] array, int target) {
		int count = 0;

		for (int i = 0; i<array.length; i++) {
			if (array[i]==target) {
				count++;
			}
		}

		return count;
	}

	public static boolean[] booleanNegative(boolean[] array) {
		boolean[] newArray = new boolean[array.length];

		for (int i = 0; i<array.length; i++) {
			newArray[i] = !array[i];
		}

		return newArray;
	}

	/**
	 * Calculates the interquartile range of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return iqr of the array
	 */
	public static double iqr(double[] array) {
		int[] keys = Sort.quicksort(array);
		double iqr = 0;

		if (array.length<2) {
			System.err.println("Error - can't calculate an IQR for an array with "+array.length+" datapoint(s)");
			return -1;
		}
		try {
			iqr = array[keys[(int)Math.ceil(array.length*0.75)]]-array[keys[(int)Math.ceil(array.length*0.25)]];
		} catch (Exception e) {
			System.err.println("Error calculating IQR");
			e.printStackTrace();
		}

		return iqr;
	}

	/**
	 * Calculates the inter quartile range of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return iqr of the array
	 */
	public static float iqr(float[] array) {
		int[] keys = Sort.quicksort(array);
		float iqr = 0;

		if (array.length<2) {
			System.err.println("Error - can't calculate an IQR for an array with "+array.length+" datapoint(s)");
			return -1;
		}
		try {
			iqr = array[keys[(int)Math.ceil(array.length*0.75)]]-array[keys[(int)Math.ceil(array.length*0.25)]];
		} catch (Exception e) {
			System.err.println("Error calculating IQR");
			e.printStackTrace();
		}

		return iqr;
	}
	
	/**
	 * Trims null values from the end of a String array
	 * 
	 * @param array
	 *            an array of Strings
	 * @return trimmed array
	 */
	public static String[] trimArray(String[] array) {
		int index;

		index = array.length;
		while (index > 0 && array[index-1] == null) {
			index--;
		}

		return index == array.length?array:subArray(array, 0, index);
	}

	/**
	 * Returns true if the String arrays are equal at all positions
	 * 
	 * @param array1
	 *            an array of Strings
	 * @param array2
	 *            an array of Strings
	 * @param caseSensitive
	 *            boolean flag
	 * @return true if arrays are equal
	 */
	public static boolean equals(String[] array1, String[] array2, boolean caseSenstitive) {
		if (array1.length != array2.length) {
			return false;
		}

		for (int i = 0; i<array1.length; i++) {
			if (caseSenstitive) {
				if (!array1[i].equals(array2[i])) {
					return false;
				}
			} else if (!array1[i].equalsIgnoreCase(array2[i])) {
				return false;
			}
        }
		
		return true;
	}
	
	/**
	 * Returns true if the String arrays are equal at all positions
	 * 
	 * @param array1
	 *            an array of Strings
	 * @param array2
	 *            an array of Strings
	 * @param caseSensitive
	 *            boolean flag
	 * @return true if arrays are equal
	 */
	public static boolean equals(int[] array1, int[] array2) {
		if (array1.length != array2.length) {
			return false;
		}

		for (int i = 0; i<array1.length; i++) {
			if (array1[i] != array2[i]) {
				return false;
			}
        }
		
		return true;
	}
	
	/**
	 * Removes NaN values from the array
	 * 
	 * @param array
	 *            an array of doubles
	 * @return scrubbed array
	 */
	public static double[] removeNaN(double[] array) {
		double[] newArray;
		boolean[] use;
		int count;
		
		use = new boolean[array.length];
		for (int i = 0; i<use.length; i++) {
			use[i] = !(array[i]+"").equals("NaN");
        }

		newArray = new double[booleanArraySum(use)];
		count = 0;
		for (int i = 0; i<use.length; i++) {
			if (use[i]) {
				newArray[count] = array[i]; 
				count++;
			}
        }
		
		return newArray;
	}
	
	/**
	 * Copies an array exactly
	 * 
	 * @param array
	 *            an array of integers
	 * @return copy of array
	 */
	public static int[] copyArray(int[] array) {
		int[] newArray;


		newArray = new int[array.length];
		for (int i = 0; i<array.length; i++) {
			newArray[i] = array[i]; 
        }
		
		return newArray;
	}
	
	/**
	 * Removes certain values from the array
	 * 
	 * @param array
	 *            an array of Strings
	 * @param thingsToRemove
	 *            list of Strings to removed from the array 
	 * @return scrubbed array
	 */
	public static String[] removeFromArray(String[] array, String[] thingsToRemove) {
		String[] newArray;
		boolean[] use;
		int count;
		
		use = new boolean[array.length];
		for (int i = 0; i<use.length; i++) {
			use[i] = ext.indexOfStr(array[i], thingsToRemove) == -1;
        }

		newArray = new String[booleanArraySum(use)];
		count = 0;
		for (int i = 0; i<use.length; i++) {
			if (use[i]) {
				newArray[count] = array[i]; 
				count++;
			}
        }
		
		return newArray;
	}

	/**
	 * Removes String at given index of the array
	 * 
	 * @param array
	 *            an array of Strings
	 * @param index
	 *            index of String to remove from the array 
	 * @return new array
	 */
	public static String[] removeFromArray(String[] array, int index) {
		String[] newArray;
		
		newArray = new String[array.length-1];
		for (int i = 0; i<index; i++) {
			newArray[i] = array[i]; 
        }
		for (int i = index+1; i<array.length; i++) {
			newArray[i-1] = array[i]; 
        }
		
		return newArray;
	}

	/**
	 * Adds specified String to the end of an array
	 * 
	 * @param array
	 *            an array of Strings
	 * @param str
	 *            String to add to the array 
	 * @return new array
	 */
	public static String[] addStrToArray(String str, String[] array) {
		return addStrToArray(str, array, array.length);
    }

	/**
	 * Adds specified String to a specified index of an array
	 * 
	 * @param array
	 *            an array of Strings
	 * @param str
	 *            String to add to the array 
	 * @param indexOfNewArray
	 *            location of the string in the new array 
	 * @return new array
	 */
	public static String[] addStrToArray(String str, String[] array, int indexOfNewStr) {
    	String[] new_array = new String[array.length+1];
    
//    	for (int i = 0; i<array.length; i++) {
    	for (int i = 0; i<new_array.length; i++) {
//    		new_array[i<=indexOfNewStr?i:i+1] = i==indexOfNewStr?str:array[i];
    		new_array[i] = i==indexOfNewStr?str:array[i>indexOfNewStr?i-1:i];
    	}
    
    	return new_array;
    }

	/**
	 * Adds specified integer to the end of an array
	 * 
	 * @param array
	 *            an array of integers
	 * @param value
	 *            integer to add to the array 
	 * @return new array
	 */
	public static int[] addIntToArray(int value, int[] array) {
    	int[] new_array = new int[array.length+1];
    
    	for (int i = 0; i<array.length; i++) {
    		new_array[i] = array[i];
    	}
    	new_array[array.length] = value;
    
    	return new_array;
    }
	
	/**
	 * Determine if the trait contained within the file is dichotomous or continuous
	 * 
	 * @param filename
	 *            filename containing the trait to be evaluated
	 * @param col
	 *            column to extract 
	 * @param allow21
	 *            allow binary trait to be 2 and 1 instead of 1 and 0 
	 * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropriate for neither 
	 */
	public static int determineType(String filename, int col, boolean allow21) {
        return determineType(filename, col, ext.MISSING_VALUES, allow21);
	}

	/**
	 * Determine if the trait contained within the file is dichotomous or continuous
	 * 
	 * @param filename
	 *            filename containing the trait to be evaluated
	 * @param col
	 *            column to extract 
	 * @param exclusions
	 *            values to exclude 
	 * @param allow21
	 *            allow binary trait to be 2 and 1 instead of 1 and 0 
	 * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropriate for neither 
	 */
	public static int determineType(String filename, int col, String[] exclusions, boolean allow21) {
        return determineType(allow21, Array.toDoubleArray(Array.removeFromArray(HashVec.loadFileToStringArray(filename, true, new int[] {col}, true), exclusions)));
	}

	/**
	 * Determine if the trait contained within the file is dichotomous or continuous
	 * 
	 * @param array
	 *            an array of doubles
	 * @param allow21
	 *            allow binary trait to be 2 and 1 instead of 1 and 0 
	 * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropraite for neither 
	 */
	public static int determineType(double[] array, boolean allow21) {
		Hashtable<String,String> hash;
		
		hash = new Hashtable<String,String>();
		for (int i = 0; hash.size() < 3 && i < array.length; i++) {
			if (!Double.isNaN(array[i])) {
				hash.put(ext.formDeci(array[i], 2), "");
			}
		}
		
		return determineType(allow21, Array.toDoubleArray(HashVec.getKeys(hash)));
	}
	
	/**
	 * Determine if the trait contained within the file is dichotomous or continuous
	 * 
	 * @param array
	 *            an array of doubles
	 * @param allow21
	 *            allow binary trait to be 2 and 1 instead of 1 and 0 
	 * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropraite for neither 
	 */
	private static int determineType(boolean allow21, double[] array) {
        if (array.length == 0) {
        	System.err.println("Error - no valid data for this trait!");
        	return -1;
        } else if (array.length == 1) {
        	System.err.println("Error - no variation for this trait!");
        	return -1;
        } else if (array.length == 2 && Array.min(array) == 0 && Array.max(array) == 1) {
        	return 0;
        } else if (allow21 && array.length == 2 && Array.min(array) == 1 && Array.max(array) == 2) {
        	return 0;
        } else if (array.length == 2) {
        	System.err.println("Error - flag was set to prevent binary trait from being anything other than 0/1 ("+array[0]+"/"+array[1]+" is not valid)");
        	return -1;
        } else if (array.length > 2) {
        	return 1;
        } else {
        	System.err.println("Unexpected error parsing phenotype class");
        	return -1;
        }
	}
	
	/**
	 * Splits an array of Strings into an array of arrays of Strings
	 * 
	 * @param array
	 *            an array of Strings
	 * @param tab
	 *            split using tabs as opposed to white spaces 
	 * @return an array of arrays of Strings 
	 */
	public static String[][] splitStrings(String[] array, boolean tab) {
		String[][] stringArrays;
		
		stringArrays = new String[array.length][];
		for (int i = 0; i < array.length; i++) {
			if (tab) {
				stringArrays[i] = array[i].split("\t", -1);
			} else {
				stringArrays[i] = array[i].split("[\\s]+");
			}
		}
		
		return stringArrays;
	}
	
	/**
	 * Transposes array to a matrix with a Splits an array of Strings into an array of arrays of Strings
	 * 
	 * @param array
	 *            an array of Strings
	 * @param tab
	 *            split using tabs as opposed to white spaces 
	 * @return an array of arrays of Strings 
	 */
	public static String[][] toMatrix(String[] array) {
		String[][] matrix;
		
		matrix = new String[array.length][1];
		for (int i = 0; i < array.length; i++) {
			matrix[i][0] = array[i];
		}
		
		return matrix;
	}
	
	/**
	 * Replace all of one value of integers with another
	 * 
	 * @param array
	 *            an array of integers
	 * @param from
	 *            value to replace 
	 * @param to
	 *            value to replace with 
	 * @return new array of integers
	 */
	public static int[] replaceAllWith(int[] array, int from, int to) {
		int[] newArray;
		
		newArray = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			if (array[i] == from) {
				newArray[i] = to;
			} else {
				newArray[i] = array[i];
			}
		}
		
		return newArray;
	}	
	
	/**
	 * Merge contents of two String arrays
	 * 
	 * @param array1
	 *            an array of Strings
	 * @param array2
	 *            an array of Strings
	 * @param numberOfMismatchesAllowed
	 *            maximum number of mismatched values allowed before the null set is returned
	 * @return merged array of Strings if successful, otherwise null
	 */
	public static String[] merge(String[] array1, String[] array2, int numberOfMismatchesAllowed) {
		String[] newArray;
		int count;
		
		if (array1.length != array2.length) {
			System.err.println("Error - mismatched number of values in the two arrays to be merged");
			return null;
		}
		
		count = 0;
		newArray = new String[array1.length];
		for (int i = 0; i < array1.length; i++) {
			if (array1[i] == null) {
				newArray[i] = array2[i];
			} else if (array2[i] == null) {
				newArray[i] = array1[i];
			} else if (array1[i].equals(array2[i])) {
				newArray[i] = array1[i];
			} else {
				newArray[i] = array1[i]+"/"+array2[i];
				count++;
			}
		}
		
		if (count > numberOfMismatchesAllowed) {
			return null;
		} else {
			return newArray;
		}
	}	

	/**
	 * Clones an array of Strings
	 * 
	 * @param array
	 *            the array of Strings to clone
	 * @return cloned array of Strings
	 */
	public static String[] clone(String[] array) {
		String[] newArray;

		newArray = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[i];
		}
		
		return newArray;
	}	

	/**
	 * Returns true if all Strings in the array represent numbers
	 * 
	 * @param array
	 *            an array of Strings
	 * @return merged array of Strings if successful, otherwise null
	 */
	public static boolean isAllNumbers(String[] array) {
		for (int i = 0; i < array.length; i++) {
			try {
				Double.parseDouble(array[i]);
			} catch (NumberFormatException nfe) {
				return false;
			}
		}

		return true;
	}	

	/**
	 * Returns the length of the largest String in the array
	 * 
	 * @param array
	 *            an array of Strings
	 * @return largest length
	 */
	public static int maxLength(String[] array) {
		int max;
		
		max = -1;
		for (int i = 0; i < array.length; i++) {
			if (array[i] != null) {
				max = Math.max(max, array[i].length());
			}
		}

		return max;
	}


	/**
	 * Converts frequency counts into proportions. So, the input looks similar to this:
	 * 				FrequencyCount
	 * 			   ---------------
	 * 		Female		152
	 * 		Male		148
	 * 
	 * The output looks like this:
	 * 				FrequencyCount
	 * 			   ---------------
	 * 		Female		50.67%
	 * 		Male		49.33%
	 * 
	 * @param counts the frequency counts in array format
	 * @return the corresponding proportion in array format
	 */
	public static double[] getProportions(int[] counts) {
		int total=0;
		double result[] = new double[counts.length];
		for (int i=0; i<counts.length; i++) {
			total += counts[i];
		}
		for (int i=0; i<counts.length; i++) {
			result[i] = (double) counts[i]/(double) total;
		}
		return result;
	}
	
	public static void main(String[] args) {
	    int[] counts = new int[10];
	    int[] trav;
	    
	    for (int i = 0; i<100000; i++) {
	    	trav = random(10, 5);
	    	for (int j = 0; j<trav.length; j++) {
	    		counts[trav[j]]++;
            }
        }
	    
	    Files.writeList(Array.toStr(counts).split("[\\s]+"), "oi.xln");

    }
}
