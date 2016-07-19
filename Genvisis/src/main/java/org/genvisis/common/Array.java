package org.genvisis.common;

import java.io.UnsupportedEncodingException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

import org.genvisis.stats.Maths;
import org.genvisis.stats.ProbDist;

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
	 * Return the minimum and maximum, respectively, in an array of floats
	 * 
	 * @param array
	 *            array of floats
	 * @return the minimum and maximum, in that order
	 */
	public static float[] minMax(float[] array) {
	    float min;
	    float max;
	    
	    if (array.length==0) {
	        System.err.println("Error - impossible to find the min of an empty array");
	        return new float[]{Float.NaN, Float.NaN};
	    }
	    min = max = array[0];
	    for (int i = 1; i<array.length; i++) {
	        if (!Float.isNaN(array[i])) {
	            if (array[i] < min) {
	                min = array[i];
	            } 
	            if (array[i] > max) {
	                max = array[i];
	            }
	        }
	    }
	    return new float[]{min, max};
	}

	/**
	 * Return the minimum and maximum, respectively, in an array of doubles
	 * 
	 * @param array
	 *            array of doubles
	 * @return the minimum and maximum, in that order
	 */
	public static double[] minMax(double[] array) {
	    double min;
	    double max;
	    
	    if (array.length==0) {
	        System.err.println("Error - impossible to find the min of an empty array");
	        return new double[]{Double.NaN, Double.NaN};
	    }
	    min = max = array[0];
	    for (int i = 1; i<array.length; i++) {
	        if (!Double.isNaN(array[i])) {
	            if (array[i] < min) {
	                min = array[i];
	            } 
	            if (array[i] > max) {
	                max = array[i];
	            }
	        }
	    }
	    return new double[]{min, max};
	}

	/**
	 * Return the minimum in an array of bytes
	 * 
	 * @param array
	 *            array of integers
	 * @return the minimum
	 */
	public static byte min(byte[] array) {
		byte min;

		if (array.length==0) {
			System.err.println("Error - impossible to find the min of an empty array");
			return Byte.MAX_VALUE;
		}
		min = array[0];
		for (int i = 1; i<array.length; i++) {
			min = (byte) Math.min(array[i], min);
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
	 * Return the maximum in an array of bytes
	 * 
	 * @param array
	 *            array of numbers
	 * @return the maximum
	 */
	public static byte max(byte[] array) {
		byte max;

		if (array.length==0) {
			System.err.println("Error - impossible to find the max of an empty array");
			return Byte.MAX_VALUE;
		}
		max = array[0];
		for (int i = 1; i<array.length; i++) {
			max = (byte) Math.max(array[i], max);
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
	 * Creates an integer array of given size and initializes each element 
	 * to 1 except for the first, which is set to zero
	 * 
	 * @param size
	 *            size of array
	 * @return array of integers initialized to the correct values
	 */
	public static int[] intArrayStandarddSkips(int size) {
		int[] arr = new int[size];
		for (int i = 0; i<size; i++) {
			arr[i] = i==0?0:1;
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
	 * Creates an integer array from the contents of a byte array
	 * 
	 * @param array
	 *            array of Strings to be converted
	 * @return array of the converted integers
	 */
	public static int[] toIntArray(byte[] array) {
		int[] arr = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			try {
				arr[i] = (array[i]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '" + array[i] + "' into an integer");
			}
		}
		return arr;
	}

	/**
	 * Creates an integer array from the contents of a double array
	 * 
	 * @param array
	 *            array of double to be converted
	 * @return array of the converted integers
	 */
	public static int[] toIntArray(double[] array) {
		int[] arr = new int[array.length];
		for (int i = 0; i<array.length; i++) {
			arr[i] = (int)array[i];
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
	public static double[][] toDoubleArrays(String[][] array, boolean NaNForMissing) {
		double[][] arr = new double[array.length][];
		for (int i = 0; i < array.length; i++) {
			arr[i] = Array.toDoubleArray(array[i], NaNForMissing);
		}
		return arr;
	}

	/**
	 * Creates an array of numbers from the contents of a string array
	 * 
	 * @param array
	 *            array of Strings to be converted
	 * 
	 * @return array of the converted numbers
	 */
	public static double[] toDoubleArray(String[] array) {
		return toDoubleArray(array, false);
	}

	public static double[] toDoubleArray(String[] array, boolean NaNForMissing) {
		double[] arr = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			try {
				arr[i] = Double.parseDouble(array[i]);
			} catch (NumberFormatException nfe) {
				if (NaNForMissing) {
					arr[i] = Double.NaN;
				} else {
					System.err.println("Error - failed to convert '" + array[i] + "' into a number");
					return null;
				}
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
	 * @return array of the converted numbers
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
	 * Creates an array of floats from the contents of a double array
	 * 
	 * @param array
	 *            array of doubles to be converted
	 * @return array of the converted numbers
	 */
	public static float[] toFloatArray(double[] array) {
		float[] arr = new float[array.length];
		for (int i = 0; i<array.length; i++) {
			try {
				arr[i] = (float)array[i];
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+array[i]+"' into a double");
			}
		}
		return arr;
	}

	/**
	 * Creates an array of floats from the contents of a byte array
	 * 
	 * @param array
	 *            array of bytes to be converted
	 * @return array of the converted numbers
	 */
	public static float[] toFloatArray(byte[] array) {
	    float[] arr = new float[array.length];
	    for (int i = 0; i<array.length; i++) {
	        try {
	            arr[i] = (float)array[i];
	        } catch (NumberFormatException nfe) {
	            System.err.println("Error - failed to convert '"+array[i]+"' into a double");
	        }
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
	
	public enum BYTE_DECODE_FORMAT{
		/**
		 * String will be converted to upper case
		 */
		UPPER_CASE, /**
		 * String will be converted to lower case
		 */
		LOWER_CASE, /**
		 * String will be left as is
		 */
		AS_IS
	}

	public static String[] decodeByteArray(byte[] b, BYTE_DECODE_FORMAT format, Logger log) {
		return decodeByteArray(b, "UTF-8", format, log);
	}
	
	public static String[] decodeByteArray(byte[] b, Logger log) {
		return decodeByteArray(b, "UTF-8", BYTE_DECODE_FORMAT.AS_IS, log);
	}

	/**
	 * @param b
	 *            each entry will be converted to a string
	 * @param charsetName
	 * @param format
	 * @param log
	 * @return
	 */
	public static String[] decodeByteArray(byte[] b, String charsetName, BYTE_DECODE_FORMAT format, Logger log) {
		String[] s = new String[b.length];
		for (int i = 0; i < s.length; i++) {
			if ((i + 1) % 2000000 == 0) {
				log.reportTimeInfo((i + 1) + " entries converted");
			}
			try {
				s[i] = new String(new byte[] { b[i] }, charsetName).toUpperCase();
			} catch (UnsupportedEncodingException e) {
				log.reportTimeError("Could not convert reference byte " + b[i] + " to string with charsetName" + charsetName);
				e.printStackTrace();
			}
		}
		return s;
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
	 * Creates a byte array from the contents of an integer array
	 * 
	 * @param array
	 *            array of floats to be converted
	 * @return array of the converted bytes
	 */
	public static byte[] toByteArray(float[] array) {
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
	 * Creates an array of byte and copies the contents of a vector of byte into it
	 * 
	 * @param v
	 * @return an array of bytes copied from a vector of byte
	 */
	public static byte[] toByteArray(Vector<Byte> v) {
		byte[] result = new byte[v.size()];
		for (int i=0; i<v.size(); i++) {
			result[i] = v.get(i);
		}
		return result;
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
	 *            an array of numbers
	 * @return sum of the array
	 */
	public static Float sum(Float[] array) {
	    Float sum = Float.valueOf(0);
	    
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
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of integers
	 * @return sum of the array
	 */
	public static int sum(byte[] array) {
		int sum = 0;

		for (int i = 0; i<array.length; i++) {
			sum += array[i];
		}

		return sum;
	}

	/**
	 * Calculates the sum of an array
	 * 
	 * @param array
	 *            an array of long
	 * @return sum of the array
	 */
	public static long sum(long[] array) {
		long sum = 0;

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
	 * @param array
	 * @param factor
	 *            multiply every value in the array by this number
	 * @return
	 */
	public static double[] multiply(final double[] array, final double factor) {
		double[] mult = new double[array.length];
		for (int i = 0; i < mult.length; i++) {
			mult[i] = array[i] * factor;
		}
		return mult;
	}


	/**
	 * @param n
	 *            number of points for the moving average
	 * @param array
	 * @param skipNaN
	 *            If true, every moving average will be composed of n points, or NaN; defaults to removing Nan for the average
	 * @return an array of moving averages with moving average of n sequential points
	 */
	public static double[] movingAverageForward(int n, double[] array, boolean skipNaN) {
		double[] ma = new double[array.length];
		ArrayList<Double> tmp = new ArrayList<Double>();
		for (int i = 0; i < array.length; i++) {
			if (!skipNaN || !Double.isNaN(array[i])) {
				tmp.add(array[i]);
				if (tmp.size() >= n) {
					ma[i] = mean(toDoubleArray(tmp), true);
					tmp.remove(0);
				} else {
					ma[i] = mean(toDoubleArray(tmp), true);
				}
			} else {
				ma[i] = Double.NaN;

			}
		}
		return ma;
	}

	/**
	 * @param n
	 *            number of points for the moving average
	 * @param array
	 * @return an array of moving averages with moving average of n sequential points
	 */
	public static double[] movingMedianForward(int n, double[] array) {
		double[] ma = new double[array.length];
//		double[] a = new double[n];
		ArrayList<Double> tmp = new ArrayList<Double>();
		for (int i = 0; i < array.length; i++) {
			tmp.add(array[i]);
			if (i >= n) {
				ma[i] = (double) median(toDoubleArray(tmp));
				tmp.remove(0);
			} else {
				ma[i] = Double.NaN;
			}
		}
		return ma;
	}
	
	public static double[] scale(double[] array) {
		return scale(array, 0);
	}
	
	
	/**
	 * @param array
	 * @param minForce
	 *            the minimum of the returned array, values scaled to have this value as min
	 * @return
	 */
	public static double[] scaleMinTo(double[] array, final double minForce) {
		double min = Array.min(Array.removeNaN(array));

		if (min < minForce) {
			min = minForce - min;
		} else {
			min = -1 * (min - minForce);
		}
		double[] scaled = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			scaled[i] = array[i] + min;
		}
		return scaled;
	}

	/**
	 * @param array
	 *            scale the values of this array between minForce and minForce +1
	 * @param minForce
	 *            the minimum value (max will be this plus 1)
	 * @return
	 */
	public static double[] scale(final double[] array, final double minForce) {
		double max = Array.max(array);
		double min = Array.min(array);
		double[] scaled = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			scaled[i] = (array[i] - min) / (max - min);
			scaled[i] += minForce;
		}
		return scaled;
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
	 * Calculates the mean of a float array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static float mean(float[] array, boolean ignoreNaN) {
		float sum;
		int count;
		
		sum = 0;
		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (!Float.isNaN(array[i]) || !ignoreNaN) {
				sum += array[i];
				count++;
			}
		}
		
		if (count == 0) {
			return Float.NaN;
		}
		
		return sum/(float)count;
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
	public static double mean(byte[] array) {
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
	 * Calculates the mean of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static Float mean(Float[] array) {
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
	public static double variance(int[] array) {
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
	public static double variance(double[] array) {
		double avg = mean(array);
		double sum = 0;

		for (int i = 0; i<array.length; i++) {
			sum += Math.pow(avg-array[i], 2);
		}

		return sum/(array.length-1);
	}

	/**
	 * Calculates the mean of an array if the sum is already known
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mean of the array
	 */
	public static double mean(double[] array, double sum) {
		return sum / array.length;
	}

	/**
	 * Calculates the sum of squares of an array if the mean is already known
	 * 
	 * @param array
	 *            an array of numbers
	 * @param avg
	 *            precomputed average of the array
	 * @return variance of the array
	 */
	public static double sumSq(double[] array, double avg) {
		double sum = 0;
		for (int i = 0; i < array.length; i++) {
			sum += Math.pow(avg - array[i], 2);
		}
		return sum;
	}

	/**
	 * Calculates the variance of an array if the sum of squares is known
	 * 
	 * @param array
	 *            an array of numbers
	 * @param avg
	 *            precomputed average of the array
	 * @return variance of the array
	 */
	public static double variance(double[] array, double sumSq) {
		return sumSq / (array.length - 1);
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
	public static double stdev(int[] array) {
		return Math.sqrt(variance(array));
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
	public static double stdev(Double[] array, boolean removeNaN) {
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
	    
	    return Math.sqrt(sum/(double)(count-1));
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
	public static Float stdev(Float[] array, boolean removeNaN) {
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
	 * Standardizes (calculates z-scores) for an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of z-scores
	 */
	public static float[] normalize(float[] array) {
		float[] newData = new float[array.length];
		float mean = Array.mean(array, true);
		float stdev = Array.stdev(array, false);

		for (int i = 0; i < newData.length; i++) {
			newData[i] = (array[i] - mean) / stdev;
		}

		return newData;
	}
	
	/**
	 * Normalizes (calculates z-scores) for an array of numbers using separate standard deviations for positive and negative numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of sign-specific z-scores
	 */
	public static double[] normalizeSigned(double[] array) {
		double[] newData;
		double mean = 0;
		double stdevPositive, stdevNegative;
		DoubleVector positives, negatives;

		negatives = new DoubleVector();
		positives = new DoubleVector();
		for (int i = 0; i<array.length; i++) {
			if (array[i] < 0) {
				negatives.add(array[i]);
				negatives.add(-1*array[i]);
			} else {
				positives.add(array[i]);
				positives.add(-1*array[i]);
			}
		}

		stdevNegative = Array.stdev(negatives.toArray());
		stdevPositive = Array.stdev(positives.toArray());
		
		newData = new double[array.length];
		for (int i = 0; i<newData.length; i++) {
			if (array[i] < 0) {
				newData[i] = (array[i]-mean)/stdevNegative;
			} else {
				newData[i] = (array[i]-mean)/stdevPositive;
			}
		}

		return newData;
	}
	
	/**
	 * Returns the quantiles of an array
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of quantiles
	 */
	public static double[] quantiles(double[] array) {
		double[] quantiles;
		int[] order;
		
		order = Sort.quicksort(array);
		quantiles = new double[array.length];
		for (int i = 0; i<quantiles.length; i++) {
			quantiles[order[i]] = ((double)i+1)/((double)quantiles.length+1);
		}

		return quantiles;
	}

	/**
	 * Returns the kurtosis of an array
	 * @param array
	 * @return
	 */
	public static double skewness(double[] array) {
		double skew = -1;
		double mean = Array.mean(array);
		double sd = Array.stdev(array);
		double m3;
		double n = array.length;

		m3 = 0;
		for (int i = 0; i<n; i++) {
			m3 += Math.pow((array[i]-mean)/sd, 3);
		}
		skew = m3 * n / (n-1) / (n-2);

		return skew;
	}

	/**
	 * Returns the kurtosis of an array
	 * @param array
	 * @return
	 */
	public static double kurtosis(double[] array) {
		double kurt = -1;
		double mean = Array.mean(array);
		double m2, s, m4s;
		double n = array.length;

		m2 = 0;
		for (int i = 0; i<n; i++) {
			m2 += Math.pow(array[i]-mean, 2);
		}
		m2 /= (n-1);
		s = Math.sqrt(m2);

		m4s = 0;
		for (int i = 0; i<array.length; i++) {
			m4s += Math.pow((array[i]-mean)/s, 4);
		}

		kurt = n*(n+1)/((n-1)*(n-2)*(n-3))*m4s-3*Math.pow(n-1, 2)/((n-2)*(n-3));

		return kurt;
	}
	
	/**
	 * Inverse-normalizes an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of inverse-normalized values
	 */
	public static double[] inverseNormalize(double[] array) {
		double[] probits, quantiles;
		
		quantiles = quantiles(array);
		probits = new double[array.length];
		for (int i = 0; i<probits.length; i++) {
			if (quantiles[i]<0.5) {
				probits[i] = ProbDist.NormDistReverse(quantiles[i]*2)*-1;
			} else {
				quantiles[i] = 1-quantiles[i];
				probits[i] = ProbDist.NormDistReverse(quantiles[i]*2)*1;
			}
		}

		return probits;
	}
	
	/**
	 * Inverse-normalizes an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return array of inverse-normalized values
	 */
	public static double[] inverseTdist(double[] array, int df) {
		double[] newValues, quantiles;
		
		quantiles = quantiles(array);
		newValues = new double[array.length];
		for (int i = 0; i<newValues.length; i++) {
			if (quantiles[i]<0.5) {
				newValues[i] = ProbDist.TDistReverse(quantiles[i]*2, df)*-1;
			} else {
				quantiles[i] = 1-quantiles[i];
				newValues[i] = ProbDist.TDistReverse(quantiles[i]*2, df)*1;
			}
		}

		return newValues;
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
	 * Determines the specified exclusive quantile of an array of numbers<br />
	 * Returns a number guaranteed to be a member of the given array.<br />
	 * This function matches Excel's QUARTILE.EXC function.
	 * 
	 * @param array
	 *            an array of numbers
	 * @param q
	 *            exclusive quantile to be determined
	 * @return specified exclusive quantile of the array
	 */
	public static double quantExclusive(int[] array, double q) {
		if (array.length == 0) return Double.NaN;

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
	 * Determines the specified exclusive quantile of an array of numbers<br />
	 * Returns a number guaranteed to be a member of the given array.<br />
	 * This function matches Excel's QUARTILE.EXC function.
	 * 
	 * @param array
	 *            an array of numbers
	 * @param q
	 *            exclusive quantile to be determined
	 * @return specified exclusive quantile of the array
	 */
	public static double quantExclusive(double[] array, double q) {
		if (array.length == 0) return Double.NaN;

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
     * Determines the specified exclusive quantile of an array of numbers<br />
     * Returns a number guaranteed to be a member of the given array.<br />
     * This function matches Excel's QUARTILE.EXC function.
     * 
     * @param array
     *            an array of numbers
     * @param q
     *            exclusive quantile to be determined
     * @return specified exclusive quantile of the array
     */
	public static float quantExclusive(float[] array, float q) {
	    if (array.length == 0) return Float.NaN;
	    
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
	 * Determines the specified quantile of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @param q
	 *            quantile to be determined
	 * @return specified quantile of the array
	 */
	public static int quantWithExtremeForTie(int[] array, double q) {
		int keys[] = Sort.quicksort(array);

		try {
			if (q>1||q<0) {
				return (0);
			} else {
				double index = (array.length+1)*q;
				if (index-(int)index==0) {
					return array[keys[(int)index-1]];
				} else if (q < 0.5) {
					return array[keys[(int)Math.floor(index)-1]];
				} else {
					return array[keys[(int)Math.ceil(index)-1]];
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
	 * Determines the specified quantiles of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @param q
	 *            quantiles to be determined
	 * @return specified quantiles of the array
	 */
	public static double[] quantsExclusive(double[] array, double[] qs) {
		int keys[] = Sort.quicksort(array);
		double[] quantiles;

		quantiles = new double[qs.length];
		for (int i = 0; i < quantiles.length; i++) {
			try {
				if (qs[i] > 1 || qs[i] < 0) {
					quantiles[i] = -1;
				} else {
					double index = (array.length + 1) * qs[i];
					if (index - (int) index == 0) {
						quantiles[i] = array[keys[(int) index - 1]];
					} else {
						quantiles[i] = (double) (qs[i] * array[keys[(int) Math.floor(index) - 1]] + (1 - qs[i]) * array[keys[(int) Math.ceil(index) - 1]]);
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
	 * Determines the median absolute difference of an array of double
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mad of the array
	 */
	public static double mad(double[] array) {
		return mad(array, 1);
	}

	
	
	/**
	 * Determines the median absolute difference of an array of double
	 * 
	 * @param array
	 *            an array of numbers
	 * @return mad of the array
	 */
	public static double mad(double[] array, double constant ) {
		double median = (quantExclusive(array, 0.50));
		double[] tmp = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			tmp[i] = Math.abs(array[i] - median);
		}
		return (quantExclusive(tmp, 0.50)) * constant;
	}
	
	/**
	 * Determines the median of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return median of the array
	 */
	public static double median(int[] array) {
		return (quantExclusive(array, 0.50));
	}
	
	/**
	 * Determines the median of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return median of the array
	 */
	public static double median(double[] array) {
		return (quantExclusive(array, 0.50));
	}
	
	/**
	 * Determines the median of an array of numbers
	 * 
	 * @param array
	 *            an array of numbers
	 * @return median of the array
	 */
	public static float median(float[] array) {
	    return (quantExclusive(array, 0.50f));
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
	 * @param a
	 * @param b
	 *            return the distance array from this value
	 * @return
	 */
	public static double[] distFrom(double[] a, double b) {
		return abs(minus(a, b));
	}

	/**
	 * @param a
	 * @param minus
	 *            subtract this value from every entry in the array
	 * @return
	 */
	public static double[] minus(double[] a, double minus) {
		double[] minusA = new double[a.length];
		for (int i = 0; i < minusA.length; i++) {
			minusA[i] = a[i] - minus;
		}
		return minusA;
	}
	
	/**
	 * @param a
	 * @return absolute value array
	 */
	public static double[] abs(double[] a) {
		double[] abs = new double[a.length];
		for (int i = 0; i < abs.length; i++) {
			abs[i] = Math.abs(a[i]);
		}
		return abs;
	}
	/**
	 * Computes x[i + lag] - x[i]
	 */
	public static double[] Diff(double[] x, int lag)// , uint differences = 1)
	{
		double[] diff = new double[x.length - lag];
		for (int i = lag, j = 0; i < x.length; i++, j++) {
			diff[j] = x[i] - x[j];
		}
		return diff;
	}

	public static double PartialSumOfPowers(double[] x, double pow, int start, int size) {
		double sp = 0.0;
		for (int i = start; i < start + size; i++) {
			sp += Math.pow(x[i], pow);
		}
		return sp;
	}

	public static int[] CumulativeSum(int[] x) {
		if (x == null || x.length <= 0) {
			return null;
		}
		int[] cumSum = new int[x.length];
		cumSum[0] = x[0];
		for (int i = 1; i < x.length; i++) {
			cumSum[i] = cumSum[i - 1] + x[i];
		}
		return cumSum;
	}
	
	/**
	 * Computes x[i + lag] - x[i]
	 */
	public static int[] Diff(int[] x, int lag)// , uint differences = 1)
	{
		int[] diff = new int[x.length - lag];
		for (int i = lag, j = 0; i < x.length; i++, j++) {
			diff[j] = x[i] - x[j];
		}
		return diff;
	}

	public static double[] InplaceSub(double[] x, double y) {
		for (int i = 0; i < x.length; i++) {
			x[i] -= y;
		}
		return x;
	}

	public static double[] InplaceAbs(double[] x) {
		for (int i = 0; i < x.length; i++) {
			x[i] = Math.abs(x[i]);
		}
		return x;
	}
	public static double WeightedSumOfSquares(double[] x, double[] w) {
		double wss = 0.0;
		for (int i = 0; i < x.length; i++) {
			double wi = (w == null) ? 1.0 : w[i];
			wss += wi * x[i] * x[i];
		}
		return wss;
	}
	
	/**
	 * Prints all values of an enum, separated by the specified delimiter
	 * @param enumValue An enum class, must be passed as ENUM.class (not the enum itself, but "&lt;enum&gt;.class")
	 * @param delimiter
	 * @return
	 */
	public static <T extends Enum<?>> String toStr(Class<T> enumValue, String delimiter) {
        T[] values = enumValue.getEnumConstants();
        String[] arr = new String[values.length];
        for (int i = 0; i < values.length; i++) {
            arr[i] = values[i].toString();
        }
        return Array.toStr(arr, null, delimiter, null);
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
		boolean commaDelimited;

		count = 0;
		commaDelimited = delimiter.equals(",");
		for (int i = 0; i<array.length; i++) {
			if (display == null || display[i]) {
				if (commaDelimited && array[i].contains(",")) {
					array[i] = "\""+array[i]+"\"";
				}
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
	 * Prints an array of objects separated by a tab
	 * 
	 * @param array
	 *            an array of objects
	 * @return String of printed objects
	 */
	public static String toStr(Object[] array) {
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
	public static String toStr(Object[] array, String delimiter) {
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
	public static String toStr(Object[] array, boolean[] display, String delimiter, String nullValue) {
		String str = "";
		int count;
		boolean commaDelimited;

		count = 0;
		commaDelimited = delimiter.equals(",");
		for (int i = 0; i<array.length; i++) {
			if (display == null || display[i]) {
				String val = array[i].toString();
				if (commaDelimited && val.contains(",")) {
					val = "\""+val+"\"";
				}
				str += (count==0?"":delimiter)+(val==null?nullValue:val);
				count++;
			}
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
	 * Breaks an array nChunks <br>
	 * Warning, one of my first times with the <T> stuff
	 * 
	 * @param array
	 *            the array
	 * @param nChunks
	 *            number of chunks
	 * @param log
	 */

	public static <T> ArrayList<T[]> splitUpArray(T[] array, int nChunks, Logger log) {
		int index = 0;
		if (array.length < nChunks) {
			log.reportError("Error - too many chunks (" + nChunks + ") for " + array.length + " things, setting to" + array.length);
			nChunks = array.length;
		}
		if (nChunks <= 0) {
			log.reportError("Error - not enough chunks (" + nChunks + ") for " + array.length + " things, setting to 1");
			nChunks = 1;
		}
		int[] chunks = Array.splitUpDistributeRemainder(array.length, nChunks, log);
		ArrayList<T[]> da = new ArrayList<T[]>();
		int start = 0;
		for (int i = 0; i < chunks.length; i++) {
			for (int j = 0; j < chunks[i]; j++) {
				index++;
			}
			T[] result = Arrays.copyOf(array, index - start);
			System.arraycopy(array, start, result, 0, result.length);
			da.add(result);
			start = index;
		}
		return da;
	}

	/**
	 * Breaks an array of strings into nChunks. This is geared toward spliting up filenames etc for batching
	 * 
	 * @param strings
	 *            the array
	 * @param nChunks
	 *            number of chunks
	 * @param log
	 */
	public static String[][] splitUpStringArray(String[] strings, int nChunks, Logger log) {
		int index = 0;
		if (strings.length < nChunks) {
			log.reportError("Error - too many chunks (" + nChunks + ") for " + strings.length + " strings, setting to" + strings.length);
			nChunks = strings.length;
		}
		if (nChunks <= 0) {
			log.reportError("Error - not enough chunks (" + nChunks + ") for " + strings.length + " strings, setting to 1");
			nChunks = 1;
		}
		int[] chunks = Array.splitUpDistributeRemainder(strings.length, nChunks, log);
		String[][] stringChunks = new String[chunks.length][];

		for (int i = 0; i < chunks.length; i++) {
			ArrayList<String> chunk = new ArrayList<String>(chunks[i]);
			for (int j = 0; j < chunks[i]; j++) {
				chunk.add(strings[index]);
				index++;
			}
			stringChunks[i] = chunk.toArray(new String[chunk.size()]);
		}
		return stringChunks;
	}
	
	/**
	 * Breaks an array of strings into nChunks with boolean representation, each boolean array has the same length as the original input
	 * 
	 * @param strings
	 *            the array
	 * @param nChunks
	 *            number of chunks
	 * @param log
	 */
	public static boolean[][] splitUpStringArrayToBoolean(String[] strings, int nChunks, Logger log) {
		String[][] stringSplits = splitUpStringArray(strings, nChunks, log);
		boolean[][] stringBoolSplits = new boolean[stringSplits.length][];
		for (int i = 0; i < stringBoolSplits.length; i++) {
			int[] indicesThisChunk = ext.indexLargeFactors(stringSplits[i], strings, true, log, true, false);
			stringBoolSplits[i] = new boolean[strings.length];
			Arrays.fill(stringBoolSplits[i], false);
			for (int j = 0; j < indicesThisChunk.length; j++) {
				stringBoolSplits[i][indicesThisChunk[j]] = true;
			}
		}
		return stringBoolSplits;
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
		int fullBinAmt = (int) Math.floor((double)total / (double)numSplits);
		for (int i = 0; i < numSplits - 1; i++) {
			splits[i] = fullBinAmt;
		}
		splits[numSplits - 1] = total - (numSplits - 1) * fullBinAmt;

		return splits;
	}

	/**
	 * Returns an array splitting a number equally, and distributes the remainder
	 *
	 * @param total
	 *            number to be split into groups
	 * @param numSplits
	 *            number of groups to split total into
	 * @return array of the numbers for each group
	 */
	public static int[] splitUpDistributeRemainder(int total, int numSplits, Logger log) {
		int[] splits = new int[numSplits];
		// TODO, could also redistribute if remainder is small
		for (int i = 0; i < numSplits - 1; i++) {
			splits[i] = (int) Math.floor((double) total / (double) numSplits);
		}
		int remainder = total - (numSplits - 1) * (int) Math.floor((double) total / (double) numSplits);
		if (numSplits > 1 && splits[numSplits - 2] < remainder) {
			splits[numSplits - 1] = splits[numSplits - 2];
			remainder -= splits[numSplits - 1];
			for (int i = 0; i < remainder; i++) {
				splits[i]++;
			}
		} else {
			splits[numSplits - 1] = remainder;
		}
		if (Array.sum(splits) != total) {
			log.reportError("Internal Error - could not properly split up " + total + " into " + numSplits);
			splits = null;
		}
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
	 * Creates an array of Strings and copies the contents of a String[][] into it
	 * 
	 * @param matrix
	 * 				matrix of String
	 * @param delimiter
	 * 				delimiter to use in the concatenated result 
	 * @return an array of Strings from the matrix
	 */
	public static String[] toStringArray(String[][] matrix, String delimiter) {
		String[] array = new String[matrix.length];

		for (int i = 0; i < array.length; i++) {
			array[i] = Array.toStr(matrix[i], delimiter);
		}

		return array;
	}

	/**
	 * Creates an array of Strings and copies the contents of a Vector into it in the specified order
	 * 
	 * @param v
	 *            vector of Strings
	 * @param order
	 *            order of elements
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
	 * Creates an array of Strings and copies the contents of an array of long
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
	 * Creates an array of Strings and copies the contents of an array of double
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
	 * Creates an array of Strings and copies the contents of an array of float
	 * 
	 * @param array
	 *            array of float
	 * @return an array of the converted Strings
	 */
	public static String[] toStringArray(float[] array) {
		String[] new_array;
		
		new_array = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			new_array[i] = array[i]+"";
		}
		
		return new_array;
	}
	
	/**
	 * Creates an array of Strings and copies the contents of a Hashtbable into it (in the correct order)
	 * 
	 * @param hash
	 *            Hashtable of Strings as keys, and their index position as the value
	 * @return an array of Strings from the Hashtable
	 */
	public static String[] toStringArray(Hashtable<String, Integer> hash) {
		Enumeration<String> enumer;
		String[] array;
		String trav;

		enumer = hash.keys();
		array = new String[hash.size()];
		while (enumer.hasMoreElements()) {
			trav = enumer.nextElement();
			array[hash.get(trav).intValue()] = trav;
		}
		
		return array;
	}	
	
	/**
	 * Creates an array of Strings and copies the contents of a Hashtbable into it (in the correct order)
	 * 
	 * @param hash
	 *            Hashtable of Strings as keys, and their index position as the value
	 * @return an array of Strings from the Hashtable
	 */
	public static int[] toIntArray(Hashtable<Integer, Integer> hash) {
		Enumeration<Integer> enumer;
		int[] array;
		int trav;

		enumer = hash.keys();
		array = new int[hash.size()];
		while (enumer.hasMoreElements()) {
			trav = enumer.nextElement();
			array[hash.get(trav).intValue()] = trav;
		}
		
		return array;
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
				hash.put(array[i], array[i]);
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
	
	
	/**
	 * convert a string array into a boolean representation <br>
	 * masks will all be set to false
	 */
	public static BooleanClassifier classifyStringsToBoolean(String[] toClassify, String[] masks) {
		String[] uniqs = unique(toClassify);
		ArrayList<String> uniqNoMask = new ArrayList<String>();
		for (int i = 0; i < uniqs.length; i++) {
			if (ext.indexOfStr(uniqs[i], masks) < 0) {
				uniqNoMask.add(uniqs[i]);
			}
		}
		uniqs = Array.toStringArray(uniqNoMask);
		boolean[][] classified = new boolean[uniqs.length][];
		for (int i = 0; i < classified.length; i++) {
			classified[i] = booleanArray(toClassify.length, false);
		}
		for (int i = 0; i < toClassify.length; i++) {
			int index = ext.indexOfStr(toClassify[i], uniqs);
			if (index >= 0 && ext.indexOfStr(toClassify[i], masks) < 0) {
				classified[index][i] = true;
			}
		}
		return new BooleanClassifier(classified, uniqs);
	}

	public static class BooleanClassifier {
		private boolean[][] classified;
		private String[] titles;

		public boolean[][] getClassified() {
			return classified;
		}

		public BooleanClassifier(boolean[][] classified, String[] titles) {
			super();
			this.classified = classified;
			this.titles = titles;
		}

		public String[] getTitles() {
			return titles;
		}

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
	 *            an array of int
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
	 * Tries to finds the first instance of a given character within an array and
	 * returns either the index or -1 if not found
	 * 
	 * @param array
	 *            an array of char
	 * @param target
	 *            the character to find
	 * @return the index or -1 if not found
	 */
	public static int indexOfChar(char[] array, char target) {
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
	 * Tries to find the instance in a sorted array where all values up to, but not including, that index are less than a given maximum target
	 * <p>
	 * For example, calling {@link Array#indexOfLastMinByte} using (new byte[] {0,1,2,24,25}, 23), would return 3 (all values up to index 3 are less than 23);
	 * 
	 * @param array
	 *            an array of bytes
	 * @param maxByte
	 *            the number to find
	 * 
	 * @return the index, or -1 if all values were greater than the maximum, or the array's length if all values were less than the maximum
	 */
	public static int indexOfFirstMaxByte(byte[] array, byte maxByte) {
		boolean found = false;
		int max = -1;
		for (int i = 0; i < array.length; i++) {
			if (array[i] < maxByte) {
				max = i;
				found = true;
			} else if (found) {
				return max + 1; // can stop here, and set max to the next index
			}
		}
		if (found) {
			return array.length;//all values were less than max
		} else {
			return max; //all values were greater than max
		}
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
	public static byte[] subArray(byte[] array, int start, int stopBefore) {
		byte[] arr;

		if (start < 0 || stopBefore > array.length || stopBefore <= start) {
			System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore + ") indicies for an array");
		}
		arr = new byte[stopBefore - start];
		for (int i = start; i < stopBefore; i++) {
			arr[i - start] = array[i];
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
	
	@SuppressWarnings("unchecked")
	public static <T> T[] subArray(T[] array, boolean[] use) {
		T[] subarray;
		int count;
		
		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n="+array.length+") and the boolean subset (n="+use.length+")");
			return null;
		}
		
		count = 0;
		subarray = (T[]) java.lang.reflect.Array.newInstance(array[0].getClass(), booleanArraySum(use));// new T[booleanArraySum(use)];
		
		for (int i = 0; i < array.length; i++) {
			if (use[i]) {
				subarray[count] = array[i];
				count++;
			}
		}
		
		return subarray;
	}
	
	/**
	 * Creates a new array using only the int values at indices defined by the Integer array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static int[] subArray(int[] array, int[] use) {
		int[] subarray = new int[use.length];
		int currentIndex = 0;
		try {

			for (int i = 0; i < use.length; i++) {
				currentIndex = use[i];
				subarray[i] = array[use[i]];
			}
		} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
			System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n=" + array.length + ")");
			return null;
		}
		return subarray;
	}
	
	
	
	@SuppressWarnings("unchecked")
	public static <T> T[] subArray(T[] array, int[] use) {
		T[] subarray;
		subarray = (T[]) java.lang.reflect.Array.newInstance(array[0].getClass(), use.length);
		int currentIndex = 0;
		try {
			for (int i = 0; i < use.length; i++) {
				currentIndex = use[i];
				subarray[i] = array[use[i]];
			}
		} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
			System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n=" + array.length + ")");
			return null;
		}
		return subarray;
	}
	
	/**
	 * Creates a new array using only the byte values at indices defined by the Integer array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static byte[] subArray(byte[] array, int[] use) {
		byte[] subarray = new byte[use.length];
		int currentIndex = 0;
		try {

			for (int i = 0; i < use.length; i++) {
				currentIndex = use[i];
				subarray[i] = array[use[i]];
			}
		} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
			System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n=" + array.length + ")");
			return null;
		}
		return subarray;
	}

	/**
	 * Creates a new array using only the float values at indices defined by the Integer array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static float[] subArray(float[] array, int[] use) {
		float[] subarray = new float[use.length];
		int currentIndex = 0;
		try {

			for (int i = 0; i < use.length; i++) {
				currentIndex = use[i];
				subarray[i] = array[use[i]];
			}
		} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
			System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n=" + array.length + ")");
			return null;
		}
		return subarray;
	}
	
	/**
	 * Creates a new array using only the strings at indices with a true in the boolean array
	 * 
	 * @param array
	 *            an array of Strings
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static int[] subArray(int[] array, boolean[] use) {
		int[] subarray;
		int count;

		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n="+array.length+") and the boolean subset (n="+use.length+")");
			return null;
		}
		
		count = 0;
		subarray = new int[booleanArraySum(use)];
		for (int i = 0; i<array.length; i++) {
			if (use[i]) {
				subarray[count] = array[i];
				count++;
			}
		}

		return subarray;
	}
	
	/**
	 * Creates a new array using only the bytes at indices with a true in the boolean array
	 * 
	 * @param array
	 *            an array of byte
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static byte[] subArray(byte[] array, boolean[] use) {
		byte[] subarray;
		int count;

		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n="+array.length+") and the boolean subset (n="+use.length+")");
			return null;
		}
		
		count = 0;
		subarray = new byte[booleanArraySum(use)];
		for (int i = 0; i<array.length; i++) {
			if (use[i]) {
				subarray[count] = array[i];
				count++;
			}
		}

		return subarray;
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
		return subArray(array, start, array.length);
	}

	/**
	 * Creates a new array using only the double values at indices with a true in the boolean array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static double[] subArray(double[] array, boolean[] use) {
		double[] subarray;
		int count;

		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n="+array.length+") and the boolean subset (n="+use.length+")");
			return null;
		}
		
		count = 0;
		subarray = new double[booleanArraySum(use)];
		for (int i = 0; i<array.length; i++) {
			if (use[i]) {
				subarray[count] = array[i];
				count++;
			}
		}

		return subarray;
	}
	
	/**
	 * Creates a new array using only the double values at indices with a true in the boolean array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static double[][] subArray(double[][] array, boolean[] use) {
		double[][] subarray;
		int count;

		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n=" + array.length + ") and the boolean subset (n=" + use.length + ")");
			return null;
		}

		count = 0;
		subarray = new double[booleanArraySum(use)][];
		for (int i = 0; i < array.length; i++) {
			if (use[i]) {
				subarray[count] = array[i];
				count++;
			}
		}

		return subarray;
	}
	
	/**
	 * Creates a new array using only the double values at indices defined by the Integer array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static double[] subArray(double[] array, int[] use) {
		double[] subarray = new double[use.length];
		try {
			for (int i = 0; i < use.length; i++) {
				subarray[i] = array[use[i]];
			}
		} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
			System.err.println("Error - out of bounds index for subset");
			return null;
		}
		return subarray;
	}
	
	/**
	 * Creates a new array using only the boolean values at indices defined by the Integer array
	 * 
	 * @param array
	 *            an array of double
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static boolean[] subArray(boolean[] array, int[] use) {
		boolean[] subarray = new boolean[use.length];
		try {
			for (int i = 0; i < use.length; i++) {
				subarray[i] = array[use[i]];
			}
		} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
			System.err.println("Error - out of bounds index for subset");
			return null;
		}
		return subarray;
	}

	/**
	 * Creates a new array using only the indices between start and stop
	 * 
	 * @param array
	 *            an array of double[]
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static double[][] subArray(double[][] array, int start, int stopBefore) {
		double[][] arr;

		if (start < 0 || stopBefore > array.length || stopBefore <= start) {
			System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore + ") indicies for an array");
		}
		arr = new double[stopBefore - start][];
		for (int i = start; i < stopBefore; i++) {
			arr[i - start] = array[i];
		}

		return arr;
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
	 *            an array of boolean
	 * @param start
	 *            first index to use
	 * @param stop
	 *            last index to use
	 * @return the subset of the original array
	 */
	public static boolean[] subArray(boolean[] array, int start, int stopBefore) {
		boolean[] arr;
		if (start<0||stopBefore>array.length||stopBefore<=start) {
			System.err.println("Error - invalid start ("+start+") and stopBefore ("+stopBefore+") indicies for an array");
		}
		arr = new boolean[stopBefore-start];
		for (int i = start; i<stopBefore; i++) {
			arr[i-start] = array[i];
		}
		return arr;
	}
	
	/**
	 * Creates a new array using only the float values at indices with a true in the boolean array
	 * 
	 * @param array
	 *            an array of float
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static float[] subArray(float[] array, boolean[] use) {
		float[] subarray;
		int count;

		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n="+array.length+") and the boolean subset (n="+use.length+")");
			return null;
		}
		
		count = 0;
		subarray = new float[booleanArraySum(use)];
		for (int i = 0; i<array.length; i++) {
			if (use[i]) {
				subarray[count] = array[i];
				count++;
			}
		}

		return subarray;
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
	 * Creates a new array using only the strings at indices with a true in the boolean array
	 * 
	 * @param array
	 *            an array of Strings
	 * @param use
	 *            indices to use
	 * @return the subset of the original array
	 */
	public static String[] subArray(String[] array, boolean[] use) {
		String[] strs;
		int count;

		if (array.length != use.length) {
			System.err.println("Error - mismatched array lengths for the aray (n="+array.length+") and the boolean subset (n="+use.length+")");
			return null;
		}
		
		count = 0;
		strs = new String[booleanArraySum(use)];
		for (int i = 0; i<array.length; i++) {
			if (use[i]) {
				strs[count] = array[i];
				count++;
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
	
	/**
	 * Takes two boolean arrays and sets the first array to the boolean AND of both arrays for a given element.
	 * 
	 * @param aRet First array, and altered array
	 * @param b Second array
	 * @return True if the operation succeeds, false if the operation fails
	 */
	public static boolean booleanArrayAndInPlace(boolean[] aRet, boolean[] b) {
	    if (aRet.length != b.length) {
	        return false;
	    }
	    for (int i = 0; i < aRet.length; i++) {
	        aRet[i] = aRet[i] && b[i];
	    }
	    return true;
	}
	
    public static boolean[] booleanArrayAnd(boolean[] aRet, boolean[] b) {
        if (aRet.length != b.length) {
            return new boolean[0];
        }
        boolean[] ret = new boolean[aRet.length];
        for (int i = 0; i < aRet.length; i++) {
            ret[i] = aRet[i] && b[i];
        }
        return ret;
    }
	
	public static <T extends Comparable<T>> int binarySearch(ArrayList<T[]> list, T[] value, int keyIndex, boolean exact) {
		return binarySearch(list, value, keyIndex, 0, list.size() - 1, exact);
	}
	
	public static <T extends Comparable<T>> int binarySearch(ArrayList<T[]> list, T[] value, int keyIndex, int low, int high, boolean exact) {
		int mid;

		while (low <= high) {
			mid = low + (high - low) / 2;
			if (mid >= list.size()) {
				if (exact) {
					return -9;
				} else {
					return list.size();
				}
			}
			if (list.get(mid)[keyIndex].compareTo(value[keyIndex]) > 0) {
				high = mid - 1;
			} else if (list.get(mid)[keyIndex].compareTo(value[keyIndex]) < 0) {
				low = mid + 1;
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

	public static int binarySearch(String[] array, String value, boolean exact) {
		return binarySearch(array, value, 0, array.length-1, exact);
	}
	
	public static int binarySearch(String[] array, String value, int low, int high, boolean exact) {
		int mid;

		while (low <= high) {
			mid = low + (high - low) / 2;
			if (mid >= array.length) {
				if (exact) {
					return -9;
				} else {
					return array.length;
				}
			}
			if (array[mid].compareTo(value) > 0) {
				high = mid - 1;
			} else if (array[mid].compareTo(value) < 0) {
				low = mid + 1;
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
	
	public static int binarySearch(int[] array, int value, boolean exact) {
		return binarySearch(array, value, 0, array.length-1, exact);
	}

	public static int binarySearch(int[] array, int value, int low, int high, boolean exact) {
		int mid;

		while (low<=high) {
//			System.out.println(array[low]+"\t"+array[high]);
			mid = low+(high-low)/2;
			if (mid >= array.length) {
				if (exact) {
					return -9;
				} else {
					return array.length;
				}
			}
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
	 * Calculates the interquartile range of an array (exclusive)
	 * 
	 * @param array
	 *            an array of numbers
	 * @return iqr of the array
	 */
	public static double iqrExclusive(double[] array) {
		int[] keys = Sort.quicksort(array);

		if (array.length<2) {
			System.err.println("Error - can't calculate an IQR for an array with "+array.length+" datapoint(s)");
			return -1;
		}

		double iqr = 0;
		try {
			iqr = array[keys[(int)Math.floor(array.length*0.75)]]-array[keys[(int)Math.floor(array.length*0.25)]];
		} catch (Exception e) {
			System.err.println("Error calculating IQR");
			e.printStackTrace();
		}
		return iqr;
	}

	/**
	 * Calculates the interquartile range of an array (exclusive)
	 * 
	 * @param array
	 *            an array of numbers
	 * @return iqr of the array
	 */
	public static float iqrExclusive(float[] array) {
		int[] keys = Sort.quicksort(array);
		float iqr = 0;

		if (array.length<2) {
			System.err.println("Error - can't calculate an IQR for an array with "+array.length+" datapoint(s)");
			return -1;
		}
		try {
			iqr = array[keys[(int)Math.floor(array.length*0.75)]]-array[keys[(int)Math.floor(array.length*0.25)]];
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
	 * Returns true if the byte arrays are equal at all positions
	 * 
	 * @param array1
	 *            an array of byte
	 * @param array2
	 *            an array of byte
	 * @param caseSensitive
	 *            boolean flag
	 * @return true if arrays are equal
	 */
	public static boolean equals(byte[] array1, byte[] array2) {
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
	 * @param array
	 * @param sigFigs
	 *            all entries in the array will be rounded the this many sig figs
	 * @return the rounded array
	 */
	public static double[] round(double[] array, int sigFigs) {
		double[] rounded = new double[array.length];
		for (int i = 0; i < rounded.length; i++) {
			rounded[i] = Maths.roundDouble(array[i], sigFigs);
		}
		return rounded;
	}
	
	/**
	 * Removes NaN values from the array
	 * 
	 * @param array
	 *            an array of doubles
	 * @return scrubbed array
	 */
	public static double[] removeNaN(double[] array) {
		boolean[] use;

		use = new boolean[array.length];
		for (int i = 0; i < use.length; i++) {
			use[i] = !Double.isNaN(array[i]);
		}

		return subArray(array, use);
	}

	/**
	 * Removes NaN values from the array
	 * 
	 * @param array
	 *            an array of doubles
	 * @return scrubbed array
	 */
	public static float[] removeNaN(float[] array) {
		boolean[] use;

		use = new boolean[array.length];
		for (int i = 0; i < use.length; i++) {
			use[i] = !Float.isNaN(array[i]);
		}

		return subArray(array, use);
	}
	
	/**
	 * Removes non-finite values from the array
	 * 
	 * @param array
	 *            an array of doubles
	 * @return scrubbed array
	 */
	public static double[] removeNonFinites(double[] array) {
		boolean[] use;

		use = new boolean[array.length];
		for (int i = 0; i < use.length; i++) {
			use[i] = Double.isFinite(array[i]);
		}

		return subArray(array, use);
	}
	
	/**
	 * Removes non-finite values from the array
	 * 
	 * @param array
	 *            an array of floats
	 * @return scrubbed array
	 */
	public static float[] removeNonFinites(float[] array) {
		boolean[] use;

		use = new boolean[array.length];
		for (int i = 0; i < use.length; i++) {
			use[i] = Float.isFinite(array[i]);
		}

		return subArray(array, use);
	}
	
	/**
	 * Removes all instances of a specified value from an array
	 * 
	 * @param array
	 *            an array of bytes
	 * @param valueToRemove
	 *            value to remove from the array
	 * @return filtered array
	 */
	public static byte[] removeAllValues(byte[] array, byte valueToRemove) {
		boolean[] use;
		
		use = new boolean[array.length];
		for (int i = 0; i<use.length; i++) {
			use[i] = array[i] != valueToRemove;
        }

		return subArray(array, use);
	}
	
	/**
	 * Removes all instances of a specified value from an array
	 * 
	 * @param array
	 *            an array of int
	 * @param valueToRemove
	 *            value to remove from the array
	 * @return filtered array
	 */
	public static int[] removeAllValues(int[] array, int valueToRemove) {
		boolean[] use;

		use = new boolean[array.length];
		for (int i = 0; i < use.length; i++) {
			use[i] = array[i] != valueToRemove;
		}

		return subArray(array, use);
	}

	public static double[] getValuesBetween(double[] array, double min, double max) {
		return getValuesBetween(array, min, max, false);
	}

	public static double[] getValuesBetween(double[] array, double min, double max, boolean gteLte) {
		ArrayList<Double> tmp = new ArrayList<Double>();
		for (int i = 0; i < array.length; i++) {
			if (!Double.isNaN(array[i]) && (array[i] > min && array[i] < max || (gteLte && array[i] >= min && array[i] <= max))) {
				tmp.add(array[i]);
			}
		}
		return Array.toDoubleArray(tmp);
	}

	public static float[] getValuesBetween(float[] array, double min, double max, boolean gteLte) {
		ArrayList<Float> tmp = new ArrayList<Float>();
		for (int i = 0; i < array.length; i++) {
			if (!Double.isNaN(array[i]) && (array[i] > min && array[i] < max || (gteLte && array[i] >= min && array[i] <= max))) {
				tmp.add(array[i]);
			}
		}
		return Array.toFloatArray(tmp);
	}

	public static int[] getValuesBetween(int[] array, int min, int max, boolean gteLte) {
		ArrayList<Integer> tmp = new ArrayList<Integer>();
		for (int i = 0; i < array.length; i++) {
			if ((array[i] > min && array[i] < max) || (gteLte && array[i] >= min && array[i] <= max)) {
				tmp.add(array[i]);
			}
		}
		return Array.toIntArray(tmp);
	}

	/**
	 * Reverses the entries of the array, first becomes last, etc
	 */
	public static String[] reverse(String[] forward) {
		String[] reverse = new String[forward.length];
		int index = forward.length - 1;
		for (int i = 0; i < reverse.length; i++) {
			reverse[i] = forward[index];
			index--;
		}
		return reverse;
	}

	/**
	 * Reverses the entries of the array, first becomes last, etc
	 */
	public static int[] reverse(int[] forward) {
		int[] reverse = new int[forward.length];
		int index = forward.length - 1;
		for (int i = 0; i < reverse.length; i++) {
			reverse[i] = forward[index];
			index--;
		}
		return reverse;
	}
	
	/**
	 * Reverses the entries of the array, first becomes last, etc
	 */
	public static boolean[] reverse(boolean[] forward) {
		boolean[] reverse = new boolean[forward.length];
		int index = forward.length - 1;
		for (int i = 0; i < reverse.length; i++) {
			reverse[i] = forward[index];
			index--;
		}
		return reverse;
	}

	/**
	 * Reverses the entries of the array, first becomes last, etc
	 */
	public static double[] reverse(double[] forward) {
		double[] reverse = new double[forward.length];
		int index = forward.length - 1;
		for (int i = 0; i < reverse.length; i++) {
			reverse[i] = forward[index];
			index--;
		}
		return reverse;
	}
	
	/**
	 * Extract the element at the given index in each sub-array and return a single flat array.
	 * 
	 * @param srcArr Source array of type T[][]
	 * @param index Index of elements in subArrays
	 * @return null if <code>srcArr</code> is null or empty, <br />or an array of type T[], with the same length as <code>srcArr</code>, containing only the elements at <code>index</code> in each sub-array of <code>srcArr</code>
	 */
    public static <T> T[] extract(T[][] srcArr, int index) {
        if (srcArr == null) return null;
        if (srcArr.length == 0) return null;
        Class<?> arrClz = srcArr[0][0].getClass();
        @SuppressWarnings("unchecked")
        T[] arr = (T[])java.lang.reflect.Array.newInstance(arrClz, srcArr.length);//new Object[srcArr.length];
        for (int i = 0; i < srcArr.length; i++) {
            arr[i] = srcArr[i][index];
        }
        return arr;
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
		return addIntToArray(value, array, array.length);
	}

	/**
	 * Adds specified integer to a specified index of an array
	 * 
	 * @param array
	 *            an array of integers
	 * @param value
	 *            integer to add to the array 
	 * @param indexOfNewArray
	 *            location of the integer in the new array 
	 * @return new array
	 */
	public static int[] addIntToArray(int value, int[] array, int indexOfNewStr) {
    	int[] new_array = new int[array.length+1];
    
    	for (int i = 0; i<new_array.length; i++) {
    		new_array[i] = i==indexOfNewStr?value:array[i>indexOfNewStr?i-1:i];
    	}
    
    	return new_array;
    }
	
	/**
	 * Adds specified double to the end of an array
	 * 
	 * @param array
	 *            an array of double
	 * @param value
	 *            double to add to the array 
	 * @return new array
	 */
	public static double[] addDoubleToArray(double value, double[] array) {
		return addDoubleToArray(value, array, array.length);
	}

	/**
	 * Adds specified double to a specified index of an array
	 * 
	 * @param array
	 *            an array of double
	 * @param value
	 *            double to add to the array 
	 * @param indexOfNewArray
	 *            location of the double in the new array 
	 * @return new array
	 */
	public static double[] addDoubleToArray(double value, double[] array, int indexOfNewStr) {
		double[] new_array = new double[array.length+1];
    
    	for (int i = 0; i<new_array.length; i++) {
    		new_array[i] = i==indexOfNewStr?value:array[i>indexOfNewStr?i-1:i];
    	}
    
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
	 * Clones an array of int
	 * 
	 * @param array
	 *            the array of int to clone
	 * @return cloned array of int
	 */
	public static int[] clone(int[] array) {
		int[] newArray;

		newArray = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[i];
		}
		
		return newArray;
	}	

	/**
	 * Clones an array of double
	 * 
	 * @param array
	 *            the array of double to clone
	 * @return cloned array of double
	 */
	public static double[] clone(double[] array) {
		double[] newArray;

		newArray = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[i];
		}
		
		return newArray;
	}	

	/**
	 * Clones an array of float
	 * 
	 * @param array
	 *            the array of float to clone
	 * @return cloned array of float
	 */
	public static float[] clone(float[] array) {
		float[] newArray;

		newArray = new float[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[i];
		}
		
		return newArray;
	}	

	/**
	 * Clones an array of boolean
	 * 
	 * @param array
	 *            the array of boolean to clone
	 * @return cloned array of boolean
	 */
	public static boolean[] clone(boolean[] array) {
		boolean[] newArray;

		newArray = new boolean[array.length];
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
	
	public static String[] addPrefixSuffixToArray(String[] array, String prefix, String suffix) {
		String[] newArray;
		
		newArray = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = (prefix==null?"":prefix)+array[i]+(suffix==null?"":suffix);
		}
		
		return newArray;
	}
	
	public static boolean isBimodal(double[] array) {
		return isBimodal(array, 0.5, 100);
	}

	public static boolean isBimodal(double[] array, double percentDropInPeak, int numBins) {
		return isBimodal(array, percentDropInPeak, (max(array)-min(array))/(double)numBins);
	}
	
	public static boolean isBimodal(double[] array, double percentDropInPeak, double binSize) {
		int numBins;
		int[] freqBinCounts,freqBinCountsSmooth;
		double minValue, maxFreq, localMinFreq;

		numBins = (int) ((Array.max(array)-Array.min(array))/binSize);
		minValue = Array.min(array);
		freqBinCounts = new int[numBins];
		for (int i=0; i<array.length; i++) {
			freqBinCounts[(int) (Math.floor(array[i]-minValue)/binSize)]++;
		}

		//smoothing
		freqBinCountsSmooth = new int[freqBinCounts.length];
		for (int i=1; i<numBins-1; i++) {
			freqBinCountsSmooth[i]=(freqBinCounts[i-1]+freqBinCounts[i]+freqBinCounts[i+1])/3;
		}
		freqBinCountsSmooth[0]=(freqBinCounts[0]+freqBinCounts[1])/2;
		freqBinCountsSmooth[numBins-1]=(freqBinCounts[numBins-2]+freqBinCounts[numBins-1])/2;
		
		maxFreq = Double.NEGATIVE_INFINITY;
		localMinFreq = Double.POSITIVE_INFINITY;
		for (int i=0; i<numBins; i++) {
			if (freqBinCountsSmooth[i]>maxFreq) {
				maxFreq = freqBinCountsSmooth[i];
			} else if (freqBinCountsSmooth[i]<(maxFreq*percentDropInPeak) && freqBinCountsSmooth[i]<localMinFreq) {
				localMinFreq = freqBinCountsSmooth[i];
			} else if (freqBinCountsSmooth[i]>=(maxFreq*percentDropInPeak)) {
				return true;
			}
		}
		return false;
	}

	public static boolean isMultimodal(double[] array, double proportionOfLastPeakRequiredForNewLocalMinima, double proportionOfGlobalMaxRequiredForLocalMaxima, double binSize) {
		return getLocalModes(array, proportionOfLastPeakRequiredForNewLocalMinima, proportionOfGlobalMaxRequiredForLocalMaxima, binSize, true).length > 1;
	}
	
	public static double[] getLocalModes(double[] array, double proportionOfLastPeakRequiredForNewLocalMinima, double proportionOfGlobalMaxRequiredForLocalMaxima) {
		return getLocalModes(array, proportionOfLastPeakRequiredForNewLocalMinima, proportionOfGlobalMaxRequiredForLocalMaxima, (max(array)-min(array))/40, true);
	}
	
	public static double[] getLocalModes(double[] array, double proportionOfLastPeakRequiredForNewLocalMinima, double proportionOfGlobalMaxRequiredForLocalMaxima, double binSize, boolean sensitiveToSmallNumbers) {
		int numBins;
		int[] freqBinCounts;
		double[] freqBinCountsSmooth;
		double minValue;
		int[] indicesOfLocalMaxima;
		double[] modes;
		
		numBins = (int) ((Array.max(array)-Array.min(array))/binSize)+1;
		minValue = Array.min(array);
		freqBinCounts = new int[numBins];
		for (int i=0; i<array.length; i++) {
			freqBinCounts[(int)Math.floor((array[i]-minValue)/binSize)]++;
		}

		//smoothing
		freqBinCountsSmooth = new double[freqBinCounts.length];
		for (int i=1; i<numBins-1; i++) {
			freqBinCountsSmooth[i]=(freqBinCounts[i-1]+freqBinCounts[i]+freqBinCounts[i+1])/3;
		}
		if (freqBinCounts.length >= 2) {
			freqBinCountsSmooth[0]=(freqBinCounts[0]+freqBinCounts[1])/2;
			freqBinCountsSmooth[numBins-1]=(freqBinCounts[numBins-2]+freqBinCounts[numBins-1])/2;
		}

		if (sensitiveToSmallNumbers) {
			proportionOfGlobalMaxRequiredForLocalMaxima = Math.max(proportionOfGlobalMaxRequiredForLocalMaxima, Math.min(0.50, proportionOfGlobalMaxRequiredForLocalMaxima*proportionOfGlobalMaxRequiredForLocalMaxima*300/array.length));
			if (array.length < 50) {
//				System.out.println(array.length+"\t"+proportionOfGlobalMaxRequiredForLocalMaxima);
			}
		}
		
		indicesOfLocalMaxima = getIndicesOfLocalMaxima(freqBinCountsSmooth, proportionOfLastPeakRequiredForNewLocalMinima, proportionOfGlobalMaxRequiredForLocalMaxima);
		modes = new double[indicesOfLocalMaxima.length];
		for (int i = 0; i < modes.length; i++) {
			modes[i] = minValue+indicesOfLocalMaxima[i]*binSize+binSize/2;
		}
		return modes;
	}

	public static double[] smooth(double[] array, int numOfPositionsInOneDirection) {
		double[] smoothed;
		
		smoothed = new double[array.length];

		return smoothed;
	}
	
	public static int[] getIndicesOfLocalMaxima(double[] array, double proportionOfLastPeakRequiredForNewLocalMinima, double proportionOfGlobalMaxRequiredForLocalMaxima) {
		double globalMax, localMax, localMin;
		int indexOfLocalMax;
		IntVector indicesOfMaxima;
		
		indicesOfMaxima = new IntVector();
		
		globalMax = max(array);
		if (globalMax == min(array)) {
			globalMax++;
		}
		
		localMax = Double.NEGATIVE_INFINITY;
		localMin = Double.POSITIVE_INFINITY;
		indexOfLocalMax = -1;
		for (int i=0; i<array.length; i++) {
			if (array[i]>localMax) {
				localMax = array[i];
				indexOfLocalMax = i;
			}
			if (localMax >= globalMax*proportionOfGlobalMaxRequiredForLocalMaxima && array[i]<=(localMax*proportionOfLastPeakRequiredForNewLocalMinima)) {
//				System.out.println("localMax="+localMax+" at index "+indexOfLocalMax);
				localMin = array[i];
				indicesOfMaxima.add(indexOfLocalMax);
				localMax = Double.NEGATIVE_INFINITY;
			}
			if (array[i]>=(globalMax*proportionOfGlobalMaxRequiredForLocalMaxima)) {
				if (localMin != Double.POSITIVE_INFINITY) {
//					System.out.println("localMin="+localMin);
				}
				localMin = Double.POSITIVE_INFINITY;
			}
		}
		if (localMin != Double.POSITIVE_INFINITY && localMax >= globalMax*proportionOfGlobalMaxRequiredForLocalMaxima && indexOfLocalMax != -1) {
//			System.out.println("localMax="+localMax+" at index "+indexOfLocalMax);
			indicesOfMaxima.add(indexOfLocalMax);
		}
		
		return indicesOfMaxima.toArray();
	}
	
	public static boolean[] indicesToBooleanArray(int[] rowsToKeep, int sizeOfArray) {
		boolean[] array;
		
		array = booleanArray(sizeOfArray, false);
		for (int i = 0; i < rowsToKeep.length; i++) {
			array[rowsToKeep[i]] = true;
		}
		
		return array;
	}
	
	public static int[] booleanArrayToIndices(boolean[] keeps) {
		
		int[] indices = new int[booleanArraySum(keeps)];
		int i = 0;
		for (int j = 0; j < keeps.length; j++) {
			if (keeps[j]) indices[i++] = j;
		}
		return indices;
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

	/**
	 * Creates an array of char and copies the contents of an array of String
	 * 
	 * @param array
	 *            array of String
	 * @return an array of the converted String
	 */
	public static char[] toCharArray(String[] array) {
		char[] newArray;
		
		newArray = new char[array.length];
		for (int i = 0; i < newArray.length; i++) {
			if (array[i].length() != 1) {
				System.err.println("Error - cannot convert string to char since it is longer than 1 byte: "+array[i]);
			}
			newArray[i] = array[i].charAt(0);
		}

		return newArray;
	}

	/**
	 * Transposes a List<List<>> i.e. a 2 dimensional list
	 *
	 * @param table
	 *            : the 2D list to be transposed
	 * @param <T>
	 *            : generic template
	 * @param table : the 2D list to be transposed
	 * @param <T> : generic template
	 * @return a {@link List<List>>} which is transposed
	 */
	public static <T> List<List<T>> transpose(List<List<T>> table) {
		List<List<T>> ret = new ArrayList<List<T>>();
		final int N = table.get(0).size();
		for (int i = 0; i < N; i++) {
			List<T> col = new ArrayList<T>();
			for (List<T> row : table) {
				col.add(row.get(i));
			}
			ret.add(col);
		}
		return ret;
	}

	
	
	public static double[] concatDubs(double[] first, double[] other) {
		int totalLength = first.length+other.length;
	
		double[] result = new double[totalLength];
		int index =0;
		for (int i = 0; i < first.length; i++) {
			result[index] =first[i];
			index++;
		}
		for (int i = 0; i < other.length; i++) {
			result[index] =other[i];
			index++;
		}
		return result;
	}

	/**
	 * Remove common elements from the front and/or back of strings Ex
	 * 
	 * [acatcat3cat, acatcat4cat] -> [3,4] <br>
	 * Not supremely tested though...
	 * 
	 * @param front
	 * @param back
	 * @return
	 */
	public static String[] untag(final String[] toUniq, boolean front, boolean back) {
		String[] uniq = new String[toUniq.length];
		boolean findingstart = true;
		int startIndex = 0;
		for (int i = 0; i < uniq.length; i++) {
			uniq[i] = toUniq[i];
		}

		if (front) {
			while (findingstart && startIndex < toUniq[0].length()) {
				String start = toUniq[0].charAt(startIndex) + "";
				for (int i = 0; i < uniq.length; i++) {
					if (startIndex >= toUniq[i].length() || !(toUniq[i].charAt(startIndex) + "").equals(start)) {
						findingstart = false;
						startIndex--;
						break;
					}
				}
				startIndex++;
			}
		}
		boolean findingstop = true;
		int stopIndex = 1;
		if (back) {
			while (findingstop && toUniq[0].length() > stopIndex + 1) {
				String stop = toUniq[0].charAt(toUniq[0].length() - stopIndex) + "";
				for (int i = 0; i < uniq.length; i++) {
					if (toUniq[i].length() <= stopIndex || !(toUniq[i].charAt(toUniq[i].length() - stopIndex) + "").equals(stop)) {
						stopIndex--;
						findingstop = false;
						break;
					}
				}
				stopIndex++;
			}
		}
		for (int i = 0; i < uniq.length; i++) {
			uniq[i] = toUniq[i].substring(startIndex, toUniq[i].length() - Math.max(0, stopIndex - 1));
		}

		return uniq;
	}
	
	
	
	/**
	 * @param array
	 * @param front
	 *            tag on to the front of each entry
	 * @param back
	 *            tag on to the back of each entry
	 * @return
	 */
	public static String[] tagOn(String[] array, String front, String back) {
		String[] arrayTagged = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			String tmp = array[i];
			if (front != null) {
				tmp = front + tmp;
			}
			if (back != null) {
				tmp = tmp + back;
			}
			arrayTagged[i] = tmp;

		}
		return arrayTagged;
	}

	/**
	 * Function to concatenate an arbitrary number of Arrays.
	 *
	 * @param first
	 *            the first array
	 * @param rest
	 *            rest all arrays
	 * @param <T>
	 *            generic: can take any type of object
	 * @return a flat array containing the elements of all given arrays
	 */
	public static <T> T[] concatAll(T[] first, T[]... rest) {
		int totalLength = first.length;
		for (T[] array : rest) {
			totalLength += array.length;
		}
		T[] result = Arrays.copyOf(first, totalLength);
		int offset = first.length;
		for (T[] array : rest) {
			System.arraycopy(array, 0, result, offset, array.length);
			offset += array.length;
		}
		return result;
	}
	
    public static int[][] nCr_indices(int n, int r) {
    	Vector<int[]> v;
        boolean done;
        int[] res, indices;
        
        res = new int[r];
        for (int i = 0; i < res.length; i++) {
            res[i] = i+1;
        }
        
        done = false;
        v = new Vector<int[]>();
    	indices = new int[r];
        while (!done) {
        	for (int i = 0; i < r; i++) {
        		indices[i] = res[i]-1;
			}
        	v.add(indices.clone());
            done = getNext(res, n, r);
        }
        
        return Matrix.toMatrix(v);
    }

    public static boolean getNext(int[] num, int n, int r) {
        int target = r - 1;
        num[target]++;
        if (num[target] > ((n - (r - target)) + 1)) {
            // Carry the One
            while (num[target] > ((n - (r - target)))) {
                target--;
                if (target < 0) {
                    break;
                }
            }
            if (target < 0) {
                return true;
            }
            num[target]++;
            for (int i = target + 1; i < num.length; i++) {
                num[i] = num[i - 1] + 1;
            }
        }
        return false;
    }

	public static double lambda(double[] pvals) {
		return ProbDist.ChiDistReverse(Array.median(pvals), 1)/ProbDist.ChiDistReverse(0.50, 1);
	}
	
	/**
	 * Creates an array of float and copies the contents of an ArrayList of float into it
	 * 
	 * @param al
	 * @return an array of floats copied from a ArrayList of floats
	 */
	public static float[] toFloatArray(ArrayList<Float> al) {
		float[] result = new float[al.size()];
		for (int i = 0; i < al.size(); i++) {
			result[i] = al.get(i);
		}
		return result;
	}

	/**
	 * Creates an array of double and copies the contents of an ArrayList of double into it
	 * 
	 * @param al
	 * @return an array of doubles copied from a ArrayList of doubles
	 */
	public static double[] toDoubleArray(ArrayList<Double> al) {
		double[] result = new double[al.size()];
		for (int i = 0; i < al.size(); i++) {
			result[i] = al.get(i);
		}
		return result;
	}
	
	/**
	 * Creates an array of int and copies the contents of an ArrayList of int into it
	 * 
	 * @param al
	 * @return an array of int copied from a ArrayList of int
	 */
	public static int[] toIntArray(ArrayList<Integer> al) {
		int[] result = new int[al.size()];
		for (int i = 0; i < al.size(); i++) {
			result[i] = al.get(i);
		}
		return result;
	}

	/**
	 * Creates an array of byte and copies the contents of an ArrayList of byte into it
	 * 
	 * @param al
	 * @return an array of byte copied from a ArrayList of byte
	 */
	public static byte[] toByteArray(ArrayList<Byte> al) {
		byte[] result = new byte[al.size()];
		for (int i = 0; i < al.size(); i++) {
			result[i] = al.get(i);
		}
		return result;
	}
	
	/**
	 * Creates an array of boolean and copies the contents of an ArrayList of Boolean into it
	 * 
	 * @param al
	 * @return an array of boolean copied from an ArrayList of Boolean
	 */
	public static boolean[] toBooleanArray(ArrayList<Boolean> al) {
		boolean[] result = new boolean[al.size()];
		for (int i = 0; i < al.size(); i++) {
			result[i] = al.get(i);
		}
		return result;
	}
	
	/**
	 * Creates an array of boolean and copies the contents of a Vector of Boolean into it
	 * 
	 * @param v
	 * @return an array of boolean copied from a Vector of Boolean
	 */
	public static boolean[] toBooleanArray(Vector<Boolean> v) {
		boolean[] result = new boolean[v.size()];
		for (int i = 0; i < v.size(); i++) {
			result[i] = v.get(i);
		}
		return result;
	}

	/**
	 * Creates an array of Integer and copies the contents of an ArrayList of Integer into it
	 * 
	 * @param al
	 * @return an array of Integer copied from a ArrayList of Integer
	 */
	public static int[] toIntegerArray(ArrayList<Integer> al) {
		int[] result = new int[al.size()];
		for (int i = 0; i < al.size(); i++) {
			result[i] = al.get(i);
		}
		return result;
	}
	
	/**
	 * Takes the log base 2 of every element in the array a
	 * 
	 * @param a
	 *            an array of doubles
	 * @return the array containing the log base 2 value of each element in the array
	 */
	public static double[] log2(double[] a) {
		double[] log2A = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			log2A[i] = Maths.log2(a[i]);
		}
		return log2A;
	}

	public static boolean containsMissingValue(double[] array) {
		for (int i = 0; i < array.length; i++) {
			if (!ext.isValidDouble(array[i]+"")) {
				return true;
			}
		}
		return false;
	}
	
	public static String[] removeMissingValues(String[] array) {
	    ArrayList<String> valid = new ArrayList<String>();
	    for (String str : array) {
	        if (str == null || !ext.isValidDouble(str)) continue;
	        valid.add(str);
	    }
	    return valid.toArray(new String[valid.size()]);
	}
	
    public static String[] combine(String[] array1, String[] array2) {
        String[] newArray = Arrays.copyOf(array1, array1.length + array2.length);
        for (int i = array1.length; i < array1.length + array2.length; i++) {
            newArray[i] = array2[i - array1.length];
        }
        return newArray;
    }
    
	public static void main(String[] args) {
	    double[] data = {11.8, 0.93, 1.76, 14, 16.5, 17.1, 32.5, 33.4, 16.8, 21.5, 13.1, 22.2, 22.2, 16, 16.2};
//        float[] data = {11.8f, 0.93f, 1.76f, 14, 16.5f, 17.1f, 32.5f, 33.4f, 16.8f, 21.5f, 13.1f, 22.2f, 22.2f, 16, 16.2f};
        
        System.out.println(Array.toStr(quantiles(data)));
        
	    
//	    double alleleFreq = 0.2;
//	    double stdev = 0.12;
//	    double[] array = new double[10000];
//	    for (int i = 0; i<array.length; i++) {
//	    	array[i] = 0;
//	    	for (int j = 0; j < 2; j++) {
//		    	array[i] += Math.random() < alleleFreq?0.5:0;
//			}
//	    	array[i] += (Math.random()<0.5?-1:1)*ProbDist.NormDistReverse(Math.random())*stdev;
//        }
//	    System.out.println(Array.toStr(getLocalModes(array, 0.1, 0.15, 0.01, false)));
//	    
//	    Files.writeList(Array.toStringArray(array), "oi.xln");

    }

}
