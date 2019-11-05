package org.pankratzlab.common;

import java.io.UnsupportedEncodingException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.pankratzlab.common.stats.Maths;
import org.pankratzlab.common.stats.ProbDist;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;
import com.google.common.primitives.Ints;

public class ArrayUtils {

  private static final int MINI = -999;
  private static final double MIND = Double.NaN;
  private static final int MINF = -999;
  private static final int MINB = Byte.MAX_VALUE;

  private static boolean badLength(int length) {
    return length == 0;
  }

  /**
   * Return the minimum in an array of integers
   *
   * @see {@link Ints#min(int...)}
   * @param array array of integers
   * @return the minimum
   */
  public static int min(int[] array) {
    if (badLength(array.length)) {
      return MINI;
    }

    return Ints.min(array);
  }

  /**
   * Return the minimum in an array of doubles
   *
   * @see {@link Doubles#min(double...)}
   * @param array array of doubles
   * @return the minimum
   */
  public static double min(double[] array) {
    if (badLength(array.length)) {
      return MIND;
    }
    return Doubles.min(array);
  }

  /**
   * Return the minimum, non-{@link Float#NaN} value.
   *
   * @see {@link Floats#min(float...)}
   * @param array array of floats
   * @return the minimum
   */
  public static float min(float[] array) {
    if (badLength(array.length)) {
      return MINF;
    }
    float min = Float.POSITIVE_INFINITY;
    for (float element : array) {
      if (Float.isNaN(element)) {
        continue;
      }
      min = Math.min(min, element);
    }
    return min;
  }

  /**
   * Return the minimum and maximum non-{@link Float#NaN} values, respectively, in an array of
   * floats
   *
   * @param array array of floats
   * @return the minimum and maximum, in that order
   */
  public static float[] minMax(float[] array) {
    if (badLength(array.length)) {
      return new float[] {MINF, MINF};
    }
    float min = Float.POSITIVE_INFINITY;
    float max = Float.NEGATIVE_INFINITY;
    for (float element : array) {
      if (Float.isNaN(element)) {
        continue;
      }
      min = Math.min(min, element);
      max = Math.max(max, element);
    }
    return new float[] {min, max};
  }

  /**
   * Return the minimum and maximum non-{@link Double#NaN} values, respectively, in an array of
   * doubles
   *
   * @param array array of doubles
   * @return the minimum and maximum, in that order
   */
  public static double[] minMax(double[] array) {
    if (badLength(array.length)) {
      return new double[] {MIND, MIND};
    }
    double min = Double.POSITIVE_INFINITY;
    double max = Double.NEGATIVE_INFINITY;
    for (double element : array) {
      if (Double.isNaN(element)) {
        continue;
      }
      min = Math.min(min, element);
      max = Math.max(max, element);
    }
    return new double[] {min, max};
  }

  /**
   * Return the minimum in an array of bytes
   *
   * @param array array of integers
   * @return the minimum
   */
  public static byte min(byte[] array) {
    if (badLength(array.length)) {
      return MINB;
    }
    byte min = Byte.MAX_VALUE;
    for (int i = 1; i < array.length; i++) {
      min = (byte) Math.min(array[i], min);
    }
    return min;
  }

  /**
   * Return the index of the minimum non-{@link Float#NaN} value in an array of floats
   *
   * @param array array of floats
   * @return the index of the minimum or -1 if an empty array is given
   */
  public static int minIndex(float[] array) {
    float min;
    int index = -1;

    if (badLength(array.length)) {
      return -1;
    }
    min = Float.POSITIVE_INFINITY;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i]) && Float.compare(array[i], min) < 0) {
        min = array[i];
        index = i;
      }
    }
    return index;
  }

  /**
   * Return the index of the minimum in an array of doubles
   *
   * @param array array of doubles
   * @return the minimum
   */
  public static int minIndex(double[] array) {
    return index(array, true);
  }

  /**
   * Return the index of the maximum in an array of integers
   *
   * @param array array of integers
   * @return the minimum
   */
  public static int maxIndex(double[] array) {
    return index(array, false);
  }

  private static int index(double[] array, boolean findMin) {
    if (badLength(array.length)) {
      return -1;
    }
    double v = findMin ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
    int index = -1;

    for (int i = 0; i < array.length; i++) {
      if (Double.isNaN(array[i])) {
        continue;
      }
      int c = Double.compare(array[i], v);
      if (findMin ? c < 0 : c > 0) {
        v = array[i];
        index = i;
      }
    }
    return index;
  }

  /**
   * Return the maximum in an array of integers
   *
   * @see {@link Ints#max(int...)}
   * @param array array of integers
   * @return the maximum
   */
  public static int max(int[] array) {
    if (badLength(array.length)) {
      return MINI;
    }

    return Ints.max(array);
  }

  /**
   * Return the maximum in an array of floats
   *
   * @see {@link Floats#max(float...)}
   * @param array array of floats
   * @return the minimum
   */
  public static float max(float[] array) {
    if (badLength(array.length)) {
      return MINF;
    }

    return Floats.max(array);
  }

  /**
   * Return the maximum in an array of numbers
   *
   * @param array array of numbers
   * @return the maximum
   */
  public static double max(double[] array) {
    if (badLength(array.length)) {
      return MIND;
    }
    return Doubles.max(array);
  }

  /**
   * Return the maximum in an array of numbers, dropping any NaNs
   *
   * @param array array of numbers
   * @return the maximum
   */
  public static double maxDropNaN(double[] array) {
    if (badLength(array.length)) {
      return MIND;
    }
    double max = array[0];
    int ind = 1;
    while (Double.isNaN(max)) {
      max = array[ind];
      ind++;
    }
    if (ind == array.length && Double.isNaN(max)) {
      return Double.NaN;
    }
    for (int i = ind; i < array.length; i++) {
      if (!Double.isNaN(array[i])) {
        max = Math.max(array[i], max);
      }
    }
    return max;
  }

  /**
   * Return the maximum in an array of numbers, dropping any NaNs
   *
   * @param array array of numbers
   * @return the maximum
   */
  public static double minDropNaN(double[] array) {
    if (badLength(array.length)) {
      return MIND;
    }
    double min = array[0];
    int ind = 1;
    while (Double.isNaN(min)) {
      min = array[ind];
      ind++;
    }
    if (ind == array.length && Double.isNaN(min)) {
      return Double.NaN;
    }
    for (int i = ind; i < array.length; i++) {
      if (!Double.isNaN(array[i])) {
        min = Math.min(array[i], min);
      }
    }
    return min;
  }

  /**
   * Return the maximum in an array of bytes
   *
   * @param array array of numbers
   * @return the maximum
   */
  public static byte max(byte[] array) {

    if (badLength(array.length)) {
      return MINB;
    }
    byte min = Byte.MAX_VALUE;
    for (int i = 1; i < array.length; i++) {
      min = (byte) Math.max(array[i], min);
    }
    return min;
  }

  /**
   * Creates an integer array of given size and initializes values to their respective indices
   *
   * @param size size of array
   * @return array of integers with values initialized to their respective indices
   */
  public static int[] arrayOfIndices(int size) {
    return arrayOfIndices(size, 0);
  }

  /**
   * Creates an integer array of given size and initializes values to their respective indices plus
   * the given offset
   *
   * @param size size of array
   * @param start starting offset
   * @return array of integers with values initialized to their respective indices
   */
  public static int[] arrayOfIndices(int size, int start) {
    int[] arr = new int[size];
    for (int i = 0; i < size; i++) {
      arr[i] = i + start;
    }
    return arr;
  }

  /**
   * Creates an integer array of given size and initializes each element with the given value
   *
   * @param size size of array
   * @param initValue initial value of each element
   * @return array of integers initialized to the given value
   */
  public static int[] intArray(int size, int initValue) {
    int[] arr = new int[size];
    Arrays.fill(arr, initValue);
    return arr;
  }

  /**
   * Creates an integer array of given size and initializes each element to 1 except for the first,
   * which is set to zero
   *
   * @param size size of array
   * @return array of integers initialized to the correct values
   */
  public static int[] intArrayStandarddSkips(int size) {
    int[] arr = new int[size];
    Arrays.fill(arr, 1);
    arr[0] = 0;
    return arr;
  }

  /**
   * Creates a long array of given size and initializes each element with the given value
   *
   * @param size size of array
   * @param initValue initial value of each element
   * @return array of longs initialized to the given value
   */
  public static long[] longArray(int size, long initValue) {
    long[] arr = new long[size];
    Arrays.fill(arr, initValue);
    return arr;
  }

  /**
   * Creates an integer array from the contents of a string array
   *
   * @param array array of Strings to be converted
   * @return array of the converted integers
   */
  public static int[] toIntArray(String[] array) {
    int[] arr = new int[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = Integer.parseInt(array[i]);
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into an integer");
      }
    }
    return arr;
  }

  /**
   * Creates an integer array from the contents of a string list
   *
   * @param array array of Strings to be converted
   * @return array of the converted integers
   */
  public static int[] toIntArray(List<String> list) {
    int[] arr = new int[list.size()];
    for (int i = 0; i < list.size(); i++) {
      try {
        arr[i] = Integer.parseInt(list.get(i));
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + list.get(i) + "' into an integer");
      }
    }
    return arr;
  }

  /**
   * Creates an integer array from the contents of a byte array
   *
   * @param array array of Strings to be converted
   * @return array of the converted integers
   */
  public static int[] toIntArray(byte[] array) {
    int[] arr = new int[array.length];
    for (int i = 0; i < array.length; i++) {
      arr[i] = array[i];
    }
    return arr;
  }

  /**
   * Creates an integer array from the contents of a double array
   *
   * @param array array of double to be converted
   * @return array of the converted integers
   */
  public static int[] toIntArray(double[] array) {
    int[] arr = new int[array.length];
    for (int i = 0; i < array.length; i++) {
      arr[i] = (int) array[i];
    }
    return arr;
  }

  /**
   * Creates an array of numbers from the contents of a string array
   *
   * @param array array of Strings to be converted
   * @return array of the converted numbers
   */
  public static double[][] toDoubleArrays(String[][] array, boolean NaNForMissing) {
    double[][] arr = new double[array.length][];
    for (int i = 0; i < array.length; i++) {
      arr[i] = ArrayUtils.toDoubleArray(array[i], NaNForMissing);
    }
    return arr;
  }

  /**
   * Creates an array of numbers from the contents of a string array
   *
   * @param array array of Strings to be converted
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
   * @param array array of floats to be converted
   * @return array of the converted doubles
   */
  public static double[] toDoubleArray(float[] array) {
    double[] arr = new double[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = array[i];
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into a double");
      }
    }
    return arr;
  }

  /**
   * Creates an array of doubles from the contents of a byte array
   *
   * @param array array of floats to be converted
   * @return array of the converted numbers
   */
  public static double[] toDoubleArray(byte[] array) {
    double[] arr = new double[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = array[i];
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into a double");
      }
    }
    return arr;
  }

  /**
   * Creates an array of doubles and copies the contents of an int array into it
   *
   * @param array an array of integers
   * @return an array of doubles copied from an array of integers
   */
  public static double[] toDoubleArray(int[] array) {
    double[] arr = new double[array.length];
    for (int i = 0; i < array.length; i++) {
      arr[i] = array[i];
    }
    return arr;
  }

  /**
   * Creates an array of doubles and copies the contents of a long array into it
   *
   * @param array an array of longs
   * @return an array of doubles copied from an array of longs
   */
  public static double[] toDoubleArray(long[] array) {
    double[] arr = new double[array.length];
    for (int i = 0; i < array.length; i++) {
      arr[i] = array[i];
    }
    return arr;
  }

  /**
   * @param array a float array
   * @return a {@link DoubleStream} containing the elements of the array
   */
  public static DoubleStream toDoubleStream(float[] array) {
    return IntStream.range(0, array.length).mapToDouble(i -> array[i]);
  }

  /**
   * Creates an array of floats from the contents of a double array
   *
   * @param array array of doubles to be converted
   * @return array of the converted numbers
   */
  public static float[] toFloatArray(double[] array) {
    float[] arr = new float[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = (float) array[i];
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into a double");
      }
    }
    return arr;
  }

  /**
   * Creates an array of floats from the contents of a byte array
   *
   * @param array array of bytes to be converted
   * @return array of the converted numbers
   */
  public static float[] toFloatArray(byte[] array) {
    float[] arr = new float[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = array[i];
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into a double");
      }
    }
    return arr;
  }

  /**
   * Creates a double array of given size and initializes each element with the given value
   *
   * @param size size of array
   * @param initValue initial value of each element
   * @return array of numbers initialized to the given value
   */
  public static double[] doubleArray(int size, double initValue) {
    double[] arr = new double[size];
    Arrays.fill(arr, initValue);
    return arr;
  }

  /**
   * Creates a float array of given size and initializes each element with the given value
   *
   * @param size size of array
   * @param initValue initial value of each element
   * @return array of floats initialized to the given value
   */
  public static float[] floatArray(int size, float initValue) {
    float[] arr = new float[size];
    Arrays.fill(arr, initValue);
    return arr;
  }

  /**
   * Creates a byte array of given size and initializes each element with the given value
   *
   * @param size size of array
   * @param initValue initial value of each element
   * @return array of bytes initialized to the given value
   */
  public static byte[] byteArray(int size, byte initValue) {
    byte[] arr = new byte[size];
    Arrays.fill(arr, initValue);
    return arr;
  }

  /**
   * Create a boolean array and make all states an initial value
   *
   * @param int size of array
   * @return the boolean array
   */
  public static boolean[] booleanArray(int size, boolean initValue) {
    boolean[] array = new boolean[size];

    // Default value is false
    if (initValue) {
      Arrays.fill(array, true);
    }

    return array;
  }

  public enum BYTE_DECODE_FORMAT {
    /**
     * String will be converted to upper case
     */
    UPPER_CASE,
    /**
     * String will be converted to lower case
     */
    LOWER_CASE,
    /**
     * String will be left as is
     */
    AS_IS
  }

  public static String[] decodeByteArray(byte[] b, BYTE_DECODE_FORMAT format, Logger log) {
    return decodeByteArray(b, ext.UTF_8, format, log);
  }

  public static String[] decodeByteArray(byte[] b, Logger log) {
    return decodeByteArray(b, ext.UTF_8, BYTE_DECODE_FORMAT.AS_IS, log);
  }

  /**
   * @param b each entry will be converted to a string
   * @param charsetName
   * @param format
   * @param log
   * @return
   */
  public static String[] decodeByteArray(byte[] b, String charsetName, BYTE_DECODE_FORMAT format,
                                         Logger log) {
    String[] s = new String[b.length];
    for (int i = 0; i < s.length; i++) {
      if ((i + 1) % 2000000 == 0) {
        log.reportTimeInfo((i + 1) + " entries converted");
      }
      try {
        s[i] = new String(new byte[] {b[i]}, charsetName).toUpperCase();
      } catch (UnsupportedEncodingException e) {
        log.reportError("Could not convert reference byte " + b[i] + " to string with charsetName"
                        + charsetName);
        e.printStackTrace();
      }
    }
    return s;
  }

  /**
   * Creates a byte array from the contents of an integer array
   *
   * @param array array of ints to be converted
   * @return array of the converted bytes
   */
  public static byte[] toByteArray(int[] array) {
    byte[] arr = new byte[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = (byte) array[i];
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into a byte");
      }
    }
    return arr;
  }

  public static byte[] toByteArray(String[] array) {
    byte[] arr = new byte[array.length];
    for (int i = 0; i < array.length; i++) {
      try {
        arr[i] = Byte.valueOf(array[i]);
      } catch (NumberFormatException nfe) {
        System.err.println("Error - failed to convert '" + array[i] + "' into a byte");
      }
    }
    return arr;
  }

  /**
   * Creates a byte array from the contents of an integer array
   *
   * @param array array of floats to be converted
   * @return array of the converted bytes
   */
  public static byte[] toByteArray(float[] array) {
    byte[] arr = new byte[array.length];
    for (int i = 0; i < array.length; i++) {
      arr[i] = (byte) array[i];
    }
    return arr;
  }

  /**
   * Creates a String array of given size and initializes each element with the given value
   *
   * @param size size of array
   * @param initValue initial value of each element
   * @return array of Strings initialized to the given value
   */
  public static String[] stringArray(int size, String initValue) {
    String[] array = new String[size];
    for (int i = 0; i < size; i++) {
      array[i] = initValue;
    }
    return array;
  }

  /**
   * Creates a String array of given size and initializes each element with the "base value"
   * +(index+1)
   *
   * @param size size of array
   * @param base base value of each element
   * @return array of Strings initialized to the given indexed values
   */
  public static String[] stringArraySequence(int size, String base) {
    return stringArraySequence(size, base, "");
  }

  /**
   * Creates a String array of given size and initializes each element with the
   * "prefix"+(index+1)+"suffix"
   *
   * @param size size of array
   * @param prefix prefix for each element
   * @param suffix suffix for each element
   * @return array of Strings initialized to the given indexed values
   */
  public static String[] stringArraySequence(int size, String prefix, String suffix) {
    String[] array = new String[size];
    for (int i = 0; i < size; i++) {
      array[i] = prefix + (i + 1) + suffix;
    }
    return array;
  }

  /**
   * Creates an array of blank Strings
   *
   * @param size size of array
   * @return array of blank Strings
   */
  public static String[] stringArray(int size) {
    return stringArray(size, "");
  }

  /**
   * Creates an integer array of given size and initializes values by randomly sampling zero to size
   * without replacement
   *
   * @param size size of array
   * @return array of random indices for the given size
   */
  public static int[] random(int size) {
    return random(size, size, false);
  }

  /**
   * Creates an integer array of a set number of selections, and initializes values by randomly
   * sampling zero to size with replacement
   *
   * @param size range of random values
   * @param numSelections number of random values to select
   * @return array of the specified number of selections from a distribution of the given size
   */
  public static int[] randomWithReplacement(int size, int numSelections) {
    return random(size, numSelections, true);
  }

  /**
   * Creates an integer array of a set number of selections, and initializes values by randomly
   * sampling zero to size without replacement
   *
   * @param size range of random values
   * @param numSelections number of random values to select
   * @return array of the specified number of selections from a distribution of the given size
   */
  public static int[] random(int size, int numSelections) {
    return random(size, numSelections, false);
  }

  /**
   * Helper method for generating random arrays of randomized integer indices. Allows configuration
   * of distribution range, number of draws, and whether values can be repeated or not.
   */
  private static int[] random(int distSize, int numSelections, boolean withReplacement) {
    int[] selections;
    Random r = new Random();
    if (withReplacement) {
      // We can just take a random number from the distribution for each position
      selections = new int[numSelections];
      for (int i = 0; i < numSelections; i++) {
        selections[i] = r.nextInt(selections.length);
      }
    } else {
      // Without replacement we start with an array of indices and do a fisher-yates shuffle
      selections = arrayOfIndices(distSize);
      for (int i = 0; i < selections.length; i++) {
        int idx = r.nextInt(selections.length);
        int tmp = selections[i];
        selections[i] = selections[idx];
        selections[idx] = tmp;
      }
    }

    return selections;
  }

  public static double sum(double[] array, boolean[] incl) {
    double sum = 0;
    for (int i = 0; i < array.length; i++) {
      if (incl[i]) sum += array[i];
    }
    return sum;
  }

  /**
   * Calculates the sum of an array
   * <p>
   * TODO replace with Streams in Java 8
   * </p>
   *
   * @param array an array of numbers
   * @return sum of the array
   */
  public static double sum(double[] array) {
    double sum = 0;

    for (double element : array) {
      // if (Double.isNaN(array[i])) {
      // System.err.println("Are you sure you want to sum: "+array[i]);
      // }
      sum += element;
    }

    return sum;
  }

  /**
   * Calculates the sum of all values in an array specified by the indices in a second array
   *
   * @param array an array of numbers
   * @param order an array of indices to sum
   * @return sum of the array at the given indices
   */
  public static double sumInOrder(double[] array, int[] order) {
    double sum = 0;

    for (int i = 0; i < array.length; i++) {
      sum += array[order[i]];
    }

    return sum;
  }

  /**
   * Calculates the sum of an array
   *
   * @param array an array of numbers
   * @return sum of the array
   */
  public static double sumExactUsingBigDecimal(double[] array) {
    BigDecimal sum = BigDecimal.ZERO;

    for (double element : array) {
      sum = sum.add(BigDecimal.valueOf(element));
    }

    return sum.doubleValue();
  }

  /**
   * Calculates the sum of an array
   * <p>
   * TODO replace with Streams in Java 8
   * </p>
   *
   * @param array an array of numbers
   * @return sum of the array
   */
  public static float sum(float[] array) {
    float sum = 0;

    for (float element : array) {
      sum += element;
    }

    return sum;
  }

  /**
   * Calculates the sum of an array
   * <p>
   * TODO replace with Streams in Java 8
   * </p>
   *
   * @param array an array of numbers
   * @return sum of the array
   */
  public static Float sum(Float[] array) {
    Float sum = Float.valueOf(0);

    for (Float element : array) {
      sum += element;
    }

    return sum;
  }

  /**
   * Calculates the sum of an array
   * <p>
   * TODO replace with Streams in Java 8
   * </p>
   *
   * @param array an array of integers
   * @return sum of the array
   */
  public static int sum(int[] array) {
    int sum = 0;

    for (int element : array) {
      sum += element;
    }

    return sum;
  }

  /**
   * Calculates the sum of an array
   * <p>
   * TODO replace with Streams in Java 8
   * </p>
   *
   * @param array an array of integers
   * @return sum of the array
   */
  public static int sum(byte[] array) {
    int sum = 0;

    for (byte element : array) {
      sum += element;
    }

    return sum;
  }

  /**
   * Calculates the sum of an array
   * <p>
   * TODO replace with Streams in Java 8
   * </p>
   *
   * @param array an array of long
   * @return sum of the array
   */
  public static long sum(long[] array) {
    long sum = 0;

    for (long element : array) {
      sum += element;
    }

    return sum;
  }

  public static <T extends Number> double sum(Collection<T> collection) {
    double s = 0.0;
    for (T val : collection) {
      s += val.doubleValue();
    }
    return s;
  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static double mean(double[] array) {
    return sum(array) / array.length;
  }

  /**
   * @param array
   * @param factor multiply every value in the array by this number
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
   * Note: This method simply calls {@link Maths#movingAverageForward(int, Collection, boolean)}
   * with a List view of the array and converts the result to a double[], if array results are not
   * required, use Maths#movingAverageForward(int, Collection, boolean)}
   * 
   * @param n number of points for the moving average
   * @param array
   * @param skipNaN If true, every moving average will be composed of n points, or NaN; defaults to
   *          removing Nan for the average
   * @return an array of moving averages with moving average of n sequential points
   */
  public static double[] movingAverageForward(int n, double[] array, boolean skipNaN) {
    return Doubles.toArray(Maths.movingAverageForward(n, Doubles.asList(array), skipNaN));
  }

  /**
   * @param n number of points for the moving average
   * @param array
   * @return an array of moving averages with moving average of n sequential points
   */
  public static double[] movingMedianForward(int n, double[] array) {
    double[] ma = new double[array.length];
    // double[] a = new double[n];
    ArrayList<Double> tmp = new ArrayList<>();
    for (int i = 0; i < array.length; i++) {
      tmp.add(array[i]);
      if (i >= n) {
        ma[i] = median(Doubles.toArray(tmp));
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
   * @param minForce the minimum of the returned array, values scaled to have this value as min
   * @return
   */
  public static double[] scaleMinTo(double[] array, final double minForce) {
    double min = ArrayUtils.min(ArrayUtils.removeNaN(array));

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
   * @param array scale the values of this array between minForce and minForce +1
   * @param minForce the minimum value (max will be this plus 1)
   * @return
   */
  public static double[] scale(final double[] array, final double minForce) {
    double max = ArrayUtils.max(array);
    double min = ArrayUtils.min(array);
    double[] scaled = new double[array.length];
    for (int i = 0; i < array.length; i++) {
      scaled[i] = (array[i] - min) / (max - min);
      scaled[i] += minForce;
    }
    return scaled;
  }

  /**
   * Calculates the mean of an Iterable
   * 
   * @param values values to calculate mean of
   * @param ignoreNaN true to skip NaN values
   * @return mean of the values
   */
  public static double mean(Iterable<? extends Number> values, boolean ignoreNaN) {
    double sum = 0.0;
    int count = 0;
    for (Number val : values) {
      double dblVal = val.doubleValue();
      if (!Double.isNaN(dblVal) || !ignoreNaN) {
        sum += dblVal;
        count++;
      }
    }

    if (count == 0) {
      return Double.NaN;
    }

    return sum / count;

  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static double mean(double[] array, boolean ignoreNaN) {
    double sum;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i]) || !ignoreNaN) {
        sum += array[i];
        count++;
      }
    }

    if (count == 0) {
      return Double.NaN;
    }

    return sum / count;
  }

  /**
   * Calculates the mean of a float array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static float mean(float[] array, boolean ignoreNaN) {
    float sum;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i]) || !ignoreNaN) {
        sum += array[i];
        count++;
      }
    }

    if (count == 0) {
      return Float.NaN;
    }

    return sum / count;
  }

  /**
   * Calculates the mean distance from a value of all values in a float array
   *
   * @param array an array of numbers
   * @param val value to measure distance from
   * @return mean of the array
   */
  public static float meanDist(float[] array, float val, boolean ignoreNaN) {
    float sum;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i]) || !ignoreNaN) {
        sum += Math.abs(array[i] - val);
        count++;
      }
    }

    if (count == 0) {
      return Float.NaN;
    }

    return sum / count;
  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static double mean(int[] array) {
    return (double) sum(array) / array.length;
  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static double mean(byte[] array) {
    return (double) sum(array) / array.length;
  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static double meanIf(double[] array, double[] filter, double filterValue,
                              boolean ignoreNaN) {
    double sum;
    int count;

    if (array.length != filter.length) {
      System.err.println("Error - filter size does not match array size");
      return Double.NEGATIVE_INFINITY;
    }

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (filter[i] == filterValue && (!ignoreNaN || !Double.isNaN(array[i]))) {
        sum += array[i];
        count++;
      }
    }

    return sum / count;
  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static float mean(float[] array) {
    return sum(array) / array.length;
  }

  /**
   * Calculates the mean of an array
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static Float mean(Float[] array) {
    return sum(array) / array.length;
  }

  public static <T extends Number> double mean(Collection<T> collection) {
    if (collection.isEmpty()) {
      return 0.0;
    }
    return sum(collection) / collection.size();
  }

  /**
   * Calculates the variance of an array
   *
   * @param array an array of numbers
   * @return variance of the array
   */
  public static double variance(int[] array) {
    double avg = mean(array);
    double sum = 0;

    for (int element : array) {
      sum += Math.pow(avg - element, 2);
    }

    return sum / (array.length - 1);
  }

  /**
   * Calculates the variance of an array
   *
   * @param array an array of numbers
   * @return variance of the array
   */
  public static double variance(double[] array) {
    double avg = mean(array);
    double sum = 0;

    for (double element : array) {
      sum += Math.pow(avg - element, 2);
    }

    return sum / (array.length - 1);
  }

  public static double variance(Collection<? extends Number> values) {
    double avg = mean(values);
    double sum = 0;

    for (Number element : values) {
      double dblElement = element.doubleValue();
      sum += Math.pow(avg - dblElement, 2);
    }
    return sum / (values.size() - 1);
  }

  /**
   * Calculates the variance of an array, dropping the NaNs in the array.
   *
   * @param array an array of numbers
   * @return variance of the array
   */
  public static double varianceDropNaN(double[] array) {
    double avg = mean(array, true);
    double sum = 0;
    double cnt = 0;

    for (double element : array) {
      if (Double.isNaN(element)) continue;
      sum += Math.pow(avg - element, 2);
      cnt++;
    }

    return sum / (cnt - 1);
  }

  /**
   * Calculates the mean of an array if the sum is already known
   *
   * @param array an array of numbers
   * @return mean of the array
   */
  public static double mean(double[] array, double sum) {
    return sum / array.length;
  }

  /**
   * Calculates the sum of squares of an array if the mean is already known
   *
   * @param array an array of numbers
   * @param avg precomputed average of the array
   * @return variance of the array
   */
  public static double sumSq(double[] array, double avg) {
    double sum = 0;
    for (double element : array) {
      sum += Math.pow(avg - element, 2);
    }
    return sum;
  }

  /**
   * Calculates the variance of an array if the sum of squares is known
   *
   * @param array an array of numbers
   * @param avg precomputed average of the array
   * @return variance of the array
   */
  public static double variance(double[] array, double sumSq) {
    return sumSq / (array.length - 1);
  }

  /**
   * Calculates the variance of an array
   *
   * @param array an array of numbers
   * @return variance of the array
   */
  public static float variance(float[] array) {
    double sum, avg; // allows for larger arrays

    sum = 0;
    for (float element : array) {
      sum += element;
    }
    avg = sum / array.length;

    for (float element : array) {
      sum += Math.pow(avg - element, 2);
    }

    return (float) (sum / (array.length - 1));
  }

  /**
   * Calculates the standard deviation of an array
   *
   * @param array an array of numbers
   * @return standard deviation of the array
   */
  public static double stdev(int[] array) {
    return Math.sqrt(variance(array));
  }

  /**
   * Calculates the standard deviation of an array
   *
   * @param array an array of numbers
   * @return standard deviation of the array
   */
  public static double stdev(double[] array) {
    return Math.sqrt(varianceDropNaN(array));
  }

  /**
   * Calculates the standard deviation of a Collection
   * 
   * @param values
   * @return standard deviation of values
   */
  public static double stdev(Collection<? extends Number> values) {
    return Math.sqrt(variance(values));
  }

  /**
   * Calculates the standard deviation of an array
   *
   * @param array an array of numbers
   * @param removeNaN remove any value that is not a number
   * @return standard deviation of the [filtered] array
   */
  public static float stdev(float[] array, boolean removeNaN) {
    double sum, avg;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i]) || !removeNaN) {
        sum += array[i];
        count++;
      }
    }
    avg = (sum / count);

    sum = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i])) {
        sum += Math.pow(avg - array[i], 2);
      }
    }

    return (float) Math.sqrt(sum / (count - 1));
  }

  /**
   * Calculates the standard deviation of an array
   *
   * @param array an array of numbers
   * @param removeNaN remove any value that is not a number
   * @return standard deviation of the [filtered] array
   */
  public static double stdev(Double[] array, boolean removeNaN) {
    double sum, avg;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i]) || !removeNaN) {
        sum += array[i];
        count++;
      }
    }
    avg = (sum / count);

    sum = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i])) {
        sum += Math.pow(avg - array[i], 2);
      }
    }

    return Math.sqrt(sum / (count - 1));
  }

  /**
   * Calculates the standard deviation of an array
   *
   * @param array an array of numbers
   * @param removeNaN remove any value that is not a number
   * @return standard deviation of the [filtered] array
   */
  public static Float stdev(Float[] array, boolean removeNaN) {
    double sum, avg;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i]) || !removeNaN) {
        sum += array[i];
        count++;
      }
    }
    avg = (sum / count);

    sum = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Float.isNaN(array[i])) {
        sum += Math.pow(avg - array[i], 2);
      }
    }

    return (float) Math.sqrt(sum / (count - 1));
  }

  /**
   * Calculates the standard deviation of an array
   *
   * @param array an array of numbers
   * @param removeNaN remove any value that is not a number
   * @return standard deviation of the [filtered] array
   */
  public static double stdev(double[] array, boolean removeNaN) {
    double sum, avg;
    int count;

    sum = 0;
    count = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i]) || !removeNaN) {
        sum += array[i];
        count++;
      }
    }
    avg = (sum / count);

    sum = 0;
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i])) {
        sum += Math.pow(avg - array[i], 2);
      }
    }

    return (float) Math.sqrt(sum / (count - 1));
  }

  /**
   * Calculates the standard deviation of an Iterable
   *
   * @param values
   * @param removeNaN remove any value that is not a number
   * @return standard deviation of the [filtered] values
   */
  public static double stdev(Iterable<? extends Number> values, boolean removeNaN) {

    double sum = 0;
    int count = 0;
    for (Number val : values) {
      double dblVal = val.doubleValue();
      if (!Double.isNaN(dblVal) || !removeNaN) {
        sum += dblVal;
        count++;
      }
    }
    if (count == 0) return Double.NaN;
    double avg = (sum / count);

    sum = 0;
    for (Number val : values) {
      double dblVal = val.doubleValue();
      if (!Double.isNaN(dblVal)) {
        sum += Math.pow(avg - dblVal, 2);
      }
    }

    return Math.sqrt(sum / (count - 1));
  }

  /**
   * Normalizes (calculates z-scores) for an array of numbers
   *
   * @param array an array of numbers
   * @return array of z-scores
   */
  public static double[] normalize(double[] array) {
    double[] newData = new double[array.length];
    double mean = ArrayUtils.mean(array);
    double stdev = ArrayUtils.stdev(array);

    for (int i = 0; i < newData.length; i++) {
      newData[i] = (array[i] - mean) / stdev;
    }

    return newData;
  }

  /**
   * Standardizes (calculates z-scores) for an array of numbers
   *
   * @param array an array of numbers
   * @return array of z-scores
   */
  public static float[] normalize(float[] array) {
    float[] newData = new float[array.length];
    float mean = ArrayUtils.mean(array, true);
    float stdev = ArrayUtils.stdev(array, false);

    for (int i = 0; i < newData.length; i++) {
      newData[i] = (array[i] - mean) / stdev;
    }

    return newData;
  }

  /**
   * Normalizes (calculates z-scores) for an array of numbers using separate standard deviations for
   * positive and negative numbers
   *
   * @param array an array of numbers
   * @return array of sign-specific z-scores
   */
  public static double[] normalizeSigned(double[] array) {
    double[] newData;
    double mean = 0;
    double stdevPositive, stdevNegative;
    DoubleVector positives, negatives;

    negatives = new DoubleVector();
    positives = new DoubleVector();
    for (double element : array) {
      if (element < 0) {
        negatives.add(element);
        negatives.add(-1 * element);
      } else {
        positives.add(element);
        positives.add(-1 * element);
      }
    }

    stdevNegative = ArrayUtils.stdev(Doubles.toArray(negatives));
    stdevPositive = ArrayUtils.stdev(Doubles.toArray(positives));

    newData = new double[array.length];
    for (int i = 0; i < newData.length; i++) {
      if (array[i] < 0) {
        newData[i] = (array[i] - mean) / stdevNegative;
      } else {
        newData[i] = (array[i] - mean) / stdevPositive;
      }
    }

    return newData;
  }

  /**
   * Returns the quantiles of an array
   *
   * @param array an array of numbers
   * @return array of quantiles
   */
  public static double[] quantiles(double[] array) {
    double[] quantiles;
    int[] order;

    order = Sort.getSortedIndices(array);

    quantiles = new double[array.length];
    for (int i = 0; i < quantiles.length; i++) {
      quantiles[order[i]] = ((double) i + 1) / ((double) quantiles.length + 1);
    }

    return quantiles;
  }

  /**
   * Returns the kurtosis of an array
   *
   * @param array
   * @return
   */
  public static double skewness(double[] array) {
    double skew = -1;
    double mean = ArrayUtils.mean(array);
    double sd = ArrayUtils.stdev(array);
    double m3;
    double n = array.length;

    m3 = 0;
    for (int i = 0; i < n; i++) {
      m3 += Math.pow((array[i] - mean) / sd, 3);
    }
    skew = m3 * n / (n - 1) / (n - 2);

    return skew;
  }

  /**
   * Returns the kurtosis of an array
   *
   * @param array
   * @return
   */
  public static double kurtosis(double[] array) {
    double kurt = -1;
    double mean = ArrayUtils.mean(array);
    double m2, s, m4s;
    double n = array.length;

    m2 = 0;
    for (int i = 0; i < n; i++) {
      m2 += Math.pow(array[i] - mean, 2);
    }
    m2 /= (n - 1);
    s = Math.sqrt(m2);

    m4s = 0;
    for (double element : array) {
      m4s += Math.pow((element - mean) / s, 4);
    }

    kurt = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * m4s
           - 3 * Math.pow(n - 1, 2) / ((n - 2) * (n - 3));

    return kurt;
  }

  /**
   * Inverse-normalizes an array of numbers
   *
   * @param array an array of numbers
   * @return array of inverse-normalized values
   */
  public static double[] inverseNormalize(double[] array) {
    double[] probits, quantiles;

    quantiles = quantiles(array);
    probits = new double[array.length];
    for (int i = 0; i < probits.length; i++) {
      if (quantiles[i] < 0.5) {
        probits[i] = ProbDist.NormDistReverse(quantiles[i] * 2) * -1;
      } else {
        quantiles[i] = 1 - quantiles[i];
        probits[i] = ProbDist.NormDistReverse(quantiles[i] * 2) * 1;
      }
    }

    return probits;
  }

  /**
   * Inverse-normalizes an array of numbers
   *
   * @param array an array of numbers
   * @return array of inverse-normalized values
   */
  public static double[] inverseTdist(double[] array, int df) {
    double[] newValues, quantiles;

    quantiles = quantiles(array);
    newValues = new double[array.length];
    for (int i = 0; i < newValues.length; i++) {
      if (quantiles[i] < 0.5) {
        newValues[i] = ProbDist.TDistReverse(quantiles[i] * 2, df) * -1;
      } else {
        quantiles[i] = 1 - quantiles[i];
        newValues[i] = ProbDist.TDistReverse(quantiles[i] * 2, df) * 1;
      }
    }

    return newValues;
  }

  /**
   * Returns the bootstrapped median of an array
   *
   * @param array an array of numbers
   * @param numReps number of replicates to perform
   * @return array of the median, the 2.5 and the 97.5 bootstrap percentile intervals
   */
  public static double[] bootstrap(double[] array, int numReps, boolean verbose) {
    double[] results = new double[3];
    double[] replicates = new double[numReps];
    int progress = 0;
    double sum;

    if (verbose) {
      System.out.print("  Bootstrap rep 0");
    }
    for (int i = 0; i < numReps; i++) {
      if (verbose && i + 1 == progress + numReps / 10) {
        for (int j = 0; j < (progress + "").length(); j++) {
          System.out.print("\b");
        }
        progress += numReps / 10;
        System.out.print(progress);
      }

      sum = 0;
      for (double element : array) {
        sum += array[(int) (Math.random() * array.length)];
      }
      replicates[i] = sum / array.length;
    }
    if (verbose) {
      System.out.println();
    }

    Arrays.sort(replicates);
    results[0] = replicates[(int) (numReps * 0.5)];
    results[1] = replicates[(int) (numReps * 0.025)];
    results[2] = replicates[(int) (numReps * 0.975)];

    return results;
  }

  /**
   * Determines the specified exclusive quantile of an array of numbers<br />
   * Returns a number guaranteed to be a member of the given array.<br />
   * This function matches Excel's QUARTILE.EXC function.
   *
   * @param array an array of numbers
   * @param q exclusive quantile to be determined
   * @return specified exclusive quantile of the array
   */
  public static double quantExclusive(int[] array, double q) {
    if (array.length == 0) {
      return Double.NaN;
    }

    int[] keys = Sort.getSortedIndices(array);

    try {
      if (q > 1 || q < 0) {
        return (0);
      } else {
        double index = (array.length + 1) * q;
        if (index - (int) index == 0) {
          return array[keys[(int) index - 1]];
        } else {
          return q * array[keys[(int) Math.floor(index) - 1]]
                 + (1 - q) * array[keys[(int) Math.ceil(index) - 1]];
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
   * @param array an array of numbers
   * @param q exclusive quantile to be determined
   * @return specified exclusive quantile of the array
   */
  public static double quantExclusive(double[] array, double q) {
    if (array.length == 0) {
      return Double.NaN;
    }

    int[] keys = Sort.getSortedIndices(array);

    try {
      if (q > 1 || q < 0) {
        return (0);
      } else {
        double index = (array.length + 1) * q;
        if (index - (int) index == 0) {
          return array[keys[(int) index - 1]];
        } else {
          return q * array[keys[(int) Math.floor(index) - 1]]
                 + (1 - q) * array[keys[(int) Math.ceil(index) - 1]];
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
   * @param array an array of numbers
   * @param q exclusive quantile to be determined
   * @return specified exclusive quantile of the array
   */
  public static float quantExclusive(float[] array, float q) {
    if (array.length == 0) {
      return Float.NaN;
    }

    int[] keys = Sort.getSortedIndices(array);

    try {
      if (q > 1 || q < 0) {
        return 0f;
      } else {
        double index = (array.length + 1) * q;
        if (index - (int) index == 0) {
          return array[keys[(int) index - 1]];
        } else {
          return q * array[keys[(int) Math.floor(index) - 1]]
                 + (1 - q) * array[keys[(int) Math.ceil(index) - 1]];
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
   * @param array an array of numbers
   * @param q quantile to be determined
   * @return specified quantile of the array
   */
  public static int quantWithExtremeForTie(int[] array, double q) {
    int keys[] = Sort.getSortedIndices(array);

    try {
      if (q > 1 || q < 0) {
        return (0);
      } else {
        double index = (array.length + 1) * q;
        if (index - (int) index == 0) {
          return array[keys[(int) index - 1]];
        } else if (q < 0.5) {
          return array[keys[(int) Math.floor(index) - 1]];
        } else {
          return array[keys[(int) Math.ceil(index) - 1]];
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
   * @param array an array of numbers
   * @param q quantiles to be determined
   * @return specified quantiles of the array
   */
  public static float[] quants(float[] array, double[] qs) {
    int keys[] = Sort.getSortedIndices(array);
    float[] quantiles;

    quantiles = new float[qs.length];
    for (int i = 0; i < quantiles.length; i++) {
      try {
        if (qs[i] > 1 || qs[i] < 0) {
          quantiles[i] = -1;
        } else {
          double index = (array.length + 1) * qs[i];
          if (index - (int) index == 0) {
            quantiles[i] = array[keys[(int) index - 1]];
          } else {
            quantiles[i] = (float) (qs[i] * array[keys[(int) Math.floor(index) - 1]]
                                    + (1 - qs[i]) * array[keys[(int) Math.ceil(index) - 1]]);
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
   * @param array an array of numbers
   * @param q quantiles to be determined
   * @return specified quantiles of the array
   */
  public static double[] quants(double[] array, double[] qs) {
    int keys[] = Sort.getSortedIndices(array);
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
            quantiles[i] = (float) (qs[i] * array[keys[(int) Math.floor(index) - 1]]
                                    + (1 - qs[i]) * array[keys[(int) Math.ceil(index) - 1]]);
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
   * @param array an array of numbers
   * @param q quantiles to be determined
   * @return specified quantiles of the array
   */
  public static double[] quantsExclusive(double[] array, double[] qs) {
    int keys[] = Sort.getSortedIndices(array);
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
            quantiles[i] = qs[i] * array[keys[(int) Math.floor(index) - 1]]
                           + (1 - qs[i]) * array[keys[(int) Math.ceil(index) - 1]];
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
   * Determines the median absolute difference of an array of floats
   *
   * @param array an array of numbers
   * @return mad of the array
   */
  public static double mad(float[] array) {
    double median = median(array);
    return mad(array, median);
  }

  public static double mad(float[] values, double median) {
    return mad(toDoubleStream(values), median, values.length);
  }

  /**
   * Determines the median absolute difference of an array of floats
   *
   * @param array an array of numbers
   * @param constant factor to multiply result by
   * @return mad of the array
   */
  public static double madFactor(float[] array, double constant) {
    return mad(array) * constant;
  }

  /**
   * Determines the median absolute difference of an array of double
   *
   * @param array an array of numbers
   * @return mad of the array
   */
  public static double mad(double[] array) {
    double median = median(array);
    return mad(array, median);
  }

  /**
   * As {@link #mad(double[])} but we know the initial array is already sorted
   */
  public static double madSorted(double[] array) {
    double median = medianSorted(array);
    return mad(array, median);
  }

  /**
   * Determines the median absolute difference of an array of double
   *
   * @param array an array of numbers
   * @param constant factor to multiply result by
   * @return mad of the array
   */
  public static double madFactor(double[] array, double constant) {
    double median = median(array);
    return madFactor(array, constant, median);
  }

  /**
   * Determines the median absolute difference of an array of double with a known median
   *
   * @param array an array of numbers
   * @param constant factor to multiply result by
   * @param median the median of the given array
   * @return mad of the array
   */
  public static double madFactor(double[] array, double constant, double median) {
    return mad(array, median) * constant;
  }

  public static double mad(double[] values, double median) {
    return mad(Arrays.stream(values), median, values.length);
  }

  /**
   * Determines the median absolute difference of a Collection of Numbers
   * 
   * @param values
   * @return
   */
  public static double mad(Collection<? extends Number> values) {
    double median = median(values);
    return mad(values, median);
  }

  /**
   * Determines the median absolute difference of a Collection of Numbers with a known median
   *
   * @param values values to find MAD of
   * @param median the median of the given values
   * @return MAD of the values
   */
  public static double mad(Collection<? extends Number> values, double median) {
    return mad(values.stream().mapToDouble(Number::doubleValue), median, values.size());
  }

  /**
   * Determines the median absolute difference of a {@link DoubleStream} with a known median and
   * size
   *
   * @param values values to find MAD of
   * @param median the median of the given values
   * @param size the size of values
   * @return MAD of the values
   */
  public static double mad(DoubleStream values, double median, int size) {
    return median(values.map(v -> Math.abs(v - median)), size);
  }

  /**
   * Determines the median of an array of numbers
   *
   * @param array an array of numbers
   * @return median of the array
   */
  public static double median(int[] array) {
    return median(Arrays.stream(array).asDoubleStream(), array.length);
  }

  /**
   * Determines the median of an array of numbers
   *
   * @param array an array of numbers
   * @return median of the array
   */
  public static double median(double[] array) {
    return median(Arrays.stream(array), array.length);
  }

  /**
   * Determines the median of an array of numbers
   *
   * @param array an array of numbers
   * @param dropNaN true to exclude NaN values from the array
   * @return median of the array
   */
  public static double median(double[] array, boolean dropNaN) {
    return dropNaN ? median(removeNonFinites(array)) : median(array);
  }

  /**
   * Determines the median of an array of numbers
   *
   * @param array an array of numbers
   * @return median of the array
   */
  public static double median(float[] array) {
    return median(toDoubleStream(array), array.length);
  }

  /**
   * @param list a sorted {@link List} of {@link Number}s
   * @return median of collection
   */
  public static double medianSorted(List<? extends Number> list) {
    int size = list.size();
    if (size <= 0) return Double.NaN;
    final int midpoint = size / 2;
    double median;
    if (size % 2 == 0) {
      median = (list.get(midpoint - 1).doubleValue() + list.get(midpoint).doubleValue()) / 2.0;
    } else {
      median = list.get(midpoint).doubleValue();
    }
    return median;
  }

  /**
   * @param collection a sorted {@link Collection} of {@link Number}s
   * @return median of collection
   */
  public static double medianSorted(Collection<? extends Number> collection) {
    if (collection instanceof List) return medianSorted((List<? extends Number>) collection);
    return medianSorted(collection.stream(), collection.size());
  }

  /**
   * @param collection a {@link Collections} of {@link Number}s
   * @return median of collection
   */
  public static double median(Collection<? extends Number> collection) {
    return median(collection.stream(), collection.size());
  }

  /**
   * @param stream a {@link Stream} of {@link Number}s with known size
   * @param size size of stream
   * @return median of stream
   */
  public static double median(Stream<? extends Number> stream, int size) {
    return median(stream.mapToDouble(Number::doubleValue), size);
  }

  /**
   * @param stream a {@link DoubleStream} with known size
   * @param size size of stream
   * @return median of stream
   */
  public static double median(DoubleStream stream, int size) {
    return medianSorted(stream.sorted(), size);
  }

  /**
   * @param array a sorted array of doubles
   * @return median of array
   */
  public static double medianSorted(double[] array) {
    return medianSorted(Doubles.asList(array));
  }

  /**
   * @param stream a sorted {@link Stream} of {@link Number}s with known size
   * @param size size of stream
   * @return median of stream
   */
  public static double medianSorted(Stream<? extends Number> stream, int size) {
    return medianSorted(stream.mapToDouble(Number::doubleValue), size);
  }

  /**
   * @param stream a sorted {@link DoubleStream} with known size
   * @param size size of stream
   * @return median of stream
   */
  public static double medianSorted(DoubleStream stream, int size) {
    if (size <= 0) return Double.NaN;
    final int midpoint = size / 2;
    double median;
    if (size % 2 == 0) {
      median = stream.skip(midpoint - 1).limit(2).summaryStatistics().getAverage();
    } else {
      median = stream.skip(midpoint).findFirst().getAsDouble();
    }
    return median;
  }

  /**
   * Prints an array of objects separated by a tab
   *
   * @param array an array of objects
   * @return String of printed objects
   */
  public static String toStr(String[] array) {
    return toStr(array, null, "\t", null);
  }

  /**
   * @see #toStr(Collection, String, String)
   */
  public static String toStr(Collection<?> c) {
    return toStr(c, "\t");
  }

  /**
   * @see #toStr(Collection, String, String)
   */
  public static String toStr(Collection<?> c, String delim) {
    return toStr(c, delim, null);
  }

  /**
   * Returns a Collection of Objects as a String separated by the specified delimiter
   *
   * @param collection a Collection of Strings
   * @param delimiter String delimiter
   * @param nullValue value to use in place of nulls
   * @return String of printed objects
   */
  public static String toStr(final Collection<?> collection, String delimiter, String nullValue) {
    List<String> cleanList = Lists.newArrayListWithCapacity(collection.size());
    boolean commaDelimited;

    commaDelimited = delimiter.equals(",");
    for (Object element : collection) {
      String elementString = element == null ? nullValue : element.toString();
      if (commaDelimited && elementString.contains(",")) {
        elementString = "\"" + elementString + "\"";
      }
      cleanList.add(elementString);
    }

    return Joiner.on(delimiter).join(cleanList);
  }

  /**
   * @param a
   * @param b return the distance array from this value
   * @return
   */
  public static double[] distFrom(double[] a, double b) {
    return abs(minus(a, b));
  }

  /**
   * @param a
   * @param minus subtract this value from every entry in the array
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
   *
   * @param enumValue An enum class, must be passed as ENUM.class (not the enum itself, but
   *          "&lt;enum&gt;.class")
   * @param delimiter
   * @return
   */
  public static <T extends Enum<?>> String toStr(Class<T> enumValue, String delimiter) {
    T[] values = enumValue.getEnumConstants();
    String[] arr = new String[values.length];
    for (int i = 0; i < values.length; i++) {
      arr[i] = values[i].toString();
    }
    return ArrayUtils.toStr(arr, null, delimiter, null);
  }

  /**
   * Prints an array of objects separated by the specified delimiter
   *
   * @param array an array of objects
   * @param delimiter String delimiter
   * @return String of printed objects
   */
  public static String toStr(String[] array, String delimiter) {
    return toStr(array, null, delimiter, null);
  }

  /**
   * Prints an array of objects separated by the specified delimiter
   *
   * @param array an array of objects
   * @param delimiter String delimiter
   * @return String of printed objects
   */
  public static String toStr(String[] array, boolean[] display, String delimiter,
                             String nullValue) {
    return toStr((Object[]) array, display, delimiter, nullValue);
  }

  /**
   * Prints the first element of each array in an array of arrays
   *
   * @param array an array of arrays
   * @param delimiter String delimiter
   * @return String of printed objects
   */
  public static String toStr(String[][] array, String delimiter) {
    String str = "";

    for (int i = 0; i < array.length; i++) {
      str += (i == 0 ? "" : delimiter) + array[i][0];
    }

    return str;
  }

  /**
   * Prints an array of objects separated by a tab
   *
   * @param array an array of objects
   * @return String of printed objects
   */
  public static String toStr(Object[] array) {
    return toStr(array, null, "\t", null);
  }

  /**
   * Prints an array of objects separated by the specified delimiter
   *
   * @param array an array of objects
   * @param delimiter String delimiter
   * @return String of printed objects
   */
  public static String toStr(Object[] array, String delimiter) {
    return toStr(array, null, delimiter, null);
  }

  /**
   * @param array an array of objects
   * @param delimiter String delimiter
   * @param nullValue String to use in place of null
   * @return String of printed objects
   */
  public static String toStr(Object[] array, String delimiter, String nullValue) {
    return toStr(array, null, delimiter, nullValue);
  }

  /**
   * Prints an array of objects separated by the specified delimiter
   *
   * @param array an array of objects
   * @param display boolean array indicating which values to print
   * @param delimiter String delimiter
   * @param nullValue String to use in place of null
   * @return String of printed objects
   */
  public static String toStr(Object[] array, boolean[] display, String delimiter,
                             String nullValue) {
    boolean commaDelimited = delimiter.equals(",");
    List<String> cleanList = Lists.newArrayList();

    if (nullValue == null) {
      nullValue = "null";
    }

    for (int i = 0; i < array.length; i++) {
      if (display == null || display[i]) {
        String val = array[i] == null ? nullValue : array[i].toString();
        if (commaDelimited && val.contains(",")) {
          val = "\"" + val + "\"";
        }
        cleanList.add(val);
      }
    }

    return Joiner.on(delimiter).join(cleanList);
  }

  /**
   * Prints an array of integers separated by a tab
   *
   * @param array an array of integers
   * @return String of printed integers
   */
  public static String toStr(int[] array) {
    return toStr(array, "\t");
  }

  /**
   * Prints an array of integers separated by the specified delimiter
   *
   * @param array an array of integers
   * @param delimiter String delimiter
   * @return String of printed integers
   */
  public static String toStr(int[] array, String delimiter) {
    String str = "";

    for (int i = 0; i < array.length; i++) {
      str += (i == 0 ? "" : delimiter) + array[i];
    }

    return str;
  }

  /**
   * Prints an array of booleans separated by the specified delimiter
   *
   * @param array an array of booleans
   * @param delimiter String delimiter
   * @return String of printed integers
   */
  public static String toStr(boolean[] array, String delimiter) {
    String str = "";

    for (int i = 0; i < array.length; i++) {
      str += (i == 0 ? "" : delimiter) + array[i];
    }

    return str;
  }

  /**
   * Prints an array of bytes separated by the specified delimiter
   *
   * @param array an array of bytes
   * @param delimiter String delimiter
   * @return String of printed bytes
   */
  public static String toStr(byte[] array, String delimiter) {
    String str = "";

    for (int i = 0; i < array.length; i++) {
      str += (i == 0 ? "" : delimiter) + array[i];
    }

    return str;
  }

  /**
   * Prints an array of numbers with as many sigfigs as necessary, each separated by a tab
   *
   * @param array an array of numbers
   * @return String of printed numbers
   */
  public static String toStr(double[] array) {
    return toStr(array, "\t");
  }

  /**
   * Prints an array of numbers with as many sigfigs as necessary, each separated by a given
   * delimiter
   *
   * @param array an array of numbers
   * @return String of printed numbers
   */
  public static String toStr(double[] array, String delim) {
    return toStr(array, -1, -1, delim);
  }

  /**
   * Prints an array of numbers separated by the specified delimiter
   *
   * @param array an array of numbers
   * @param sigfigs number of significant digits
   * @param delimiter String delimiter
   * @return String of printed numbers
   */
  public static String toStr(double[] array, int minSigFigs, int maxSigFigs, String delimiter) {
    StringBuilder str = new StringBuilder("");

    for (int i = 0; i < array.length; i++) {
      str.append((i == 0 ? "" : delimiter)
                 + (maxSigFigs == -1 ? ext.formDeci(array[i], 10)
                                     : ext.formDeci(array[i], minSigFigs, maxSigFigs)));
    }

    return str.toString();
  }

  /**
   * Prints an array of numbers with as many sigfigs as necessary, each separated by a tab
   *
   * @param array an array of numbers
   * @return String of printed numbers
   */
  public static String toStr(float[] array) {
    return toStr(array, -1, -1, "\t");
  }

  /**
   * Prints an array of numbers separated by the specified delimiter
   *
   * @param array an array of numbers
   * @param sigfigs number of significant digits
   * @param delimiter String delimiter
   * @return String of printed numbers
   */
  public static String toStr(float[] array, int minSigFigs, int maxSigFigs, String delimiter) {
    String str = "";

    for (int i = 0; i < array.length; i++) {
      str += (i == 0 ? "" : delimiter)
             + (maxSigFigs == -1 ? ext.formDeci(array[i], 10)
                                 : ext.formDeci(array[i], minSigFigs, maxSigFigs));
    }

    return str;
  }

  /**
   * Breaks an array nChunks <br>
   * Warning, one of my first times with the <T> stuff
   *
   * @param array the array
   * @param nChunks number of chunks
   * @param log
   */

  public static <T> List<T[]> splitUpArray(T[] array, int nChunks, Logger log) {
    int index = 0;
    if (array.length < nChunks) {
      log.reportError("Error - too many chunks (" + nChunks + ") for " + array.length
                      + " things, setting to " + array.length);
      nChunks = array.length;
    }
    if (nChunks <= 0) {
      log.reportError("Error - not enough chunks (" + nChunks + ") for " + array.length
                      + " things, setting to 1");
      nChunks = 1;
    }
    int[] chunks = ArrayUtils.splitUpDistributeRemainder(array.length, nChunks, log);
    ArrayList<T[]> da = new ArrayList<>();
    int start = 0;
    for (int chunk : chunks) {
      for (int j = 0; j < chunk; j++) {
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
   * Breaks an array of strings into nChunks. This is geared toward splitting up filenames etc for
   * batching
   *
   * @param strings the array
   * @param nChunks number of chunks
   * @param log
   */
  public static String[][] splitUpStringArray(String[] strings, int nChunks, Logger log) {
    return splitUpArray(strings, nChunks, log).toArray(new String[][] {});
  }

  /**
   * Breaks an array of strings into nChunks with boolean representation, each boolean array has the
   * same length as the original input
   *
   * @param strings the array
   * @param nChunks number of chunks
   * @param log
   */
  public static boolean[][] splitUpStringArrayToBoolean(String[] strings, int nChunks, Logger log) {
    String[][] stringSplits = splitUpStringArray(strings, nChunks, log);
    boolean[][] stringBoolSplits = new boolean[stringSplits.length][];
    for (int i = 0; i < stringBoolSplits.length; i++) {
      int[] indicesThisChunk = ext.indexLargeFactors(stringSplits[i], strings, true, log, true);
      stringBoolSplits[i] = new boolean[strings.length];
      Arrays.fill(stringBoolSplits[i], false);
      for (int j = 0; j < indicesThisChunk.length; j++) {
        stringBoolSplits[i][indicesThisChunk[j]] = true;
      }
    }
    return stringBoolSplits;
  }

  /**
   * Returns an array splitting a number as equally as possible into different amounts
   *
   * @param total number to be split into groups
   * @param numSplits number of groups to split total into
   * @return array of the numbers for each group
   */
  public static int[] splitUp(int total, int numSplits) {
    int[] splits = new int[numSplits];
    int fullBinAmt = (int) Math.floor((double) total / (double) numSplits);
    for (int i = 0; i < numSplits - 1; i++) {
      splits[i] = fullBinAmt;
    }
    splits[numSplits - 1] = total - (numSplits - 1) * fullBinAmt;

    return splits;
  }

  /**
   * Returns an array splitting a number equally, and distributes the remainder
   *
   * @param total number to be split into groups
   * @param numSplits number of groups to split total into
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
    if (ArrayUtils.sum(splits) != total) {
      log.reportError("Internal Error - could not properly split up " + total + " into "
                      + numSplits);
      splits = null;
    }
    return splits;
  }

  public static int[][] splitUpIntoBinsOfIndices(String[] values, Set<String> drops, int binSizeMax,
                                                 Logger log) {
    List<int[]> batches = new ArrayList<>();
    if (drops == null || drops.isEmpty()) {
      // break into batches normally
      int[] batch = new int[Math.min(binSizeMax, values.length)];
      for (int i = 0; i < values.length; i++) {
        if (i > 0 && i % binSizeMax == 0) {
          batches.add(batch);
          batch = new int[Math.min(binSizeMax, values.length - i)];
        }
        batch[i % binSizeMax] = i;
      }
      batches.add(batch);
    } else {
      List<Integer> orphans = new ArrayList<>();
      for (int i = 0; i < values.length; i++) {
        if (drops.contains(values[i])) {
          continue;
        } else {
          int s = i;
          int e = i;
          for (int j = i + 1; j < values.length; j++) {
            if (drops.contains(values[j])) {
              break;
            }
            e++;
          }
          while (e - s + 1 >= binSizeMax) {
            int[] batch = arrayOfIndices(binSizeMax, s);
            batches.add(batch);
            s += binSizeMax;
          }
          for (int j = s; j <= e; j++) {
            orphans.add(j);
          }
          i = e + 1; // skip next, as we know it's in "complete"
        }
      }
      while (orphans.size() > binSizeMax) {
        int[] batch = new int[binSizeMax];
        for (int i = 0; i < batch.length; i++) {
          batch[i] = orphans.remove(0);
        }
        batches.add(batch);
      }
      if (!orphans.isEmpty()) {
        int[] batch = new int[orphans.size()];
        for (int i = 0; i < batch.length; i++) {
          batch[i] = orphans.get(i);
        }
        batches.add(batch);
      }
    }

    return batches.toArray(new int[batches.size()][]);
  }

  /**
   * Creates an array of Strings and copies the contents of a {@link Collection} into it
   *
   * @param collection {@link Collection} of Strings
   * @return an array of Strings from the {@link Collection}
   */
  public static String[] toStringArray(Collection<String> collection) {
    return collection.toArray(new String[collection.size()]);
  }

  /**
   * Creates an array of Strings and copies the contents of a String[][] into it
   *
   * @param matrix matrix of String
   * @param delimiter delimiter to use in the concatenated result
   * @return an array of Strings from the matrix
   */
  public static String[] toStringArray(String[][] matrix, String delimiter) {
    String[] array = new String[matrix.length];

    for (int i = 0; i < array.length; i++) {
      array[i] = ArrayUtils.toStr(matrix[i], delimiter);
    }

    return array;
  }

  /**
   * Creates an array and copies the contents of a List into it in the specified order
   * 
   * @param <T>
   * @param list
   * @param order desired order of elements
   * @return an ordered array from the List
   */
  public static <T> T[] toStringArray(List<T> list, int[] order) {
    if (order.length != list.size()) {
      System.err.println("Error - order does not have the same number of elements (n="
                         + order.length + ") as the List (n=" + list.size() + ")");
      return null;
    }
    @SuppressWarnings("unchecked")
    T[] array = (T[]) new Object[list.size()];

    for (int i = 0; i < array.length; i++) {
      array[i] = list.get(order[i]);
    }
    return array;
  }

  /**
   * Creates an array of Strings and copies the contents of an ArrayList into it
   *
   * @param v vector of Strings
   * @return an array of Strings from the ArrayList
   */
  public static String[] toStringArray(List<String> al) {
    return al.toArray(new String[al.size()]);
  }

  /**
   * Creates an array of Strings and copies the contents of an array of booleans
   *
   * @param array array of boolean
   * @param onesAndZeros array of boolean
   * @return an array of the converted Strings
   */
  public static String[] toStringArray(boolean[] array, boolean onesAndZeros) {
    String[] new_array;

    new_array = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      new_array[i] = onesAndZeros ? (array[i] ? "1" : "0") : array[i] + "";
    }

    return new_array;
  }

  /**
   * Creates an array of Strings and copies the contents of an array of int
   *
   * @param array array of int
   * @return an array of the converted Strings
   */
  public static String[] toStringArray(int[] array) {
    String[] new_array;

    new_array = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      new_array[i] = array[i] + "";
    }

    return new_array;
  }

  /**
   * Creates an array of Strings and copies the contents of an array of long
   *
   * @param array array of long
   * @return an array of the converted Strings
   */
  public static String[] toStringArray(long[] array) {
    String[] new_array;

    new_array = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      new_array[i] = array[i] + "";
    }

    return new_array;
  }

  /**
   * Creates an array of Strings and copies the contents of an array of double
   *
   * @param array array of double
   * @return an array of the converted Strings
   */
  public static String[] toStringArray(double[] array) {
    String[] new_array;

    new_array = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      new_array[i] = array[i] + "";
    }

    return new_array;
  }

  /**
   * Creates an array of Strings and copies the contents of an array of float
   *
   * @param array array of float
   * @return an array of the converted Strings
   */
  public static String[] toStringArray(float[] array) {
    String[] new_array;

    new_array = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      new_array[i] = array[i] + "";
    }

    return new_array;
  }

  /**
   * Creates an array and copies the Keys of a Map into it according to the order specified by the
   * Values
   *
   * @param map Map with intended index position as the Value
   * @return array of Keys from map indexed by Values
   */
  @SuppressWarnings("unchecked")
  public static <T> T[] mapToValueSortedArray(Map<T, Integer> map) {
    T[] array = null;
    for (Map.Entry<T, Integer> entry : map.entrySet()) {
      if (array == null) {
        array = (T[]) java.lang.reflect.Array.newInstance(entry.getKey().getClass(), map.size());
      }
      array[entry.getValue()] = entry.getKey();
    }
    return array;
  }

  /**
   * Creates a Vector and copies the conents of a String array into it
   *
   * @param array an array of Strings
   * @return a vector of Strings
   */
  public static Vector<String> toStringVector(String[] array) {
    Vector<String> v = new Vector<>();
    for (String element : array) {
      v.add(element);
    }
    return v;
  }

  public static <T> List<T> toList(T[] array) {
    return Lists.newArrayList(array);
  }

  /**
   * Returns an array of the unique Strings
   *
   * @param array an array of Strings
   * @return array of the unique Strings
   */
  public static String[] unique(String[] array) {
    Set<String> filter = new LinkedHashSet<>();

    for (int i = 0; i < array.length; i++) {
      filter.add(array[i]);
    }

    return filter.toArray(new String[filter.size()]);
  }

  public static <T extends Number> List<T> unique(List<T> list) {
    Set<T> filter = new LinkedHashSet<>(list);

    return new ArrayList<>(filter);
  }

  // /**
  // * Returns an array minus any entries with a null or value listed as invalid
  // *
  // * @param array
  // * an array of Strings
  // * @param thoseToBeRemoved
  // * an array of Strings representing invalid values
  // * @return array of valid Strings
  // */
  // public static String[] removeThese(String[] array, String[] thoseToBeRemoved) {
  // Vector<String> v = new Vector<String>();
  //
  // for (int i = 0; i<array.length; i++) {
  // if (array[i]!=null&&ext.indexOfStr(array[i], thoseToBeRemoved)==-1) {
  // v.add(array[i]);
  // }
  // }
  //
  // return Array.toStringArray(v);
  // }

  /**
   * Inserts the specified String into an array the specified position
   *
   * @param str String to be inserted
   * @param array an array of Strings
   * @param pos position for str to be inserted
   * @return altered array of Strings
   */
  public static String[] insertStringAt(String str, String[] array, int pos) {
    if (pos < 0 || pos > array.length) {
      throw new ArrayIndexOutOfBoundsException(pos);
    }
    List<String> list = toList(array);
    list.add(pos, str);

    return ArrayUtils.toStringArray(list);
  }

  /**
   * Populate an array of doubles by converting the indicated members of a String array
   *
   * @param line array of Strings
   * @param indices indices of the Strings to be converted to doubles
   * @return the resulting double array
   */
  public static double[] populateDoubles(String[] line, int[] indices) {
    double[] array = new double[indices.length];

    for (int i = 0; i < indices.length; i++) {
      array[i] = line[indices[i]].equals(".") ? -999 : Double.parseDouble(line[indices[i]]);
    }

    return array;
  }

  /**
   * convert a string array into a boolean representation <br>
   * masks will all be set to false
   */
  public static BooleanClassifier classifyStringsToBoolean(String[] toClassify, String[] masks) {
    String[] uniqs = unique(toClassify);
    ArrayList<String> uniqNoMask = new ArrayList<>();
    for (String uniq : uniqs) {
      if (ext.indexOfStr(uniq, masks) < 0) {
        uniqNoMask.add(uniq);
      }
    }
    uniqs = ArrayUtils.toStringArray(uniqNoMask);
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

    private final boolean[][] classified;
    private final String[] titles;

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

  // /**
  // * Create an array of initialized IntVectors
  // *
  // * @param int
  // * size of array
  // * @return the array of IntVectors
  // */
  // public static IntVector[] intVectorArray(int size) {
  // IntVector[] array = new IntVector[size];
  //
  // for (int i = 0; i<size; i++) {
  // array[i] = new IntVector();
  // }
  //
  // return array;
  // }
  //
  // /**
  // * Create an array of initialized DoubleVectors
  // *
  // * @param int
  // * size of array
  // * @return the array of Vectors
  // */
  // public static DoubleVector[] doubleVectorArray(int size) {
  // DoubleVector[] array = new DoubleVector[size];
  //
  // for (int i = 0; i<size; i++) {
  // array[i] = new DoubleVector();
  // }
  //
  // return array;
  // }

  /**
   * Find all instances of the values found in the first column of the matrix and replace them with
   * the values from the second column
   *
   * @param array an array of integers
   * @param swaps a matrix of old and new values
   * @return the same array but with the swapped values
   */
  public static int[] findReplace(int[] array, int[][] swaps) {
    int[] newArray = array.clone();

    for (int i = 0; i < array.length; i++) {
      for (int[] swap : swaps) {
        if (array[i] == swap[0]) {
          newArray[i] = swap[1];
        }
      }
    }

    return newArray;
  }

  /**
   * Tries to find the instance in a sorted array where all values up to, but not including, that
   * index are less than a given maximum target
   * <p>
   * For example, calling {@link ArrayUtils#indexOfLastMinByte} using (new byte[] {0,1,2,24,25},
   * 23), would return 3 (all values up to index 3 are less than 23);
   *
   * @param array an array of bytes
   * @param maxByte the number to find
   * @return the index, or -1 if all values were greater than the maximum, or the array's length if
   *         all values were less than the maximum
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
      return array.length;// all values were less than max
    } else {
      return max; // all values were greater than max
    }
  }

  /**
   * Creates a new array using only the indices at and after start
   *
   * @param array an array
   * @param start first index to use (inclusive)
   * @return the subset of the original array
   */
  public static <T> T[] subArray(T[] array, int start) {
    return subArray(array, start, array.length);
  }

  /**
   * Creates a new array using only the indices at and after start and before end
   *
   * @param array an array
   * @param start first index to use (inclusive)
   * @param stop index to stop at (non-inclusive)
   * @return the subset of the original array
   */
  public static <T> T[] subArray(T[] array, int start, int end) {
    return Arrays.copyOfRange(array, start, end);
  }

  @SuppressWarnings("unchecked")
  public static <T> T[] subArray(T[] array, boolean[] use) {
    T[] subarray;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
      return null;
    }

    count = 0;
    subarray = (T[]) java.lang.reflect.Array.newInstance(array[0].getClass(), booleanArraySum(use));// new
                                                                                                    // T[booleanArraySum(use)];

    for (int i = 0; i < array.length; i++) {
      if (use[i]) {
        subarray[count] = array[i];
        count++;
      }
    }

    return subarray;
  }

  /**
   * Creates a new array using only the indices at and after start
   *
   * @param array an array of ints
   * @param start first index to use
   * @return the subset of the original array
   */
  public static int[] subArray(int[] array, int start) {
    return Arrays.copyOfRange(array, start, array.length);
  }

  /**
   * Creates a new array using only the int values at indices defined by the Integer array
   *
   * @param array an array of double
   * @param use indices to use
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
      System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n="
                         + array.length + ")");
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
      System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n="
                         + array.length + ")");
      return null;
    }
    return subarray;
  }

  /**
   * Creates a new array using only the byte values at indices defined by the Integer array
   *
   * @param array an array of double
   * @param use indices to use
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
      System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n="
                         + array.length + ")");
      return null;
    }
    return subarray;
  }

  /**
   * Creates a new array using only the float values at indices defined by the Integer array
   *
   * @param array an array of double
   * @param use indices to use
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
      System.err.println("Error - out of bounds index (" + currentIndex + ") for subsetting (n="
                         + array.length + ")");
      return null;
    }
    return subarray;
  }

  /**
   * Creates a new array using only the strings at indices with a true in the boolean array
   *
   * @param array an array of Strings
   * @param use indices to use
   * @return the subset of the original array
   */
  public static int[] subArray(int[] array, boolean[] use) {
    int[] subarray;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
      return null;
    }

    count = 0;
    subarray = new int[booleanArraySum(use)];
    for (int i = 0; i < array.length; i++) {
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
   * @param array an array of byte
   * @param use indices to use
   * @return the subset of the original array
   */
  public static byte[] subArray(byte[] array, boolean[] use) {
    byte[] subarray;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
      return null;
    }

    count = 0;
    subarray = new byte[booleanArraySum(use)];
    for (int i = 0; i < array.length; i++) {
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
   * @param array an array of doubles
   * @param start first index to use
   * @param stop last index to use
   * @return the subset of the original array
   */
  public static double[] subArray(double[] array, int start, int stopBefore) {
    double[] arr;

    if (start < 0 || stopBefore > array.length || stopBefore <= start) {
      System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore
                         + ") indicies for an array");
    }
    arr = new double[stopBefore - start];
    for (int i = start; i < stopBefore; i++) {
      arr[i - start] = array[i];
    }

    return arr;
  }

  /**
   * Creates a new array using only the indices at and after start
   *
   * @param array an array of doubles
   * @param start first index to use
   * @param stop last index to use
   * @return the subset of the original array
   */
  public static double[] subArray(double[] array, int start) {
    return subArray(array, start, array.length);
  }

  /**
   * Creates a new array using only the double values at indices with a true in the boolean array
   *
   * @param array an array of double
   * @param use indices to use
   * @return the subset of the original array
   */
  public static double[] subArray(double[] array, boolean[] use) {
    double[] subarray;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
      return null;
    }

    count = 0;
    subarray = new double[booleanArraySum(use)];
    for (int i = 0; i < array.length; i++) {
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
   * @param array an array of double
   * @param use indices to use
   * @return the subset of the original array
   */
  public static double[][] subArray(double[][] array, boolean[] use) {
    double[][] subarray;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
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
   * @param array an array of double
   * @param use indices to use
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
   * @param array an array of double
   * @param use indices to use
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
   * @param array an array of double[]
   * @param start first index to use
   * @param stop last index to use
   * @return the subset of the original array
   */
  public static double[][] subArray(double[][] array, int start, int stopBefore) {
    double[][] arr;

    if (start < 0 || stopBefore > array.length || stopBefore <= start) {
      System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore
                         + ") indicies for an array");
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
   * @param array an array of floats
   * @param start first index to use
   * @param stop last index to use
   * @return the subset of the original array
   */
  public static float[] subArray(float[] array, int start, int stopBefore) {
    float[] arr;

    if (start < 0 || stopBefore > array.length || stopBefore <= start) {
      System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore
                         + ") indicies for an array");
    }
    arr = new float[stopBefore - start];
    for (int i = start; i < stopBefore; i++) {
      arr[i - start] = array[i];
    }

    return arr;
  }

  /**
   * Creates a new array using only the indices between start and stop
   *
   * @param array an array of boolean
   * @param start first index to use
   * @param stop last index to use
   * @return the subset of the original array
   */
  public static boolean[] subArray(boolean[] array, int start, int stopBefore) {
    boolean[] arr;
    if (start < 0 || stopBefore > array.length || stopBefore <= start) {
      System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore
                         + ") indicies for an array");
    }
    arr = new boolean[stopBefore - start];
    for (int i = start; i < stopBefore; i++) {
      arr[i - start] = array[i];
    }
    return arr;
  }

  /**
   * As {@link #subArrayInRange(float[], boolean[], float, float)}, using all samples by default.
   */
  public static float[] subArrayInRange(float[] array, float min, float max) {
    return subArrayInRange(array, null, min, max);
  }

  /**
   * Creates a new array using only the float values within the given range. Prunes out
   * {@link Float#NaN} as well.
   *
   * @param array base array to filter
   * @param use indices of array to use
   * @param min smallest value to include in array
   * @param max largest value to include in array
   * @return an array of values filtered to the specified range
   */
  public static float[] subArrayInRange(float[] array, boolean[] use, float min, float max) {
    boolean[] samples = new boolean[array.length];
    for (int i = 0; i < array.length; i++) {
      samples[i] = (use == null || use[i]) && !Float.isNaN(array[i])
                   && (Float.compare(array[i], min) >= 0 && Float.compare(array[i], max) <= 0);
    }
    return subArray(array, samples);
  }

  /**
   * Creates a new array using only the float values at indices with a true in the boolean array
   *
   * @param array an array of float
   * @param use indices to use
   * @return the subset of the original array
   */
  public static float[] subArray(float[] array, boolean[] use) {
    float[] subarray;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
      return null;
    }

    count = 0;
    subarray = new float[booleanArraySum(use)];
    for (int i = 0; i < array.length; i++) {
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
   * @param array an array of Strings
   * @param start first index to use
   * @param stop last index to use (exclusive)
   * @return the subset of the original array
   */
  public static String[] subArray(String[] array, int start, int stopBefore) {
    String[] arr;

    if (start < 0 || stopBefore > array.length || stopBefore < start) {
      System.err.println("Error - invalid start (" + start + ") and stopBefore (" + stopBefore
                         + ") indicies for an array of size " + array.length);
    }
    arr = new String[stopBefore - start];
    for (int i = start; i < stopBefore; i++) {
      arr[i - start] = array[i];
    }

    return arr;
  }

  /**
   * Creates a new array using only the indices at and after start
   *
   * @param array an array of Strings
   * @param start first index to use
   * @return the subset of the original array
   */
  public static String[] subArray(String[] array, int start) {
    return subArray(array, start, array.length);
  }

  /**
   * Creates a new array using only the indices provided
   *
   * @param array an array of Strings
   * @param indices indices to use
   * @return the subset of the original array
   */
  public static String[] subArray(String[] array, int[] indices) {
    return subArray(array, replaceAllWith(indices, -1, -2), null);
  }

  /**
   * Creates a new array using only the indices provided
   *
   * @param array an array of Strings
   * @param indices indices to use
   * @param missingValue if an index equals -1, then missingValue will be inserted
   * @return the subset of the original array
   */
  public static String[] subArray(String[] array, int[] indices, String missingValue) {
    String[] strs = new String[indices.length];

    if (indices.length == 0) {
      return strs;
    }
    if (min(indices) < -1) {
      System.err.println("Error missing an index; subarray will fail");
      System.exit(1);
    }
    for (int i = 0; i < strs.length; i++) {
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
   * @param array an array of Strings
   * @param use indices to use
   * @return the subset of the original array
   */
  public static String[] subArray(String[] array, boolean[] use) {
    String[] strs;
    int count;

    if (array.length != use.length) {
      System.err.println("Error - mismatched array lengths for the aray (n=" + array.length
                         + ") and the boolean subset (n=" + use.length + ")");
      return null;
    }

    count = 0;
    strs = new String[booleanArraySum(use)];
    for (int i = 0; i < array.length; i++) {
      if (use[i]) {
        strs[count] = array[i];
        count++;
      }
    }

    return strs;
  }

  /**
   * Create a sub-list using only indicies with {@code true} values in the given use array.
   */
  public static <T> List<T> subList(List<T> listIn, boolean[] use) {
    List<T> listOut = new ArrayList<>();
    for (int i = 0; i < listIn.size(); i++) {
      if (use[i]) {
        listOut.add(listIn.get(i));
      }
    }
    return listOut;
  }

  /**
   * Increments a boolean array
   *
   * @param array an array of boolean
   */
  public static void incBoolean(boolean[] array) {
    for (int i = array.length - 1; i >= 0; i--) {
      if (!array[i]) {
        array[i] = true;
        for (int j = i + 1; j < array.length; j++) {
          array[j] = false;
        }
        return;
      }
    }

  }

  /**
   * Provides a frequency distribution of an array of Strings
   *
   * @param array an array of Strings
   * @return the unique elements and their counts
   */
  public static String[][] frequency(String[] array) {
    Hashtable<String, String> hash = new Hashtable<>();
    int count;
    String[] keys;
    String[][] summary;

    for (String element : array) {
      if (hash.containsKey(element)) {
        count = Integer.parseInt(hash.get(element));
      } else {
        count = 0;
      }
      count++;
      hash.put(element, count + "");
    }

    keys = HashVec.getKeys(hash);
    summary = new String[keys.length][2];
    for (int i = 0; i < summary.length; i++) {
      summary[i][0] = keys[i];
      summary[i][1] = hash.get(keys[i]);
    }

    return summary;
  }

  public static int booleanArraySum(boolean[] array) {
    int count = 0;

    for (boolean element : array) {
      if (element) {
        count++;
      }
    }

    return count;
  }

  /**
   * Takes two boolean arrays and sets the first array to the boolean AND of both arrays for a given
   * element.
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

  public static <T extends Comparable<T>> int binarySearch(List<T[]> list, T[] value, int keyIndex,
                                                           boolean exact) {
    return binarySearch(list, value, keyIndex, 0, list.size() - 1, exact);
  }

  public static <T extends Comparable<T>> int binarySearch(List<T[]> list, T[] value, int keyIndex,
                                                           int low, int high, boolean exact) {
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
    return binarySearch(array, value, 0, array.length - 1, exact);
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
    return binarySearch(array, value, 0, array.length - 1, exact);
  }

  public static int binarySearch(int[] array, int value, int low, int high, boolean exact) {
    int mid;

    while (low <= high) {
      // System.out.println(array[low]+"\t"+array[high]);
      mid = low + (high - low) / 2;
      if (mid >= array.length) {
        if (exact) {
          return -9;
        } else {
          return array.length;
        }
      }
      if (array[mid] > value) {
        high = mid - 1;
      } else if (array[mid] < value) {
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

  public static String[] booleanArrayToStringArray(boolean[] array) {
    String[] strs = new String[array.length];

    for (int i = 0; i < array.length; i++) {
      strs[i] = array[i] ? "1" : "0";
    }

    return strs;
  }

  public static int countIf(String[] array, String target) {
    int count = 0;

    for (String element : array) {
      if (element.equals(target)) {
        count++;
      }
    }

    return count;
  }

  public static int countIf(int[] array, int target) {
    int count = 0;

    for (int element : array) {
      if (element == target) {
        count++;
      }
    }

    return count;
  }

  public static boolean[] booleanNegative(boolean[] array) {
    boolean[] newArray = new boolean[array.length];

    for (int i = 0; i < array.length; i++) {
      newArray[i] = !array[i];
    }

    return newArray;
  }

  /**
   * Calculates the interquartile range of an array (exclusive)
   *
   * @param array an array of numbers
   * @return iqr of the array
   */
  public static double iqrExclusive(double[] array) {
    if (array.length < 2) {
      System.err.println("Error - can't calculate an IQR for an array with " + array.length
                         + " datapoint(s)");
      return -1;
    }

    int[] keys = Sort.getSortedIndices(array);

    double iqr = 0;
    try {
      iqr = array[keys[(int) Math.floor(array.length * 0.75)]]
            - array[keys[(int) Math.floor(array.length * 0.25)]];
    } catch (Exception e) {
      System.err.println("Error calculating IQR");
      e.printStackTrace();
    }
    return iqr;
  }

  /**
   * Calculates the interquartile range of an array (exclusive)
   *
   * @param array an array of numbers
   * @return iqr of the array
   */
  public static float iqrExclusive(float[] array) {
    if (array.length < 2) {
      System.err.println("Error - can't calculate an IQR for an array with " + array.length
                         + " datapoint(s)");
      return -1;
    }

    int[] keys = Sort.getSortedIndices(array);
    float iqr = 0;

    try {
      iqr = array[keys[(int) Math.floor(array.length * 0.75)]]
            - array[keys[(int) Math.floor(array.length * 0.25)]];
    } catch (Exception e) {
      System.err.println("Error calculating IQR");
      e.printStackTrace();
    }

    return iqr;
  }

  /**
   * Trims null values from the end of an array
   * 
   * @param <T>
   * @param array an array of Strings
   * @return trimmed array
   */
  public static <T> T[] trimArray(T[] array) {
    int index;

    index = array.length;
    while (index > 0 && array[index - 1] == null) {
      index--;
    }

    return index == array.length ? array : subArray(array, 0, index);
  }

  /**
   * Returns true if the String arrays are equal at all positions
   *
   * @param array1 an array of Strings
   * @param array2 an array of Strings
   * @param caseSensitive boolean flag
   * @return true if arrays are equal
   */
  public static boolean equals(String[] array1, String[] array2, boolean caseSenstitive) {
    if (array1.length != array2.length) {
      return false;
    }

    for (int i = 0; i < array1.length; i++) {
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
   * @param array1 an array of Strings
   * @param array2 an array of Strings
   * @param caseSensitive boolean flag
   * @return true if arrays are equal
   */
  public static boolean equals(int[] array1, int[] array2) {
    if (array1.length != array2.length) {
      return false;
    }

    for (int i = 0; i < array1.length; i++) {
      if (array1[i] != array2[i]) {
        return false;
      }
    }

    return true;
  }

  /**
   * Returns true if the byte arrays are equal at all positions
   *
   * @param array1 an array of byte
   * @param array2 an array of byte
   * @param caseSensitive boolean flag
   * @return true if arrays are equal
   */
  public static boolean equals(byte[] array1, byte[] array2) {
    if (array1.length != array2.length) {
      return false;
    }

    for (int i = 0; i < array1.length; i++) {
      if (array1[i] != array2[i]) {
        return false;
      }
    }

    return true;
  }

  /**
   * @param array
   * @param sigFigs all entries in the array will be rounded the this many sig figs
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
   * @param array an array of doubles
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
   * @param array an array of doubles
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
   * @param array an array of doubles
   * @return scrubbed array
   */
  public static double[] removeNonFinites(double[] array) {
    return Arrays.stream(array).filter(Double::isFinite).toArray();
  }

  /**
   * Removes non-finite values from the array
   *
   * @param array an array of floats
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
   * @param array an array of bytes
   * @param valueToRemove value to remove from the array
   * @return filtered array
   */
  public static byte[] removeAllValues(byte[] array, byte valueToRemove) {
    boolean[] use;

    use = new boolean[array.length];
    for (int i = 0; i < use.length; i++) {
      use[i] = array[i] != valueToRemove;
    }

    return subArray(array, use);
  }

  /**
   * Removes all instances of a specified value from an array
   *
   * @param array an array of int
   * @param valueToRemove value to remove from the array
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
    ArrayList<Double> tmp = new ArrayList<>();
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i])
          && (array[i] > min && array[i] < max || (gteLte && array[i] >= min && array[i] <= max))) {
        tmp.add(array[i]);
      }
    }
    return Doubles.toArray(tmp);
  }

  /**
   * Create a new array with any non-finite value replaced with {@link Float#NaN}.
   *
   * @param array input array
   * @return a copy of the input array with all non-finite values set to NaN.
   */
  public static float[] replaceNonFinites(float[] array) {
    float[] ret = new float[array.length];
    for (int i = 0; i < array.length; i++) {
      if (Float.isFinite(array[i])) {
        ret[i] = array[i];
      } else {
        ret[i] = Float.NaN;
      }
    }
    return ret;
  }

  public static float[] getValuesBetween(float[] array, double min, double max, boolean gteLte) {
    ArrayList<Float> tmp = new ArrayList<>();
    for (int i = 0; i < array.length; i++) {
      if (!Double.isNaN(array[i])
          && (array[i] > min && array[i] < max || (gteLte && array[i] >= min && array[i] <= max))) {
        tmp.add(array[i]);
      }
    }
    return Floats.toArray(tmp);
  }

  public static int[] getValuesBetween(int[] array, int min, int max, boolean gteLte) {
    ArrayList<Integer> tmp = new ArrayList<>();
    for (int element : array) {
      if ((element > min && element < max) || (gteLte && element >= min && element <= max)) {
        tmp.add(element);
      }
    }
    return Ints.toArray(tmp);
  }

  /**
   * Reverses the entries of the array, first becomes last, etc
   */
  public static String[] reverse(String[] forward) {
    String[] reverse = new String[forward.length];
    return reverse(forward, reverse);
  }

  /**
   * As {@link #reverse(String[])} but happens in place
   */
  public static String[] reverseInPlace(String[] forward) {
    return reverse(forward, forward);
  }

  private static String[] reverse(String[] forward, String[] reverse) {
    for (int i = 0; i < (reverse.length + 1) / 2; i++) {
      String t = forward[i];
      reverse[i] = forward[forward.length - 1 - i];
      reverse[forward.length - 1 - i] = t;
    }
    return reverse;
  }

  /**
   * Creates a new array holding the reverse entries of the given array
   */
  public static int[] reverse(int[] forward) {
    int[] reverse = new int[forward.length];
    return reverse(forward, reverse);
  }

  /**
   * As {@link #reverse(int[])} but happens in place
   */
  public static int[] reverseInPlace(int[] forward) {
    return reverse(forward, forward);
  }

  private static int[] reverse(int[] forward, int[] reverse) {
    for (int i = 0; i < (reverse.length + 1) / 2; i++) {
      int t = forward[i];
      reverse[i] = forward[forward.length - 1 - i];
      reverse[forward.length - 1 - i] = t;
    }
    return reverse;
  }

  /**
   * Creates a new array holding the reverse entries of the given array
   */
  public static long[] reverse(long[] forward) {
    long[] reverse = new long[forward.length];
    return reverse(forward, reverse);
  }

  /**
   * As {@link #reverse(long[])} but happens in place
   */
  public static long[] reverseInPlace(long[] forward) {
    return reverse(forward, forward);
  }

  private static long[] reverse(long[] forward, long[] reverse) {
    for (int i = 0; i < (reverse.length + 1) / 2; i++) {
      long t = forward[i];
      reverse[i] = forward[forward.length - 1 - i];
      reverse[forward.length - 1 - i] = t;
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
   * Creates a new array holding the reverse entries of the given array
   */
  public static double[] reverse(double[] forward) {
    double[] reverse = new double[forward.length];
    return reverse(forward, reverse);
  }

  /**
   * As {@link #reverse(double[])} but happens in place
   */
  public static double[] reverseInPlace(double[] forward) {
    return reverse(forward, forward);
  }

  private static double[] reverse(double[] forward, double[] reverse) {
    for (int i = 0; i < (reverse.length + 1) / 2; i++) {
      double t = forward[i];
      reverse[i] = forward[forward.length - 1 - i];
      reverse[forward.length - 1 - i] = t;
    }
    return reverse;
  }

  /**
   * Extract the element at the given index in each sub-array and return a single flat array.
   *
   * @param srcArr Source array of type T[][]
   * @param index Index of elements in subArrays
   * @return null if <code>srcArr</code> is null or empty, <br />
   *         or an array of type T[], with the same length as <code>srcArr</code>, containing only
   *         the elements at <code>index</code> in each sub-array of <code>srcArr</code>
   */
  public static <T> T[] extract(T[][] srcArr, int index) {
    if (srcArr == null) {
      return null;
    }
    if (srcArr.length == 0) {
      return null;
    }
    Class<?> arrClz = srcArr[0][0].getClass();
    @SuppressWarnings("unchecked")
    T[] arr = (T[]) java.lang.reflect.Array.newInstance(arrClz, srcArr.length);// new
                                                                               // Object[srcArr.length];
    for (int i = 0; i < srcArr.length; i++) {
      arr[i] = srcArr[i][index];
    }
    return arr;
  }

  /**
   * Removes certain values from the array
   *
   * @param array an array of Strings
   * @param thingsToRemove list of Strings to removed from the array
   * @return scrubbed array
   */
  public static String[] removeFromArray(String[] array, String[] thingsToRemove) {
    String[] newArray;
    boolean[] use;
    int count;

    use = new boolean[array.length];
    for (int i = 0; i < use.length; i++) {
      use[i] = ext.indexOfStr(array[i], thingsToRemove) == -1;
    }

    newArray = new String[booleanArraySum(use)];
    count = 0;
    for (int i = 0; i < use.length; i++) {
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
   * @param array an array of Strings
   * @param index index of String to remove from the array
   * @return new array
   */
  public static String[] removeFromArray(String[] array, int index) {
    String[] newArray;

    newArray = new String[array.length - 1];
    for (int i = 0; i < index; i++) {
      newArray[i] = array[i];
    }
    for (int i = index + 1; i < array.length; i++) {
      newArray[i - 1] = array[i];
    }

    return newArray;
  }

  /**
   * Add entries to the end of an array
   * 
   * @param array
   * @param additions
   * @return a new array consisting of the entries in array with the entries in additions appended
   *         to the end
   */
  @SafeVarargs
  public static <T> T[] appendToArray(T[] array, T... additions) {
    T[] newArray = Arrays.copyOf(array, array.length + additions.length);
    for (int i = 0; i < additions.length; i++) {
      newArray[i + array.length] = additions[i];
    }
    return newArray;
  }

  /**
   * add a single entry to the end of an array
   * 
   * @param array
   * @param addition
   * @return a new array consisting of the entries in array with addition appended to the end
   */
  public static <T> T[] appendToArray(T[] array, T addition) {
    T[] newArray = Arrays.copyOf(array, array.length + 1);
    newArray[array.length] = addition;
    return newArray;
  }

  /**
   * Adds specified String to the end of an array
   *
   * @param array an array of Strings
   * @param str String to add to the array
   * @return new array
   */
  public static String[] addStrToArray(String str, String[] array) {
    return addStrToArray(str, array, array.length);
  }

  /**
   * Adds specified String to a specified index of an array
   *
   * @param array an array of Strings
   * @param str String to add to the array
   * @param indexOfNewArray location of the string in the new array
   * @return new array
   */
  public static String[] addStrToArray(String str, String[] array, int indexOfNewStr) {
    String[] new_array = new String[array.length + 1];

    // for (int i = 0; i<array.length; i++) {
    for (int i = 0; i < new_array.length; i++) {
      // new_array[i<=indexOfNewStr?i:i+1] = i==indexOfNewStr?str:array[i];
      new_array[i] = i == indexOfNewStr ? str : array[i > indexOfNewStr ? i - 1 : i];
    }

    return new_array;
  }

  /**
   * Adds specified integer to the end of an array
   *
   * @param array an array of integers
   * @param value integer to add to the array
   * @return new array
   */
  public static int[] addIntToArray(int value, int[] array) {
    return addIntToArray(value, array, array.length);
  }

  /**
   * Adds specified integer to a specified index of an array
   *
   * @param array an array of integers
   * @param value integer to add to the array
   * @param indexOfNewArray location of the integer in the new array
   * @return new array
   */
  public static int[] addIntToArray(int value, int[] array, int indexOfNewStr) {
    int[] new_array = new int[array.length + 1];

    for (int i = 0; i < new_array.length; i++) {
      new_array[i] = i == indexOfNewStr ? value : array[i > indexOfNewStr ? i - 1 : i];
    }

    return new_array;
  }

  /**
   * Adds specified double to the end of an array
   *
   * @param array an array of double
   * @param value double to add to the array
   * @return new array
   */
  public static double[] addDoubleToArray(double value, double[] array) {
    return addDoubleToArray(value, array, array.length);
  }

  /**
   * Adds specified double to a specified index of an array
   *
   * @param array an array of double
   * @param value double to add to the array
   * @param indexOfNewArray location of the double in the new array
   * @return new array
   */
  public static double[] addDoubleToArray(double value, double[] array, int indexOfNewStr) {
    double[] new_array = new double[array.length + 1];

    for (int i = 0; i < new_array.length; i++) {
      new_array[i] = i == indexOfNewStr ? value : array[i > indexOfNewStr ? i - 1 : i];
    }

    return new_array;
  }

  /**
   * Determine if the trait contained within the file is dichotomous or continuous
   *
   * @param filename filename containing the trait to be evaluated
   * @param col column to extract
   * @param allow21 allow binary trait to be 2 and 1 instead of 1 and 0
   * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropriate for
   *         neither
   */
  public static int determineType(String filename, int col, boolean allow21) {
    return determineType(filename, col, ext.MISSING_VALUES, allow21);
  }

  /**
   * Determine if the trait contained within the file is dichotomous or continuous
   *
   * @param filename filename containing the trait to be evaluated
   * @param col column to extract
   * @param exclusions values to exclude
   * @param allow21 allow binary trait to be 2 and 1 instead of 1 and 0
   * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropriate for
   *         neither
   */
  public static int determineType(String filename, int col, String[] exclusions, boolean allow21) {
    return determineType(allow21,
                         ArrayUtils.toDoubleArray(ArrayUtils.removeFromArray(HashVec.loadFileToStringArray(filename,
                                                                                                           true,
                                                                                                           new int[] {col},
                                                                                                           true),
                                                                             exclusions)));
  }

  /**
   * Determine if the trait contained within the file is dichotomous or continuous
   *
   * @param array an array of doubles
   * @param allow21 allow binary trait to be 2 and 1 instead of 1 and 0
   * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropraite for
   *         neither
   */
  public static int determineType(double[] array, boolean allow21) {
    List<BigDecimal> bds = new ArrayList<>();
    // Find the first two non-NaN doubles
    for (int i = 0; bds.size() < 3 && i < array.length; i++) {
      if (!Double.isNaN(array[i])) {
        BigDecimal bd = BigDecimal.valueOf(array[i]).round(new MathContext(2));
        bds.add(bd);
      }
    }

    Collections.sort(bds);

    return determineType(allow21,
                         new double[] {bds.get(0).doubleValue(), bds.get(1).doubleValue()});
  }

  /**
   * Determine if the trait contained within the file is dichotomous or continuous
   *
   * @param array an array of doubles
   * @param allow21 allow binary trait to be 2 and 1 instead of 1 and 0
   * @return 0 if appropriate for logistic, 1 if appropriate for linear, -1 if appropraite for
   *         neither
   */
  private static int determineType(boolean allow21, double[] array) {
    if (array.length == 0) {
      System.err.println("Error - no valid data for this trait!");
      return -1;
    } else if (array.length == 1) {
      System.err.println("Error - no variation for this trait!");
      return -1;
    } else if (array.length == 2 && ArrayUtils.min(array) == 0 && ArrayUtils.max(array) == 1) {
      return 0;
    } else if (allow21 && array.length == 2 && ArrayUtils.min(array) == 1
               && ArrayUtils.max(array) == 2) {
      return 0;
    } else if (array.length == 2) {
      System.err.println("Error - flag was set to prevent binary trait from being anything other than 0/1 ("
                         + array[0] + "/" + array[1] + " is not valid)");
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
   * @param array an array of Strings
   * @param tab split using tabs as opposed to white spaces
   * @return an array of arrays of Strings
   */
  public static String[][] splitStrings(String[] array, boolean tab) {
    String[][] stringArrays;

    stringArrays = new String[array.length][];
    for (int i = 0; i < array.length; i++) {
      if (tab) {
        stringArrays[i] = array[i].split("\t", -1);
      } else {
        stringArrays[i] = array[i].split(PSF.Regex.GREEDY_WHITESPACE);
      }
    }

    return stringArrays;
  }

  /**
   * Transposes array to a matrix with a Splits an array of Strings into an array of arrays of
   * Strings
   *
   * @param array an array of Strings
   * @param tab split using tabs as opposed to white spaces
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
   * @param array an array of integers
   * @param from value to replace
   * @param to value to replace with
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
   * @param array1 an array of Strings
   * @param array2 an array of Strings
   * @param numberOfMismatchesAllowed maximum number of mismatched values allowed before the null
   *          set is returned
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
        newArray[i] = array1[i] + "/" + array2[i];
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
   * Returns true if all Strings in the array represent numbers
   *
   * @param array an array of Strings
   * @return merged array of Strings if successful, otherwise null
   */
  public static boolean isAllNumbers(String[] array) {
    for (String element : array) {
      try {
        Double.parseDouble(element);
      } catch (NumberFormatException nfe) {
        return false;
      }
    }

    return true;
  }

  /**
   * Returns the length of the largest String in the array
   *
   * @param array an array of Strings
   * @return largest length
   */
  public static int maxLength(String[] array) {
    int max;

    max = -1;
    for (String element : array) {
      if (element != null) {
        max = Math.max(max, element.length());
      }
    }

    return max;
  }

  public static String[] addPrefixToArray(String[] array, String prefix) {
    return addPrefixSuffixToArray(array, prefix, "");
  }

  public static String[] addSuffixToArray(String[] array, String suffix) {
    return addPrefixSuffixToArray(array, "", suffix);
  }

  public static String[] addPrefixSuffixToArray(String[] array, String prefix, String suffix) {
    String[] newArray;

    newArray = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      newArray[i] = (prefix == null ? "" : prefix) + array[i] + (suffix == null ? "" : suffix);
    }

    return newArray;
  }

  public static boolean isBimodal(double[] array) {
    return isBimodal(array, 0.5, 100);
  }

  public static boolean isBimodal(double[] array, double percentDropInPeak, int numBins) {
    return isBimodal(array, percentDropInPeak, (max(array) - min(array)) / numBins);
  }

  public static boolean isBimodal(double[] array, double percentDropInPeak, double binSize) {
    int numBins;
    int[] freqBinCounts, freqBinCountsSmooth;
    double minValue, maxFreq, localMinFreq;

    numBins = (int) ((ArrayUtils.max(array) - ArrayUtils.min(array)) / binSize);
    minValue = ArrayUtils.min(array);
    freqBinCounts = new int[numBins];
    for (int i = 0; i < array.length; i++) {
      freqBinCounts[(int) (Math.floor(array[i] - minValue) / binSize)]++;
    }

    // smoothing
    freqBinCountsSmooth = new int[freqBinCounts.length];
    for (int i = 1; i < numBins - 1; i++) {
      freqBinCountsSmooth[i] = (freqBinCounts[i - 1] + freqBinCounts[i] + freqBinCounts[i + 1]) / 3;
    }
    freqBinCountsSmooth[0] = (freqBinCounts[0] + freqBinCounts[1]) / 2;
    freqBinCountsSmooth[numBins - 1] = (freqBinCounts[numBins - 2] + freqBinCounts[numBins - 1])
                                       / 2;

    maxFreq = Double.NEGATIVE_INFINITY;
    localMinFreq = Double.POSITIVE_INFINITY;
    for (int i = 0; i < numBins; i++) {
      if (freqBinCountsSmooth[i] > maxFreq) {
        maxFreq = freqBinCountsSmooth[i];
      } else if (freqBinCountsSmooth[i] < (maxFreq * percentDropInPeak)
                 && freqBinCountsSmooth[i] < localMinFreq) {
        localMinFreq = freqBinCountsSmooth[i];
      } else if (freqBinCountsSmooth[i] >= (maxFreq * percentDropInPeak)) {
        return true;
      }
    }
    return false;
  }

  public static boolean isMultimodal(double[] array,
                                     double proportionOfLastPeakRequiredForNewLocalMinima,
                                     double proportionOfGlobalMaxRequiredForLocalMaxima,
                                     double binSize) {
    return getLocalModes(array, proportionOfLastPeakRequiredForNewLocalMinima,
                         proportionOfGlobalMaxRequiredForLocalMaxima, binSize, true).length > 1;
  }

  public static double[] getLocalModes(double[] array,
                                       double proportionOfLastPeakRequiredForNewLocalMinima,
                                       double proportionOfGlobalMaxRequiredForLocalMaxima) {
    return getLocalModes(array, proportionOfLastPeakRequiredForNewLocalMinima,
                         proportionOfGlobalMaxRequiredForLocalMaxima,
                         (max(array) - min(array)) / 40, true);
  }

  public static double[] getLocalModes(double[] array,
                                       double proportionOfLastPeakRequiredForNewLocalMinima,
                                       double proportionOfGlobalMaxRequiredForLocalMaxima,
                                       double binSize, boolean sensitiveToSmallNumbers) {
    int numBins;
    int[] freqBinCounts;
    double[] freqBinCountsSmooth;
    double minValue;
    int[] indicesOfLocalMaxima;
    double[] modes;

    numBins = (int) ((ArrayUtils.max(array) - ArrayUtils.min(array)) / binSize) + 1;
    minValue = ArrayUtils.min(array);
    freqBinCounts = new int[numBins];
    for (int i = 0; i < array.length; i++) {
      freqBinCounts[(int) Math.floor((array[i] - minValue) / binSize)]++;
    }

    // smoothing
    freqBinCountsSmooth = new double[freqBinCounts.length];
    for (int i = 1; i < numBins - 1; i++) {
      freqBinCountsSmooth[i] = (freqBinCounts[i - 1] + freqBinCounts[i] + freqBinCounts[i + 1]) / 3;
    }
    if (freqBinCounts.length >= 2) {
      freqBinCountsSmooth[0] = (freqBinCounts[0] + freqBinCounts[1]) / 2;
      freqBinCountsSmooth[numBins - 1] = (freqBinCounts[numBins - 2] + freqBinCounts[numBins - 1])
                                         / 2;
    }

    if (sensitiveToSmallNumbers) {
      proportionOfGlobalMaxRequiredForLocalMaxima = Math.max(proportionOfGlobalMaxRequiredForLocalMaxima,
                                                             Math.min(0.50,
                                                                      proportionOfGlobalMaxRequiredForLocalMaxima
                                                                            * proportionOfGlobalMaxRequiredForLocalMaxima
                                                                            * 300 / array.length));
      if (array.length < 50) {
        // System.out.println(array.length+"\t"+proportionOfGlobalMaxRequiredForLocalMaxima);
      }
    }

    indicesOfLocalMaxima = getIndicesOfLocalMaxima(freqBinCountsSmooth,
                                                   proportionOfLastPeakRequiredForNewLocalMinima,
                                                   proportionOfGlobalMaxRequiredForLocalMaxima);
    modes = new double[indicesOfLocalMaxima.length];
    for (int i = 0; i < modes.length; i++) {
      modes[i] = minValue + indicesOfLocalMaxima[i] * binSize + binSize / 2;
    }
    return modes;
  }

  public static double[] smooth(double[] array, int numOfPositionsInOneDirection) {
    double[] smoothed;

    smoothed = new double[array.length];

    return smoothed;
  }

  public static int[] getIndicesOfLocalMaxima(double[] array,
                                              double proportionOfLastPeakRequiredForNewLocalMinima,
                                              double proportionOfGlobalMaxRequiredForLocalMaxima) {
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
    for (int i = 0; i < array.length; i++) {
      if (array[i] > localMax) {
        localMax = array[i];
        indexOfLocalMax = i;
      }
      if (localMax >= globalMax * proportionOfGlobalMaxRequiredForLocalMaxima
          && array[i] <= (localMax * proportionOfLastPeakRequiredForNewLocalMinima)) {
        // System.out.println("localMax="+localMax+" at index "+indexOfLocalMax);
        localMin = array[i];
        indicesOfMaxima.add(indexOfLocalMax);
        localMax = Double.NEGATIVE_INFINITY;
      }
      if (array[i] >= (globalMax * proportionOfGlobalMaxRequiredForLocalMaxima)) {
        if (localMin != Double.POSITIVE_INFINITY) {
          // System.out.println("localMin="+localMin);
        }
        localMin = Double.POSITIVE_INFINITY;
      }
    }
    if (localMin != Double.POSITIVE_INFINITY
        && localMax >= globalMax * proportionOfGlobalMaxRequiredForLocalMaxima
        && indexOfLocalMax != -1) {
      // System.out.println("localMax="+localMax+" at index "+indexOfLocalMax);
      indicesOfMaxima.add(indexOfLocalMax);
    }

    return Ints.toArray(indicesOfMaxima);
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
      if (keeps[j]) {
        indices[i++] = j;
      }
    }
    return indices;
  }

  /**
   * Converts frequency counts into proportions. So, the input looks similar to this: FrequencyCount
   * --------------- Female 152 Male 148 The output looks like this: FrequencyCount ---------------
   * Female 50.67% Male 49.33%
   *
   * @param counts the frequency counts in array format
   * @return the corresponding proportion in array format
   */
  public static double[] getProportions(int[] counts) {
    int total = 0;
    double result[] = new double[counts.length];
    for (int count : counts) {
      total += count;
    }
    for (int i = 0; i < counts.length; i++) {
      result[i] = (double) counts[i] / (double) total;
    }
    return result;
  }

  /**
   * Creates an array of char and copies the contents of an array of String
   *
   * @param array array of String
   * @return an array of the converted String
   */
  public static char[] toCharArray(String[] array) {
    char[] newArray;

    newArray = new char[array.length];
    for (int i = 0; i < newArray.length; i++) {
      if (array[i].length() != 1) {
        System.err.println("Error - cannot convert string to char since it is longer than 1 byte: "
                           + array[i]);
      }
      newArray[i] = array[i].charAt(0);
    }

    return newArray;
  }

  /**
   * Transposes a List<List<>> i.e. a 2 dimensional list
   *
   * @param table : the 2D list to be transposed
   * @param <T> : generic template
   * @param table : the 2D list to be transposed
   * @param <T> : generic template
   * @return a {@link List<List>>} which is transposed
   */
  public static <T> List<List<T>> transpose(List<List<T>> table) {
    List<List<T>> ret = new ArrayList<>();
    final int N = table.get(0).size();
    for (int i = 0; i < N; i++) {
      List<T> col = new ArrayList<>();
      for (List<T> row : table) {
        col.add(row.get(i));
      }
      ret.add(col);
    }
    return ret;
  }

  public static double[] concatDubs(double[] first, double[] other) {
    int totalLength = first.length + other.length;

    double[] result = new double[totalLength];
    int index = 0;
    for (double element : first) {
      result[index] = element;
      index++;
    }
    for (double element : other) {
      result[index] = element;
      index++;
    }
    return result;
  }

  /**
   * Remove common elements from the front and/or back of strings Ex [acatcat3cat, acatcat4cat] ->
   * [3,4] <br>
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
          if (startIndex >= toUniq[i].length()
              || !(toUniq[i].charAt(startIndex) + "").equals(start)) {
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
          if (toUniq[i].length() <= stopIndex
              || !(toUniq[i].charAt(toUniq[i].length() - stopIndex) + "").equals(stop)) {
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
   * @param front tag on to the front of each entry
   * @param back tag on to the back of each entry
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
   * @param first the first array
   * @param rest rest all arrays
   * @param <T> generic: can take any type of object
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
      res[i] = i + 1;
    }

    done = false;
    v = new Vector<>();
    indices = new int[r];
    while (!done) {
      for (int i = 0; i < r; i++) {
        indices[i] = res[i] - 1;
      }
      v.add(indices.clone());
      done = getNext(res, n, r);
    }

    return Matrix.toMatrix(v);
  }

  public static boolean[] getFinite(double[] arr) {
    boolean[] fin = new boolean[arr.length];
    for (int i = 0; i < arr.length; i++) {
      fin[i] = Double.isFinite(arr[i]);
    }
    return fin;
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
    return ProbDist.ChiDistReverse(ArrayUtils.median(pvals), 1) / ProbDist.ChiDistReverse(0.50, 1);
  }

  /**
   * Takes the log base 2 of every element in the array a
   *
   * @param a an array of doubles
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
      if (!ext.isValidDouble(array[i] + "")) {
        return true;
      }
    }
    return false;
  }

  public static String[] removeMissingValues(String[] array) {
    ArrayList<String> valid = new ArrayList<>();
    for (String str : array) {
      if (str == null || !ext.isValidDouble(str)) {
        continue;
      }
      valid.add(str);
    }
    return valid.toArray(new String[valid.size()]);
  }

  public static <T> T[] combine(T[] array1, T[] array2) {
    T[] newArray = Arrays.copyOf(array1, array1.length + array2.length);
    for (int i = array1.length; i < array1.length + array2.length; i++) {
      newArray[i] = array2[i - array1.length];
    }
    return newArray;
  }

  /**
   * Appends the second matrix to the starting matrix Assumes the second matrix has at least as many
   * rows as the starting matrix Will not add rows to accommodate a longer second matrix
   * 
   * @param start Starting matrix
   * @param toAppend Matrix to append
   * @return
   */
  public static String[][] append(String[][] start, String[][] toAppend) {
    String[][] m = new String[start.length][start[0].length + toAppend[0].length];

    for (int i = 0; i < m.length; i++) {
      for (int j = 0; j < m[0].length; j++) {
        if (start[0].length <= j)
          m[i][j] = toAppend[i][j - start[0].length];
        else
          m[i][j] = start[i][j];
      }
    }

    return m;
  }

  /**
   * Helper method to avoid modification of original data structure
   *
   * @return a sorted copy of the given array
   */
  public static double[] sortedCopy(double[] array) {
    double[] sorted = Arrays.copyOf(array, array.length);
    Arrays.sort(sorted);
    return sorted;
  }

  /**
   * As {@link ArrayUtils#sortedCopy(double[])} for float arrays.
   */
  public static float[] sortedCopy(float[] array) {
    float[] sorted = Arrays.copyOf(array, array.length);
    Arrays.sort(sorted);
    return sorted;
  }

  /**
   * As {@link ArrayUtils#sortedCopy(double[])} for object arrays.
   */
  public static <T> T[] sortedCopy(T[] array) {
    T[] sorted = Arrays.copyOf(array, array.length);
    Arrays.sort(sorted);
    return sorted;
  }

  public static String[] sortedCopyAlphanum(String[] array) {
    String[] sorted = Arrays.copyOf(array, array.length);
    Arrays.sort(sorted, new Comparator<String>() {

      final Pattern p = Pattern.compile("^\\d+");

      @Override
      public int compare(String object1, String object2) {
        Matcher m = p.matcher(object1);
        Integer number1 = null;
        if (!m.find()) {
          return object1.compareToIgnoreCase(object2);
        } else {
          Integer number2 = null;
          number1 = Integer.parseInt(m.group());
          m = p.matcher(object2);
          if (!m.find()) {
            return object1.compareToIgnoreCase(object2);
          } else {
            number2 = Integer.parseInt(m.group());
            int comparison = number1.compareTo(number2);
            if (comparison != 0) {
              return comparison;
            } else {
              return object1.compareToIgnoreCase(object2);
            }
          }
        }
      }
    });
    return sorted;
  }

  /**
   * @param source to copy
   * @return a deep copy of the input array, copying only the int values to new 1st and 2nd
   *         dimension arrays
   */
  public static int[][] deepCopy(int[][] source) {
    int[][] copy = new int[source.length][];
    for (int i = 0; i < source.length; i++) {
      copy[i] = source[i].clone();
    }
    return copy;
  }

  /**
   * @param source to copy
   * @return a deep copy of the input array, copying only the int values to new 1st and 2nd
   *         dimension arrays
   */
  public static int[][][] deepCopy(int[][][] source) {
    int[][][] copy = new int[source.length][][];
    for (int i = 0; i < source.length; i++) {
      copy[i] = deepCopy(source[i]);
    }
    return copy;
  }

  /**
   * @param source to copy
   * @return a deep copy of the input array, copying only the double values to new 1st and 2nd
   *         dimension arrays
   */
  public static double[][] deepCopy(double[][] source) {
    double[][] copy = new double[source.length][];
    for (int i = 0; i < source.length; i++) {
      copy[i] = source[i].clone();
    }
    return copy;
  }

  /**
   * @param source to copy
   * @return a deep copy of the input array, copying only the double values to new 1st and 2nd
   *         dimension arrays
   */
  public static double[][][] deepCopy(double[][][] source) {
    double[][][] copy = new double[source.length][][];
    for (int i = 0; i < source.length; i++) {
      copy[i] = deepCopy(source[i]);
    }
    return copy;
  }

  public static int[] booleanArrayRunLengths(boolean[] samplesToLoadAfterIndex) {
    ArrayList<Integer> runs = new ArrayList<>();
    int run = 0;
    for (boolean b : samplesToLoadAfterIndex) {
      if (b) {
        run++;
      } else {
        runs.add(run);
        run = 0;
      }
    }
    runs.add(run);
    return Ints.toArray(runs);
  }

  /**
   * Preserves input array iteration order
   * 
   * @param values
   * @return
   */
  public static <T> ImmutableMap<T, Integer> immutableIndexMap(T[] values) {
    ImmutableMap.Builder<T, Integer> builder = new ImmutableMap.Builder<>();
    for (int i = 0; i < values.length; i++) {
      builder.put(values[i], i);
    }
    return builder.build();
  }

  public static <T> Map<T, Integer> indexMap(T[] values) {
    Map<T, Integer> map = new HashMap<>();
    for (int i = 0; i < values.length; i++) {
      map.put(values[i], i);
    }
    return map;
  }

  public static float[][] computeZScores(float[][] regionMedianValues) {
    float[][] scores = new float[regionMedianValues.length][];
    for (int i = 0; i < regionMedianValues.length; i++) {
      scores[i] = computeZScores(regionMedianValues[i]);
    }
    return scores;
  }

  public static float[] computeZScores(float[] regionMedianValues) {
    float mean = mean(regionMedianValues);
    float sd = stdev(regionMedianValues, false);
    float[] scores = new float[regionMedianValues.length];
    for (int i = 0; i < regionMedianValues.length; i++) {
      scores[i] = (regionMedianValues[i] - mean) / sd;
    }
    return scores;
  }

}
