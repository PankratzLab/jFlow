package org.genvisis.common;

import java.util.Vector;

import com.google.common.primitives.Ints;

/**
 * Static utility class for working with {@link PrimitiveVector}s.
 */
public final class Vectors {

  /**
   * Helper method to ensure an order array is the same size as a given vector.
   */
  private static boolean checkLength(Vector<?> v, int[] order) {
    if (v.size() != order.length) {
      Logger log = new Logger();
      log.reportError("Error - order array does not have the same number of elements ("
                      + order.length + ") as the target (" + v.size() + ")");
      return false;
    }
    return true;
  }

  /**
   * Create an array of initialized, empty vectors of the given type.
   *
   * @param type {@link PrimitiveVector} to use as the array type.
   * @param numVectors size of the array
   * @return An array of vectors, with each element initialized to an empty vector.
   */
  public static <T extends PrimitiveVector> T[] initializedArray(Class<T> type, int numVectors) {
    Logger log = new Logger();
    @SuppressWarnings("unchecked")
    T[] vectors = (T[]) java.lang.reflect.Array.newInstance(type, numVectors);
    for (int i = 0; i < vectors.length; i++) {
      try {
        vectors[i] = type.newInstance();
      } catch (InstantiationException e) {
        log.reportException(e);
      } catch (IllegalAccessException e) {
        log.reportException(e);
      }
    }
    return vectors;
  }

  /**
   * Convert an {@link IntVector} to an {@code int[]}, using the values of the {@code order} array
   * as indices into the vector.
   * <p>
   * FIXME: this will need to be rewritten for every primitive type. Instead of using primitive
   * arrays directly, we should use array classes directly backed by primitive arrays.
   * </p>
   *
   * @param v Base vector to convert
   * @param order Each element {@code i} of the output array will have the value
   *        {@code Vector.get(order[i])}.
   * @return An ordered array, or {@code null} if the order array is not the same size as the base
   *         vector.
   */
  public static int[] orderedArray(IntVector v, int[] order) {
    if (!checkLength(v, order)) {
      return null;
    }

    int[] a = new int[order.length];

    for (int i = 0; i < a.length; i++) {
      a[i] = v.get(order[i]);
    }

    return a;
  }

  /**
   * Convert an array of vectors to a 2-dimensional primitive array
   * <p>
   * FIXME: this will need to be rewritten for every primitive type. Instead of using primitive
   * arrays directly, we should use array classes directly backed by primitive arrays.
   * </p>
   *
   * @param va Array of vectors
   * @return A 2-D array representation of the given vector array.
   */
  public static int[][] toMatrix(IntVector[] va) {
    int[][] matrix = new int[va.length][];
    for (int i = 0; i < matrix.length; i++) {
      matrix[i] = Ints.toArray(va[i]);
    }
    return matrix;
  }

  private Vectors() {
    // prevent instantiation of static utility class
  }
}
