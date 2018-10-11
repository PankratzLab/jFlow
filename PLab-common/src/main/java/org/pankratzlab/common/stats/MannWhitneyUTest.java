package org.pankratzlab.common.stats;

import java.util.Arrays;
import java.util.List;

public class MannWhitneyUTest {

  private static final double TIE_SCORE = 0.5;
  private static final double WIN_SCORE = 1;

  public static double mannWhitneyU(List<Double> pheno, List<Integer> aIndices,
                                    List<Integer> bIndices) {
    double[] a = new double[aIndices.size()];
    double[] b = new double[bIndices.size()];

    int index = 0;
    for (Integer k : aIndices) {
      a[index++] = pheno.get(k);
    }
    index = 0;
    for (Integer k : bIndices) {
      b[index++] = pheno.get(k);
    }

    return mannWhitneyU(a, b);
  }

  public static double mannWhitneyU(double[] a, double[] b) {
    Arrays.sort(a);
    Arrays.sort(b);
    return mannWhitneyUSorted(a, b);
  }

  public static double mannWhitneyUSorted(double[] aSorted, double[] bSorted) {
    double u1 = mannWhitneyUHelper(aSorted, bSorted);
    double u2 = (aSorted.length * bSorted.length) - u1;
    return Math.max(u1, u2);
  }

  /**
   * Compute the Mann-Whitney U1 score
   */
  private static double mannWhitneyUHelper(double[] aSorted, double[] bSorted) {
    double aSum = 0.0;

    int aIndex = 0;
    int bIndex = 0;

    // We want to sum the rank of all the A's
    while (aIndex < aSorted.length) {
      double a = aSorted[aIndex];

      int c = 0;
      int equalA = 0;
      int equalB = 0;

      // If we haven't exceeded all the B's, we have to see if A or B is smaller
      if (bIndex < bSorted.length) {
        double b = bSorted[bIndex];

        c = Double.compare(a, b);

        if (c >= 0) {
          // We're going to consume from the B array, so we want to know how many sequential 
          // elements of B have this same value
          equalB = countEqual(bSorted, bIndex);
        }
      }

      if (c <= 0) {
        // We're going to consume from the B array, either because all the B's are consumed or the
        // leading A(s) is(are) smaller than the leading B(s).
        // So we count how many A's have the same current value
        equalA = countEqual(aSorted, aIndex);

        // We're scoring rank so we have to consider how many elements are we consuming from the A
        // and B arrays, and then calculate the average of their sequential rankings.
        // e.g. if the first two A's have the same value as the first three B's, then the rankings
        // are 1, 2, 3, 4, 5
        int combinedCount = equalA + equalB;
        int combinedRank = 0;
        for (int i = 0; i < combinedCount; i++) {
          combinedRank += (1 + aIndex + bIndex + i);
        }

        // The score is the average of the combined ranks
        double combinedScore = ((double) combinedRank) / combinedCount;

        // And the A sum gets shares of this score equal to the number of A values that were used
        aSum += (equalA * combinedScore);
      }

      // Finally, we advance the A and B positions for the values consumed tihs pass
      aIndex += equalA;
      bIndex += equalB;
    }

    double u = aSum - (((aSorted.length) * (aSorted.length + 1)) / 2);

    return u;
  }

  /**
   * @return The number of leading items in the given array with the same value
   */
  private static int countEqual(double[] arr, int start) {
    int i = start + 1;
    for (; i < arr.length; i++) {
      if (Double.compare(arr[i], arr[start]) != 0) {
        break;
      }
    }
    return i - start;
  }

  /**
   * Implementation of method 1 for computing U
   */
  private static double mannWhitneyUMethod1(double[] aSorted, double[] bSorted) {
    double score = 0;
    int aIndex = 0;
    int bIndex = 0;

    while (aIndex < aSorted.length) {
      double aTest = aSorted[aIndex];

      // Advance through B values that are smaller than the current a value
      while (bIndex < bSorted.length && Double.compare(aTest, bSorted[bIndex]) > 0) {
        bIndex++;
      }

      // If there are no B's >= the current A, then all remaining A's scored 0
      if (bIndex >= bSorted.length) {
        break;
      }

      // Count how many B's are equal to the current A and update score appropriately
      int bEqualIndex = bIndex;
      while (bEqualIndex < bSorted.length && Double.compare(aTest, bSorted[bEqualIndex]) == 0) {
        bEqualIndex++;
      }

      score += (TIE_SCORE * (bEqualIndex - bIndex));

      // All remaining B's that weren't ties were beaten by the current A, so we update the score
      score += (WIN_SCORE * (bSorted.length - bEqualIndex));

      aIndex++;
    }

    return score;
  }

  public static double calculatePValue(double Umin, int n1, int n2) {
    long prod = (long) n1 * n2;

    // http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U#Normal_approximation
    double mean = prod / 2.0;
    double se = prod * (n1 + n2 + 1) / 12.0;

    // assumes number of ties is small
    double z = (Umin - mean) / Math.sqrt(se);

    return ProbDist.NormDist(z);
  }
}
