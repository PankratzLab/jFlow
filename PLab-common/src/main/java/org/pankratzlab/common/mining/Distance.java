package org.pankratzlab.common.mining;

import org.pankratzlab.common.ArrayUtils;

public class Distance {

  public static final double euclidean(int[] p1, int[] p2) {
    return euclidean(ArrayUtils.toDoubleArray(p1), ArrayUtils.toDoubleArray(p2));
  }

  public static final double euclidean(double[] p1, double[] p2) {
    double dist = 0;

    if (p1.length != p2.length) {
      System.err.println("Error - points have different numbers of dimensions");
    }

    for (int i = 0; i < p1.length; i++) {
      dist += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }

    return Math.sqrt(dist);
  }
}