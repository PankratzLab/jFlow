package org.genvisis.stats;

import static org.junit.Assert.assertEquals;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import org.junit.Test;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.shared.stats.MannWhitneyUTest;
import org.pankratzlab.shared.stats.ProbDist;
import com.google.common.collect.ImmutableList;

public class testMeanWhitney {

  @Test
  public void testSimple() {
    assertEquals(25.0,
                 MannWhitneyUTest.mannWhitneyUSorted(new double[] {1.0, 3.0, 3.0, 3.0, 3.0, 3.0},
                                                     new double[] {2.0, 2.0, 2.0, 2.0, 2.0, 4.0}),
                 0.000005);
  }

  public void testIndexed() {
    List<Double> v = ImmutableList.of(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    List<Integer> i1 = ImmutableList.of(0, 2, 4);
    List<Integer> i2 = ImmutableList.of(1, 3, 5);
    MannWhitneyUTest.mannWhitneyU(v, i1, i2);
  }

  public void testComplex() {
    List<Double> pheno = new ArrayList<>();
    List<Integer> aIndices = new ArrayList<>();
    List<Integer> bIndices = new ArrayList<>();

    Random rand = new Random(12345678);

    for (int i = 0; i < 600; i++) {
      pheno.add(rand.nextDouble());
      if (rand.nextBoolean()) {
        aIndices.add(i);
      } else {
        bIndices.add(i);
      }
    }

    long t = System.currentTimeMillis();

    for (int i = 0; i < 100; i++) {
      Collections.shuffle(pheno, rand);

      // run permutation for this shuffled phenotype value list
      for (int j = 0; j < 14000; j++) {
        MannWhitneyMR.mannWhitneyU(pheno, aIndices, bIndices);
        //                MannWhitneyUStat.mannWhitneyU(pheno, aIndices, bIndices);
      }
    }

    System.out.println(System.currentTimeMillis() - t);
  }

  private static class MannWhitneyMR {

    public static double mannWhitneyU(double[] a, double[] b) {
      double u = 0.0;

      Arrays.sort(a);
      Arrays.sort(b);

      double[] rankA = new double[a.length];
      double[] rankB = new double[b.length];

      int rank = 1;
      int bIndex = 0;
      int aIndex = 0;
      while (aIndex < a.length || bIndex < b.length) {
        if (bIndex >= b.length || (aIndex < a.length && a[aIndex] < b[bIndex])) {
          double rankAvg = rank++;
          int endA = aIndex + 1;
          while (endA < a.length && a[aIndex] == a[endA]) {
            rankAvg += rank;
            rank++;
            endA++;
          }
          rankAvg /= endA - aIndex;
          for (int i = aIndex; i < endA; i++) {
            rankA[i] = rankAvg;
          }
          aIndex = endA;
        } else if (aIndex >= a.length || b[bIndex] < a[aIndex]) {
          double rankAvg = rank++;
          int endB = bIndex + 1;
          while (endB < b.length && b[bIndex] == b[endB]) {
            rankAvg += rank;
            rank++;
            endB++;
          }
          rankAvg /= endB - bIndex;
          for (int i = bIndex; i < endB; i++) {
            rankB[i] = rankAvg;
          }
          bIndex = endB;
        } else if (a[aIndex] == b[bIndex]) {
          double rankAvg = rank + ++rank;
          rank++;
          int endA = aIndex + 1;
          int endB = bIndex + 1;
          int countTies = 2;
          while (endA < a.length && a[aIndex] == a[endA]) {
            rankAvg += rank;
            rank++;
            countTies++;
            endA++;
          }
          while (endB < b.length && b[bIndex] == b[endB]) {
            rankAvg += rank;
            rank++;
            countTies++;
            endB++;
          }
          rankAvg /= countTies;

          for (int i = aIndex; i < endA; i++) {
            rankA[i] = rankAvg;
          }
          for (int i = bIndex; i < endB; i++) {
            rankB[i] = rankAvg;
          }

          aIndex = endA;
          bIndex = endB;
        }
      }
      double u1 = ArrayUtils.sum(rankA) - (a.length * (a.length + 1) / 2);
      double u2 = ArrayUtils.sum(rankB) - (b.length * (b.length + 1) / 2);

      return Math.max(u1, u2);
    }

    public static void mannWhitneyU(List<Double> pheno, List<Integer> aIndices,
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

      MannWhitneyMR.mannWhitneyU(a, b);
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

}
