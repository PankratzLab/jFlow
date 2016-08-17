package org.genvisis.stats;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.common.Array;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;

public class Stats {

  public static double cdf(OpdfGaussian dist, ObservationReal n) {
    return new NormalDistribution(dist.mean(), dist.variance()).cumulativeProbability(n.value);
  }

  public static double FishersExact(double a, double b, double c, double d, boolean oneTailed) {
    double p = 0;
    double[] LogFactorial = new double[(int) (a + b + c + d + 1)];
    LogFactorial[0] = 0;

    for (int i = 1; i <= a + b + c + d; i++) {
      LogFactorial[i] = LogFactorial[i - 1] + Math.log(i);
    }
    if (a * d > b * c) {
      a = a + b;
      b = a - b;
      a = a - b;

      c = c + d;
      d = c - d;
      c = c - d;
    }

    if (a > d) {
      a = a + d;
      d = a - d;
      a = a - d;
    }
    if (b > c) {
      b = b + c;
      c = b - c;
      b = b - c;
    }

    double a_org = a;
    double p_sum = 0;
    double p_1 = p = FisherSub(LogFactorial, (int) a, (int) b, (int) c, (int) d);

    while (a >= 0) {
      p_sum += p;
      if (a == 0) {
        break;
      }
      a--;
      b++;
      c++;
      d--;
      p = FisherSub(LogFactorial, (int) a, (int) b, (int) c, (int) d);
    }

    if (oneTailed) {
      p_sum = 1e-8 * Math.round(p_sum * 1e8);
      return p_sum;
    } else {
      a = b;
      b = 0;
      c = c - a;
      d = d + a;

      p = FisherSub(LogFactorial, (int) a, (int) b, (int) c, (int) d);
      while (p < p_1) {
        if (a == a_org) {
          break;
        }
        p_sum += p;
        a--;
        b++;
        c++;
        d--;
        p = FisherSub(LogFactorial, (int) a, (int) b, (int) c, (int) d);
      }
      p_sum = 1e-8 * Math.round(p_sum * 1e8);
      return p_sum;
    }
  }

  public static double FisherSub(double[] LogFactorial, int a, int b, int c, int d) {
    double p = 0;

    p = LogFactorial[a + b] + LogFactorial[c + d] + LogFactorial[a + c] + LogFactorial[b + d]
        - LogFactorial[a + b + c + d] - LogFactorial[a] - LogFactorial[b] - LogFactorial[c]
        - LogFactorial[d];

    return Math.exp(p);
  }

  public static double PearsonChiSquare(double a, double b, double c, double d) {
    double pr22 = (a + c) / (a + b + c + d);
    double aExp = pr22 * (a + b);
    double bExp = a + b - aExp;
    double cExp = pr22 * (c + d);
    double dExp = c + d - cExp;

    return Math.pow(a - aExp, 2.0) / aExp + Math.pow(b - bExp, 2.0) / bExp
           + Math.pow(c - cExp, 2.0) / cExp + Math.pow(d - dExp, 2.0) / dExp;
  }

  public static double ttestOneSample(double mean, double stdev, int n, double expected) {
    return (mean - expected) / (stdev / Math.sqrt(n));
  }

  public static double ttestOneSample(double[] array, double expected) {
    return ttestOneSample(Array.mean(array), Array.stdev(array), array.length, expected);
  }

  public static double ztest(double pHat, double stdev, double expected) {
    return (pHat - expected) / stdev;
  }

  public static double ztest(double p1, int n1, double p2, int n2) {
    double phat, sp, z;

    phat = (p1 + p2) / 2;
    sp = Math.sqrt(phat * (1 - phat) * ((1 / (double) n1) + (1 / (double) n2)));
    z = (p1 - p2) / sp;

    return z;
  }

  public static double ztest(int x1, int n1, int x2, int n2) {
    double phat, sp, z;

    phat = (double) (x1 + x2) / (double) (n1 + n2);
    sp = Math.sqrt(phat * (1 - phat) * ((1 / (double) n1) + (1 / (double) n2)));
    z = (((double) x1 / (double) n1) - ((double) x2 / (double) n2)) / sp;

    return z;
  }

  /**
   * The code has been moved to common.Array
   */
  // public static double kurtosis(double[] array) {
  // double kurt = -1;
  // double mean = Array.mean(array);
  // double m2, s, m4s;
  // double n = array.length;
  //
  // m2 = 0;
  // for (int i = 0; i<n; i++) {
  // m2 += Math.pow(array[i]-mean, 2);
  // }
  // m2 /= (n-1);
  // s = Math.sqrt(m2);
  //
  // m4s = 0;
  // for (int i = 0; i<array.length; i++) {
  // m4s += Math.pow((array[i]-mean)/s, 4);
  // }
  //
  // kurt = n*(n+1)/((n-1)*(n-2)*(n-3))*m4s-3*Math.pow(n-1, 2)/((n-2)*(n-3));
  //
  // return kurt;
  // }

  // alternate method
  // private static double kurtosisCloseButNoCigar(double[] array) {
  // double kurt = -1;
  // double mean = Array.mean(array);
  // double m4, m2;
  //
  // m4 = 0;
  // for (int i = 0; i<array.length; i++) {
  // m4 += Math.pow(array[i]-mean, 4);
  // }
  // m4 /= array.length;
  //
  // m2 = 0;
  // for (int i = 0; i<array.length; i++) {
  // m2 += Math.pow(array[i]-mean, 2);
  // }
  // m2 /= array.length;
  //
  // kurt = m4 / Math.pow(m2, 2) - 3;
  //
  // return kurt;
  // }
}
