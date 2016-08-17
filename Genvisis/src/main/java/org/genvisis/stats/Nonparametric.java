package org.genvisis.stats;

import org.genvisis.common.Array;
import org.genvisis.common.CountHash;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.mining.Transformations;

public class Nonparametric {


  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "Nonparametric.dat";
    String logfile = null;
    Logger log;

    String usage = "\n" + "stats.Nonparametric requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      System.out.println(runWilcoxonRankSumTest(
          Array.toDoubleArray(
              HashVec.loadFileToStringArray("file1.txt", false, new int[] {0}, false)),
          Array.toDoubleArray(
              HashVec.loadFileToStringArray("file2.txt", false, new int[] {0}, false)),
          log));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  // non paramteric form of ANOVA, one-way analysis of variance
  public static double[] runKruskalWallis() {
    // TODO Auto-generated method stub

    return null;
  }

  // nonparametric form of the independent T-Test, does not assume the traits are normally
  // distributed
  public static double runWilcoxonRankSumTest(double[] group1, double[] group2, Logger log) { // aka
                                                                                              // Mann-Whitney
                                                                                              // U
                                                                                              // test,
                                                                                              // Mann-Whitney-Wilcoxon
                                                                                              // test,
                                                                                              // Wilcoxon-MannWhitney
                                                                                              // test
    double[] merged;
    double[] ranks, ranks1, ranks2;
    int n1, n2, n;
    double u, mu, su, z, p;
    CountHash ch;
    double extra;
    int[] counts;

    n1 = group1.length;
    n2 = group2.length;
    n = n1 + n2;

    ch = new CountHash();
    merged = new double[n];
    for (int i = 0; i < n; i++) {
      if (i < n1) {
        merged[i] = group1[i];
      } else {
        merged[i] = group2[i - n1];
      }
      ch.add(merged[i] + "");
    }

    ranks = Transformations.rankTransform(merged);
    ranks1 = Array.subArray(ranks, 0, n1);
    ranks2 = Array.subArray(ranks, n1);

    // System.out.println("mean rank for group1: "+Array.mean(ranks1));
    // System.out.println("mean rank for group2: "+Array.mean(ranks2));
    //
    // System.out.println("rank sum for group1: "+Array.sum(ranks1));
    // System.out.println("rank sum for group2: "+Array.sum(ranks2));

    u = Array.sum(ranks1) - n1 * (n1 + 1) / 2;
    if (u < n1 * n2 - u) {
      ranks = ranks1;
    } else {
      u = n1 * n2 - u;
      ranks = ranks2;
    }

    // System.out.println("U: "+u);

    // mu = n1*(n1+n2+1)/2; // wikipedia entry is wrong
    mu = n1 * n2 / 2;

    // su = Math.sqrt(n1*n2*(n1+n2+1)/12); // assumes there are no ties
    // System.out.println("stdev1: "+su);
    // z = (u - mu) / su;
    // System.out.println("z1: "+z);
    // p = ProbDist.NormDist(z);
    // System.out.println("p1: "+p);

    extra = 0;
    counts = ch.getCounts();
    for (int count : counts) {
      if (count > 1) {
        extra += (Math.pow(count, 3) - count) / 12;
      }
    }

    // System.out.println("extra tail: "+extra);

    extra = (Math.pow(n, 3) - n) / 12 - extra;

    // System.out.println("extra: "+extra);

    su = Math.sqrt(extra * n1 * n2 / (n * (n - 1))); // assumes there are no ties

    // System.out.println("stdev2: "+su);

    z = (u - mu) / su;

    // System.out.println("z2: "+z);

    p = ProbDist.NormDist(z);

    // System.out.println("p2: "+p);

    return p;
  }

  public static double runWilcoxonRankSumTest(int[] groupings, double[] values, Logger log) {
    return runWilcoxonRankSumTest(Ttest.splitOut(groupings, values, Array.min(groupings)),
        Ttest.splitOut(groupings, values, Array.min(groupings) + 1), log);
  }

  // nonparametric form of the Paired T-Test, does not assume the traits are normally distributed
  public static double[] runWilcoxonSignedRankTest() {
    // TODO Auto-generated method stub

    return null;
  }
}
