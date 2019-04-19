package org.pankratzlab.common;

import java.io.IOException;

import org.pankratzlab.common.stats.ContingencyTable;
import org.pankratzlab.common.stats.Maths;
import org.pankratzlab.common.stats.ProbDist;

public class AlleleFreq {

  public static double calcFrequency(int[] genotypes) {
    return calcFrequency(ArrayUtils.toDoubleArray(genotypes));
  }

  public static double calcFrequency(double[] genotypes) {
    if (genotypes.length != 3) {
      System.err.println("Error - can't compute an allele frequency from " + genotypes.length
                         + " classes of genotypes");
    }
    return calcFrequency(genotypes[0], genotypes[1], genotypes[2]);
  }

  public static double calcFrequency(double pp, double pq, double qq) {
    double total = pp + pq + qq;

    if (total == 0) {
      return Double.NaN;
    } else {
      return (pp * 2 + pq) / (2 * total);
    }
  }

  public static double calcMAF(int pp, int pq, int qq) {
    double total = pp + pq + qq;
    if (total == 0) {
      return Double.NaN;
    } else {
      double minor = Maths.min(pp, qq);
      return (minor * 2 + pq) / (2 * total);
    }
  }

  public static double computeHeterozygosity(int[] counts) {
    double p = (double) (counts[0] * 2 + counts[1]) / (double) (ArrayUtils.sum(counts) * 2);
    return 1 - p * p - (1 - p) * (1 - p);
  }

  public static double HWE(int pp, int pq, int qq) {
    return HWE((double) pp, (double) pq, (double) qq);
  }

  public static double HWE(int pp, int pq, int pr, int qq, int qr, int rr) {
    return HWE((double) pp, (double) pq, (double) pr, (double) qq, (double) qr, (double) rr);
  }

  public static double[][] getObsExp(double pp, double pq, double qq) {
    double total = pp + pq + qq;
    double p = (pp * 2 + pq) / (2 * total);
    double q = (qq * 2 + pq) / (2 * total);
    double[] observed = new double[3];
    double[] expected = new double[3];

    observed[0] = pp;
    observed[1] = pq;
    observed[2] = qq;

    expected[0] = p * p * total;
    expected[1] = 2 * p * q * total;
    expected[2] = q * q * total;

    return new double[][] {observed, expected};
  }

  public static double[][] getObsExp(double pp, double pq, double pr, double qq, double qr,
                                     double rr) {
    double total = pp + pq + pr + qq + qr + rr;
    double p = (pp * 2 + pq + pr) / (2 * total);
    double q = (qq * 2 + pq + qr) / (2 * total);
    double r = (rr * 2 + pr + qr) / (2 * total);
    double[] observed = new double[6];
    double[] expected = new double[6];

    observed[0] = pp;
    observed[1] = pq;
    observed[2] = pr;
    observed[3] = qq;
    observed[4] = qr;
    observed[5] = rr;

    expected[0] = p * p * total;
    expected[1] = 2 * p * q * total;
    expected[2] = 2 * p * r * total;
    expected[3] = q * q * total;
    expected[4] = 2 * q * r * total;
    expected[5] = r * r * total;

    return new double[][] {observed, expected};
  }

  public static double HWE(double pp, double pq, double qq) {
    return HWE(getObsExp(pp, pq, qq));
  }

  public static double HWE(double pp, double pq, double pr, double qq, double qr, double rr) {
    return HWE(getObsExp(pp, pq, pr, qq, qr, rr));
  }

  public static double HWE(double[][] obsExp) {
    double chi;
    chi = 0;
    for (int i = 0; i < obsExp[0].length; i++) {
      chi += (obsExp[0][i] - obsExp[1][i]) * (obsExp[0][i] - obsExp[1][i]) / obsExp[1][i];
    }

    return chi;
  }

  public static double HWEsig(double[] genotypes) {
    if (genotypes.length == 3) {
      return HWEsig(genotypes[0], genotypes[1], genotypes[2]);
    }
    if (genotypes.length == 6) {
      return HWEsig(genotypes[0], genotypes[1], genotypes[2], genotypes[3], genotypes[4],
                    genotypes[5]);
    }
    System.err.println("Error - can't compute Hardy Weinberg from " + genotypes.length
                       + " classes of genotypes");
    return Double.NaN;
  }

  public static double HWEsig(int[] genotypes) {
    return HWEsig(ArrayUtils.toDoubleArray(genotypes));
  }

  public static double HWEsig(int pp, int pq, int qq) {
    return ProbDist.ChiDist(HWE((double) pp, (double) pq, (double) qq), 1);
  }

  public static double HWEsig(int pp, int pq, int pr, int qq, int qr, int rr) {
    return ProbDist.ChiDist(HWE((double) pp, (double) pq, (double) pr, (double) qq, (double) qr,
                                (double) rr),
                            1);
  }

  public static double HWEsig(double pp, double pq, double qq) {
    return ProbDist.ChiDist(HWE(pp, pq, qq), 1);
  }

  public static double HWEsig(double pp, double pq, double pr, double qq, double qr, double rr) {
    return ProbDist.ChiDist(HWE(pp, pq, pr, qq, qr, rr), 3);
  }

  public static double[] HetExcess(int pp, int pq, int qq) {
    return HetExcess((double) pp, (double) pq, (double) qq);
  }

  public static double[] HetExcess(double pp, double pq, double qq) {
    double p, q;
    double total = pp + pq + qq;
    double[][] calculations = new double[2][3];
    double[] results = new double[2];

    calculations[0][0] = pp / total * 100;
    calculations[0][1] = pq / total * 100;
    calculations[0][2] = qq / total * 100;

    p = (pp * 2 + pq) / (2 * total);
    q = (qq * 2 + pq) / (2 * total);

    calculations[1][0] = p * p * 100;
    calculations[1][1] = 2 * p * q * 100;
    calculations[1][2] = q * q * 100;

    for (int i = 0; i < 3; i++) {
      results[0] += Math.abs(calculations[0][i] - calculations[1][i]);
    }
    results[0] /= 100;
    results[1] = ProbDist.ChiDist(ContingencyTable.ChiSquare(calculations, false, false), 1);

    return results;
  }

  public static void main(String[] args) throws IOException {
    System.err.println("This is much more useful when called from a program...");
  }

  public static String getHWETableHTML(int[] values) {
    if (values.length != 3 && values.length != 6) {
      System.err.println("Error - can't compute Hardy Weinberg from " + values.length
                         + " classes of alleleCounts");
      return null;
    }
    double[][] vals = values.length == 3 ? getObsExp(values[0], values[1], values[2])
                                         : getObsExp(values[0], values[1], values[2], values[3],
                                                     values[4], values[5]);
    StringBuilder output = new StringBuilder("<html><table border=\"1\"><tr><td></td>");
    String[] colLabels = values.length == 3 ? new String[] {"PP", "PQ", "QQ"}
                                            : new String[] {"PP", "PQ", "PR", "QQ", "QR", "RR"};

    for (String col : colLabels) {
      output.append("<td>");
      output.append(col);
      output.append("</td>");
    }
    output.append("</tr>");

    output.append("<tr>");
    output.append("<td>Obs:</td>");
    for (double v : vals[0]) {
      output.append("<td align=\"center\">");
      output.append(v);
      output.append("</td>");
    }
    output.append("</tr>");
    output.append("<tr>");
    output.append("<td>Exp:</td>");
    for (double v : vals[1]) {
      output.append("<td align=\"center\">");
      output.append(ext.prettyP(v));
      output.append("</td>");
    }
    output.append("</tr>");
    output.append("</table></html>");

    return output.toString();
  }

}
