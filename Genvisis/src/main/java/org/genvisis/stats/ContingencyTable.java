// these bad boys assume that the rows represent the explanatory variable levels, and the columns
// response levels
// ergo data[explanatory][response]

package org.genvisis.stats;

import java.io.IOException;

import org.genvisis.common.Array;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class ContingencyTable {
  public static final double DEFAULT_ALPHA = 0.05;
  public static final int DEFAULT_SIGFIGS = 3;

  public static double stderr(double[][] data) {
    double p1, p2, n1, n2;

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    p1 = data[0][0] / (data[0][0] + data[0][1]);
    p2 = data[1][0] / (data[1][0] + data[1][1]);
    n1 = data[0][0] + data[0][1];
    n2 = data[1][0] + data[1][1];

    return Math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
  }

  public static double diffProportions(double[][] data) {
    double p1, p2;

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    p1 = data[0][0] / (data[0][0] + data[0][1]);
    p2 = data[1][0] / (data[1][0] + data[1][1]);

    return p1 - p2;
  }

  public static double[] diffCI(double[][] data) {
    return diffCI(data, DEFAULT_ALPHA);
  }

  public static double[] diffCI(double[][] data, double alpha) {
    double diff, z, se;

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    diff = diffProportions(data);
    z = ProbDist.NormDistReverse(alpha);
    se = stderr(data);

    return new double[] {diff - z * se, diff + z * se};
  }

  public static double relativeRisk(double[][] data) {
    double p1, p2;

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    p1 = data[0][0] / (data[0][0] + data[0][1]);
    p2 = data[1][0] / (data[1][0] + data[1][1]);

    return p1 / p2;
  }

  public static double[] relativeRiskCI(double[][] data) {
    return relativeRiskCI(data, DEFAULT_ALPHA);
  }

  public static double[] relativeRiskCI(double[][] data, double alpha) { // for
    // large
    // samples
    // only
    double logRR, p1, p2, n1, n2, se, z;

    p1 = data[0][0] / (data[0][0] + data[0][1]);
    p2 = data[1][0] / (data[1][0] + data[1][1]);
    n1 = Array.sum(data[0]);
    n2 = Array.sum(data[1]);

    logRR = Math.log(relativeRisk(data));
    z = ProbDist.NormDistReverse(alpha);
    se = Math.sqrt((1 - p1) / (n1 * p1) + (1 - p2) / (n2 * p2));

    return new double[] {Math.exp(logRR - z * se), Math.exp(logRR + z * se)};
  }

  public static double oddsRatio(double[][] data) {
    double p1, p2, odds1, odds2;

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    p1 = data[0][0] / (data[0][0] + data[0][1]);
    p2 = data[1][0] / (data[1][0] + data[1][1]);
    odds1 = p1 / (1 - p1);
    odds2 = p2 / (1 - p2);

    return odds1 / odds2;
  }

  public static double oddsRatioSE(double[][] data) {
    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    return Math.sqrt(1 / data[0][0] + 1 / data[0][1] + 1 / data[1][0] + 1 / data[1][1]);
  }

  public static double[] oddsRatioCI(double[][] data) {
    return oddsRatioCI(data, DEFAULT_ALPHA);
  }

  public static double[] oddsRatioCI(double[][] data, double alpha) {
    double logOR, z, se;

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - assuming a 2x2 table");
      System.exit(1);
    }

    logOR = Math.log(oddsRatio(data));
    z = ProbDist.NormDistReverse(alpha);
    se = oddsRatioSE(data);

    return new double[] {Math.exp(logOR - z * se), Math.exp(logOR + z * se)};
  }

  public static double[] computeRowSums(double[][] data) {
    double[] rowSums = new double[data.length];

    for (int i = 0; i < data.length; i++) {
      if (data[i].length != data[0].length) {
        System.err.println("Error - can't do a ChiSquare on a matrix with different numbers of columns");
      }
      rowSums[i] = Array.sum(data[i]);
    }

    return rowSums;
  }

  public static double[] computeColSums(double[][] data) {
    double[] colSums = new double[data[0].length];

    for (int j = 0; j < data[0].length; j++) {
      colSums[j] = Array.sum(Matrix.extractColumn(data, j));
    }

    return colSums;
  }

  public static double[][] computeExpecteds(double[][] data) {
    double[] rowSums, colSums;
    double[][] expecteds;
    double total;

    rowSums = computeRowSums(data);
    colSums = computeColSums(data);
    total = Array.sum(rowSums);

    expecteds = new double[data.length][data[0].length];
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        expecteds[i][j] = rowSums[i] * colSums[j] / total;
      }
    }

    return expecteds;
  }

  public static double ChiSquare(int[][] counts) {
    return ChiSquare(Matrix.toDoubleArrays(counts), true, true);
  }

  public static double ChiSquare(int[][] counts, boolean verbose) {
    return ChiSquare(Matrix.toDoubleArrays(counts), verbose, true);
  }

  public static double ChiSquare(double[][] data) {
    return ChiSquare(data, true, true);
  }

  public static double ChiSquare(double[][] data, boolean verbose,
                                 boolean possibleYatesCorrection) {
    double[][] expecteds;
    double chi;
    int numCellsLT5, numCells;
    double minExpected = 5;
    double diff;
    boolean useYates;

    expecteds = computeExpecteds(data);
    numCellsLT5 = 0;
    useYates = false;
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        if (expecteds[i][j] < 5) {
          numCellsLT5++;
          if (expecteds[i][j] < minExpected) {
            minExpected = expecteds[i][j];
          }
        }
      }
    }

    numCells = data.length * data[0].length;
    if ((double) numCellsLT5 / (double) numCells >= 0.20) {
      if (verbose) {
        System.out.println("FYI, " + numCellsLT5 + " cells ("
                           + ext.formDeci((double) numCellsLT5 / (double) numCells, 1, true)
                           + "%) have expected count less than 5. The minimum expected count is "
                           + ext.formDeci(minExpected, 2, true)
                           + ". Might consider using Fisher's exact test"
                           + (possibleYatesCorrection ? "; in the meantime using a Yates continuity correction"
                                                      : ""));
      }
      useYates = true;
    }

    chi = 0;
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        diff = Math.abs((data[i][j] - expecteds[i][j]))
               - (useYates && possibleYatesCorrection ? 0.5 : 0);
        chi += diff * diff / expecteds[i][j];
      }
    }

    return chi;
  }

  public static double likelihoodRatioStatistic(int[][] counts) {
    return likelihoodRatioStatistic(Matrix.toDoubleArrays(counts));
  }

  public static double likelihoodRatioStatistic(double[][] data) {
    double[][] expecteds;
    double llr;

    llr = 0;
    expecteds = computeExpecteds(data);
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        llr += data[i][j] * Math.log(data[i][j] / expecteds[i][j]);
      }
    }
    llr *= 2;

    return llr;
  }

  public static double[][] residuals(int[][] counts) {
    return residuals(Matrix.toDoubleArrays(counts));
  }

  public static double[][] residuals(double[][] data) {
    double[] rowSums, colSums;
    double[][] expecteds, residuals;
    double total;

    rowSums = computeRowSums(data);
    colSums = computeColSums(data);
    total = Array.sum(rowSums);

    residuals = new double[data.length][data[0].length];
    expecteds = computeExpecteds(data);
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        residuals[i][j] =
            (data[i][j] - expecteds[i][j]) / Math.sqrt(expecteds[i][j] * (1 - rowSums[i] / total)
                                                       * (1 - colSums[j] / total));
      }
    }

    return residuals;
  }

  public static double linearCorrelation(double[][] data) {
    return linearCorrelation(data, Array.toDoubleArray(Array.arrayOfIndices(data.length)),
                             Array.toDoubleArray(Array.arrayOfIndices(data[0].length)));
  }

  public static double linearCorrelation(double[][] data, double[] rowScores, double[] colScores) {
    double[] rowSums, colSums;
    double total, uBar, vBar, uSD, vSD;
    double covar;

    rowSums = computeRowSums(data);
    colSums = computeColSums(data);
    total = Array.sum(rowSums);

    if (rowScores.length != rowSums.length) {
      System.err.println("Error - number of rowScores does not jive with the number of rows");
    }
    if (colScores.length != colSums.length) {
      System.err.println("Error - number of colScores does not jive with the number of columns");
    }

    uBar = 0;
    for (int i = 0; i < rowSums.length; i++) {
      uBar += rowSums[i] / total * rowScores[i];
    }

    vBar = 0;
    for (int j = 0; j < colSums.length; j++) {
      vBar += colSums[j] / total * colScores[j];
    }

    covar = 0;
    for (int i = 0; i < rowSums.length; i++) {
      for (int j = 0; j < colSums.length; j++) {
        covar += (rowScores[i] - uBar) * (colScores[j] - vBar) * data[i][j] / total;
      }
    }

    uSD = 0;
    for (int i = 0; i < rowSums.length; i++) {
      uSD += Math.pow((rowScores[i] - uBar), 2) * rowSums[i] / total;
    }

    vSD = 0;
    for (int j = 0; j < colSums.length; j++) {
      vSD += Math.pow((colScores[j] - vBar), 2) * colSums[j] / total;
    }

    return covar / Math.sqrt(uSD * vSD);
  }

  public static double linearTrendStatistic(double[][] data) {
    return linearTrendStatistic(data, Array.toDoubleArray(Array.arrayOfIndices(data.length)),
                                Array.toDoubleArray(Array.arrayOfIndices(data[0].length)));
  }

  public static double linearTrendStatistic(double[][] data, double[] rowScores,
                                            double[] colScores) {
    double r = linearCorrelation(data, rowScores, colScores);
    double n = Matrix.sum(data);

    return (n - 1) * Math.pow(r, 2);
  }

  public static double[] midRanks(double[] array) {
    double[] midRanks = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      midRanks[i] = 0;
      for (int j = 0; j < i; j++) {
        midRanks[i] += array[j];
      }
      midRanks[i] += (array[i] + 1) / 2;
    }

    return midRanks;
  }

  public static double binomialCoeff(int a, int b) {
    double d;

    d = BinomialDistribution.fact(a)
        / (BinomialDistribution.fact(b) * BinomialDistribution.fact(a - b));

    return d;
  }

  public static double probabilityBC(int n11, int n, int n1_, int n2_, int n_1) {
    double d;

    d = binomialCoeff(n1_, n11) * binomialCoeff(n2_, n_1 - n11) / binomialCoeff(n, n_1);

    return d;
  }

  public static double FishersExact(int[][] data) {
    return FishersExact(data, false, false);
  }

  public static double FishersExact(int[][] data, boolean twosided, boolean midvalue) {
    double d, d2, exact, temp;
    int n11, n, n1_, n2_, n_1, max;

    n11 = data[0][0];
    n1_ = data[0][0] + data[0][1];
    n2_ = data[1][0] + data[1][1];
    n_1 = data[0][0] + data[1][0];
    n = n1_ + n2_;
    max = n1_ < n_1 ? n1_ : n_1;

    d = 0;
    for (int i = n11; i <= max; i++) {
      d += probabilityBC(i, n, n1_, n2_, n_1);
      if (i == n11 && midvalue) {
        d /= 2;
      }
    }

    d2 = 0;
    if (twosided) {
      exact = probabilityBC(n11, n, n1_, n2_, n_1);
      if (midvalue) {
        d = exact / 2;
        for (int i = 0; i <= max; i++) {
          temp = probabilityBC(i, n, n1_, n2_, n_1);
          if (temp < exact) {
            d2 += probabilityBC(i, n, n1_, n2_, n_1) * 2;
          }
        }
      } else {
        for (int i = 0; i <= max - n11; i++) {
          temp = probabilityBC(i, n, n1_, n2_, n_1);
          if (temp <= exact) {
            d2 += probabilityBC(i, n, n1_, n2_, n_1);
          }
        }
      }
    }

    return d + d2;
  }

  public static double phiCoefficient(int[][] data) {
    double chisq = ChiSquare(data);
    double n = Matrix.sum(data);

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - a phi coefficient can only be calculated for a 2x2 table");
    }

    return Math.sqrt(chisq / (n));
  }

  public static double tetrachoricCorrelation(int[][] data) {
    double chisq = ChiSquare(data);
    double n = Matrix.sum(data);

    if (data.length != 2 || data[0].length != 2) {
      System.err.println("Error - a tetrachoric correlation can only be calculated for a 2x2 table");
    }

    return Math.sqrt(chisq / (chisq + n));
  }

  public static double tetrachoricCorrelationAdjusted(int[][] data) {
    return ContingencyTable.tetrachoricCorrelation(data) / Math.sqrt(0.5);
  }

  public static void demo() {
    int[][] iData;
    double[][] data;
    double[] ci;
    double[][] residuals;

    // data = new int[][] {{ 189, 10845 },
    // { 104, 10933 }};
    // data = new int[][] {{ 172, 90 },
    // { 173, 346 }};
    // data = new int[][] {{ 509, 116 },
    // { 398, 104 }};
    // data = new double[][] {{ 19.5, 132.5 }, // Question 2.33
    // { 11.5, 52.5 }};
    // data = new double[][] {{ 0.5, 9.5 }, // Question 2.33
    // { 6.5, 97.5 }};
    // data = new double[][] {{ 19.5, 141.5 }, // Question 2.33b
    // { 17.5, 149.5 }};
    // data = new double[][] {{ 60.57435793, 502.4256421 }, // GBA 2+2
    // variants, unbiased
    // { 21, 346 }};
    // data = new double[][] {{ 20.5, 542.5 }, // GBA 2 mutations, unbiased
    // { 4, 363 }};
    // data = new double[][] {{ 85.167, 884.834 }, // Rep1-263 allele count,
    // no Parkin
    // { 43, 679 }};
    // data = new double[][] {{ 79, 701 }, // Rep1-263 allele count, no
    // Parkin
    // { 77, 989 }};
    // data = new double[][] {{ 87.5, 882.5 }, // Rep1-263 allele count, no
    // Parkin
    // { 43, 685 }};
    // data = new double[][] {{ 6, 31}, // Marder's GBA - Early versus Late
    // { 2, 59}};
    // data = new double[][] { {200.9, 294.0833333}, // Our GBA - Early versus Late
    // {38.51666667, 29.5}};
    // data = new double[][] {{ 39, 15}, // Random data test
    // { 25, 21}};
    // data = new double[][] {{ 3, 1}, // Random data test
    // { 1, 3}};
    data = new double[][] {{3, 500}, // Max spread to get significance for CNVs
                           {1, 6000}};

    System.out.println("For the following table: ");
    System.out.println("                  Response");
    System.out.println("              Yes          No");
    System.out.println("Group 1" + ext.formStr(data[0][0] + "", 10) + "\t"
                       + ext.formStr(data[0][1] + "", 10));
    System.out.println("Group 2" + ext.formStr(data[1][0] + "", 10) + "\t"
                       + ext.formStr(data[1][1] + "", 10));

    System.out.println("p1 = " + ext.formDeci(data[0][0] / (data[0][0] + data[0][1]), 4, true));
    System.out.println("p2 = " + ext.formDeci(data[1][0] / (data[1][0] + data[1][1]), 4, true));
    System.out.println();
    ci = diffCI(data);
    System.out.println("The difference of proportions is "
                       + ext.formDeci(diffProportions(data), 4, true) + " ("
                       + ext.formDeci(ci[0], 4, true) + ", " + ext.formDeci(ci[1], 4, true) + ")");
    System.out.println();
    ci = relativeRiskCI(data);
    System.out.println("The relative risk is " + ext.formDeci(relativeRisk(data), 4, true) + " ("
                       + ext.formDeci(ci[0], 4, true) + ", " + ext.formDeci(ci[1], 4, true) + ")");
    System.out.println();
    ci = oddsRatioCI(data);
    System.out.println("The odds ratio is " + ext.formDeci(oddsRatio(data), 4, true) + " ("
                       + ext.formDeci(ci[0], 4, true) + ", " + ext.formDeci(ci[1], 4, true) + ")");
    System.out.println();

    // data = new int[][] {{ 762, 327, 468},
    // { 484, 239, 477}};
    // data = new int[][] {{ 509, 116 },
    // { 398, 104 }};
    // data = new int[][] {{ 871, 444, 873}, // Question 2.19
    // { 302, 80, 43}};
    // data = new int[][] {{ 871, 444}, // Question 2.19c
    // { 302, 80}};
    // data = new int[][] {{ 1315, 873}, // Question 2.19c
    // { 382, 43}};
    // data = new double[][] {{ 9, 44, 13, 10}, // Question 2.27
    // { 11, 52, 23, 22},
    // { 9, 41, 12, 27}};
    // data = new double[][] {{ 266.417, 618.417, 85.167},
    // { 207.000, 472.000, 43.000}}; // chi for Rep1 , no Parkin
    // data = new double[][] {{ 79, 185, 516 }, // Rep1-263 allele count, no
    // Parkin
    // { 77, 243, 746 }};

    System.out.println("The Pearson (score) chi-square statistic is "
                       + ext.formDeci(ChiSquare(data), 4, true) + " (p="
                       + ProbDist.ChiDist(ChiSquare(data), (data.length - 1) * (data[0].length - 1))
                       + ")");
    System.out.println("The Likelihood-Ratio chi-square statistic is "
                       + ext.formDeci(likelihoodRatioStatistic(data), 4, true) + " (p="
                       + ProbDist.ChiDist(likelihoodRatioStatistic(data),
                                          (data.length - 1) * (data[0].length - 1))
                       + ")");
    System.out.println();
    System.out.println("The residuals are:");
    residuals = residuals(data);
    for (double[] residual : residuals) {
      for (int j = 0; j < residuals[0].length; j++) {
        System.out.print((j == 0 ? "" : "\t") + ext.formDeci(residual[j], 3));
      }
      System.out.println();
    }
    System.out.println();

    // data = new int[][] {{ 17066, 14464, 788, 126, 37},
    // { 48, 38, 5, 1, 1}};
    data = new double[][] {{266.417, 618.417, 85.167}, {207.000, 472.000, 43.000}}; // chi
                                                                                    // for
                                                                                    // Rep1
                                                                                    // , no
                                                                                    // Parkin
    // rowScores = new double[] {0, 1};
    // colScores = new double[] {0, 0.5, 1.5, 4, 7};
    //
    // System.out.println("Using provided row and column scores:");
    // System.out.println("The linear correlation is
    // "+ext.formDeci(linearCorrelation(Array.doubleArrays(data), rowScores,
    // colScores), 5));
    // System.out.println("The linear trend test statistic is
    // "+ext.formDeci(linearTrendStatistic(Array.doubleArrays(data),
    // rowScores, colScores), 3) +"
    // (p="+ext.prettyP(ProbDist.ChiDist(linearTrendStatistic(Array.doubleArrays(data),
    // rowScores, colScores), 1))+")");
    System.out.println("Using equally spaced ranks:");
    System.out.println("The linear correlation is " + ext.formDeci(linearCorrelation(data), 5));
    System.out.println("The linear trend test statistic is "
                       + ext.formDeci(linearTrendStatistic(data), 3) + " (p="
                       + ext.prettyP(ProbDist.ChiDist(linearTrendStatistic(data), 1)) + ")");
    System.out.println("Using midranks:");
    System.out.println(Array.toStr(midRanks(computeColSums(data))));
    System.out.println("The linear correlation is "
                       + ext.formDeci(linearCorrelation(data, midRanks(computeRowSums(data)),
                                                        midRanks(computeColSums(data))),
                                      5));
    System.out.println("The linear trend test statistic is "
                       + ext.formDeci(linearTrendStatistic(data, midRanks(computeRowSums(data)),
                                                           midRanks(computeColSums(data))),
                                      3)
                       + " (p="
                       + ext.prettyP(ProbDist.ChiDist(linearTrendStatistic(data,
                                                                           midRanks(computeRowSums(data)),
                                                                           midRanks(computeColSums(data))),
                                                      1))
                       + ")");
    System.out.println();

    // data = new int[][] {{ 3, 1},
    // { 1, 3}};
    // data = new int[][] {{ 7, 8}, // Question 2.29
    // { 0, 15}};
    // iData = new int[][] {{ 21, 2}, // Question 2.31
    // { 15, 3}};
    // iData = new int[][] {{ 6, 4}, // Question 3.10a
    // { 2, 8}};
    // iData = new int[][] {{ 4, 3}, // Question 3.10b
    // { 1, 5}};
    // iData = new int[][] {{ 5, 3}, // Question 3.10c
    // { 3, 6}};
    // iData = new int[][] {{ 15, 12}, // Question 3.10all
    // { 6, 19}};
    // iData = new int[][] {{ 6, 31}, // Marder's GBA - Early versus Late
    // { 2, 59}};
    // iData = new int[][] {{ 39, 15}, // Random data test
    // { 25, 21}};
    iData = new int[][] {{8, 456}, // Our L444P variant data (GBA)
                         {0, 344}};

    System.out.println("Exact, one-sided p=" + FishersExact(iData));
    System.out.println("Exact, two-sided p=" + FishersExact(iData, true, false));
    System.out.println("Midvalue, one-sided p=" + FishersExact(iData, false, true));
    System.out.println("Midvalue, two-sided p=" + FishersExact(iData, true, true));
    System.out.println("Odds-ratio: " + oddsRatio(Matrix.toDoubleArrays(iData)));

  }

  public static double ChiSquareOptimizedSignificance(double[][] data, boolean twosided,
                                                      boolean verbose) {
    double[][] expecteds;
    int numCellsLT5, numCells;
    double minExpected = 5;
    boolean fishy = false;

    numCellsLT5 = 0;
    expecteds = computeExpecteds(data);
    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        if (expecteds[i][j] < 5) {
          numCellsLT5++;
          if (expecteds[i][j] < minExpected) {
            minExpected = expecteds[i][j];
          }
        }
      }
    }

    numCells = data.length * data[0].length;
    if ((double) numCellsLT5 / (double) numCells >= 0.20) {
      System.out.print(numCellsLT5 + " cell" + (numCellsLT5 > 1 ? "s have" : "  has ")
                       + " an expected count less than 5 (min=" + ext.formDeci(minExpected, 2, true)
                       + "); ");
      if (numCells == 4) {
        System.out.println("using Fisher's exact");
        fishy = true;
      } else {
        System.out.println("however, can't use Fisher's exact because the table is larger than 2x2");
      }
    }
    if (fishy) {
      return FishersExact(Matrix.toIntArraysRounded(data), twosided, false);
    } else {
      return ProbDist.ChiDist(ChiSquare(data), (data.length - 1) * (data[0].length - 1));
    }
  }

  public static void main(String[] args) throws IOException {
    demo();
  }
}
