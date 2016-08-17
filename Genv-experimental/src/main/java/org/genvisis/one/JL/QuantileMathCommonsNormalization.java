package org.genvisis.one.JL;

import java.util.ArrayList;

import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.genvisis.common.Array;

import com.google.common.primitives.Doubles;

public class QuantileMathCommonsNormalization {
  private static final RankingAlgorithm COV_RANKER_TIE =
      new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);

  /**
   * Quantile normalize data where missing values are allowed. the missing values are replaced with
   * the median or mean of the other probes in the sample If chosen keep median value instead of the
   * NA value.
   *
   * NB: For now NA is defined as "-999"!!
   *
   * @param rawData matrix containing expression/methylation data
   * @param useMean use mean for guessing the NA value if false use median
   * @param retainNA retain the NA values, put NA values back after normalization
   */
  public static void QuantileNormAdressingNaValuesBeforeQN(double[][] rawData, boolean useMedian,
                                                           boolean retainNA) {
    boolean[][] wasNA = new boolean[rawData.length][rawData[1].length];

    for (int s = 0; s < rawData[1].length; ++s) {

      ArrayList<Double> nonNAvalues = new ArrayList<Double>();

      boolean needsReplacement = false;

      for (int p = 0; p < rawData.length; ++p) {
        if (Double.isNaN(rawData[p][s])) {
          needsReplacement = true;
          wasNA[p][s] = true;
        } else {
          wasNA[p][s] = false;
          nonNAvalues.add(rawData[p][s]);
        }
      }

      if (needsReplacement) {
        double replacementValue;
        if (useMedian) {
          replacementValue = Array.median(Doubles.toArray(nonNAvalues));
        } else {
          replacementValue = Array.mean(Doubles.toArray(nonNAvalues));

        }

        for (int p = 0; p < rawData.length; ++p) {
          if (wasNA[p][s]) {
            rawData[p][s] = replacementValue;
          }
        }
      }
    }

    quantilenormalize(rawData);

    if (retainNA) {
      for (int s = 0; s < rawData[1].length; ++s) {
        for (int p = 0; p < rawData.length; ++p) {
          if (wasNA[p][s]) {
            rawData[p][s] = Double.NaN;
          }
        }
      }
    }

  }

  /**
   * Quantile normalize a double[][] double[probes][sample]
   *
   * @param rawData matrix containing expression/methylation data
   */
  public static void quantilenormalize(double[][] rawData) {
    System.out.println("\nPerforming quantile normalization:");
    // Calculate the average expression, when per sample all raw expression levels have been
    // ordered:


    int probeCount = rawData.length;
    int sampleCount = rawData[probeCount - 1].length;

    double[] rankedMean = new double[probeCount];
    for (int sampleID = 0; sampleID < sampleCount; sampleID++) {
      double[] x = new double[probeCount];

      for (int probeID = 0; probeID < probeCount; probeID++) {
        x[probeID] = rawData[probeID][sampleID];
      }
      java.util.Arrays.sort(x);
      for (int probeID = 0; probeID < probeCount; probeID++) {
        rankedMean[probeID] += x[probeID];
      }
    }

    for (int probeID = 0; probeID < probeCount; probeID++) {
      rankedMean[probeID] /= sampleCount;
    }

    double[] rankedMeanClasses = new double[probeCount - 1];

    for (int probeID = 0; probeID < (probeCount - 1); probeID++) {
      rankedMeanClasses[probeID] = ((rankedMean[probeID] + rankedMean[probeID + 1]) / 2);
    }

    // Iterate through each sample:
    for (int s = 0; s < sampleCount; s++) {
      double[] probes = new double[probeCount];
      for (int p = 0; p < probeCount; p++) {
        probes[p] = rawData[p][s];
      }
      double[] probesRanked = COV_RANKER_TIE.rank(probes);

      double[] probesQuantileNormalized = new double[probeCount];
      for (int p = 0; p < probeCount; p++) {

        if ((probesRanked[p] % 1) != 0) {
          probesQuantileNormalized[p] = rankedMeanClasses[(int) Math.floor((probesRanked[p] - 1))];
        } else {
          probesQuantileNormalized[p] = rankedMean[(int) (probesRanked[p] - 1)];
        }

        rawData[p][s] = probesQuantileNormalized[p];
      }
      // double[] probesRankedAfterQQNorm = rda.rank(probesQuantileNormalized, false);

    }
  }
}
