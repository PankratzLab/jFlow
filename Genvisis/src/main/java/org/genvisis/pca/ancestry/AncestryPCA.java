/**
 * 
 */
package org.genvisis.pca.ancestry;

import java.util.HashMap;
import java.util.Map;
import org.genvisis.common.Logger;
import org.genvisis.common.matrix.MatrixDataLoading;
import org.genvisis.common.matrix.NamedRealMatrix;
import org.genvisis.common.matrix.SVD;
import com.google.common.math.Stats;
import com.google.common.math.StatsAccumulator;

/**
 * Runs PCA that focuses on ancestry PCS (prepares genotype data for ancestry PCA)
 */
public class AncestryPCA {

  private final SVD svd;

  private final Map<String, Stats> statMap;

  private AncestryPCA(SVD svd, Map<String, Stats> statMap) {
    this.svd = svd;
    this.statMap = statMap;
  }

  /**
   * @param loader {@link MatrixDataLoading} that provides genotypes in double format (0,1,2,NaN)
   * @param numComponents number of components that will be stored
   * @param log
   * @return {@link SVD} holding results
   */
  public static SVD generatePCs(MatrixDataLoading loader, int numComponents, Logger log) {
    NamedRealMatrix m = loader.getData();

    Map<String, Stats> statMap = generateStats(m, log);
    normalizeGenotypeData(m, statMap, log);
    return computePCA(m, numComponents, log);
  }

  /**
   * @param svd {@link SVD} that will be used to extrapolate the data
   * @param loader {@link MatrixDataLoading} that will load the matrix
   * @param log
   * @return {@link NamedRealMatrix} holding extrapolated PCs
   */
  public static NamedRealMatrix extrapolatePCs(SVD svd, MatrixDataLoading loader, Logger log) {
    NamedRealMatrix m = loader.getData();

    Map<String, Stats> statMap = generateStats(m, log);
    normalizeGenotypeData(m, statMap, log);
    return svd.getExtraploatedPCs(m, log);
  }

  /**
   * @param m compute {@link Stats} for each marker (row)
   * @param log
   * @return Map from row name to {@link Stats}
   */
  private static Map<String, Stats> generateStats(NamedRealMatrix m, Logger log) {
    Map<String, Stats> statMap = new HashMap<>();
    log.reportTimeInfo("Preparing genotype data");
    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {
      StatsAccumulator statsAccumulator = new StatsAccumulator();
      for (int column = 0; column < m.getDenseMatrix().numCols; column++) {
        double val = m.getDenseMatrix().get(row, column);
        if (valid(val)) {
          if (isNonMissing(val)) {
            statsAccumulator.add(val);
          }
        } else {
          throw new IllegalArgumentException("Invalid  value at marker "
                                             + m.getIndexRowMap().get(row) + " and sample "
                                             + m.getIndexColumnMap().get(column));
        }
      }
      Stats stats = statsAccumulator.snapshot();
      statMap.put(m.getIndexRowMap().get(row), stats);
    }
    return statMap;
  }

  /**
   * Prepares genotype data for PCA (see https://www.nature.com/articles/ng1847). Can also be used
   * to prepare data for extrapolating PCs
   * 
   * @param m {@link NamedRealMatrix} with genotypes to prepare
   * @param statMap a map from row names to the marker {@link Stats} required for normalization
   * @param log
   * @return
   */
  private static Map<String, Stats> normalizeGenotypeData(NamedRealMatrix m,
                                                          Map<String, Stats> statMap, Logger log) {
    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {

      Stats stats = statMap.get(m.getIndexRowMap().get(row));

      double norm = getNormalizationFactor(stats);
      for (int column = 0; column < m.getDenseMatrix().numCols; column++) {
        double val = m.getDenseMatrix().get(row, column);
        if (isNonMissing(val)) {
          double mc = val - stats.mean();
          mc /= norm;
          m.getDenseMatrix().set(row, column, mc);
        } else {
          m.getDenseMatrix().set(row, column, 0);
        }
      }
    }
    log.reportTimeInfo("Finished preparing genotype data");
    return statMap;

  }

  private static double getNormalizationFactor(Stats stats) {
    double possible = (double) (2 + 2 * stats.count());
    double posteriorEstimateOfAlleleFrequency = (1 + stats.sum()) / possible;
    return Math.sqrt(posteriorEstimateOfAlleleFrequency * (1 - posteriorEstimateOfAlleleFrequency));
  }

  private static boolean isNonMissing(double val) {
    return Double.isFinite(val) && val >= 0;
  }

  private static boolean valid(double val) {
    // since these are genotypes, we are using strict equality
    return Double.isNaN(val) || val == 2 || val == 1 || val == 0;
  }

  /**
   * @param m {@link NamedRealMatrix} that has been prepared with
   *          {@link AncestryPCA#generateStats(NamedRealMatrix, Logger)}
   * @param log
   * @return {@link SVD}
   */
  private static SVD computePCA(NamedRealMatrix m, int numComponents, Logger log) {

    return SVD.computeSVD(m, numComponents, log);
  }
}
