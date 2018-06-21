/**
 * 
 */
package org.genvisis.pca.ancestry;

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

  private AncestryPCA() {

  }

  /**
   * @param loader {@link MatrixDataLoading} that provides genotypes in double format (0,1,2,NaN)
   * @param numComponents number of components that will be stored
   * @param log
   * @return {@link SVD} holding results
   */
  public static SVD generatePCs(MatrixDataLoading loader, int numComponents, Logger log) {
    NamedRealMatrix m = loader.getData();

    normalizeGenotypeData(m, log);
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

    normalizeGenotypeData(m, log);
    return svd.getExtraploatedPCs(m, log);
  }

  /**
   * Prepares genotype data for PCA (see https://www.nature.com/articles/ng1847). Can also be used
   * to prepare data for extrapolating PCs
   * 
   * @param m {@link NamedRealMatrix} with genotypes to prepare
   * @param log
   * @return
   */
  static void normalizeGenotypeData(NamedRealMatrix m, Logger log) {

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
      double possible = (double) (2 + 2 * stats.count());
      double pi = (1 + stats.sum()) / possible;
      double norm = Math.sqrt(pi * (1 - pi));
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
   *          {@link AncestryPCA#normalizeGenotypeData(NamedRealMatrix, Logger)}
   * @param log
   * @return {@link SVD}
   */
  private static SVD computePCA(NamedRealMatrix m, int numComponents, Logger log) {

    return SVD.computeSVD(m, numComponents, log);
  }
}
