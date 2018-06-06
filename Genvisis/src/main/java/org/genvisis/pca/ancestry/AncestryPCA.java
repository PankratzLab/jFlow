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
   * @param log
   * @return {@link SVD} holding results
   */
  public static SVD generatePCs(MatrixDataLoading loader, Logger log) {
    NamedRealMatrix m = loader.getData();

    normalizeGenotypeData(m, log);
    return computePCA(m, log);
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
    //    Let gij be a matrix of genotypes for SNP i and individual j, where i = 1 to M 
    //        and j = 1 to N. We subtract the row mean μi = (Σjgij)/N from each entry in row i 
    //        to obtain a matrix with row sums equal to 0; missing entries are excluded from
    //        the computation of μi and are subsequently set to 0. We then normalize row i by 
    //        dividing each entry by √(pi(1 − pi)), where pi is a posterior estimate of the unobserved
    //        underlying allele frequency of SNP i defined by pi = (1 + Σjgij)/(2 + 2N), 
    //        with missing entries excluded from the computation.

    for (int row = 0; row < m.getM().getRowDimension(); row++) {
      StatsAccumulator statsAccumulator = new StatsAccumulator();
      for (int column = 0; column < m.getM().getColumnDimension(); column++) {
        double val = m.getM().getEntry(row, column);
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
      double possible = (double) (2 + 2 * m.getM().getColumnDimension());
      double pi = (1 + stats.sum()) / possible;
      double norm = Math.sqrt(pi * (1 - pi));
      for (int column = 0; column < m.getM().getColumnDimension(); column++) {
        double val = m.getM().getEntry(row, column);
        if (isNonMissing(val)) {
          double mc = val - stats.mean();
          mc /= norm;
          m.getM().setEntry(row, column, mc);
        } else {
          m.getM().setEntry(row, column, 0);
        }
      }
    }
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
  private static SVD computePCA(NamedRealMatrix m, Logger log) {

    return SVD.computeSVD(m, log);
  }
}
