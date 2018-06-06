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
 * Runs PCA starting from a {@link NamedRealMatrix}
 */
public class AncestryPCA {

  private AncestryPCA() {

  }

  public static SVD generatePCs(MatrixDataLoading loader, Logger log) {
    NamedRealMatrix m = loader.getData();

    prepareForPCA(m, log);
    return computePCA(m, log);
  }

  /**
   * Prepares genotype data for PCA (see https://www.nature.com/articles/ng1847)
   * 
   * @param m
   * @param log
   * @return
   */
  public static void prepareForPCA(NamedRealMatrix m, Logger log) {
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
          if (Double.isFinite(val)) {
            statsAccumulator.add(val);
          }
        } else {
          throw new IllegalArgumentException("Invalid  value at marker "
                                             + m.getIndexRowMap().get(row) + " and sample "
                                             + m.getIndexColumnMap().get(column));
        }
      }
      Stats stats = statsAccumulator.snapshot();
      for (int column = 0; column < m.getM().getColumnDimension(); column++) {
        double val = m.getM().getEntry(row, column);
        if (Double.isFinite(val)) {
          m.getM().setEntry(row, column, 0);
        } else {
          m.getM().setEntry(row, column, val - stats.mean());
        }
      }
    }
  }

  private static boolean valid(double val) {
    // since these are genotypes, we are using strict equality
    return Double.isNaN(val) || val == 2 || val == 1 || val == 0;
  }

  private static SVD computePCA(NamedRealMatrix m, Logger log) {

    return SVD.computeSVD(m, log);
  }
}
