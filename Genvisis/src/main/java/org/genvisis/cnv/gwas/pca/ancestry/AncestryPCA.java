/**
 * 
 */
package org.genvisis.cnv.gwas.pca.ancestry;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import org.genvisis.cnv.ejml.matrix.MatrixDataLoading;
import org.genvisis.cnv.ejml.matrix.NamedRealMatrix;
import org.genvisis.cnv.ejml.matrix.SVD;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import com.google.common.math.Stats;
import com.google.common.math.StatsAccumulator;

/**
 * Runs PCA that focuses on ancestry PCS (prepares genotype data for ancestry PCA)
 */
public class AncestryPCA implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  private final SVD svd;

  private final Map<String, Stats> statMap;

  private AncestryPCA(SVD svd, Map<String, Stats> statMap) {
    this.svd = svd;
    this.statMap = statMap;
  }

  /**
   * @return the underlying {@link SVD}
   */
  public SVD getSvd() {
    return svd;
  }

  /**
   * @param loader {@link MatrixDataLoading} that provides genotypes in double format (0,1,2,NaN)
   * @param numComponents number of components that will be stored
   * @param log
   * @return {@link AncestryPCA} holding the {@link SVD} and required {@link Stats} for
   *         normalization
   */
  public static AncestryPCA generatePCs(MatrixDataLoading loader, int numComponents, Logger log) {
    NamedRealMatrix m = loader.getData();

    Map<String, Stats> statMap = generateStats(m, log);
    normalizeGenotypeData(m, statMap, log);
    SVD svd = computePCA(m, numComponents, log);

    return new AncestryPCA(svd, statMap);
  }

  /**
   * @param svd {@link AncestryPCA} that will be used to extrapolate the data
   * @param loader {@link MatrixDataLoading} that will load the matrix
   * @param log
   * @return {@link NamedRealMatrix} holding extrapolated PCs
   */
  public static NamedRealMatrix extrapolatePCs(AncestryPCA ancestryPCA, MatrixDataLoading loader,
                                               Logger log) {
    NamedRealMatrix m = loader.getData();

    normalizeGenotypeData(m, ancestryPCA.statMap, log);
    return ancestryPCA.svd.getExtraploatedPCs(m, log);
  }

  public void writeSerial(String serFile) {
    SerializedFiles.writeSerial(this, serFile, true);
  }

  /**
   * @param serFile
   * @return
   */
  public static AncestryPCA readSerial(String serFile) {
    return (AncestryPCA) SerializedFiles.readSerial(serFile);
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
        checkValid(val, row, column);
        if (isNonMissing(val)) {
          statsAccumulator.add(val);
        }
      }
      Stats stats = statsAccumulator.snapshot();
      statMap.put(m.getNameForRowIndex(row), stats);
    }
    return statMap;
  }

  /**
   * Prepares genotype data for PCA (see https://www.nature.com/articles/ng1847). Can also be used
   * to prepare data for extrapolating PCs Note that any data not found in the statMap will be
   * ignored(and not normalized).
   * 
   * @param m {@link NamedRealMatrix} with genotypes to prepare
   * @param statMap a map from row names to the marker {@link Stats} required for normalization
   * @param log
   * @return
   */
  private static Map<String, Stats> normalizeGenotypeData(NamedRealMatrix m,
                                                          Map<String, Stats> statMap, Logger log) {

    if (!statMap.keySet().containsAll(m.getRowMap().keySet())) {
      log.reportTimeWarning(m.getRowMap().keySet().size() - statMap.keySet().size()
                            + " markers in the input data were not found in the data "
                            + "used to generate normalization stats, these markers will be ignored. "
                            + "This is typically not a problem");
    }
    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {
      String rowKey = m.getNameForRowIndex(row);
      if (statMap.containsKey(rowKey)) {
        Stats stats = statMap.get(rowKey);

        double norm = getNormalizationFactor(stats);
        for (int column = 0; column < m.getDenseMatrix().numCols; column++) {
          double val = m.getDenseMatrix().get(row, column);
          checkValid(val, row, column);
          if (isNonMissing(val)) {
            double mc = val - stats.mean();
            mc /= norm;
            m.getDenseMatrix().set(row, column, mc);
          } else {
            m.getDenseMatrix().set(row, column, 0);
          }

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

  private static void checkValid(double val, int row, int column) {
    // since these are genotypes, we are using strict equality
    if (!(Double.isNaN(val) || val == 2 || val == 1 || val == 0)) {
      throw new IllegalArgumentException("Invalid  value at row " + row + " and column " + column);
    }

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
