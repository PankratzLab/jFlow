/**
 * 
 */
package org.genvisis.common;

import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.data.DenseMatrix64F;
import org.genvisis.stats.Maths;

/**
 * Common methods for operating on {@link RealMatrix}
 */
public class RealMatrixUtils {

  private RealMatrixUtils() {

  }

  /**
   * Attaches column and row names to the matrix
   */
  public static class NamedRealMatrix {

    private final String[] rowNames;
    private final String[] columnNames;
    private final RealMatrix m;

    /**
     * @param rowNames
     * @param columnNames
     * @param m
     */
    public NamedRealMatrix(String[] rowNames, String[] columnNames, RealMatrix m) {
      super();
      if (rowNames == null || rowNames.length != m.getRowDimension()) {
        throw new IllegalArgumentException("Mismatched row and row names length");
      }
      if (columnNames == null || columnNames.length != m.getColumnDimension()) {
        throw new IllegalArgumentException("Mismatched column and column names length");
      }
      this.rowNames = rowNames;
      this.columnNames = columnNames;
      this.m = m;
    }

    /**
     * @return the rowNames
     */
    public String[] getRowNames() {
      return rowNames;
    }

    /**
     * @return the columnNames
     */
    public String[] getColumnNames() {
      return columnNames;
    }

    /**
     * @return the m
     */
    public RealMatrix getM() {
      return m;
    }

  }

  /**
   * @param m Convert this apache {@link RealMatrix} to an EJML style {@link DenseMatrix64F}
   * @return {@link DenseMatrix64F}
   */
  public static DenseMatrix64F toDenseMatrix64F(RealMatrix m) {
    DenseMatrix64F dm = new DenseMatrix64F(1, 1);
    dm.reshape(m.getRowDimension(), m.getColumnDimension());
    for (int row = 0; row < m.getRowDimension(); row++) {
      for (int column = 0; column < m.getColumnDimension(); column++) {
        dm.add(row, column, m.getEntry(row, column));
      }
    }
    return dm;
  }

  /**
   * compute fold-change (by column) , and then center the matrix so each row has median of 0;
   * 
   * @param m an {@link RealMatrix} that has been FC-ed by column and centered by row
   */
  public static void foldChangeAndCenter(RealMatrix m) {
    double[] medians = new double[m.getColumnDimension()]; //In genvisis, samples are typically columns

    //    convert columns to log2 fold-change from median
    for (int column = 0; column < m.getColumnDimension(); column++) {
      double[] tmp = new double[m.getRowDimension()];//In genvisis, markers are typically rows
      for (int row = 0; row < m.getRowDimension(); row++) {
        tmp[row] += m.getEntry(row, column);
      }
      medians[column] = ArrayUtils.median(tmp);
    }
    for (int row = 0; row < m.getRowDimension(); row++) {
      for (int column = 0; column < m.getColumnDimension(); column++) {
        double entry = m.getEntry(row, column);
        if (entry > 0) {
          double standard = Maths.log2(entry / medians[column]);
          m.setEntry(row, column, standard);
        } else {
          m.setEntry(row, column, 0);
        }
      }
    }
    // center rows to median of 0

    centerRowsToMedian(m);
  }

  /**
   * @param m Center the rows of this {@link RealMatrix} to a median of 0
   */
  public static void centerRowsToMedian(RealMatrix m) {
    for (int row = 0; row < m.getRowDimension(); row++) {
      double[] tmp = m.getRow(row);
      double median = ArrayUtils.median(tmp);
      for (int j = 0; j < m.getColumnDimension(); j++) {
        m.setEntry(row, j, tmp[j] - median);
      }
    }
  }

  /**
   * @param m scale the columns of this matrix to a mean of 0 and SD of 1
   */
  public static void scaleAndCenterColumns(RealMatrix m) {
    double[] sds = new double[m.getColumnDimension()];
    double[] mean = new double[m.getColumnDimension()];

    for (int column = 0; column < m.getColumnDimension(); column++) {
      double[] tmp = new double[m.getRowDimension()];
      for (int i = 0; i < m.getRowDimension(); i++) {
        tmp[i] += m.getEntry(i, column);
      }
      mean[column] = ArrayUtils.mean(tmp);
      sds[column] = ArrayUtils.stdev(tmp);
    }
    for (int row = 0; row < m.getRowDimension(); row++) {
      for (int column = 0; column < m.getColumnDimension(); column++) {
        double standard = m.getEntry(row, column) - mean[column];
        standard /= sds[column];
        m.setEntry(row, column, standard);
      }
    }
  }
}
