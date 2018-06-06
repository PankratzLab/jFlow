/**
 * 
 */
package org.genvisis.common.matrix;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Attaches column and row names to the matrix
 */
public class NamedRealMatrix {

  //TODO conv to map String->index
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
