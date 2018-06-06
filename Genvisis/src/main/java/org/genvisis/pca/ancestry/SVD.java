/**
 * 
 */
package org.genvisis.pca.ancestry;

import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * 
 *
 */
public class SVD {

  private static final long serialVersionUID = 5517545544997813383L;
  private final RealMatrix m;
  private final String[] colNames;
  private final String[] rowNames;

  //  With V,W, and original data M we can always compute U
  private RealMatrix v;
  private DiagonalMatrix w;

  /**
   * @param m
   * @param colNames
   * @param rowNames
   * @param v
   * @param w
   */
  public SVD(RealMatrix m, String[] colNames, String[] rowNames, RealMatrix v, DiagonalMatrix w) {
    super();
    this.m = m;
    this.colNames = colNames;
    this.rowNames = rowNames;
    this.v = v;
    this.w = w;
  }


}
