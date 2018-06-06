/**
 * 
 */
package org.genvisis.common.matrix;

import org.apache.commons.math3.linear.RealMatrix;

/**
 *
 *
 */
public interface MatrixDataLoading {

  /**
   * Loads the date (typically genotypes that will enter the ancestry PCA)
   * 
   * @return {@link RealMatrix}
   */
  NamedRealMatrix getData();

}
