/**
 * 
 */
package org.genvisis.cnv.analysis.pca.ancestry;

import org.apache.commons.math3.linear.RealMatrix;
import org.genvisis.common.RealMatrixUtils.NamedRealMatrix;

/**
 *
 *
 */
interface AncestryDataLoading {

  /**
   * Loads the date (typically genotypes that will enter the ancestry PCA)
   * 
   * @return {@link RealMatrix}
   */
  NamedRealMatrix getData();

}
