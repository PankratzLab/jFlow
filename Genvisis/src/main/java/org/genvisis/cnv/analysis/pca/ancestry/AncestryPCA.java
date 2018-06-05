/**
 * 
 */
package org.genvisis.cnv.analysis.pca.ancestry;

import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.SingularOps;
import org.genvisis.common.Logger;
import org.genvisis.common.RealMatrixUtils;
import org.genvisis.common.RealMatrixUtils.NamedRealMatrix;

/**
 * Runs PCA starting from a {@link NamedRealMatrix}
 */
public class AncestryPCA {

  private AncestryPCA() {

  }

  public static void computePCA(NamedRealMatrix m, Logger log) {
    DenseMatrix64F a = RealMatrixUtils.toDenseMatrix64F(m.getM());

    log.reportTimeInfo("Computing EJML PCs");
    SvdImplicitQrDecompose_D64 svd = new SvdImplicitQrDecompose_D64(false, false, true, false);
    svd.decompose(a);
    log.reportTimeInfo("Finished Computing EJML PCs");

    log.reportTimeInfo("finished computing SVD base");

    DenseMatrix64F tv = svd.getV(null, true);

    DenseMatrix64F tmpW = svd.getW(null);
    SingularOps.descendingOrder(null, false, tmpW, tv, true);
    int numSingular = Math.min(tmpW.numRows, tmpW.numCols);
    double[] singularValues = new double[numSingular];
    for (int i = 0; i < numSingular; i++) {
      singularValues[i] = tmpW.get(i, i);
    }
    DiagonalMatrix w = new DiagonalMatrix(singularValues);
    RealMatrix v = MatrixUtils.createRealMatrix(tv.numRows, tv.numCols);
    for (int row = 0; row < tv.numRows; row++) {
      for (int col = 0; col < tv.numCols; col++) {
        v.addToEntry(row, col, tv.get(row, col));
      }
    }
  }

}