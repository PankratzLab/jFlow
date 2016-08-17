package org.genvisis.one.JL;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.genvisis.common.Logger;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.LeastSquares.LS_TYPE;

public class LinRegTest {

  public static void main(String[] args) {
    run();

  }

  private static void run() {

    Logger log = new Logger();
    int numSamps = 10000;
    int numVar = 100;
    double[] y = new double[numSamps];
    double[][] x2 = new double[numSamps][numVar];
    for (int i = 0; i < 10000; i++) {
      y[i] = Math.random();
      for (int j = 0; j < x2[i].length; j++) {
        x2[i][j] = Math.random();
      }

    }

    long time = System.currentTimeMillis();
    LeastSquares lsols = new LeastSquares(y, x2, true, true, LS_TYPE.QR_DECOMP);
    log.reportTimeElapsed(LS_TYPE.QR_DECOMP.toString() + ": ", time);
    System.out.println(lsols.getRsquare());

    time = System.currentTimeMillis();
    LeastSquares ls = new LeastSquares(y, x2, true, true, LS_TYPE.REGULAR);
    log.reportTimeElapsed(LS_TYPE.REGULAR.toString() + ": ", time);
    System.out.println(ls.getRsquare());
    time = System.currentTimeMillis();
    LeastSquares lssvd = new LeastSquares(y, x2, true, true, LS_TYPE.SVD);
    log.reportTimeElapsed(LS_TYPE.SVD.toString() + ": ", time);
    System.out.println(lssvd.getRsquare());
    time = System.currentTimeMillis();



    OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
    ols.newSampleData(y, x2);

    // ols.estimateRegressionParameters();
    System.out.println(ols.calculateRSquared());
    log.reportTimeElapsed(LS_TYPE.QR_DECOMP.toString() + " 2 ", time);

    for (int i = 0; i < ls.getBetas().length; i++) {
      if (i < 10) {
        // System.out.println(LS_TYPE.REGULAR + ": " + ls.getBetas()[i] + "\t" + LS_TYPE.SVD + ": "
        // + lssvd.getBetas()[i] + "\tOLS: " + lssvd.getBetas()[i]);
        // System.out.println(LS_TYPE.REGULAR + ": " + ls.getSigs()[i] + "\t" + LS_TYPE.SVD + ": " +
        // lssvd.getSigs()[i] + "\tOLS: " + lssvd.getSigs()[i]);
      }

    }

  }
}
