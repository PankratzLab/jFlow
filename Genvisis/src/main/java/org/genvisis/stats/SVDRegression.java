package org.genvisis.stats;

import java.util.Arrays;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;

/**
 * A class to compute betas and (partially)standard error of betas using singular value
 * decomposition
 * <p>
 * This is faster when using many (aprox. 150) independent variables
 * <p>
 * Methods are described here: http://www.jstor.org/stable/2684086
 *
 */
public class SVDRegression {
  private int numComponents;
  private final double[] deps;
  private double[] W;// num components
  private double[] a_hat;
  private double[] betas;
  private double[][] invP_Diagonal;
  private final double[][] indepsT;
  private double[][] U;// num components by num inds
  private double[][] V;// num components by num independent variables (full rank)
  private final boolean verbose;

  private PrincipalComponentsCompute principalComponentsCompute;
  private final Logger log;

  /**
   * @param deps
   * @param new_indeps organized as indeps[sample0][indep0...]
   * @param log Warning - no data checks are done and it is assumed that this is called after a
   *        regression model's QC steps
   */

  public SVDRegression(double[] deps, double[][] new_indeps, boolean verbose, Logger log) {
    super();
    this.deps = deps;
    indepsT = addConstant(Matrix.transpose(new_indeps));// transpose is neccesary;
    this.log = log;
    this.verbose = verbose;
    // for testing against the paper
    // this.deps =new double[]{41.38,31.01,37.41,50.05,39.17,38.86,46.14,44.47};
  }

  /**
   * Computes everything
   * <p>
   * Warning - no data checks are done and it is assumed that this is called after a regression
   * model's QC steps
   */
  public void svdRegression() {
    numComponents = indepsT.length;// number of independent variables, plus the constant = full rank
    betas = new double[numComponents];
    principalComponentsCompute =
        PrincipalComponentsCompute.getPrincipalComponents(numComponents, false, indepsT, verbose,
                                                          log);
    extractU();
    extractW();
    computeV();
    computeAhat();
    computeBetas();
    computeInvPDiagonal();
  }

  /**
   * @return a matrix with the diagonal (and only the diagonal) equalivalent to inverse p
   */
  public double[][] getInvP_Diagonal_Equivalent() {
    return invP_Diagonal;
  }

  public double[] getBetas() {
    return betas;
  }

  private void extractU() {// technically U transposed
    U = PrincipalComponentsCompute.getPCs(principalComponentsCompute, numComponents, verbose, log);
  }

  private void extractW() {
    W = principalComponentsCompute.getSingularValues();// W -singular values - must be sorted in
                                                       // descending order or else things blow up
  }

  private void computeV() {// components are the rows of the V matrix
    V = new double[numComponents][indepsT.length];
    for (int i = 0; i < indepsT.length; i++) {
      for (int j = 0; j < numComponents; j++) {
        V[i][j] = getVSingle(W[j], indepsT[i], U[j]);
      }
    }
  }

  /**
   * Equation 12 from the paper above
   */
  private void computeAhat() {// the "predictive linker" between dependent and independent variables
    a_hat = new double[numComponents];
    for (int i = 0; i < numComponents; i++) {
      a_hat[i] = 0;
      for (int j = 0; j < deps.length; j++) {
        a_hat[i] += deps[j] * U[i][j];
      }
    }
  }

  /**
   * Equation 18/19 from paper above
   */
  private void computeBetas() {// should be identical within precision to the normal equations
    for (int i = 0; i < numComponents; i++) {
      betas[i] = 0;
      for (int j = 0; j < numComponents; j++) {
        betas[i] += V[i][j] * (a_hat[j] / W[j]);
      }
    }
  }

  /**
   * Computes the first part of equation 20a from the paper above, does not multiply by sigma
   * squared to
   * <p>
   * be equivalent to the invP diagonal
   */
  private void computeInvPDiagonal() {
    invP_Diagonal = new double[numComponents][numComponents];
    for (int i = 0; i < numComponents; i++) {
      invP_Diagonal[i][i] = 0;
      for (int j = 0; j < numComponents; j++) {
        invP_Diagonal[i][i] += Math.pow(V[i][j], 2) / Math.pow(W[j], 2);
      }
    }
  }

  /**
   * Adds the constant term to a variable dominant indepsT
   */
  private static double[][] addConstant(double[][] indepsT) {
    double[][] tmpIndeps = new double[indepsT.length + 1][indepsT[0].length];
    Arrays.fill(tmpIndeps[0], 1);
    for (int i = 0; i < indepsT.length; i++) {
      tmpIndeps[i + 1] = indepsT[i];
    }
    return tmpIndeps;
  }

  /**
   * Compute a V matrix entry at a given point, data[] is variable dominant
   */
  private static double getVSingle(double singularValue, double[] data, double[] basis) {
    double sum = 0;
    for (int i = 0; i < data.length; i++) {
      sum += data[i] * basis[i];
    }
    return sum / singularValue;
  }

}
