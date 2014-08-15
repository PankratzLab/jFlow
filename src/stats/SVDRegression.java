package stats;

import java.util.Arrays;

import cnv.analysis.pca.PrincipalComponentsCompute;
import common.Logger;
import common.Matrix;

/**
 * A class to compute betas and (partially)standard error of betas using singular value decomposition
 * <p>
 * This is faster when using many (aprox. 150) independent variables
 * <p>
 * Methods are described here: http://www.jstor.org/stable/2684086
 *
 */
public class SVDRegression {
	private int numComponents;
	private double[] deps;
	private double[] W;// num components
	private double[] a_hat;
	private double[] betas;
	private double[][] invP_Diagonal;
	private double[][] indepsT;
	private double[][] U;// num components by num inds
	private double[][] V;// num components by num independent variables (full rank)

	private PrincipalComponentsCompute principalComponentsCompute;
	private Logger log;

	/**
	 * @param deps
	 * @param new_indeps
	 *            organized as indeps[sample0][indep0...]
	 * @param log
	 *            Warning - no data checks are done and it is assumed that this is called after a regression model's QC steps
	 */

	public SVDRegression(double[] deps, double[][] new_indeps, Logger log) {
		super();
		this.deps = deps;
		this.indepsT = addConstant(Matrix.transpose(new_indeps));// transpose is neccesary;
		this.log = log;
		// for testing against the paper
		// this.deps =new double[]{41.38,31.01,37.41,50.05,39.17,38.86,46.14,44.47};
	}

	/**
	 * Computes everything
	 * <p>
	 * Warning - no data checks are done and it is assumed that this is called after a regression model's QC steps
	 */
	public void svdRegression() {
		this.numComponents = indepsT.length;// number of independent variables, plus the constant = full rank
		this.betas = new double[numComponents];
		this.principalComponentsCompute = PrincipalComponentsCompute.getPrincipalComponents(numComponents, false, indepsT, log);
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
		this.U = PrincipalComponentsCompute.getPCs(principalComponentsCompute, numComponents, log);
	}

	private void extractW() {
		this.W = principalComponentsCompute.getSingularValues();// W -singular values - must be sorted in descending order or else things blow up
	}

	private void computeV() {// components are the rows of the V matrix
		this.V = new double[numComponents][indepsT.length];
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
		this.a_hat = new double[numComponents];
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
	 * Computes the first part of equation 20a from the paper above, does not multiply by sigma squared to
	 * <p>
	 * be equivalent to the invP diagonal
	 */
	private void computeInvPDiagonal() {
		this.invP_Diagonal = new double[numComponents][numComponents];
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
