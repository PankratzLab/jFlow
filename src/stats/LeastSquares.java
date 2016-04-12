package stats;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import common.Array;
import common.Matrix;
import common.ext;

public class LeastSquares extends RegressionModel {
	private double[][] X;
	private double[][] Y;
	private double meanY;
	private double[][] invP;
	private LS_TYPE lType;

	public enum LS_TYPE {
		/**
		 * The good ole regression that has been well tested
		 */
		REGULAR,
		/**
		 * Uses singular value decomposition to get the lin reg. <br>
		 * Basically only faster than {@link LS_TYPE#REGULAR} due to QR decomposition in EJML <br>
		 * Extra baggage from full SVD though
		 */
		SVD,
		/**
		 * math.commons Implementation
		 */
		OLS;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public LeastSquares(Vector vDeps, Vector vIndeps) { // deps = Vector of int/double as String, indeps = Vector of double[]
		this(vDeps, vIndeps, false, true);
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public LeastSquares(Vector vDeps, Vector vIndeps, LS_TYPE lType) { // deps = Vector of int/double as String, indeps = Vector of double[]
		this(vDeps, vIndeps, false, false, lType);
	}
	
	@SuppressWarnings({ "rawtypes" })
	public LeastSquares(Vector<String> vDeps, Vector vIndeps, boolean bypassDataCheck, boolean verbose) {
		this(processDeps(vDeps), processIndeps(vIndeps), bypassDataCheck, verbose);
	}
	
	@SuppressWarnings({ "rawtypes" })
	public LeastSquares(Vector<String> vDeps, Vector vIndeps, boolean bypassDataCheck, boolean verbose, LS_TYPE lType) {
		this(processDeps(vDeps), processIndeps(vIndeps), bypassDataCheck, verbose, lType);
	}
	
	public LeastSquares(double[] new_deps, double[][] new_indeps) {
		this(new_deps, new_indeps, false, true);
	}

	public LeastSquares(double[] new_deps, double[] new_indeps) {
		this(new_deps, Matrix.toMatrix(new_indeps), false, true);
	}

	public LeastSquares(int[] iDeps, int[][] iIndeps) {
		this(Array.toDoubleArray(iDeps), Matrix.toDoubleArrays(iIndeps));
	}

	public LeastSquares(double[] new_deps, double[][] new_indeps, boolean bypassDataCheck, boolean verbose) {
		this(new_deps, new_indeps, null, bypassDataCheck, verbose);
	}
	
	public LeastSquares(double[] new_deps, double[][] new_indeps, boolean bypassDataCheck, boolean verbose, LS_TYPE lType) {
		this(new_deps, new_indeps, null, bypassDataCheck, verbose, lType);
	}
	
	public LeastSquares(double[] new_deps, double[][] new_indeps, String[] indepVariableNames, boolean bypassDataCheck, boolean verbose) {
		this(new_deps, new_indeps, null, bypassDataCheck, verbose, LS_TYPE.REGULAR);
	}

	public LeastSquares(double[] new_deps, double[][] new_indeps, String[] indepVariableNames, boolean bypassDataCheck, boolean verbose, LS_TYPE lType) {
		this.deps = new_deps;
		this.indeps = new_indeps;
		this.verbose = verbose;
		this.analysisFailed = false;
		this.logistic = false;
		this.lType = lType;

		if (deps.length!=indeps.length) {
			System.err.println("Error - mismatched number of records: "+deps.length+" dependent elements and "+indeps.length+" independent elements");
			fail();
			return;
		}
		
		if (new_indeps.length > 0) {
			M = new_indeps[0].length;
			varNames = new String[M+1];
			varNames[0] = "Constant";
			for (int i = 1; i<M+1; i++) {
				varNames[i] = "Indep "+i;
			}
			maxNameSize = (M+1)<10?8:7+((M+1)+"").length();
			
			if (indepVariableNames != null) {
				setVarNames(indepVariableNames);
			}
		}


		if (!bypassDataCheck) {
			checkForMissingData();
			analysisFailed = !dataQC();
		}

		calc();
	}

	public double Determinant(double[][] matrix) {
		int tms = matrix.length;
		double det = 1;
		double f1 = 0;
		double temp = 0;
		int v = 1;
		int iDF = 1;

		for (int col = 0; col<tms-1; col++) {
			for (int row = col+1; row<tms; row++) {
				v = 1;

				outahere: while (matrix[col][col]==0) {
					if (col+v>=tms) {
						iDF = 0;
						break outahere;
					} else {
						for (int c = 0; c<tms; c++) {
							temp = matrix[col][c];
							matrix[col][c] = matrix[col+v][c];
							matrix[col+v][c] = temp;
						}
						v++;
						iDF = iDF*-1;
					}
				}

				if (matrix[col][col]!=0) {
					f1 = (-1)*matrix[row][col]/matrix[col][col];
					for (int i = col; i<tms; i++) {
						matrix[row][i] = f1*matrix[col][i]+matrix[row][i];
					}
				}
			}
		}

		for (int i = 0; i<tms; i++) {
			det = det*matrix[i][i];
		}

		det = det*iDF; // adjust w/ determinant factor

		return det;
	}

	public double[][] Inverse(double[][] a) {
		// Formula used to Calculate Inverse:
		// inv(A) = 1/det(A) * adj(A)
		int tms = a.length;
		double m[][] = new double[tms][tms];
		double mm[][] = Adjoint(a);
		double dd = 1/Determinant(a);

		for (int i = 0; i<tms; i++)
			for (int j = 0; j<tms; j++) {
				m[i][j] = dd*mm[i][j];
			}

		return m;
	}

	public double[][] Adjoint(double[][] a) {
		int tms = a.length;
		double m[][] = new double[tms][tms];
		int ii, jj, ia, ja;
		double det;

		for (int i = 0; i<tms; i++)
			for (int j = 0; j<tms; j++) {
				ia = ja = 0;

				double ap[][] = new double[tms-1][tms-1];

				for (ii = 0; ii<tms; ii++) {
					for (jj = 0; jj<tms; jj++) {

						if ((ii!=i)&&(jj!=j)) {
							ap[ia][ja] = a[ii][jj];
							ja++;
						}

					}
					if ((ii!=i)&&(jj!=j)) {
						ia++;
					}
					ja = 0;
				}

				det = Determinant(ap);
				m[i][j] = (double)Math.pow(-1, i+j)*det;
			}

		m = Transpose(m);

		return m;
	}

	public static double[][] Transpose(double[][] a) {
		double m[][] = new double[a[0].length][a.length];

		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				m[j][i] = a[i][j];
		return m;
	}

	public double[][] MultiplyMatrix(double[][] a, double[][] b) {
		double matrix[][] = new double[a.length][b[0].length];

		if (a[0].length!=b.length) {
			System.err.println("Matrices incompatible for multiplication");
			fail();
		}

		for (int i = 0; i<a.length; i++) {
			for (int j = 0; j<b[i].length; j++)
				matrix[i][j] = 0;
		}

		for (int i = 0; i<matrix.length; i++) {
			for (int j = 0; j<matrix[i].length; j++) {
				matrix[i][j] = calculateRowColumnProduct(a, i, b, j);
			}
		}

		return matrix;
	}

	public double calculateRowColumnProduct(double[][] A, int row, double[][] B, int col) {
		double product = 0;

		for (int i = 0; i<A[row].length; i++)
			product += A[row][i]*B[i][col];

		return product;
	}

	public void linregr() {
		if (N<M+1) {
			System.err.println("Error - your need more data points than you have variables");
			fail();
		}

		invP = Inverse(MultiplyMatrix(X, Transpose(X)));

		betas = MultiplyMatrix(MultiplyMatrix(Y, Transpose(X)), invP)[0];
	}

	public void calc() {
		double SSy, SSerr, S2;

		logistic = false;
		N = deps.length;
		if (N==0) {
			if (verbose) {
				System.err.println("Error - cannot perform a least squares regression with zero indivudals!!");
			}
			fail();
			return;
		}
		M = indeps[0].length;

		X = new double[M+1][N];
		Y = new double[][] {deps};

		meanY = 0;
		for (int i = 0; i<N; i++) {
			// if (i > 856) {
			// double hi = Y[0][i];
			// System.out.println(hi);
			// }

			meanY += Y[0][i];
			// if ((Y[0][i]+"").equals("NaN")) {
			// System.err.println("Error - record #"+(i+1)+" is NaN");
			// }
			X[0][i] = 1;
			for (int j = 1; j<=M; j++) {
				X[j][i] = indeps[i][j-1];
			}
		}
		meanY /= N;

		if (!analysisFailed) {
			maxNameSize = (M + 1) < 10 ? 8 : 7 + ((M + 1) + "").length();

			switch (lType) {
			case OLS:
				break;
			case REGULAR:
				linregr();
				break;
			case SVD:
				SVDRegression svdRegress = new SVDRegression(deps, indeps, verbose, log);
				svdRegress.svdRegression();
				betas = svdRegress.getBetas();
				invP = svdRegress.getInvP_Diagonal_Equivalent();// Warning! this is not a full invP matrix, only the diagonal elements are computed

				break;
			default:
				throw new IllegalArgumentException("Invalid regression type " + lType);

			}

//			meanRes = 0;
			predicteds = new double[N];
			for (int i = 0; i<N; i++) {
				predicteds[i] = betas[0];
				for (int j = 1; j<=M; j++)
					predicteds[i] += betas[j]*X[j][i];
//				meanRes += predicteds[i];
			}
//			meanRes /= N;

			residuals = new double[N];
			SSy = SSerr = 0;
			for (int i = 0; i<N; i++) {
				SSy += (Y[0][i]-meanY)*(Y[0][i]-meanY);
				residuals[i] = Y[0][i]-predicteds[i];
				SSerr += Math.pow(residuals[i], 2);
			}
			S2 = SSerr/(N-M-1);
			SEofBs = new double[M+1];
			stats = new double[M+1];
			sigs = new double[M+1];
			for (int i = 0; i<M+1; i++) {
				SEofBs[i] = Math.sqrt(S2)*Math.sqrt(invP[i][i]);
				stats[i] = betas[i]/SEofBs[i];
				sigs[i] = ProbDist.TDist(stats[i], N-M-1);
				if (sigs[i]<0) {
					if (sigs[i]<-1E-10) {
						System.out.println("Negative p-value: "+stats[i]+"\t"+(N-M-1)+"\t"+ProbDist.TDist(stats[i], N-M-1));
					}
					sigs[i] = 0;
				}
			}

			Rsquare = 1-SSerr/SSy;
			overall = (Rsquare/M)/((1-Rsquare)/(N-M-1)); // F
			// stat
			overallSig = ProbDist.FDist(overall, M, N-M-1);
			// System.out.println(betas[1]+"\t"+SEofBs[1]);
		} else {
			fail();
		}
	}

	public void setVarNames(String[] names, int maxOverride) {
		setVarNames(names);
		maxNameSize = maxOverride;
	}

	public double getF() {
		return overall;
	}

	public double getFsig() {
		return overallSig;
	}

	public double[] getTs() {
		return stats;
	}

	public String getEquation() {
		String output = "y = "+ext.formDeci(betas[0], sigfigs);

		for (int i = 1; i<=M; i++)
			output += " + "+ext.formDeci(betas[i], sigfigs)+"x"+i;

		return output;
	}

	public String getSummary() {
		String str = "";
		String eol;

		if (analysisFailed) {
			return "Did not run";
		}
		
		eol = System.getProperty("os.name").startsWith("Windows")?"\r\n":"\n";
		if (onePer) {
			str += "One per family was permuted "+numPermutations+" times"+eol;
			str += "Statistics were bootstrapped "+numBootReps+" times"+eol;
			str += ""+eol;
			str += "Number of independent observations: "+Array.unique(famIDs).length+""+eol;
			str += ""+eol;
		} else {
			str += "n="+N+" Fstat = "+ext.formDeci(overall, 3)+", Sig. "+ext.formDeci(overallSig, 3, true)+" and R2= "+ext.formDeci(Rsquare, 4)+""+eol;
			str += ""+eol;
		}
		str += "Coefficients:"+eol;
		str += ext.formStr("Model", maxNameSize, true)+"\t   Beta\t StdErr\t      t\t   Sig."+eol;
		str += modelSummary()+""+eol;

		return str;
	}

	public String modelSummary() {
		String str = "";
		String eol;

		eol = System.getProperty("os.name").startsWith("Windows")?"\r\n":"\n";
		for (int i = 0; i<betas.length; i++) {
			str += ext.formStr(varNames[i], maxNameSize, true)+"\t"+ext.formStr(ext.formDeci(betas[i], sigfigs, true), 7)+"\t"+ext.formStr(ext.formDeci(SEofBs[i], sigfigs, true), 7)+"\t"+ext.formStr((stats[i]>100?ext.formSciNot(stats[i], 1, true):ext.formDeci(stats[i], sigfigs, true)), 7)+"\t"+ext.formStr(ext.formDeci(sigs[i], sigfigs, true), 7)+"\t"+ext.prettyP(sigs[i])+"\t=TDIST("+Math.abs(stats[i])+","+((onePer?Array.unique(famIDs).length:N)-2)+",2)"+eol;
		}

		return str;
	}
	
	public double[][] getEffectsAndConfidenceIntervals() {
		double[][] statsAndConfidenceIntervals;
		
		statsAndConfidenceIntervals = new double[betas.length][3];
		for (int i = 0; i<betas.length; i++) {
			statsAndConfidenceIntervals[i][0] = betas[i];
			statsAndConfidenceIntervals[i][1] = betas[i]-1.96*SEofBs[i];
			statsAndConfidenceIntervals[i][2] = betas[i]+1.96*SEofBs[i];
        }
		
		return statsAndConfidenceIntervals;
	}

	public void dumpData(String filename) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(filename));
			writer.print("Dep");
			for (int i = 1; i<varNames.length; i++) {
				writer.print("\t"+varNames[i]);
			}
			writer.println();
			for (int i = 0; i<N; i++) {
				writer.print(Y[0][i]);
				for (int j = 1; j<M+1; j++) {
					writer.print("\t"+X[j][i]);
				}
				writer.println();
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing dump file: "+filename);
		}
	}
	
	public void destroy() {
		X = null;
		Y = null;
		invP = null;
	}
	
	
}
