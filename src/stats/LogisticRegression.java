package stats;

import java.util.*;

import common.*;

public class LogisticRegression extends RegressionModel {
	private int sY0, sY1;
	private String output;
	private double overallDF;
	private double[][] odds_ratios;
	private double logLikeFinal;
	private double logLikeNull;
	private int offset;
	private double CSRsquare;

	@SuppressWarnings("rawtypes")
	public LogisticRegression(Vector iDeps, Vector iIndeps) {
		this(iDeps, iIndeps, false, true);
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public LogisticRegression(Vector vDeps, Vector vIndeps, boolean bypassDataChecks, boolean verbose) {
		this(processDeps(vDeps), processIndeps(vIndeps), bypassDataChecks, verbose);
	}

	public LogisticRegression(int[] deps, double[][] indeps) {
		this(Array.toDoubleArray(deps), indeps);
	}

	public LogisticRegression(int[] deps, int[][] indeps) {
		this(Array.toDoubleArray(deps), Matrix.toDoubleArrays(indeps));
	}

	public LogisticRegression(int[] deps, int[][] indeps, boolean bypassDataChecks, boolean verbose) {
		this(Array.toDoubleArray(deps), Matrix.toDoubleArrays(indeps), bypassDataChecks, verbose);
	}

	public LogisticRegression(int[] deps, double[] indeps) {
		this(Array.toDoubleArray(deps), Matrix.toMatrix(indeps), false, true);
	}

	public LogisticRegression(double[] deps, double[] indeps) {
		this(deps, Matrix.toMatrix(indeps), false, true);
	}

	public LogisticRegression(double[] deps, double[][] indeps) {
		this(deps, indeps, false, true);
	}

	public LogisticRegression(double[] deps, double[][] indeps, boolean bypassDataChecks, boolean verbose) {
		this(deps, indeps, null, bypassDataChecks, verbose);
	}
	
	public LogisticRegression(double[] deps, double[][] indeps, String[] indepVariableNames, boolean bypassDataChecks, boolean verbose) {
		this.deps = deps;
		this.indeps = indeps;
		this.verbose = verbose;
		this.analysisFailed = false;
		this.logistic = true;
		
		varNames = new String[indeps[0].length+1];
		
		if (indeps.length > 0) {
			M = indeps[0].length;
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
		
		if (!bypassDataChecks) {
			checkForMissingData();
			analysisFailed = !dataQC();
		}

		compute();
	}

	public int ix(int j, int k, int nCols) {
		return j*nCols+k;
	}

	// detects offset, but is never used
//	public boolean dataChecksOut() {
//		IntVector iv = new IntVector();
//
//		for (int i = 0; i<deps.length; i++) {
//			if (!iv.contains((int)deps[i])) {
//				iv.add((int)deps[i]);
//			}
//		}
//		if (iv.size()==1) {
//			System.err.println("Error in logistic regression - Dependent variables must be 2 consecutive integers, not just '"+iv.elementAt(0)+"')");
//			return false;
//		} else if (iv.size()>2||Math.abs(iv.elementAt(0)-iv.elementAt(1))!=1) {
//			System.err.println("Error in logistic regression - Dependent variables must be 2 and no more than 2 consecutive integers");
//			return false;
//		}
//
//		offset = iv.elementAt(0)<iv.elementAt(1)?iv.elementAt(0):iv.elementAt(1);
//		if (offset!=0) {
//			for (int i = 0; i<deps.length; i++) {
//				deps[i] = deps[i]-offset;
//			}
//		}
//
//		return true;
//	}

	public void compute() {
		N = indeps.length; // number of cases
		if (N==0) {
			System.err.println("Error - cannot perform logistic regression with zero individuals!!");
			fail();
			return;
		}
		M = indeps[0].length;

		double x;
		double q;

		int nP = M+1;
		int nP1 = nP+1;
		sY0 = 0; // sum of YOs, i.e. number of unaffecteds
		sY1 = 0; // sum of Y1s, i.e. number of affecteds
		int sC = 0; // sum of YOs + Y1s, i.e. total number of cases

		double[] X = new double[N*(M+1)]; // arrays should be zeroed out
		int[] Y0 = new int[N];
		int[] Y1 = new int[N];
		double[] xM = new double[M+1]; // mean
		double[] xSD = new double[M+1]; // std dev
		double[] Par = new double[nP];
		double[] SEP = new double[nP];
		double[] Arr = new double[nP*nP1];

		for (int i = 0; i<N; i++) {
			X[ix(i, 0, M+1)] = 1;

			for (int j = 1; j<=M; j++) {
				X[ix(i, j, M+1)] = indeps[i][j-1];
			}
			// x = iDeps[i];
			x = deps[i];
			if (x==0) {
				Y0[i] = 1;
				sY0 = sY0+1;
			} else {
				Y1[i] = 1;
				sY1 = sY1+1;
			}
			sC = sC+(Y0[i]+Y1[i]);
			for (int j = 1; j<=M; j++) {
				x = X[ix(i, j, M+1)];
				xM[j] = xM[j]+(Y0[i]+Y1[i])*x;
				xSD[j] = xSD[j]+(Y0[i]+Y1[i])*x*x;
			}
		}

		output = ext.formDeci(sY0, 3)+"\t"+ext.formDeci(sY1, 3);

		for (int j = 1; j<=M; j++) {
			xM[j] = xM[j]/sC;
			xSD[j] = xSD[j]/sC;
			xSD[j] = Math.sqrt(Math.abs(xSD[j]-xM[j]*xM[j]));
			output += "\t"+ext.formDeci(xM[j]*sC, 3);
		}
		xM[0] = 0;
		xSD[0] = 1;

		for (int i = 0; i<N; i++) {
			for (int j = 1; j<=M; j++) {
				X[ix(i, j, M+1)] = ((X[ix(i, j, M+1)]-xM[j])/xSD[j]);
			}
		}

		Par[0] = Math.log((double)sY1/(double)sY0);
		for (int j = 1; j<=M; j++) {
			Par[j] = 0;
		}

		double LnV = 0;
		double Ln1mV = 0;

		double LLp = 2E10;
		double LL = 1E10;
		double LLn = 1E10;

		int count = 0;
		while (Math.abs(LLp-LL)>0.0000001 && count < 10000) {
			count++;
			LLp = LL;
			LL = 0;
			for (int j = 0; j<=M; j++) {
				for (int k = j; k<=M+1; k++) {
					Arr[ix(j, k, M+2)] = 0;
				}
			}

			for (int i = 0; i<N; i++) {
				double v = Par[0];
				for (int j = 1; j<=M; j++) {
					v = v+Par[j]*X[ix(i, j, M+1)];
				}
				if (v>15) {
					LnV = -Math.exp(-v);
					Ln1mV = -v;
					q = Math.exp(-v);
				} else {
					if (v<-15) {
						LnV = v;
						Ln1mV = -Math.exp(v);
						q = Math.exp(v);
					} else {
						v = 1/(1+Math.exp(-v));
						LnV = Math.log(v);
						Ln1mV = Math.log(1-v);
						q = v*(1-v);
					}
				}
				LL = LL-2*Y1[i]*LnV-2*Y0[i]*Ln1mV;
				for (int j = 0; j<=M; j++) {
					double xij = X[ix(i, j, M+1)];
					Arr[ix(j, M+1, M+2)] = Arr[ix(j, M+1, M+2)]+xij*(Y1[i]*(1-v)+Y0[i]*(-v));
					for (int k = j; k<=M; k++) {
						Arr[ix(j, k, M+2)] = Arr[ix(j, k, M+2)]+xij*X[ix(i, k, M+1)]*q*(Y0[i]+Y1[i]);
					}
				}
			}

			if (LLp==1e+10) {
				LLn = LL;
			}

			for (int j = 1; j<=M; j++) {
				for (int k = 0; k<j; k++) {
					Arr[ix(j, k, M+2)] = Arr[ix(k, j, M+2)];
				}
			}

			for (int i = 0; i<=M; i++) {
				double s = Arr[ix(i, i, M+2)];
				Arr[ix(i, i, M+2)] = 1;
				for (int k = 0; k<=M+1; k++) {
					Arr[ix(i, k, M+2)] = Arr[ix(i, k, M+2)]/s;
				}
				for (int j = 0; j<=M; j++) {
					if (i!=j) {
						s = Arr[ix(j, i, M+2)];
						Arr[ix(j, i, M+2)] = 0;
						for (int k = 0; k<=M+1; k++) {
							Arr[ix(j, k, M+2)] = Arr[ix(j, k, M+2)]-s*Arr[ix(i, k, M+2)];
						}
					}
				}
			}

			for (int j = 0; j<=M; j++) {
				Par[j] = Par[j]+Arr[ix(j, M+1, M+2)];
			}

		}
		logLikeNull = LLn;
		logLikeFinal = LL;

		double CSq = LLn-LL;
		overall = CSq;
		overallDF = M;
		overallSig = ProbDist.ChiDist(CSq, M);
		CSRsquare = 1-Math.exp(-1*(LLn-LL)/N);
		Rsquare = CSRsquare/(1-Math.exp(-1*LLn/N)); // Nagelkerke R
		// Square

		sigs = new double[nP];
		stats = new double[nP];
		SEofBs = new double[nP];
		betas = new double[nP];
		for (int j = 1; j<=M; j++) {
			Par[j] = Par[j]/xSD[j];
			SEP[j] = Math.sqrt(Arr[ix(j, j, nP+1)])/xSD[j];
			Par[0] = Par[0]-Par[j]*xM[j];
			if (count==10000) {
				sigs[j] = 99999999;
			} else {
				sigs[j] = ProbDist.NormDist(Math.abs(Par[j]/SEP[j]));
			}
			stats[j] = Math.pow(Par[j]/SEP[j], 2);
			SEofBs[j] = SEP[j];
			betas[j] = Par[j];
		}
		betas[0] = Par[0];

		computeOddsRatios();
		for (int j = 1; j<=M; j++) {
			output += "\t"+ext.formDeci(odds_ratios[j][0], 3, true)+"\t"+ext.formDeci(odds_ratios[j][1], 3, true)+"\t"+ext.formDeci(odds_ratios[j][2], 3, true);
		}

		predicteds = new double[N];
		residuals = new double[N];
		for (int i = 0; i<N; i++) {
			x = betas[0];
			for (int j = 0; j<M; j++) {
				x += indeps[i][j]*betas[j+1];
			}
			predicteds[i] = Math.exp(x)/(1+Math.exp(x));
			// residuals[i] = iDeps[i] - predicteds[i];
			residuals[i] = deps[i]-predicteds[i];
		}
	}

	public void computeOddsRatios() {
		odds_ratios = new double[M+1][3];

		for (int j = 1; j<=M; j++) {
			odds_ratios[j][0] = Math.exp(betas[j]);
			odds_ratios[j][1] = Math.exp(betas[j]-1.96*SEofBs[j]);
			odds_ratios[j][2] = Math.exp(betas[j]+1.96*SEofBs[j]);
		}
	}
	
	public int getY0() {
		return sY0;
	}

	public int getY1() {
		return sY1;
	}
	
	// output returns "num_unaffected num_affected
	// num_vars*(mean_var*num_records) num_vars*(OR lower_bound upper_bound)
	public String getOutput() {
		return output;
	}

	public double[] getProbabilities() {
		return predicteds;
	}

	public double[] getWalds() {
		return stats;
	}

	public double[][] getOddsRatios() {
		return odds_ratios;
	}

	public double[][] getEffectsAndConfidenceIntervals() {
		return odds_ratios;
	}

	public String getSummary() {
		String str = "";
		String eol;

		if (analysisFailed) {
			return "Did not run";
		}

		eol = Files.isWindows()?"\r\n":"\n";
		if (onePer) {
			str += "One per family was permuted "+numPermutations+" times"+eol;
			str += "Statistics were bootstrapped "+numBootReps+" times"+eol;
			str += ""+eol;
			str += "Average of "+ext.formDeci((double)logCounts[0]/numPermutations, 2)+" cases with Y="+offset+""+eol;
			str += "Average of "+ext.formDeci((double)logCounts[1]/numPermutations, 2)+" cases with Y="+(offset+1)+""+eol;
			str += "Number of independent observations: "+Array.unique(famIDs).length+""+eol;
			str += ""+eol;
		} else {
			str += sY0+" cases with Y="+offset+""+eol;
			str += sY1+" cases with Y="+(offset+1)+""+eol;
			str += "Total "+(sY0+sY1)+""+eol;
			str += ""+eol;
			str += "-2 Log likelihood = "+ext.formDeci(logLikeNull, 3)+" (Null)"+eol;
			str += "-2 Log likelihood = "+ext.formDeci(logLikeFinal, 3)+" (Converged)"+eol;
			str += "ChiSquare = "+ext.formDeci(overall, 3)+", df = "+overallDF+", p = "+ext.formDeci(overallSig, 3, true)+""+eol;
			str += "Cox & Snell R-square = "+ext.formDeci(CSRsquare, 3, true)+", Nagelkerke R-square = "+ext.formDeci(Rsquare, 4, true)+""+eol;
			str += ""+eol;
		}

		str += "Coefficients:"+eol;
		str += ext.formStr("Model", maxNameSize, true)+"\t   Beta\t StdErr\t   Wald\t   Sig.\t  O.R."+eol;
		str += modelSummary()+""+eol;

		if (!onePer) {
			// str += ""+delimiter;
			// for (int i = 0; i<5; i++) {
			// str += residuals[i]+""+delimiter;
			// }
		}

		return str;
	}

	public String modelSummary() {
		String str = "";
		String eol;
		
		if (analysisFailed) {
			return "Did not run";
		}

		eol = Files.isWindows()?"\r\n":"\n";
		for (int i = 1; i<betas.length; i++) {
			str += ext.formStr(varNames[i], maxNameSize, true)+"\t"+ext.formStr(ext.formDeci(betas[i], sigfigs, true), 7)+"\t"+ext.formStr(ext.formDeci(SEofBs[i], sigfigs, true), 7)+"\t"+ext.formStr(ext.formDeci(stats[i], sigfigs, true), 7)+"\t"+ext.formStr(ext.formDeci(sigs[i], sigfigs, true), 7)+"\t  "+ext.formDeci(odds_ratios[i][0], sigfigs, true)+" ("+ext.formDeci(odds_ratios[i][1], sigfigs, true)+", "+ext.formDeci(odds_ratios[i][2], sigfigs, true)+")"+eol;
			str += ext.formStr(varNames[i], maxNameSize, true)+"\t"+ext.formStr(ext.formDeci(betas[i], sigfigs, true), 7)+"\t"+ext.formStr(ext.formDeci(SEofBs[i], sigfigs, true), 7)+"\t"+ext.formStr((stats[i]>100?ext.formSciNot(stats[i], 1, true):ext.formDeci(stats[i], sigfigs, true)), 7)+"\t"+ext.formStr(ext.formDeci(sigs[i], sigfigs, true), 7)+"\t"+ext.prettyP(sigs[i])+"\t==CHIDIST("+Math.abs(stats[i])+",1) (untested)"+eol;
		}
		str += ext.formStr(varNames[0], maxNameSize, true)+"\t"+ext.formStr(ext.formDeci(betas[0], sigfigs, true), 7)+eol;

		return str;
	}

	public double[][] getAccuracy() {
		int guess;
		double max = -1;
		DoubleVector dv = new DoubleVector();
		int[] keys;
		double[][] space;
		int n = deps.length;

		double[][] accs = new double[2][4];

		int numAffs = 0;

		for (int i = 0; i<n; i++) {
			if (!dv.contains(predicteds[i])) {
				dv.add(predicteds[i]);
			}
			if ((int)deps[i]==1) {
				numAffs++;
			}
		}
		keys = Sort.quicksort(dv);
		space = new double[dv.size()-1][4];
		for (int i = 0; i<dv.size()-1; i++) {
			space[i][0] = dv.elementAt(keys[i])+(dv.elementAt(keys[i+1])-dv.elementAt(keys[i]))/2;
			for (int j = 0; j<n; j++) {
				guess = predicteds[j]>=space[i][0]?1:0;
				if ((int)deps[j]==guess) {
					space[i][1]++;
					space[i][2+(int)deps[j]]++;
				}
			}
			space[i][1] /= (double)n;
			if (space[i][1]>max) {
				max = space[i][1];
				accs[0] = space[i];
			}
			space[i][2] /= (double)(n-numAffs);
			space[i][3] /= (double)numAffs;
			if (accs[1][1]==0&&space[i][2]>space[i][3]) {
				accs[1] = space[i];
			}
		}
		return accs;
	}

	public String getAccuracySummary() {
		String str = "";
		double[][] accs = getAccuracy();
		String eol;
		
		eol = Files.isWindows()?"\r\n":"\n";

		str += "\t\tcut pt\taccuracy\tspecificity\tsensitivity"+eol;
		str += "Maxed:\t"+ext.formDeci(accs[0][0], 4, true)+"\t"+ext.formDeci(accs[0][1], 3)+"\t\t"+ext.formDeci(accs[0][2], 3)+"\t\t"+ext.formDeci(accs[0][3], 3)+eol;
		str += "Cross:\t"+ext.formDeci(accs[1][0], 4, true)+"\t"+ext.formDeci(accs[1][1], 3)+"\t\t"+ext.formDeci(accs[1][2], 3)+"\t\t"+ext.formDeci(accs[1][3], 3)+eol;

		return str;
	}
}
