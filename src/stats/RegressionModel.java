package stats;

import java.io.*;
import java.util.*;

import common.Array;
import common.HashVec;
import common.IntVector;
import common.Logger;
import common.ext;

public abstract class RegressionModel {
	protected boolean verbose;
	protected int M;
	protected int N;
	protected double[] deps;
	protected double[][] indeps;
	protected double[] betas;
	protected double[] SEofBs;
	protected double[] stats; // T statistic for linear, Wald statistic for logistic
	protected double[] sigs;
	protected double overall; // F statistic for linear, chi-square(?) for logistic
	protected double overallSig;
	protected double[] predicteds;
	protected double[] residuals;
	protected double Rsquare;
	protected String[] varNames;
	protected int maxNameSize;
	protected boolean logistic;
	protected boolean onePer;
	protected int numPermutations;
	protected int numBootReps;
	protected int[] logCounts;
	protected String[] famIDs;
	protected int failures;
	protected boolean analysisFailed;
	
	public RegressionModel() {
		verbose = true;
	}

	public abstract String modelSummary();
	
	public abstract String getSummary();

	public abstract double[][] getEffectsAndConfidenceIntervals();

	public boolean isLogistic() {
		return logistic;
	}
	
	public double[] getFinalDependentVariables() {
		return deps;
	}

	public double[][] getFinalIndependentVariables() {
		return indeps;
	}

	public double[] getBetas() {
		return betas;
	}

	public double[] getSEofBs() {
		return SEofBs;
	}

	public double[] getStats() {
		return stats;
	}

	public double[] getSigs() {
		return sigs;
	}

	public double getOverall() {
		return overall;
	}

	public double getOverallSig() {
		return overallSig;
	}

	public double[] getPredicteds() {
		return predicteds;
	}

	public double[] getResiduals() {
		return residuals;
	}

	public double getRsquare() {
		return Rsquare;
	}

	public String[] getVarNames() {
		return varNames;
	}

	public int getNumFailures() {
		return failures;
	}
	
	public boolean analysisFailed() {
		for (int i = 0; i<betas.length; i++) {
			if (Double.isNaN(betas[i])) {
				return true;
			}
		}
		return analysisFailed;
	}
	
	protected void fail() {
		betas = new double[M+1];
		SEofBs = new double[M+1];
		stats = new double[M+1];
		sigs = new double[M+1];
		predicteds = new double[N];
		residuals = new double[N];
		for (int i = 0; i<sigs.length; i++) {
			betas[i] = SEofBs[i] = stats[i] = sigs[i] = Double.NaN;
		}
		for (int i = 0; i<N; i++) {
			predicteds[i] = residuals[i] = Double.NaN;
		}

		Rsquare = overall = overallSig = Double.NaN;
		analysisFailed = true;
	}
	

	public void setVarNames(String[] names, int maxOverride) {
		setVarNames(names);
		maxNameSize = maxOverride;
	}

	public void setVarNames(String[] names) {
		if (names.length!=M) {
			System.err.println("Error naming independent variables: "+M+" variables, and only "+names.length+" names");
			return;
		}
		varNames = new String[M+1];
		varNames[0] = "Constant";
		maxNameSize = 8;
		for (int i = 0; i<M; i++) {
			varNames[i+1] = names[i];
			if (names[i].length()>maxNameSize) {
				maxNameSize = names[i].length();
			}
		}
	}

	public void dumpData(String filename) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(filename));
			writer.print("Dep");
			for (int i = 1; i<varNames.length; i++) {
				writer.print("\t"+varNames[i]);
			}
			writer.println(famIDs==null?"":"\tFamID");
			for (int i = 0; i<N; i++) {
				writer.print(ext.formDeci(deps[i], 5));
				for (int j = 0; j<M; j++) {
					writer.print("\t"+indeps[i][j]);
				}
				writer.println(famIDs==null?"":"\t"+famIDs[i]);
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing dump file: "+filename);
		}
	}

	public void onePerFamily(String[] famIDs, int numReps, int bootReps) {
		double[][] statPerms = new double[stats.length][numReps];
		double[][] betaPerms = new double[betas.length][numReps];
		onePer = true;
		String[] fams = Array.unique(famIDs);
		Hashtable<String,IntVector> hash = new Hashtable<String,IntVector>();
		IntVector iv;
		double[] rDeps = new double[fams.length];
		double[][] rIndeps = new double[fams.length][];
		int index;
		RegressionModel model;
		double[] rBetas, rStats;
		String progress = "0";

		this.famIDs = famIDs;

		analysisFailed = true;
		numPermutations = numReps;
		numBootReps = bootReps;
		logCounts = new int[2];

		for (int i = 0; i<fams.length; i++) {
			hash.put(fams[i], new IntVector());
		}

		overall = Double.NaN;
		overallSig = Double.NaN;
		Rsquare = Double.NaN;
		for (int i = 0; i<famIDs.length; i++) {
			predicteds[i] = Double.NaN;
			residuals[i] = Double.NaN;
			((IntVector)hash.get(famIDs[i])).add(i);
		}

		if (verbose) {
			System.out.print("Rep 0");
		}
		failures = 0;
		for (int i = 0; i<numPermutations; i++) {
			if (verbose&&(i+1)%(numReps/10)==0) {
				for (int j = 0; j<progress.length(); j++) {
					System.out.print("\b");
				}
				progress = (i+1)+(failures>0?", "+failures+" failures":"");
				System.out.print(progress);
			}

			do {
				for (int j = 0; j<fams.length; j++) {
					iv = (IntVector)hash.get(fams[j]);
					index = iv.elementAt((int)(Math.random()*iv.size()));
					rDeps[j] = deps[index];
					rIndeps[j] = indeps[index];
				}

				if (logistic) {
					model = new LogisticRegression(rDeps, rIndeps, true, true);
					logCounts[0] += ((LogisticRegression)model).getY0();
					logCounts[1] += ((LogisticRegression)model).getY1();
				} else {
					model = new LeastSquares(rDeps, rIndeps);
				}
				rBetas = model.getBetas();
				rStats = model.getStats();
				failures += model.analysisFailed()?1:0;
				for (int j = 1; j<=M; j++) {
					betaPerms[j][i] = rBetas[j];
					statPerms[j][i] = rStats[j];
				}
			} while (model.analysisFailed());
		}
		if (verbose) {
			System.out.println("");
		}

		for (int i = 1; i<=M; i++) {
			betas[i] = Array.bootstrap(betaPerms[i], numBootReps, verbose)[0];
			stats[i] = Array.bootstrap(statPerms[i], numBootReps, verbose)[0];
			sigs[i] = logistic?ProbDist.ChiDist(stats[i], 1):ProbDist.TDist(stats[i], N-M-1);
			SEofBs[i] = Double.NaN;
		}
		if (logistic) {
			((LogisticRegression)this).computeOddsRatios();
		}

		analysisFailed = false;
	}
	
	public void checkForMissingData() {
		boolean[] use = new boolean[deps.length];
		int countNoDeps = 0;
		int countNoIndeps = 0;
		int count;
		double[] newDeps;
		double[][] newIndeps;

		for (int i = 0; i<use.length; i++) {
			use[i] = true;
			if ((deps[i]+"").equals("NaN")) {
				countNoDeps++;
				use[i] = false;
			}
			for (int j = 0; j<indeps[i].length; j++) {
				if (use[i]&&(indeps[i]+"").equals("NaN")) {
					countNoIndeps++;
					use[i] = false;
				}
			}
		}

		count = Array.booleanArraySum(use);
		newDeps = new double[count];
		newIndeps = new double[count][];
		count = 0;
		for (int i = 0; i<use.length; i++) {
			if (use[i]) {
				newDeps[count] = deps[i];
				newIndeps[count] = indeps[i];
				count++;
			}
		}

		if (verbose) {
			if (countNoDeps>0) {
				System.err.println(countNoDeps+" individuals were dropped because they lacked the dependent variable");
			}
			if (countNoIndeps>0) {
				System.err.println(countNoIndeps+" individuals were dropped because they lacked at least one independent variable");
			}
		}

		deps = newDeps;
		indeps = newIndeps;
	}
	
	public boolean dataQC() {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		double[][] oldIndeps;
		double diff, offset;
		String[] keys;
		
		if (deps==null) {
			System.err.println("Array with dependent variables is null");
			return false;
		}

		if (indeps==null) {
			System.err.println("Array with independent variables is null");
			return false;
		}

		for (int i = 0; i<deps.length && hash.size()<3; i++) {
			if (!hash.containsKey(deps[i]+"")) {
				hash.put(deps[i]+"", "");
			}
		}
		if (hash.size()<2) {
			System.err.println("No variance in the dependent variable");
			return false;
		}
		if (logistic) {
			if (hash.size()>2) {
				System.err.println("Error in logistic regression - Dependent variables must be 2 and no more than 2 consecutive integers");
				return false;
			}
			keys = HashVec.getKeys(hash);
			diff = Double.parseDouble(keys[0])-Double.parseDouble(keys[1]);
			if (Math.abs(diff) != 1) {
				System.err.println("Error in logistic regression - Dependent variables must be 2 and no more than 2 consecutive integers");
				return false;
			}

			offset = diff<0?Double.parseDouble(keys[0]):Double.parseDouble(keys[1]);
			if (offset!=0) {
				for (int i = 0; i<deps.length; i++) {
					deps[i] = deps[i]-offset;
				}
			}
		}
		
		for (int j = 0; j<indeps[0].length; j++) {
			hash.clear();
			for (int i = 0; i<indeps.length && hash.size()<2; i++) {
				if (!hash.containsKey(indeps[i][j]+"")) {
					hash.put(indeps[i][j]+"", "");
				}
			}
			if (hash.size()<2) {
				if (verbose) {
					System.err.println("No variance in independent variable number "+(j+1)+"; collapsing and ignoring");
				}
				oldIndeps = indeps.clone();
				indeps = new double[oldIndeps.length][oldIndeps[0].length-1];
				for (int k = 0; k<j; k++) {
					for (int i = 0; i<oldIndeps.length; i++) {
						indeps[i][k] = oldIndeps[i][k];
                    }
                }
				for (int k = j+1; k<oldIndeps[0].length; k++) {
					for (int i = 0; i<oldIndeps.length; i++) {
						indeps[i][k-1] = oldIndeps[i][k];
                    }
                }
				j--;
			}
		}

		return indeps[0].length>0?true:false;
	}

	public static double[] processDeps(Vector<String> vDeps) {
		double[] deps = new double[vDeps.size()];

		for (int i = 0; i<deps.length; i++) {
			deps[i] = Double.parseDouble(vDeps.elementAt(i));
		}

		return deps;
	}
	
	@SuppressWarnings({ "rawtypes" })
	public static double[][] processIndeps(Vector vIndeps) {
		double[][] indeps = new double[vIndeps.size()][];
		int[] intarray = {0};
		
		for (int i = 0; i<indeps.length; i++) {
			if (vIndeps.elementAt(i).getClass()==intarray.getClass()) {
				indeps[i] = Array.toDoubleArray((int[])vIndeps.elementAt(i));
			} else {
				indeps[i] = (double[])vIndeps.elementAt(i);
			}
		}
		
		return indeps;
	}
	
	public static RegressionModel determineAppropriate(double[] deps, double[] indeps, boolean bypassDataCheck, boolean verbose) {
		double[][] newIndeps = new double[indeps.length][1];
		
		for (int i = 0; i<indeps.length; i++) {
			newIndeps[i][0] = indeps[i]; 
        }
		
		return determineAppropriate(deps, newIndeps, bypassDataCheck, verbose);
	}
	
	public static RegressionModel determineAppropriate(double[] deps, double[][] indeps, boolean bypassDataCheck, boolean verbose) {
		Vector<String> depCount = new Vector<String>();
		
		for (int i = 0; i<deps.length && depCount.size()<=2; i++) {
			if (!(deps[i]+"").equals("NaN")) {
				HashVec.addIfAbsent(deps[i]+"", depCount);
			}
        }
		
		return depCount.size()==2?new LogisticRegression(deps, indeps, bypassDataCheck, verbose):new LeastSquares(deps, indeps, bypassDataCheck, verbose);		
	}
	
	public static boolean isBinaryTrait(String[] pheno, Logger log) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] keys;
		double diff;

		hash = new Hashtable<String,String>();
		if (pheno==null) {
			log.reportError("Error - null array");
			return false;
		}

		for (int i = 0; i<pheno.length && hash.size()<4; i++) {
			if (!hash.containsKey(pheno[i]+"") && !pheno[i].equalsIgnoreCase("NaN") && !pheno[i].equalsIgnoreCase("NA") && !pheno[i].equals(".")) {
				hash.put(pheno[i]+"", "");
			}
		}
		
		if (hash.size()<2) {
			log.reportError("No variance in the dependent variable");
			return false;
		}
		
		if (hash.size() == 3 && hash.containsKey("0")) {
			log.report("Warning - phenotype would be binary, if it weren't for the zero values; not flagged for logistic");
		}

		if (hash.size() == 2) {
			keys = HashVec.getKeys(hash);
			diff = Double.parseDouble(keys[0])-Double.parseDouble(keys[1]);
			if (Math.abs(diff) != 1) {
				System.err.println("Error - only two values, but they are not consecutive integers; not flagged for logistic");
				return false;
			}
		}
		
		return hash.size() == 2;		
	}	
}
