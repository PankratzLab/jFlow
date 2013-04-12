package stats;

import java.io.*;
import java.util.*;

import common.IntVector;
import common.ext;

public class Stepwise {
	public static double ENTRY_PROB = 0.05;
//	public static double ENTRY_PROB = 0.10;
	public static double REMOVAL_PROB = 0.10;
	private int M, N;
	private Vector<double[]> Xs;
	private Vector<String> Ys;
	private boolean logistic;
	private String[] varNames;
	private int maxNameSize;
	private Vector<IntVector> increments;
	private RegressionModel finalModel;

	public Stepwise(Vector<String> deps, Vector<double[]> indeps) {
		Xs = indeps;
		Ys = deps;
		run();
	}

	public Stepwise(double[] new_deps, double[][] new_indeps) {
		Xs = new Vector<double[]>();
		Ys = new Vector<String>();

		for (int i = 0; i<new_deps.length; i++) {
			Xs.add(new_indeps[i]);
			Ys.add(new_deps[i]+"");
		}

		run();
	}

	public Stepwise(int[] new_deps, int[][] new_indeps) {
		Xs = new Vector<double[]>();
		Ys = new Vector<String>();
		double[] data;

		for (int i = 0; i<new_deps.length; i++) {
			data = new double[new_indeps[i].length];
			for (int j = 0; j<new_indeps[i].length; j++) {
				data[j] = new_indeps[i][j];
			}
			Xs.add(data);
			Ys.add(new_deps[i]+"");
		}

		run();
	}

	public static RegVectors procFile(String filename) {
		BufferedReader reader = null;
		String[] line;
		Vector<String> deps = new Vector<String>();
		Vector<double[]> indeps = new Vector<double[]>();
		double[] data;
		boolean names = false;
		Vector<String> indNames = new Vector<String>();

		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.mark(5000);
			line = reader.readLine().split("[\\s]+");
			for (int i = 1; i<line.length; i++) {
				try {
					indNames.add(line[i]);
					Double.parseDouble(line[i]);
				} catch (NumberFormatException nfe) {
					names = true;
				}
			}
			if (!names) {
				reader.reset();
			}
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				deps.add(line[0]);
				data = new double[line.length-1];
				for (int i = 1; i<line.length; i++) {
					data[i-1] = Double.valueOf(line[i]).doubleValue();
				}
				indeps.add(data);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(1);
		}

		return (names?new RegVectors(deps, indeps, indNames):new RegVectors(deps, indeps));
	}

	public void run() {
		RegressionModel model;
		IntVector in = new IntVector();
		IntVector out = new IntVector();
		Vector<String> v = new Vector<String>();
		double[][] pvals;
		int bestModel;
		boolean done;
		double lowestP, highestRsq;

		if (Xs.size()!=Ys.size()) {
			System.err.println("Error using Vectors for stepwise regression; "+Ys.size()+" dependent elements and "+Xs.size()+" independent elements");
			System.exit(11);
		}

		N = Ys.size();
		for (int i = 0; i<N&&v.size()<3; i++) {
			if (!v.contains(Ys.elementAt(i))) {
				v.add(Ys.elementAt(i));
			}
		}

		if (v.size()==1) {
			System.err.println("Error - dependent variable has zero variation ("+v.elementAt(0)+")");
			System.exit(5);
		}
		logistic = (v.size()==2);
		System.out.println("Performing stepwise "+(logistic?"logistic":"linear")+" regression");

		M = Xs.firstElement().length;

		varNames = new String[M];
		for (int i = 0; i<M; i++) {
			varNames[i] = "Indep "+(i+1);
		}
		maxNameSize = M<10?8:7+(M+"").length();

		for (int i = 0; i<M; i++) {
			out.add(i);
		}

		done = false;
		increments = new Vector<IntVector>();
		while (!done) {
			pvals = new double[out.size()][];
			bestModel = -1;
			highestRsq = 0;
			lowestP = 1;
			for (int i = 0; i<out.size(); i++) {
				in.add(out.popFirst());
				if (logistic) {
					model = new LogisticRegression(Ys, travXs(in));
				} else {
					model = new LeastSquares(Ys, travXs(in));
				}
				if (!model.analysisFailed()&&model.getSigs().length==in.size()+1) {
					pvals[i] = model.getSigs();
					if (pvals[i][pvals[i].length-1]<lowestP) {
						lowestP = pvals[i][pvals[i].length-1];
					}
					if (model.getRsquare()>highestRsq) {
						highestRsq = model.getRsquare();
						bestModel = i;
					}
				}
				out.add(in.popLast());
			}
			if (lowestP<ENTRY_PROB) {
				in.add(out.popAt(bestModel));
				increments.add(in.clone());
				for (int i = in.size(); i>=1; i--) {
					if (pvals[bestModel][i]>REMOVAL_PROB) {
						in.removeElementAt(i-1);
						done = true;
					}
					if (done) {
						increments.add(in.clone());
						done = false;
					}
				}
			} else {
				done = true;
			}
		}

		if (increments.size()>0) {
			if (logistic) {
				finalModel = new LogisticRegression(Ys, travXs(in));
			} else {
				finalModel = new LeastSquares(Ys, travXs(in));
			}
		}
	}

	public Vector<double[]> travXs(IntVector ins) {
		Vector<double[]> v = new Vector<double[]>();
		double[] data, newdata;

		for (int i = 0; i<N; i++) {
			newdata = new double[ins.size()];
			data = Xs.elementAt(i);
			for (int j = 0; j<ins.size(); j++) {
				newdata[j] = data[ins.elementAt(j)];
			}
			v.add(newdata);
		}

		return v;
	}

	public void setVarNames(Vector<String> vNames) {
		String[] names = new String[vNames.size()];
		for (int i = 0; i<names.length; i++) {
			names[i] = vNames.elementAt(i);
		}
		setVarNames(names);
	}

	public void setVarNames(String[] names) {
		if (names.length!=M) {
			System.err.println("Error naming independent variables: "+M+" variables, and only "+names.length+" names");
			return;
		}
		varNames = names;
		maxNameSize = 8;
		for (int i = 0; i<M; i++) {
			if (names[i].length()>maxNameSize) {
				maxNameSize = names[i].length();
			}
		}
	}

	public String getSummary() {
		RegressionModel model;
		String Rsum = " Model\t"+(logistic?" ChiSq":"    F")+"\t   Sig\t R-square\n";
		String ModelSum = ext.formStr("Variable", maxNameSize, true)+"\t   Beta\t StdErr\t      T\t    Sig\n";
		IntVector ins;
		String[] travNames;

		for (int i = 0; i<increments.size(); i++) {
			ins = increments.elementAt(i);
			travNames = new String[ins.size()];
			for (int j = 0; j<ins.size(); j++) {
				travNames[j] = varNames[ins.elementAt(j)];
			}
			if (logistic) {
				model = new LogisticRegression(Ys, travXs(ins));
			} else {
				model = new LeastSquares(Ys, travXs(ins));
			}
			model.setVarNames(travNames, maxNameSize);
			Rsum += ext.formStr(i+1+"", 4)+"\t"+(Double.isInfinite(model.getOverall())?"    .":ext.formStr(ext.formDeci(model.getOverall(), 1, true), 7))+"\t  "+ext.formDeci(model.getOverallSig(), 3, true)+"\t  "+ext.formDeci(model.getRsquare(), 3, true)+"\n";
			ModelSum += "------ Model "+(i+1)+(i<10?" -":" ")+"---------------------------\n"+model.modelSummary();
		}

		return Rsum+"\n"+ModelSum;
	}

	public String getAccuracySummary() {
		return increments.size()==0?"no variables in final model":(logistic?((LogisticRegression)finalModel).getAccuracySummary():"");
	}

	public double[] getFinalPredicteds() {
		return increments.size()==0?new double[N]:finalModel.getPredicteds();
	}

	public double[] getFinalResiduals() {
		return increments.size()==0?new double[N]:finalModel.getResiduals();
	}

	public String getFinalNames() {
		IntVector ins;
		String finalNames = "";

		if (increments.size()==0) {
			return "no variables in final model";
		}

		ins = increments.lastElement();
		for (int j = 0; j<ins.size(); j++) {
			finalNames += (j==0?"":"\t")+varNames[ins.elementAt(j)];
		}

		return finalNames;
	}

	public void dumpData(String filename) {
		double[] data;

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(filename));
			writer.print("Dep");
			for (int i = 0; i<varNames.length; i++) {
				writer.print("\t"+varNames[i]);
			}
			writer.println();
			for (int i = 0; i<N; i++) {
				writer.print(Ys.elementAt(i));
				data = Xs.elementAt(i);
				for (int j = 0; j<M; j++) {
					writer.print("\t"+data[j]);
				}
				writer.println();
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing dump file: "+filename);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "vars.txt";
		Stepwise sw;

		String usage = "\n"+"park.stepwise requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			RegVectors rvs = procFile(filename);
			sw = new Stepwise(rvs.getDeps(), rvs.getIndeps());
			if (rvs.getVarNames()!=null) {
				sw.setVarNames(rvs.getVarNames());
			}
			System.out.println(sw.getSummary());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class RegVectors {
	private Vector<String> deps;

	private Vector<double[]> indeps;

	private Vector<String> varNames;

	public RegVectors(Vector<String> deps, Vector<double[]> indeps) {
		this.deps = deps;
		this.indeps = indeps;
		this.varNames = null;
	}

	public RegVectors(Vector<String> deps, Vector<double[]> indeps, Vector<String> varNames) {
		this.deps = deps;
		this.indeps = indeps;
		this.varNames = varNames;
	}

	public Vector<String> getDeps() {
		return deps;
	}

	public Vector<double[]> getIndeps() {
		return indeps;
	}

	public Vector<String> getVarNames() {
		return varNames;
	}
}
