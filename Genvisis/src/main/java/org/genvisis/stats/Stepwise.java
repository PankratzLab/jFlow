package org.genvisis.stats;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.Producer;

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

	public Stepwise(Vector<String> deps, Vector<double[]> indeps, boolean bonferroni) {
		Xs = indeps;
		Ys = deps;
		run(Integer.MAX_VALUE,bonferroni, 1);
	}

	public Stepwise(double[] new_deps, double[][] new_indeps, int svdRegressionSwitch, boolean bonferroniEntry, int numThreads) {
		Xs = new Vector<double[]>();
		Ys = new Vector<String>();

		for (int i = 0; i<new_deps.length; i++) {
			Xs.add(new_indeps[i]);
			Ys.add(new_deps[i]+"");
		}

		run(svdRegressionSwitch, bonferroniEntry, numThreads);
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

		run(Integer.MAX_VALUE, false, 1);
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

	/**
	 * @param svdRegressionSwitch
	 *            when the number of indeps in the model is greater than this number, svd regression will be used
	 * @param numThreads
	 */
	public void run(int svdRegressionSwitch, boolean bonferroniEntry, int numThreads) {
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
			if(numThreads>1){
				RegressionProducer producer = new RegressionProducer(in, out, logistic, Ys, Xs, N, svdRegressionSwitch);
				WorkerTrain<RegressionModel> train = new WorkerTrain<RegressionModel>(producer, numThreads, 2, new Logger());
				int index = 0;
				int currentSize = in.size()+1;

				while (train.hasNext()) {
					RegressionModel model2 = train.next();
					if (!model2.analysisFailed() && model2.getSigs().length == currentSize + 1) {
						pvals[index] = model2.getSigs();

						if (pvals[index][pvals[index].length - 1] < lowestP) {
							lowestP = pvals[index][pvals[index].length - 1];
						}
						if (model2.getRsquare() > highestRsq) {
							highestRsq = model2.getRsquare();
							bestModel = index;
						}
					}
					index++;
				}

			} else {
			
				for (int i = 0; i < out.size(); i++) {
					in.add(out.popFirst());
					if (logistic) {
						model = new LogisticRegression(Ys, travXs(N, Xs, in));
					} else {
						model = new LeastSquares(Ys, travXs(N, Xs, in));
					}
					if (!model.analysisFailed() && model.getSigs().length == in.size() + 1) {
						pvals[i] = model.getSigs();
						if (pvals[i][pvals[i].length - 1] < lowestP) {
							lowestP = pvals[i][pvals[i].length - 1];
						}
						if (model.getRsquare() > highestRsq) {
							highestRsq = model.getRsquare();
							bestModel = i;
						}
					}
					out.add(in.popLast());
				}
			}
			if (lowestP < (bonferroniEntry ? ENTRY_PROB / M : ENTRY_PROB)) {
				in.add(out.popAt(bestModel));
				System.out.println(ext.getTime() + "\t" + in.size() + " independant variables added to the model, lowest p-value = " + lowestP);
				System.out.println(ext.getTime() + "\t" + in.size() + " independant variables added to the model, highest Rsquare = " + highestRsq);
				if (bonferroniEntry) {
					System.out.println(ext.getTime() + "bonf=" + (ENTRY_PROB / M));
				}
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
				finalModel = new LogisticRegression(Ys, travXs(N,Xs,in));
			} else {
				Vector<double[]> currentXs = travXs(N, Xs, in);
				finalModel = new LeastSquares(Ys, currentXs, false, currentXs.size() > svdRegressionSwitch);
			}
			
		}
	}
	private static class RegressionProducer implements Producer<RegressionModel>{
		private IntVector in ;
		private IntVector out;
		private boolean logistic;
		private Vector<String> Ys;
		private Vector<double[]> Xs;
		private int N;
		private int svdRegressionSwitch;

		private int index;

		public RegressionProducer(IntVector in, IntVector out, boolean logistic, Vector<String> ys, Vector<double[]> xs, int n, int svdRegressionSwitch) {
			super();
			this.in = in;
			this.out = out;
			this.logistic = logistic;
			this.Ys = ys;
			this.Xs = xs;
			this.N = n;
			this.index =0;
			this.svdRegressionSwitch = svdRegressionSwitch;
		}

		@Override
		public boolean hasNext() {
			return index < out.size();
		}

		@Override
		public Callable<RegressionModel> next() {
			in.add(out.popFirst());
			Vector<double[]> currentXs = travXs(N, Xs, in);
			RegressionWorker worker = new RegressionWorker(Ys, currentXs, logistic, svdRegressionSwitch);
			out.add(in.popLast());
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub
			
		}
	}
	
	private static class RegressionWorker implements Callable<RegressionModel>{
		private Vector<double[]> currentXs;
		private  Vector<String> Ys;
		private boolean logistic;
		private int svdRegressionSwitch;
		
		public RegressionWorker(Vector<String> ys,Vector<double[]> currentXs,  boolean logistic,int svdRegressionSwitch) {
			super();
			this.currentXs = currentXs;
			this.Ys = ys;
			this.logistic = logistic;
			this.svdRegressionSwitch = svdRegressionSwitch;
		}

		@Override
		public RegressionModel call() throws Exception {
			if (logistic) {
				return new LogisticRegression(Ys, currentXs);
			} else {
				return new LeastSquares(Ys, currentXs, false, currentXs.size() > svdRegressionSwitch);
			}
		}

	}


	public static Vector<double[]> travXs(int N, Vector<double[]> Xs, IntVector ins) {
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
		String line_ending = System.getProperty("os.name").startsWith("Windows")?"\r\n":"\n";

		RegressionModel model;
		String Rsum = " Model\t"+(logistic?" ChiSq":"    F")+"\t   Sig\t R-square"+line_ending;
		String ModelSum = ext.formStr("Variable", maxNameSize, true)+"\t   Beta\t StdErr\t      T\t    Sig"+line_ending;
		IntVector ins;
		String[] travNames;

		for (int i = 0; i<increments.size(); i++) {
			ins = increments.elementAt(i);
			travNames = new String[ins.size()];
			for (int j = 0; j<ins.size(); j++) {
				travNames[j] = varNames[ins.elementAt(j)];
			}
			if (logistic) {
				model = new LogisticRegression(Ys, travXs(N, Xs, ins));
			} else {
				model = new LeastSquares(Ys, travXs(N, Xs, ins));
			}
			model.setVarNames(travNames, maxNameSize);
			Rsum += ext.formStr(i+1+"", 4)+"\t"+(Double.isInfinite(model.getOverall())?"    .":ext.formStr(ext.formDeci(model.getOverall(), 1, true), 7))+"\t  "+ext.formDeci(model.getOverallSig(), 3, true)+"\t  "+ext.formDeci(model.getRsquare(), 3, true)+line_ending;
			ModelSum += "------ Model "+(i+1)+(i<10?" -":" ")+"---------------------------"+line_ending+model.modelSummary();
		}

		return Rsum+line_ending+ModelSum;
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

	public static class StepWiseSummary {
		private double[] sigs;
		private double[] stats;
		private int[] orderOfOriginal;

		public StepWiseSummary(double[] sigs, double[] stats, int[] orderOfOriginal) {
			super();
			this.sigs = sigs;
			this.stats = stats;
			this.orderOfOriginal = orderOfOriginal;
		}

		public double[] getSigs() {
			return sigs;
		}

		public double[] getStats() {
			return stats;
		}

		public int[] getOrderOfOriginal() {
			return orderOfOriginal;
		}

	}
	
	public StepWiseSummary getStepWiseSummary(int svdRegressionSwitch,int numThreads){
		
		if (increments.size() > 0) {
			IntVector in = new IntVector();
			IntVector out = increments.lastElement();

			double[] sigs = new double[out.size()];
			double[] stats = new double[out.size()];
			int[] orderOfOriginal = out.toArray();

			RegressionProducer producer = new RegressionProducer(in, out, logistic, Ys, Xs, N, svdRegressionSwitch);
			WorkerTrain<RegressionModel> train = new WorkerTrain<RegressionModel>(producer, numThreads, 2, new Logger());
			int index = 0;

			while (train.hasNext()) {
				RegressionModel model = train.next();
				double[] modelSigs = model.getSigs();
				sigs[index] = modelSigs[modelSigs.length - 1];
				stats[index] = model.getRsquare();
				index++;
			}
			return new StepWiseSummary(sigs, stats, orderOfOriginal);
		} else {
			return null;
		}
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
		boolean bonferroni = false;

		String usage = "\n"+
				"stats.stepwise requires 0-1 arguments\n"+
				"   (1) filename (i.e. file="+filename+" (default)\n"+
				"   (2) Bonferroni threhsold for entry instead of nominal (i.e. bonferroni="+bonferroni+" (default)\n"+
				"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bonferroni=")) {
				bonferroni = ext.parseBooleanArg(args[i]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			RegVectors rvs = procFile(filename);
			sw = new Stepwise(rvs.getDeps(), rvs.getIndeps(), bonferroni);
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
