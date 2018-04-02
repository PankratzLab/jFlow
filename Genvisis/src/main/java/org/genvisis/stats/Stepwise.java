package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import com.google.common.primitives.Ints;

public class Stepwise {

  public static double ENTRY_PROB = 0.05;
  // public static double ENTRY_PROB = 0.10;
  public static double REMOVAL_PROB = 0.10;
  private int M, N;
  private final List<double[]> Xs;
  private final List<String> Ys;
  private boolean logistic;
  private String[] varNames;
  private int maxNameSize;
  private Vector<IntVector> increments;
  private RegressionModel finalModel;

  public Stepwise(List<String> list, List<double[]> list2, boolean bonferroni) {
    Xs = list2;
    Ys = list;
    run(Integer.MAX_VALUE, bonferroni, 1);
  }

  public Stepwise(double[] new_deps, double[][] new_indeps, int svdRegressionSwitch,
                  boolean bonferroniEntry, int numThreads) {
    Xs = new Vector<>();
    Ys = new Vector<>();

    for (int i = 0; i < new_deps.length; i++) {
      Xs.add(new_indeps[i]);
      Ys.add(new_deps[i] + "");
    }

    run(svdRegressionSwitch, bonferroniEntry, numThreads);
  }

  public Stepwise(int[] new_deps, int[][] new_indeps) {
    Xs = new Vector<>();
    Ys = new Vector<>();
    double[] data;

    for (int i = 0; i < new_deps.length; i++) {
      data = new double[new_indeps[i].length];
      for (int j = 0; j < new_indeps[i].length; j++) {
        data[j] = new_indeps[i][j];
      }
      Xs.add(data);
      Ys.add(new_deps[i] + "");
    }

    run(Integer.MAX_VALUE, false, 1);
  }

  public static RegVectors procFile(String filename) {
    BufferedReader reader = null;
    String[] line;
    Vector<String> deps = new Vector<>();
    Vector<double[]> indeps = new Vector<>();
    double[] data;
    boolean names = false;
    Vector<String> indNames = new Vector<>();

    try {
      reader = new BufferedReader(new FileReader(filename));
      reader.mark(5000);
      line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 1; i < line.length; i++) {
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
        line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
        deps.add(line[0]);
        data = new double[line.length - 1];
        for (int i = 1; i < line.length; i++) {
          data[i - 1] = Double.valueOf(line[i]).doubleValue();
        }
        indeps.add(data);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(1);
    }

    return (names ? new RegVectors(deps, indeps, indNames) : new RegVectors(deps, indeps));
  }

  /**
   * @param svdRegressionSwitch when the number of indeps in the model is greater than this number,
   *          svd regression will be used
   * @param numThreads
   */
  public void run(int svdRegressionSwitch, boolean bonferroniEntry, int numThreads) {
    RegressionModel model;
    IntVector in = new IntVector();
    IntVector out = new IntVector();
    Vector<String> v = new Vector<>();
    double[][] pvals;
    int bestModel;
    boolean done;
    double lowestP, highestRsq;

    if (Xs.size() != Ys.size()) {
      System.err.println("Error using Vectors for stepwise regression; " + Ys.size()
                         + " dependent elements and " + Xs.size() + " independent elements");
      System.exit(11);
    }

    N = Ys.size();
    for (int i = 0; i < N && v.size() < 3; i++) {
      if (!v.contains(Ys.get(i))) {
        v.add(Ys.get(i));
      }
    }

    if (v.size() == 1) {
      System.err.println("Error - dependent variable has zero variation (" + v.elementAt(0) + ")");
      System.exit(5);
    }
    logistic = (v.size() == 2);
    System.out.println("Performing stepwise " + (logistic ? "logistic" : "linear") + " regression");

    M = Xs.get(0).length;

    varNames = new String[M];
    for (int i = 0; i < M; i++) {
      varNames[i] = "Indep " + (i + 1);
    }
    maxNameSize = M < 10 ? 8 : 7 + (M + "").length();

    for (int i = 0; i < M; i++) {
      out.add(i);
    }

    done = false;
    increments = new Vector<>();
    while (!done) {
      pvals = new double[out.size()][];
      bestModel = -1;
      highestRsq = 0;
      lowestP = 1;
      if (numThreads > 1) {
        RegressionProducer producer = new RegressionProducer(in, out, logistic, Ys, Xs, N,
                                                             svdRegressionSwitch);
        try (WorkerTrain<RegressionModel> train = new WorkerTrain<>(producer,
                                                                                   numThreads, 2,
                                                                                   new Logger())) {
          int index = 0;
          int currentSize = in.size() + 1;

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
        }

      } else {

        for (int i = 0; i < out.size(); i++) {
          in.add(out.remove(0));
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
          out.add(in.remove(in.size() - 1));
        }
      }
      if (lowestP < (bonferroniEntry ? ENTRY_PROB / M : ENTRY_PROB)) {
        in.add(out.remove(bestModel));
        System.out.println(ext.getTime() + "\t" + in.size()
                           + " independant variables added to the model, lowest p-value = "
                           + lowestP);
        System.out.println(ext.getTime() + "\t" + in.size()
                           + " independant variables added to the model, highest Rsquare = "
                           + highestRsq);
        if (bonferroniEntry) {
          System.out.println(ext.getTime() + "bonf=" + (ENTRY_PROB / M));
        }
        increments.add(in.clone());
        for (int i = in.size(); i >= 1; i--) {
          if (pvals[bestModel][i] > REMOVAL_PROB) {
            in.removeElementAt(i - 1);
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

    if (increments.size() > 0) {
      if (logistic) {
        finalModel = new LogisticRegression(Ys, travXs(N, Xs, in));
      } else {
        Vector<double[]> currentXs = travXs(N, Xs, in);
        finalModel = new LeastSquares(Ys, currentXs, false, currentXs.size() > svdRegressionSwitch);
      }

    }
  }

  private static class RegressionProducer extends AbstractProducer<RegressionModel> {

    private final IntVector in;
    private final IntVector out;
    private final boolean logistic;
    private final List<String> Ys;
    private final List<double[]> Xs;
    private final int N;
    private final int svdRegressionSwitch;

    private int index;

    public RegressionProducer(IntVector in, IntVector out, boolean logistic, List<String> ys2,
                              List<double[]> xs2, int n, int svdRegressionSwitch) {
      super();
      this.in = in;
      this.out = out;
      this.logistic = logistic;
      Ys = ys2;
      Xs = xs2;
      N = n;
      index = 0;
      this.svdRegressionSwitch = svdRegressionSwitch;
    }

    @Override
    public boolean hasNext() {
      return index < out.size();
    }

    @Override
    public Callable<RegressionModel> next() {
      in.add(out.remove(0));
      Vector<double[]> currentXs = travXs(N, Xs, in);
      RegressionWorker worker = new RegressionWorker(Ys, currentXs, logistic, svdRegressionSwitch);
      out.add(in.remove(in.size() - 1));
      index++;
      return worker;
    }
  }

  private static class RegressionWorker implements Callable<RegressionModel> {

    private final List<double[]> currentXs;
    private final List<String> Ys;
    private final boolean logistic;
    private final int svdRegressionSwitch;

    public RegressionWorker(List<String> ys2, List<double[]> currentXs, boolean logistic,
                            int svdRegressionSwitch) {
      super();
      this.currentXs = currentXs;
      Ys = ys2;
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

  public static Vector<double[]> travXs(int N, List<double[]> xs2, IntVector ins) {
    Vector<double[]> v = new Vector<>();
    double[] data, newdata;

    for (int i = 0; i < N; i++) {
      newdata = new double[ins.size()];
      data = xs2.get(i);
      for (int j = 0; j < ins.size(); j++) {
        newdata[j] = data[ins.elementAt(j)];
      }
      v.add(newdata);
    }

    return v;
  }

  public void setVarNames(List<String> vNames) {
    String[] names = new String[vNames.size()];
    for (int i = 0; i < names.length; i++) {
      names[i] = vNames.get(i);
    }
    setVarNames(names);
  }

  public void setVarNames(String[] names) {
    if (names.length != M) {
      System.err.println("Error naming independent variables: " + M + " variables, and only "
                         + names.length + " names");
      return;
    }
    varNames = names;
    maxNameSize = 8;
    for (int i = 0; i < M; i++) {
      if (names[i].length() > maxNameSize) {
        maxNameSize = names[i].length();
      }
    }
  }

  public String getSummary() {
    String line_ending = System.getProperty("os.name").startsWith("Windows") ? "\r\n" : "\n";

    RegressionModel model;
    String Rsum = " Model\t" + (logistic ? " ChiSq" : "    F") + "\t   Sig\t R-square"
                  + line_ending;
    String ModelSum = ext.formStr("Variable", maxNameSize, true)
                      + "\t   Beta\t StdErr\t      T\t    Sig" + line_ending;
    IntVector ins;
    String[] travNames;

    for (int i = 0; i < increments.size(); i++) {
      ins = increments.elementAt(i);
      travNames = new String[ins.size()];
      for (int j = 0; j < ins.size(); j++) {
        travNames[j] = varNames[ins.elementAt(j)];
      }
      if (logistic) {
        model = new LogisticRegression(Ys, travXs(N, Xs, ins));
      } else {
        model = new LeastSquares(Ys, travXs(N, Xs, ins));
      }
      model.setVarNames(travNames, maxNameSize);
      Rsum += ext.formStr(i + 1 + "", 4) + "\t"
              + (!Double.isFinite(model.getOverall()) ? "    ."
                                                      : ext.formStr(ext.formDeci(model.getOverall(),
                                                                                 1, true),
                                                                    7))
              + "\t  " + ext.formDeci(model.getOverallSig(), 3, true) + "\t  "
              + ext.formDeci(model.getRsquare(), 3, true) + line_ending;
      ModelSum += "------ Model " + (i + 1) + (i < 10 ? " -" : " ") + "---------------------------"
                  + line_ending + model.modelSummary();
    }

    return Rsum + line_ending + ModelSum;
  }

  public String getAccuracySummary() {
    return increments.size() == 0 ? "no variables in final model"
                                  : (logistic ? ((LogisticRegression) finalModel).getAccuracySummary()
                                              : "");
  }

  public double[] getFinalPredicteds() {
    return increments.size() == 0 ? new double[N] : finalModel.getPredicteds();
  }

  public double[] getFinalResiduals() {
    return increments.size() == 0 ? new double[N] : finalModel.getResiduals();
  }

  public static class StepWiseSummary {

    private final double[] sigs;
    private final double[] stats;
    private final int[] orderOfOriginal;

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

  public StepWiseSummary getStepWiseSummary(int svdRegressionSwitch, int numThreads) {

    if (increments.size() > 0) {
      IntVector in = new IntVector();
      IntVector out = increments.lastElement();

      double[] sigs = new double[out.size()];
      double[] stats = new double[out.size()];
      int[] orderOfOriginal = Ints.toArray(out);

      RegressionProducer producer = new RegressionProducer(in, out, logistic, Ys, Xs, N,
                                                           svdRegressionSwitch);
      try (WorkerTrain<RegressionModel> train = new WorkerTrain<>(producer,
                                                                                 numThreads, 2,
                                                                                 new Logger())) {
        int index = 0;

        while (train.hasNext()) {
          RegressionModel model = train.next();
          double[] modelSigs = model.getSigs();
          sigs[index] = modelSigs[modelSigs.length - 1];
          stats[index] = model.getRsquare();
          index++;
        }
      }
      return new StepWiseSummary(sigs, stats, orderOfOriginal);
    } else {
      return null;
    }
  }

  public String getFinalNames() {
    IntVector ins;
    String finalNames = "";

    if (increments.size() == 0) {
      return "no variables in final model";
    }

    ins = increments.lastElement();
    for (int j = 0; j < ins.size(); j++) {
      finalNames += (j == 0 ? "" : "\t") + varNames[ins.elementAt(j)];
    }

    return finalNames;
  }

  public void dumpData(String filename) {
    double[] data;

    try {
      PrintWriter writer = Files.openAppropriateWriter(filename);
      writer.print("Dep");
      for (String varName : varNames) {
        writer.print("\t" + varName);
      }
      writer.println();
      for (int i = 0; i < N; i++) {
        writer.print(Ys.get(i));
        data = Xs.get(i);
        for (int j = 0; j < M; j++) {
          writer.print("\t" + data[j]);
        }
        writer.println();
      }
      writer.close();
    } catch (IOException ioe) {
      System.err.println("Error writing dump file: " + filename);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "vars.txt";
    Stepwise sw;
    boolean bonferroni = false;

    String usage = "\n" + "stats.stepwise requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default)\n"
                   + "   (2) Bonferroni threhsold for entry instead of nominal (i.e. bonferroni="
                   + bonferroni + " (default)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("bonferroni=")) {
        bonferroni = ext.parseBooleanArg(arg);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      RegVectors rvs = procFile(filename);
      sw = new Stepwise(rvs.getDeps(), rvs.getIndeps(), bonferroni);
      if (rvs.getVarNames() != null) {
        sw.setVarNames(rvs.getVarNames());
      }
      System.out.println(sw.getSummary());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}

class RegVectors {

  private final List<String> deps;

  private final List<double[]> indeps;

  private final List<String> varNames;

  public RegVectors(List<String> deps, List<double[]> indeps) {
    this.deps = deps;
    this.indeps = indeps;
    varNames = null;
  }

  public RegVectors(List<String> deps, List<double[]> indeps, List<String> varNames) {
    this.deps = deps;
    this.indeps = indeps;
    this.varNames = varNames;
  }

  public List<String> getDeps() {
    return deps;
  }

  public List<double[]> getIndeps() {
    return indeps;
  }

  public List<String> getVarNames() {
    return varNames;
  }
}
