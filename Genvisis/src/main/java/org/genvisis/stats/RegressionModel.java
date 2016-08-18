package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public abstract class RegressionModel {
  protected boolean verbose;
  protected Logger log;
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
  protected int sigfigs;

  public RegressionModel() {
    verbose = true;
    log = new Logger();
    sigfigs = 3;
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
    for (double beta : betas) {
      if (Double.isNaN(beta)) {
        return true;
      }
    }
    return analysisFailed;
  }

  protected void fail() {
    betas = new double[M + 1];
    SEofBs = new double[M + 1];
    stats = new double[M + 1];
    sigs = new double[M + 1];
    predicteds = new double[N];
    residuals = new double[N];
    for (int i = 0; i < sigs.length; i++) {
      betas[i] = SEofBs[i] = stats[i] = sigs[i] = Double.NaN;
    }
    for (int i = 0; i < N; i++) {
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
    if (names.length != M) {
      System.err.println("Error naming independent variables: " + M + " variables, and "
                         + names.length + " names");
      return;
    }
    varNames = new String[M + 1];
    varNames[0] = "Constant";
    maxNameSize = 8;
    for (int i = 0; i < M; i++) {
      varNames[i + 1] = names[i];
      if (names[i].length() > maxNameSize) {
        maxNameSize = names[i].length();
      }
    }
  }

  public void setSigFigs(int sigfigs) {
    this.sigfigs = sigfigs;
  }

  public void dumpData(String filename) {
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(filename));
      writer.print("Dep");
      for (int i = 1; i < varNames.length; i++) {
        writer.print("\t" + varNames[i]);
      }
      writer.println(famIDs == null ? "" : "\tFamID");
      for (int i = 0; i < N; i++) {
        writer.print(ext.formDeci(deps[i], 5));
        for (int j = 0; j < M; j++) {
          writer.print("\t" + indeps[i][j]);
        }
        writer.println(famIDs == null ? "" : "\t" + famIDs[i]);
      }
      writer.close();
    } catch (IOException ioe) {
      log.reportError("Error writing dump file: " + filename);
    }
  }

  public void onePerFamily(String[] famIDs, int numReps, int bootReps) {
    double[][] statPerms = new double[stats.length][numReps];
    double[][] betaPerms = new double[betas.length][numReps];
    onePer = true;
    String[] fams = Array.unique(famIDs);
    Hashtable<String, IntVector> hash = new Hashtable<String, IntVector>();
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

    for (String fam : fams) {
      hash.put(fam, new IntVector());
    }

    overall = Double.NaN;
    overallSig = Double.NaN;
    Rsquare = Double.NaN;
    for (int i = 0; i < famIDs.length; i++) {
      predicteds[i] = Double.NaN;
      residuals[i] = Double.NaN;
      hash.get(famIDs[i]).add(i);
    }

    if (verbose) {
      log.reportError("Rep 0", false, true);
    }
    failures = 0;
    for (int i = 0; i < numPermutations; i++) {
      if (verbose && (i + 1) % (numReps / 10) == 0) {
        for (int j = 0; j < progress.length(); j++) {
          log.reportError("\b", false, true);
        }
        progress = (i + 1) + (failures > 0 ? ", " + failures + " failures" : "");
        log.report(progress);
      }

      do {
        for (int j = 0; j < fams.length; j++) {
          iv = hash.get(fams[j]);
          index = iv.elementAt((int) (Math.random() * iv.size()));
          rDeps[j] = deps[index];
          rIndeps[j] = indeps[index];
        }

        if (logistic) {
          model = new LogisticRegression(rDeps, rIndeps, true, true);
          logCounts[0] += ((LogisticRegression) model).getY0();
          logCounts[1] += ((LogisticRegression) model).getY1();
        } else {
          model = new LeastSquares(rDeps, rIndeps);
        }
        rBetas = model.getBetas();
        rStats = model.getStats();
        failures += model.analysisFailed() ? 1 : 0;
        for (int j = 1; j <= M; j++) {
          betaPerms[j][i] = rBetas[j];
          statPerms[j][i] = rStats[j];
        }
      } while (model.analysisFailed());
    }
    if (verbose) {
      log.report("");
    }

    for (int i = 1; i <= M; i++) {
      betas[i] = Array.bootstrap(betaPerms[i], numBootReps, verbose)[0];
      stats[i] = Array.bootstrap(statPerms[i], numBootReps, verbose)[0];
      sigs[i] = logistic ? ProbDist.ChiDist(stats[i], 1) : ProbDist.TDist(stats[i], N - M - 1);
      SEofBs[i] = Double.NaN;
    }
    if (logistic) {
      ((LogisticRegression) this).computeOddsRatios();
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

    for (int i = 0; i < use.length; i++) {
      use[i] = true;
      if ((deps[i] + "").equals("NaN")) {
        countNoDeps++;
        use[i] = false;
      }
      for (int j = 0; j < indeps[i].length; j++) {
        if (use[i] && (indeps[i][j] + "").equals("NaN")) {
          countNoIndeps++;
          use[i] = false;
        }
      }
    }

    count = Array.booleanArraySum(use);
    newDeps = new double[count];
    newIndeps = new double[count][];
    count = 0;
    for (int i = 0; i < use.length; i++) {
      if (use[i]) {
        newDeps[count] = deps[i];
        newIndeps[count] = indeps[i];
        count++;
      }
    }

    if (verbose) {
      if (countNoDeps > 0) {
        log.reportError(countNoDeps
                        + " individuals were dropped because they lacked the dependent variable");
      }
      if (countNoIndeps > 0) {
        log.reportError(countNoIndeps
                        + " individuals were dropped because they lacked at least one independent variable");
      }
    }

    deps = newDeps;
    indeps = newIndeps;
  }

  public boolean dataQC() {
    Hashtable<String, String> hash = new Hashtable<String, String>();
    double[][] oldIndeps;
    double diff, offset;
    String[] keys;
    Vector<String> newVariableNames;

    if (deps == null) {
      log.reportError("Array with dependent variables is null");
      return false;
    }

    if (indeps == null) {
      log.reportError("Array with independent variables is null");
      return false;
    }

    for (int i = 0; i < deps.length && hash.size() < 3; i++) {
      if (!hash.containsKey(deps[i] + "")) {
        hash.put(deps[i] + "", "");
      }
    }
    if (hash.size() < 2) {
      if (verbose) {
        log.reportError("No variance in the dependent variable");
      }
      return false;
    }
    if (logistic) {
      if (hash.size() > 2) {
        log.reportError("Error in logistic regression - Dependent variables must be 2 and no more than 2 consecutive integers");
        return false;
      }
      keys = HashVec.getKeys(hash);
      diff = Double.parseDouble(keys[0]) - Double.parseDouble(keys[1]);
      if (Math.abs(diff) != 1) {
        log.reportError("Error in logistic regression - Dependent variables must be 2 and no more than 2 consecutive integers");
        return false;
      }

      offset = diff < 0 ? Double.parseDouble(keys[0]) : Double.parseDouble(keys[1]);
      if (offset != 0) {
        for (int i = 0; i < deps.length; i++) {
          deps[i] = deps[i] - offset;
        }
      }
    }

    newVariableNames = new Vector<String>();
    newVariableNames.add(varNames[0]);
    for (int j = 0; j < indeps[0].length; j++) {
      hash.clear();
      for (int i = 0; i < indeps.length && hash.size() < 2; i++) {
        if (!hash.containsKey(indeps[i][j] + "")) {
          hash.put(indeps[i][j] + "", "");
        }
      }
      if (hash.size() < 2) {
        if (verbose) {
          log.reportError("No variance in independent variable number " + (j + 1)
                          + "; collapsing and ignoring");
        }
        oldIndeps = indeps.clone();
        indeps = new double[oldIndeps.length][oldIndeps[0].length - 1];
        for (int k = 0; k < j; k++) {
          for (int i = 0; i < oldIndeps.length; i++) {
            indeps[i][k] = oldIndeps[i][k];
          }
        }
        for (int k = j + 1; k < oldIndeps[0].length; k++) {
          for (int i = 0; i < oldIndeps.length; i++) {
            indeps[i][k - 1] = oldIndeps[i][k];
          }
        }
        varNames = Array.removeFromArray(varNames, j + 1);
        j--;
      } else {
        // newVariableNames.add(varNames[j+1]);
      }
    }
    // varNames = Array.toStringArray(newVariableNames);

    return indeps[0].length > 0 ? true : false;
  }

  public static double[] processDeps(Vector<String> vDeps) {
    double[] deps = new double[vDeps.size()];

    for (int i = 0; i < deps.length; i++) {
      deps[i] = Double.parseDouble(vDeps.elementAt(i));
    }

    return deps;
  }

  // TODO These may eventually replace the others, but only after getting rid of the ambiguous
  // erasures
  // public static double[] processStringDeps(Vector<String> vDeps) {
  // double[] deps = new double[vDeps.size()];
  //
  // for (int i = 0; i<deps.length; i++) {
  // deps[i] = Double.parseDouble(vDeps.elementAt(i));
  // }
  //
  // return deps;
  // }
  //
  // public static double[][] processIntIndeps(Vector<int[]> vIndeps) {
  // double[][] indeps = new double[vIndeps.size()][];
  //
  // for (int i = 0; i<indeps.length; i++) {
  // indeps[i] = Array.toDoubleArray(vIndeps.elementAt(i));
  // }
  //
  // return indeps;
  // }
  //
  // public static double[][] processDoubleIndeps(Vector<double[]> vIndeps) {
  // double[][] indeps = new double[vIndeps.size()][];
  //
  // for (int i = 0; i<indeps.length; i++) {
  // indeps[i] = vIndeps.elementAt(i);
  // }
  //
  // return indeps;

  @SuppressWarnings({"rawtypes"})
  public static double[][] processIndeps(Vector vIndeps) {
    double[][] indeps = new double[vIndeps.size()][];
    int[] intarray = {0};

    for (int i = 0; i < indeps.length; i++) {
      if (vIndeps.elementAt(i).getClass() == intarray.getClass()) {
        indeps[i] = Array.toDoubleArray((int[]) vIndeps.elementAt(i));
      } else {
        indeps[i] = (double[]) vIndeps.elementAt(i);
      }
    }

    return indeps;
  }

  public static RegressionModel determineAppropriate(double[] deps, double[] indeps,
                                                     boolean bypassDataCheck, boolean verbose) {
    double[][] newIndeps = new double[indeps.length][1];

    for (int i = 0; i < indeps.length; i++) {
      newIndeps[i][0] = indeps[i];
    }

    return determineAppropriate(deps, newIndeps, bypassDataCheck, verbose);
  }

  public static RegressionModel determineAppropriate(double[] deps, double[][] indeps,
                                                     boolean bypassDataCheck, boolean verbose) {
    Vector<String> depCount = new Vector<String>();

    for (int i = 0; i < deps.length && depCount.size() <= 2; i++) {
      if (!(deps[i] + "").equals("NaN")) {
        HashVec.addIfAbsent(deps[i] + "", depCount);
      }
    }

    return depCount.size() == 2 ? new LogisticRegression(deps, indeps, bypassDataCheck, verbose)
                                : new LeastSquares(deps, indeps, bypassDataCheck, verbose);
  }

  public static boolean isBinaryTrait(String[] pheno, Logger log) {
    Hashtable<String, String> hash = new Hashtable<String, String>();
    String[] keys;
    double diff;

    hash = new Hashtable<String, String>();
    if (pheno == null) {
      log.reportError("Error - null array");
      return false;
    }

    for (int i = 0; i < pheno.length && hash.size() < 4; i++) {
      if (!hash.containsKey(pheno[i] + "") && !ext.isMissingValue(pheno[i])) {
        hash.put(pheno[i] + "", "");
      }
    }

    if (hash.size() < 2) {
      log.reportError("No variance in the dependent variable");
      return false;
    }

    if (hash.size() == 3 && hash.containsKey("0")) {
      log.report("Warning - phenotype would be binary, if it weren't for the zero values; not flagged for logistic");
    }

    if (hash.size() == 2) {
      keys = HashVec.getKeys(hash);
      diff = Double.parseDouble(keys[0]) - Double.parseDouble(keys[1]);
      if (Math.abs(diff) != 1) {
        log.reportError("Error - only two values, but they are not consecutive integers; not flagged for logistic");
        return false;
      }
    }

    return hash.size() == 2;
  }

  public static String[] getIDsWithCompleteData(String filename, boolean alsoUseFamID, Logger log) {
    BufferedReader reader;
    String[] line;
    Vector<String> v;
    String delimiter;
    int numElements;
    boolean use;

    v = new Vector<String>();
    try {
      reader = Files.getAppropriateReader(filename);
      delimiter = Files.determineDelimiter(filename, log);
      line = reader.readLine().trim().split(delimiter);
      numElements = line.length;
      if (ext.indexOfStr(line[0], ext.COMMON_IDS) == -1) {
        log.report("Warning - unexpected id for phenotype file " + filename + ": " + line[0]);
      }
      while (reader.ready()) {
        line = reader.readLine().split(delimiter);
        if (line.length != numElements) {
          System.err.println("Error - mismatched number of elements for ID '" + line[0]
                             + "' (expecting " + numElements + ", found " + line.length
                             + "); check delimiter or for trailing whitespace");
        }
        use = true;
        for (int i = 1; i < line.length; i++) {
          if (ext.isMissingValue(line[i]) || !ext.isValidDouble(line[i])) {
            use = false;
          }
        }
        if (use) {
          v.add(line[0]);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    return Array.toStringArray(v);
  }

  public static boolean[] getRowsWithCompleteData(double[] deps, double[][] indeps, Logger log) {
    boolean[] use;

    if (deps == null || indeps == null || deps.length != indeps.length) {
      log.reportError("Error - cannot determine rows with compelte data since the deps length and the indeps length are not equal");
      return new boolean[] {};
    }

    use = Array.booleanArray(deps.length, true);
    for (int i = 0; i < deps.length; i++) {
      if (Double.isNaN(deps[i])) {
        use[i] = false;
      }
      for (int j = 0; j < indeps[i].length; j++) {
        if (Double.isNaN(indeps[i][j])) {
          use[i] = false;
        }
      }
    }

    return use;
  }

  public static boolean[] getRowsWithCompleteData(String[] deps, String[][] indeps, Logger log) {
    return getRowsWithCompleteData(Array.toDoubleArray(deps, true),
                                   Array.toDoubleArrays(indeps, true), log);
  }
}
