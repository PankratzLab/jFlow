package org.genvisis.mining;

import com.google.common.primitives.Doubles;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class BoxCox {
  public static final int LOG_LIKELIHOOD = 1;
  public static final int KURTOSIS = 2;
  public static final double INIT_LOW = -10;
  public static final double INIT_HIGH = 10;
  public static final int NUM_STEPS = 10;
  public static final int SIGFIGS = 5;
  public static final double KURT_CI = 0.02;

  private final double[] data;
  private final double maxLL_Lambda;
  private final double minKurt_Lambda;
  private final double N;
  private double sumLN;
  private double shift;
  private final Logger log;

  public BoxCox(double[] datas, Logger log) {
    double low, high, convergence_limit;
    double min;

    data = datas.clone();
    this.log = log;
    N = data.length;

    shift = 0;
    min = Array.min(data);
    if (min < 0) {
      log.reportError("Warning - the original data was not suitable for transformation -- the minimum ("
                      + min + ") cannot be less than 1.0; therefore added " + (Math.abs(min) + 1)
                      + " to every element");
      shift = Math.abs(min) + 1;
      for (int i = 0; i < N; i++) {
        data[i] += shift;
      }
    }

    // if (Array.min(data)<1) {
    // log.reportError("Minimum is less than 1, adding 1.0 to all values");
    // for (int i = 0; i<N; i++) {
    // data[i]++;
    // }
    // shift = 1;
    // }

    sumLN = 0;
    for (int i = 0; i < N; i++) {
      sumLN += Math.log(data[i]);
    }

    low = INIT_LOW;
    high = INIT_HIGH;
    convergence_limit = Math.pow(0.1, SIGFIGS + 1);

    maxLL_Lambda = findLambda(low, high, convergence_limit, LOG_LIKELIHOOD);
    minKurt_Lambda =
        findLambda(maxLL_Lambda - KURT_CI, maxLL_Lambda + KURT_CI, convergence_limit, KURTOSIS);
  }

  private double findLambda(double low, double high, double convergence_limit, int method) {
    double step, value, maxVal;
    int maxPos;
    double lam;

    while (high - low > convergence_limit) {
      // System.out.println(ext.formDeci(low, SIGFIGS) +"\t"+
      // ext.formDeci(high, SIGFIGS));
      maxPos = -1;
      switch (method) {
        case LOG_LIKELIHOOD:
          maxVal = Double.NEGATIVE_INFINITY;
          break;
        case KURTOSIS:
          maxVal = Double.POSITIVE_INFINITY;
          break;
        default:
          maxVal = 0;
          log.reportError("Error - '" + method + "' is an invalid evaluation method");
          System.exit(1);
      }
      step = (high - low) / NUM_STEPS;
      for (int i = 0; i < NUM_STEPS; i++) {
        lam = low + i * step;
        value = evaluate(data, lam, method);
        switch (method) {
          case LOG_LIKELIHOOD:
            if (value > maxVal) {
              maxPos = i;
              maxVal = value;
            }
            break;
          case KURTOSIS:
            if (Math.abs(value) < Math.abs(maxVal)) {
              maxPos = i;
              maxVal = value;
            }
        }

      }
      if (maxPos == -1) {
        log.reportError("Error - yo, you got serious problems of non-convergence");
      } else if (maxPos == 0) {
        high = low + step;
        low = low - NUM_STEPS * step * 2;
      } else if (maxPos == NUM_STEPS - 1) {
        low = high - step;
        high = high + NUM_STEPS * step * 2;
      } else {
        low = low + (maxPos - 1) * step;
        high = high - (10 - maxPos - 1) * step;
      }
    }

    return Double.parseDouble(ext.formDeci((high + low) / 2, SIGFIGS));

  }

  private double evaluate(double[] data, double lambda, int method) {
    switch (method) {
      case LOG_LIKELIHOOD:
        return -1 * (N - 1) / 2 * Math.log(Array.variance(transform(data, lambda)))
               + (lambda - 1) * (N - 1) / N * sumLN;
      case KURTOSIS:
        return Array.kurtosis(transform(data, lambda));
    }
    return -1;
  }

  public boolean isShifted() {
    return shift == 1;
  }

  public double getLambda_MaxLL() {
    return maxLL_Lambda;
  }

  public double[] getTransform_MaxLL() {
    return transform(data, getLambda_MaxLL());
  }

  public double getLambda_MinKurt() {
    return minKurt_Lambda;
  }

  public double[] getTransform_MinKurt() {
    return transform(data, getLambda_MinKurt());
  }

  public static double[] transform(double[] data, double lambda) {
    double[] tr = new double[data.length];

    if (lambda == 0) {
      for (int i = 0; i < data.length; i++) {
        tr[i] = Math.log(data[i]);
      }
    } else {
      for (int i = 0; i < data.length; i++) {
        tr[i] = (Math.pow(data[i], lambda) - 1) / lambda;
      }
    }

    return tr;
  }

  // public double lookUpValue_MaxLL(double value) {
  // if (maxLL_Lambda==0) {
  // return Math.log(value+shift);
  // } else {
  // return (Math.pow(value+shift, maxLL_Lambda)-1)/maxLL_Lambda;
  // }
  // }
  //
  // public double lookUpValue_MinKurt(double value) {
  // if (minKurt_Lambda==0) {
  // return Math.log(value+shift);
  // } else {
  // return (Math.pow(value+shift, minKurt_Lambda)-1)/minKurt_Lambda;
  // }
  // }

  public static void procFile(String filename) throws IOException {
    BufferedReader reader;
    DoubleVector dv = new DoubleVector();
    String temp;
    BoxCox bc;
    double[] data;

    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        temp = reader.readLine();
        try {
          dv.add(Double.parseDouble(temp));
        } catch (NumberFormatException nfe) {
          System.err.println("Error - '" + temp + "' is not a valid number");
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

    data = Doubles.toArray(dv);
    bc = new BoxCox(data, new Logger());
    System.out.println("Maximizing Log-Likelihood...");
    System.out.println("Optimal lambda: " + bc.getLambda_MaxLL());
    System.out.println("Kurtosis: " + Array.kurtosis(bc.getTransform_MaxLL()));
    System.out.println("");
    System.out.println("Minimizing kurtosis...");
    System.out.println("Optimal lambda: " + bc.getLambda_MinKurt());
    System.out.println("Kurtosis: " + Array.kurtosis(bc.getTransform_MinKurt()));
    System.out.println("");

  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "BoxCox.dat";

    String usage = "\n" + "park.BoxCox requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      procFile(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
