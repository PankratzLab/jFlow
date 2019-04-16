package org.pankratzlab.common.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CountHash;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

/**
 * This class will keep a running tab of the first two statistical moments (mean and standard
 * deviation), as well as the number of outliers seen while sequentially adding numbers, without
 * actually keeping any of those numbers in memory. This is useful when it would be impractical to
 * keep or load all numbers into memory.
 */
public class Moments {

  private double first;
  private double sum;
  private double variance;
  private double lowerThrehsold;
  private double upperThrehsold;
  private double n;
  private double min;
  private double max;
  private int numLowOutliers;
  private int numHighOutliers;
  private int numNaNs;

  public Moments() {
    this(Double.NaN, Double.NaN);
  }

  /**
   * This class will keep a running tab of the first two statistical moments (mean and standard
   * deviation), as well as the number of outliers seen while sequentially adding numbers, without
   * actually keeping any of those numbers in memory. This is useful when it would be impractical to
   * keep or load all numbers into memory.
   * 
   * @param lowerThrehsold the lower threshold at which to count something an outlier
   * @param upperThrehsold the upper threshold at which to count something an outlier
   */
  public Moments(double lowerThrehsold, double upperThrehsold) {
    variance = sum = n = numLowOutliers = numHighOutliers = numNaNs = 0;
    this.min = Double.MAX_VALUE;
    this.max = Double.MIN_VALUE;
    this.lowerThrehsold = lowerThrehsold;
    this.upperThrehsold = upperThrehsold;
  }

  /**
   * Adds a new number to the mix
   * 
   * @param d the new number
   */
  public void addNum(double d) {
    if (Double.isNaN(d)) {
      numNaNs++;
      return;
    }

    if (n == 0) {
      first = d;
    } else if (n == 1) {
      variance = 2 * Math.pow((d - (first + d) / 2), 2);
    } else {
      variance = ((n - 1) * variance + (d - (sum + d) / (n + 1)) * (d - sum / n)) / (n);
    }
    if (!Double.isNaN(lowerThrehsold) && d < lowerThrehsold) {
      numLowOutliers++;
    }
    if (!Double.isNaN(upperThrehsold) && d > upperThrehsold) {
      numHighOutliers++;
    }
    if (d > max) {
      max = d;
    }
    if (d < min) {
      min = d;
    }
    sum += d;
    n++;
  }

  public double getN() {
    return n;
  }

  public double getStandard() {
    return Math.sqrt(variance);
  }

  public int getNumLowOutliers() {
    return numLowOutliers;
  }

  public int getNumHighOutliers() {
    return numHighOutliers;
  }

  public int getNumNaNs() {
    return numNaNs;
  }

  public double getMin() {
    return min;
  }

  public double getMax() {
    return max;
  }

  /**
   * Reports the number of numbers added, the mean and standard deviation and the number of outliers
   * seen (if thresholds were defined).
   */
  public String report() {
    StringBuilder sb = new StringBuilder();

    sb.append("n=" + (int) n);
    sb.append("\tmean=" + (sum / (double) n));
    sb.append("\tsd=" + Math.sqrt(variance));
    sb.append("\tmin=" + min);
    sb.append("\tmax=" + max);
    if (!Double.isNaN(lowerThrehsold)) {
      sb.append("\tnumLowOutliers=" + numLowOutliers + " ("
                + ext.prettyP((double) numLowOutliers / n) + ")");
    }
    if (!Double.isNaN(upperThrehsold)) {
      sb.append("\tnumHighOutliers=" + numHighOutliers + " ("
                + ext.prettyP((double) numHighOutliers / n) + ")");
    }
    sb.append("\tnumNaNs=" + numNaNs + " (" + ext.prettyP((double) numNaNs / n) + ")");

    return sb.toString();
  }

  /**
   * Reads in a file and keeps track of the moments for the numbers within the indices defined.
   * 
   * @param filename file with the matrix of numbers
   * @param indices the indices of the columns in the file that will be tracked
   * @param lowerThrehsold the lower threshold at which to count something an outlier
   * @param upperThrehsold the upper threshold at which to count something an outlier
   */
  private static void summarizeFile(String filename, int[] indices, double lowerThreshold,
                                    double upperThreshold) {
    BufferedReader reader;
    String temp;
    String[] line;
    Moments[] moments;
    double num;
    CountHash countHash = new CountHash();

    Logger log = new Logger(filename + ".log");

    try {
      reader = Files.getAppropriateReader(filename);
      moments = new Moments[indices.length];
      for (int i = 0; i < moments.length; i++) {
        moments[i] = new Moments(lowerThreshold, upperThreshold);
      }
      // skip header and any comments
      do {
        temp = reader.readLine();
      } while (temp.startsWith("#"));
      while ((temp = reader.readLine()) != null) {
        if (temp.startsWith("#")) continue;
        line = temp.trim().split("[\\s]+");
        for (int i = 0; i < moments.length; i++) {
          num = Double.parseDouble(line[indices[i]]);
          if (Double.isNaN(num)) {
            countHash.add(ArrayUtils.toStr(line));
          }
          moments[i].addNum(Double.parseDouble(line[indices[i]]));
        }
      }
      for (int i = 0; i < moments.length; i++) {
        log.report(moments[i].report());
      }
      String[] values = countHash.getValues();
      int[] counts = countHash.getCounts();
      if (counts.length > 0) {
        log.report("The following lines produced NaNs at their given frequency:");
      }
      for (int i = 0; i < values.length; i++) {
        log.report(values[i] + " (n=" + counts[i] + ")");
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + filename + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      log.reportException(ioe);
      return;
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "Moments.dat";
    int[] indices = new int[] {0};
    double lowerThreshold = Double.NaN;
    double upperThreshold = Double.NaN;

    StringBuilder usage = new StringBuilder();
    usage.append("\n");
    usage.append("org.genvisis.stats.Moments requires 0-1 arguments\n");
    usage.append("   (1) input filename (i.e. file=" + filename + " (default))\n");
    usage.append("   (2) indices within file to summarize (i.e. indices=0,1,5 (default is 0))\n");
    usage.append("   (3) (optional) lower bound for defining outliers (i.e. lowerThreshold=-32 (default is none))\n");
    usage.append("   (4) (optional) upper bound for defining outliers (i.e. upperThreshold=32 (default is none))\n");

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage.toString());
        System.exit(1);
      } else if (args[i].startsWith("file=")) {
        filename = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith("indices=")) {
        indices = ext.parseIntArrayArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith("lowerThreshold=")) {
        lowerThreshold = ext.parseDoubleArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith("upperThreshold=")) {
        upperThreshold = ext.parseDoubleArg(args[i]);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + args[i]);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage.toString());
      System.exit(1);
    }
    try {
      summarizeFile(filename, indices, lowerThreshold, upperThreshold);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
