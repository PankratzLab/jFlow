// k-nearest neighbors algorithm
package org.genvisis.mining;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

public class kNN {
  public static final int[] PARTITION_DEFAULTS = {60, 40, 0};
  public static final boolean NORMALIZE_DEFAULT = true;
  public static final int MAX_K = 20;
  public static final double BINARY_CUTOFF_DEFAULT = 0.5;
  public static final double TARGET_TARGET_DEFAULT = 1.0;
  public static final boolean WEIGHTED_DEFAULT = true;
  private double[] targets;
  private double[][] predictors;
  private int[][] partitions;
  private int bestK;
  private double[][] errorRates;
  private boolean quant;
  private int startK;
  private int stopK;
  private String output;
  private boolean normalize;
  private double cutoff;
  private double target_target;
  private final boolean weighted;
  private boolean run;

  public kNN() {
    startK = -1;
    stopK = -1;
    output = "kNN-report.out";
    normalize = NORMALIZE_DEFAULT;
    cutoff = BINARY_CUTOFF_DEFAULT;
    target_target = TARGET_TARGET_DEFAULT;
    weighted = WEIGHTED_DEFAULT;
    run = false;
  }

  public kNN(double[] newTargets, double[][] newPredictors,
             int[][] newPartitions) throws IOException {
    this();
    targets = newTargets;
    predictors = newPredictors;
    partitions = newPartitions;
  }

  public kNN(double[] newTargets, double[][] newPredictors) throws IOException {
    this();
    targets = newTargets;
    predictors = newPredictors;
    partitions = createPartitions(targets.length);
  }

  public kNN calc() {
    if (startK == -1) {
      startK = 1;
    }
    if (stopK == -1) {
      stopK = partitions[0].length < MAX_K ? partitions[0].length : MAX_K;
    }
    if (normalize) {
      predictors = Transformations.transform(predictors, Transformations.NORMALIZE);
    }
    if (partitions[1].length == 0) {
      bestK = -1;
    } else {
      errorRates = new double[stopK + 1][3];
      double minError = Double.POSITIVE_INFINITY;
      Vector<String> v = new Vector<String>();
      int count = 0;

      while (v.size() < 3 && count < targets.length) {
        HashVec.addIfAbsent(targets[count++] + "", v);
      }
      if (v.size() > 2) {
        quant = true;
      } else {
        for (int i = 0; i < targets.length; i++) {
          targets[i] = targets[i] == target_target ? 1 : 0;
        }
      }

      for (int i = startK; i <= stopK; i++) {
        for (int j = 0; j < 3; j++) {
          errorRates[i][j] = calcError(targetsFromPartition(j),
                                       score(predictorsFromPartition(j), i, false));
        }
        if (errorRates[i][1] < minError) {
          minError = errorRates[i][1];
          bestK = i;
        }
      }
    }

    run = true;
    return this;
  }

  public double[] targetsFromPartition(int part) {
    double[] deps = new double[partitions[part].length];

    for (int i = 0; i < partitions[part].length; i++) {
      deps[i] = targets[partitions[part][i]];
    }

    return deps;
  }

  public double[][] predictorsFromPartition(int part) {
    double[][] indeps = new double[partitions[part].length][];

    for (int i = 0; i < partitions[part].length; i++) {
      indeps[i] = predictors[partitions[part][i]];
    }

    return indeps;
  }

  public double calcError(double[] actuals, double[] predicteds) {
    double err = 0;

    if (quant) {
      err = PredictionEvaluation.calculateError(actuals, predicteds);
    } else {
      for (int i = 0; i < predicteds.length; i++) {
        if (actuals[i] != predicteds[i]) {
          err += 1;
        }
      }
      err /= predicteds.length;
    }

    return err;
  }

  public void setNormalize(boolean aFlag) {
    normalize = aFlag;
  }

  public void setCutoff(double value) {
    cutoff = value;
  }

  public void setTarget(double value) {
    target_target = value;
  }

  public void setPartitions(int[] newParts) {
    partitions = createPartitions(targets.length, newParts);
  }

  public void setMaxK(int value) {
    stopK = value;
  }

  public void setK(int value) {
    startK = stopK = value;
  }

  public double[] score(double[][] indeps, int K, boolean raw) {
    double[] predicteds = new double[indeps.length];
    int[] keys;
    double[] distances;
    double sum;
    // double[] weights = new double[K];

    for (int i = 0; i < indeps.length; i++) {
      distances = new double[partitions[0].length];
      for (int j = 0; j < partitions[0].length; j++) {
        for (int k = 0; k < predictors[0].length; k++) {
          distances[j] += Math.pow(indeps[i][k] - predictors[partitions[0][j]][k], 2);
        }
        distances[j] = Math.sqrt(distances[j]);
      }
      keys = Sort.quicksort(distances);

      sum = 0;
      for (int j = 0; j < K; j++) {
        distances[keys[j]] = distances[keys[j]] == 0 ? 0.00000001 : distances[keys[j]];
        sum += 1 / distances[keys[j]];
      }

      for (int j = 0; j < K; j++) {
        predicteds[i] += (quant && weighted ? (1 / distances[keys[j]]) / sum : 1 / (double) K)
                         * targets[partitions[0][keys[j]]];
      }

      if (!raw && !quant) {
        predicteds[i] = predicteds[i] >= cutoff ? 1 : 0;
      }
    }

    return predicteds;
  }

  public void report() {
    PrintWriter writer;
    double[] predicteds, actuals;
    int[][] counts;

    if (!run) {
      calc();
    }

    try {
      writer = new PrintWriter(new FileWriter(output));
      writer.println("k-nearest neighbors report:");
      writer.println("k\tT Error\tV Error" + (partitions[2].length == 0 ? "" : "\tTest Error"));
      for (int i = startK; i <= stopK; i++) {
        writer.println(i + "\t" + ext.formDeci(errorRates[i][0], 4, true) + "\t"
                       + ext.formDeci(errorRates[i][1], 4, true)
                       + (partitions[2].length == 0 ? ""
                                                    : "\t"
                                                      + ext.formDeci(errorRates[i][2], 4, true))
                       + (i == bestK ? "\t<-- Best K" : ""));
      }
      writer.println();
      writer.println();
      if (!quant) {
        for (int set = 0; set < 3; set++) {
          if (partitions[set].length > 0) {
            actuals = targetsFromPartition(set);
            predicteds = score(predictorsFromPartition(set), bestK, false);
            counts = new int[2][2];
            for (int i = 0; i < actuals.length; i++) {
              counts[actuals[i] == 1.0 ? 0 : 1][predicteds[i] == 1.0 ? 0 : 1]++;
            }

            writer.println("Cut off Prob.Val. for Success (Updatable)			" + cutoff);
            writer.println((set == 0 ? "Training" : (set == 0 ? "Validation" : "Test"))
                           + " Data scoring - Summary Report (for k=" + bestK + ")");
            writer.println("\nClassification Confusion Matrix");
            writer.println("\tPredicted Class");
            writer.println("Actual\t1\t0");
            writer.println("1\t\t" + counts[0][0] + "\t" + counts[0][1]);
            writer.println("0\t\t" + counts[1][0] + "\t" + counts[1][1]);
            writer.println();
            writer.println("\nError Report");
            writer.println("Class\t# Cases\t# Errors\t% Error");
            writer.println("1\t" + (counts[0][0] + counts[0][1]) + "\t\t" + counts[0][1] + "\t\t"
                           + ext.formDeci(((double) counts[0][1]
                                           / (double) (counts[0][0] + counts[0][1]))
                                          * 100, 2, true)
                           + "%");
            writer.println("0\t" + (counts[1][0] + counts[1][1]) + "\t\t" + counts[1][0] + "\t\t"
                           + ext.formDeci(((double) counts[1][0]
                                           / (double) (counts[1][0] + counts[1][1]))
                                          * 100, 2, true)
                           + "%");
            writer.println();

          }
        }
      }

      actuals = targetsFromPartition(1);
      predicteds = score(predictorsFromPartition(1), bestK, true);
      writer.println("Classification of Validation Data (for k=" + bestK + ")");
      writer.println(quant ? "Actual\tPredicted" : "Predicted\tActual\tprobability");
      for (int i = 0; i < predicteds.length; i++) {
        writer.println((quant ? "" : (predicteds[i] >= cutoff ? 1 : 0) + "\t") + actuals[i] + "\t"
                       + predicteds[i]);
      }

      writer.close();
    } catch (IOException ioe) {
      System.err.println("Error writing kNN report");
    }
  }

  public kNN(String filename) {
    this();

    BufferedReader reader;
    String[] header, line;
    boolean manualPartition;
    Vector<double[]> records = new Vector<double[]>();
    double[] dataline;
    DoubleVector dependents = new DoubleVector();
    IntVector[] partitionIndices = new IntVector[] {new IntVector(), new IntVector(),
                                                    new IntVector()};
    int count = 0;

    try {
      reader = new BufferedReader(new FileReader(filename));
      output = (filename.indexOf(".") >= 0 ? filename.substring(0, filename.lastIndexOf("."))
                                           : filename)
               + "-report.out";
      header = reader.readLine().split("\t", -1);
      System.out.println("Assuming the last column (" + header[header.length - 1]
                         + ") is the dependent variable");
      if (manualPartition = header[0].equals("Partition")) {
        System.out.println("Assuming the first column (" + header[0]
                           + ") will be used to partition the data");
      } else {
        System.out.println("Assuming there is no provided 'Partition' column; using defaults ("
                           + PARTITION_DEFAULTS[0] + ":" + PARTITION_DEFAULTS[1] + ":"
                           + PARTITION_DEFAULTS[2] + ")");
      }
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        if (line.length != header.length) {
          System.err.println("Error - number of columns for row " + (count + 1)
                             + " does not match that for header");
          System.exit(3);
        }
        dataline = new double[header.length - (manualPartition ? 2 : 1)];

        for (int i = 0; i < line.length; i++) {
          try {
            if (i == line.length - 1) {
              dependents.add(Double.parseDouble(line[i]));
            } else if (manualPartition && i == 0) {
              if (line[0].length() != 1
                  || (!line[0].equals("T") && !line[0].equals("V") && !line[0].equals("S"))) {
                System.err.println("Error - '" + line[0]
                                   + "' is not a  valid partition setting (use T for Training, V for Validation, and S for test set)");
                System.exit(4);
              }
              partitionIndices[line[0].equals("T") ? 0 : (line[0].equals("V") ? 1 : 2)].add(count);
            } else {
              dataline[i - (manualPartition ? 1 : 0)] = Double.parseDouble(line[i]);
            }
          } catch (NumberFormatException nfe) {
            System.err.println("Error - '" + line[i]
                               + "' is not a valid integer, and we're not currently not set up to process missing fields");
            System.exit(5);
          }
        }
        records.add(dataline);
        count++;
      }
      reader.close();

      targets = Doubles.toArray(dependents);

      predictors = new double[count][];
      for (int i = 0; i < count; i++) {
        predictors[i] = records.elementAt(i);
      }

      if (manualPartition) {
        partitions = new int[3][];
        for (int i = 0; i < partitions.length; i++) {
          partitions[i] = Ints.toArray(partitionIndices[i]);
        }
      } else {
        partitions = createPartitions(records.size());
      }
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static int[][] createPartitions(int n) {
    return createPartitions(n, PARTITION_DEFAULTS);
  }

  public static int[][] createPartitions(int n, int[] percentages) {
    int[][] partitions = new int[3][];
    int[] randIndices;
    int count;

    if (percentages.length < 1 || percentages.length > 3) {
      System.err.println("Error - partition percentages array has to have a size within the range of 1-3");
      System.exit(1);
    }
    while (percentages.length < 3) {
      percentages = Array.addIntToArray(0, percentages);
    }
    if (Array.sum(percentages) != 100) {
      System.err.println("Error - partition percentages do not sum to 100");
      System.exit(2);
    }
    partitions[0] = new int[(int) Math.round((double) percentages[0] * n / 100)];
    partitions[1] = new int[(int) Math.round((double) (percentages[1]) * n / 100)];
    partitions[2] = new int[n - partitions[0].length - partitions[1].length];

    randIndices = Array.random(n);
    count = 0;
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < partitions[j].length; i++) {
        partitions[j][i] = randIndices[count];
        count++;
      }
    }

    return partitions;
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "kNNdata5.txt";
    kNN knn;

    String usage = "\n" + "park.kNN requires 0-1 arguments\n" + "   (1) filename (i.e. file="
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
      knn = new kNN(filename);
      // knn.setK(5);
      knn.calc();
      knn.report();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
