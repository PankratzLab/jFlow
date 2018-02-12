package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class Bootstrap {

  public Bootstrap(String filename) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;

    int count = 0;
    double[] array, bs;
    int column = 1, reps = 1000;
    double mean;

    try {
      reader = new BufferedReader(new FileReader(filename));
      line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
      try {
        if (line.length != 2) {
          reader.close();
          throw new Exception();
        }
        column = Integer.parseInt(line[0]);
        reps = Integer.parseInt(line[1]);
      } catch (Exception e) {
        System.err.println("Error - first line of file must contain exactly two integers: column to bootstrap and number of replicates");
        System.exit(1);
      }
      while (reader.ready()) {
        reader.readLine();
        count++;
      }
      System.out.println("Using all " + count + " records.");
      System.out.println("Bootstrapping column " + column + " a total of " + reps + " times.");
      reader.close();

      reader = new BufferedReader(new FileReader(filename));
      reader.readLine();
      array = new double[count];
      count = 0;
      while (reader.ready()) {
        line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
        array[count] = Double.parseDouble(line[column - 1]);
        count++;
      }
      mean = ArrayUtils.mean(array);
      System.out.println("mean: " + mean + " (chi-square, 1 df p-value: "
                         + ext.formDeci(ProbDist.ChiDist(mean, 1), 4, true) + ")");
      bs = ArrayUtils.bootstrap(array, reps, true);
      System.out.println("bootstrapped mean (95% CI):\t" + ext.formDeci(bs[0], 4, true) + " ("
                         + ext.formDeci(bs[1], 4, true) + "," + ext.formDeci(bs[2], 4, true) + ")");
      System.out.println("corresponding chi-square, 1 df p-values:\t"
                         + ext.formDeci(ProbDist.ChiDist(bs[0], 1), 4, true) + " ("
                         + ext.formDeci(ProbDist.ChiDist(bs[1], 1), 4, true) + ","
                         + ext.formDeci(ProbDist.ChiDist(bs[2], 1), 4, true) + ")");
      reader.close();

      writer = Files.openAppropriateWriter(filename + "-bootsrapped.out");
      writer.println("Using all " + count + " records.");
      writer.println("Bootstrapping column " + column + " a total of " + reps + " times.");
      writer.println();
      writer.println("mean: " + ArrayUtils.mean(array));
      writer.println("bootstrapped mean (95% CI):\t" + ext.formDeci(bs[0], 4, true) + " ("
                     + ext.formDeci(bs[1], 4, true) + "," + ext.formDeci(bs[2], 4, true) + ")");
      writer.println("corresponding chi-square, 1 df p-values:\t"
                     + ext.formDeci(ProbDist.ChiDist(bs[0], 1), 4, true) + " ("
                     + ext.formDeci(ProbDist.ChiDist(bs[1], 1), 4, true) + ","
                     + ext.formDeci(ProbDist.ChiDist(bs[2], 1), 4, true) + ")");
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    // String filename = "Bootstrap.dat";
    String filename = "Survival_Rep1-ALL_VPD.500K_summary.out";

    String usage = "\n" + "park.Bootstrap requires 0-1 arguments\n" + "   (1) filename (i.e. file="
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
      new Bootstrap(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
