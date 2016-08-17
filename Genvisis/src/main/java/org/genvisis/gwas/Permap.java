package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.ext;
import org.genvisis.mining.Transformations;

public class Permap {
  public static final String DIR =
      "C:\\Documents and Settings\\npankrat\\My Documents\\PD-singleton\\";

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "permap.dat";

    String usage = "\n" + "park.permap requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n" + "";

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
      makePermapFile(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void makePermapFile(String filename) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    Vector<String> v = new Vector<String>();
    String[] peeps;
    double[] values;
    double mean, stdev;
    int count;

    try {
      reader = new BufferedReader(new FileReader(DIR + "pdXY.fam"));
      while (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        v.add(line[1]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "pdXY.fam" + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "pdXY.fam" + "\"");
      System.exit(2);
    }

    peeps = Array.toStringArray(v);
    try {
      // reader = new BufferedReader(new
      // FileReader(DIR+"plink.mdist.missing"));
      // writer = new PrintWriter(new
      // FileWriter(DIR+"plink_similarity.txt"));
      // writer.println("TITLE= "+"plink.mdist.missing");
      reader = new BufferedReader(new FileReader(DIR + "plink.mibs"));
      writer = new PrintWriter(new FileWriter(DIR + "plink.mibs.txt"));
      writer.println("TITLE= " + "plink.mibs");
      writer.println("NOBJECTS= " + peeps.length);
      writer.println();
      writer.println("SIMILARITYLIST");

      count = 0;
      values = new double[peeps.length * (peeps.length - 1) / 2];
      for (int i = 0; i < peeps.length; i++) {
        line = reader.readLine().split("[\\s]+");
        for (int j = 0; j < i; j++) {
          values[count++] = Double.parseDouble(line[j]);
        }
      }

      mean = Array.mean(values);
      stdev = Array.stdev(values);
      System.out.println(values.length);

      System.out.println("mean: " + mean);
      System.out.println("std: " + stdev);
      System.out.println("range: " + Array.min(values) + " - " + Array.max(values));

      values = Transformations.standardizeRange(values);
      mean = Array.mean(values);
      stdev = Array.stdev(values);

      System.out.println("mean: " + mean);
      System.out.println("std: " + stdev);
      System.out.println("range: " + Array.min(values) + " - " + Array.max(values));

      count = 0;
      for (int i = 0; i < peeps.length; i++) {
        writer.print(peeps[i]);
        for (int j = 0; j < i; j++) {
          writer.print("\t" + ext.formDeci(values[count++], 5));
        }
        writer.println("\t" + 1);
      }
      writer.println(1);
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + DIR + "plink.mdist.missing"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + DIR + "plink.mdist.missing" + "\"");
      System.exit(2);
    }

  }
}
