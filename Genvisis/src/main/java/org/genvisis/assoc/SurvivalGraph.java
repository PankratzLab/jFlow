// requires a file with 3 columns and no header: FamID, CLASS, AGE

package org.genvisis.assoc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class SurvivalGraph {
  public static void graph(String filename) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    String temp;
    Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
    Vector<String> v, classV, ageV;
    String[] fams, classes, ages;
    double[][] counts;
    int[] keys;
    double[] sums, freqs;

    classV = new Vector<String>();
    ageV = new Vector<String>();
    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.split("\t");
        if (line.length != 3) {
          System.err.println("Error - line does not have exactly 3 tab-delimited columns:\n"
                             + temp);
        }
        if (!line[0].equals(".") && !line[1].equals(".") && !line[2].equals(".")) {
          if (hash.containsKey(line[0])) {
            v = hash.get(line[0]);
          } else {
            hash.put(line[0], v = new Vector<String>());
          }
          v.add(line[1] + "\t" + line[2]);
          HashVec.addIfAbsent(line[1], classV, false);
          HashVec.addIfAbsent(line[2], ageV, false);
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

    classes = Array.toStringArray(classV);
    ages = Array.toStringArray(ageV);
    keys = Sort.quicksort(Array.toDoubleArray(ages));
    fams = HashVec.getKeys(hash);
    counts = new double[classes.length][ages.length];
    for (String fam : fams) {
      v = hash.get(fam);
      for (int j = 0; j < v.size(); j++) {
        line = (v.elementAt(j)).split("\t");
        counts[ext.indexOfStr(line[0], classes)][ext.indexOfStr(line[1], ages)] += 1.0 / v.size();
      }
    }

    sums = new double[counts.length];
    for (int i = 0; i < counts.length; i++) {
      sums[i] = Array.sum(counts[i]);
    }

    try {
      writer = new PrintWriter(new FileWriter(filename + "-survival.xls"));
      freqs = Array.doubleArray(classes.length, 1);
      writer.println("\t" + Array.toStr(classes));
      writer.println("0\t" + Array.toStr(Array.stringArray(classes.length, "1")));
      for (int j = 0; j < ages.length; j++) {
        writer.print(ages[keys[j]]);
        for (int i = 0; i < classes.length; i++) {
          freqs[i] -= counts[i][keys[j]] / sums[i];
          writer.print("\t" + ext.formDeci(freqs[i], 5, true));
        }
        writer.println();
      }
      writer.close();
    } catch (IOException ioe) {
      System.err.println("Error writing final file");
      ioe.printStackTrace();
    }

  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    // String filename = "test_survival-1.dat";
    // String filename = "graph_Rep1-263_VPD.dat";
    // String filename = "graph_Rep1-ALL_VPD.dat";
    // String filename = "graph_Rep1-263_VPD_noParkin.dat";
    // String filename = "graph_Rep1-ALL_VPD_noParkin.dat";
    String filename = "GBA_survivial_Graph.dat";

    String usage = "\n" + "assoc.SurvivalGraph requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default)\n" + "";

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
      graph(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
