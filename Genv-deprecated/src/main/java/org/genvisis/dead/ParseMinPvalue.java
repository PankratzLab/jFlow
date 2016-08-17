// -Xms1024M -Xmx1024M
package org.genvisis.dead;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

import com.google.common.primitives.Ints;

public class ParseMinPvalue {
  // public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\merged\\results\\";
  public static final String DEFAULT_DIR =
      "C:\\Documents and Settings\\npankrat\\My Documents\\ADNI\\results\\";

  // public static final String[] METHODS = {"linear", "logistic"};
  public static final String[] MODELS = {"Additive", "Dominant", "Recessive", "Genotypic"};
  // public static final String[] MODELS = {"Additive", "Recessive"};
  // public static final String[] MODELS = {"Additive"};

  public static final String[] METHODS = {"logistic"};

  public static final String TAG = "p-value";

  public static final int TOP_N = 100;

  public static final boolean USE_ASSOC = false;

  public static void parseMinPvalue(String dir, String filename) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, marks;
    Hashtable<String, String> hash = new Hashtable<String, String>();
    Vector<String> markers = new Vector<String>();
    IntVector iv;
    int[] indices, keys;
    double min, d;
    int count;
    double[] ps;

    System.out.println("Loading marker names...");
    if (new File(dir + "markers.dat").exists()) {
      markers = HashVec.loadFileToVec(dir + "markers.dat", true, true, false);
    } else {
      markers = HashVec.loadFileToVec("C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\"
                                      + "markers.dat", true, true, false);
    }
    if (USE_ASSOC) {
      System.out.println("Loading results of assoc test...");
      hash = HashVec.loadFileToHashString(dir + "plink.assoc", 1, new int[] {8}, "\t", true);
    } else {
      System.out.println("Ignoring the assoc test...");
    }
    for (int i = 0; i < METHODS.length; i++) {
      try {
        System.out.println("Loading " + METHODS[i] + " results...");
        reader = new BufferedReader(new FileReader(dir + METHODS[i] + ".xls"));
        line = reader.readLine().split("[\\s]+");
        iv = new IntVector();
        for (int j = 0; j < line.length; j++) {
          if (line[j].equals(TAG) && ext.indexOfStr(line[j - 2], MODELS) != -1) {
            iv.add(j);
          }
        }
        indices = Ints.toArray(iv);
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          if (USE_ASSOC) {
            min = Double.parseDouble(hash.get(line[0]));
          } else {
            min = 1;
          }
          for (int j = 0; j < indices.length; j++) {
            if (!line[indices[j]].equals(".")) {
              d = Double.parseDouble(line[indices[j]]);
              if (d < min) {
                min = d;
              }
            }
          }
          hash.put(line[0], min + "");
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + METHODS[i] + ".xls"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + METHODS[i] + ".xls" + "\"");
        System.exit(2);
      }
    }

    writer = Files.getWriter(dir + filename);
    writer.println("SNP\tmin p-value\t-log10 p");
    count = 0;
    for (int i = 0; i < markers.size(); i++) {
      if (hash.containsKey(markers.elementAt(i))) {
        writer.println(markers.elementAt(i) + "\t" + hash.get(markers.elementAt(i)) + "\t"
                       + -1 * Math.log10(Double.parseDouble(hash.get(markers.elementAt(i)))));
        count++;
      }
    }
    writer.close();

    marks = new String[count];
    ps = new double[count];
    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      reader.readLine();
      for (int i = 0; i < count; i++) {
        line = reader.readLine().split("[\\s]+");
        marks[i] = line[0];
        ps[i] = Double.parseDouble(line[1]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    hash = HashVec.loadFileToHashString(dir + "plink.assoc", 1, new int[] {0}, "\t", false);
    try {
      keys = Sort.quicksort(ps);
      count = 0;
      writer = new PrintWriter(new FileWriter(dir + "Top" + TOP_N + ".txt"));
      for (int i = 0; count < TOP_N; i++) {
        if (Integer.parseInt(hash.get(marks[keys[i]])) < 23) {
          writer.println(marks[keys[i]]);
          count++;
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error determining top " + TOP_N + " list.");
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = DEFAULT_DIR;
    String filename = "min_pvalue.dat";

    String usage = "\\n" + "park.gwa.ParseMinPvalue requires 0-2 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("file=")) {
        filename = args[i].split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      parseMinPvalue(dir, filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
