package org.genvisis.kaput;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class genoOnIBD {
  public genoOnIBD(String mutation, String dumpfile, String pos, String pre) throws IOException {
    BufferedReader reader;
    PrintWriter writer;
    String[] line = null, pair;
    String temp, trav, prev;
    Hashtable<String, String> hash = new Hashtable<String, String>();
    Hashtable<String, Vector<String>> preIDs;
    boolean done;
    Vector<String> ids, founders;
    String geno1, geno2;
    if (!new File(mutation).exists()) {
      System.err.println("Error - could not find " + mutation + " in current directory");
      System.exit(2);
    }
    reader = new BufferedReader(new FileReader(mutation));
    reader.readLine();
    while (reader.ready()) {
      line = reader.readLine().split("[\\s]+");
      hash.put(line[1] + "\t" + line[2],
               (Integer.valueOf(line[3]).intValue() + Integer.valueOf(line[4]).intValue() - 2)
                                         + "");
    }
    reader.close();

    if (!new File(pre).exists()) {
      System.err.println("Error - could not find " + pre + " in current directory");
      System.exit(2);
    }
    reader = new BufferedReader(new FileReader(pre));
    preIDs = new Hashtable<String, Vector<String>>();
    ids = new Vector<String>();
    founders = new Vector<String>();
    prev = "";
    done = false;
    while (!done) {
      if (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        trav = line[0];
      } else {
        done = true;
        trav = "";
      }
      if (!trav.equals(prev)) {
        if (founders.size() != 2 && prev != "") {
          System.err.println("Error - there should only be sibpairs in file " + pre
                             + "; problem encountered with family " + prev);
          System.exit(3);
        }
        preIDs.put(prev, ids);
        ids = new Vector<String>();
        founders.removeAllElements();
        prev = trav;
      }
      if (line[2].equals("0") && line[3].equals("0")) {
        if (!founders.contains(line[1])) {
          founders.add(line[1]);
        }
      } else {
        ids.add(line[1]);
        if (!founders.contains(line[2])) {
          founders.add(line[2]);
        }
        if (!founders.contains(line[3])) {
          founders.add(line[3]);
        }
      }
    }
    reader.close();

    if (!new File(dumpfile).exists()) {
      System.err.println("Error - could not find " + dumpfile + " in current directory");
      System.exit(2);
    }
    reader = new BufferedReader(new FileReader(dumpfile));
    writer = new PrintWriter(new FileWriter("regressed_" + dumpfile));

    while (reader.ready()) {
      line = reader.readLine().split("[\\s]+");
      if (Math.abs(Double.valueOf(line[0]).doubleValue()
                   - Double.valueOf(pos).doubleValue()) < 0.5) {
        ids = preIDs.get(line[1]);
        pair = line[2].split("-");
        for (int i = 0; i <= 1; i++) {
          if (Integer.valueOf(pair[i]).intValue() <= ids.size()) {
            pair[i] = ids.elementAt(Integer.valueOf(pair[i]).intValue() - 1);
          } else {
            System.err.println("Error - Individual " + line[1] + "-" + pair[i]
                               + " was not found in " + pre);
          }
        }
        geno1 = hash.get(line[1] + "\t" + pair[0]);
        geno2 = hash.get(line[1] + "\t" + pair[1]);
        temp = ".";
        // if (Integer.valueOf(geno1).intValue() < 0 ||
        // Integer.valueOf(geno2).intValue() < 0) {
        // temp = ".";
        // }
        if (Integer.valueOf(geno1).intValue() == 0 && Integer.valueOf(geno2).intValue() == 0) {
          temp = "0";
        }
        if (Integer.valueOf(geno1).intValue() > 0 && Integer.valueOf(geno2).intValue() > 0) {
          temp = "1";
        }
        if (Integer.valueOf(geno1).intValue() == 2 && Integer.valueOf(geno2).intValue() == 2) {
          temp = "2";
        }
        writer.print(line[1] + "\t" + pair[0] + "\t" + geno1 + "\t" + line[1] + "\t" + pair[1]
                     + "\t" + geno2 + "\t" + line[3] + "\t" + line[4] + "\t" + line[5]);
        writer.println("\t" + temp + "\t"
                       + ext.formDeci(Double.valueOf(line[4]).doubleValue() / 2
                                      + Double.valueOf(line[5]).doubleValue(), 6, true));
      }
    }

    reader.close();
    writer.close();

    // VectorI dependent = new DenseVector(dep);
    // MatrixI independent = new RowArrayMatrix(indep, true);
    // LinearRegressionI reg = new ReverseLinear(dependent, independent);

    // Print coefficients
    // VectorI coef = reg.getCoefficients();
    // System.out.println("Intercept Coefficient = "+coef.elementAt(2));
    // System.out.println("Column-0 Coefficient = "+coef.elementAt(0));
    // System.out.println("Column-1 Coefficient = "+coef.elementAt(1));

    // for (int i=0; i<numIDs; i++) {
    // writer.println(dep[i] - (coef.elementAt(2) +
    // coef.elementAt(0)*indep[i][0] + coef.elementAt(1)*indep[i][1]));
    // }
  }

  public static void main(String[] args) throws IOException {
    String mutation = "BRI3.dat", dump = "mibds02.dat", pos = "233", pre = "re_chrom02.pre";
    String usage = "\n" + "park.genoOnIBD requires 3 arguments:\n"
                   + "   (1) a chromosome#.dat-like file with genotypes of the putative mutation (i.e. mut="
                   + mutation + " (default))\n" + "   (2) a mapmaker dumpIBD file (i.e. dump="
                   + dump + " (default))\n" + "   (3) position to regress (i.e. pos=" + pos
                   + " (default))\n"
                   + "   (4) the (re)chrom##.pre file used to create the dump (i.e. pre=" + pre
                   + " (default))\n" + "";
    int numArgs = args.length;

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("mut=")) {
        mutation = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("dump=")) {
        dump = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("pos=")) {
        pos = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("pre=")) {
        pre = args[i].split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    System.out.println("Regressing the counts of " + mutation + " on the IBD estimates of " + dump
                       + " at position " + pos + " using the IDs from " + pre);
    try {
      new genoOnIBD(mutation, dump, pos, pre);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
