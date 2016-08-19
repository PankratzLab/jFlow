package org.genvisis.gwas;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class Haplotypes {

  // expects 4 columns: FID, IID, Hap1, Hap2
  public static void countsPerPerson(String filename) {
    PrintWriter writer;
    String[] keys;
    Hashtable<String, String> hash = new Hashtable<String, String>();

    String[][] data = HashVec.loadFileToStringMatrix(filename, true, new int[] {0, 1, 2, 3}, false);
    hash = new Hashtable<String, String>();
    for (String[] element : data) {
      for (int j = 0; j < 2; j++) {
        hash.put(element[2 + j], "");
      }
    }
    try {
      keys = HashVec.getKeys(hash);
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_counts.xln"));
      writer.println("FID\tIID\t" + Array.toStr(keys));
      int[] counts;
      for (int i = 0; i < data.length; i++) {
        counts = new int[keys.length];
        for (int j = 0; j < 2; j++) {
          counts[ext.indexOfStr(data[i][2 + j], keys)]++;
        }
        writer.println(data[i][0] + "\t" + data[i][1] + "\t" + Array.toStr(counts));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename, false) + "_counts.xln");
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename =
                    "D:\\umn\\Myron\\CARe\\Pathway analysis\\results\\haplo.phase-CAGCA_parsed.xln";

    String usage = "\n" + "gwas.Haplotypes requires 0-1 arguments\n" + "   (1) filename (i.e. file="
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
      countsPerPerson(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
