package org.genvisis.gwas;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class Haploview {
  public static final String[] LD_HEADER = {"L1", "L2", "D'", "r^2"};

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "ld_info.txt";
    String order = "order.txt";
    int column = 2;

    String usage = "\n" + "gwas.Haploview requires 0-1 arguments\n"
        + "   (1) name of file contianing pairwise LD information (i.e. file=" + filename
        + " (default))\n" + "   (2) name of file contianing the order to be presented (i.e. order="
        + order + " (default))\n"
        + "   (3) column number containing the data to be presented (i.e. col=" + column
        + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("order=")) {
        order = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("col=")) {
        column = ext.parseIntArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // filename = "D:\\tWork\\Consortium\\hg19\\phenos\\SNCA\\LD\\hg18\\ld_info.txt";
    // order = "D:\\tWork\\Consortium\\hg19\\phenos\\SNCA\\LD\\hg18\\order.txt";
    // column = 2;
    // column = 3;

    try {
      parseLD(filename, order, column);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseLD(String filename, String order, int column) {
    Hashtable<String, String> hash;
    String[] markerNames;
    PrintWriter writer;

    markerNames = HashVec.loadFileToStringArray(order, false, new int[] {0}, false);
    hash = HashVec.loadFileToHashString(filename, new int[] {0, 1}, new int[] {column}, false, "\t",
        true, false, false);
    HashVec.mergeHash2IntoHash1(hash, HashVec.loadFileToHashString(filename, new int[] {1, 0},
        new int[] {column}, false, "\t", true, false, false));

    try {
      writer =
          new PrintWriter(new FileWriter(ext.rootOf(order, false) + "_LD_col" + column + ".xln"));
      writer.println("D'");
      for (String markerName : markerNames) {
        writer.print("\t" + markerName);
      }
      writer.println();
      for (int i = 0; i < markerNames.length; i++) {
        writer.print(markerNames[i]);
        for (int j = 0; j < markerNames.length; j++) {
          writer.print("\t");
          if (j > i) {
            writer.print(hash.get(markerNames[i] + "\t" + markerNames[j]));
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      System.err
          .println("Error writing to " + ext.rootOf(order, false) + "_LD_col" + column + ".xln");
      e.printStackTrace();
    }
  }
}
