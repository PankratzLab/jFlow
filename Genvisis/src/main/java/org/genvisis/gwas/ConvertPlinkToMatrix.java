package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Logger;
import org.genvisis.filesys.SnpMarkerSet;

public class ConvertPlinkToMatrix {
  public static void convert(String dir, String root) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String[] markerNames;

    markerNames = new SnpMarkerSet(dir + root + ".map", false, new Logger()).getMarkerNames();
    try {
      reader = new BufferedReader(new FileReader(dir + root + ".ped"));
      writer = new PrintWriter(new FileWriter(dir + root + ".matrix"));
      writer.print("IID");
      for (String markerName : markerNames) {
        writer.print("\t" + markerName);
      }
      writer.println();
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.print(line[1]);
        for (int i = 0; i < markerNames.length; i++) {
          writer.print("\t" + line[6 + 2 * i + 0] + line[6 + 2 * i + 1]);
        }
        writer.println();
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - file \"" + root + ".ped" + "\" not found in current directory");
      fnfe.printStackTrace();
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + root + ".ped" + "\"");
      ioe.printStackTrace();
      System.exit(2);
    }
  }

  public static void main(String[] args) {
    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\SequencingProjectWithCIDR\\SNP_comparisons\\";
    // String root = "48snps";
    String root = "768snps";

    try {
      convert(dir, root);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
