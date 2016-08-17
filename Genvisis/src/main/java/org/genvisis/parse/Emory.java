// -Xms1024M -Xmx1024M
package org.genvisis.parse;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;

public class Emory {
  public static final String[] HEADER_STARTER = {"IID", "LabID", "SEX", "Case"};

  public static void main(String[] args) {
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\UMN\\Myron\\ExcisionPathway\\00src\\";
    // String filename = "map2.dat";
    // String filename = "map1.dat";
    // String filename = "cpru.dat";
    // String filename = "cpru2.dat";

    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\UMN\\Myron\\CARDIA\\00src\\";
    // String filename = "yalta.dat";

    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\UMN\\Myron\\Indian_Diabetes\\00src\\";
    // String filename = "firstSet.dat";
    // String filename = "second11.dat";

    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Myron\\ExcisionPathway\\finalScores\\00src\\";
    String filename = "pooled.dat";

    try {
      parse(dir, filename, new Logger());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse(String dir, String filename, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String[] markerNames = null;
    SnpMarkerSet map;

    try {
      System.out.println("Writing pedigree file...");
      reader = new BufferedReader(new FileReader(dir + filename));
      writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename) + ".ped"));
      line = reader.readLine().trim().split("[\\s]+");
      ext.checkHeader(Array.subArray(line, 0, 4), HEADER_STARTER, false, true);
      markerNames = Array.subArray(line, 4);
      while (reader.ready()) {
        line = reader.readLine().split("\t", -1);
        if (line.length != markerNames.length + 4) {
          System.err.println("Error - mismatched number of columns for record '" + line[0] + "'");
        }
        writer.print(line[0] + "\t" + line[1] + "\t0\t0\t" + line[2] + "\t"
            + (line[3].equals(".") ? "0" : (Integer.parseInt(line[3]) + 1)));
        for (int i = 0; i < markerNames.length; i++) {
          if (line[4 + i].equals("")) {
            writer.print("\t0\t0");
          } else {
            writer.print("\t" + line[4 + i].charAt(0) + "\t"
                + line[4 + i].charAt(line[4 + i].length() == 2 ? 1 : 0));
          }
        }
        writer.println();
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    map = new SnpMarkerSet(markerNames);
    map.parseSNPlocations();
    map.writeToFile(dir + ext.rootOf(filename) + ".map", SnpMarkerSet.PLINK_MAP_FORMAT, log);
  }
}
