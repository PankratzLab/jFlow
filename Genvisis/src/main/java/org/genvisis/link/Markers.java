package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class Markers {
  public static final String DEFAULT_MARKER_DATABASE = "marker.database";

  public static final String[] ALT_LOCS =
      {"C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\park\\", "/home/npankrat/"};

  public static int[] order(String[] markers, boolean shouldBeOnSameChromosome) {
    return order(markers, DEFAULT_MARKER_DATABASE, shouldBeOnSameChromosome);
  }

  public static int[] order(String[] markers, String databaseFile,
      boolean shouldBeOnSameChromosome) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int[] chrs = new int[markers.length];
    double[] distances = Array.doubleArray(markers.length, -1);
    DecimalFormat myFormatter = new DecimalFormat("##0.00");
    int[] keys = new int[markers.length];
    int index, chr = -1;
    double diff;
    String output;

    try {
      reader = Files.getReader(databaseFile, ALT_LOCS);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line[0].equals("")) {
        } else if (line[0].startsWith("chromosome")) {
          chr = Integer.parseInt(line[0].substring(10));
        } else {
          index = ext.indexOfStr(line[0], markers);
          if (index != -1) {
            chrs[index] = chr;
            distances[index] = Double.parseDouble(line[chr == 23 ? 3 : 2]);
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + databaseFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + databaseFile + "\"");
      System.exit(2);
    }

    try {
      if (shouldBeOnSameChromosome) {
        keys = Sort.quicksort(distances);
        chr = chrs[keys[0]];
        output = markers[keys[0]];
        writer = new PrintWriter(new FileWriter("markerMap.dat"));
        for (int i = 0; i < keys.length; i++) {
          if (distances[keys[i]] == -1) {
            System.err
                .println("Error - marker '" + markers[i] + "' could not be found in the database");
          }
          if (chrs[keys[i]] != chr) {
            System.err.println("Error - all markers were not on the same chromosome ('"
                + markers[keys[i]] + "' for instance)");
          }
          if (i > 0) {
            diff = distances[keys[i]] - distances[keys[i - 1]];
            writer.print(myFormatter.format(diff) + " ");
            output += "\n     " + myFormatter.format(diff) + "\n" + markers[keys[i]] + "    "
                + distances[keys[i]];
          }
        }
        writer.println();
        writer.println(output);
        writer.close();
      } else {
        System.err.println(
            "Sorry, ordering markers on different chromosomes has not been implemented yet. Have at it!");
      }
    } catch (Exception e) {
      System.err.println("Error writing markerMap.dat");
      e.printStackTrace();
    }

    return keys;
  }
}
