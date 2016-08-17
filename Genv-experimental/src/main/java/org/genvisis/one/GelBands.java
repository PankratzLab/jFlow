package org.genvisis.one;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.ext;
import org.genvisis.stats.Maths;

public class GelBands {
  public static final String LADDER_HEADER = "Ladder header";
  public static final String LADDER_TAG = "100 bp ladder";
  public static final String[] SUFFIXES = {"X", "Y"};

  public static final int REFERENCE_FRAGMENT_SIZE = 1329;
  public static final int REFERENCE_REPEATS_SIZE = 896;
  public static final int SIZE_OF_UNIT = 32;

  public static void call(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, cell, values;
    Vector<String[]> v = new Vector<String[]>();
    // int count;
    boolean betweenLadders;
    int[] ladder;
    int[][][] ladderCoordinates;

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(filename + "_calls.xln"));
      betweenLadders = false;
      // count = 0;
      ladder = null;
      ladderCoordinates = new int[2][][];
      while (reader.ready()) {
        line = reader.readLine().trim().split("\\t", -1);
        if (line[0].equalsIgnoreCase(LADDER_TAG)) {
          if (betweenLadders) {
            ladderCoordinates[1] = parseLadderCoordinates(ladder, line);
            for (int i = 0; i < v.size(); i++) {
              writer.print(Array.toStr(v.elementAt(i)));
              values = translateIntoAlleles(ladder, ladderCoordinates, v.elementAt(i));
              writer.print("\t" + Array.toStr(values));
              values = translateAllelesIntoRepeats(values);
              writer.print("\t" + Array.toStr(values));
              writer.println();
            }
            v.removeAllElements();
            betweenLadders = false;
          } else {
            ladderCoordinates[0] = parseLadderCoordinates(ladder, line);
            betweenLadders = true;
            // count++;
          }
          writer.println(Array.toStr(line));
        } else if (betweenLadders) {
          v.add(line);
        } else if (line[0].equalsIgnoreCase(LADDER_HEADER)) {
          ladder = Array.intArray((line.length - 1) / 2, -1);
          for (int i = 0; i < ladder.length; i++) {
            for (int j = 0; j < 2; j++) {
              cell = line[i * 2 + 1 + j].split("[\\s]+");
              try {
                if (ladder[i] == -1) {
                  ladder[i] = Integer.parseInt(cell[0]);
                } else if (ladder[i] != Integer.parseInt(cell[0])) {
                  System.err.println("Error - mismatched ladder sizes ('" + line[i * 2 + 1 + 0]
                                     + "' and '" + line[i * 2 + 1 + 1]
                                     + "'); expecting the same size with X and Y suffixes");
                }
                if (!cell[1].equals(SUFFIXES[j])) {
                  System.err.println("Error - mismatched ladder format; expecting the same size with X and Y suffixes");
                }
              } catch (Exception e) {
                System.err.println("Error parsing ladder size: " + Array.toStr(cell));
                e.printStackTrace();
              }
            }
          }
          writer.println(Array.toStr(line));
        } else {
          System.out.println("ignoring line: " + Array.toStr(line));
          writer.println(Array.toStr(line));
        }

      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static int linearTranslation(int[] ladder, int[] customLadderCoordinates, double x) {
    int rung;
    double[] slopeAndIntercept;
    boolean inverse;

    inverse = customLadderCoordinates[0] - customLadderCoordinates[ladder.length - 1] > 0;

    rung = 1;
    if (inverse) {
      while (rung < customLadderCoordinates.length - 1 && x < customLadderCoordinates[rung]) {
        rung++;
      }
    } else {
      while (rung < customLadderCoordinates.length - 1 && x > customLadderCoordinates[rung]) {
        rung++;
      }
    }

    slopeAndIntercept = Maths.slopeAndIntercept(customLadderCoordinates[rung - 1], ladder[rung - 1],
                                                customLadderCoordinates[rung], ladder[rung]);

    return (int) (slopeAndIntercept[0] * x + slopeAndIntercept[1]);
  }

  public static void main(String[] args) {
    // int[] ladder = new int[] {900, 1000, 1200, 1517};
    // int[] ladderCoor = new int[] {188, 154, 88, 14};
    // double y = 88;
    //
    // System.out.println(linearTranslation(ladder, ladderCoor, y));
    //
    // System.exit(1);

    String filename =
        "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\DOCK5\\Bill's replication\\points.txt";

    try {
      call(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static int[][] parseLadderCoordinates(int[] ladder, String[] values) {
    int[][] ladderCoordinates;

    if (!values[0].equalsIgnoreCase(LADDER_TAG)) {
      System.err.println("Error - invalid ladder tag: '" + values[0] + "'");
      System.exit(1);
    }
    if (values.length != ladder.length * 2 + 1) {
      System.err.println("Error - invalid number of ladder coordinates (need X and Y for each rung defined in the ladder header)");
      System.exit(1);
    }

    ladderCoordinates = new int[ladder.length][2];
    for (int i = 0; i < ladder.length; i++) {
      for (int j = 0; j < 2; j++) {
        try {
          ladderCoordinates[i][j] = Integer.parseInt(values[1 + i * 2 + j]);
        } catch (Exception e) {
          System.err.println("Error parsing " + SUFFIXES[j] + " coordinate for " + ladder[i]
                             + " bp ladder rung : " + values[1 + i * 2 + j]);
          e.printStackTrace();
        }
      }
    }

    return ladderCoordinates;
  }

  public static String[] translateAllelesIntoRepeats(String[] alleles) {
    String[] repeats;

    repeats = new String[alleles.length];
    for (int i = 0; i < alleles.length; i++) {
      if (alleles[i].equals(".")) {
        repeats[i] = ".";
      } else {
        repeats[i] = ext.formDeci(
                                  (double) (Integer.parseInt(alleles[i])
                                            - (REFERENCE_FRAGMENT_SIZE - REFERENCE_REPEATS_SIZE))
                                  / (double) SIZE_OF_UNIT, 2);
      }
    }

    return repeats;
  }

  public static String[] translateIntoAlleles(int[] ladder, int[][][] ladderCoordinates,
                                              String[] values) {
    String[] alleles;
    int[] ys;
    int[] customLadderCoordinates;
    double[] slopeAndIntercept;
    int numPoints;
    double meanX, count, d;

    d = (values.length - 1.0) / 2.0;
    if (d - Math.floor(d) > 0.0001) {
      System.err.println("Error - odd number of coordinates after identifier '" + values[0]
                         + "'; translation requires pairs of coordinates");
    }

    numPoints = (int) Math.floor(d);
    ys = new int[numPoints];
    meanX = 0;
    count = 0;
    for (int i = 0; i < numPoints; i++) {
      try {
        if (values[i * 2 + 1 + 0].equals(".") || values[i * 2 + 1 + 1].equals(".")) {
          ys[i] = Integer.MIN_VALUE;
        } else {
          ys[i] = Integer.parseInt(values[i * 2 + 1 + 1]);

          meanX += Integer.parseInt(values[i * 2 + 1 + 0]);
          count++;
        }
      } catch (Exception e) {
        System.err.println("Error parsing point " + (i + 1) + " for sample " + values[0] + ": "
                           + values[i * 2 + 1 + 0] + "," + values[i * 2 + 1 + 1]);
      }
    }
    meanX /= count;

    customLadderCoordinates = new int[ladder.length];
    for (int i = 0; i < ladder.length; i++) {
      slopeAndIntercept =
          Maths.slopeAndIntercept(ladderCoordinates[0][i][0], ladderCoordinates[0][i][1],
                                  ladderCoordinates[1][i][0], ladderCoordinates[1][i][1]);
      customLadderCoordinates[i] = (int) (slopeAndIntercept[0] * meanX + slopeAndIntercept[1]);
    }

    alleles = new String[numPoints];
    for (int i = 0; i < numPoints; i++) {
      if (ys[i] == Integer.MIN_VALUE) {
        alleles[i] = ".";
      } else {
        alleles[i] = linearTranslation(ladder, customLadderCoordinates, ys[i]) + "";
      }
    }
    //
    // try {
    // writer = new PrintWriter(new FileWriter("ladders.xln", true));
    // for (int i = 0; i<ladder.length; i++) {
    // writer.println(customLadderCoordinates[i]+"\t"+ladder[i]);
    // }
    // writer.println();
    // writer.close();
    // } catch (Exception e) {
    // System.err.println("Error writing to "+"ladders.xln");
    // e.printStackTrace();
    // }

    return alleles;
  }
}
