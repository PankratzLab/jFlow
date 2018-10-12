package org.genvisis.seq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

public class Mapability {

  public static final String[] REQS = {"Chr", "Position", "totalreads", "notduplicated",
                                       "brokenmatepairs", "mapqualitygtzero", "notproperlypaired"};

  public static void parseBeds(String filename) {
    PrintWriter writer;
    String[] header;
    String[][] data;
    int[] indices;
    int pos;
    double value, sum, mean, max;

    header = Files.getHeaderOfFile(filename, "\t", new Logger());
    data = HashVec.loadFileToStringMatrix(filename, true, ArrayUtils.arrayOfIndices(header.length));
    indices = ext.indexFactors(REQS, header, false);

    max = -1;
    sum = 0;
    for (String[] element : data) {
      value = Double.parseDouble(element[indices[3]]);
      sum += value;
      if (value > max) {
        max = value;
      }
    }
    mean = sum / data.length;

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + "_coverage.bed");
      writer.println("track name=totalUniqueReads description=\"Coverage relative to mean (mean=500)\" useScore=1");
      for (String[] element : data) {
        pos = Integer.parseInt(element[indices[1]]);
        writer.println("chr" + element[indices[0]] + " " + (pos - 50) + " " + (pos + 50) + " " + pos
                       + " " + (int) Math.min(Double.parseDouble(element[indices[3]]) / mean * 500,
                                              1000));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_coverage.bed");
      e.printStackTrace();
    }

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + "_brokenMatePairs.bed");
      writer.println("track name=percentBrokenMatePairs description=\"Percent of unique reads with broken mate pairs (scaled from 0-1000) \" useScore=1");
      for (String[] element : data) {
        pos = Integer.parseInt(element[indices[1]]);
        writer.println("chr" + element[indices[0]] + " " + (pos - 50) + " " + (pos + 50) + " " + pos
                       + " " + (int) (Double.parseDouble(element[indices[4]])
                                      / Double.parseDouble(element[indices[3]]) * 1000));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_coverage.bed");
      e.printStackTrace();
    }

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + "_mapQualityZero.bed");
      writer.println("track name=percentMapQualityZero description=\"Percent of unique reads with map quality of zero (scaled from 0-1000) \" useScore=1");
      for (String[] element : data) {
        pos = Integer.parseInt(element[indices[1]]);
        writer.println("chr" + element[indices[0]] + " " + (pos - 50) + " " + (pos + 50) + " " + pos
                       + " "
                       + (int) ((Double.parseDouble(element[indices[3]])
                                 - Double.parseDouble(element[indices[5]]))
                                / Double.parseDouble(element[indices[3]]) * 1000));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_coverage.bed");
      e.printStackTrace();
    }

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + "_notProperlyPaired.bed");
      writer.println("track name=percentNotProperlyPaired description=\"Percent of unique reads that are not properly paired (scaled from 0-1000) \" useScore=1");
      for (String[] element : data) {
        pos = Integer.parseInt(element[indices[1]]);
        writer.println("chr" + element[indices[0]] + " " + (pos - 50) + " " + (pos + 50) + " " + pos
                       + " " + (int) (Double.parseDouble(element[indices[6]])
                                      / Double.parseDouble(element[indices[3]]) * 1000));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_coverage.bed");
      e.printStackTrace();
    }

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false)
                                           + "_compositePoorQualty.bed");
      writer.println("track name=CompositePoorQuality description=\"Composite Poor Quality = MapQualityZero + NotProperlyPaired\" useScore=1");
      for (String[] element : data) {
        pos = Integer.parseInt(element[indices[1]]);
        writer.println("chr" + element[indices[0]] + " " + (pos - 50) + " " + (pos + 50) + " " + pos
                       + " "
                       + (int) ((Double.parseDouble(element[indices[3]])
                                 - Double.parseDouble(element[indices[5]])
                                 + Double.parseDouble(element[indices[6]]))
                                / Double.parseDouble(element[indices[3]]) * 1000));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_coverage.bed");
      e.printStackTrace();
    }
  }

  public static void assignPass(String filename, String[] passes) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String[][] data;
    Segment[][] segs;
    Segment seg;
    boolean overlaps;

    segs = new Segment[passes.length][];
    for (int i = 0; i < passes.length; i++) {
      data = HashVec.loadFileToStringMatrix(passes[i], false, new int[] {0, 1, 2});
      segs[i] = new Segment[data.length];
      for (int j = 0; j < data.length; j++) {
        segs[i][j] = new Segment("chr" + data[j][0] + ":" + data[j][1] + "-" + data[j][2]);
      }
    }
    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = Files.openAppropriateWriter(filename + "_positions.xln");
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      ext.checkHeader(line, new String[] {"Chr", "Position"}, new int[] {0, 1}, false, new Logger(),
                      true);
      for (int i = 0; i < passes.length; i++) {
        line = ArrayUtils.insertStringAt(ext.rootOf(passes[i]), line, 2 + i);
      }
      writer.println(ArrayUtils.toStr(line));
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        seg = new Segment("chr" + line[0] + ":" + line[1] + "-" + (Integer.parseInt(line[1]) + 1));
        for (int i = 0; i < passes.length; i++) {
          overlaps = false;
          for (int j = 0; j < segs[i].length; j++) {
            if (seg.overlaps(segs[i][j])) {
              overlaps = true;
            }
          }
          line = ArrayUtils.insertStringAt(overlaps ? "1" : "0", line, 2 + i);
        }
        writer.println(ArrayUtils.toStr(line));
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "D:\\tWork\\SequencingProjectWithCIDR\\Coverage\\all.out";
    String filename = "positions.txt";
    boolean parse = false;
    String[] passes = new String[] {"1stPass.bed", "2ndPass.bed"};

    String usage = "\n" + "seq.Mapability requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n"
                   + "   (2) parse beds (i.e. -parse (not the default))\n" + " OR\n"
                   + "   (2) assign passes to file (i.e. passes=passOne.bed,passTwo.bed (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-parse")) {
        parse = true;
        numArgs--;
      } else if (arg.startsWith("passes=")) {
        passes = arg.split("=")[1].split(",");
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (parse) {
        parseBeds(filename);
      } else if (passes != null) {
        assignPass(filename, passes);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
