package org.genvisis.widgets;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Date;
import java.util.LinkedList;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class peakAt {
  public static final int DEFAULT_NUM_LINES = 1000;

  public static void peak(String filename) throws Elision {
    String temp;
    int numLines;
    boolean columnsNotLines = false;
    boolean tail, counting;
    String outputFilename;

    temp = filename.substring(0, filename.length() - (".peakAt").length());
    tail = filename.indexOf(".tail.") > 0;
    temp = ext.replaceAllWith(temp, ".tail", "");
    counting = filename.indexOf(".count.") > 0 || filename.indexOf(".wc.") > 0;

    if (temp.indexOf(".") > 0) {
      try {
        if (temp.substring(temp.lastIndexOf(".") + 1).startsWith("cols")) {
          columnsNotLines = true;
          temp = temp.substring(0, temp.lastIndexOf("."));
        } else {
          columnsNotLines = false;
        }
      } catch (Exception e) {
        columnsNotLines = false;
      }

      try {
        numLines = Integer.parseInt(temp.substring(temp.lastIndexOf(".") + 1));
        temp = temp.substring(0, temp.lastIndexOf("."));
      } catch (NumberFormatException nfe) {
        numLines = DEFAULT_NUM_LINES;
      }
    } else {
      numLines = DEFAULT_NUM_LINES;
    }

    outputFilename =
        "PeakAt_" + (tail ? "last" : "first") + numLines + (columnsNotLines ? "column" : "line")
                     + (numLines == 1 ? "" : "s") + "_" + ext.removeDirectoryInfo(temp);

    actualPeak(filename, outputFilename, numLines, tail, counting, columnsNotLines);
  }

  public static void actualPeak(String filename, String outputFilename, int numLines, boolean tail,
                                boolean counting, boolean columnsNotLines) throws Elision {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    LinkedList<String> ll;
    int count;
    long time;
    InputStreamReader isReader;

    try {
      writer = Files.getAppropriateWriter(ext.parseDirectoryOfFile(filename)
                                          + (outputFilename.endsWith(".gz") ? outputFilename.substring(0,
                                                                                                       outputFilename.length()
                                                                                                          - 3)
                                                                            : outputFilename));
      if (counting) {
        System.out.println("Counting the number of rows.");
        time = new Date().getTime();
        count = Files.countLines(filename, 0);
        writer.println(count);
        writer.println("Counted " + count + " rows in " + ext.getTimeElapsed(time));
      } else {
        isReader = null;
        if (outputFilename.endsWith(".gz")) {
          isReader = new InputStreamReader(new GZIPInputStream(new FileInputStream(filename)));
        } else if (outputFilename.endsWith(".zip")) {
          isReader = new InputStreamReader(new ZipInputStream(new FileInputStream(filename)));
        } else {
          isReader = new FileReader(filename);
        }

        reader = new BufferedReader(isReader);
        if (columnsNotLines) {
          System.out.println("Taking the first " + numLines + " columns of all rows.");
          while (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
            for (int i = 0; i < Math.min(numLines, line.length); i++) {
              writer.print((i == 0 ? "" : "\t") + line[i]);
            }
            writer.println();
            writer.flush();
          }
        } else if (tail) {
          ll = new LinkedList<String>();
          while (reader.ready()) {
            ll.addLast(reader.readLine());
            if (ll.size() > numLines) {
              ll.removeFirst();
            }
          }
          while (!ll.isEmpty()) {
            writer.println(ll.removeFirst());
          }
        } else {
          for (int i = 0; i < numLines && reader.ready(); i++) {
            writer.println(reader.readLine());
          }
        }
        reader.close();
      }
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      throw new Elision();
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      throw new Elision();
    }

  }

  public static void fromParameters(String filename, Logger log) throws Elision {
    Vector<String> params;
    String trav;
    String inputFilename = "very_large_file.txt";
    String outputFilename = "fraction_of_large_file.txt";
    int numLines = 10000;
    boolean tail = false;
    boolean counting = false;
    boolean columnsNotLines = false;

    params =
        Files.parseControlFile(filename, "peakat",
                               new String[] {"input=" + inputFilename, "output=" + outputFilename,
                                             "numLines=" + numLines, "tailNotHead=" + tail,
                                             "countLines=" + counting,
                                             "columnsNotLines=" + columnsNotLines},
                               log);

    if (params != null) {
      for (int i = 0; i < params.size(); i++) {
        trav = params.elementAt(i);
        if (trav.startsWith("input=")) {
          inputFilename = ext.parseStringArg(trav, null);
        } else if (trav.startsWith("output=")) {
          outputFilename = ext.parseStringArg(trav, null);
        } else if (trav.startsWith("numLines=")) {
          numLines = ext.parseIntArg(trav);
        } else if (trav.startsWith("tailNotHead=")) {
          tail = ext.parseBooleanArg(trav);
        } else if (trav.startsWith("countLines=")) {
          counting = ext.parseBooleanArg(trav);
        } else if (trav.startsWith("columnsNotLines=")) {
          columnsNotLines = ext.parseBooleanArg(trav);
        }
      }

      actualPeak(inputFilename, outputFilename, numLines, tail, counting, columnsNotLines);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "peakAt.dat";

    String usage = "\n" + "widgets.peakAt requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n" + "";

    try {
      for (String arg : args) {
        if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
          System.err.println(usage);
          throw new Elision();
        } else {
          filename = arg;
          numArgs--;
        }
      }
      if (numArgs != 0) {
        System.err.println(usage);
        throw new Elision();
      }

      peak(filename);
    } catch (Exception e) {
      e.printStackTrace();
      ext.waitForResponse();
    }
  }
}
