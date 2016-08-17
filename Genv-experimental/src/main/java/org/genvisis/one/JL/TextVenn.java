package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * @author lane0212 Does a text type venn ish comparison for overlapping keys...
 */
public class TextVenn {

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "D:/data/logan/OSv2_seq/RegNovo/Reg_superComp.txt";
    String logfile = null;
    Logger log;

    String usage = "\n" + "one.JL.Overlapper requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      overlapIt(filename, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static String overlapIt(String inputFile, Logger log) {
    String[] header = Files.getHeaderOfFile(inputFile, log);

    log.reportTimeInfo("Found " + header.length + " categories in " + inputFile);
    log.reportTimeInfo(Array.toStr(header));
    ArrayList<HashSet<String>> cats = new ArrayList<HashSet<String>>();
    for (String element : header) {
      cats.add(new HashSet<String>());
    }
    try {
      BufferedReader reader = Files.getAppropriateReader(inputFile);

      reader.readLine();
      while (reader.ready()) {
        String[] line = reader.readLine().split("\t");
        for (int i = 0; i < line.length; i++) {
          if (!line[i].equals("")) {
            cats.get(i).add(line[i].trim());
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + inputFile + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + inputFile + "\"");
      return null;
    }
    for (int i = 0; i < cats.size(); i++) {
      log.reportTimeInfo("cat " + header[i] + " n= " + cats.get(i).size());
    }

    ArrayList<String> comp = new ArrayList<String>();
    ArrayList<String> compSummary = new ArrayList<String>();
    ArrayList<ArrayList<String>> overlaps = new ArrayList<ArrayList<String>>();
    int maxOverlap = 0;
    for (int i = 0; i < cats.size(); i++) {
      for (int j = i + 1; j < cats.size(); j++) {
        comp.add(header[i] + "(n=" + cats.get(i).size() + ") vs " + header[j] + " (n= "
                 + cats.get(j).size() + ")");
        int overlap = 0;
        ArrayList<String> tmp = new ArrayList<String>();
        for (String hit : cats.get(i)) {
          if (cats.get(j).contains(hit)) {
            overlap++;
            tmp.add(hit);

          }
        }
        if (overlap > maxOverlap) {
          maxOverlap = overlap;
        }
        overlaps.add(tmp);
        compSummary.add(overlap + "");
        System.out.println(overlap);
      }
    }
    String output = ext.addToRoot(inputFile, ".overlapSummary.txt");
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(output));
      writer.println(Array.toStr(comp.toArray(new String[comp.size()])));
      writer.println(Array.toStr(compSummary.toArray(new String[compSummary.size()])));
      for (int i = 0; i < maxOverlap; i++) {
        for (int j = 0; j < overlaps.size(); j++) {
          if (overlaps.get(j).size() > 0) {
            writer.print((j == 0 ? "" : "\t") + overlaps.get(j).remove(0));
          } else {
            writer.print((j == 0 ? "" : "\t") + "");
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + output);
      log.reportException(e);
    }

    return null;
  }

}
