package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

import com.google.common.primitives.Ints;

public class ConsolidateTaxons {
  public static void consolidate(String rootDir) {
    String[] dirs = Files.listDirectories(rootDir, false);
    Logger log = new Logger();
    ArrayList<Hashtable<String, String[]>> taxa = new ArrayList<Hashtable<String, String[]>>();
    HashSet<String> allTaxa = new HashSet<String>();
    String output = rootDir + "taxaSummary.txt";
    log.reportTimeInfo("NUMDIRS = " + dirs.length);
    for (String dir : dirs) {
      String[] contams = Files.listFullPaths(rootDir + dir + "/", "", false);
      log.reportTimeInfo("Current directory " + rootDir + dir + " Number of files "
                         + contams.length);

      for (String contam : contams) {

        if (contam.endsWith("contam")) {
          log.reportTimeInfo("Consolidating " + contam);
          Hashtable<String, String[]> current = new Hashtable<String, String[]>();
          try {
            BufferedReader reader = Files.getAppropriateReader(contam);
            while (reader.ready()) {
              String[] line = reader.readLine().trim().split("\t");
              allTaxa.add(line[0]);
              current.put(line[0], Array.subArray(line, 1, line.length));

              // System.out.println(Array.toStr(current.get(line[0])));
            }
            reader.close();
            taxa.add(current);
          } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \"" + contam + "\" not found in current directory");
            return;
          } catch (IOException ioe) {
            log.reportError("Error reading file \"" + contam + "\"");
            return;
          }
        }
      }
    }

    try {
      PrintWriter writer = new PrintWriter(new FileWriter(output));
      writer.print("Taxa");
      for (int i = 0; i < taxa.size(); i++) {
        String[] files = taxa.get(i).get("Taxa");
        for (String file : files) {
          writer.print("\t" + ext.rootOf(file));
        }
      }
      writer.println();
      for (String ataxa : allTaxa) {
        if (!ataxa.equals("Taxa")) {
          ArrayList<Integer> counts = new ArrayList<Integer>();
          for (int i = 0; i < taxa.size(); i++) {
            int[] blankCounts = new int[taxa.get(i).get("Taxa").length];
            Arrays.fill(blankCounts, 0);
            String[] blanks = new String[taxa.get(i).get("Taxa").length];
            Arrays.fill(blanks, "0");
            if (taxa.get(i).containsKey(ataxa)) {
              String[] tmp = taxa.get(i).get(ataxa);
              for (String element : tmp) {
                counts.add(Integer.parseInt(element));
              }
            } else {
              for (int blankCount : blankCounts) {
                counts.add(blankCount);

              }
            }
          }
          int[] allCounts = Ints.toArray(counts);
          if (Array.max(allCounts) > 100) {
            writer.println(ataxa + "\t" + Array.toStr(allCounts));
          }
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + output);
      log.reportException(e);
    }
  }

  public static void main(String[] args) {
    // int numArgs = args.length;
    String rootDir = "/home/tsaim/shared/Project_Tsai_Project_021/contamination/";
    consolidate(rootDir);
  }
}
