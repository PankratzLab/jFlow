package org.genvisis.dead;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.park.*;

public class AgeAtDeath {
  public static final String[] HEADER =
      {"FamNo", "IndNo", "Mother", "MOB", "DOB", "YOB", "MOD", "DOD", "YOD", "Dead", "Father",
       "MOB", "DOB", "YOB", "MOD", "DOD", "YOD", "Dead"};
  public static final String DEFAULT_FILE =
      "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\mtDNA\\matchingParentalAAD\\parents_age_at_death-020210.txt";
  public static final double DIFF_AGE_THRESHOLD = 10.0;
  public static final int[] FLIP = {1, 0};

  public static void parse(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> dx;
    double[] ages;
    String[] dates;
    boolean use;
    String fatherID, motherID;
    int affected;

    dx = tools.getBestPDdx();

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(filename + ".out"));
      writer.println("FID\tIID\tMother_Affected\tMother_AAD\tFather_Affected\tFather_AAD\tUse");
      line = reader.readLine().trim().split("[\\s]+");
      ext.checkHeader(line, new String[] {"Mother", "Father"}, false);
      line = reader.readLine().trim().split("[\\s]+");
      ext.checkHeader(line, HEADER, false);
      while (reader.ready()) {
        line = reader.readLine().split("\\t", -1);
        ages = new double[2];
        use = true;
        for (int i = 0; i < ages.length; i++) {
          if (line[2 + i * 8 + 3].equals("")
              || (line[2 + i * 8 + 6].equals("") && line[2 + i * 8 + 7].equals("TRUE"))) {
            ages[i] = -1;
          } else {
            dates = new String[line[2 + i * 8 + 7].equals("TRUE") ? 2 : 1];
            for (int j = 0; j < dates.length; j++) {
              dates[j] =
                  (line[2 + i * 8 + 1 + j * 3 + 0].equals("") ? "6"
                                                              : line[2 + i * 8 + 1 + j * 3 + 0])
                         + "/"
                         + (line[2 + i * 8 + 1 + j * 3 + 1].equals("") ? "15"
                                                                       : line[2 + i * 8 + 1 + j * 3
                                                                              + 1])
                         + "/" + line[2 + i * 8 + 1 + j * 3 + 2];
            }
            ages[i] =
                ext.calcDays(ext.parseDate(dates[0]),
                             dates.length == 1 ? new Date() : ext.parseDate(dates[1]))
                      / 365.25;
          }
        }
        motherID = tools.getUniqueID(line[0], line[2]);
        fatherID = tools.getUniqueID(line[0], line[10]);
        if (tools.isAffected(dx, motherID) && tools.isAffected(dx, fatherID)) {
          System.err.println("Error - both parents are affected in family " + line[0]);
          affected = -1;
        } else if (tools.isAffected(dx, motherID)) {
          affected = 0;
        } else if (tools.isAffected(dx, fatherID)) {
          affected = 1;
        } else {
          System.err.println("Error - neither parent is affected in family " + line[0]);
          affected = -1;
        }

        if (affected < 0 || Array.min(ages) < 0
            || ages[affected] - ages[FLIP[affected]] > DIFF_AGE_THRESHOLD) {
          use = false;
        }
        writer.println(line[0] + "\t" + line[1] + "\t" + (tools.isAffected(dx, motherID) ? 1 : 0)
                       + "\t" + ages[0] + "\t" + (tools.isAffected(dx, fatherID) ? 1 : 0) + "\t"
                       + ages[1] + "\t" + (use ? 1 : 0));
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
    String filename = DEFAULT_FILE;

    String usage = "\n" + "dead.AgeAtDeath requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n" + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("file=")) {
        filename = args[i].split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      parse(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
