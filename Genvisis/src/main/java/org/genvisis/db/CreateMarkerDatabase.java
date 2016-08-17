// accessed via a .crf file
package org.genvisis.db;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class CreateMarkerDatabase {
  public static final String[] ALT_LOCS =
      {"C:\\Documents and Settings\\npankrat\\My Documents\\1_CRFdb\\genotypes\\"};
  public static final String[] FIELDS =
      {"!UniqueID", "!FamID", "!IndID", "!Source", ".AOO", "!AgeAtExam", ".Dx", ".VPD", "!Affected",
          "!Male", "!Caucasian", ".AffParent", "0parkin", "1noKnownHomozygousParkinMutation",
          "0G2019S", "0AnyLRRK2", "0CausativeMutation", "!Use", "!Comment"};
  public static final String[] DATABASES = {
      "C:\\Documents and Settings\\npankrat\\My Documents\\1_CRFdb\\crf_db.dat",
      "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Global PD files\\control.db.xln"};
  public static final String[] POST_HOCS = {"Age:AOO|AgeAtExam|AgeAtExam|AgeAtExam:Use=1",
      "VPD_Analysis:VPD|1|1|1:Use=1", "ALL_Analysis:Affected|1|1|1:Use=1"};
  public static final String UNIQUEID = "UniqueID";

  public static void create(String outfile, Vector<String> infoV) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, header, dataline, postHoc, requirement, options;
    String trav, markerName;
    Hashtable<String, String[][]> hash = new Hashtable<String, String[][]>();
    String[][] data;
    int[] indices;
    String[] fields;
    boolean[] required;
    String[] defaults;
    boolean flag;
    int index;
    String[][] info;

    info = new String[infoV.size()][];
    for (int i = 0; i < info.length; i++) {
      info[i] = infoV.elementAt(i).trim().split("[\\s]+");
      markerName = info[i][0];
      try {
        reader = Files.getReader(markerName, ALT_LOCS);
        header = reader.readLine().trim().split(",");
        if (header.length == 1) {
          System.err
              .println("Error - CreateMarkerDatabase currently requires comma delimited files");
          System.exit(1);
        }
        if (info[i].length == 1) {
          System.out.println("Using all columns for " + markerName + " since none were specified");
          index = ext.indexOfStr(UNIQUEID, header);
          if (index == -1) {
            System.err.println(
                "Error - " + UNIQUEID + " was not found as a column header in " + markerName);
          }
          info[i] = new String[header.length + 1];
          info[i][0] = markerName;
          info[i][1] = UNIQUEID;
          for (int j = 2; j < info[i].length; j++) {
            info[i][j] = j - 2 < index ? header[j - 2] : header[j - 1];
          }
        }
        indices = new int[info[i].length - 1];
        for (int j = 0; j < indices.length; j++) {
          trav = info[i][j + 1].contains("|")
              ? info[i][j + 1].substring(0, info[i][j + 1].indexOf("|")) : info[i][j + 1];
          indices[j] = ext.indexOfStr(trav, header);
          if (indices[j] == -1) {
            System.err
                .println("Error - '" + trav + "' is not a column header in file " + markerName);
            System.exit(1);
          }
        }
        while (reader.ready()) {
          line = reader.readLine().trim().split(",");
          if (line.length == 1) {
            System.err.println("Error - expecting " + markerName + " to be a comma delimited file");
            System.exit(1);
          }
          dataline = new String[indices.length - 1];
          for (int j = 0; j < dataline.length; j++) {
            dataline[j] = line[indices[j + 1]];
          }
          if (hash.containsKey(line[indices[0]])) {
            data = hash.get(line[indices[0]]);
          } else {
            hash.put(line[indices[0]], data = new String[info.length][]);
          }
          if (data[i] != null) {
            System.err.println(
                "Error - duplicate data in file '" + markerName + "' for " + line[indices[0]]);
          }
          data[i] = dataline;
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println(fnfe.getMessage());
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + markerName + "\"");
        System.exit(2);
      }
    }

    fields = new String[FIELDS.length];
    required = new boolean[FIELDS.length];
    defaults = new String[FIELDS.length];
    for (int j = 0; j < FIELDS.length; j++) {
      fields[j] = FIELDS[j].substring(1);
      defaults[j] = FIELDS[j].substring(0, 1);
      required[j] = defaults[j].equals("!");
    }

    try {
      writer = new PrintWriter(new FileWriter(outfile));
      writer.print("UniqueID");
      for (String[] element : info) {
        markerName = element[0].substring(0, element[0].lastIndexOf("."));
        for (int j = 2; j < element.length; j++) {
          trav =
              element[j].contains("|") ? element[j].substring(element[j].indexOf("|")) : element[j];
          writer.print("\t" + (info.length == 1 ? "" : markerName + "_") + trav);
        }
      }
      for (int j = 1; j < fields.length; j++) {
        writer.print("\t" + fields[j]);
      }
      for (String element : POST_HOCS) {
        writer.print("\t" + element.split(":")[0]);
      }
      writer.println();

      for (int i = 0; i < DATABASES.length; i++) {
        try {
          reader = new BufferedReader(new FileReader(DATABASES[i]));
          header = reader.readLine().trim().split("\\t", -1);
          indices = new int[fields.length];
          flag = false;
          for (int j = 0; j < fields.length; j++) {
            indices[j] = ext.indexOfStr(fields[j], header);
            if (required[j] && indices[j] == -1) {
              System.err.println("Error - '" + fields[j] + "' is not a column header in file "
                  + DATABASES[i] + " but is declared as being required");
              flag = true;
            }
          }
          if (flag) {
            System.exit(1);
          }
          while (reader.ready()) {
            line = reader.readLine().trim().split("\\t", -1);
            if (hash.containsKey(line[indices[0]])) {
              data = hash.get(line[indices[0]]);
              if (data.length == 0) {
                System.err.println("Error - multiple refrences to '" + line[indices[0]] + "' in "
                    + DATABASES[i] + " and/or earlier databases");
              } else {
                hash.put(line[indices[0]], new String[0][0]);
                flag = false;
                writer.print(line[indices[0]]);
                for (String[] element : data) {
                  if (element == null) {
                    flag = true;
                    writer.print("\t" + Array.toStr(Array.stringArray(info[i].length - 2, ".")));
                  } else {
                    writer.print("\t" + Array.toStr(element));
                  }
                }
                if (flag) {
                  System.err.println("Warning - '" + line[indices[0]]
                      + "' is missing data for at least one of the markers.");
                }
                for (int j = 1; j < FIELDS.length; j++) {
                  writer.print("\t" + (indices[j] == -1
                      || (line[indices[j]].equals(".") && !FIELDS[j].substring(0, 1).equals("!"))
                          ? FIELDS[j].substring(0, 1) : line[indices[j]]));
                }
                for (String element : POST_HOCS) {
                  postHoc = element.trim().split(":");
                  requirement = postHoc[2].trim().split("=");
                  if (ext.indexOfStr(requirement[0], header) == -1) {
                    System.err
                        .println("Error - '" + requirement[0] + "' is a requirement for posthoc "
                            + postHoc[0] + ", but is not in header of " + DATABASES[i]);
                    System.exit(1);
                  }
                  if (!line[ext.indexOfStr(requirement[0], header)].equals(requirement[1])) {
                    writer.print("\t.");
                  } else {
                    options = postHoc[1].trim().split("\\|");
                    if (ext.indexOfStr(options[i], header) == -1) {
                      writer.print("\t" + options[i]);
                    } else {
                      writer.print("\t" + line[ext.indexOfStr(options[i], header)]);
                    }
                  }
                }
                writer.println();
              }
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + DATABASES[i] + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + DATABASES[i] + "\"");
          System.exit(2);
        }
      }

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + outfile);
      e.printStackTrace();
    }
  }

  public static void createFromParameters(String filename, Logger log) {
    Vector<String> paramV;

    paramV = Files.parseControlFile(filename, "db",
        new String[] {"file1.csv ColumnHeaderForUniqueID ColumnHeader1 ColumnHeader3 ColumnHeader5",
            "file2.csv ColumnHeaderForUniqueID ColumnHeader2 ColumnHeader3	ColumnHeader9",
            "file3.csv ColumnHeaderForUniqueID ColumnHeader9", "..."},
        log);
    if (paramV != null) {
      create(ext.rootOf(filename) + ".xln", paramV);
    }
  }
}
