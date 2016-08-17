package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.ext;

public class prioritize {
  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "anyExtras.xls";
    String database = tools.CRF_DIR + "crf_db.dat";
    String variables = "FamID:^G2019S;>OtherLRRK2;^parkin;^VPD";

    String usage = "\n" + "park.prioritize requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new prioritize(filename, database, variables);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public prioritize(String filename, String databases, String variables) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, l2;
    String temp, trav;
    Hashtable<String, String[][]> hash = new Hashtable<String, String[][]>();
    String[] DBs = databases.split("\\|");
    String lookup = "";
    int[] lookupIndices = new int[DBs.length + 1];
    int[][] variableIndices = new int[DBs.length][];
    char[][] variableOperators = new char[DBs.length][];
    String[][] data;

    for (int db = 0; db < DBs.length; db++) {
      try {
        reader = new BufferedReader(new FileReader(DBs[db]));
        line = reader.readLine().split("\t", -1);
        if (lookup.equals("")) {
          lookup = (variables.split("\\|")[db]).split(":")[0];
        } else if (!lookup.equals((variables.split("\\|")[db]).split(":")[0])) {
          System.err.println("Error - lookup values need to be the same in all files and databases");
          System.exit(1);
        }
        lookupIndices[db] = ext.indexOfStr(lookup, line);
        if (lookupIndices[db] == -1) {
          System.err.println("Error - lookup variable '" + lookup + "' was not found in database '"
                             + DBs[db] + "'");
          System.exit(1);
        }
        l2 = (variables.split("\\|")[db]).split(":")[1].split(";");
        variableOperators[db] = new char[l2.length];
        variableIndices[db] = new int[l2.length];
        for (int i = 0; i < l2.length; i++) {
          variableOperators[db][i] = l2[i].substring(0, 1).charAt(0);
          variableIndices[db][i] = ext.indexOfStr(l2[i].substring(1), line);
          if (variableIndices[db][i] == -1) {
            System.err.println("Error - variable '" + l2[i].substring(1)
                               + "' was not found in database '" + DBs[db]
                               + "' (do you need an operator)");
            System.exit(1);
          }
        }
        while (reader.ready()) {
          line = reader.readLine().split("\t", -1);
          trav = line[lookupIndices[db]];
          if (hash.containsKey(trav)) {
            data = hash.get(trav);
          } else {
            data = new String[DBs.length][];
            hash.put(trav, data);
          }
          if (data[db] == null) {
            data[db] = Array.stringArray(variableIndices[db].length);
          }
          for (int i = 0; i < data[db].length; i++) {
            temp = line[variableIndices[db][i]];
            switch (variableOperators[db][i]) {
              case '^':
                if (temp.equals(".")) {
                  temp = "-1";
                }
                try {
                  if (data[db][i].equals("")
                      || Integer.parseInt(temp) > Integer.parseInt(data[db][i])) {
                    data[db][i] = temp;
                  }
                } catch (NumberFormatException nfe) {
                  System.err.println("Error - couldn't parse variable: '" + temp + "'");
                  data[db][i] = "999";
                }
                break;
              case '=':
                data[db][i] = temp;
                break;
              case '+':
                if (data[db][i].equals("")) {
                  data[db][i] = "0";
                }
                if (temp.equals(".")) {
                  temp = "0";
                }
                data[db][i] = "" + (Integer.parseInt(temp) + Integer.parseInt(data[db][i]));
                break;
              case '>':
                if (temp.equals(".")) {
                  temp = "";
                }
                if (temp.length() > data[db][i].length()) {
                  data[db][i] = temp;
                }
                break;
              default:
                System.err.println("Error - unknown operator '" + variableOperators[db][i] + "'");
                System.exit(1);
            }
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error - could not find " + DBs[db] + " in current directory");
        System.exit(2);
      } catch (IOException ioe) {
        System.err.println("Error parsing " + DBs[db] + "");
        System.exit(3);
      }
    }

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(filename + "-meshed.xls"));
      line = reader.readLine().split("\t", -1);
      lookupIndices[DBs.length] = ext.indexOfStr(lookup, line);
      if (lookupIndices[DBs.length] == -1) {
        System.err.println("Error - lookup value for the databases ('" + lookup
                           + "') does not match a column in the file " + filename);
        System.exit(1);
      }
      writer.print(Array.toStr(line));
      for (int i = 0; i < DBs.length; i++) {
        line = (variables.split("\\|")[i]).split(":")[1].split(";");
        for (int j = 0; j < line.length; j++) {
          line[j] = line[j].substring(1);
        }
        writer.print("\t" + Array.toStr(line));
      }
      writer.println();
      while (reader.ready()) {
        line = reader.readLine().split("\t", -1);
        writer.print(Array.toStr(line));
        data = hash.get(line[lookupIndices[DBs.length]]);
        if (data == null) {
          System.err.println("Error - no data for '" + line[lookupIndices[DBs.length]] + "'");
        }
        for (int i = 0; i < DBs.length; i++) {
          for (int j = 0; j < variableIndices[i].length; j++) {
            writer.print("\t" + (data == null ? "XXX"
                                              : (data[i][j].equals("-1")
                                                 || data[i][j].equals("") ? "." : data[i][j])));
          }
        }
        writer.println();
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - could not find " + filename + " in current directory");
      System.exit(2);
    } catch (IOException ioe) {
      System.err.println("Error parsing " + filename + "");
      System.exit(3);
    }

  }
}
