package org.genvisis.assoc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;

public class MergeDuplicateCalls {
  public static void main(String[] args) {
    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\APOE\\";
    String filename = "snps.ped";

    merge(dir, filename);
  }

  public static void merge(String dir, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] trav, line, keys;
    String id;
    Hashtable<String, String[]> hash;
    int numFields;

    hash = new Hashtable<String, String[]>();
    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      numFields = -1;
      while (reader.ready()) {
        trav = reader.readLine().trim().split("[\\s]+");
        if (trav.length == 1) {
          if (!trav[0].equals("")) {
            System.err.println("Error - what am I supposed to do with: " + trav[0]);
          }
        } else {
          id = trav[0] + "-" + trav[1];
          if (numFields == -1) {
            numFields = trav.length;
          } else if (trav.length != numFields) {
            System.err.println("Error - different number of fields for: " + id);
            System.exit(1);
          }
          if (hash.containsKey(id)) {
            line = hash.get(id);
            for (int i = 2; i < numFields; i++) {
              if (line[i].equals("0") || line[i].equals(".")) {
                line[i] = trav[i];
              } else if (!trav[i].equals("0") && !line[i].equals(trav[i])) {
                System.err.println("Error - mismatch in column " + (i + 1) + " for " + id + ": "
                                   + line[i] + " versus " + trav[i]);
              }
            }
          } else {
            hash.put(id, trav);
          }
        }
      }
      reader.close();

    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + filename + "_marged.xln"));
      keys = HashVec.getKeys(hash);
      for (String key : keys) {
        writer.println(Array.toStr(hash.get(key)));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + filename + "_marged.xln");
      e.printStackTrace();
    }
  }
}
