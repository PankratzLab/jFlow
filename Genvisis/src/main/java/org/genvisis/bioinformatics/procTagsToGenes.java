package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;

public class procTagsToGenes {
  public static final String[] DEFAULT_FILES =
      {"knownToGNfAtlas2.prn", "knownToGnf1h.prn", "knownToU133.prn", "knownToU133Plus2.prn",
       "knownToU95.prn"};

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    Vector<String> v = null;

    String usage = "\n" + "park.procTagsToGenes requires 0+ arguments\n" + "   filenames (i.e. "
                   + Array.toStr(DEFAULT_FILES, " ") + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else {
        if (v == null) {
          v = new Vector<String>();
        }
        v.add(arg);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      runProcTagsToGenes(v == null ? DEFAULT_FILES : Array.toStringArray(v));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void runProcTagsToGenes(String[] filenames) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, genes;
    Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
    Vector<String> v;

    for (String filename : filenames) {
      try {
        reader = new BufferedReader(new FileReader(filename));
        reader.readLine();
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          HashVec.addToHashVec(hash, line[0], line[1], true);
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + filename + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + filename + "\"");
        System.exit(2);
      }
    }
    try {
      genes = HashVec.getKeys(hash);
      writer = new PrintWriter(new FileWriter("knownToAll.prn"));
      writer.println("#name\tvalue");
      for (String gene : genes) {
        v = hash.get(gene);
        for (int j = 0; j < v.size(); j++) {
          writer.println(gene + "\t" + v.elementAt(j));
        }
      }
      writer.close();
    } catch (IOException ioe) {
      System.err.println("Error writing " + "knownToAll.prn");
      ioe.printStackTrace();
    }

  }
}
