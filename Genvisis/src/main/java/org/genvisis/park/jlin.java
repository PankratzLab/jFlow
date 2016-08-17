package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

public class jlin {
  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "progeni_haploview.ped";
    // String filename = "negneuros.ped";
    String mapfile = "map.info";

    String usage = "\n" + "park.jlin requires 2 arguments:\n" + "   (1) a pre file (i.e. file="
                   + filename + " (default))\n" + "   (2) a Haploview-style map file (i.e. map="
                   + mapfile + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("map=")) {
        mapfile = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (args.length == 0) {
      System.err.println("Using defaults (file=" + filename + " map=" + mapfile + ")");
    }

    try {
      new jlin(filename, mapfile);
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

  public String[][] ALLELE_LOOKUP = {{"", "A", "G"}, {"", "A", "G"}, {"", "A", "G"}};

  public jlin(String filename, String mapfile) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    Hashtable<String, Vector<int[]>> hash = new Hashtable<String, Vector<int[]>>();
    Vector<String> fams = new Vector<String>();
    Vector<int[]> v;
    String temp;
    int[] data;
    int score;

    if (!new File(filename).exists()) {
      System.err.println("Error - could not find " + filename + " in current directory");
      System.exit(2);
    }
    reader = new BufferedReader(new FileReader(filename));
    reader.readLine();
    while (reader.ready()) {
      line = reader.readLine().split("[\\s]+");
      data = new int[line.length - 6];
      if (data.length % 2 != 0) {
        System.err.println("Error - Expecting 6 columns (including affection status) followed by an even number of alleles");
        System.exit(4);
      }
      score = 0;
      for (int i = 0; i < line.length - 6; i++) {
        data[i] = Integer.valueOf(line[i + 6]).intValue();
        if (data[i] != 0) {
          score++;
        }
      }
      if (score > 1 && line[5].equals("2")) {
        // if (score == 6 && line[5].equals("2")) {
        if (hash.containsKey(line[0])) {
          v = hash.get(line[0]);
        } else {
          hash.put(line[0], v = new Vector<int[]>());
          fams.add(line[0]);
        }
        v.add(data);
      }
    }
    reader.close();

    writer = new PrintWriter(new FileWriter(filename + "-jlin.csv"));
    if (!new File(mapfile).exists()) {
      System.err.println("Error - could not find " + mapfile + " in current directory");
      System.exit(2);
    }
    reader = new BufferedReader(new FileReader(mapfile));
    temp = "";
    score = 0;
    while (reader.ready()) {
      line = reader.readLine().split("[\\s]+");
      writer.print((score == 0 ? "" : ",") + line[0]);
      temp += (score == 0 ? "" : ",") + line[1];
      score++;
    }
    writer.println();
    writer.println(temp);
    reader.close();

    for (int i = 0; i < fams.size(); i++) {
      v = hash.get(fams.elementAt(i));
      // writer.println(translate(v.elementAt(0), ALLELE_LOOKUP));
      for (int j = 0; j < v.size(); j++) {
        writer.println(translate(v.elementAt(j), ALLELE_LOOKUP));
      }
    }

    writer.close();
  }

  public String translate(int[] data, String[][] lookup) throws IOException {
    String str = "";
    for (int i = 0; i < data.length / 2; i++) {
      str += (i == 0 ? "" : ",") + lookup[i][data[i * 2]] + lookup[i][data[i * 2 + 1]];
    }
    return str;
  }
}
