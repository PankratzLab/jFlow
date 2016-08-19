package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.AlleleFreq;
import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class extractMarkers {
  public static final String BASE_DIR = "extractedMarkers";

  public extractMarkers(String filename) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    String temp, chrome;
    Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
    int count;
    String[] chrs, markerNames, markersPicked;
    int[] indices;
    double sum;
    int sigfigs;
    int[][] genotypeCounts;
    String dir;

    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        if (line.length != 2) {
          System.err.println("Error - expecting 2 columns: chr#\tmarkerName");
          System.exit(1);
        }
        try {
          if (Integer.parseInt(line[0]) < 1 || Integer.parseInt(line[0]) > 25) {
            System.err.println("Error - expecting chromosome to be between 1 and 25");
            System.exit(1);
          }
        } catch (Exception e) {
          System.err.println("Error - expecting 2 columns: chr#\tmarkerName");
          System.exit(1);
        }
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

    count = 0;
    while (new File(BASE_DIR + "." + (++count)).exists()) {
      ;
    }
    dir = BASE_DIR + "." + count + "/";
    new File(dir).mkdir();

    chrs = HashVec.getKeys(hash);
    for (String chr : chrs) {
      chrome = ext.chrome(Integer.parseInt(chr));
      markerNames = new LinkageMap("map" + chrome + ".dat").getMarkerNames();
      markersPicked = Array.toStringArray(hash.get(chr));
      indices = ext.indexFactors(markersPicked, markerNames, false, false);

      if (Array.min(indices) < 0) {
        System.err.println("Error - the following markers were not found in " + "map" + chrome
                           + ".dat");
        for (int j = 0; j < indices.length; j++) {
          if (indices[j] < 0) {
            System.err.print(markersPicked[j] + " ");
          }
        }
        System.err.println();
        System.exit(1);
      }

      try {
        reader = new BufferedReader(new FileReader("map" + chrome + ".dat"));
        writer = new PrintWriter(new FileWriter(dir + "map" + chrome + ".dat"));

        temp = reader.readLine().trim();
        writer.println((markersPicked.length + 1)
                       + temp.substring(temp.split("[\\s]+")[0].length()));
        writer.println(reader.readLine());
        reader.readLine();
        writer.println(Array.toStr(Array.stringArraySequence(markersPicked.length + 1, ""), " "));
        for (int j = 0; j < (chrome.equals("23") ? 5 : 4); j++) {
          writer.println(reader.readLine());
        }

        for (int j = 0; j < markerNames.length; j++) {
          if (ext.indexOfInt(j, indices) != -1) {
            writer.println(reader.readLine());
            writer.println(reader.readLine());
          } else {
            reader.readLine();
            reader.readLine();
          }
        }
        writer.println(reader.readLine());

        line = reader.readLine().split("[\\s]+");
        indices = Sort.putInOrder(indices);
        sum = Double.parseDouble(line[0]) * -1;
        sigfigs = 0;
        for (int j = 0; j < line.length - 3; j++) {
          sum += Double.parseDouble(line[j]);
          if (line[j].indexOf(".") >= 0
              && line[j].substring(line[j].indexOf(".") + 1).length() > sigfigs) {
            sigfigs = line[j].substring(line[j].indexOf(".") + 1).length();
          }
          if (ext.indexOfInt(j, indices) != -1) {
            if (j == indices[0]) {
              writer.print(line[0]);
            } else {
              writer.print(" " + ext.formDeci(sum, sigfigs, true));
            }
            sum = 0;
          }
        }
        writer.println("  " + line[line.length - 3] + " " + line[line.length - 2] + " "
                       + line[line.length - 1]);
        writer.println(reader.readLine());

        reader.close();
        writer.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + "map" + chrome + ".dat"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + "map" + chrome + ".dat" + "\"");
        System.exit(2);
      }

      try {
        reader = new BufferedReader(new FileReader("re_chrom" + chrome + ".pre"));
        writer = new PrintWriter(new FileWriter(dir + "re_chrom" + chrome + ".pre"));
        genotypeCounts = new int[indices.length][3];
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          writer.print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4]
                       + "\t" + line[5]);
          for (int j = 0; j < markerNames.length; j++) {
            if (ext.indexOfInt(j, indices) != -1) {
              writer.print("\t" + line[6 + 2 * j + 0] + "\t" + line[6 + 2 * j + 1]);
              count = Integer.parseInt(line[6 + 2 * j + 0]) + Integer.parseInt(line[6 + 2 * j + 1])
                      - 2;
              if (count >= 0) {
                genotypeCounts[ext.indexOfInt(j, indices)][count]++;
              }
            }
          }
          writer.println();
        }
        writer.close();
        reader.close();

        writer = new PrintWriter(new FileWriter(dir + "heterozygosity" + chrome + ".xls"));
        for (int j = 0; j < indices.length; j++) {
          writer.println(markerNames[indices[j]] + "\t"
                         + ext.formDeci(AlleleFreq.computeHeterozygosity(genotypeCounts[j]), 3));
          if (Math.abs(AlleleFreq.computeHeterozygosity(genotypeCounts[j])) < 0.00001) {
            System.out.println(markerNames[indices[j]] + " is uninformative");
          }
        }
        writer.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + "re_chrom" + chrome + ".pre"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + "re_chrom" + chrome + ".pre" + "\"");
        System.exit(2);
      } catch (Exception e) {
        e.printStackTrace();
        System.exit(2);
      }
      Cyrillic.format(dir + "re_chrom" + chrome + ".pre", dir + "map" + chrome + ".dat",
                      dir + "cyrill" + chrome + ".pre", dir + "cyrill" + chrome + ".dat",
                      Integer.parseInt(chr));
    }
    System.out.println("Files with the extracted markers were successfully created in '" + dir
                       + "'");
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "extractThese.txt";

    String usage = "\n" + "park.extractMarkers requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

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
      new extractMarkers(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
