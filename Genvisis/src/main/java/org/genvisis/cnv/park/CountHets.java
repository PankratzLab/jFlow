package org.genvisis.cnv.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.ext;

class CNVregion {
  private final String sample;

  private final String famID;

  private final String indID;

  private final String chr;

  private final int start;

  private final int stop;

  private final int type;

  private int countHet;

  private int countCalled;

  private int baf15;

  private int baf25;

  private int countValidBAF;

  private int countMarkers;

  public CNVregion(String[] line) {
    sample = line[0];
    famID = line[1];
    indID = line[2];
    chr = line[3];
    start = Integer.parseInt(line[4]);
    stop = Integer.parseInt(line[5]);
    type = Integer.parseInt(line[6]);
  }

  public void checkPosition(String travChr, int travPosition, String call, String baf) {
    if (travChr.equals(chr) && travPosition >= start && travPosition <= stop) {
      if (call.equals("AB")) {
        countHet++;
      }
      if (call.equals("AA") || call.equals("AB") || call.equals("BB")) {
        countCalled++;
      }
      try {
        if (Double.parseDouble(baf) >= 0.15 && Double.parseDouble(baf) <= 0.85) {
          baf15++;
        }
        if (Double.parseDouble(baf) >= 0.25 && Double.parseDouble(baf) <= 0.75) {
          baf25++;
        }
        countValidBAF++;
      } catch (NumberFormatException nfe) {
      }
      countMarkers++;
    }
  }

  public String getSampleID() {
    return sample;
  }

  @Override
  public String toString() {
    return sample + "\t" + famID + "\t" + indID + "\t" + chr + "\t" + start + "\t" + stop + "\t"
        + type + "\t" + countMarkers + "\t"
        + ext.formDeci((double) countHet / (double) countCalled, 3) + "\t"
        + ext.formDeci((double) baf15 / (double) countValidBAF, 3) + "\t"
        + ext.formDeci((double) baf25 / (double) countValidBAF, 3); // countHet+"\t"+countCalled+"\t"+
  }
}


public class CountHets {
  public static final String WINDOWS_DIRECTORY =
      "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\allMarkers\\chr14\\";
  public static final String LINUX_DIRECTORY = "/home/npankrat/penncnv/allMarkers/cnvs/";
  public static final String[] EXPECTED_HEADER =
      {"Sample", "FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"};
  public static final String[] CNV_SUFFIXES =
      {"Name", "Chr", "Position", ".GType", ".Log R Ratio", ".B Allele Freq"};

  public static void count(String dir, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String trav;
    Hashtable<String, Vector<CNVregion>> hash = new Hashtable<String, Vector<CNVregion>>();
    Vector<String> hashKeys = new Vector<String>();
    Vector<CNVregion> v = new Vector<CNVregion>();
    CNVregion region;
    String cnvDirectory;
    boolean problem;

    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      ext.checkHeader(reader.readLine().trim().split("\t"), EXPECTED_HEADER, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length != EXPECTED_HEADER.length) {
          System.err.println("Error - mismatched number of columns in " + dir + filename);
          System.exit(1);
        }
        region = new CNVregion(line);
        trav = region.getSampleID();
        if (hash.containsKey(trav)) {
          v = hash.get(trav);
        } else {
          hash.put(trav, v = new Vector<CNVregion>());
          hashKeys.add(trav);
        }
        v.add(region);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    if (new File(WINDOWS_DIRECTORY).exists()) {
      cnvDirectory = WINDOWS_DIRECTORY;
    } else if (new File(LINUX_DIRECTORY).exists()) {
      cnvDirectory = LINUX_DIRECTORY;
    } else {
      cnvDirectory = "";
      System.err.println("Error - could not resolve directory to find raw cnv data");
      System.err.println(WINDOWS_DIRECTORY);
      System.err.println(LINUX_DIRECTORY);
      System.exit(1);
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + filename + "_hetCounts.xln"));
      writer.println(
          Array.toStr(EXPECTED_HEADER) + "\t#Markers\t%CalledHet\t%0.15<BAF,0.85\t%0.25<BAF,0.75");
      for (int i = 0; i < hashKeys.size(); i++) {
        v = hash.get(hashKeys.elementAt(i));
        try {
          reader = new BufferedReader(new FileReader(cnvDirectory + hashKeys.elementAt(i)));
          problem = false;
          line = reader.readLine().trim().split("\t");
          for (int j = 0; j < CNV_SUFFIXES.length; j++) {
            if (!line[j].endsWith(CNV_SUFFIXES[j])) {
              System.err.println("Error - expecting column " + (j + 1) + " to end in "
                  + CNV_SUFFIXES[j] + " (found " + line[j] + ")");
              problem = true;
            }
          }
          if (problem) {
            System.exit(1);
          }
          while (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
            for (int j = 0; j < v.size(); j++) {
              v.elementAt(j).checkPosition(line[1], Integer.parseInt(line[2]), line[3], line[5]);
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + cnvDirectory + hashKeys.elementAt(i)
              + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + cnvDirectory + hashKeys.elementAt(i) + "\"");
          System.exit(2);
        }
        for (int j = 0; j < v.size(); j++) {
          writer.println(v.elementAt(j).toString());
        }
        writer.flush();
      }

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to '" + dir + filename + "_hetCounts.xln" + "'");
      e.printStackTrace();
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\allMarkers\\";
    String filename = "chr14.cnv";

    String usage = "\\n" + "park.cnv.CountHets requires 0-1 arguments\n"
        + "   (1) directory (i.e. dir=" + dir + " (default; use dir=./ for pwd))\n"
        + "   (2) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
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
      count(dir, filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
