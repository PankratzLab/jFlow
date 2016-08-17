package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.ByteVector;
import org.genvisis.common.CountVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

import com.google.common.primitives.Bytes;

public class Relpair {
  public static final int MAX_NUM_DIGITS_IN_ALLELE = 2;

  public static final String[][] FILTER_PAIRS =
      {{"UN", "CO"}, {"UN", "AV"}, {"CO", "UN"}, {"CO", "AV"}};

  public static void createFiles(String dir) {
    BufferedReader[] readers;
    PrintWriter writer = null;
    String[] line;
    int maxfamIDsize = 1, maxindIDsize = 1;
    CountVector cv = new CountVector();
    String[] fams;
    int[] indsPerFam;
    ByteVector found, missing;
    byte[] chrs;
    int[] numMarkers;
    LinkageMap map;
    String[] markerNames;
    double[][] alleleFreqs;
    double[] cumulativePositions;

    found = new ByteVector();
    missing = new ByteVector();
    for (byte chr = 1; chr <= 23; chr++) {
      if (new File(dir + "re_chrom" + ext.chrome(chr) + ".pre").exists()) {
        found.add(chr);
      } else {
        missing.add(chr);
      }
    }
    chrs = Bytes.toArray(found);
    if (missing.size() > 0) {
      System.err.println("Warning - Files not found for chromosome"
          + (missing.size() > 1 ? "s" : "") + ": " + Array.toStr(missing.toArray(), " "));
    }
    numMarkers = new int[chrs.length];
    readers = new BufferedReader[chrs.length];
    for (int i = 0; i < chrs.length; i++) {
      try {
        readers[i] =
            new BufferedReader(new FileReader(dir + "re_chrom" + ext.chrome(chrs[i]) + ".pre"));
      } catch (FileNotFoundException fnfe) {
        fnfe.printStackTrace();
      }
    }
    try {
      while (readers[0].ready()) {
        line = readers[0].readLine().trim().split("[\\s]+");
        if (line.length > 1) {
          maxfamIDsize = Math.max(maxfamIDsize, line[0].length());
          maxindIDsize = Math.max(maxindIDsize, line[1].length());
          cv.add(line[0]);
        }
      }
      readers[0].close();
      readers[0] =
          new BufferedReader(new FileReader(dir + "re_chrom" + ext.chrome(chrs[0]) + ".pre"));
    } catch (IOException ioe) {
      System.err
          .println("Error parsing first file: " + dir + "re_chrom" + ext.chrome(chrs[0]) + ".pre");
      ioe.printStackTrace();
    }
    fams = cv.getValues();
    indsPerFam = cv.getCounts();

    try {
      writer = new PrintWriter(new FileWriter(dir + "relpair.loc"));
      for (int chrI = 0; chrI < chrs.length; chrI++) {
        try {
          map = new LinkageMap(dir, chrs[chrI]);
          markerNames = map.getMarkerNames();
          alleleFreqs = map.getAlleleFreqs();
          cumulativePositions = map.getCumulativePositions(true);
          for (int i = 0; i < markerNames.length; i++) {
            writer.println(ext
                .formStr((markerNames[i].length() < 10 ? markerNames[i]
                    : markerNames[i].substring(0, 8) + "-"), 10, true)
                + ((chrs[chrI] < 23) ? "AUTOSOME" : "X-LINKED") + " "
                + ext.formStr(alleleFreqs[i].length + "", 2, true)
                + ext.formStr(chrs[chrI] + "", 4, false) + "    "
                + ext.formDeci(cumulativePositions[i], 3, true));
            for (int j = 0; j < alleleFreqs[i].length; j++) {
              writer.println(ext.formStr((j + 1) + "", 8, true)
                  + ext.formStr(ext.formDeci(alleleFreqs[i][j], 6, true), 8, true));
            }
          }
          numMarkers[chrI] = markerNames.length;
        } catch (Exception e) {
          System.err.println("Error parsing map for chromosome " + chrs[chrI]);
          e.printStackTrace();
        }
      }
      writer.close();
    } catch (IOException ex) {
      System.err.println("Error: could not write to " + dir + "relpair.loc");
      System.exit(1);
    }

    System.out.println("Using " + Array.sum(numMarkers) + " markers.");

    try {
      writer = new PrintWriter(new FileWriter(dir + "relpair.ped"));
      System.out.print("Processing families");

      writer.println("(I2,1X,A" + maxfamIDsize + ")");
      writer.println("(3A" + (maxindIDsize + 1) + ",2A1,A" + (MAX_NUM_DIGITS_IN_ALLELE * 2 + 1)
          + "," + Array.sum(numMarkers) + "(1X,A" + (MAX_NUM_DIGITS_IN_ALLELE * 2 + 1) + "))");

      for (int fam = 0; fam < indsPerFam.length; fam++) {
        writer.println(ext.formStr(indsPerFam[fam] + "", 2, false) + " " + fams[fam]);
        for (int ind = 0; ind < indsPerFam[fam]; ind++) {
          for (int i = 0; i < chrs.length; i++) {
            try {
              line = readers[i].readLine().trim().split("[\\s]+");
              if (line.length > 1) {
                if (i == 0) {
                  writer.print(ext.formStr(line[1], maxindIDsize + 1, true));
                  if (line[2].equals("0") || line[3].equals("0")) {
                    line[2] = " ";
                    line[3] = " ";
                  }
                  writer.print(ext.formStr(line[2], maxindIDsize + 1, true)
                      + ext.formStr(line[3], maxindIDsize + 1, true));
                  writer.print((line[4].equals("1")) ? "M" : "F");
                }

                for (int j = 0; j < numMarkers[i]; j++) {
                  if (line.length < 6 + j * 2 + 2) {
                    System.err.println("Error - what's up with chromosome " + chrs[i] + " for "
                        + line[0] + "," + line[1]);
                  }
                  if (line[6 + j * 2 + 0].equals("0") || line[6 + j * 2 + 1].equals("0")) {
                    writer.print("      ");
                  } else {
                    writer.print(" " + ext.formStr(line[6 + j * 2 + 0] + "/" + line[6 + j * 2 + 1],
                        MAX_NUM_DIGITS_IN_ALLELE * 2 + 1, true));
                  }
                }
              }
            } catch (Exception ex2) {
              System.err
                  .println("Error in processing the re_chrom" + ext.chrome(chrs[i]) + ".pre files");
              ex2.printStackTrace();
            }
          }
          writer.println();
        }
        System.out.print(".");
      }
      writer.close();
    } catch (IOException ex) {
      System.err.println("Error: could not write to relpair.ped");
      System.exit(1);
    }

    for (int i = 0; i < chrs.length; i++) {
      try {
        readers[i].close();
      } catch (IOException ioe) {
        ioe.printStackTrace();
      }
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + "relpair.ctl"));
      writer.println("relpair.loc");
      writer.println("relpair.ped");
      writer.println("relpair.out");
      writer.println("all"); // can be changed to "family" if you don't
      // want to check individuals between
      // families
      writer.println("n");
      writer.println("n");
      writer.println("F");
      writer.println("M");
      writer.println("50");
      writer.println("0.03");
      writer.println("1");
      writer.println("1000.0"); // for v2.0
      writer.close();
    } catch (IOException ex) {
      System.err.println("Error: could not write to relpair.ctl");
      System.exit(1);
    }
  }

  public static boolean filterPairs(String first, String second) {
    for (String[] element : FILTER_PAIRS) {
      if (element[0].equalsIgnoreCase(first) && element[1].equalsIgnoreCase(second)) {
        return false;
      }
    }
    return true;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "";

    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\";
    // boolean create = true;
    // boolean parse = false;

    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\Family structure
    // issues\\relpair\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\";
    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\Family structure issues\\fifthPass\\";
    boolean create = false;
    boolean parse = true;

    // String genome = "../../plink.genome";
    String genome = "plink.genome";
    boolean filter = true;

    String usage =
        "\n" + "link.Relpair requires 1-4 arguments (make sure you use -Xms1024M -Xmx1024M for parsing)\n"
            + "   (1) create files using all available re_chrom##.pre files (i.e. -createFiles (not the default))\n"
            + "   (2) parse results files (i.e. -parseResults (not the default))\n"
            + "   (3) include data from a plink genome file when parsing (i.e. genome=plink.genome (not the default))\n"
            + "   (4) filter out distantly related pairs (i.e. -filter (not the default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("-createFiles")) {
        create = true;
        numArgs--;
      } else if (arg.startsWith("-parseResults")) {
        parse = true;
        numArgs--;
      } else if (arg.startsWith("genome=")) {
        genome = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-filter")) {
        filter = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (!new File(dir).exists()) {
      System.err.println("Error - directory '" + dir + "' not found; trying current directory");
      dir = "";
    }
    if (genome.equals("") || !new File(dir + genome).exists()) {
      genome = null;
    }
    try {
      if (create) {
        createFiles(dir);
      }
      if (parse) {
        parseFiles(dir, genome, filter);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseFiles(String dir, String genome, boolean filter) {
    BufferedReader reader;
    PrintWriter withinWriter = null, betweenWriter = null;
    String[] line;
    String temp;
    Hashtable<String, String> within, between;
    String[] keys;

    within = new Hashtable<String, String>();
    between = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader(dir + "relpair.out"));
      do {
        temp = reader.readLine();
      } while (!temp.contains("STRONGEST DISCREPANCIES WITHIN FAMILIES"));
      for (int i = 0; i < 5; i++) {
        reader.readLine();
      }
      if (genome == null) {
        withinWriter = Files.getWriter(dir + "within.xln");
      }
      do {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length > 1 && (!filter || filterPairs(line[3], line[4]))) {
          if (line[5].equals(">")) {
            line[5] = "10^6";
            line[6] = "";
          }
          if (genome == null) {
            withinWriter.println(Array.toStr(line));
          } else {
            within.put(line[0] + "\t" + line[1] + "\t" + line[2],
                line[3] + "\t" + line[4] + "\t" + line[5]);
          }
        }
      } while (line.length > 1);

      if (genome == null) {
        withinWriter.close();
        betweenWriter = Files.getWriter(dir + "between.xln");
      }

      do {
        temp = reader.readLine();
      } while (!temp.contains("STRONGEST DISCREPANCIES ACROSS FAMILIES"));
      for (int i = 0; i < 5; i++) {
        reader.readLine();
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (!line[0].equals("NONE.") && (!filter || filterPairs(line[4], line[5]))) {
          if (line[6].equals(">")) {
            line[6] = "10^6";
          }
          if (genome == null) {
            betweenWriter.println(Array.toStr(line));
          } else {
            between.put(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3],
                line[4] + "\t" + line[5] + "\t" + line[6]);
          }
        }
      }
      if (genome == null) {
        betweenWriter.close();
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err
          .println("Error: file \"" + dir + "relpair.out" + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + "relpair.out" + "\"");
      System.exit(2);
    }

    if (genome != null) {
      withinWriter = Files.getWriter(dir + "within.xln");
      withinWriter.println("PED\tID1\tID2\tEZ\tZ0\tZ1\tZ2\tPI_HAT\tPUT\tINF\tRATIO");
      betweenWriter = Files.getWriter(dir + "between.xln");
      betweenWriter.println("PED1\tID1\tPED2\tID2\tZ0\tZ1\tZ2\tPI_HAT\tPUT\tINF\tRATIO");
      try {
        reader = new BufferedReader(new FileReader(dir + genome));
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");

          if (within.containsKey(line[0] + "\t" + line[1] + "\t" + line[3])) {
            withinWriter.println(line[0] + "\t" + line[1] + "\t" + line[3] + "\t" + line[5] + "\t"
                + line[6] + "\t" + line[7] + "\t" + line[8] + "\t" + line[9] + "\t"
                + within.remove(line[0] + "\t" + line[1] + "\t" + line[3]));
            withinWriter.flush();
          }

          if (between.containsKey(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3])) {
            betweenWriter.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t"
                + line[6] + "\t" + line[7] + "\t" + line[8] + "\t" + line[9] + "\t"
                + between.remove(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3]));
            betweenWriter.flush();
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error - file \"" + dir + genome + "\" not found in current directory");
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + genome + "\"");
        ioe.printStackTrace();
      }
      keys = HashVec.getKeys(within);
      for (String key : keys) {
        withinWriter
            .println(key + "\t" + Array.toStr(Array.stringArray(5, ".")) + "\t" + within.get(key));
      }
      keys = HashVec.getKeys(between);
      for (String key : keys) {
        betweenWriter
            .println(key + "\t" + Array.toStr(Array.stringArray(4, ".")) + "\t" + within.get(key));
      }

      withinWriter.close();
      betweenWriter.close();
    }
  }
}
