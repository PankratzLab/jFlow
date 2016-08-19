package org.genvisis.cnv.park;

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
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.stats.Maths;

public class CompareCalls {
  public static final String DEFAULT_ROOT =
                                          "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\allCalls\\";

  // public static final String[] DEFAULT_FILES = {"conf.cnv", "allMarkers.cnv"};
  // public static final String[] DEFAULT_FILES = {"conf_100kb_5SNP_10.0.cnv",
  // "allMarkers_100kb_5SNP_10.0.cnv"};
  public static final String[] DEFAULT_FILES = {"conf_100kb_5SNP_10.0.cnv",
                                                "conf_100kb_20SNP_10.0.cnv"};

  public static void compare(String rootDir, String[] files) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, inds;
    Hashtable<String, Hashtable<String, Vector<CNVariant>>> hash =
                                                                 new Hashtable<String, Hashtable<String, Vector<CNVariant>>>();
    Hashtable<String, Vector<CNVariant>> source = new Hashtable<String, Vector<CNVariant>>();
    CNVariant[][] cnvs;
    Vector<CNVariant> v = new Vector<CNVariant>();
    int match;
    int[] counts;
    int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(files.length, 2);

    for (int i = 0; i < files.length; i++) {
      try {
        reader = new BufferedReader(new FileReader(rootDir + files[i]));
        if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER,
                             false)) {
          reader.close();
          return;
        }
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          if (hash.containsKey(line[0] + "\t" + line[1])) {
            source = hash.get(line[0] + "\t" + line[1]);
          } else {
            hash.put(line[0] + "\t" + line[1], source = new Hashtable<String, Vector<CNVariant>>());
          }
          if (source.containsKey(i + "")) {
            v = source.get(i + "");
          } else {
            source.put(i + "", v = new Vector<CNVariant>());
          }
          v.add(new CNVariant(line, i));
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + rootDir + files[i]
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + rootDir + files[i] + "\"");
        System.exit(2);
      }
    }

    inds = HashVec.getKeys(hash);
    for (int[] allPossibleCombination : allPossibleCombinations) {
      try {
        writer = new PrintWriter(new FileWriter(rootDir + "Compare "
                                                + ext.rootOf(files[allPossibleCombination[0]])
                                                + " and "
                                                + ext.rootOf(files[allPossibleCombination[1]])
                                                + ".xln"));
        writer.println("FID\tIID\tTotal" + ext.rootOf(files[allPossibleCombination[0]]) + "\tTotal"
                       + ext.rootOf(files[allPossibleCombination[1]]) + "\tUnique"
                       + ext.rootOf(files[allPossibleCombination[0]]) + "\tUnique"
                       + ext.rootOf(files[allPossibleCombination[1]]) + "\tOverlaps\tExactMatches");
        for (String ind : inds) {
          cnvs = new CNVariant[][] {
                                    CNVariant.toCNVariantArray(hash.get(ind)
                                                                   .get(allPossibleCombination[0]
                                                                        + "")),
                                    CNVariant.toCNVariantArray(hash.get(ind).get(
                                                                                 allPossibleCombination[1] + ""))};
          counts = new int[4];
          if (cnvs[0].length == 0) {
            System.err.println("Error - " + ind + " not found in "
                               + files[allPossibleCombination[0]]);
          }
          if (cnvs[1].length == 0) {
            System.err.println("Error - " + ind + " not found in "
                               + files[allPossibleCombination[1]]);
          }
          for (int a = 0; a < cnvs[0].length; a++) {
            match = 0;
            for (int b = 0; b < cnvs[1].length; b++) {
              if (cnvs[0][a].equals(cnvs[1][b])) {
                match = 3;
                cnvs[1][b].setSource(99);
              } else if (match < 2 && cnvs[0][a].overlaps(cnvs[1][b])) {
                match = 2;
                cnvs[1][b].setSource(99);
              }
            }
            counts[match]++;
          }
          for (int b = 0; b < cnvs[1].length; b++) {
            match = 1;
            for (int a = 0; a < cnvs[0].length; a++) {
              if (cnvs[1][b].getSource() != 99 && cnvs[1][b].equals(cnvs[0][a])) {
                match = 3;
              } else if (match < 2 && cnvs[1][b].getSource() != 99
                         && cnvs[1][b].overlaps(cnvs[0][a])) {
                match = 2;
              }
            }
            if (cnvs[1][b].getSource() != 99) {
              counts[match]++;
            }
          }
          writer.println(ind + "\t" + cnvs[0].length + "\t" + cnvs[1].length + "\t"
                         + Array.toStr(counts));
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error comparing " + files[allPossibleCombination[0]] + " and "
                           + files[allPossibleCombination[1]]);
        e.printStackTrace();
      }
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String rootDirectory = DEFAULT_ROOT;
    String[] files = DEFAULT_FILES;

    String usage = "\\n" + "park.cnv.ComparePlinkResults requires 0-1 arguments\n"
                   + "   (1) directory (i.e. dir=" + rootDirectory + " (default))\n"
                   + "   (2) files to be compared (i.e. files=" + Array.toStr(files, ",")
                   + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        rootDirectory = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("files=")) {
        files = arg.split("=")[1].split(",");
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      compare(rootDirectory, files);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
