package org.genvisis.gwas;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.StringVector;
import org.genvisis.common.Vectors;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialIntMatrix;
import org.genvisis.filesys.SerialStringMatrix;
import org.genvisis.parse.GenParser;

public class Fst {
  private final double[] p; // allele frequency of the second allele
  private final double[] f; // inbreeding coefficient (higher indicates inbreeding, negative
                            // indicates isolate breaking (i.e. Wahlund effect))
  private final double[] hetObs;
  private final double[] hetExp;
  private double hI; // observed heterozygosities in individuals in subpopulations
  private double hS; // expected heterozygosities in subpopulations
  private final double hT; // expected heterozygosities for overall total population
  private final double fis;
  private final double fst;
  private final double fit;

  public Fst(int[][] counts, boolean unweightedBySampleSize) {
    double pBar;
    int[] sums;

    p = new double[counts.length];
    hetObs = new double[counts.length];
    hetExp = new double[counts.length];
    f = new double[counts.length];
    sums = new int[counts.length];
    for (int i = 0; i < counts.length; i++) {
      sums[i] = Array.sum(counts[i]);
      p[i] = (double) (counts[i][2] * 2 + counts[i][1]) / (double) (2 * sums[i]);
      hetObs[i] = (double) counts[i][1] / (double) sums[i];
      hetExp[i] = 1 - (p[i] * p[i] + (1 - p[i]) * (1 - p[i]));
      f[i] = (hetExp[i] - hetObs[i]) / hetExp[i];
    }

    if (unweightedBySampleSize) {
      pBar = Array.mean(p);
      hI = Array.mean(hetObs);
      hS = Array.mean(hetExp);
    } else {
      pBar = 0;
      hI = 0;
      hS = 0;
      for (int i = 0; i < counts.length; i++) {
        pBar += p[i] * sums[i];
        hI += hetObs[i] * sums[i];
        hS += hetExp[i] * sums[i];
      }
      pBar /= Array.sum(sums);
      hI /= Array.sum(sums);
      hS /= Array.sum(sums);
    }

    hT = 1 - (pBar * pBar + (1 - pBar) * (1 - pBar));

    fis = (hS - hI) / hS;
    fst = (hT - hS) / hT;
    fit = (hT - hI) / hT;
  }

  public double[] getP() {
    return p;
  }

  public double[] getF() {
    return f;
  }

  public double getHI() {
    return hI;
  }

  public double getHS() {
    return hS;
  }

  public double getHT() {
    return hT;
  }

  public double getFis() {
    return fis;
  }

  public double getFst() {
    return fst;
  }

  public double getFit() {
    return fit;
  }

  public static double[] calcPs(int[][] counts) {
    double[] p;

    p = new double[counts.length];
    for (int i = 0; i < counts.length; i++) {
      p[i] = (double) (counts[i][2] * 2 + counts[i][1]) / (double) (2 * Array.sum(counts[i]));
    }

    return p;
  }

  // if all you want is Fst, then this is more than twice as fast for large numbers of markers
  public static double calcFst(int[][] counts, boolean unweightedBySampleSize) {
    double[] p;
    double[] hetExp;
    double hS;
    double hT;
    double pBar;
    int[] sums;

    p = new double[counts.length];
    hetExp = new double[counts.length];
    sums = new int[counts.length];
    for (int i = 0; i < counts.length; i++) {
      sums[i] = Array.sum(counts[i]);
      p[i] = (double) (counts[i][2] * 2 + counts[i][1]) / (double) (2 * sums[i]);
      hetExp[i] = 1 - (p[i] * p[i] + (1 - p[i]) * (1 - p[i]));
    }

    if (unweightedBySampleSize) {
      pBar = Array.mean(p);
      hS = Array.mean(hetExp);
    } else {
      pBar = 0;
      hS = 0;
      for (int i = 0; i < counts.length; i++) {
        pBar += p[i] * sums[i];
        hS += hetExp[i] * sums[i];
      }
      pBar /= Array.sum(sums);
      hS /= Array.sum(sums);
    }

    hT = 1 - (pBar * pBar + (1 - pBar) * (1 - pBar));

    return (hT - hS) / hT;
  }

  public static String[][] prepFiles(String filename, String plinkRoot) {
    String[] line, keys, array, files;
    String trav;
    Hashtable<String, Vector<String>> hash;
    Vector<String> v;
    int count;
    String dir, root, list;
    GenParser parser;
    String[][] markerInfo;
    int[][] alleleCounts;
    boolean reverse;

    dir = ext.parseDirectoryOfFile(plinkRoot);
    root = ext.removeDirectoryInfo(plinkRoot);

    if (new File(dir + root + ".bim.ser").exists()) {
      markerInfo = SerialStringMatrix.load(dir + root + ".bim.ser", false).getMatrix();
    } else if (!new File(dir + root + ".bim").exists()) {
      System.err.println("Error - algorithm requires a '" + root
                         + ".bim' file to ensure marker order and allelic strand");
      return null;
    } else {
      markerInfo = HashVec.loadFileToStringMatrix(dir + root + ".bim", false,
                                                  new int[] {1, 0, 3, 4, 5}, false);
      new SerialStringMatrix(markerInfo).serialize(dir + root + ".bim.ser");
      files = Files.list(dir, ".hwe.ser", false);
      for (String file : files) {
        count = 1;
        do {
          trav = file;
          for (int j = 0; j < count; j++) {
            trav += "_";
          }
        } while (new File(dir + trav).exists());
        new File(dir + file).renameTo(new File(dir + trav));
      }
    }

    hash = HashVec.loadFileToHashVec(filename, 2, new int[] {0, 1}, "\t", false, false);
    keys = HashVec.getKeys(hash);
    for (int i = 0; i < keys.length; i++) {
      v = hash.get(keys[i]);
      if (v.size() == 1) {
        if (!v.elementAt(0).equals("FID\tIID")) {
          System.err.println("Error - found population subgroup (" + keys[i]
                             + ") with only one member: " + v.elementAt(0)
                             + "; this will not be run");
          return null;
        }
      } else {
        if (v.size() < 20) {
          System.err.println("Warning - subgroup " + keys[i] + " only has " + v.size()
                             + " members");
        }
        list = keys[i] + ".list";
        v.insertElementAt("FID\tIID", 0);
        array = Array.toStringArray(v);
        if (!new File(dir + list).exists()
            || !Array.equals(array,
                             HashVec.loadFileToStringArray(dir + list, false, new int[] {0, 1},
                                                           false),
                             false)
            || !new File(dir + keys[i] + ".hwe.ser").exists()) {
          Files.writeList(array, dir + list);
          CmdLine.run("plink --bfile " + root + " --keep " + keys[i] + ".list --hardy --out "
                      + keys[i], dir);
          line = new String[] {dir + keys[i] + ".hwe", "!2=ALL", "1", "3", "4", "5"};
          parser = new GenParser(line, null);
          ext.checkHeader(parser.getColumnNames(), new String[] {"SNP", "A1", "A2", "GENO"}, true);
          count = 0;
          alleleCounts = new int[markerInfo.length][];
          while (parser.ready()) {
            line = parser.nextLine();
            if (line != null) {
              if (!line[0].equals(markerInfo[count][0])) {
                System.err.println("Error - mismatched marker order for " + keys[i] + ".hwe"
                                   + "; expecting " + markerInfo[count][0]
                                   + " from .bim file, found " + line[0]);
                return null;
              } else if (line[1].equals(markerInfo[count][4])
                         && line[2].equals(markerInfo[count][3])) {
                reverse = true;
              } else if (!line[1].equals(markerInfo[count][3])
                         || !line[2].equals(markerInfo[count][4])) {
                System.err.println("Error - mismatched alleles in " + keys[i] + ".hwe"
                                   + " for marker " + line[0] + "; expecting "
                                   + markerInfo[count][3] + "/" + markerInfo[count][4] + ", found "
                                   + line[1] + "/" + line[2]);
                return null;
              } else {
                reverse = false;
              }
              alleleCounts[count] =
                                  Array.toIntArray(Sort.putInOrder(line[3].split("/"),
                                                                   reverse ? new int[] {2, 1, 0}
                                                                           : new int[] {0, 1, 2}));
              count++;
            }
          }
          new SerialIntMatrix(alleleCounts).serialize(dir + keys[i] + ".hwe.ser");
        } else {
          alleleCounts = SerialIntMatrix.load(dir + keys[i] + ".hwe.ser", false).getMatrix();
        }
      }
    }

    return markerInfo;
  }

  public static void runAnalyses(String filename, String plinkRoot, String[] incl) {
    PrintWriter writer;
    String[] files;
    int count;
    String[][] markerInfo;
    int[][] counts;
    int[][][] alleleCounts;
    String dir;

    prepFiles(filename, plinkRoot);

    dir = ext.parseDirectoryOfFile(filename);
    if (incl == null) {
      System.out.println("No specific ethnicities provided, using all in Fst calculations:");
      files = Files.list(dir, ".hwe.ser", false);
      incl = new String[files.length];
      for (int i = 0; i < files.length; i++) {
        incl[i] = files[i].substring(0, files[i].indexOf(".hwe.ser"));
      }
      System.out.println(ext.listWithCommas(incl, false));
    } else {
      count = 0;
      for (int i = 0; i < incl.length; i++) {
        if (!new File(dir + incl[i] + ".hwe.ser").exists()) {
          System.err.println("Error - " + incl[i] + ".hwe.ser not found in directory");
          count++;
        }
      }
      if (count > 0) {
        System.err.println("Error - could not find " + count + " of the " + incl.length
                           + " necessary files to perform calculations");
        return;
      }
    }

    if (incl.length == 0) {
      System.err.println("Error - cannot compute Fst with no valid populations provided...");
      return;
    } else if (incl.length == 1) {
      System.err.println("Error - cannot compute Fst using a single (" + incl[0] + ") population");
      return;
    }

    markerInfo = prepFiles(filename, plinkRoot);
    if (markerInfo == null) {
      System.err.println("Error - file prep failed");
      return;
    }

    alleleCounts = new int[incl.length][][];
    for (int i = 0; i < incl.length; i++) {
      alleleCounts[i] = SerialIntMatrix.load(dir + incl[i] + ".hwe.ser", false).getMatrix();
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + Array.toStr(incl, "-") + "_Fst.xln"));
      writer.println("Marker\tFst_equal\tFst_weighted\t" + Array.toStr(incl));
      counts = new int[incl.length][];
      for (int i = 0; i < markerInfo.length; i++) {
        for (int j = 0; j < alleleCounts.length; j++) {
          counts[j] = alleleCounts[j][i];
        }
        writer.println(markerInfo[i][0] + "\t" + calcFst(counts, true) + "\t"
                       + calcFst(counts, false) + "\t" + Array.toStr(calcPs(counts)));

      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename, false) + "_Fst.xln");
      e.printStackTrace();
    }
  }

  public static void splitIndividualsIntoNgroups(String filename, int n) {
    String[] keys, inds;
    Hashtable<String, Vector<String>> hash;
    int count;
    StringVector[] fileContents;

    ext.checkHeader(Files.getHeaderOfFile(filename, "[\\s]+", new Logger()),
                    new String[] {"FID", "IID"}, new int[] {0, 1}, false, new Logger(), true);
    hash = HashVec.loadFileToHashVec(filename, 2, new int[] {0, 1}, "\t", true, false);
    keys = HashVec.getKeys(hash);

    count = 0;
    fileContents = Vectors.initializedArray(StringVector.class, n);
    for (String key : keys) {
      inds = Array.toStringArray(hash.get(key));
      for (String ind : inds) {
        fileContents[count % n].add(ind + "\t" + key);
        count++;
      }
    }

    for (int i = 0; i < fileContents.length; i++) {
      String[] s = new String[fileContents[i].size()];
      Files.writeList(fileContents[i].toArray(s),
                      ext.rootOf(filename, false) + "." + (i + 1) + ".dat");
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String plinkRoot = "plink";
    String filename = "";
    String[] incl = null;
    int n = -1;

    String usage = "\n" + "gwas.Fst requires 0-1 arguments\n"
                   + "   (1) filename with ancestry (i.e. file=" + filename + " (default))\n"
                   + "   (2) ethnicities to include in calculations (i.e. incl=Ashk,Brit,Ital (default all defined populations))\n"
                   + "   (3) PLINK root (i.e. root=" + plinkRoot + " (default))\n" + "  OR\n"
                   + "   (2) number of files to evenly split ethnicites into (i.e. split=2 (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("root=")) {
        plinkRoot = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("incl=")) {
        incl = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("split=")) {
        n = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // filename = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Consortium\\Fst\\testAshk\\ashkbrit.txt";
    // plinkRoot = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Consortium\\Fst\\testAshk\\plink";

    // filename = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Consortium\\Fst\\EuropeanClasses.txt";
    // plinkRoot = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Consortium\\Fst\\plink";
    // n = 2;

    filename =
             "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\Fst\\Discovery\\EuropeanClasses.1clean.dat";
    plinkRoot =
              "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\Fst\\Discovery\\plink";
    // incl = new String[] {"Ashk", "British"};
    incl = new String[] {"Ashk", "Italian"};

    try {
      if (n > 0) {
        splitIndividualsIntoNgroups(filename, n);
      } else if (!filename.equals("")) {
        runAnalyses(filename, plinkRoot, incl);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
