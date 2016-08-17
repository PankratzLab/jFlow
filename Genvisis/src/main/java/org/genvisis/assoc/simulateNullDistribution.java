// You just need a database and the indices of the relevant columns, pedigree.pre is made for you
// make sure to trim it to only contain the individuals you want to permute (i.e. VPD, no caustitive
// mutation, no missing data)
// Note: this algorithm does not currently require specifying gender, so any individuals not found
// in ninfo2 will default to being female
// currently have to run allegro simulate.opt by hand
// currently PD specific (required ninfo2)
package org.genvisis.assoc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.link.LinkageMap;
import org.genvisis.link.TrimFam;
import org.genvisis.park.tools;
import org.genvisis.stats.LogisticRegression;
import org.genvisis.stats.RegressionModel;

public class simulateNullDistribution {
  public static final int NUM_REPS = 10000;
  public static final int REP_STEP = 100;
  public static final boolean FREQS_FROM_CONTROLS_ONLY = false;
  public static final int FAM_REPS_DEFAULT = 10000;
  public static final int BOOT_REPS_DEFAULT = 1000;

  public static Hashtable<String, Vector<String>> cloneHash(
      Hashtable<String, Vector<String>> hash) {
    Hashtable<String, Vector<String>> clone = new Hashtable<String, Vector<String>>();
    String[] keys = HashVec.getKeys(hash);

    for (String key : keys) {
      clone.put(key, hash.get(key));
    }

    return clone;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "Rep1.db.xls";
    String filename = "SNCA3_database_VPD_noCM_Cauc.xln";
    String outputfile = "pedigree.pre";
    int numReps = NUM_REPS;
    boolean controlFreqs = FREQS_FROM_CONTROLS_ONLY;
    int famIDcol = 2;
    int indIDcol = 3;
    // int affCol = 26;
    int affCol = 28;
    // int alleleCol = 14;
    int alleleCol = 6;
    // String targets = "263,259";
    String targets = "2";

    String usage = "\n" + "park.simulateNullDistribution requires 0-7 arguments\n"
        + "   (1) input filename (i.e. file=" + filename + " (default)\n"
        + "   (2) output filename (i.e. out=" + filename + " (default)\n"
        + "   (3) number of replicates to create (i.e. reps=" + numReps + " (default)\n"
        + "   (4) column of FamIDs (i.e. fam=" + famIDcol + " (default)\n"
        + "   (5) column of IndIDs (i.e. ind=" + indIDcol + " (default)\n"
        + "   (6) column of affection status (i.e. aff=" + affCol + " (default)\n"
        + "   (7) column of first allele (i.e. allele=" + alleleCol + " (default)\n"
        + "   (8) target alleles (separated by commas) (i.e. target=" + alleleCol + " (default)\n"
        + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outputfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("reps=")) {
        numReps = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("fam=")) {
        famIDcol = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("ind=")) {
        indIDcol = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("aff=")) {
        affCol = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("allele=")) {
        alleleCol = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("targets=")) {
        targets = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new simulateNullDistribution(filename, outputfile, numReps, controlFreqs, famIDcol - 1,
          indIDcol - 1, affCol - 1, alleleCol - 1, targets);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public simulateNullDistribution(String filename, String outputfile, int reps,
      boolean controlFreqs, int famCol, int indCol, int affCol, int allCol, String targets) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line = null, trav, fams;
    String temp;
    Hashtable<String, Hashtable<String, String>> phenos =
        new Hashtable<String, Hashtable<String, String>>();
    Hashtable<String, String> affs;
    Hashtable<String, Vector<String>> vips = new Hashtable<String, Vector<String>>(), backup;
    Hashtable<String, Vector<String>> genos = new Hashtable<String, Vector<String>>();
    Vector<String> v;
    IntVector alleles = new IntVector();
    DoubleVector counts = new DoubleVector();
    int numParticipants, allele;
    Vector<String> pre = new Vector<String>();
    String prev;
    TrimFam tf;
    int numCols;
    int[] keys;
    boolean done;
    int[] targetIndices;
    int[][][] indeps;
    int[] deps;
    int count;
    RegressionModel model;
    double[][] alleleFreqs;
    double[] freqs;

    try {
      reader = new BufferedReader(new FileReader(filename));
      numCols = reader.readLine().split("\t").length;
      if (numCols < Array.max(new int[] {famCol, indCol, affCol, allCol})) {
        System.err.println("Error - there are fewer columns than the max column specified");
        System.exit(1);
      }
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        if (line.length != numCols) {
          System.err.println("Inconsistent number of columns");
        }
        if (!line[allCol + 0].equals(".") && !line[allCol + 0].equals("0")
            && !line[allCol + 1].equals(".") && !line[allCol + 1].equals("0")) {
          HashVec.addToHashVec(vips, line[famCol], line[indCol], false);
          HashVec.addToHashVec(genos, line[famCol], line[allCol + 0], false);
          HashVec.addToHashVec(genos, line[famCol], line[allCol + 1], false);
          HashVec.addToHashHash(phenos, line[famCol], line[indCol], line[affCol]);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
    try {
      reader = tools.getNinfoReader(2);
      writer = new PrintWriter(new FileWriter(outputfile));
      prev = "";
      done = false;
      // backup = vips.clone();
      backup = cloneHash(vips);
      while (!done) {
        if (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
        } else {
          prev = "";
          done = true;
        }
        if (done || !line[0].equals(prev)) {
          if (vips.containsKey(prev)) {
            System.err.println("ERROR! YOU HAVE YET TO CONVERT TO THE NEW INPUT FOR TrimFam!");
            // tf = new TrimFam(pre, vips.get(prev), false, true);
            tf = new TrimFam(pre, false, false, true, TrimFam.SCORE_99_NAMING_SCHEME, 0, false,
                false, new Logger());
            v = tf.getExtendedFamilyInformation();
            affs = phenos.get(prev);
            for (int i = 0; i < v.size(); i++) {
              trav = v.elementAt(i).split("[\\s]+");
              writer
                  .println(
                      Array.toStr(trav, "\t") + "\t"
                          + (affs.containsKey(trav[1])
                              ? (Integer.parseInt(affs.get(trav[1])) + 1) + "\t1\t1" : "0\t0\t0")
                          + "");
            }
          }
          pre.removeAllElements();
          vips.remove(prev);
          if (controlFreqs) {
            genos.remove(prev);
          }
        }
        pre.add(line[0] + "\t" + line[1] + "\t" + line[4] + "\t" + line[5] + "\t"
            + (line[2].equals("F") ? "2" : (line[2].equals("M") ? "1" : "0")));
        prev = line[0];
      }
      fams = HashVec.getKeys(vips);
      for (String fam : fams) {
        v = vips.get(fam);
        affs = phenos.get(fam);
        if (v.size() > 1) {
          System.err
              .println("Error - more than one individual found in a family not found in ninfo2: "
                  + fam + "\n" + "        The following indivuals will be listed as unrelated: "
                  + Array.toStr(Array.toStringArray(v), ",") + "");
        }
        for (int j = 0; j < v.size(); j++) {
          writer
              .println(fam + "\t" + v.elementAt(j) + "\t0\t0\t2\t"
                  + (affs.containsKey(v.elementAt(j))
                      ? (Integer.parseInt(affs.get(v.elementAt(j))) + 1) + "\t1\t1" : "0\t9\t9")
                  + "");
        }

      }
      vips = backup;

      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + outputfile + "\" is otherwise occupied");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println(
          "Error reading file \"" + filename + "\" or writing file \"" + outputfile + "\"");
      System.exit(2);
    }

    fams = HashVec.getKeys(genos);
    numParticipants = 0;
    for (String fam : fams) {
      v = genos.get(fam);
      numParticipants += v.size() / 2;
      for (int j = 0; j < v.size(); j++) {
        allele = Integer.parseInt(v.elementAt(j));
        int index = alleles.indexOf(allele);
        if (index >= 0) {
          counts.set(index, counts.get(index) + (2.0 / v.size()));
        } else {
          alleles.add(allele);
          counts.add(2.0 / v.size());
        }
      }

    }
    System.out.println("Allele frequencies were permuted one per family (" + numParticipants
        + " from " + fams.length + " families)");
    keys = Sort.quicksort(alleles, Sort.ASCENDING);
    line = targets.split(",");
    targetIndices = new int[line.length];
    for (int i = 0; i < line.length; i++) {
      targetIndices[i] = keys[alleles.indexOf(Integer.parseInt(line[i]))];
    }

    freqs = new double[alleles.size()];
    for (int j = 0; j < alleles.size(); j++) {
      freqs[j] = counts.elementAt(keys[j]) / ((double) fams.length * 2);
    }
    alleleFreqs = new double[REP_STEP][];
    for (int i = 0; i < REP_STEP; i++) {
      alleleFreqs[i] = freqs;
    }
    new LinkageMap(1, Array.stringArraySequence(REP_STEP, "Rep"), alleleFreqs,
        Array.doubleArray(REP_STEP, 50), false, false).createFile("map.dat");

    try {
      writer = new PrintWriter(new FileWriter("simulate.opt"));
      writer.println("% Read input in LINKAGE style format:\n" + "PREFILE " + outputfile + "\n"
          + "DATFILE map.dat\n\n" + "% Simulate stroke reconstruction pedigrees\n"
          + "SIMULATE het:1\n\n" + "% Other options:\n" + "MAXMEMORY 100");
      writer.close();
    } catch (IOException ioe) {
      System.err.println("Problem writing map.dat");
    }

    for (int i = 0; i < reps / REP_STEP; i++) {
      try {
        Runtime.getRuntime().exec("allegro simulate.opt").waitFor();
      } catch (Exception e) {
        System.err.println("Error - running allegro at the command prompt");
      }

      try {
        reader = new BufferedReader(new FileReader(outputfile + "." + ext.formNum(i + 1, 3)));
        writer = new PrintWriter(new FileWriter("tr_" + outputfile + ".1"));
        while (reader.ready()) {
          temp = reader.readLine();
          for (int j = 0; j < temp.length(); j++) {
            writer.print(temp.charAt(j) == 32 ? "\t"
                : (temp.charAt(j) < 10 ? (char) (temp.charAt(j) + 48) : temp.charAt(j)));
          }
          writer.println();
        }
        reader.close();
        writer.close();
      } catch (FileNotFoundException fnfe) {
        System.err
            .println("Error: file \"" + outputfile + ".1" + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + outputfile + ".1" + "\"");
        System.exit(2);
      }

      count = 0;
      fams = new String[numParticipants];
      deps = new int[numParticipants];
      indeps = new int[REP_STEP][numParticipants][targetIndices.length];
      try {
        reader = new BufferedReader(new FileReader("tr_" + outputfile + ".1"));
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          v = vips.get(line[0]);
          affs = phenos.get(line[0]);
          if (v.contains(line[1])) {
            fams[count] = line[0];
            deps[count] = Integer.parseInt(affs.get(line[1]));
            for (int rep = 0; rep < REP_STEP; rep++) {
              for (int j = 0; j < targetIndices.length; j++) {
                for (int k = 0; k < 2; k++) {
                  if (Integer.parseInt(line[6 + rep * 2 + k]) == targetIndices[j] + 1) {
                    indeps[rep][count][j]++;
                  }
                }
              }
            }
            count++;
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println(
            "Error: file \"" + "tr_" + outputfile + ".1" + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + "tr_" + outputfile + ".1" + "\"");
        System.exit(2);
      }

      try {
        writer = new PrintWriter(new FileWriter((i + 1) + ".out"));
        for (int j = 0; j < REP_STEP; j++) {
          model = new LogisticRegression(deps, indeps[j], false, false);
          model.onePerFamily(fams, FAM_REPS_DEFAULT, BOOT_REPS_DEFAULT);
          writer.println(Array.toStr(model.getStats()) + "\t" + Array.toStr(model.getBetas()));
          writer.flush();
        }
        writer.close();
      } catch (IOException ioe) {
        System.err.println("Error writing " + i + ".out");
      }
    }
  }
}
