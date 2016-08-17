package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.stats.McNemarsTest;

public class MitoDNA {
  public static final String DEFAULT_DIR =
      "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\mtDNA\\";
  // public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My
  // Documents\\mtDNA\\";
  public static final String DEFAULT_REFS = "referenceHaps.dat";
  public static final String DEFAULT_SNPS = "db.xln";
  // public static final String DEFAULT_SNPS = "errs.xln";
  public static final String DEFAULT_PAIRS = "pairs.dat";

  public static final String[] ANALYSIS_VARS =
      {"IJK", "JTUK", "JTIWX", "mito1719_dom", "mito4580_dom", "mito7028_dom", "mito8251_dom",
          "mito9055_dom", "mito10398_dom", "mito10400_dom", "mito12308_dom", "mito13368_dom",
          "mito13708_dom", "mito16391_dom", "H", "I", "J", "K", "M", "T", "U", "V", "W", "X",
          "Other", "I_v_H", "J_v_H", "K_v_H", "M_v_H", "T_v_H", "U_v_H", "V_v_H", "W_v_H", "X_v_H"};

  public static void hapit(String dir, String referenceHaps, String snps) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp, trav, genos;
    Vector<String> names = new Vector<String>();
    Vector<String[]> v = new Vector<String[]>();
    int count, missing;
    String[] refSNPs = null;
    String[][] refs;
    int[] indices;
    boolean prob;
    boolean[] poss;
    String[] trans;

    try {
      reader = new BufferedReader(new FileReader(dir + referenceHaps));
      line = reader.readLine().trim().split("[\\s]+");
      if (!line[0].toLowerCase().contains("haplo")) {
        System.err.println("Error - I do not think that word means what you think it means");
      }
      refSNPs = Array.subArray(line, 1);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length != refSNPs.length + 1) {
          System.err.println("Error - mismatched number of SNPs in reference haplotypes");
          reader.close();
          return;
        }
        names.add(line[0]);
        v.add(Array.subArray(line, 1));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err
          .println("Error: file \"" + dir + referenceHaps + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + referenceHaps + "\"");
      System.exit(2);
    }

    refs = Matrix.toStringArrays(v);
    System.out.println(
        "Found " + refs.length + " reference haplotypes using " + refSNPs.length + " SNPs each");

    try {
      reader = new BufferedReader(new FileReader(dir + snps));
      temp = reader.readLine();
      line = temp.trim().split("[\\s]+");
      indices = Array.intArray(refSNPs.length, -1);
      prob = false;
      for (int i = 0; i < line.length; i++) {
        for (int j = 0; j < refSNPs.length; j++) {
          if (line[i].contains(refSNPs[j])) {
            if (indices[j] == -1) {
              indices[j] = i;
            } else {
              System.err.println(
                  "Error - more than one column contains '" + refSNPs[j] + "'; please advise");
              prob = true;
            }
          }
        }
      }
      if (prob) {
        System.exit(1);
      }
      writer = new PrintWriter(new FileWriter(dir + ext.rootOf(snps) + "_haps.xln"));
      writer.print("Haplotype\t#missing\tmissingAny");
      for (String refSNP : refSNPs) {
        writer.print("\t" + refSNP + "_dom");
      }
      writer.print("\t" + Array.toStr(Array.toStringArray(names)));
      for (int i = 0; i < names.size(); i++) {
        writer.print("\t" + names.elementAt(i) + "_v_" + names.elementAt(0));
      }
      writer.println("\t" + temp);
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.trim().split("\t", -1);
        count = 0;
        missing = 0;
        genos = "";
        prob = false;
        for (int i = 0; i < refSNPs.length; i++) {
          trav = indices[i] == -1 ? "." : line[indices[i]];
          if (!trav.equals(".") && !trav.equals("1")) {
            if (trav.length() != 2) {
              writer.println("! INVALID GENOTYPE (" + trav + ")\t" + missing + "\t"
                  + (missing > 0 ? 1 : 0) + genos + "\t"
                  + Array.toStr(Array.stringArray(refs.length * 2, ".")) + "\t" + temp);
              prob = true;
              genos += "\t?";
            } else if (trav.charAt(0) != trav.charAt(1)) {
              count++;
              genos += "\t.";
            } else {
              genos += "\t" + (trav.charAt(0) < 69 ? "0" : "1");
            }
          } else {
            genos += "\t.";
            if (indices[i] >= 0) {
              missing++;
            }
          }
        }
        if (!prob) {
          if (count > 0) {
            writer.println("! " + count + " HET" + (count > 1 ? "S" : "") + "\t" + missing + "\t"
                + (missing > 0 ? 1 : 0) + genos + "\t"
                + Array.toStr(Array.stringArray(refs.length * 2, ".")) + "\t" + temp);
          } else {
            poss = Array.booleanArray(refs.length, true);
            for (int i = 0; i < refSNPs.length; i++) {
              trav = indices[i] == -1 ? "." : line[indices[i]];
              for (int j = 0; j < refs.length; j++) {
                if (!trav.equals(".") && !trav.equals("1") && !refs[j][i].equals(".")
                    && !trav.substring(1).equals(refs[j][i])) {
                  poss[j] = false;
                }
              }
            }
            count = Array.booleanArraySum(poss);
            if (count == 0) {
              writer.println("Other\t" + missing + "\t" + (missing > 0 ? 1 : 0) + genos + "\t"
                  + Array.toStr(Array.stringArray(refs.length, "0")) + "\t"
                  + Array.toStr(Array.stringArray(refs.length, ".")) + "\t" + temp);
            } else {
              trav = "";
              for (int i = 0; i < poss.length; i++) {
                if (poss[i]) {
                  trav += (trav.equals("") ? "" : " / ") + names.elementAt(i);
                }
              }
              trans = count == 1 ? Array.booleanArrayToStringArray(poss)
                  : Array.stringArray(refs.length, ".");
              writer.print((count > 1 ? count + " " : "") + trav + "\t" + missing + "\t"
                  + (missing > 0 ? 1 : 0) + genos + "\t" + Array.toStr(trans));
              for (int i = 0; i < poss.length; i++) {
                if (trans[i].equals("0") && !poss[0]) {
                  trans[i] = ".";
                }
              }
              writer.println("\t" + Array.toStr(trans) + "\t" + temp);
            }
          }
        }
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + snps + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + snps + "\"");
      System.exit(2);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = DEFAULT_DIR;
    String refs = DEFAULT_REFS;
    String snps = DEFAULT_SNPS;
    String pairs = DEFAULT_PAIRS;
    boolean hapit = false;
    boolean testit = true;

    String usage = "\\n" + "park.MitoDNA requires 0-1 arguments\n" + "   (1) directory (i.e. dir="
        + dir + " (default))\n" + "   (2) reference haplotypes (i.e. refs=" + refs + " (default))\n"
        + "   (3) SNP database (i.e. dir=" + snps + " (default))\n"
        + "   (4) create haplotype database (i.e. -hapit (not the default))\n"
        + "   (5) file with pairs of UniqueIDs (i.e. pairs=" + pairs + " (default))\n"
        + "   (6) test haplotype database (i.e. -testit (not the default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("refs=")) {
        refs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("snps=")) {
        snps = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-hapit")) {
        hapit = true;
        numArgs--;
      } else if (arg.startsWith("-testit")) {
        testit = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (hapit) {
        hapit(dir, refs, snps);
      }
      if (testit) {
        McNemarsTest.batch(dir, pairs, ext.rootOf(snps) + "_haps.xln", ANALYSIS_VARS);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
