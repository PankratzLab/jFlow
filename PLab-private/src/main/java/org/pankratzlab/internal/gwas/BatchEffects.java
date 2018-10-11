package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class BatchEffects {

  // public static final String DIR = "C:\\Documents and
  // Settings\\npankrat\\My Documents\\gwas\\batchEffects\\";
  public static final String DIR = "";

  public static final double[] THRESHOLDS = {1E-2, 1E-3, 1E-4, 1E-5, 1E-10, 1E-15, 1E-20};

  public static void setupBatchEffects(String pedfile, String dnas) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    String temp;
    Hashtable<String, String> hash = new Hashtable<>();
    Vector<String> v = new Vector<>();
    int numBatches;
    int[] set;

    try {
      reader = new BufferedReader(new FileReader(dnas));
      reader.readLine();
      while (reader.ready()) {
        temp = reader.readLine();
        if (temp.indexOf("\t") == -1) {
          line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        } else {
          line = temp.split("\t", -1);
        }
        hash.put(line[0], line[1]);
        HashVec.addIfAbsent(line[1], v);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dnas + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dnas + "\"");
      System.exit(2);
    }
    numBatches = v.size();

    try {
      reader = new BufferedReader(new FileReader(pedfile));
      writer = Files.openAppropriateWriter("batchListing.txt");
      writer.println("FID\tIID\t" + ArrayUtils.toStr(ArrayUtils.toStringArray(v)));
      while (reader.ready()) {
        line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
        set = ArrayUtils.intArray(numBatches, 1);
        set[v.indexOf(hash.get(line[6]))] = 2;
        writer.println(line[0] + "\t" + line[1] + "\t" + ArrayUtils.toStr(set));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + pedfile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + pedfile + "\"");
      System.exit(2);
    }

    try {
      writer = Files.openAppropriateWriter("batchEffects");
      for (int i = 1; i <= numBatches; i++) {
        // writer.println("plink --bfile pd_gwas --pheno
        // batchListing.txt --mpheno "+i+" --freq --maf 0.01 --geno 0.20
        // --mind 0.5 --out freq."+i);
        // writer.println("plink --bfile pd_gwas --pheno
        // batchListing.txt --mpheno "+i+" --assoc --out assoc."+i);
        // writer.println("plink --bfile pd_gwas --keep keeps.txt
        // --pheno batchListing.txt --mpheno "+i+" --test-missing --maf
        // 0 --geno 1 --mind 1 --out missing."+i);

        writer.println("plink --bfile pd_gwas --keep noWGAs.txt --pheno batchListing.txt --mpheno "
                       + i + " --assoc --out assoc." + i);
        // writer.println("plink --bfile pd_gwas --keep noWGAs.txt
        // --pheno batchListing.txt --mpheno "+i+" --test-missing --maf
        // 0 --geno 1 --mind 1 --out missing."+i);

      }
      writer.close();
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

  public static void parseBatchEffects(String dir) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, bins = null;
    int[][] missCounts, sigCounts;
    double d;

    System.out.println(ext.getTime());
    try {
      reader = new BufferedReader(new FileReader("batchListing.txt"));
      line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
      bins = ArrayUtils.subArray(line, 2);
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "batchListing.txt"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "batchListing.txt" + "\"");
      System.exit(2);
    }

    sigCounts = new int[bins.length][THRESHOLDS.length + 1];
    for (int i = 1; i <= bins.length; i++) {
      try {
        System.out.println("assoc." + i + ".assoc");
        reader = Files.getReader("assoc." + i + ".assoc", dir);
        reader.readLine();
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          if (!line[8].equals("NA") && Integer.parseInt(line[0]) < 23) {
            d = Double.parseDouble(line[8]);
            for (int j = 0; j < THRESHOLDS.length; j++) {
              if (d < THRESHOLDS[j]) {
                sigCounts[i - 1][j]++;
              }
            }
            sigCounts[i - 1][THRESHOLDS.length]++;
          }
        }
        reader.close();
      } catch (Exception e) {}
    }

    missCounts = new int[bins.length][THRESHOLDS.length + 1];
    for (int i = 1; i <= bins.length; i++) {
      try {
        System.out.println("missing." + i + ".missing");
        reader = Files.getReader("missing." + i + ".missing", dir);
        reader.readLine();
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          if (!line[4].equals("NA") && Integer.parseInt(line[0]) < 23) {
            d = Double.parseDouble(line[4]);
            for (int j = 0; j < THRESHOLDS.length; j++) {
              if (d < THRESHOLDS[j]) {
                missCounts[i - 1][j]++;
              }
            }
            missCounts[i - 1][THRESHOLDS.length]++;
          }
        }
        reader.close();
      } catch (Exception e) {}
    }

    try {
      writer = Files.openAppropriateWriter("batchEffectsCounts.xln");
      for (int i = 0; i < bins.length; i++) {
        writer.println(bins[i] + "\tASSOC\t" + ArrayUtils.toStr(sigCounts[i]));
        writer.println(bins[i] + "\tMISS\t" + ArrayUtils.toStr(missCounts[i]));
      }
      writer.close();
    } catch (Exception e) {}
    System.out.println(ext.getTime());
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String pedfile = "pedfile.txt";
    String dnas = "dnas.txt";

    String usage = "\\n" + "park.gwa.BatchEffects requires 0-1 arguments\n"
                   + "   (1) pedigree file (i.e. ped=" + pedfile + " (default))\n"
                   + "   (2) list of DNAs and their batch (i.e. dnas=" + dnas + " (default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("ped=")) {
        pedfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("dnas=")) {
        dnas = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      // setupBatchEffects(pedfile, dnas);
      parseBatchEffects(DIR);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}