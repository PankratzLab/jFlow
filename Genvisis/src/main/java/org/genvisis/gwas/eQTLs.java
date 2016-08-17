package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.stats.Correlation;

public class eQTLs {
  public static final double DEFAULT_TRANS_THRESHOLD = 0.0001;
  public static final int DEFAULT_CIS_DISTANCE = 50000;
  public static final double DEFAULT_CIS_THRESHOLD = 0.001;

  public static final String[][] REQS_WITH_STDERR =
      {{"SNP", "Marker", "Name", "MarkerName", "name"}, {"A1", "Allele1", "REF"},
       {"Effect", "beta", "beta_SNP_add"}, {"P", "pval", "p-val", "p-value"},
       {"SE", "StdErr", "sebeta_SNP_add"}};
  public static final String[][] REQS_WITHOUT_STDERR =
      {{"SNP", "Marker", "Name", "MarkerName", "name"}, {"A1", "Allele1", "REF"},
       {"Effect", "beta", "beta_SNP_add"}, {"P", "pval", "p-val", "p-value"}};

  public static final String DEFAULT_TRANSCRIPT_BED =
      "D:/Myron/eQTLs/GEO_DataSet_GSE9703/GPL5188-122_withData.bed";

  public static void generateScores(String dir, String phenoFile, double transThreshold,
                                    int cisDistance, double cisThreshold, boolean divideByStderr) {
    // Hashtable<String,String> hash;
    String[] phenos;
    Logger log;
    byte chr;
    int start, stop;
    Segment probeSeg;
    double[][] pairs;

    log = new Logger(dir + ext.rootOf(phenoFile) + "_weights.log");
    phenos = Array.subArray(Files.getHeaderOfFile(dir + phenoFile, "\t", log), 2);

    // hash = HashVec.loadFileToHashString(DEFAULT_TRANSCRIPT_BED, new int[] {3}, new int[] {0, 1,
    // 2}, false, "\t", true, false, false);
    chr = 19;
    start = 10242517;
    stop = 10258291;
    probeSeg = new Segment(chr, start, stop);

    for (String pheno : phenos) {
      parseWeights(dir + pheno + "_SE1.out", dir + "allSNPs.dat", dir + pheno + "_meta.weights",
                   transThreshold, probeSeg, cisDistance, cisThreshold, false);
      // CmdLine.run("plink --bfile full --pheno plink_pheno.dat --keep whites.txt --score
      // "+phenos[i]+"_meta.weights --out "+phenos[i]+"_meta", dir);
      CmdLine.run("plink --bfile CEU --pheno pheno.dat --pheno-name " + pheno + " --score " + pheno
                  + "_meta.weights --out " + pheno + "_meta", dir);
      pairs = new double[][] {
                              Array.toDoubleArray(HashVec.loadFileToStringArray(dir + pheno
                                                                                + "_meta.profile",
                                                                                true, new int[] {2},
                                                                                false)),
                              Array.toDoubleArray(HashVec.loadFileToStringArray(dir + pheno
                                                                                + "_meta.profile", true, new int[] {5}, false))};
      pairs = Matrix.transpose(Matrix.removeRowsWithNaN(Matrix.transpose(pairs)));
      System.out.println(pheno + "_meta\t" + Array.toStr(Correlation.Pearson(pairs)));

      parseWeights(dir + "CEU_" + pheno + ".se.metal", dir + "allSNPs.dat",
                   dir + pheno + "_CEU.weights", transThreshold, probeSeg, cisDistance,
                   cisThreshold, false);
      // CmdLine.run("plink --bfile full --pheno plink_pheno.dat --keep whites.txt --score
      // "+phenos[i]+"_CEU.weights --out "+phenos[i]+"_CEU", dir);
      CmdLine.run("plink --bfile CEU --pheno pheno.dat --pheno-name " + pheno + " --score " + pheno
                  + "_CEU.weights --out " + pheno + "_CEU", dir);
      pairs = new double[][] {
                              Array.toDoubleArray(HashVec.loadFileToStringArray(dir + pheno
                                                                                + "_CEU.profile",
                                                                                true, new int[] {2},
                                                                                false)),
                              Array.toDoubleArray(HashVec.loadFileToStringArray(dir + pheno
                                                                                + "_CEU.profile", true, new int[] {5}, false))};
      pairs = Matrix.transpose(Matrix.removeRowsWithNaN(Matrix.transpose(pairs)));
      System.out.println(pheno + "_CEU\t" + Array.toStr(Correlation.Pearson(pairs)));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "pheno.dat";
    String dir = "";

    String usage = "\n" + "gwas.eQTLs requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("pheno=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/region/";
    // dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/chr1_crp/";
    // dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/chr1_selp/";
    // dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/IBCvariants/";
    // dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/filteredGenome/";
    dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/sigExtremes/";
    // dir = "D:/Myron/eQTLs/GEO_DataSet_GSE9703/HapMap2/icam1/";
    try {
      // run(dir, filename);
      generateScores(dir, filename, DEFAULT_TRANS_THRESHOLD, DEFAULT_CIS_DISTANCE,
                     DEFAULT_CIS_THRESHOLD, true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse(String filename) {
    for (int chr = 1; chr <= 22; chr++) {
      System.out.println(Files.getRunString() + " gwas.Minimac extract=snps/EUR.chr" + chr
                         + ".snps hapFile=hap/EUR/unzipped/EUR.chr" + chr
                         + ".hap mapFile=map/EUR.chr" + chr + ".map");
    }
  }

  public static void parseWeights(String filename, String mapFile, String outfile,
                                  double transThreshold, Segment probeSeg, int cisDistance,
                                  double cisThreshold, boolean divideByStderr) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, header;
    String[][] reqs;
    Hashtable<String, Segment> segs;
    Segment region;
    Logger log;
    int[] indices;
    double pval;

    log = new Logger(ext.rootOf(filename, false) + "_parseWeights.log");
    log.report("Parsing '" + filename + "'");
    if (mapFile != null) {
      segs = new SnpMarkerSet(mapFile, SnpMarkerSet.GENERIC_FORMAT_IGNORE_FIRST_LINE, false,
                              log).getSegments();
    } else {
      segs = null;
    }

    region = new Segment(probeSeg.getChr(), probeSeg.getStart() - cisDistance,
                         probeSeg.getStop() + cisDistance);

    try {
      header = Files.getHeaderOfFile(filename, "\t", log);
      if (divideByStderr) {
        reqs = REQS_WITH_STDERR;
      } else {
        reqs = REQS_WITHOUT_STDERR;
      }
      indices = ext.indexFactors(reqs, header, false, true, true, log, false);

      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(outfile));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (!ext.isMissingValue(line[indices[3]])) {
          pval = Double.parseDouble(line[indices[3]]);
          if (pval < transThreshold || (pval < cisThreshold && segs != null
                                        && region.overlaps(segs.get(line[indices[0]])))) {
            writer.println(line[indices[0]] + "\t" + line[indices[1]].toUpperCase() + "\t"
                           + (divideByStderr ? ext.formDeci(Double.parseDouble(line[indices[2]])
                                                            / Double.parseDouble(line[indices[4]]), 5)
                                             : line[indices[2]]));
          }
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }


  public static void run(String dir, String phenoFile) {
    String[] phenos, parameters;
    Logger log;

    log = new Logger(dir + ext.rootOf(phenoFile) + ".log");
    phenos = Array.subArray(Files.getHeaderOfFile(dir + phenoFile, "\t", log), 2);
    //// CmdLine.run("plink --bfile plink --keep ../CEU.txt --freq --out CEU", dir);
    // CmdLine.run("plink --bfile CEU --freq --out CEU", dir);
    //// CmdLine.run("plink --bfile plink --keep ../YRI.txt --freq --out YRI", dir);
    // CmdLine.run("plink --bfile YRI --freq --out YRI", dir);

    if (!new File(dir + "allSNPs.dat").exists()) {
      SnpMarkerSet.merge(new SnpMarkerSet[] {new SnpMarkerSet(dir + "CEU.bim"),
                                             new SnpMarkerSet(dir + "YRI.bim")})
                  .writeToFile(dir + "allSNPs.dat", SnpMarkerSet.GENERIC_FORMAT_ANNOTATED, log);
    }

    // int count = 1;
    parameters = new String[phenos.length * 3];
    for (int i = 0; i < phenos.length; i++) {
      // CmdLine.run("plink --bfile plink --keep ../CEU.txt --pheno "+phenoFile+" --pheno-name
      // "+phenos[i]+" --linear --ci 0.95 --sex --perm --out CEU_"+phenos[i], dir);
      // CmdLine.run("plink --bfile CEU --pheno "+phenoFile+" --pheno-name "+phenos[i]+" --linear
      // --ci 0.95 --sex --out CEU_"+phenos[i], dir);
      // Files.writeList(new String[] {"plink --bfile CEU --pheno "+phenoFile+" --pheno-name
      // "+phenos[i]+" --linear --ci 0.95 --sex --out CEU_"+phenos[i]}, dir+(count++)+".bat");
      Metal.convertPlinkResults(dir, "CEU_" + phenos[i] + ".assoc.linear", "ADD", "linear",
                                "CEU.frq", true, true, dir + "CEU_" + phenos[i] + ".se.metal",
                                false);
      // CmdLine.run("plink --bfile plink --keep ../YRI.txt --pheno "+phenoFile+" --pheno-name
      // "+phenos[i]+" --linear --ci 0.95 --sex --perm --out YRI_"+phenos[i], dir);
      // CmdLine.run("plink --bfile YRI --pheno "+phenoFile+" --pheno-name "+phenos[i]+" --linear
      // --ci 0.95 --sex --out YRI_"+phenos[i], dir);
      // Files.writeList(new String[] {"plink --bfile YRI --pheno "+phenoFile+" --pheno-name
      // "+phenos[i]+" --linear --ci 0.95 --sex --out YRI_"+phenos[i]}, dir+(count++)+".bat");
      Metal.convertPlinkResults(dir, "YRI_" + phenos[i] + ".assoc.linear", "ADD", "linear",
                                "YRI.frq", true, true, dir + "YRI_" + phenos[i] + ".se.metal",
                                false);
      Metal.metaAnalyze(dir,
                        new String[] {"CEU_" + phenos[i] + ".se.metal",
                                      "YRI_" + phenos[i] + ".se.metal"},
                        phenos[i] + "_SE", true, log);

      parameters[i * 3 + 0] = dir + phenos[i] + "_SE1.out" + " 0 5=Meta_" + phenos[i] + "_pval";
      parameters[i * 3 + 1] =
          dir + "CEU_" + phenos[i] + ".se.metal" + " 0 5=CEU_" + phenos[i] + "_pval";
      parameters[i * 3 + 2] =
          dir + "YRI_" + phenos[i] + ".se.metal" + " 0 5=YRI_" + phenos[i] + "_pval";
    }
    // Files.combineWithLessMemory(HashVec.loadFileToStringArray(dir+"allSNPs.dat", false, new int[]
    // {0}, false), parameters, null, "MarkerName", dir+"results.xln", log, true, true, false,
    // false, false);
  }
}
