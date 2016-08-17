// -Xms1024M -Xmx1024M
package org.genvisis.nrss;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Sort;
import org.genvisis.common.Vectors;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

public class NrsHap {
  public static final double INDEX_THRESHOLD = 0.001;

  public static final double INCLUSION_THRESHOLD = 0.01;

  public static final int WINDOW = 150000;

  public static final String[][] MODELS =
      {{"0.001", "0.01", "150000"}, {"0.001", "0.01", "50000"}, {"0.001", "0.01", "25000"},
       {"0.0001", "0.01", "250000"}, {"0.0001", "0.01", "150000"}, {"0.0001", "0.01", "25000"},};

  // public static final String[][] MODELS = {{"0.01", "0.1", "150000"}};

  // public static final int EST_NUM_MARKERS_IN_LARGEST_CHR = 25000;
  public static final int EST_NUM_MARKERS_IN_LARGEST_CHR = 40000;

  public static double computeStatistic(String chr, String[] markerNames, double[] pvalues,
                                        double threshold) {
    return computeStatistic(chr, markerNames, pvalues, 0, pvalues.length - 1, threshold, "");
  }

  public static double computeStatistic(String chr, String[] markerNames, double[] pvalues,
                                        int lowIndex, int highIndex, double threshold,
                                        String extraStuff) {
    int[] keys;
    double[] pval_window;
    Vector<String> snps = new Vector<String>();
    BufferedReader reader;
    String[] line;

    pval_window = new double[highIndex - lowIndex + 1];
    for (int i = 0; i < pval_window.length; i++) {
      pval_window[i] = pvalues[lowIndex + i];
    }
    keys = Sort.quicksort(pval_window);

    for (int i = 0; i < pval_window.length; i++) {
      if (pval_window[keys[i]] < threshold) {
        snps.add(markerNames[lowIndex + keys[i]]);
      }
    }
    String temp = "plink --bfile chr" + chr + extraStuff + " --hap-assoc --snps "
                  + Array.toStr(Array.toStringArray(snps), ",") + " --hap-window " + snps.size();
    // System.out.println(temp);
    CmdLine.run(temp, ".");

    try {
      reader = new BufferedReader(new FileReader("plink.assoc.hap"));
      if (!reader.readLine().trim().split("[\\s]+")[1].equals("HAPLOTYPE")) {
        System.err.println("Error - could not find the result of the haplotype test for "
                           + Array.toStr(Array.toStringArray(snps), ","));
      }
      line = reader.readLine().trim().split("[\\s]+");
      reader.close();
      new File("plink.assoc.hap").delete();
      return Double.parseDouble(line[4]);
    } catch (Exception e) {
      e.printStackTrace();
    }

    return 0;
  }

  public static void formDistributions(int column, int levelsDeep) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp;
    int count;
    File[] files, dirs;
    DoubleVector[][] dvs;
    DoubleVector trav;
    int[] keys;
    String col;

    dirs = new File(".").listFiles(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return file.isDirectory() && filename.startsWith("index");
      }
    });

    dvs = new DoubleVector[dirs.length][];
    for (int i = 0; i < dirs.length; i++) {
      files = dirs[i].listFiles(new FilenameFilter() {
        @Override
        public boolean accept(File file, String filename) {
          return filename.endsWith("parsed_cluster.xln");
        }
      });
      dvs[i] = Vectors.initializedArray(DoubleVector.class, levelsDeep);
      for (File file : files) {
        trav = new DoubleVector();
        try {
          reader = new BufferedReader(new FileReader(file));
          reader.readLine();
          while (reader.ready()) {
            line = reader.readLine().split("[\\s]+");
            trav.add(Double.parseDouble(line[column]));
          }
          reader.close();
          keys = Sort.quicksort(trav, Sort.DESCENDING);
          for (int k = 0; k < levelsDeep && trav.size() > k; k++) {
            dvs[i][k].add(trav.elementAt(keys[k]));
          }
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + file.getName()
                             + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + file.getName() + "\"");
          System.exit(2);
        }
      }
    }

    try {
      writer = new PrintWriter(new FileWriter("NRSS_distributions.xln"));
      for (int i = 0; i < dirs.length; i++) {
        writer.print((i == 0 ? "" : "\t") + dirs[i].getName()
                     + Array.toStr(Array.stringArray(levelsDeep), "\t"));
      }
      writer.println();
      for (int i = 0; i < dirs.length; i++) {
        writer.print((i == 0 ? "" : "\t")
                     + Array.toStr(Array.stringArraySequence(levelsDeep, ""), "\t"));
      }
      writer.println();
      writer.println(Array.toStr(Array.stringArray(dirs.length * levelsDeep, "0.5")));
      temp = "";
      for (int i = 0; i < dirs.length; i++) {
        for (int j = 0; j < levelsDeep; j++) {
          col = ext.getExcelColumn(i * levelsDeep + j);
          temp += (i == 0 && j == 0 ? "" : "\t") + "=COUNTIF(" + col + "5:" + col
                  + (Math.max(4 + dvs[i][j].size(), 5)) + ", \">\"&" + col + "3)/COUNT(" + col
                  + "5:" + col + (Math.max(4 + dvs[i][j].size(), 5)) + ")";
        }
      }

      count = 0;
      while (!temp.equals(Array.toStr(Array.stringArray(dirs.length * levelsDeep)))) {
        writer.println(temp);
        temp = "";
        for (int i = 0; i < dirs.length; i++) {
          for (int j = 0; j < levelsDeep; j++) {
            if (dvs[i][j].size() > count) {
              temp += (i == 0 && j == 0 ? "" : "\t") + dvs[i][j].elementAt(count);
            } else {
              temp += (i == 0 && j == 0 ? "" : "\t");
            }
          }
        }
        count++;
      }

      writer.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "logistic.xls";
    String filename = "plink.assoc.logistic.flipped";
    // int procMap = -1;
    int simulate = -1;
    boolean parseSims = false;
    boolean procSims = false;
    boolean formDist = false;
    int col = 8;
    // int col = 6;

    String usage = "\\n" + "park.gwa.NrsHap requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n"
                   + "   (1) column of p-values (i.e. col=" + col + " (default))\n"
                   + "   (2) proc map for Haploview (i.e. procMap=4 for chromosome 4 (not the default))\n"
                   + "   (4) simulate null distribution (i.e. sim=100 to do 100 replicates (not the default))\n"
                   + "   (3) parse simulations (i.e. -parseSims (not the default))\n"
                   + "   (3) process simulation results (i.e. -procSims (not the default))\n"
                   + "   (3) form distributions (i.e. -formDist (not the default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("col=")) {
        col = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
        // } else if (args[i].startsWith("procMap=")) {
        // procMap = Integer.parseInt(args[i].split("=")[1]);
        // numArgs--;
      } else if (arg.startsWith("sim=")) {
        simulate = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-parseSims")) {
        parseSims = true;
        numArgs--;
      } else if (arg.startsWith("-procSims")) {
        procSims = true;
        numArgs--;
      } else if (arg.startsWith("-formDist")) {
        formDist = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (simulate > 0) {
        // simulateNull("pd_gwas.fam", simulate);
        simulateNull("plink.fam", simulate);
      } else if (parseSims) {
        parseSimulations("logistic", 8, 0, 2, 1);
      } else if (procSims) {
        // procSimulations("logistic.parsed", 8, 0, 2, 1);
        procSimulations("logistic.parsed", 3, 1, 2, 0);
        // procSimulations(".final", 3, 1, 2, 0);
        // procSimulations(".included", 3, 1, 2, 0);
      } else if (formDist) {
        formDistributions(3, 5);
      } else {
        // procRegion("C:\\Documents and Settings\\npankrat\\My
        // Documents\\gwas\\postHocs\\nonredundency\\SNCA.results.txt",
        // 3, "C:\\Documents and Settings\\npankrat\\My
        // Documents\\gwas\\postHocs\\nonredundency\\SNCA.recode.ped.LD");
        // procFile("chr4.xls", 6);
        // procFile("SNCA.results.txt", 3);
        // procFile("logistic.xls", 6);
        // procFile("3plus.prn", 6);
        // procFile("sim.1.additive.assoc.logistic.parsed", 3);
        // procFile("sim.1.additive.assoc.logistic.parsed", 3);
        procFile(filename, col);

      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static double minusLog(double p) {
    return -1 * Math.log10(p);
  }

  public static void parseSimulations(final String suffix, int p_column, int chr_column,
                                      int pos_column, int markerName_column) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String filename, outfile;

    File[] files = new File(".").listFiles(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith(suffix);
      }
    });
    System.out.println(ext.getTime() + "\tFound " + files.length + " files with suffix '" + suffix
                       + "' to parse");

    for (File file : files) {
      filename = file.getName();
      outfile = filename + ".parsed";

      try {
        reader = new BufferedReader(new FileReader(filename));
        writer = new PrintWriter(new FileWriter(outfile));
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          writer.println(line[markerName_column] + "\t" + line[chr_column] + "\t" + line[pos_column]
                         + "\t" + (line[p_column].equals("NA") ? "." : line[p_column]));
        }
        reader.close();
        writer.close();
      } catch (Exception e) {
        System.err.println("Error parsing '" + filename + "'");
      }
    }
    System.out.println(ext.getTime() + "\tDone");
  }

  public static void procFile(String filename, int p_column) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Vector<String> markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
    IntVector markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    DoubleVector pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    int count, lowIndex, highIndex;
    IntVector indexVector = new IntVector();
    String chr;
    double d;
    String[] markerNames;
    int[] positions;
    double[] pvals, stats;
    int[] indexSNPs, markerCounts;
    boolean done;
    IntVector cluster;
    double maxStat;
    int maxSNP;
    double minPvalue;

    try {
      System.out.println(ext.getTime());
      writer = new PrintWriter(new FileWriter("nrsHap_results.xln"));
      writer.println("Marker\tChr\tPostition\tNRS Statistic\tNumMarkers\tUCSC coordinates\tMin p-value");
      reader = new BufferedReader(new FileReader(filename));
      reader.readLine();
      line = reader.readLine().trim().split("[\\s]+");
      chr = "";
      count = 0;
      done = false;
      while (!done) {
        // System.out.println(line[0]+"\t"+line[1]);
        if (chr.equals("")) {
          chr = line[1];
          count = 0;
        } else {
          if (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
          } else {
            done = true;
            line[1] = "done";
          }
        }

        if (line[1].equals(chr)) {
          if (!line[p_column].equals(".") && !line[p_column].equals("NA")) {
            d = Double.parseDouble(line[p_column]);
            pvalueVector.add(d);
            if (d < INDEX_THRESHOLD) {
              indexVector.add(count);
            }
            markerVector.add(line[0]);
            markerLocations.add(Integer.parseInt(line[2]));
            count++;
          }
        } else {
          markerNames = Array.toStringArray(markerVector);
          positions = Ints.toArray(markerLocations);
          pvals = Doubles.toArray(pvalueVector);
          indexSNPs = Ints.toArray(indexVector);
          markerLocations.clear();
          indexVector.clear();
          markerVector.clear();
          pvalueVector.clear();

          if (!new File("chr" + chr + ".bed").exists()) {
            CmdLine.run("plink --bfile plink --chr " + chr + " --make-bed --out chr" + chr, ".");
          }
          stats = new double[indexSNPs.length];
          markerCounts = new int[indexSNPs.length];
          System.out.println("Analyzing " + indexSNPs.length + " index SNPs on chromosome " + chr);
          cluster = new IntVector();
          for (int i = 0; i < indexSNPs.length; i++) {
            lowIndex = indexSNPs[i];
            minPvalue = 1;
            while (lowIndex >= 0 && positions[lowIndex] > positions[indexSNPs[i]] - WINDOW) {
              if (pvals[lowIndex] < minPvalue) {
                minPvalue = pvals[lowIndex];
              }
              lowIndex--;
            }
            lowIndex++;
            highIndex = indexSNPs[i];
            while (highIndex < positions.length
                   && positions[highIndex] < positions[indexSNPs[i]] + WINDOW) {
              if (pvals[highIndex] < minPvalue) {
                minPvalue = pvals[highIndex];
              }
              highIndex++;
            }
            highIndex--;
            stats[i] = computeStatistic(chr, markerNames, pvals, lowIndex, highIndex,
                                        INCLUSION_THRESHOLD, "");
            markerCounts[i] = highIndex - lowIndex + 1;
            // System.out.println((i+1)+")
            // "+markerNames[indexSNPs[i]]+"\t"+positions[indexSNPs[i]]+"\t"+ext.formDeci(stats[i],
            // 2)+" (using "+(highIndex-lowIndex+1)+" markers)");

            if (i == 0) {
              cluster.add(i);
            } else if (positions[indexSNPs[i]] - positions[indexSNPs[i - 1]] < WINDOW) {
              cluster.add(i);
            } else {
              maxStat = -1;
              maxSNP = -1;
              for (int j = 0; j < cluster.size(); j++) {
                if (stats[cluster.elementAt(j)] > maxStat) {
                  maxSNP = cluster.elementAt(j);
                  maxStat = stats[maxSNP];
                }
              }
              writer.println(markerNames[indexSNPs[maxSNP]] + "\t" + chr + "\t"
                             + positions[indexSNPs[maxSNP]] + "\t" + ext.formDeci(stats[maxSNP], 2)
                             + "\t" + markerCounts[maxSNP] + "\tchr" + chr + ":"
                             + (positions[indexSNPs[maxSNP]] - WINDOW) + "-"
                             + (positions[indexSNPs[maxSNP]] + WINDOW) + "\t" + minPvalue);
              writer.flush();
              cluster.clear();
              cluster.add(i);
            }

            if (i == indexSNPs.length - 1) {
              maxStat = -1;
              maxSNP = -1;
              for (int j = 0; j < cluster.size(); j++) {
                if (stats[cluster.elementAt(j)] > maxStat) {
                  maxSNP = cluster.elementAt(j);
                  maxStat = stats[maxSNP];
                }
              }
              writer.println(markerNames[indexSNPs[maxSNP]] + "\t" + chr + "\t"
                             + positions[indexSNPs[maxSNP]] + "\t" + ext.formDeci(stats[maxSNP], 2)
                             + "\t" + markerCounts[maxSNP] + "\tchr" + chr + ":"
                             + (positions[indexSNPs[maxSNP]] - WINDOW) + "-"
                             + (positions[indexSNPs[maxSNP]] + WINDOW) + "\t" + minPvalue);
              writer.flush();
            }

          }

          markerNames = null;
          positions = null;
          pvals = null;
          indexSNPs = null;
          markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
          markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
          pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
          indexVector = new IntVector();
          chr = "";
        }

      }
      reader.close();
      writer.close();
      System.out.println(ext.getTime());
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

  }

  public static void procRegion(String results, int p_column) {
    BufferedReader reader = null;
    String[] line;
    Vector<String> markerVector = new Vector<String>();
    DoubleVector pvalueVector = new DoubleVector();
    String chr = null;

    try {
      reader = new BufferedReader(new FileReader(results));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        markerVector.add(line[0]);
        if (chr == null) {
          chr = line[1];
        } else if (!line[1].equals(chr)) {
          System.err.println("Error - different chromosomes present in file to be processed: " + chr
                             + " and " + line[1]);
        }
        pvalueVector.add(Double.parseDouble(line[p_column]));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + results + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + results + "\"");
      System.exit(2);
    }
    System.out.println("Stat is "
                       + ext.formDeci(computeStatistic(chr, Array.toStringArray(markerVector),
                                                       Doubles.toArray(pvalueVector),
                                                       INCLUSION_THRESHOLD),
                                      2));
  }

  public static void procSimulations(final String suffix, int p_column, int chr_column,
                                     int pos_column, int markerName_column) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Vector<String> markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
    IntVector markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    DoubleVector pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    int count, lowIndex, highIndex;
    IntVector indexVector = new IntVector();
    // String chr;
    double d;
    String[] markerNames;
    int[] positions;
    double[] pvals, stats;
    int[] indexSNPs, markerCounts;
    IntVector cluster;
    double maxStat;
    int maxSNP;
    String filename, outfile;
    String[] models;
    double[] indexThresholds, inclusionThresholds;
    double maxThreshold;
    int[] windowSizes;
    boolean firstOne;

    File[] files = new File(".").listFiles(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith(suffix);
      }
    });
    System.out.println("Found " + files.length + " files with suffix '" + suffix
                       + "' to compute from");

    maxThreshold = 0;
    models = new String[MODELS.length];
    indexThresholds = new double[MODELS.length];
    inclusionThresholds = new double[MODELS.length];
    windowSizes = new int[MODELS.length];
    for (int t = 0; t < MODELS.length; t++) {
      models[t] = "index" + MODELS[t][0] + "_incl" + MODELS[t][1] + "_"
                  + (Integer.parseInt(MODELS[t][2]) / 1000) + "kb/";
      new File(models[t]).mkdirs();
      indexThresholds[t] = Double.parseDouble(MODELS[t][0]);
      inclusionThresholds[t] = Double.parseDouble(MODELS[t][1]);
      windowSizes[t] = Integer.parseInt(MODELS[t][2]);
      if (indexThresholds[t] > maxThreshold) {
        maxThreshold = indexThresholds[t];
      }
    }
    System.out.println("Found " + MODELS.length + " thresholds to consider (min index pvalue of "
                       + maxThreshold + ")");

    // for (int chr_target = 1; chr_target <= 23; chr_target++) {
    for (int chr_target = 1; chr_target <= 22; chr_target++) {
      System.out.println(ext.getTime() + "\tParsing chromosome " + chr_target);
      if (!new File("chr" + chr_target + ".bed").exists()) {
        CmdLine.run("plink --bfile plink --chr " + chr_target + " --make-bed --out chr"
                    + chr_target, ".");
      }
      for (File file : files) {
        filename = file.getName();
        outfile = filename + "_cluster.xln";
        // System.out.println(filename);
        try {
          reader = new BufferedReader(new FileReader(filename));
          reader.readLine();
          count = 0;
          while (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
            if (line[chr_column].equals(chr_target + "")) {
              if (!line[p_column].equals(".") && !line[p_column].equals("NA")) {
                d = Double.parseDouble(line[p_column]);
                pvalueVector.add(d);
                if (d < maxThreshold) {
                  indexVector.add(count);
                }
                markerVector.add(line[markerName_column]);
                markerLocations.add(Integer.parseInt(line[pos_column]));
                count++;
              }
            }
          }
          reader.close();
          markerNames = Array.toStringArray(markerVector);
          positions = Ints.toArray(markerLocations);
          pvals = Doubles.toArray(pvalueVector);
          indexSNPs = Ints.toArray(indexVector);
          markerLocations.clear();
          indexVector.clear();
          markerVector.clear();
          pvalueVector.clear();

          for (int t = 0; t < MODELS.length; t++) {
            writer = new PrintWriter(new FileWriter(models[t] + outfile, chr_target > 1));
            if (chr_target == 1) {
              writer.println("Marker\tChr\tPostition\tNRS Statistic\tNumMarkers\tUCSC coordinates");
            }

            stats = new double[indexSNPs.length];
            markerCounts = new int[indexSNPs.length];
            //
            // System.out.println("Analyzing "+indexSNPs.length+"
            // index SNPs on chromosome "+chr_target);

            cluster = new IntVector();
            firstOne = true;
            for (int i = 0; i < indexSNPs.length; i++) {
              if (indexSNPs[i] == -1) {
                System.err.println("Error - model " + t + " indexSNP " + i);
              }
              if (pvals[indexSNPs[i]] < indexThresholds[t]) {
                lowIndex = indexSNPs[i];
                while (lowIndex >= 0
                       && positions[lowIndex] > positions[indexSNPs[i]] - windowSizes[t]) {
                  lowIndex--;
                }
                lowIndex++;
                highIndex = indexSNPs[i];
                while (highIndex < positions.length
                       && positions[highIndex] < positions[indexSNPs[i]] + windowSizes[t]) {
                  highIndex++;
                }
                highIndex--;
                stats[i] = computeStatistic(chr_target + "", markerNames, pvals, lowIndex,
                                            highIndex, inclusionThresholds[t],
                                            " --pheno sim." + filename.split("\\.")[1] + ".dat");
                markerCounts[i] = highIndex - lowIndex + 1;
                //
                // System.out.println((i+1)+")
                // "+markerNames[indexSNPs[i]]+"\t"+positions[indexSNPs[i]]+"\t"+ext.formDeci(stats[i],
                // 2)+" (using "+(highIndex-lowIndex+1)+"
                // markers)");

                if (firstOne) {
                  cluster.add(i);
                  firstOne = false;
                } else if (positions[indexSNPs[i]] - positions[indexSNPs[i - 1]] < windowSizes[t]) {
                  cluster.add(i);
                } else {
                  maxStat = 0;
                  maxSNP = -1;
                  for (int j = 0; j < cluster.size(); j++) {
                    if (stats[cluster.elementAt(j)] > maxStat) {
                      maxSNP = cluster.elementAt(j);
                      maxStat = stats[maxSNP];
                    }
                  }
                  writer.println(markerNames[indexSNPs[maxSNP]] + "\t" + chr_target + "\t"
                                 + positions[indexSNPs[maxSNP]] + "\t"
                                 + ext.formDeci(stats[maxSNP], 2) + "\t" + markerCounts[maxSNP]
                                 + "\tchr" + chr_target + ":"
                                 + (positions[indexSNPs[maxSNP]] - windowSizes[t]) + "-"
                                 + (positions[indexSNPs[maxSNP]] + windowSizes[t]));
                  writer.flush();
                  cluster.clear();
                  cluster.add(i);
                }

                if (i == indexSNPs.length - 1) {
                  maxStat = 0;
                  maxSNP = -1;
                  for (int j = 0; j < cluster.size(); j++) {
                    if (stats[cluster.elementAt(j)] > maxStat) {
                      maxSNP = cluster.elementAt(j);
                      maxStat = stats[maxSNP];
                    }
                  }
                  writer.println(markerNames[indexSNPs[maxSNP]] + "\t" + chr_target + "\t"
                                 + positions[indexSNPs[maxSNP]] + "\t"
                                 + ext.formDeci(stats[maxSNP], 2) + "\t" + markerCounts[maxSNP]
                                 + "\tchr" + chr_target + ":"
                                 + (positions[indexSNPs[maxSNP]] - windowSizes[t]) + "-"
                                 + (positions[indexSNPs[maxSNP]] + windowSizes[t]));
                  writer.flush();
                }

              }
            }

            writer.close();
          }

          markerNames = null;
          positions = null;
          pvals = null;
          indexSNPs = null;
          markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
          markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
          pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
          indexVector = new IntVector();

        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + filename + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + filename + "\"");
          System.exit(2);
        }
      }
    }
  }

  public static void simulateNull(String pedfile, int replicates) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Vector<String> v = new Vector<String>();
    String[][] data;
    int[] order;
    int count;

    try {
      reader = new BufferedReader(new FileReader(pedfile));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        v.add(line[0] + "\t" + line[1] + "\t" + line[5]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - could not find " + pedfile + " in current directory");
      System.exit(2);
    } catch (IOException ioe) {
      System.err.println("Error parsing " + pedfile + "");
      System.exit(3);
    }

    data = new String[v.size()][];
    for (int i = 0; i < v.size(); i++) {
      data[i] = v.elementAt(i).split("[\\s]+");
    }

    for (int i = 1; i <= replicates; i++) {
      try {
        writer = new PrintWriter(new FileWriter("sim." + i + ".dat"));
        writer.println("FID\tIID\tAff");
        order = Array.random(data.length);
        for (int j = 0; j < data.length; j++) {
          writer.println(data[j][0] + "\t" + data[j][1] + "\t" + data[order[j]][2]);
        }
        writer.close();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }

    try {
      writer = null;
      count = 0;
      for (int i = 0; i < replicates; i++) {
        if (i % 250 == 0) {
          if (writer != null) {
            writer.close();
            Files.chmod("batchSims." + count);
          }
          count++;
          writer = new PrintWriter(new FileWriter("batchSims." + count));
        }

        // writer.println("plink --bfile pd_gwas --pheno
        // sim."+(i+1)+".dat --logistic --out sim."+(i+1)+".additive");
        writer.println("plink --bfile plink --pheno sim." + (i + 1) + ".dat --logistic --out sim."
                       + (i + 1) + ".additive");
      }
      writer.close();
      Files.chmod("batchSims." + count);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
