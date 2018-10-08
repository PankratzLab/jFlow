// -Xms1024M -Xmx1024M
package org.genvisis.nrss;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;
import org.genvisis.seq.biostatistics.MapSNPsAndGenes;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.DoubleVector;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.Vectors;
import org.pankratzlab.common.ext;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

public class OriginalNrss {

  public static final double INDEX_THRESHOLD = 0.001;
  public static final double INCLUSION_THRESHOLD = 0.01;
  public static final int DEFAULT_WINDOW = 150000;

  // public static final String[][] MODELS = {{"0.001", "0.01", "150000"},
  // {"0.001", "0.01", "50000"},
  // {"0.001", "0.01", "25000"},
  // {"0.0001", "0.01", "250000"},
  // {"0.0001", "0.01", "150000"},
  // {"0.0001", "0.01", "25000"},
  // };

  // public static final String[][] MODELS = {{"0.01", "0.1", "150000"}};

  public static final String[][] MODELS = {{"0.001", "0.01", "150000"}};

  public static final int EST_NUM_MARKERS_IN_LARGEST_CHR = 25000;

  public static final String DEFAULT_LD_LOC = "";

  // public static final String DEFAULT_LD_LOC = "C:\\Documents and
  // Settings\\npankrat\\My Documents\\gwas\\udall\\LD\\";

  public static void procFile(String dir, String filename, int p_column, int variate_column,
                              int window, String outfile) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash;
    Vector<String> markerVector = new Vector<>(EST_NUM_MARKERS_IN_LARGEST_CHR);
    IntVector markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    DoubleVector pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    DoubleVector varVector = variate_column == -1 ? null
                                                  : new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
    int count, lowIndex, highIndex;
    IntVector indexVector = new IntVector();
    String chr;
    double d;
    String[] markerNames;
    int[] positions;
    int[][] markerPositions;
    double[] pvals, vars = null;
    double[][] stats;
    int[] indexSNPs, markerCounts;
    boolean done;
    IntVector cluster;
    double maxStat;
    int maxSNP;
    double minPvalue;
    Vector<String> v = new Vector<>();
    String[] genes;

    try {
      System.out.println(ext.getTime());
      reader = new BufferedReader(new FileReader(dir + filename));
      reader.readLine();
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      chr = "";
      count = 0;
      done = false;
      while (!done) {
        if (chr.equals("")) {
          chr = line[1];
          count = 0;
        } else {
          if (reader.ready()) {
            line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          } else {
            done = true;
            line[1] = "done";
          }
        }

        if (line[1].equals(chr)) {
          if (!line[p_column].equals(".") && !line[p_column].equals("NA")) {
            d = Double.parseDouble(line[p_column]);
            pvalueVector.add(d);
            if (varVector != null) {
              varVector.add(Double.parseDouble(line[variate_column]));
            }
            if (d < INDEX_THRESHOLD) {
              indexVector.add(count);
            }
            markerVector.add(line[0]);
            markerLocations.add(Integer.parseInt(line[2]));
            count++;
          }
        } else {
          markerNames = ArrayUtils.toStringArray(markerVector);
          positions = Ints.toArray(markerLocations);
          pvals = Doubles.toArray(pvalueVector);
          indexSNPs = Ints.toArray(indexVector);
          markerLocations.clear();
          indexVector.clear();
          markerVector.clear();
          pvalueVector.clear();
          if (varVector != null) {
            vars = Doubles.toArray(varVector);
            varVector.clear();
          }

          stats = new double[indexSNPs.length][];
          markerCounts = new int[indexSNPs.length];
          System.out.println("Analyzing " + indexSNPs.length + " index SNPs on chromosome " + chr);
          if (indexSNPs.length > 0) {
            // hash = new Hashtable<String, String>();
            hash = loadLD(DEFAULT_LD_LOC + "chr" + chr + ".recode.ped.LD");
          } else {
            hash = null;
          }
          cluster = new IntVector();

          for (int i = 0; i < indexSNPs.length; i++) {
            lowIndex = indexSNPs[i];
            minPvalue = 1;
            while (lowIndex >= 0 && positions[lowIndex] > positions[indexSNPs[i]] - window) {
              if (pvals[lowIndex] < minPvalue) {
                minPvalue = pvals[lowIndex];
              }
              lowIndex--;
            }
            lowIndex++;
            highIndex = indexSNPs[i];
            while (highIndex < positions.length
                   && positions[highIndex] < positions[indexSNPs[i]] + window) {
              if (pvals[highIndex] < minPvalue) {
                minPvalue = pvals[highIndex];
              }
              highIndex++;
            }
            highIndex--;
            stats[i] = computeStatistic(markerNames, pvals, lowIndex, highIndex, hash,
                                        INCLUSION_THRESHOLD, vars);
            markerCounts[i] = highIndex - lowIndex + 1;
            // System.out.println((i+1)+")
            // "+markerNames[indexSNPs[i]]+"\t"+positions[indexSNPs[i]]+"\t"+ext.formDeci(stats[i],
            // 2)+" (using "+(highIndex-lowIndex+1)+" markers)");

            if (i == 0) {
              cluster.add(i);
            } else if (positions[indexSNPs[i]] - positions[indexSNPs[i - 1]] < window) {
              cluster.add(i);
            } else {
              maxStat = 0;
              maxSNP = -1;
              for (int j = 0; j < cluster.size(); j++) {
                if (stats[cluster.elementAt(j)][0] > maxStat) {
                  maxSNP = cluster.elementAt(j);
                  maxStat = stats[maxSNP][0];
                }
              }
              v.add(markerNames[indexSNPs[maxSNP]] + "\t" + chr + "\t"
                    + positions[indexSNPs[maxSNP]] + "\t" + ext.formDeci(stats[maxSNP][0], 2) + "\t"
                    + markerCounts[maxSNP] + "\t"
                    + ext.formDeci(stats[maxSNP][0] * Math.sqrt(markerCounts[maxSNP]), 2) + "\tchr"
                    + chr + ":" + (positions[indexSNPs[maxSNP]] - window) + "-"
                    + (positions[indexSNPs[maxSNP]] + window) + "\t" + minPvalue
                    + (variate_column == -1 ? ""
                                            : "\t"
                                              + ArrayUtils.toStr(ArrayUtils.subArray(stats[maxSNP],
                                                                                     1,
                                                                                     stats[maxSNP].length))));
              cluster.clear();
              cluster.add(i);
            }

            if (i == indexSNPs.length - 1) {
              maxStat = 0;
              maxSNP = -1;
              for (int j = 0; j < cluster.size(); j++) {
                if (stats[cluster.elementAt(j)][0] > maxStat) {
                  maxSNP = cluster.elementAt(j);
                  maxStat = stats[maxSNP][0];
                }
              }
              v.add(markerNames[indexSNPs[maxSNP]] + "\t" + chr + "\t"
                    + positions[indexSNPs[maxSNP]] + "\t" + ext.formDeci(stats[maxSNP][0], 2) + "\t"
                    + markerCounts[maxSNP] + "\t"
                    + ext.formDeci(stats[maxSNP][0] * Math.sqrt(markerCounts[maxSNP]), 2) + "\tchr"
                    + chr + ":" + (positions[indexSNPs[maxSNP]] - window) + "-"
                    + (positions[indexSNPs[maxSNP]] + window) + "\t" + minPvalue
                    + (variate_column == -1 ? ""
                                            : "\t"
                                              + ArrayUtils.toStr(ArrayUtils.subArray(stats[maxSNP],
                                                                                     1,
                                                                                     stats[maxSNP].length))));
            }
          }

          markerNames = null;
          positions = null;
          pvals = null;
          indexSNPs = null;
          markerVector = new Vector<>(EST_NUM_MARKERS_IN_LARGEST_CHR);
          markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
          pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
          indexVector = new IntVector();
          chr = "";
        }
      }
      reader.close();

      markerPositions = new int[v.size()][2];
      for (int i = 0; i < v.size(); i++) {
        line = v.elementAt(i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
        markerPositions[i][0] = Integer.parseInt(line[1]);
        markerPositions[i][1] = Integer.parseInt(line[2]);
      }
      genes = MapSNPsAndGenes.mapSNPsToGenesLoosely(markerPositions, window, 36, new Logger());

      writer = Files.openAppropriateWriter(dir + outfile);
      writer.println("Marker\tChr\tPostition\tGenes(s)\tWeightedStatistic\tNumMarkers\tUnweightedStatistic\tUCSC coordinates\tMin p-value"
                     + (variate_column == -1 ? "" : "\tvars"));
      for (int i = 0; i < v.size(); i++) {
        line = v.elementAt(i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
        writer.println(ArrayUtils.toStr(ArrayUtils.subArray(line, 0, 3)) + "\t" + genes[i] + "\t"
                       + ArrayUtils.toStr(ArrayUtils.subArray(line, 3)));
      }
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

  public static Hashtable<String, String> loadLD(String filename) {
    BufferedReader reader;
    String[] line;
    Hashtable<String, String> hash = new Hashtable<>();

    try {
      reader = new BufferedReader(new FileReader(filename));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        hash.put(line[0] + ":" + line[1], line[4]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    return hash;
  }

  public static void procRegion(String results, int p_column, String ld) {
    BufferedReader reader = null;
    String[] line;
    Hashtable<String, String> hash;
    Vector<String> markerVector = new Vector<>();
    DoubleVector pvalueVector = new DoubleVector();

    hash = loadLD(ld);

    try {
      reader = new BufferedReader(new FileReader(results));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        markerVector.add(line[0]);
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
                       + ext.formDeci(computeStatistic(ArrayUtils.toStringArray(markerVector),
                                                       Doubles.toArray(pvalueVector), hash,
                                                       INCLUSION_THRESHOLD)[0],
                                      2));
  }

  public static double[] computeStatistic(String[] markerNames, double[] pvalues,
                                          Hashtable<String, String> hash, double threshold) {
    return computeStatistic(markerNames, pvalues, 0, pvalues.length - 1, hash, threshold, null);
  }

  public static double[] computeStatistic(String[] markerNames, double[] pvalues, int lowIndex,
                                          int highIndex, Hashtable<String, String> hash,
                                          double threshold, double[] vars) {
    double[] stats = new double[5];
    int[] keys;
    double maxR2, r2;
    double[] pval_window, var_window = null;
    int count;

    if (vars != null) {
      var_window = new double[highIndex - lowIndex + 1];
      for (int i = 0; i < var_window.length; i++) {
        var_window[i] = vars[lowIndex + i];
      }
    }

    pval_window = new double[highIndex - lowIndex + 1];
    for (int i = 0; i < pval_window.length; i++) {
      pval_window[i] = pvalues[lowIndex + i];
    }
    keys = Sort.getSortedIndices(pval_window);

    count = 0;
    stats[2] = 1;
    for (int i = 0; i < pval_window.length; i++) {
      if (pval_window[keys[i]] < threshold) {
        maxR2 = 0;
        for (int j = 0; j < i; j++) {
          if (hash.containsKey(markerNames[lowIndex + keys[i]] + ":"
                               + markerNames[lowIndex + keys[j]])) {
            r2 = Double.parseDouble(hash.get(markerNames[lowIndex + keys[i]] + ":"
                                             + markerNames[lowIndex + keys[j]]));
          } else if (hash.containsKey(markerNames[lowIndex + keys[j]] + ":"
                                      + markerNames[lowIndex + keys[i]])) {
            r2 = Double.parseDouble(hash.get(markerNames[lowIndex + keys[j]] + ":"
                                             + markerNames[lowIndex + keys[i]]));
          } else {
            // System.err.println("Error - no LD computed between
            // "+markerNames[lowIndex+keys[j]]+" and
            // "+markerNames[lowIndex+keys[i]]);
            r2 = 0;
          }
          if (r2 > maxR2) {
            maxR2 = r2;
          }
        }
        stats[0] += (1 - maxR2) * minusLog(pval_window[keys[i]]);
        if (vars != null) {
          stats[1] += var_window[keys[i]];
          stats[2] = Math.min(var_window[keys[i]], stats[2]);
          if (var_window[keys[i]] < 0.05) {
            stats[3]++;
          }
          if (var_window[keys[i]] < 0.10) {
            stats[4]++;
          }
          count++;
        }
      }
    }

    stats[0] /= Math.sqrt(pval_window.length);
    if (vars != null) {
      stats[1] /= count;
      stats[3] /= count;
      stats[4] /= count;
    }

    return stats;
  }

  public static double minusLog(double p) {
    return -1 * Math.log10(p);
  }

  public static void batchHaploview(int numBatches) {
    String commands = "";

    // commands += "plink --bfile pd_gwas --chr [%0] --out chr[%0] --recode
    // --remove WGAs.txt\n";
    commands += "plink --bfile plink --chr [%0] --out chr[%0].recode --recode\n\n";
    commands += Files.getRunString() + " nrss.Nrss procMap=[%0]\n";
    commands += "java -Dsun.java2d.noddraw=true -Xmx1024m -classpath /home/npankrat/Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -log -pedfile chr[%0].recode.ped -info chr[%0].map -skipcheck -hwcutoff 0 -maxDistance 500 -dprime\n\n";

    Files.batchIt("haplo", null, numBatches, commands, ArrayUtils.stringArraySequence(23, ""));
  }

  public static void procMapForHaploview(String in, String out) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    try {
      reader = new BufferedReader(new FileReader(in));
      writer = Files.openAppropriateWriter(out);
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        writer.println(line[1] + "\t" + line[3]);
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + in + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + in + "\"");
      System.exit(2);
    }
  }

  public static void simulateNull(String pedfile, int replicates) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Vector<String> v = new Vector<>();
    String[][] data;
    int[] order;
    int count;

    try {
      reader = new BufferedReader(new FileReader(pedfile));
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
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
      data[i] = v.elementAt(i).split(PSF.Regex.GREEDY_WHITESPACE);
    }

    for (int i = 1; i <= replicates; i++) {
      try {
        writer = Files.openAppropriateWriter("sim." + i + ".dat");
        writer.println("FID\tIID\tAff");
        order = ArrayUtils.random(data.length);
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
          writer = Files.openAppropriateWriter("batchSims." + count);
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
        writer = Files.openAppropriateWriter(outfile);
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
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

  public static void procSimulations(final String suffix, int p_column, int chr_column,
                                     int pos_column, int markerName_column) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash;
    Vector<String> markerVector = new Vector<>(EST_NUM_MARKERS_IN_LARGEST_CHR);
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

    for (int chr_target = 1; chr_target <= 23; chr_target++) {
      System.out.println(ext.getTime() + "\tLoading LD information for chromosome " + chr_target);
      hash = loadLD(DEFAULT_LD_LOC + "chr" + chr_target + ".recode.ped.LD");
      // hash = new Hashtable<String, String>();

      System.out.println(ext.getTime() + "\tparsing chromosome " + chr_target);

      for (File file : files) {
        filename = file.getName();
        outfile = filename + "_cluster.xln";
        // System.out.println(filename);

        try {
          reader = new BufferedReader(new FileReader(filename));
          reader.readLine();
          count = 0;
          while (reader.ready()) {
            line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
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
          markerNames = ArrayUtils.toStringArray(markerVector);
          positions = Ints.toArray(markerLocations);
          pvals = Doubles.toArray(pvalueVector);
          indexSNPs = Ints.toArray(indexVector);
          markerLocations = null;
          indexVector = null;
          markerVector = null;
          pvalueVector = null;

          for (int t = 0; t < MODELS.length; t++) {
            writer = Files.openAppropriateWriter(models[t] + outfile, chr_target > 1);
            if (chr_target == 1) {
              writer.println("Marker\tChr\tPostition\tWeightedStatistic\tNumMarkers\tUnweightedStatistic\tUCSC coordinates");
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
                stats[i] = computeStatistic(markerNames, pvals, lowIndex, highIndex, hash,
                                            inclusionThresholds[t], null)[0];
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
                                 + "\t"
                                 + ext.formDeci(stats[maxSNP] * Math.sqrt(markerCounts[maxSNP]), 2)
                                 + "\tchr" + chr_target + ":"
                                 + (positions[indexSNPs[maxSNP]] - windowSizes[t]) + "-"
                                 + (positions[indexSNPs[maxSNP]] + windowSizes[t]));
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
                                 + "\t"
                                 + ext.formDeci(stats[maxSNP] * Math.sqrt(markerCounts[maxSNP]), 2)
                                 + "\tchr" + chr_target + ":"
                                 + (positions[indexSNPs[maxSNP]] - windowSizes[t]) + "-"
                                 + (positions[indexSNPs[maxSNP]] + windowSizes[t]));
                }

              }
            }

            writer.close();
          }

          markerNames = null;
          positions = null;
          pvals = null;
          indexSNPs = null;
          markerVector = new Vector<>(EST_NUM_MARKERS_IN_LARGEST_CHR);
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

  public static void formDistributions(int column, int levelsDeep) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp;
    int count;
    File[] files, dirs;
    DoubleVector[][] dvs;
    DoubleVector trav;
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
            line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
            trav.add(Double.parseDouble(line[column]));
          }
          reader.close();
          Collections.sort(trav, Collections.reverseOrder());
          for (int k = 0; k < levelsDeep && trav.size() > k; k++) {
            dvs[i][k].add(trav.elementAt(k));
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
      writer = Files.openAppropriateWriter("NRSS_distributions.xln");
      for (int i = 0; i < dirs.length; i++) {
        writer.print((i == 0 ? "" : "\t") + dirs[i].getName()
                     + ArrayUtils.toStr(ArrayUtils.stringArray(levelsDeep), "\t"));
      }
      writer.println();
      for (int i = 0; i < dirs.length; i++) {
        writer.print((i == 0 ? "" : "\t")
                     + ArrayUtils.toStr(ArrayUtils.stringArraySequence(levelsDeep, ""), "\t"));
      }
      writer.println();
      writer.println(ArrayUtils.toStr(ArrayUtils.stringArray(dirs.length * levelsDeep, "0.5")));
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
      while (!temp.equals(ArrayUtils.toStr(ArrayUtils.stringArray(dirs.length * levelsDeep)))) {
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
    int batch = 0;
    int procMap = -1;
    int simulate = -1;
    boolean parseSims = false;
    boolean procSims = false;
    boolean formDist = false;

    // String filename = "plink.assoc.logistic.flipped";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\_GAW16\\results\\goodMarkers\\replicate\\";
    // String filename = "additive.included";
    // int col = 8;

    // String dir = "";
    // String filename = "additive.mafs.prn";
    // int col = 3;
    // int var = 4;

    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\gwas\\udall\\";
    // String filename = "MetalResults.xln";
    // int col = 7;
    // int var = -1;

    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Renewal\\progress report\\2q NRSS\\";
    String filename = "SFH2q_final_perm.assoc.logistic.xln";
    int col = 4;
    int var = -1;

    String usage = "\\n" + "park.gwa.Nrss requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n" + "   (1) column of p-values (i.e. col=" + col
                   + " (default))\n" + "   (1) column of variate (i.e. var=" + var + " (default))\n"
                   + "   (2) proc map for Haploview (i.e. procMap=4 for chromosome 4 (not the default))\n"
                   + "   (3) batch haploview LD computation, using N files (i.e. batch=N (not the default))\n"
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
      } else if (arg.startsWith("var=")) {
        var = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("procMap=")) {
        procMap = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("sim=")) {
        simulate = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("batch=")) {
        batch = Integer.parseInt(arg.split("=")[1]);
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
      if (batch > 0) {
        batchHaploview(batch);
      } else if (procMap > 0) {
        procMapForHaploview("chr" + procMap + ".recode.map", "chr" + procMap + ".map");
      } else if (simulate > 0) {
        // simulateNull("pd_gwas.fam", simulate);
        simulateNull("plink.fam", simulate);
      } else if (parseSims) {
        parseSimulations("logistic", 8, 0, 2, 1);
      } else if (procSims) {
        // procSimulations("logistic.parsed", 8, 0, 2, 1);
        // procSimulations("logistic.parsed", 3, 1, 2, 0);
        // procSimulations(".final", 3, 1, 2, 0);
        procSimulations(".included", 3, 1, 2, 0);
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

        // procFile(dir, filename, 7, var, DEFAULT_WINDOW,
        // "IU_Add_Nrss.xln");
        // procFile(dir, filename, 10, var, DEFAULT_WINDOW,
        // "Udall_Add_Nrss.xln");
        // procFile(dir, filename, 11, var, DEFAULT_WINDOW,
        // "Meta_Add_Nrss.xln");
        // procFile(dir, filename, 13, var, DEFAULT_WINDOW,
        // "Meta_Dom_Nrss.xln");
        // procFile(dir, filename, 15, var, DEFAULT_WINDOW,
        // "Meta_Rec_Nrss.xln");

        procFile(dir, filename, col, var, DEFAULT_WINDOW, "2q_Add.xln");
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
