package org.genvisis.expression;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.ext;
import org.genvisis.stats.Maths;
import org.genvisis.stats.ProbDist;
import org.genvisis.stats.Stats;

import com.google.common.primitives.Doubles;

public class Allelotyping {
  public static final String[] HEADER = {"Sample", "SNP", "Freq"};
  public static final String SWAP_INFO = "alleles.xln";

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD1\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD1corr\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD2\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD2corr\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\merge\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\mergeCorr\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD1fixedAgain\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD2fixedAgain\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD1RawAgain\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD2RawAgain\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\Econs1\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD3\\Plate1\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD3\\Plate2\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\Alleleotypes\\PD3\\";
    String dir = "D:\\tWork\\Expression\\Alleleotypes\\Econs2\\\\";
    String filename = "Allelotypes.dat";

    String usage = "\\n" + "expression.Allelotyping requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
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
      parse(dir, filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse(String dir, String filename) {
    BufferedReader reader;
    PrintWriter writer, summary;
    String[] line, indKeys, regionKeys, snpKeys;
    Hashtable<String, Hashtable<String, Hashtable<String, DoubleVector>>> snps =
        new Hashtable<String, Hashtable<String, Hashtable<String, DoubleVector>>>();
    Hashtable<String, Hashtable<String, DoubleVector>> regions;
    Hashtable<String, DoubleVector> inds;
    DoubleVector values, regionMeans, regionStdevs, allObs, allObsAroundZero, regionLog2ratios;
    IntVector regionCounts;
    String ind, region, snpName;
    Hashtable<String, String> lookup =
        HashVec.loadFileToHashString(dir + "snpGeneLookup.txt", false);
    double mean, stdev, z, t, val, ciLow, ciHigh;
    double[] log2ratios;
    Hashtable<String, String> swaps;
    String swapCheck;
    int countAbove, countBelow;

    swaps = HashVec.loadFileToHashString(dir + SWAP_INFO, 0, new int[] {1, 2, 3}, "\t", true);

    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      ext.checkHeader(reader.readLine().trim().split("[\\s]+"), HEADER, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        ind = line[0].substring(0, line[0].length() - 2);
        region = line[0].substring(line[0].length() - 2, line[0].length() - 1);
        snpName = line[1];
        if (snps.containsKey(snpName)) {
          regions = snps.get(snpName);
        } else {
          snps.put(snpName, regions = new Hashtable<String, Hashtable<String, DoubleVector>>());
        }
        if (regions.containsKey(region)) {
          inds = regions.get(region);
        } else {
          regions.put(region, inds = new Hashtable<String, DoubleVector>());
        }
        if (inds.containsKey(ind)) {
          values = inds.get(ind);
        } else {
          inds.put(ind, values = new DoubleVector());
        }
        if (line.length > 2 && Math.abs(Double.parseDouble(line[2])) > 0.0001) {
          swapCheck = swaps.get(snpName).split("[\\s]+")[0];
          if (swapCheck == null) {
            System.err.println("Error - trying to swap check for an unlisted SNP: " + snpName);
          } else if (swapCheck.equals("yes")) {
            values.add(1 - Double.parseDouble(line[2]));
          } else if (swapCheck.equals("no")) {
            values.add(Double.parseDouble(line[2]));
          } else {
            System.err.println("Error - invalid swap designation (expecting yes/no): " + swapCheck);
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + "allelotypes.xln"));
      writer.println("Gene\tSNP\tAllele\tRegion\tIndividual\tN\tMean\tStdev\tpval\tlog2ratio\tlog2ratio_stdev\tpval\tratio\tRegionMean\tRegionStdev\tZ-score\tpval from Zdist\tT\tratio\tpval from Tdist");
      summary = new PrintWriter(new FileWriter(dir + "allelotype_summary.xln"));
      // summary.println("gene\tsnp\tN\tproportion\tstd dev\tt-value df\tregion A\tregion B\tregion
      // C");
      // summary.println("Gene\tSNP\tAllele\tN\tProportion\tStdev\tt-value\tdf\tregion A
      // ratio\tp-value\tcountAbove\tcountBelow\tregion B
      // ratio\tp-value\tcountAbove\tcountBelow\tregion C ratio\tp-value\tcountAbove\tcountBelow");
      summary.println("Gene\tSNP\tAllele\tN\tProportion\tStdev\tt-value\tdf\tregion A ratio\tp-value\t#Above\t#Below\tregion B ratio\tp-value\t#Above\t#Below\tregion C ratio\tp-value\t#Above\t#Below");

      snpKeys = HashVec.getKeys(snps);
      for (String snpKey : snpKeys) {
        regions = snps.get(snpKey);
        regionKeys = HashVec.getKeys(regions);
        for (int j = 0; j < regionKeys.length; j++) {
          inds = regions.get(regionKeys[j]);
          indKeys = HashVec.getKeys(inds);
          regionMeans = new DoubleVector();
          regionStdevs = new DoubleVector();
          regionCounts = new IntVector();
          allObs = new DoubleVector();
          allObsAroundZero = new DoubleVector();
          regionLog2ratios = new DoubleVector();
          countAbove = 0;
          countBelow = 0;
          for (int k = 0; k < indKeys.length; k++) {
            values = inds.get(indKeys[k]);
            regionMeans.add(Array.mean(Doubles.toArray(values)));
            regionStdevs.add(Array.stdev(Doubles.toArray(values)));
            regionCounts.add(values.size());
            log2ratios = new double[values.size()];
            for (int l = 0; l < values.size(); l++) {
              val = values.elementAt(l);
              allObs.add(val);
              allObsAroundZero.add(val - 0.5);
              log2ratios[l] = Maths.log2(val / (1 - val));
            }
            regionLog2ratios.add(Array.mean(log2ratios));
            writer.println((j == 0
                            && k == 0 ? lookup.get(snpKey) + "\t" + snpKey + "\t"
                                        + (swaps.get(snpKey)
                                                .split("[\\s]+")[0].equals("yes") ? swaps.get(snpKey)
                                                                                         .split("[\\s]+")[2]
                                                                                  : swaps.get(snpKey)
                                                                                         .split("[\\s]+")[1])
                                      : "\t\t")
                           + "\t" + (k == 0 ? regionKeys[j] : "") + "\t" + indKeys[k] + "\t"
                           + values.size() + "\t" + Array.mean(Doubles.toArray(values)) + "\t"
                           + Array.stdev(Doubles.toArray(values)) + "\t=NORMSDIST(-1*ABS("
                           + Array.mean(Doubles.toArray(values)) + "-0.5)/"
                           + Array.stdev(Doubles.toArray(values)) + ")*2" + "\t"
                           + Array.mean(log2ratios) + "\t" + Array.stdev(log2ratios)
                           + "\t=NORMSDIST(-1*ABS(" + Array.mean(log2ratios) + ")/"
                           + Array.stdev(log2ratios) + ")*2" + "\t=POWER(2,"
                           + Array.mean(log2ratios) + ")");
            if (ProbDist.NormDist(Array.mean(log2ratios) / Array.stdev(log2ratios)) < 0.05) {
              if (Array.mean(log2ratios) > 0) {
                countAbove++;
              } else {
                countBelow++;
              }
            }

            // writer.println(lookup.get(snpKeys[i])+"\t"+snpKeys[i]+"\t"+regionKeys[j]+"\t"+indKeys[k]+"\t"+values.size()+"\t"+ext.formDeci(Array.mean(values.toArray()),
            // 3, true)+"\t"+ext.formDeci(Array.stdev(values.toArray()), 3, true));
          }
          mean = Array.mean(Doubles.toArray(regionMeans));
          stdev = Array.stdev(Doubles.toArray(regionMeans));
          z = Stats.ztest(mean, stdev, 0.5);
          t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0.5);
          writer.println("\t\t\t\t\t\t\t\t\t\t\t\t\t" + mean + "\t" + stdev + "\t" + z + "\t"
                         + ProbDist.NormDist(z) + "\t" + t + "\t\t"
                         + ProbDist.TDist(t, indKeys.length - 1)
                         + "\tOne observation (mean) per sample");
          if (j == 0) {
            summary.print(lookup.get(snpKey) + "\t" + snpKey + "\t"
                          + (swaps.get(snpKey)
                                  .split("[\\s]+")[0].equals("yes") ? swaps.get(snpKey)
                                                                           .split("[\\s]+")[2]
                                                                    : swaps.get(snpKey)
                                                                           .split("[\\s]+")[1])
                          + "\t" + indKeys.length + "\t" + mean + "\t" + stdev + "\t" + t + "\t"
                          + (indKeys.length - 1));
          }
          // summary.print("\t"+ext.prettyP(ProbDist.TDist(t, indKeys.length-1)));

          mean = Array.mean(Doubles.toArray(regionLog2ratios));
          stdev = Array.stdev(Doubles.toArray(regionLog2ratios));
          // ciLow = mean-1.96*stdev/Math.sqrt(regionLog2ratios.size()); // for z-test
          // ciHigh = mean+1.96*stdev/Math.sqrt(regionLog2ratios.size());
          t = ProbDist.TDistReverse(0.05, regionLog2ratios.size() - 1);
          ciLow = mean - t * stdev / Math.sqrt(regionLog2ratios.size());
          ciHigh = mean + t * stdev / Math.sqrt(regionLog2ratios.size());
          z = Stats.ztest(mean, stdev, 0);
          t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0);

          writer.println("\t\t\t\t\t\t\t\t\t\t\t\t\t" + mean + "\t" + stdev + "\t" + z + "\t"
                         + ProbDist.NormDist(z) + "\t" + t + "\t"
                         + ext.formDeci(Math.pow(2, mean), 3, true) + " ("
                         + ext.formDeci(Math.pow(2, ciLow), 3, true) + ", "
                         + ext.formDeci(Math.pow(2, ciHigh), 3, true) + ")" + "\t"
                         + ProbDist.TDist(t, indKeys.length - 1)
                         + "\tOne observation (ratio) per sample");

          summary.print("\t" + ext.formDeci(Math.pow(2, mean), 3, true) + " ("
                        + ext.formDeci(Math.pow(2, ciLow), 3, true) + ", "
                        + ext.formDeci(Math.pow(2, ciHigh), 3, true) + ")" + "\t"
                        + ext.prettyP(ProbDist.TDist(t, indKeys.length - 1)));
          summary.print("\t" + countAbove + "\t" + countBelow);

          mean = Array.mean(Doubles.toArray(regionMeans));
          stdev = Array.stdev(Doubles.toArray(allObs));
          z = Stats.ztest(mean, stdev, 0.5);
          t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0.5);
          writer.println("\t\t\t\t\t\t\t\t\t\t\t\t\t" + mean + "\t" + stdev + "\t" + z + "\t"
                         + ProbDist.NormDist(z) + "\t" + t + "\t\t"
                         + ProbDist.TDist(t, indKeys.length - 1) + "\tAll observations");

          stdev = pool(regionStdevs, regionCounts);
          z = Stats.ztest(mean, stdev, 0.5);
          t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0.5);
          writer.println("\t\t\t\t\t\t\t\t\t\t\t\t\t" + mean + "\t" + stdev + "\t" + z + "\t"
                         + ProbDist.NormDist(z) + "\t" + t + "\t\t"
                         + ProbDist.TDist(t, indKeys.length - 1)
                         + "\tWith pooled variance estimator");
        }
        summary.println();
      }
      writer.close();
      summary.close();
    } catch (Exception e) {
      System.err.println("Error writing to allelptypes.out in " + dir);
      e.printStackTrace();
    }
  }

  // pooled two sample variance estimator
  public static double pool(DoubleVector stdevs, IntVector counts) {
    double sumVars = 0;
    int totalN = 0;

    for (int i = 0; i < stdevs.size(); i++) {
      sumVars += stdevs.elementAt(i) * stdevs.elementAt(i) * (counts.elementAt(i) - 1);
      totalN += counts.elementAt(i) - 1;
    }

    return Math.sqrt(sumVars / totalN);
  }

}
