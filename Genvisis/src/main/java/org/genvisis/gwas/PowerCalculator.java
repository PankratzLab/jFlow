package org.genvisis.gwas;

import java.util.Hashtable;
import java.util.Map;

import org.genvisis.common.Array;
import org.genvisis.common.Internat;
import org.genvisis.common.ext;

public class PowerCalculator {
  // public static final double[] MAFs = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
  // 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
  public static final double[] MAFs =
      {0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
  public static final double[] RELATIVE_RISKS =
      {1.10, 1.20, 1.30, 1.40, 1.60, 1.80, 2.00, 2.2, 2.4, 2.6};
  // public static final double[] MAFs = {0.20};
  public static final String[] FORMATTING_TO_REMOVE = {"<em><font color=\"navy\">", "</font></em>"};
  public static final double RR_INCREMENT = 0.01;

  public static double getRelativeRiskAtEightyPercentPower(double prevalence, double maf,
                                                           int numCases, int numControls,
                                                           double alpha,
                                                           boolean unselected) throws Exception {
    boolean found;
    int index, prev;
    int[] array;
    double rr;


    found = false;
    index = 20;
    array = Array.intArray(10000, -1);
    while (!found) {
      rr = 1 + index * RR_INCREMENT;
      if (array[index] == -1) {
        array[index] =
            getSampleSize(prevalence, rr, maf, numCases, numControls, alpha, unselected, false);
        if (array[index] == -9) {
          array[index] =
              getSampleSize(prevalence, rr, maf, numCases, numControls, alpha, unselected, true);
        }
        // System.err.println("array["+index+"]: "+array[index]);
      } else if (array[index] == -9) {
        return -9;
      } else if (array[index] <= numCases && array[index - 1] >= numCases) {
        return rr;
      } else if (array[index] > numCases) {
        prev = index;
        do {
          prev++;
        } while (prev - index < 20 && array[prev] == -1);
        index = (int) Math.ceil((prev - index) / 2.0) + index;
      } else {
        prev = index;
        do {
          prev--;
          // System.err.println("Error - how did this happen at "+prev+": "+array[prev] +"<"+
          // numCases);
        } while (array[prev] == -1 && prev > 0);
        index = (int) Math.floor((index - prev) / 2.0) + prev;
      }
    }

    return -7;
  }

  public static int getSampleSize(double prevalence, double relativeRisk, double maf, int numCases,
                                  int numControls, double alpha, boolean unselected,
                                  boolean dominance) {
    String[] results, line, cells;
    double ccratio;
    String trav;
    double[] powers;
    int[] sampleSizes;

    ccratio = (double) numControls / (double) numCases;

    if (unselected) {
      System.err.println("Warning - are your controls really unselected?");
    }

    if (alpha < 1E-8) {
      System.err.println("Error - cannot set alpha to less than 1E-8; truncating");
      alpha = 1E-8;
    }

    Map<String, String> data = new Hashtable<String, String>();
    data.put("fA", ext.prettyP(maf));
    data.put("k", ext.formDeci(prevalence, 3));
    data.put("rAa", ext.formDeci(relativeRisk, 3));
    // data.put("rAA", ext.formDeci(relativeRisk*relativeRisk, 3));
    if (dominance) {
      data.put("rAA", ext.formDeci(0.99 / prevalence, 3));
    } else {
      data.put("rAA", ext.formDeci(relativeRisk + relativeRisk - 1, 3));
    }
    data.put("dprime", "1.0");
    data.put("m1", ext.prettyP(maf));
    data.put("n", numCases + "");
    data.put("ccratio", ext.formDeci(ccratio, 4, false));
    data.put("alpha", ext.prettyP(alpha, 2, 100, 2, true));
    data.put("power", "0.80");
    if (unselected) {
      data.put("unsel", "TRUE");
    }

    try {
      results =
          Internat.doSubmit("http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi", data, 1000);
    } catch (Exception e) {
      System.err.println("Error - failed to connect to website");
      e.printStackTrace();
      System.exit(1);
      return -999;
    }

    if (results[0].contains("allelic")) {
      trav = results[0].substring(results[0].indexOf("allelic"));
      for (String element : FORMATTING_TO_REMOVE) {
        trav = ext.replaceAllWith(trav, element, "");
      }
      line = trav.split("<tr>");

      powers = new double[5];
      sampleSizes = new int[5];
      for (int i = 0; i < 5; i++) {
        cells = line[i + 2].trim().split("<td>");
        for (int j = 0; j < 3; j++) {
          powers[i] = Double.parseDouble(cells[2].substring(0, cells[2].indexOf("<")).trim());
          sampleSizes[i] = Integer.parseInt(cells[3].substring(0, cells[3].indexOf("<")).trim());
        }
        // System.out.println(powers[i]+"\t"+sampleSizes[i]);
      }
      // System.out.println("\t\t\t\t\t\t\t"+relativeRisk+"\t"+powers[4]+"\t"+sampleSizes[4]);
      return sampleSizes[4];
    }

    return -9;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "PowerCalculator.dat";

    String usage = "\n" + "gwas.PowerCalculator requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      // rangeOfMaf(0.15, 0.01, 200, 200, 6, false); // diabetes
      // rangeOfMaf(0.03, 0.01, 250, 250, 100, false); // HB
      // rangeOfMaf(0.001, 0.01, 252, 871, 500000, false); // OS
      // rangeOfMaf(0.001, 0.01, 21, 120, 1000, false); // OS
      // rangeOfMaf(0.001, 0.01, 273, 991, 500000, false); // OS
      // rangeOfMaf(0.01, 0.01, 400, 400, 100000, false); // Diabetes
      // rangeOfMaf(0.15, 0.01, 882, 838, 118, false); // Indian follow up v1 all samples
      // rangeOfMaf(0.15, 0.01, 441, 438, 118, false); // Indian follow up v1 early onset cases
      // rangeOfMaf(0.001, 0.01, 285, 1500, 500, false); // Logan's ALL grant v1
      // rangeOfMaf(0.001, 0.01, 1500, 1500, 1000, false); // Logan's ALL grant v1
      // rangeOfMaf(0.001, 0.01, 1350, 1350, 1000, false); // Logan's ALL grant v2
      // rangeOfMaf(0.001, 0.01, 365, 1350, 500, false); // Logan's ALL grant v2
      // rangeOfMaf(0.001, 0.01, 985, 1350, 500, false); // Logan's ALL grant v2
      // rangeOfMaf(0.001, 0.01, 1350, 1350, 1, false); // Logan's ALL grant v2
      // rangeOfMaf(0.001, 0.01, 365, 1350, 1, false); // Logan's ALL grant v2
      // rangeOfMaf(0.001, 0.01, 985, 1350, 1, false); // Logan's ALL grant v2
      // rangeOfMaf(0.001, 0.01, 1700, 1700, 500, false); // Heather's OSCA2 grant
      // rangeOfMaf(0.001, 0.01, 1700, 1700, 2291, false); // Heather's OSCA2 grant
      // rangeOfMaf(0.001, 0.01, 1700, 1700, 2000, false); // Heather's OSCA2 grant, aim 1 discovery
      // rangeOfMaf(0.001, 0.01, 1500, 1500, 1, false); // Heather's OSCA2 grant, aim 1 replication
      // rangeOfMaf(0.001, 0.01, 1500, 1500, 5, false); // Heather's OSCA2 grant, aim 1 replication
      // rangeOfMaf(0.5, 0.01, 178, 178, 3, false); // Heather's OSCA2 grant, aim 2a
      // rangeOfMaf(0.19, 0.01, 266, 1400-266, 3, false); // Heather's OSCA2 grant, aim 2b i
      // rangeOfMaf(0.19, 0.01, (int)(1230*0.19), (int)(1230*(1-0.19)), 3, false); // Heather's
      // OSCA2 grant, aim 2b ii
      // rangeOfMaf(0.001, 0.01, 900, 800, 3, false); // Heather's OSCA2 grant, aim 3
      // rangeOfMaf(0.001, 0.01, 465, 1119, 66, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 434, 1119, 66, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 465+434, 1119, 66, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 465, 1119, 5043, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 434, 1119, 5043, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 465+434, 1119, 5043, false); // Poynter's MDS-AML grant, MDS
      rangeOfMaf(0.001, 0.01, 465, 1119, 1000000, false); // Poynter's MDS-AML grant, MDS
      rangeOfMaf(0.001, 0.01, 434, 1119, 1000000, false); // Poynter's MDS-AML grant, MDS
      rangeOfMaf(0.001, 0.01, 465 + 434, 1119, 1000000, false); // Poynter's MDS-AML grant, MDS



      // rangeOfRelativeRisk(0.15, 200, false);
      // getSampleSize();
      // getSampleSizeForASetOfPairings("D:/Myron/Indian_Diabetes/SequencingPilot/power.input");
      // getSampleSizeForASetOfPairings("D:/Myron/Indian_Diabetes/SequencingPilot/population.input");
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void rangeOfMaf(double prevalence, double relativeRiskIncrement, int numCases,
                                int numControls, int numTests,
                                boolean unselected) throws Exception {
    double alpha, rr;

    alpha = 0.05 / numTests;
    System.out.println("Prevalence = " + ext.formDeci(prevalence * 100, 2) + "%");
    System.out.println("n cases = " + numCases);
    System.out.println("n controls = " + numControls);
    System.out.println("n tests = " + numTests + " (alpha=" + ext.prettyP(alpha) + ")");
    System.out.println();
    System.out.println("MinorAlleleFreq\tRelativeRisk @80% power");
    for (double maf : MAFs) {
      // System.out.println(MAFs[mafIndex]+"\t"+getSampleSize(prevalence, 1.6, MAFs[mafIndex],
      // numCases, numControls, alpha, false));
      rr = getRelativeRiskAtEightyPercentPower(prevalence, maf, numCases, numControls, alpha,
                                               unselected);
      System.out.println(maf + "\t" + (rr == -9 ? "failed" : ext.formDeci(rr, 2)));
    }
  }

  public static void rangeOfRelativeRisk(double prevalence, int numTests,
                                         boolean unselected) throws Exception {
    double alpha;

    alpha = 0.05 / numTests;
    System.out.println("Prevalence = " + ext.formDeci(prevalence * 100, 2) + "%");
    System.out.println("cells = n cases = n controls (for total sample size, multiply by 2)");
    System.out.println("n tests = " + numTests + " (alpha=" + ext.prettyP(alpha) + ")");
    System.out.println();
    System.out.println("\tRelativeRisk");
    System.out.print("MinorAlleleFreq");
    for (double element : RELATIVE_RISKS) {
      System.out.print("\t" + ext.formDeci(element, 2, true));
    }
    System.out.println();
    for (double maf : MAFs) {
      System.out.print(maf);
      for (double element : RELATIVE_RISKS) {
        System.out.print("\t"
                         + getSampleSize(prevalence, element, maf, 1, 1, alpha, unselected, false));
      }
      System.out.println();
    }
  }
}
