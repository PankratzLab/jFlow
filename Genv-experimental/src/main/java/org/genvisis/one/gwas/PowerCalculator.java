package org.genvisis.one.gwas;

import java.util.Hashtable;
import java.util.Map;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Internat;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.ProbDist;
import org.pankratzlab.common.stats.Stats;

public class PowerCalculator {

  // public static final double[] MAFs = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
  // 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
  // common and rare
  public static final double[] MAFs = {
                                       // 0.001, 0.005, 0.01, 0.02, 0.03, 0.04,
                                       0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
  // just common
  // public static final double[] MAFs = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
  // 0.50};
  public static final double[] RELATIVE_RISKS = {1.10, 1.20, 1.30, 1.40, 1.60, 1.80, 2.00, 2.2, 2.4,
                                                 2.6};
  // public static final double[] MAFs = {0.20};
  public static final String[] FORMATTING_TO_REMOVE = {"<em><font color=\"navy\">", "</font></em>"};
  public static final double RR_INCREMENT = 0.01;

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

    Map<String, String> data = new Hashtable<>();
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
      // results = Internat.doSubmit("http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi", data,
      // 1000);
      results = Internat.doSubmit("http://zzz.bwh.harvard.edu/cgi-bin/cc2k.cgi", data, 1000);

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
    array = ArrayUtils.intArray(10000, -1);
    while (!found) {
      rr = 1 + index * RR_INCREMENT;
      if (array[index] == -1) {
        array[index] = getSampleSize(prevalence, rr, maf, numCases, numControls, alpha, unselected,
                                     false);
        if (array[index] == -9) {
          array[index] = getSampleSize(prevalence, rr, maf, numCases, numControls, alpha,
                                       unselected, true);
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

  // simulateSomaticDistribution(210, 0.33, 20000, 9);
  private static double simulateSomaticDistribution(int numCases, double targetPercentage,
                                                    int numGeneTests,
                                                    int numExpectedSomaticPerExome, int numReps,
                                                    boolean verbose) {

    double emp2 = 1;
    int target = (int) (numCases * targetPercentage);
    int numExceedingTarget = 0;
    for (int i = 0; i < numReps; i++) {
      int[] geneCounts = new int[numGeneTests];
      boolean exceeds = false;
      for (int j = 0; j < numCases; j++) {
        int numActualSomatic = numExpectedSomaticPerExome
                               + (int) ((Math.random() < 0.50 ? -1 : 1)
                                        * ProbDist.NormDistReverse(Math.random()) * 2);
        // System.out.println(numActualSomatic);
        for (int k = 0; k < numActualSomatic; k++) {
          geneCounts[(int) (Math.random() * numGeneTests)]++;
        }
      }
      for (int j = 0; j < geneCounts.length; j++) {
        if (geneCounts[j] >= target) {
          exceeds = true;
        }
      }
      if (exceeds) {
        numExceedingTarget++;
      }
    }

    emp2 = ((double) numExceedingTarget + 1) / ((double) numReps + 1);
    if (verbose) {
      System.out.println(numExceedingTarget + " reps exceeded target case count of " + target
                         + " (p=" + ext.prettyP(emp2) + ") which corresponds to "
                         + ext.formDeci(targetPercentage * 100, 1) + "% of cases");
    }

    return emp2;
  }

  private static void powerSomaticDistribution(int numCases, double targetPercentage,
                                               int numGeneTests, int numExpectedSomaticPerExome,
                                               int numReps, int numPowerReps) {

    int numExceeding = 0;
    for (int i = 0; i < numPowerReps; i++) {
      System.out.print(".");
      double actualNumberOfCarriersInRep = 0;
      for (int j = 0; j < numCases; j++) {
        if (Math.random() < targetPercentage) {
          actualNumberOfCarriersInRep++;
        }
      }
      if (simulateSomaticDistribution(numCases,
                                      ((double) actualNumberOfCarriersInRep / (double) numCases),
                                      numGeneTests, numExpectedSomaticPerExome, numReps,
                                      false) < 0.05) {
        numExceeding++;
      }
    }

    int target = (int) (numCases * targetPercentage);
    double power = ((double) numExceeding + 1) / ((double) numPowerReps + 1);
    System.out.println("\nWe have " + ext.formDeci(power * 100, 1)
                       + "% power to detect a target case count of " + target
                       + " or greater which corresponds to "
                       + ext.formDeci(targetPercentage * 100, 1) + "% of cases");

  }

  // simulateSomaticDistribution(210, 0.33, 20000, 9);
  private static double simulateSomaticCaseControl(int numCases, int numControls,
                                                   int numCaseCarriers, int numControlCarriers,
                                                   int numGeneTests, int numExpectedSomaticPerExome,
                                                   int numReps, boolean verbose) {

    double chiToBeat = Stats.PearsonChiSquare(numCaseCarriers, numCases - numCaseCarriers,
                                              numControlCarriers, numControls - numControlCarriers);

    double emp2 = 1;
    int numExceedingTarget = 0;
    int numMeetingCriterion = 0;
    for (int i = 0; i < numReps; i++) {
      int[][] geneCounts = new int[numGeneTests][2];
      boolean exceeds = false;
      for (int j = 0; j < numCases + numControls; j++) {
        int numActualSomatic = numExpectedSomaticPerExome
                               + (int) ((Math.random() < 0.50 ? -1 : 1)
                                        * ProbDist.NormDistReverse(Math.random()) * 2);
        // System.out.println(numActualSomatic);
        for (int k = 0; k < numActualSomatic; k++) {
          geneCounts[(int) (Math.random() * numGeneTests)][j < numCases ? 0 : 1]++;
        }
      }
      for (int j = 0; j < geneCounts.length; j++) {
        if ((double) geneCounts[j][0] / (double) numCases >= (double) geneCounts[j][0]
                                                             / (double) numControls
            && Stats.PearsonChiSquare(geneCounts[j][0], numCases - geneCounts[j][0],
                                      geneCounts[j][1],
                                      numControls - geneCounts[j][0]) >= chiToBeat) {
          exceeds = true;
        }
        // if (geneCounts[j][0] + geneCounts[j][1] >= 4) {
        // numMeetingCriterion++;
        // }
      }
      if (exceeds) {
        numExceedingTarget++;
      }
      // System.out.println(numMeetingCriterion);
    }

    emp2 = ((double) numExceedingTarget + 1) / ((double) numReps + 1);
    if (verbose) {
      System.out.println(numExceedingTarget + " of " + numReps
                         + "  reps exceeded target chi-square value of " + chiToBeat + " (p="
                         + ext.prettyP(emp2) + ") which corresponds to " + " (an odds ratio of "
                         + ext.formDeci(((double) numCaseCarriers / (double) numCases)
                                        / ((double) numControlCarriers / (double) numControls), 2)
                         + ")");
    }

    return emp2;
  }

  private static void powerSomaticCaseControl(int numCases, int numControls,
                                              double targetPercentageInCases,
                                              double targetPercentageInControls, int numGeneTests,
                                              int numExpectedSomaticPerExome, int numReps,
                                              int numPowerReps) {

    int numExceeding = 0;
    for (int i = 0; i < numPowerReps; i++) {
      System.out.print(".");

      int actualNumberOfCarrierCasesInRep = 0;
      for (int j = 0; j < numCases; j++) {
        if (Math.random() < targetPercentageInCases) {
          actualNumberOfCarrierCasesInRep++;
        }
      }

      int actualNumberOfCarrierControlsInRep = 0;
      for (int j = 0; j < numControls; j++) {
        if (Math.random() < targetPercentageInControls) {
          actualNumberOfCarrierControlsInRep++;
        }
      }

      if (simulateSomaticCaseControl(numCases, numControls, actualNumberOfCarrierCasesInRep,
                                     actualNumberOfCarrierControlsInRep, numGeneTests,
                                     numExpectedSomaticPerExome, numReps, false) < 0.05) {
        numExceeding++;
      }
    }

    double or = targetPercentageInCases / targetPercentageInControls;
    double power = ((double) numExceeding + 1) / ((double) numPowerReps + 1);
    double overallFreq = (targetPercentageInCases * (double) numCases
                          + targetPercentageInControls * (double) numControls)
                         / (double) (numCases + numControls);
    System.out.println();
    System.out.println("We have " + ext.formDeci(power * 100, 1)
                       + "% power to detect an odds ratio of " + ext.formDeci(or, 2, true)
                       + " or greater when the overall frequency of GCT cases with the variant is ("
                       + ext.prettyP(overallFreq) + ")");
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
      // rangeOfMaf(0.001, 0.01, 465, 1119, 1000000, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 434, 1119, 1000000, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 465 + 434, 1119, 1000000, false); // Poynter's MDS-AML grant, MDS
      // rangeOfMaf(0.001, 0.01, 875 + 700, 875 + 700, 1000000, false); // Poynter's MDS-AML grant,
      // MDS

      // rangeOfMaf(0.001, 0.01, 101, 406, 300000, false); // Poynter's renewal on ototoxicity
      // rangeOfMaf(0.001, 0.01, 38, 152, 500000, false); // Poynter's renewal on ototoxicity

      // rangeOfMaf(0.001, 0.01, 867 + 700 + 800, 867 + 700 + 800, 1000000, false); // Poynter's GCT
      // R21
      // rangeOfMaf(0.001, 0.01, 867 + 700 + 800, 867 + 700 + 800, 50000, false); // Poynter's GCT
      // R21
      // rangeOfMaf(0.001, 0.01, 2386, 2386 + 1000, 1000000, false); // Poynter's GCT resubmission
      // rangeOfMaf(0.001, 0.01, 479, 479, 1, false); // Poynter's GCT replication
      // rangeOfMaf(0.001, 0.01, 1000, 1000, 1000000, false); // Poynter's GCT replication

      // rangeOfMaf(0.001, 0.01, 163, 656, 300000, false); // Poynter's renewal on ototoxicity -
      // discovery
      // rangeOfMaf(0.001, 0.01, 271, 855, 300000, false); // Poynter's renewal on ototoxicity -
      // replication
      // rangeOfMaf(0.001, 0.01, 163 + 271, 656 + 855, 300000, false); // Poynter's renewal on
      // ototoxicity - combined

      // rangeOfMaf(0.001, 0.01, 163, 656, 20, false); // Poynter's renewal on ototoxicity -
      // discovery
      // rangeOfMaf(0.001, 0.01, 271, 855, 20, false); // Poynter's renewal on ototoxicity -
      // replication
      // rangeOfMaf(0.001, 0.01, 163 + 271, 656 + 855, 20, false); // Poynter's renewal on
      // ototoxicity - combined

      // rangeOfMaf(0.001, 0.01, 163 + 271, 656 + 855, 16000, false); // Poynter's renewal on
      // ototoxicity - combined, gene-based burden and PrediXcan

      rangeOfMaf(0.001, 0.01, 163 - 30, 656 - 122, 300000, false); // Poynter's renewal on
                                                                   // ototoxicity - discovery // no
                                                                   // IGCT
      rangeOfMaf(0.001, 0.01, 271, 855, 300000, false); // Poynter's renewal on ototoxicity -
                                                        // replication // no IGCT
      rangeOfMaf(0.001, 0.01, 163 + 271 - 30, 656 + 855 - 122, 300000, false); // Poynter's renewal
                                                                               // on ototoxicity -
                                                                               // combined // no
                                                                               // IGCT

      rangeOfMaf(0.001, 0.01, 163 - 30, 656 - 122, 20, false); // Poynter's renewal on ototoxicity -
                                                               // discovery // no IGCT
      rangeOfMaf(0.001, 0.01, 271, 855, 20, false); // Poynter's renewal on ototoxicity -
                                                    // replication // no IGCT
      rangeOfMaf(0.001, 0.01, 163 + 271 - 30, 656 + 855 - 122, 20, false); // Poynter's renewal on
                                                                           // ototoxicity - combined
                                                                           // // no IGCT

      rangeOfMaf(0.001, 0.01, 163 + 271 - 30, 656 + 855 - 122, 16000, false); // Poynter's renewal
                                                                              // on ototoxicity -
                                                                              // combined,
                                                                              // gene-based burden
                                                                              // and PrediXcan // no
                                                                              // IGCT

      // rangeOfMaf(0.001, 0.01, 665 , 1119 , 100000, false); // Poynter's MDS-AML resubmission,
      // discovery
      // rangeOfMaf(0.001, 0.01, 665 , 1119 , 50000, false); // Poynter's MDS-AML resubmission, meta
      // rangeOfMaf(0.001, 0.01, 1700 , 4597, 100000, false); // Poynter's MDS-AML resubmission,
      // replication
      // rangeOfMaf(0.001, 0.01, 1700 , 4597, 50000, false); // Poynter's MDS-AML resubmission,
      // replication
      // rangeOfMaf(0.001, 0.01, 665 + 1700 , 1119 + 4597, 100000, false); // Poynter's MDS-AML
      // resubmission, meta

      // rangeOfMaf(0.001, 0.01, 192, 406, 1000000, false); // Logan's Ewings grant
      // rangeOfMaf(0.001, 0.01, 598, 598, 1000000, false); // Logan's Ewings grant
      // rangeOfMaf(0.001, 0.01, 338, 717, 1000000, false); // Logan's Ewings grant

      // rangeOfMaf(0.001, 0.01, 148, 316, 1000000, false); // Logan's Ewings grant, ALSF, mets
      // rangeOfMaf(0.001, 0.01, 464, 464, 1000000, false); // Logan's Ewings grant, ALSF, gwas

      // rangeOfMaf(0.001, 0.01, 181, 385, 1000000, false); // Logan's Ewings grant, ALSF v2, mets
      // rangeOfMaf(0.001, 0.01, 566, 464, 1000000, false); // Logan's Ewings grant, ALSF v2, gwas

      // rangeOfRelativeRisk(0.15, 200, false);
      // getSampleSize();
      // getSampleSizeForASetOfPairings("D:/Myron/Indian_Diabetes/SequencingPilot/power.input");
      // getSampleSizeForASetOfPairings("D:/Myron/Indian_Diabetes/SequencingPilot/population.input");

      // simulateSomaticDistribution(210, 0.025, 20000, 9, 1000, true);
      // powerSomaticDistribution(210, 0.03, 20000, 9, 100, 100);
      // powerSomaticDistribution(285, 0.03, 500, 1, 100, 100);
      // simulateSomaticCaseControl(80, 130, 16, 12, 20000, 9, 100, true);
      // powerSomaticCaseControl(80, 130, 0.18, 0.05, 20000, 9, 100, 100);

      // rangeOfMaf(0.001, 0.01, 2305, 7251, 1000000, false); // Lindsay Males
      // rangeOfMaf(0.001, 0.01, 1889, 7251, 1000000, false); // Lindsay Females
      // rangeOfMaf(0.001, 0.01, 4191, 14502, 1000000, false); // Lindsay Combined
      // rangeOfMaf(0.001, 0.01, 846, 651, 5000, false); // Lindsay B-ALL CNV
      // rangeOfMaf(0.001, 0.01, 197, 61, 5000, false); // Lindsay T-ALL CNV

      // rangeOfMaf(0.001, 0.01, 1043, 712, 5000, false); // Lindsay Combined CNV

      // rangeOfMaf(0.001, 0.01, 1368, 7251, 1000000, false); // Lindsay Males take 2
      // rangeOfMaf(0.001, 0.01, 1119, 7251, 1000000, false); // Lindsay Females take 2

    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
