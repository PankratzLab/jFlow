package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import org.genvisis.cnv.filesys.CNVQC;
import org.genvisis.cnv.filesys.MarkerFreqs;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.stats.Maths;

public class QCIterator implements Runnable {
  private static boolean alreadyDefinedID(
      Hashtable<String, Hashtable<String, Integer>> defineCompHash, String[] line, Logger log) {
    boolean alreadyDefined = false;
    for (String element : line) {
      // same id not allowed in same line
      if (defineCompHash.containsKey(element)) {
        log.reportError("Warning - duplicate IDs were detectected in the comparision file: "
            + element + " was seen twice, only the first comparison will be used ");
        alreadyDefined = true;
      }
    }
    return alreadyDefined;
  }

  private static double[] binIt(double startVal, double stopVal, int numBins) {
    double[] values;
    if (startVal == stopVal) {
      values = new double[1];
      values[0] = startVal;
    } else {
      values = new double[numBins + 1];
      double inc = getIncrement(startVal, stopVal, numBins);
      for (int i = 0; i < numBins + 1; i++) {
        values[i] = (inc * i) + startVal;
      }
    }
    return values;
  }

  private static void checkThreadStatus(int processors, Thread[] threads) {
    boolean complete;
    complete = false;
    while (!complete) {
      complete = true;
      for (int i = 0; i < processors; i++) {
        if (threads[i].isAlive()) {
          complete = false;
        }
      }
      if (!complete) {
        try {
          Thread.sleep(1000L);
        } catch (InterruptedException ex) {
        }
      }
    }
  }

  // This is the main qc tracker, if there are more calls passing the thresholds for a given percent
  // target, the qc parameters will be updated
  // copy number specific
  private static OptimizedQCThresholds[] checkThresholds(CNVComparison cnvComp,
      OptimizedQCThresholds[] bestOptqcs, double targetConcordancePercentage,
      OptimizedQCThresholds qcIteration) {
    double[] averageCNPercent = cnvComp.getAverageCNPercent();
    int[] numCallsPassing = cnvComp.getFilteredCallsAvailable();
    for (int i = 0; i < averageCNPercent.length; i++) {
      if (averageCNPercent[i] >= targetConcordancePercentage
          && bestOptqcs[i].getCallsPassingFilter() < numCallsPassing[i]) {
        bestOptqcs[i] = new OptimizedQCThresholds(qcIteration, targetConcordancePercentage,
            averageCNPercent[i], numCallsPassing[i], cnvComp.getGoodCalls()[i],
            cnvComp.getTotalCallsAvailable()[i], cnvComp.getNumberIndsCompared(), i);
      }
    }
    return bestOptqcs;
  }

  // collects results from all threads
  private static OptimizedQCThresholds[][] collectOptqcs(int processors, QCIterator[] qcIts,
      double[] allPercents, Logger log) {
    OptimizedQCThresholds[][] allOptQcs = new OptimizedQCThresholds[allPercents.length][];
    log.report(ext.getTime() + " Collecting QC thresholds for " + allPercents.length
        + " target percents from available threads...");
    int indIndex = 0;
    int counter = 0;
    for (int i = 0; i < allPercents.length; i++) {
      counter++;
      if (counter > processors) {
        indIndex += 1;
        counter = 1;
      }

      if (qcIts[i % processors].getTargetPercentages()[indIndex] != allPercents[i]) {
        log.reportError("Error - recieved unmatched results while collecting results for "
            + Array.toStr(qcIts[i % processors].getOptqcs()[indIndex], ",") + "\t"
            + allPercents[i]);
        System.exit(1);
      } else if (qcIts[i % processors].getTargetPercentages()[indIndex] == allPercents[i]) {
        allOptQcs[i] = qcIts[i % processors].getOptqcs()[indIndex];
      }
    }
    log.report(ext.getTime() + " Sucessfully collected QC thresholds for " + allPercents.length
        + " target percents from available threads...");
    return allOptQcs;
  }

  // TODO to compare two different plink format files;
  public static void compareFiles(Project proj, String[] plinkCnvQCs, String[] SampleQCFiles,
      String compareThresholdFileName, Logger log) {
    Hashtable<Integer, Hashtable<String, CNVariantQC[]>> fileIndcnVariantQCs =
        new Hashtable<Integer, Hashtable<String, CNVariantQC[]>>();
    String[][] inds = new String[plinkCnvQCs.length][];
    int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(plinkCnvQCs.length, 2);
    // CNVComparison cnvComparison;
    for (int i = 0; i < plinkCnvQCs.length; i++) {
      CNVariantQC[] cnVariantQCs = CNVariantQC.getCNVariantQCFromPlinkFile(proj, plinkCnvQCs[i]);
      inds[i] = CNVariantQC.getIDList(cnVariantQCs, null);
      fileIndcnVariantQCs.put(i, CNVariantQC.getIndCNVQCs(inds[i], cnVariantQCs));
    }
    for (int[] allPossibleCombination : allPossibleCombinations) {
      ArrayList<CNVariantQC[]> cnvQCX = new ArrayList<CNVariantQC[]>();
      ArrayList<CNVariantQC[]> cnvQCY = new ArrayList<CNVariantQC[]>();
      int fileX = allPossibleCombination[0];
      int fileY = allPossibleCombination[1];
      String[] indsx = inds[fileX];
      String[] indsy = inds[fileY];
      for (String element : indsx) {
        for (String element2 : indsy) {
          if (fileIndcnVariantQCs.get(fileX).containsKey(element)
              && fileIndcnVariantQCs.get(fileY).containsKey(element2) && element.equals(element2)) {
            cnvQCX.add(fileIndcnVariantQCs.get(fileX).get(element));
            cnvQCY.add(fileIndcnVariantQCs.get(fileY).get(element2));
            // log.report(fileX +"\t"+indsx[j] +"\t"+fileY+"\t"+indsy[k]);
          }
        }
      }
      CNVariantQC[][] unfilteredcnvsQCs1 = cnvQCX.toArray(new CNVariantQC[cnvQCX.size()][]);
      CNVariantQC[][] unfilteredcnvsQCs2 = cnvQCY.toArray(new CNVariantQC[cnvQCY.size()][]);
      log.report("Comparing " + unfilteredcnvsQCs1.length + " individuals with "
          + unfilteredcnvsQCs2.length + " individuals");
      // public CNVComparison(CNVariantQC[][] unfilteredcnvsQCs1, CNVariantQC[][]
      // unfilteredcnvsQCs2, Hashtable<String, CNVSampleQC> cnvSampleQCHash, OptimizedQCThresholds
      // qcThresholds, int filterType, Logger log) {
      CNVComparison cnvComparison;
      if (SampleQCFiles == null || compareThresholdFileName == null) {
        cnvComparison =
            new CNVComparison(unfilteredcnvsQCs1, unfilteredcnvsQCs2, null, null, 0, log);
      } else {
        log.report("Filtering on sample QC before comparisons");
        OptimizedQCThresholds qcThresholds = OptimizedQCThresholds.loadThresholdsFromTxt(
            proj.PROJECT_DIRECTORY.getValue() + compareThresholdFileName, log);
        Hashtable<String, CNVSampleQC> cnvSampleQCHash1 =
            CNVSampleQC.getSampleQCs(proj, SampleQCFiles[allPossibleCombination[0]]);
        Hashtable<String, CNVSampleQC> cnvSampleQCHash2 =
            CNVSampleQC.getSampleQCs(proj, SampleQCFiles[allPossibleCombination[1]]);
        CNVComparison filteredOnly1 =
            new CNVComparison(unfilteredcnvsQCs1, cnvSampleQCHash1, qcThresholds, 1, log);
        CNVComparison filteredOnly2 =
            new CNVComparison(unfilteredcnvsQCs2, cnvSampleQCHash2, qcThresholds, 1, log);
        CNVariantQC[][] filteredcnvsQCs1 = filteredOnly1.getFilteredcnvQCs1();
        CNVariantQC[][] filteredcnvsQCs2 = filteredOnly2.getFilteredcnvQCs1();
        ArrayList<CNVariantQC[]> cnvQC1 = new ArrayList<CNVariantQC[]>();
        ArrayList<CNVariantQC[]> cnvQC2 = new ArrayList<CNVariantQC[]>();
        for (int j = 0; j < filteredcnvsQCs1.length; j++) {
          if (filteredcnvsQCs1[j] != null && filteredcnvsQCs2[j] != null) {
            cnvQC1.add(filteredcnvsQCs1[j]);
            cnvQC2.add(filteredcnvsQCs2[j]);
          }
        }
        if (cnvQC1.size() != cnvQC2.size()) {
          log.reportError("Error - number of individuals for comparisions do not match");
          System.exit(1);
        }
        log.report(cnvQC1.size() + "\t" + cnvQC2.size());
        filteredcnvsQCs1 = cnvQC1.toArray(new CNVariantQC[cnvQC1.size()][]);
        filteredcnvsQCs2 = cnvQC2.toArray(new CNVariantQC[cnvQC2.size()][]);
        cnvComparison = new CNVComparison(filteredcnvsQCs1, filteredcnvsQCs2, null, null, 0, log);
      }
      ArrayList<CNVariantQC> misses = cnvComparison.getMisses();
      writeMisses(proj, misses, log);
      log.report(Array.toStr(cnvComparison.getGoodCalls()));
      log.report(Array.toStr(cnvComparison.getFilteredCallsAvailable()));
      log.report(Array.toStr(cnvComparison.getAverageCNPercent()));
    }
  }

  public static void convertToQCFormat(Project proj, String plinkCnvs, String markerMAFser,
      String output, String QCsubset, int threads) {
    String[] inds;
    MarkerSet markerSet = proj.getMarkerSet();
    int[][] indices = markerSet.getIndicesByChr();
    int[] positions = markerSet.getPositions();
    String[] markerNames = markerSet.getMarkerNames();
    MarkerFreqs markerFreqs;
    Logger log;
    double[] mafs;

    log = proj.getLog();
    if (markerMAFser != null) {
      markerFreqs = MarkerFreqs.load(proj.PROJECT_DIRECTORY.getValue() + markerMAFser, false);
      mafs = markerFreqs.getMafs();
    } else {
      markerFreqs = null;
      mafs = new double[markerNames.length];
      Arrays.fill(mafs, 0.5);
      log.report(
          "Warning - a marker Frequency file was not provided, setting all MAF values to 0.5");
    }
    if (markerFreqs != null && markerFreqs.getFingerprint() != markerSet.getFingerprint()) {
      log.reportError(
          "Error - mismatched marker fingerprints in the project's marker set and the imported AlleleFrequency file ("
              + markerMAFser + "); aborting");
      System.exit(1);
    }
    log.report(
        ext.getTime() + " Retrieving cnvs from " + proj.PROJECT_DIRECTORY.getValue() + plinkCnvs);
    CNVariantQC[] cnVariantQCs = CNVariantQC.getCNVariantQCFromPlinkFile(proj, plinkCnvs);
    log.report(ext.getTime() + " Finished retrieving cnvs from " + proj.PROJECT_DIRECTORY.getValue()
        + plinkCnvs);
    if (QCsubset != null) {
      log.report(ext.getTime() + " Filtering cnvs by the QC subset file");
      inds = CNVariantQC.getIDList(cnVariantQCs,
          defineCompLists(proj.PROJECT_DIRECTORY.getValue(), QCsubset, log));
      log.report(ext.getTime()
          + " Finished filtering cnvs by the QC subset file , computing cnv QC metrics for "
          + inds.length + " individuals");
    } else {
      inds = CNVariantQC.getIDList(cnVariantQCs, null);
      log.report(
          ext.getTime() + " Using " + inds.length + " individuals to compute cnv QC metrics ");
    }
    Hashtable<String, CNVariantQC[]> allIndcnVariantQCs =
        CNVariantQC.getIndCNVQCs(inds, cnVariantQCs);
    Hashtable<String, Double> markerMAFhash = hashMAFs(markerNames, mafs);
    log.report(ext.getTime() + " Retrieving markers in cnvs and assigning MAFs");
    for (String ind : inds) {
      for (int j = 0; j < allIndcnVariantQCs.get(ind).length; j++) {
        allIndcnVariantQCs.get(ind)[j].findMarkerNamesinCNV(proj, indices, positions, markerNames);
        allIndcnVariantQCs.get(ind)[j].assignMAFs(markerMAFhash, log);
      }
    }
    CNVariantQC[][] allIndcnVariantQCsArrays = CNValidate.computeMultiThreadedValidations(proj,
        inds, allIndcnVariantQCs, markerSet, threads);
    printNewCNVariantQCFile(proj, output, inds, allIndcnVariantQCsArrays);
    log.report(ext.getTime() + " Completed QC computations");
  }

  private static void defineComparisons(
      Hashtable<String, Hashtable<String, Integer>> defineCompHash, String[] line, Logger log) {
    Hashtable<String, Integer> defined = new Hashtable<String, Integer>();
    for (int i = 0; i < line.length; i++) {
      for (int j = 0; j < line.length; j++) {
        if (i != j && line[i] == line[j]) {
          log.reportError("Warning - duplicate IDs were detectected in the comparision file: "
              + line[i] + " in column " + i + "and " + line[j] + " in column " + j);
        }
        defined.put(line[j], j);
        defineCompHash.put(line[i], defined);
      }
    }
  }

  // Defines the comparisions of interest
  private static Hashtable<String, Hashtable<String, Integer>> defineCompLists(String rootDir,
      String compFile, Logger log) {
    Hashtable<String, Hashtable<String, Integer>> defineCompHash =
        new Hashtable<String, Hashtable<String, Integer>>();
    BufferedReader reader;
    String[] line;
    int maxNumComparisions = 0;
    try {
      reader = new BufferedReader(new FileReader(rootDir + compFile));
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t");
        if (line.length > maxNumComparisions) {
          maxNumComparisions = line.length;
        }
        // same id in multiple lines not allowed
        if (alreadyDefinedID(defineCompHash, line, log)) {
          continue;
        }
        defineComparisons(defineCompHash, line, log);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + rootDir + compFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + rootDir + compFile + "\"");
      System.exit(2);
    }
    return defineCompHash;
  }

  public static void filterByComparison(Project proj, String plinkCnvQCs, String duplicatesFile,
      Logger log) {
    Hashtable<String, Hashtable<String, Integer>> defineCompHash =
        defineCompLists(proj.PROJECT_DIRECTORY.getValue(), duplicatesFile, log);
    CNVariantQC.filterCNVQCsByComparison(proj, plinkCnvQCs, defineCompHash);

  }

  public static void filterCNVsByQCThresholds(Project proj, String plinkCnvQCs, String SampleQCFile,
      String qcThresholdFileName, String output, String QCsubset, int optimizationType) {
    Logger log = proj.getLog();
    OptimizedQCThresholds qcThresholds = OptimizedQCThresholds
        .loadThresholdsFromTxt(proj.PROJECT_DIRECTORY.getValue() + qcThresholdFileName, log);
    log.report(ext.getTime() + " Loaded thresholds from " + proj.PROJECT_DIRECTORY.getValue()
        + qcThresholdFileName);
    log.report(ext.getTime() + " Retrieving sample qc data from "
        + proj.PROJECT_DIRECTORY.getValue() + SampleQCFile);
    Hashtable<String, CNVSampleQC> cnvSampleQCHash = CNVSampleQC.getSampleQCs(proj, SampleQCFile);
    log.report(
        ext.getTime() + " Loading cnvQCs from " + proj.PROJECT_DIRECTORY.getValue() + plinkCnvQCs);
    CNVariantQC[] cnVariantQCs =
        CNVQC.load(proj.PROJECT_DIRECTORY.getValue() + plinkCnvQCs, false).getCnVariantQCs();
    log.report(ext.getTime() + " Finished loading cnvQCs from " + proj.PROJECT_DIRECTORY.getValue()
        + plinkCnvQCs);
    String[] inds = CNVariantQC.getIDList(cnVariantQCs, null);
    CNVariantQC[][] unfilteredcnvsQCs = new CNVariantQC[inds.length][];
    Hashtable<String, CNVariantQC[]> allIndcnVariantQCs =
        CNVariantQC.getIndCNVQCs(inds, cnVariantQCs);
    for (int i = 0; i < inds.length; i++) {
      unfilteredcnvsQCs[i] = allIndcnVariantQCs.get(inds[i]);
    }
    log.report(ext.getTime() + " Beginning to filter " + plinkCnvQCs + " using filter type "
        + CNVComparison.QC_PARAMETERs[optimizationType]);
    CNVComparison filteredOnly =
        new CNVComparison(unfilteredcnvsQCs, cnvSampleQCHash, qcThresholds, optimizationType, log);
    log.report(ext.getTime() + " Finished Filtering " + plinkCnvQCs + " using filter type "
        + CNVComparison.QC_PARAMETERs[optimizationType]);
    CNVariantQC[][] filteredcnvsQCs = filteredOnly.getFilteredcnvQCs1();
    log.report(ext.getTime() + " Creating output files");
    summarizeFiltering(proj, inds, filteredcnvsQCs, cnvSampleQCHash, output, QCsubset);
    log.report(ext.getTime() + " Finished");

  }



  private static ArrayList<ArrayList<Double>> getcabinet(double[] percents, int processors) {
    ArrayList<ArrayList<Double>> cabinet = new ArrayList<ArrayList<Double>>();

    for (int i = 0; i < processors; i++) {
      cabinet.add(new ArrayList<Double>());
    }
    for (int i = 0; i < percents.length; i++) {
      cabinet.get(i % processors).add(percents[i]);
    }
    return cabinet;
  }

  private static double[] getDoubleThreadPercents(ArrayList<Double> doubles) {
    double[] threadPercents = new double[doubles.size()];
    for (int i = 0; i < doubles.size(); i++) {
      threadPercents[i] = doubles.get(i);
    }
    return threadPercents;
  }

  private static double getIncrement(double startVal, double stopVal, int numBins) {
    return (stopVal - startVal) / numBins;
  }

  private static OptimizedQCThresholds[] getQCIterations() {
    // TODO these will be input parameters, with defaults if not provided
    // CNV -specific
    int numMarkerStart = 20;
    int numMarkerStop = 20;
    double[] confCutoffs = binIt(0, 100, 300);
    double[] alphas = binIt(0, 100, 300);
    // less than
    double[] BAFQCcutoffs = binIt(200000, 200000, 1000);
    // greater than
    double[] twopqCutoffs = binIt(0, 0, 10);
    // less than
    double[] hetCutoffs = binIt(1, 1, 10);
    // less than
    double[] cnvLRRSTDevs = binIt(100, 100, 100);
    double[] pennconf = binIt(0, 0, 10);
    double[] bafDrifts = binIt(0.03, .03, 10);
    double[] kbSize = binIt(0, 0, 10);
    double[] kbDensity = binIt(0, 0, 10);
    // SampleSpecific
    int numSampleCNVsStart = 100;
    int numSampleCNVsStop = 100;
    double[] lrrCutoffs = binIt(0.35, .35, 20);
    double[] GCWFCutoffs = binIt(0.02, 0.02, 20);
    double[] sampleCallRates = binIt(0.96, .96, 20);

    ArrayList<OptimizedQCThresholds> qcThresholds = new ArrayList<OptimizedQCThresholds>();
    // qcThresholds.add(new OptimizedQCThresholds(alphas[i], confCutoffs[j], lrrCutoffs[k], l,
    // BAFQCcutoffs[m], twopqCutoffs[n], hetCutoffs[o], GCWFCutoffs[p], q));
    // OptimizedQCThresholds noFilter = new OptimizedQCThresholds(0.5, 1.5, 20.0, 200000, 1, 0.0,
    // 1.0, 10.0, 80000, 10.0);
    // qcThresholds.add(noFilter);
    // public OptimizedQCThresholds(double alpha, double confCutoff, double lrrCutoff, int
    // numMarkers, double BAFQCcutoff, double twopqCutoff, double hetCutoff, double GCWF, int
    // numSampleCNVs, double cnvLRRSTDev) {
    for (double alpha : alphas) {
      for (double confCutoff : confCutoffs) {
        for (double lrrCutoff : lrrCutoffs) {
          for (int l = numMarkerStart; l <= numMarkerStop; l++) {
            for (double bafqCcutoff : BAFQCcutoffs) {
              for (double twopqCutoff : twopqCutoffs) {
                for (double hetCutoff : hetCutoffs) {
                  for (double gcwfCutoff : GCWFCutoffs) {
                    for (int q = numSampleCNVsStart; q <= numSampleCNVsStop; q++) {
                      for (double cnvLRRSTDev : cnvLRRSTDevs) {
                        for (double sampleCallRate : sampleCallRates) {
                          for (double element : pennconf) {
                            for (double bafDrift : bafDrifts) {
                              for (double element2 : kbSize) {
                                for (double element3 : kbDensity) {
                                  qcThresholds.add(new OptimizedQCThresholds(alpha, confCutoff,
                                      lrrCutoff, l, bafqCcutoff, twopqCutoff, hetCutoff, gcwfCutoff,
                                      q, cnvLRRSTDev, sampleCallRate, element, bafDrift, element2,
                                      element3));
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return qcThresholds.toArray(new OptimizedQCThresholds[qcThresholds.size()]);
  }

  private static Hashtable<String, Double> hashMAFs(String[] markerNames, double[] mafs) {
    Hashtable<String, Double> markerMAFs = new Hashtable<String, Double>();
    for (int i = 0; i < markerNames.length; i++) {
      markerMAFs.put(markerNames[i], mafs[i]);
    }
    return markerMAFs;
  }

  private static QCIterator[] iteratePercents(int processors, Thread[] threads,
      ArrayList<ArrayList<Double>> cabinet, CNVariantQC[][][] cnvQCsAssigned,
      Hashtable<String, CNVSampleQC> cnvSampleQCHash, int optimizationType, Logger log) {
    QCIterator[] qcIts = new QCIterator[processors];
    for (int i = 0; i < processors; i++) {
      qcIts[i] = new QCIterator(cnvQCsAssigned, cnvSampleQCHash,
          getDoubleThreadPercents(cabinet.get(i)), optimizationType, log);
      threads[i] = new Thread(qcIts[i]);
      threads[i].start();
    }
    checkThreadStatus(processors, threads);
    return qcIts;
  }

  public static void main(String[] args) {
    String filename = null;
    int numArgs = args.length;
    String usage = "TODO";
    // String usage = "\n"+
    // "cnv.qc.QCIterator requires 0-4 arguments\n"+
    // " (0) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+"
    // (default))\n"+
    // " (1) plink CNV format .cnv file to test qc thresholds (i.e. cnvs=all_cnvs.cnv\n"+
    // " (2) a tab delimited file defining the comparisons (i.e. comp=replicates.txt (default)) "+
    // " (3) name of the output file (i.e. out=qcThresholds.txt (default)"+
    // "";
    // String filename = "C:/workspace/Genvisis/projects/ARICGenvisis_CEL_11908.properties";
    // C:\workspace\Genvisis\projects\ARICGenvisis_CEL_11908.properties
    String duplicatesFile = "rootsANDdubs.comp.txt";
    String plinkCnvs = "Not_all_gw6.cnv";
    String MarkerFreqs = "MarkerFreq.ser";
    // String plinkCnvs = "all_gw6.cnv ";
    // String beastHeights = plinkCnvs.replaceAll(".cnv", "heights");
    String plinkCnvQCs = "qcThresholds.txt.GENQC_COMP.ser";
    String logfile = null;
    String QCsubset = null;
    // String output;
    String qcThresholdFileName = "qcParmaters.txt";
    String output = "qcThresholds.txt";
    String SampleQCFile = "Sample_QC.xln";
    int optimizationType = 1;
    boolean convert = false;
    boolean filter = false;
    int threads = 8;
    // String plinkCnvs = null;
    Logger log;
    Project proj;
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("cnvs=")) {
        plinkCnvs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("comp=")) {
        duplicatesFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        output = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("cnvQCs=")) {
        plinkCnvQCs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("qcsubset=")) {
        QCsubset = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("mafs=")) {
        MarkerFreqs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("params=")) {
        qcThresholdFileName = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("type=")) {
        optimizationType = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("threads=")) {
        threads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-convert")) {
        convert = true;
        numArgs--;
      } else if (arg.startsWith("-filter")) {
        filter = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      proj = new Project(filename, false);
      log = proj.getLog();
      System.out.println("Logging progress to " + log.getFilename());
      if (convert) {
        output += ".GENQC";
        convertToQCFormat(proj, plinkCnvs, MarkerFreqs, output, QCsubset, threads);
      } else if (filter) {
        filterCNVsByQCThresholds(proj, plinkCnvQCs, SampleQCFile, qcThresholdFileName, output,
            QCsubset, optimizationType);
      } else {
        // String compareThresholdFileName = "qcParms_Compare.txt";
        // String[] plinkCnvsCompare = {"all_gw6.cnv" ,"all_JL_gw6.cnv"};
        // String[] sampleQCs ={"Sample_QC.xln","Sample_QC_JL.xln"};
        // log.report("Comparing "+ Array.toStr(plinkCnvsCompare));
        // compareFiles( proj, plinkCnvsCompare,sampleQCs,compareThresholdFileName, log);
        // output += ".summary";

        optimizeQCThresholds(proj, plinkCnvQCs, duplicatesFile, SampleQCFile, output,
            optimizationType);
      }
      Files.backup(logfile, proj.PROJECT_DIRECTORY.getValue(),
          proj.PROJECT_DIRECTORY.getValue() + "LOGBACKUP/");

    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  // preps CNV file for comparision, loads CNV heights, and sends jobs out for parameter iteration
  // at each percent target
  public static void optimizeQCThresholds(Project proj, String plinkCnvQCs, String duplicatesFile,
      String SampleQCFile, String output, int optimizationType) {
    Logger log = proj.getLog();
    int processors = Runtime.getRuntime().availableProcessors();
    double[] allPercents = binIt(.40, 1, 30);
    getcabinet(allPercents, processors);
    Thread[] threads = new Thread[processors];
    log.report(ext.getTime() + " Prepping cnvQCs in " + plinkCnvQCs + " for comparisons");
    Hashtable<String, Hashtable<String, Integer>> defineCompHash =
        defineCompLists(proj.PROJECT_DIRECTORY.getValue(), duplicatesFile, log);
    CNVariantQC[][][] cnvQCsAssigned =
        CNVariantQC.prepCNVQCsForComparison(proj, plinkCnvQCs, defineCompHash);
    log.report(ext.getTime() + " Finished prepping cnvQCs in " + plinkCnvQCs + " for "
        + cnvQCsAssigned[0].length + " comparisons");
    log.report(
        ext.getTime() + " Beginning iterations for " + allPercents.length + " target percentages");
    Hashtable<String, CNVSampleQC> cnvSampleQCHash = CNVSampleQC.getSampleQCs(proj, SampleQCFile);
    QCIterator[] qcits = iteratePercents(processors, threads, getcabinet(allPercents, processors),
        cnvQCsAssigned, cnvSampleQCHash, optimizationType, log);
    log.report(
        ext.getTime() + " Finished iterations for " + allPercents.length + " target percentages");
    OptimizedQCThresholds[][] optqcs = collectOptqcs(processors, qcits, allPercents, log);
    summarizeOptqcs(proj, optqcs, output);
  }

  private static void printNewCNVariantQCFile(Project proj, String output, String[] inds,
      CNVariantQC[][] allIndcnVariantQCsArrays) {
    ArrayList<CNVariantQC> toSerialize = new ArrayList<CNVariantQC>();
    for (int i = 0; i < inds.length; i++) {
      CNVariantQC[] indcnVariantQCs = allIndcnVariantQCsArrays[i];
      for (CNVariantQC indcnVariantQC : indcnVariantQCs) {
        toSerialize.add(indcnVariantQC);
      }
    }

    // CNVariantQC[] toTest = toSerialize.toArray(new CNVariantQC[toSerialize.size()]);

    new CNVQC(toSerialize.toArray(new CNVariantQC[toSerialize.size()]))
        .serialize(proj.PROJECT_DIRECTORY.getValue() + output + ".ser");
  }

  private static void summarizeFiltering(Project proj, String[] inds,
      CNVariantQC[][] filteredcnvsQCs, Hashtable<String, CNVSampleQC> cnvSampleQCHash,
      String output, String QCsubset) {
    PrintWriter sampleWriter, cnvWriter, famWriter;
    Hashtable<String, Hashtable<String, Integer>> defineCompHash = null;
    Logger log = proj.getLog();

    if (QCsubset != null) {
      defineCompHash = defineCompLists(proj.PROJECT_DIRECTORY.getValue(), QCsubset, log);
    }
    try {
      sampleWriter =
          new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output + ".txt"));
      famWriter =
          new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output + ".fam"));
      cnvWriter =
          new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output + ".cnv"));
      sampleWriter.println("Sample\tPassQC?0:1");
      cnvWriter.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
      for (int i = 0; i < inds.length; i++) {
        if (defineCompHash.containsKey(inds[i]) || QCsubset == null) {
          if (filteredcnvsQCs[i] != null) {
            String lookup = filteredcnvsQCs[i][0].getSourceFile();
            int sex = proj.getSampleData(0, false).getSexForIndividual(lookup);
            famWriter.println(filteredcnvsQCs[i][0].getCnVariant().getFamilyID() + "\t"
                + filteredcnvsQCs[i][0].getCnVariant().getIndividualID() + "\t0\t0\t" + sex
                + "\t1");
            sampleWriter.println(inds[i] + "\t0");
            for (int j = 0; j < filteredcnvsQCs[i].length; j++) {
              cnvWriter.println(filteredcnvsQCs[i][j].getCnVariant().toPlinkFormat());
            }
          } else {
            sampleWriter.println(inds[i] + "\t1");
          }
        }
      }
      sampleWriter.close();
      famWriter.close();
      cnvWriter.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + output);
      log.reportException(e);
      System.exit(1);
    }
  }

  private static void summarizeOptqcs(Project proj, OptimizedQCThresholds[][] optqcs,
      String output) {
    Logger log = proj.getLog();

    try {
      PrintWriter qcoutput =
          new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output, false));
      qcoutput.println(Array.toStr(OptimizedQCThresholds.OPT_QC_HEADS));
      for (OptimizedQCThresholds[] optqc : optqcs) {
        qcoutput.println(Array.toStr(OptimizedQCThresholds.OPT_QC_HEADS));
        for (int k = 0; k < optqc.length; k++) {
          qcoutput.print(optqc[k].getDisplayString());
          qcoutput.print("\n");
        }
      }
      qcoutput.close();
    } catch (IOException e) {
      log.reportError("Could Not Open " + proj.PROJECT_DIRECTORY.getValue() + output);
      log.reportException(e);
    }
  }

  private static void writeMisses(Project proj, ArrayList<CNVariantQC> misses, Logger log) {
    log.report("writing unmatched calls to " + proj.PROJECT_DIRECTORY.getValue() + "misses.cnv");
    try {
      PrintWriter missWriter =
          new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + "misses.cnv"));
      missWriter.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
      for (int i = 0; i < misses.size(); i++) {
        missWriter.println(misses.get(i).getCnVariant().toPlinkFormat());
      }
      missWriter.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + "misses.cnv");
      e.printStackTrace();
      System.exit(1);
    }
  }

  private final Hashtable<String, CNVSampleQC> cnvSampleQCHash;



  private final CNVariantQC[][][] cnvQCsAssigned;

  private final double[] targetPercentages;

  private final OptimizedQCThresholds[][] optqcs;

  private final int optimizationType;

  private final Logger log;

  public QCIterator(CNVariantQC[][][] cnvQCsAssigned,
      Hashtable<String, CNVSampleQC> cnvSampleQCHash, double[] targetPercentages,
      int optimizationType, Logger log) {
    this.cnvQCsAssigned = cnvQCsAssigned;
    this.cnvSampleQCHash = cnvSampleQCHash;
    this.targetPercentages = targetPercentages;
    optqcs = new OptimizedQCThresholds[targetPercentages.length][];
    this.optimizationType = optimizationType;
    this.log = log;
  }

  public OptimizedQCThresholds[][] getOptqcs() {
    return optqcs;
  }

  public double[] getTargetPercentages() {
    return targetPercentages;
  }

  private OptimizedQCThresholds[] iterate(Hashtable<String, CNVSampleQC> cnvSampleQCHash,
      CNVariantQC[][][] cnvQCsAssigned, double targetConcordancePercentage, int optimizationType,
      Logger log) {
    OptimizedQCThresholds[] bestOptqcs =
        OptimizedQCThresholds.getNewOptqcs(targetConcordancePercentage, 6);
    OptimizedQCThresholds[] qcThresholds = getQCIterations();
    log.report(ext.getTime() + " Iterating " + qcThresholds.length + " parameter combinations");
    for (OptimizedQCThresholds qcThreshold : qcThresholds) {
      CNVComparison cnvComp = new CNVComparison(cnvQCsAssigned[0], cnvQCsAssigned[1],
          cnvSampleQCHash, qcThreshold, optimizationType, log);
      bestOptqcs = checkThresholds(cnvComp, bestOptqcs, targetConcordancePercentage, qcThreshold);
    }
    return bestOptqcs;
  }

  @Override
  public void run() {
    for (int i = 0; i < targetPercentages.length; i++) {
      log.report(
          ext.getTime() + " Beginning iterations for target percentage " + targetPercentages[i]);
      optqcs[i] =
          iterate(cnvSampleQCHash, cnvQCsAssigned, targetPercentages[i], optimizationType, log);
      log.report(
          ext.getTime() + " Finished iterations for target percentage " + targetPercentages[i]);
    }
  }
}
