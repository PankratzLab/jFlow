package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

// class to store optimized thresholds
public class OptimizedQCThresholds {
  public static final String[] QC_Thresholds =
      {"ALPHA", "BEAST_CONF", "SAMPLE_LRR_SD", "CNV_NUM_MARKERS", "BAF_QC", "TWOPQ", "PERCENT_HET",
       "GCWF", "NUM_SAMPLE_CNVS", "CNV_LRR_SD", "SAMPLE_CALL_RATE", "PENN_CONF", "BAF_DRIFT",
       "SIZE_IN_KB", "PROBE_DENSITY_IN_KB"};
  public static final String[] OPT_QC_HEADS =
      {"targetConcordance", "actualConcordance", Array.toStr(QC_Thresholds), "goodCalls",
       "callsPassingFilter", "totalCalls", "pairsOfDuplicates", "CN"};

  public static OptimizedQCThresholds[] getNewOptqcs(double targetConcordance, int number) {
    OptimizedQCThresholds[] optqcs = new OptimizedQCThresholds[number];
    for (int i = 0; i < number; i++) {
      optqcs[i] = new OptimizedQCThresholds(targetConcordance, i);
    }
    return optqcs;
  }

  public static OptimizedQCThresholds loadThresholdsFromTxt(String QCThresholdFileName,
                                                            Logger log) {
    OptimizedQCThresholds thresholds = null;
    BufferedReader reader;
    String[] line, header;
    int[] indices;
    int index;
    try {
      reader = new BufferedReader(new FileReader(QCThresholdFileName));
      header = reader.readLine().trim().split("[\\s]+");
      indices = Array.intArray(QC_Thresholds.length, -1);
      for (int i = 0; i < header.length; i++) {
        index = ext.indexOfEndsWith(header[i], QC_Thresholds, true);
        if (index >= 0) {
          indices[index] = i;
        }
      }
      if (Array.min(indices) == -1) {
        log.reportError("Error - Need a column header ending with the following suffixes; missing at least one");
        log.reportError("        " + Array.toStr(QC_Thresholds, "  "));
        System.exit(1);
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        thresholds = parseLine(line, indices, log);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + QCThresholdFileName + "\" not found ");
      System.exit(1);
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + QCThresholdFileName + "\"");
      System.exit(2);
    }
    return thresholds;
  }

  private static OptimizedQCThresholds parseLine(String[] line, int[] indices, Logger log) {
    for (int i = 0; i < indices.length; i++) {
      if (line[indices[i]].equals("NA") && (QC_Thresholds[i].equals(QC_Thresholds[3])
                                            || QC_Thresholds[i].equals(QC_Thresholds[8]))) {
        log.reportError("Warning - not filtering with " + QC_Thresholds[i]);
        line[indices[i]] = "-1";
      } else if (line[indices[i]].equals("NA")) {
        line[indices[i]] = "NaN";
        log.reportError("Warning - not filtering with " + QC_Thresholds[i]);
      } else {
        log.report("Info - filtering with threshold " + QC_Thresholds[i] + "=" + line[indices[i]]);
      }
    }
    return new OptimizedQCThresholds(line, indices);
  }

  private double targetConcordance;
  private double actualConcordance;
  private final double alpha;
  private final double beastConfCutoff;
  private int callsPassingFilter;
  private int goodCalls;
  private int totalCalls;
  private int pairsOfDuplicates;
  private int CN;
  private final int numMarkers;
  private int numSampleCNVs;
  private final double lrrCutoff;
  private final double BAFQCcutoff;
  private final double twopqCutoff;
  private double hetCutoff;
  private double GCWF;
  private double cnvLRRSTDev;
  private double sampleCallRate;
  private double pennConf;


  private double bafDrift;

  private double kbSize;

  private double kbDensity;

  public OptimizedQCThresholds(double alpha, double beastConfCutoff, double lrrCutoff,
                               int numMarkers, double BAFQCcutoff, double twopqCutoff,
                               double hetCutoff, double GCWF, int numSampleCNVs, double cnvLRRSTDev,
                               double sampleCallRate, double pennConf, double bafDrift,
                               double kbSize, double kbDensity) {
    this.alpha = alpha;
    this.beastConfCutoff = beastConfCutoff;
    this.lrrCutoff = lrrCutoff;
    this.numMarkers = numMarkers;
    this.BAFQCcutoff = BAFQCcutoff;
    this.twopqCutoff = twopqCutoff;
    this.hetCutoff = hetCutoff;
    this.GCWF = GCWF;
    this.numSampleCNVs = numSampleCNVs;
    this.cnvLRRSTDev = cnvLRRSTDev;
    this.sampleCallRate = sampleCallRate;
    this.pennConf = pennConf;
    this.bafDrift = bafDrift;
    this.kbSize = kbSize;
    this.kbDensity = kbDensity;
  }

  public OptimizedQCThresholds(double targetConcordance, int CN) {
    this.targetConcordance = targetConcordance;
    actualConcordance = 0;
    alpha = 0;
    beastConfCutoff = 0;
    callsPassingFilter = 0;
    goodCalls = 0;
    totalCalls = 0;
    pairsOfDuplicates = 0;
    this.CN = CN;
    lrrCutoff = 0;
    numMarkers = 0;
    BAFQCcutoff = 0;
    twopqCutoff = 0;

  }

  public OptimizedQCThresholds(OptimizedQCThresholds comparedQCs, double targetConcordance,
                               double actualConcordance, int callsPassingFilter, int goodCalls,
                               int totalCalls, int pairsOfDuplicates, int CN) {
    this.targetConcordance = targetConcordance;
    this.actualConcordance = actualConcordance;
    this.callsPassingFilter = callsPassingFilter;
    this.goodCalls = goodCalls;
    this.totalCalls = totalCalls;
    this.pairsOfDuplicates = pairsOfDuplicates;
    this.CN = CN;
    alpha = comparedQCs.getAlpha();
    beastConfCutoff = comparedQCs.getBeastConfCutoff();
    lrrCutoff = comparedQCs.getLrrCutoff();
    numMarkers = comparedQCs.getNumMarkers();
    BAFQCcutoff = comparedQCs.getBAFQCcutoff();
    twopqCutoff = comparedQCs.getTwopqCutoff();
    hetCutoff = comparedQCs.getHetCutoff();
    GCWF = comparedQCs.getGCWF();
    numSampleCNVs = comparedQCs.getNumSampleCNVs();
    cnvLRRSTDev = comparedQCs.getCnvLRRSTDev();
    sampleCallRate = comparedQCs.getSampleCallRate();
    pennConf = comparedQCs.getPennConf();
    bafDrift = comparedQCs.getBafDrift();
    kbSize = comparedQCs.getKbSize();
    kbDensity = comparedQCs.getKbDensity();

  }

  public OptimizedQCThresholds(String[] line, int[] indices) {
    // "ALPHA", "BEAST_CONF", "SAMPLE_LRR_SD", "CNV_NUM_MARKERS", "BAF_QC", "TWOPQ", "PERCENT_HET",
    // "GCWF",
    // "NUM_SAMPLE_CNVS", "CNV_LRR_SD", "SAMPLE_CALL_RATE", "PENN_CONF", "BAF_DRIFT" };
    // return new OptimizedQCThresholds(Double.parseDouble(line[indices[0]]),
    // Double.parseDouble(line[indices[1]]),
    // Double.parseDouble(line[indices[2]]), Integer.parseInt(line[indices[3]]),
    // Double.parseDouble(line[indices[4]]),
    // Double.parseDouble(line[indices[5]]), Double.parseDouble(line[indices[6]]),
    // Double.parseDouble(line[indices[7]]),
    // Integer.parseInt(line[indices[8]]), Double.parseDouble(line[indices[9]]),
    // Double.parseDouble(line[indices[10]]),
    // Double.parseDouble(line[indices[11]]));
    alpha = Double.parseDouble(line[indices[0]]);
    beastConfCutoff = Double.parseDouble(line[indices[1]]);
    lrrCutoff = Double.parseDouble(line[indices[2]]);
    numMarkers = Integer.parseInt(line[indices[3]]);
    BAFQCcutoff = Double.parseDouble(line[indices[4]]);
    twopqCutoff = Double.parseDouble(line[indices[5]]);
    hetCutoff = Double.parseDouble(line[indices[6]]);
    GCWF = Double.parseDouble(line[indices[7]]);
    numSampleCNVs = Integer.parseInt(line[indices[8]]);
    cnvLRRSTDev = Double.parseDouble(line[indices[9]]);
    sampleCallRate = Double.parseDouble(line[indices[10]]);
    pennConf = Double.parseDouble(line[indices[11]]);
    bafDrift = Double.parseDouble(line[indices[12]]);
    kbSize = Double.parseDouble(line[indices[13]]);
    kbDensity = Double.parseDouble(line[indices[14]]);
  }

  public double getActualConcordance() {
    return actualConcordance;
  }

  public double getAlpha() {
    return alpha;
  }

  public double getBafDrift() {
    return bafDrift;
  }

  public double getBAFQCcutoff() {
    return BAFQCcutoff;
  }

  public double getBeastConfCutoff() {
    return beastConfCutoff;
  }

  public int getCallsPassingFilter() {
    return callsPassingFilter;
  }

  public int getCN() {
    return CN;
  }

  public double getCnvLRRSTDev() {
    return cnvLRRSTDev;
  }

  public String getDisplayString() {
    if (callsPassingFilter > 0) {
      return targetConcordance + "\t" + actualConcordance + "\t" + alpha + "\t" + beastConfCutoff
             + "\t" + lrrCutoff + "\t" + numMarkers + "\t" + BAFQCcutoff + "\t" + twopqCutoff + "\t"
             + hetCutoff + "\t" + GCWF + "\t" + numSampleCNVs + "\t" + cnvLRRSTDev + "\t"
             + sampleCallRate + "\t" + pennConf + "\t" + bafDrift + "\t" + kbSize + "\t" + kbDensity
             + "\t" + goodCalls + "\t" + callsPassingFilter + "\t" + totalCalls + "\t"
             + pairsOfDuplicates + "\t" + CN;
    } else if (CN == 5) {
      return "";
    } else {
      return "No parameters achieved target concordance of " + targetConcordance
             + " for copy number " + CN;
    }
  }

  public double getGCWF() {
    return GCWF;
  }

  public int getGoodCalls() {
    return goodCalls;
  }

  public double getHetCutoff() {
    return hetCutoff;
  }

  public double getKbDensity() {
    return kbDensity;
  }

  public double getKbSize() {
    return kbSize;
  }

  public double getLrrCutoff() {
    return lrrCutoff;
  }

  public int getNumMarkers() {
    return numMarkers;
  }

  public int getNumSampleCNVs() {
    return numSampleCNVs;
  }

  public int getPairsOfDuplicates() {
    return pairsOfDuplicates;
  }

  public double getPennConf() {
    return pennConf;
  }

  public double getSampleCallRate() {
    return sampleCallRate;
  }

  public double getTargetConcordance() {
    return targetConcordance;
  }

  public int getTotalCalls() {
    return totalCalls;
  }

  public double getTwopqCutoff() {
    return twopqCutoff;
  }

  public boolean moreCalls(int numCalls) {
    return callsPassingFilter < numCalls;
  }

  public void setCN(int cN) {
    CN = cN;
  }

}
