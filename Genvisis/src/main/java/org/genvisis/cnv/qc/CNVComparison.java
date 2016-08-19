package org.genvisis.cnv.qc;

import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;

// class to determine the percent concordance given QC parameters, currently only alpha, confidence,
// numMarkers, and LRR_SD

public class CNVComparison {
  static final String[] QC_PARAMETERs = {"0 - No Filtering ", " 1 - CNVQC and Sample QC",
                                         "2 -Sample QC only"};
  private final CNVariantQC[][] filteredcnvQCs1;
  private CNVariantQC[][] filteredcnvQCs2;
  private ArrayList<CNVariantQC> misses;;
  private final Hashtable<String, CNVSampleQC> cnvSampleQCHash;
  private final OptimizedQCThresholds qcThresholds;
  private double[] averageCNPercent;
  private int[] goodCalls;
  private int[] totalCallsAvailable;
  private int[] filteredCallsAvailable;
  private int numberIndsCompared;
  private final int filterType;
  private Logger log;

  // TODO undefined sampleQCFile, or missing QC info

  public CNVComparison(CNVariantQC[][] unfilteredcnvsQCs1, CNVariantQC[][] unfilteredcnvsQCs2,
                       Hashtable<String, CNVSampleQC> cnvSampleQCHash,
                       OptimizedQCThresholds qcThresholds, int filterType, Logger log) {
    totalCallsAvailable = getCallsAvailable(unfilteredcnvsQCs1, unfilteredcnvsQCs2);
    this.cnvSampleQCHash = cnvSampleQCHash;
    this.qcThresholds = qcThresholds;
    this.filterType = filterType;
    misses = new ArrayList<CNVariantQC>();
    filteredcnvQCs1 = getFilteredCnvsByType(unfilteredcnvsQCs1);
    filteredcnvQCs2 = getFilteredCnvsByType(unfilteredcnvsQCs2);
    filteredCallsAvailable = getCallsAvailable(filteredcnvQCs1, filteredcnvQCs2);
    goodCalls = getGoodCalls(filteredcnvQCs1, filteredcnvQCs2);
    averageCNPercent = averageCNPercents(goodCalls, filteredCallsAvailable);
    numberIndsCompared = getNumberCompared(filteredcnvQCs1, filteredcnvQCs2);

  }

  public CNVComparison(CNVariantQC[][] unfilteredcnvsQCs1,
                       Hashtable<String, CNVSampleQC> cnvSampleQCHash,
                       OptimizedQCThresholds qcThresholds, int filterType, Logger log) {
    this.cnvSampleQCHash = cnvSampleQCHash;
    this.qcThresholds = qcThresholds;
    this.filterType = filterType;
    filteredcnvQCs1 = getFilteredCnvsByType(unfilteredcnvsQCs1);

  }

  public ArrayList<CNVariantQC> getMisses() {
    return misses;
  }

  public CNVariantQC[][] getFilteredcnvQCs1() {
    return filteredcnvQCs1;
  }

  public OptimizedQCThresholds getQcThresholds() {
    return qcThresholds;
  }

  public int getNumberIndsCompared() {
    return numberIndsCompared;
  }

  public double[] getAverageCNPercent() {
    return averageCNPercent;
  }

  public int[] getGoodCalls() {
    return goodCalls;
  }

  public int[] getTotalCallsAvailable() {
    return totalCallsAvailable;
  }

  public int[] getFilteredCallsAvailable() {
    return filteredCallsAvailable;
  }

  private int getNumberCompared(CNVariantQC[][] filteredcnvQC1s, CNVariantQC[][] filteredcnvQC2s) {
    int count = 0;
    for (int i = 0; i < filteredcnvQC1s.length; i++) {
      if (filteredcnvQC1s[i] != null && filteredcnvQC2s[i] != null) {
        if (filteredcnvQC1s[i].length > 0 && filteredcnvQC2s[i].length > 0) {
          count++;
        }
      }
    }
    return count;
  }

  private int[] getGoodCalls(CNVariantQC[][] filteredcnvQCs1, CNVariantQC[][] filteredcnvQCs2) {
    int[][] counts = new int[6][5];
    if (filteredcnvQCs1.length == filteredcnvQCs2.length) {
      for (int i = 0; i < filteredcnvQCs1.length; i++) {
        if (filteredcnvQCs1[i] != null && filteredcnvQCs2[i] != null) {
          if (filteredcnvQCs1[i].length > 0 && filteredcnvQCs2[i].length > 0) {
            counts = countMatches(filteredcnvQCs1[i], filteredcnvQCs2[i], counts);
          }
        }
      }
    } else {
      log.reportError("Error - unmatched number of individuals are being compared, this should not happen");
      System.exit(1);
    }
    return goodCallCN(counts);
  }

  private double[] averageCNPercents(int[] goodCalls, int[] cnNumbers) {
    double[] averageCNPercent = new double[6];
    for (int i = 0; i < goodCalls.length; i++) {
      if (goodCalls[i] > 0) {
        averageCNPercent[i] = (double) goodCalls[i] / (double) cnNumbers[i];
      } else {
        averageCNPercent[i] = 0;
      }
    }
    return averageCNPercent;
  }

  // currently defined as 2*exact matches + significant overlap matches
  private int[] goodCallCN(int[][] counts) {
    int[] goodCalls = new int[6];
    for (int i = 0; i < counts.length; i++) {
      goodCalls[i] += 2 * counts[i][3] + counts[i][4];
    }
    return goodCalls;

  }

  private int[] getCallsAvailable(CNVariantQC[][] cnvQCs1, CNVariantQC[][] cnvsQC2) {
    int[] callsAvailable = new int[6];
    callsAvailable = getAllCalls(cnvQCs1, callsAvailable);
    callsAvailable = getAllCalls(cnvsQC2, callsAvailable);
    return callsAvailable;
  }

  private int[] getAllCalls(CNVariantQC[][] cnvQCs, int[] totalCallsAvailable) {
    for (CNVariantQC[] cnvQC : cnvQCs) {
      if ((cnvQC != null)) {
        if (cnvQC.length > 0) {
          totalCallsAvailable = getCNNumbers(cnvQC, totalCallsAvailable);
        }
      }
    }
    return totalCallsAvailable;
  }

  private CNVariantQC[][] getFilteredCnvsByType(CNVariantQC[][] unfilteredcnvQCs) {
    CNVariantQC[][] filteredcnvs = new CNVariantQC[unfilteredcnvQCs.length][];
    checkCompatability();

    for (int i = 0; i < unfilteredcnvQCs.length; i++) {
      if (unfilteredcnvQCs[i] == null) {
        log.reportError("Warning - Some individuals do not have cnvs");
      }
      if (unfilteredcnvQCs[i].length > 0) {
        if (filterType == 0) {
          filteredcnvs[i] = unfilteredcnvQCs[i];
        } else if (filterType == 1) {
          if (passesSampleQC(unfilteredcnvQCs[i][0], unfilteredcnvQCs[i].length)) {
            filteredcnvs[i] = filterCNVQC(unfilteredcnvQCs[i]);
          }
        } else if (filterType == 2) {
          if (passesSampleQC(unfilteredcnvQCs[i][0], unfilteredcnvQCs[i].length)) {
            filteredcnvs[i] = unfilteredcnvQCs[i];

          }
        } else {
          log.reportError("Error - invaled filter Type , need to use "
                          + Array.toStr(QC_PARAMETERs));
          System.exit(1);
        }
      }
    }
    return filteredcnvs;
  }

  private void checkCompatability() {
    if (cnvSampleQCHash == null && (filterType == 1 || filterType == 2)) {
      log.reportError("Error - invaled filter Type with undefined sample QC file, need to use "
                      + QC_PARAMETERs[0]);
      System.exit(1);
    }
  }

  private boolean passesSampleQC(CNVariantQC unfilteredcnvQC, int numCNVs) {
    String lookup = unfilteredcnvQC.getCnVariant().getFamilyID() + "\t"
                    + unfilteredcnvQC.getCnVariant().getIndividualID();
    boolean passSampleQC = false;
    if (cnvSampleQCHash.containsKey(lookup)) {
      CNVSampleQC cnvSampleQC = cnvSampleQCHash.get(lookup);
      if (lrr_SD(cnvSampleQC.getLrrSDev()) && gcwf(cnvSampleQC.getGCWF()) && numberCNVs(numCNVs)
          && callRate(unfilteredcnvQC.getSampleCallRate()) && bafDrift(cnvSampleQC.getBafDrift())) {
        passSampleQC = true;

      }
    } else {
      log.reportError("Error - Sample QC was defined, but Sample QC was not available for "
                      + lookup);
      System.exit(1);
    }
    return passSampleQC;
  }

  private CNVariantQC[] filterCNVQC(CNVariantQC[] unfilteredcnvQCs) {
    ArrayList<CNVariantQC> cnvQCs = new ArrayList<CNVariantQC>();

    for (CNVariantQC unfilteredcnvQC : unfilteredcnvQCs) {
      // double[] mafs = unfilteredcnvQCs[i].getMafs();
      // double[] bafs = unfilteredcnvQCs[i].getBafs();
      // double[] lrrs = unfilteredcnvQCs[i].getLrrs();
      // byte[] abGenoytypes = unfilteredcnvQCs[i].getGenotypes();
      double pennConf = unfilteredcnvQC.getCnVariant().getScore();
      int numMarkers = unfilteredcnvQC.getCnVariant().getNumMarkers();
      double height = unfilteredcnvQC.getHeight();
      double kbSize = (double) (unfilteredcnvQC.getCnVariant().getSize()) / 1000;
      double kbDensity = (numMarkers) / kbSize;
      // System.out.println(kbSize + "\t" + qcThresholds.getKbSize());
      // double sumTwopq = 0;
      // double sumBaf = 0;
      // int numPoly = 0;
      // int numHets = 0;
      // for (int j = 0; j < mafs.length; j++) {
      // if (mafs[j] > 0) {
      // sumTwopq += ((2 * mafs[j]) * (1 - mafs[j]));
      // numPoly++;
      // sumBaf += bafDistance(bafs[j]);
      // if (abGenoytypes[j] == 1) {
      // numHets++;
      // }
      // }
      // }
      if (goodScore(numMarkers, height) && checkKbSize(kbSize) && checkKbDensity(kbDensity)
          && goodNumMarkers(numMarkers) && checkPennConf(pennConf)) {
        cnvQCs.add(unfilteredcnvQC);
      }
    }
    return CNVariantQC.toCNVQCArray(cnvQCs);
  }

  private boolean checkKbDensity(double kbDensity) {
    if (Double.isNaN(qcThresholds.getKbDensity())) {
      return true;
    } else {
      return kbDensity >= qcThresholds.getKbDensity();
    }
  }

  private boolean checkKbSize(double kbSize) {
    if (Double.isNaN(qcThresholds.getKbSize())) {
      return true;
    } else {
      return kbSize >= qcThresholds.getKbSize();
    }
  }

  private boolean bafDrift(double sampleBafDrift) {
    if (Double.isNaN(qcThresholds.getBafDrift())) {
      return true;
    } else {
      return sampleBafDrift <= qcThresholds.getBafDrift();
    }
  }

  private boolean callRate(double sampleCallRate) {
    if (Double.isNaN(qcThresholds.getSampleCallRate())) {
      return true;
    } else {
      return sampleCallRate >= qcThresholds.getSampleCallRate();
    }
  }

  private boolean numberCNVs(int samplenumCNVs) {
    if (qcThresholds.getNumSampleCNVs() < 0) {
      return true;
    } else {
      return samplenumCNVs <= qcThresholds.getNumSampleCNVs();
    }
  }

  private boolean gcwf(double sampleGCWF) {
    if (Double.isNaN(qcThresholds.getGCWF())) {
      return true;
    } else {
      return Math.abs(sampleGCWF) <= qcThresholds.getGCWF();
    }
  }

  private boolean lrr_SD(double sampleLrrSD) {
    if (Double.isNaN(qcThresholds.getLrrCutoff())) {
      return true;
    } else {
      return sampleLrrSD <= qcThresholds.getLrrCutoff();
    }
  }

  private boolean goodNumMarkers(int numMarkers) {
    if (qcThresholds.getNumMarkers() < 0) {
      return true;
    } else {
      return numMarkers >= qcThresholds.getNumMarkers();
    }
  }


  private boolean goodScore(int numMarkers, double height) {
    if (Double.isNaN(qcThresholds.getAlpha()) || Double.isNaN(qcThresholds.getBeastConfCutoff())) {
      return true;
    } else {
      return score(numMarkers, qcThresholds.getAlpha(),
                   height) >= qcThresholds.getBeastConfCutoff();
    }
  }

  private double score(int numMarkers, double alpha, double cnvHeight) {
    return Math.abs(Math.pow(numMarkers, alpha) * cnvHeight);

  }

  private boolean checkPennConf(double cnvPennConf) {
    if (Double.isNaN(qcThresholds.getPennConf())) {
      return true;
    } else {
      return cnvPennConf >= qcThresholds.getPennConf();
    }
  }
  //
  // private boolean checkHets(int numHets, int numPoly) {
  // if (Double.isNaN(qcThresholds.getHetCutoff())) {
  // return true;
  // } else {
  // return (double) (numHets / numPoly) <= qcThresholds.getHetCutoff();
  // }
  // }
  //
  // private boolean aPoly(int numPoly) {
  // return numPoly > 0;
  // }
  //
  //
  // private double bafDistance(double baf) {
  // return Maths.min(baf, 1 - baf);
  // }

  private int[][] countMatches(CNVariantQC[] filteredcnvQCs1, CNVariantQC[] filteredcnvQCs2,
                               int[][] counts) {
    // CN states as defined here are ,0,1,2,3,4,5...where 2 is a total non-copy number aware count,
    // and 5 is for non matching (overlap,exact,sigoverlap)
    int match;
    int CN = 5;
    for (CNVariantQC element : filteredcnvQCs1) {
      match = 0;
      for (CNVariantQC element2 : filteredcnvQCs2) {
        if (element2.getCnVariant().getCN() == element.getCnVariant().getCN()) {
          if (element.getCnVariant().equals(element2.getCnVariant())) {
            match = 3;
            element2.getCnVariant().setSource(99);
            CN = element2.getCnVariant().getCN();
          } else if (match < 2
                     && element.getCnVariant().significantOverlap(element2.getCnVariant())) {
            match = 4;
            CN = element2.getCnVariant().getCN();

            // an overlap is not assumed in both directions
            // filteredcnvQCs2[b].getCnVariant().setSource(99);
          } else if (match < 2 && element.getCnVariant().overlaps(element2.getCnVariant())) {
            match = 2;
            CN = element2.getCnVariant().getCN();
            element2.getCnVariant().setSource(99);
          }
        }
      }
      if (match == 0) {
        misses.add(element);
      }
      counts[CN][match]++;
      counts[2][match]++;
    }
    CN = 5;
    for (CNVariantQC element : filteredcnvQCs2) {
      match = 1;
      for (CNVariantQC element2 : filteredcnvQCs1) {
        if (element.getCnVariant().getCN() == element2.getCnVariant().getCN()) {
          if (element.getCnVariant().getSource() != 99
              && element.getCnVariant().equals(element2.getCnVariant())) {
            match = 3;
            CN = element.getCnVariant().getCN();
          } else if (match < 2 && element.getCnVariant().getSource() != 99
                     && element.getCnVariant().significantOverlap(element2.getCnVariant())) {
            match = 4;
            CN = element.getCnVariant().getCN();
          } else if (match < 2 && element.getCnVariant().getSource() != 99
                     && element.getCnVariant().overlaps(element2.getCnVariant())) {
            match = 2;
            CN = element.getCnVariant().getCN();
          }
        }
      }

      if (element.getCnVariant().getSource() != 99) {
        if (match == 1) {
          misses.add(element);
        }
        counts[CN][match]++;
        counts[2][match]++;
      }
    }
    return counts;
  }

  // TODO FID
  private int[] getCNNumbers(CNVariantQC[] cnvQCs, int[] cnNumbers) {
    String IID = cnvQCs[0].getCnVariant().getIndividualID();
    for (int i = 0; i < cnvQCs.length; i++) {
      if (IID.equals(cnvQCs[i].getCnVariant().getIndividualID())) {
        cnNumbers[cnvQCs[i].getCnVariant().getCN()]++;
        cnNumbers[2]++;
      } else {
        log.reportError("Error - found unmatched IDs for individuals CNVS, this should not happen");
        System.exit(1);
      }
    }
    return cnNumbers;
  }
}
