package org.genvisis.seq.manage;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.filesys.Segment;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Stores the pileUp of a bam file for each position found (by bin)
 */
public class BamPile extends Segment implements Serializable {

  /**
   *
   */
  private static final long serialVersionUID = 1L;
  private static final String[] BASE_HEADER = new String[] {"UCSC", "REF", "NUM_REF", "NUM_ALT",
                                                            "PROP_REF", "PROP_ALT"};
  private static final String[] COUNT_HEADER = new String[] {"NUMA", "NUMG", "NUMC", "NUMT", "NUMN",
                                                             "NUMDEL", "NUMINS"};
  private static final String[] EXT_HEADER = new String[] {"MAPQ", "PHRED"};

  private final Segment bin;
  private final int[] counts;
  private String refAllele;
  private int numBasesOverlap;
  private int numBasesWithMismatch;
  private int numOverlappingReads;

  private final double[] avgMapQ;
  private final double[] avgPhread;
  private double overallAvgMapQ;
  private double overallAvgDepth;

  public static String getBampPileHeader() {
    String header = ArrayUtils.toStr(BASE_HEADER);
    for (String element : COUNT_HEADER) {
      header += "\t" + element;
    }

    for (String element : COUNT_HEADER) {
      for (String element2 : EXT_HEADER) {
        header += "\t" + element + "_" + element2;
      }
    }

    return header;
  }

  public BamPile(Segment bin) {
    super(bin.getChr(), bin.getStart(), bin.getStop());
    this.bin = bin;
    counts = new int[7];// A,G,C,T,N,Del,Ins
    avgMapQ = new double[7];
    avgPhread = new double[7];
    overallAvgMapQ = 0;
    overallAvgDepth = 0;
    refAllele = "NA";
    numBasesOverlap = 0;
    numBasesWithMismatch = 0;
    numOverlappingReads = 0;

  }

  public double[] getAvgMapQ() {
    return avgMapQ;
  }

  public int[] getCounts() {
    return counts;
  }

  public int getNumBasesOverlap() {
    return numBasesOverlap;
  }

  public int getNumBasesWithMismatch() {
    return numBasesWithMismatch;
  }

  public int getTotalDepth(boolean includeIndels, boolean includeNs) {
    boolean[] indelmask = new boolean[7];
    Arrays.fill(indelmask, true);
    if (!includeIndels) {
      indelmask[6] = false;
      indelmask[5] = false;
    }
    if (!includeNs) {
      indelmask[4] = false;
    }

    return ArrayUtils.sum(ArrayUtils.subArray(counts, indelmask));

  }

  public void setReference(ReferenceGenome referenceGenome) {

    String[] ref = referenceGenome.getSequenceFor(bin);
    refAllele = ref[0];
    if (ref.length > 1) {
      System.out.println(ArrayUtils.toStr(ref));
    }

    if (ref.length > 1) {
      for (int i = 1; i < ref.length; i++) {
        refAllele += "/" + ref[i];
      }
    }
  }

  public String getOuput(Logger log) {
    double numRef = getNumRef(log);
    double numAlt = getNumAlt(log);
    String out = "";
    out += bin.getUCSClocation();
    out += "\t" + refAllele;
    out += "\t" + numRef;
    out += "\t" + numAlt;
    out += "\t" + getPropRef(log);
    out += "\t" + getPropAlt(log);
    out += "\t" + counts[0];
    out += "\t" + counts[1];
    out += "\t" + counts[2];
    out += "\t" + counts[3];
    out += "\t" + counts[4];
    out += "\t" + counts[5];
    out += "\t" + counts[6];
    out += "\t" + ArrayUtils.toStr(avgMapQ);
    out += "\t" + ArrayUtils.toStr(avgPhread);
    return out;
  }

  public double getPropRef(Logger log) {
    double numRef = getNumRef(log);
    double numAlt = getNumAlt(log);
    double total = numRef + numAlt;
    return numRef / total;
  }

  private double getPropAlt(Logger log) {
    double numRef = getNumRef(log);
    double numAlt = getNumAlt(log);
    double total = numRef + numAlt;
    return numAlt / total;
  }

  private int[] getAltCounts(Logger log) {
    boolean[] referenceMask = new boolean[7];
    Arrays.fill(referenceMask, true);
    referenceMask[6] = false;
    referenceMask[5] = false;

    if (refAllele.equals("A")) {
      referenceMask[0] = false;
    } else if (refAllele.equals("G")) {
      referenceMask[1] = false;
    } else if (refAllele.equals("C")) {
      referenceMask[2] = false;
    } else if (refAllele.equals("T")) {
      referenceMask[3] = false;
    } else if (refAllele.equals("N")) {
      referenceMask[4] = false;
    }
    if (ArrayUtils.booleanArraySum(referenceMask) != 4) {
      log.reportError("Invalid number of alternate allele possibilities, found "
                      + ArrayUtils.booleanArraySum(referenceMask) + " with ref allele" + refAllele);

    }
    return ArrayUtils.subArray(counts, referenceMask);
  }

  public int getNumAlt(Logger log) {
    return ArrayUtils.sum(getAltCounts(log));
  }

  public boolean hasOnlyOneAlt(Logger log) {
    return ArrayUtils.countIf(getAltCounts(log), 0) == 3;// everthing else is 0
  }

  public int getNumRef(Logger log) {
    boolean[] referenceMask = new boolean[7];
    Arrays.fill(referenceMask, false);
    referenceMask[6] = false;
    referenceMask[5] = false;

    if (refAllele.equals("A")) {
      referenceMask[0] = true;
    } else if (refAllele.equals("G")) {
      referenceMask[1] = true;
    } else if (refAllele.equals("C")) {
      referenceMask[2] = true;
    } else if (refAllele.equals("T")) {
      referenceMask[3] = true;
    } else if (refAllele.equals("N")) {
      referenceMask[4] = true;
    }
    if (ArrayUtils.booleanArraySum(referenceMask) != 1) {
      log.reportError("Invalid number of reference allele possibilities");
    }
    return (ArrayUtils.sum(ArrayUtils.subArray(counts, referenceMask)));
  }

  public boolean hasAltAllele(Logger log) {
    return getNumAlt(log) > 0;
  }

  public void summarize() {

    int totalDepth = 0;
    double mapQ = 0;
    // TODO count deletions and insertions here?
    for (int i = 0; i < counts.length - 2; i++) {
      totalDepth += counts[i];
      mapQ += avgMapQ[i];
    }
    overallAvgDepth = (double) totalDepth / bin.getSize();
    if (totalDepth > 0) {
      overallAvgMapQ = mapQ / totalDepth;
    }
    for (int i = 0; i < counts.length; i++) {
      if (counts[i] > 0) {
        avgMapQ[i] = avgMapQ[i] / counts[i];
        avgPhread[i] = avgPhread[i] / counts[i];
      }
    }
  }

  public double getOverallAvgMapQ() {
    return overallAvgMapQ;
  }

  public double getOverallAvgDepth() {
    return overallAvgDepth;
  }

  public Segment getBin() {
    return bin;
  }

  public int getNumOverlappingReads() {
    return numOverlappingReads;
  }

  public void addRecord(SAMRecord samRecord, String[] refMatchedSegment, double phredFilter,
                        Logger log) {
    Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
    String[] ref = refMatchedSegment;
    Segment toPile = bin.getIntersection(samRecordSegment, log);
    int mapQ = samRecord.getMappingQuality();
    if (mapQ == 255) {
      String error = "Detected invalid mapping quality (255)";
      log.reportError(error);
      throw new IllegalArgumentException(error);
    } else {
      // System.out.println(toPile.getSize());
      String r = samRecord.getReadString();
      double[] p = SamRecordOps.getReadPhred(samRecord);
      // double[] p = new double[r.length()];
      numOverlappingReads++;
      int curRefBase = samRecord.getAlignmentStart();
      int curReadBase = 0;
      List<CigarElement> cigarEls = samRecord.getCigar().getCigarElements();
      for (CigarElement cigarEl : cigarEls) {
        if (curRefBase > toPile.getStop()) {
          break;

        }
        CigarOperator op = cigarEl.getOperator();
        for (int i = 0; i < cigarEl.getLength(); i++) {
          String base = null;
          if (curRefBase > toPile.getStop()) {
            break;

          }
          if (curRefBase >= toPile.getStart() && curRefBase <= toPile.getStop()) {
            if (p[curReadBase] > phredFilter) {
              base = r.charAt(curReadBase) + "";
            }
          }

          if (op.consumesReadBases() && op.consumesReferenceBases()) {
            if (base != null) {
              addRegBase(base, mapQ, p[curReadBase], log);
              numBasesOverlap++;// we only count aligned bases
              if (ref != null && !base.equals(ref[curRefBase - bin.getStart()])) {
                numBasesWithMismatch++;
              }
            }
            curRefBase++;
            curReadBase++;
          } else if (op.consumesReadBases()) {
            if (base != null) {
              addIns(mapQ, p[curReadBase]);
            }
            curReadBase++;
          } else if (op.consumesReferenceBases()) {
            if (base != null) {
              addDel(mapQ, p[curReadBase]);
            }
            curRefBase++;
          }
        }
      }

      if (numBasesWithMismatch > numBasesOverlap) {
        String error = "Internal Error: Found more reads with mismatching bases than total reads";
        log.reportError(error);
        throw new IllegalStateException(error);
      }
    }
  }

  private void addDel(int mapQ, double p) {
    counts[5]++;
    avgPhread[5] += p;
    avgMapQ[5] += mapQ;
  }

  private void addIns(int mapQ, double p) {
    counts[6]++;
    avgPhread[6] += p;
    avgMapQ[6] += mapQ;
  }

  private void addRegBase(String b, int mapQ, double p, Logger log) {
    if (b.equals("A")) {
      counts[0]++;
      avgPhread[0] += p;
      avgMapQ[0] += mapQ;
    } else if (b.equals("G")) {
      counts[1]++;
      avgPhread[1] += p;
      avgMapQ[1] += mapQ;
    } else if (b.equals("C")) {
      counts[2]++;
      avgPhread[2] += p;
      avgMapQ[2] += mapQ;
    } else if (b.equals("T")) {
      counts[3]++;
      avgPhread[3] += p;
      avgMapQ[3] += mapQ;
    } else if (b.equals("N")) {
      counts[4]++;
      avgPhread[4] += p;
      avgMapQ[4] += mapQ;
    } else {
      String error = "Invalid base " + b + " for regular base tracking";
      log.reportError(error);
      throw new IllegalArgumentException(error);
    }
  }

  public static void writeSerial(BamPile[] bamPiles, String filename) {
    SerializedFiles.writeSerial(bamPiles, filename, true);
  }

  public static BamPile[] readSerial(String filename, Logger log) {
    return (BamPile[]) SerializedFiles.readSerial(filename, log, false, true);
  }

}
