package org.genvisis.seq.manage;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import javax.annotation.Nonnull;
import org.genvisis.seq.ReferenceGenome;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.filesys.Segment;
import com.google.common.base.Enums;
import com.google.common.base.Optional;
import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
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
  private static final String[] COUNT_HEADER = Arrays.stream(Base.values()).map(Base::toString)
                                                     .map(b -> "NUM" + b).toArray(String[]::new);
  private static final String[] EXT_HEADER = new String[] {"MAPQ", "PHRED"};

  public static final ImmutableSet<Base> INDELS = Sets.immutableEnumSet(Base.DEL, Base.INS);
  public static final ImmutableSet<Base> NUCLEOTIDES = Arrays.stream(Base.values())
                                                             .filter(Predicates.not(INDELS::contains))
                                                             .collect(Sets.toImmutableEnumSet());

  public static enum Base {
    A, G, C, T, N, DEL, INS;
  }

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
    counts = new int[Base.values().length]; // Indexed by Base ordinals
    avgMapQ = new double[counts.length];
    avgPhread = new double[counts.length];
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

  /**
   * @return int[] of base counts, indexed by ordinal of {@link Base}
   */
  public int[] getCounts() {
    return counts;
  }

  /**
   * @param base
   * @return count of given base
   */
  public int getBaseCount(@Nonnull Base base) {
    return counts[base.ordinal()];
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
      INDELS.stream().map(Base::ordinal).forEach(i -> indelmask[i] = false);
    }
    if (!includeNs) {
      indelmask[Base.N.ordinal()] = false;
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
    for (int i = 0; i < counts.length; i++) {
      out += "\t" + counts[i];
    }
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

    Base ref = Enums.getIfPresent(Base.class, refAllele).orNull();
    if (ref != null) referenceMask[ref.ordinal()] = false;

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

    Base ref = Enums.getIfPresent(Base.class, refAllele).orNull();
    if (ref != null) referenceMask[ref.ordinal()] = true;

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
    for (Base base : NUCLEOTIDES) {
      totalDepth += counts[base.ordinal()];
      mapQ += avgMapQ[base.ordinal()];
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
    counts[Base.DEL.ordinal()]++;
    avgPhread[Base.DEL.ordinal()] += p;
    avgMapQ[Base.DEL.ordinal()] += mapQ;
  }

  private void addIns(int mapQ, double p) {
    counts[Base.INS.ordinal()]++;
    avgPhread[Base.INS.ordinal()] += p;
    avgMapQ[Base.INS.ordinal()] += mapQ;
  }

  private void addRegBase(String b, int mapQ, double p, Logger log) {
    Optional<Base> base = Enums.getIfPresent(Base.class, b);
    if (base.isPresent() && NUCLEOTIDES.contains(base.get())) {
      int i = base.get().ordinal();
      counts[i]++;
      avgPhread[i] += p;
      avgMapQ[i] += mapQ;
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
