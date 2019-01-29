package org.genvisis.seq.manage;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;
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
import com.google.common.util.concurrent.AtomicDoubleArray;
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
  private String refAllele;

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

  private final AtomicInteger numOverlappingReads_A;
  private final AtomicInteger numBasesOverlap_A;
  private final AtomicInteger numBasesWithMismatch_A;
  private final AtomicIntegerArray counts_A;
  private final AtomicDoubleArray avgMapQ_A;
  private final AtomicDoubleArray avgPhread_A;

  public BamPile(Segment bin) {
    super(bin.getChr(), bin.getStart(), bin.getStop());
    this.bin = bin;
    counts_A = new AtomicIntegerArray(Base.values().length);
    avgMapQ_A = new AtomicDoubleArray(counts_A.length());
    avgPhread_A = new AtomicDoubleArray(counts_A.length());
    numOverlappingReads_A = new AtomicInteger(0);
    numBasesOverlap_A = new AtomicInteger(0);
    numBasesWithMismatch_A = new AtomicInteger(0);
    overallAvgMapQ = 0;
    overallAvgDepth = 0;
    refAllele = "NA";
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = super.hashCode();
    result = prime * result + avgMapQ_A.hashCode();
    result = prime * result + avgPhread_A.hashCode();
    result = prime * result + ((bin == null) ? 0 : bin.hashCode());
    result = prime * result + counts_A.hashCode();
    result = prime * result + numBasesOverlap_A.hashCode();
    result = prime * result + numBasesWithMismatch_A.hashCode();
    result = prime * result + numOverlappingReads_A.hashCode();
    long temp;
    temp = Double.doubleToLongBits(overallAvgDepth);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    temp = Double.doubleToLongBits(overallAvgMapQ);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    result = prime * result + ((refAllele == null) ? 0 : refAllele.hashCode());
    return result;
  }

  private static boolean eq(AtomicDoubleArray a, AtomicDoubleArray b) {
    if (a == null && b == null) return true;
    if (a == null && b != null) return false;
    if (a != null && b == null) return false;
    if (a.length() != b.length()) return false;
    for (int i = 0; i < a.length(); i++) {
      if (a.get(i) != b.get(i)) return false;
    }
    return true;
  }

  private static boolean eq(AtomicIntegerArray a, AtomicIntegerArray b) {
    if (a == null && b == null) return true;
    if (a == null && b != null) return false;
    if (a != null && b == null) return false;
    if (a.length() != b.length()) return false;
    for (int i = 0; i < a.length(); i++) {
      if (a.get(i) != b.get(i)) return false;
    }
    return true;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (!super.equals(obj)) return false;
    if (getClass() != obj.getClass()) return false;
    BamPile other = (BamPile) obj;
    if (!eq(avgMapQ_A, other.avgMapQ_A)) return false;
    if (!eq(avgPhread_A, other.avgPhread_A)) return false;
    if (bin == null) {
      if (other.bin != null) return false;
    } else if (!bin.matches(other.bin)) return false;
    if (!eq(counts_A, other.counts_A)) return false;
    if (numBasesOverlap_A.get() != other.numBasesOverlap_A.get()) return false;
    if (numBasesWithMismatch_A.get() != other.numBasesWithMismatch_A.get()) return false;
    if (numOverlappingReads_A.get() != other.numOverlappingReads_A.get()) return false;
    if (Double.doubleToLongBits(overallAvgDepth) != Double.doubleToLongBits(other.overallAvgDepth)) return false;
    if (Double.doubleToLongBits(overallAvgMapQ) != Double.doubleToLongBits(other.overallAvgMapQ)) return false;
    if (refAllele == null) {
      if (other.refAllele != null) return false;
    } else if (!refAllele.equals(other.refAllele)) return false;
    return true;
  }

  /**
   * @param base
   * @return count of given base
   */
  public int getBaseCount(@Nonnull Base base) {
    return counts_A.get(base.ordinal());
  }

  public int getNumBasesOverlap() {
    return numBasesOverlap_A.get();
  }

  public int getNumBasesWithMismatch() {
    return numBasesWithMismatch_A.get();
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

    int sum = 0;
    for (int i = 0; i < indelmask.length; i++) {
      if (indelmask[i]) {
        sum += counts_A.get(i);
      }
    }
    return sum;
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
    for (int i = 0; i < counts_A.length(); i++) {
      out += "\t" + counts_A.get(i);
    }
    out += "\t" + toStr(avgMapQ_A);
    out += "\t" + toStr(avgPhread_A);
    return out;
  }

  private String toStr(AtomicDoubleArray arr) {
    StringBuilder sb = new StringBuilder(Double.toString(arr.get(0)));
    for (int i = 1; i < arr.length(); i++) {
      sb.append("\t").append(arr.get(i));
    }
    return sb.toString();
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

    int cnt = 0;
    int[] altArray = new int[ArrayUtils.booleanArraySum(referenceMask)];
    for (int i = 0; i < referenceMask.length; i++) {
      if (referenceMask[i]) {
        altArray[cnt++] = counts_A.get(i);
      }
    }

    return altArray;
  }

  /**
   * @param log
   * @return int[] array with counts for {A, G, C, T, N}
   */
  public int[] getAlleleCounts(Logger log) {
    int[] array = new int[5];
    array[0] = counts_A.get(0);
    array[1] = counts_A.get(1);
    array[2] = counts_A.get(2);
    array[3] = counts_A.get(3);
    array[4] = counts_A.get(4);
    return array;
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
    int sum = 0;
    for (int i = 0; i < referenceMask.length; i++) {
      if (referenceMask[i]) {
        sum += counts_A.get(i);
      }
    }
    return sum;
  }

  public boolean hasAltAllele(Logger log) {
    return getNumAlt(log) > 0;
  }

  public void summarize() {
    int totalDepth = 0;
    double mapQ = 0;
    // TODO count deletions and insertions here?
    for (Base base : NUCLEOTIDES) {
      totalDepth += counts_A.get(base.ordinal());
      mapQ += avgMapQ_A.get(base.ordinal());
    }
    overallAvgDepth = (double) totalDepth / bin.getSize();
    if (totalDepth > 0) {
      overallAvgMapQ = mapQ / totalDepth;
    }
    for (int i = 0; i < counts_A.length(); i++) {
      if (counts_A.get(i) > 0) {
        avgMapQ_A.set(i, avgMapQ_A.get(i) / counts_A.get(i));
        avgPhread_A.set(i, avgPhread_A.get(i) / counts_A.get(i));
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
    return numOverlappingReads_A.get();
  }

  public void addRecordAtomic(SAMRecord samRecord, String[] refMatchedSegment, double phredFilter,
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
      numOverlappingReads_A.incrementAndGet();
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
              addRegBaseAtomic(base, mapQ, p[curReadBase], log);
              numBasesOverlap_A.incrementAndGet();// we only count aligned bases
              if (ref != null && !base.equals(ref[curRefBase - bin.getStart()])) {
                numBasesWithMismatch_A.incrementAndGet();
              }
            }
            curRefBase++;
            curReadBase++;
          } else if (op.consumesReadBases()) {
            if (base != null) {
              addInsAtomic(mapQ, p[curReadBase]);
            }
            curReadBase++;
          } else if (op.consumesReferenceBases()) {
            if (base != null) {
              addDelAtomic(mapQ, p[curReadBase]);
            }
            curRefBase++;
          }
        }
      }

      // TODO concurrency may cause issues here:
      int currMis, currOver;
      synchronized (this) {
        currMis = numBasesWithMismatch_A.get();
        currOver = numBasesOverlap_A.get();
      }
      if (currMis > currOver) {
        String error = "Internal Error: Found more reads with mismatching bases than total reads";
        log.reportError(error);
        throw new IllegalStateException(error);
      }
    }
  }

  private void addDelAtomic(int mapQ, double p) {
    counts_A.addAndGet(Base.DEL.ordinal(), 1);
    avgPhread_A.addAndGet(Base.DEL.ordinal(), p);
    avgMapQ_A.addAndGet(Base.DEL.ordinal(), mapQ);
  }

  private void addInsAtomic(int mapQ, double p) {
    avgPhread_A.addAndGet(Base.INS.ordinal(), p);
    avgMapQ_A.addAndGet(Base.INS.ordinal(), mapQ);
  }

  private void addRegBaseAtomic(String b, int mapQ, double p, Logger log) {
    Optional<Base> base = Enums.getIfPresent(Base.class, b);
    if (base.isPresent() && NUCLEOTIDES.contains(base.get())) {
      int i = base.get().ordinal();
      counts_A.addAndGet(i, 1);
      avgPhread_A.addAndGet(i, p);
      avgMapQ_A.addAndGet(i, mapQ);
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
