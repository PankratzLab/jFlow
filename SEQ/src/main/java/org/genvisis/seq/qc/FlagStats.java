package org.genvisis.seq.qc;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;

/**
 * Store read summary stats akin to samtools flagstat
 */
public class FlagStats {

  private int total;
  private int qcPassed;
  private int secondary;
  private int supplementary;
  private int duplicates;
  private int mapped;
  private int paired;
  private int read1;
  private int read2;
  private int properlyPaired;
  private int itselfAndMateMapped;
  private int singleton;
  private int mateMappedToDiffChr;
  private int mateMappedToDiffChrMQgte5;

  public FlagStats() {
    // no need
  }

  public int getTotal() {
    return total;
  }

  public int getQcPassed() {
    return qcPassed;
  }

  public int getSecondary() {
    return secondary;
  }

  public int getSupplementary() {
    return supplementary;
  }

  public int getDuplicates() {
    return duplicates;
  }

  public int getMapped() {
    return mapped;
  }

  public int getPaired() {
    return paired;
  }

  public int getRead1() {
    return read1;
  }

  public int getRead2() {
    return read2;
  }

  public int getProperlyPaired() {
    return properlyPaired;
  }

  public int getItselfAndMateMapped() {
    return itselfAndMateMapped;
  }

  public int getSingleton() {
    return singleton;
  }

  public int getMateMappedToDiffChr() {
    return mateMappedToDiffChr;
  }

  public int getMateMappedToDiffChrMQgte5() {
    return mateMappedToDiffChrMQgte5;
  }

  public static List<String> getHeader() {
    List<String> header = new ArrayList<>();
    header.add("TOTAL");
    header.add("QC_PASSED");
    header.add("SECONDARY");
    header.add("SUPPLEMENTARY");
    header.add("DUPLICATES");
    header.add("MAPPED");
    header.add("PAIRED");
    header.add("READ1");
    header.add("READ2");
    header.add("PROPERLY_PAIRED");
    header.add("ITSELF_AND_MATE_MAPPED");
    header.add("SINGLETON");
    header.add("MATE_MAPPED_TO_DIFFERENT_CHR");
    header.add("MATE_MAPPED_TO_DIFFERENT_CHR_MQ_GTE_5");

    return header;

  }

  public List<Integer> getData() {
    List<Integer> data = new ArrayList<>();
    data.add(total);
    data.add(qcPassed);
    data.add(secondary);
    data.add(supplementary);
    data.add(duplicates);
    data.add(mapped);
    data.add(paired);
    data.add(read1);
    data.add(read2);
    data.add(properlyPaired);
    data.add(itselfAndMateMapped);
    data.add(singleton);
    data.add(mateMappedToDiffChr);
    data.add(mateMappedToDiffChrMQgte5);
    return data;

  }

  /**
   * @param toCombine Flagstats to combine
   * @return the combination of all {@link FlagStats}
   */
  public static FlagStats combine(List<FlagStats> toCombine) {
    FlagStats combo = new FlagStats();
    for (FlagStats flagStats : toCombine) {
      combo.total += flagStats.total;
      combo.qcPassed += flagStats.qcPassed;
      combo.secondary += flagStats.secondary;
      combo.supplementary += flagStats.supplementary;
      combo.duplicates += flagStats.duplicates;
      combo.mapped += flagStats.mapped;
      combo.paired += flagStats.paired;
      combo.read1 += flagStats.read1;
      combo.read2 += flagStats.read2;
      combo.properlyPaired += flagStats.properlyPaired;
      combo.itselfAndMateMapped += flagStats.itselfAndMateMapped;
      combo.singleton += flagStats.singleton;
      combo.mateMappedToDiffChr += flagStats.mateMappedToDiffChr;
      combo.mateMappedToDiffChrMQgte5 += flagStats.mateMappedToDiffChrMQgte5;

    }
    return combo;
  }

  /**
   * @param samRecord compute stats for this record
   */
  public void stat(SAMRecord samRecord) {
    total++;
    qcPassed += !samRecord.getReadFailsVendorQualityCheckFlag() ? 1 : 0;
    secondary += samRecord.getNotPrimaryAlignmentFlag() ? 1 : 0;
    supplementary += samRecord.getSupplementaryAlignmentFlag() ? 1 : 0;
    duplicates += samRecord.getDuplicateReadFlag() ? 1 : 0;
    mapped += !samRecord.getReadUnmappedFlag() ? 1 : 0;
    paired += samRecord.getReadPairedFlag() ? 1 : 0;
    read1 += samRecord.getFirstOfPairFlag() ? 1 : 0;
    read2 += !samRecord.getFirstOfPairFlag() ? 1 : 0;
    properlyPaired += samRecord.getProperPairFlag() ? 1 : 0;
    itselfAndMateMapped += (!samRecord.getReadUnmappedFlag()
                            && !samRecord.getMateUnmappedFlag()) ? 1 : 0;
    singleton += (!samRecord.getReadUnmappedFlag() && samRecord.getMateUnmappedFlag()) ? 1 : 0;
    mateMappedToDiffChr += (samRecord.getMateReferenceIndex() != samRecord.getReferenceIndex()
                            && !samRecord.getMateUnmappedFlag() && !samRecord.getReadUnmappedFlag())
                                                                                                     ? 1
                                                                                                     : 0;
    mateMappedToDiffChrMQgte5 += (samRecord.getMateReferenceIndex() != samRecord.getReferenceIndex()
                                  && !samRecord.getMateUnmappedFlag()
                                  && !samRecord.getReadUnmappedFlag()
                                  && samRecord.getMappingQuality() >= 5
                                     & samRecord.getMappingQuality() != 255) ? 1 : 0;
  }

}
