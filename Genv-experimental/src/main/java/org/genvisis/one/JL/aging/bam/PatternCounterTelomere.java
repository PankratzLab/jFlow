package org.genvisis.one.JL.aging.bam;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.one.JL.aging.bam.BamAnalysisFilters.GCRangeFilter;
import org.genvisis.one.JL.aging.bam.BamAnalysisFilters.PatternFilter;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.SeqOps.GC_COMP_METHOD;
import org.genvisis.seq.qc.FlagStats;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * primarily methods to mimic Telseq: A issues noticed in TelSeq (that are repeated here): 1. GC
 * content calculation includes Ns -
 * https://github.com/zd1/telseq/blob/35019d2291c89162a34f322680a23d1a46a87990/
 * src/Telseq/telseq.cpp#L516-L527 2. Telomeric reads are included as normalization reads 3. Total
 * reads is off by one -
 * https://github.com/zd1/telseq/blob/35019d2291c89162a34f322680a23d1a46a87990/
 * src/Telseq/telseq.cpp#L214-L215 4. No "wiggle" for motif counting 5.
 */
public class PatternCounterTelomere implements BamAnalysis {

  // https://github.com/zd1/telseq/blob/35019d2291c89162a34f322680a23d1a46a87990/src/Telseq/telseq.h#L38
  private static final double GC_TELOMERIC_LOWERBOUND = .52;
  private static final double GC_TELOMERIC_UPPERBOUND = .48;

  private static final double GENOME_LENGTH_AT_TEL_GC = 332720800;
  private static final double LENGTH_UNIT = 1000;
  private static final double TELOMERE_ENDS = 46;

  private GCRangeFilter gcTelomericRange;
  private FlagStats telomericFlagStats;
  private FlagStats normalizationFlagStats;
  private FlagStats totalFlagStats;
  private PatternFilter pattern;
  private List<SamRecordFilter> baseFilters;
  private List<SamRecordFilter> nonCountFilters;
  private String analysisName;

  /**
   * @param pattern a {@link PatternFilter} to count
   * @param gcPatternRange range of GC content to use for pattern matches
   * @param gcBinNonTelomericRanges ranges of GC content to use for non - pattern matches
   * @param exclusionSet reads overlapping any ({@link Segment#overlaps(Segment)}) will not be used
   * @param baseFilters All of these {@link SamRecordFilter}s will be applied to every read
   * @param nonCountFilters will exclusively be applied to reads that pass the baseFilter, but are
   *          not flagged by the {@link PatternFilter}
   */
  public PatternCounterTelomere(PatternFilter pattern, GCRangeFilter gcPatternRange,
                                List<SamRecordFilter> baseFilters,
                                List<SamRecordFilter> nonCountFilters) {
    super();
    this.pattern = pattern;
    this.gcTelomericRange = gcPatternRange;
    this.baseFilters = baseFilters;
    this.nonCountFilters = nonCountFilters;
    this.telomericFlagStats = new FlagStats();
    this.normalizationFlagStats = new FlagStats();
    this.totalFlagStats = new FlagStats();
  }

  // SRR1738845 test case

  @Override
  public void analyze(SAMRecord samRecord, Logger log) {
    totalFlagStats.stat(samRecord);

    if (!pattern.filterOut(samRecord)) {
      telomericFlagStats.stat(samRecord);
    }
    // Telseq uses telomeric reads for normalization :(
    boolean use = true;
    for (SamRecordFilter filter : nonCountFilters) {
      if (use) {
        use = !filter.filterOut(samRecord);
      }
    }
    if (use) {
      // https://github.com/zd1/telseq/blob/35019d2291c89162a34f322680a23d1a46a87990/src/Telseq/telseq.cpp#L288-L295
      if (!gcTelomericRange.filterOut(samRecord)) {
        normalizationFlagStats.stat(samRecord);
      }
    }

  }

  @Override
  public void summarize() {
    // sf

  }

  @Override
  public List<String> getOutput() {

    List<String> output = new ArrayList<>();
    output.add("SOURCE\t" + ArrayUtils.toStr(FlagStats.getHeader()) + "\tLENGTH_ESTIMATE");
    output.add("TELOMERIC\t" + ArrayUtils.toStr(telomericFlagStats.getData()) + "\t"
               + computeTelomereLength(telomericFlagStats, normalizationFlagStats));
    output.add("NORMALIZATION\t" + ArrayUtils.toStr(normalizationFlagStats.getData()) + "\tNaN");
    output.add("TOTAL\t" + ArrayUtils.toStr(totalFlagStats.getData()) + "\tNaN");
    return output;

  }

  @Override
  public boolean shouldAnalyze(SAMRecord samRecord, Logger log) {

    for (SamRecordFilter filter : baseFilters) {
      if (filter.filterOut(samRecord)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public void init(String bamFile, Logger log) {
    this.analysisName = BamOps.getSampleName(bamFile, log) + "tel.results.txt";
  }

  /**
   * @param minCount minimum number of times for a pattern to be seen
   * @return
   */
  public static PatternFilter getTelomericPattern(int minCount) {
    // https://github.com/zd1/telseq/blob/35019d2291c89162a34f322680a23d1a46a87990/src/Telseq/telseq.h#L33-L34
    HashSet<String> telPattern = new HashSet<>();
    telPattern.add("TTAGGG");
    telPattern.add("CCCTAA");// rev comp
    return new PatternFilter(telPattern, minCount);
  }

  @Override
  public String getRootOutputFile() {
    return analysisName;
  }

  /**
   * @param telomericFlagStats {@link FlagStats} generated for telomeric reads
   * @param normalizationFlagStats {@link FlagStats} generated for normalization reads
   * @param genomeLengthAtTelGC length of the genome at telomeric GC
   * @param lengthUnit bp length unit
   * @param telomereEnds number of telomere ends, 46 in humans
   * @return
   */
  private static double computeTelomereLength(FlagStats telomericFlagStats,
                                              FlagStats normalizationFlagStats,
                                              double genomeLengthAtTelGC, double lengthUnit,
                                              double telomereEnds) {
    double ratio = (double) telomericFlagStats.getTotal() / normalizationFlagStats.getTotal();
    double scale = genomeLengthAtTelGC / lengthUnit / telomereEnds;
    return ratio * scale;
  }

  private static double computeTelomereLength(FlagStats telomericFlagStats,
                                              FlagStats normalizationFlagStats) {
    return computeTelomereLength(telomericFlagStats, normalizationFlagStats,
                                 GENOME_LENGTH_AT_TEL_GC, LENGTH_UNIT, TELOMERE_ENDS);
  }

  public static GCRangeFilter getDefaultGCTelRange() {
    return new GCRangeFilter(GC_TELOMERIC_UPPERBOUND, GC_TELOMERIC_LOWERBOUND,
                             GC_COMP_METHOD.N_COUNTS_FOR_TOTAL);
  }

}
