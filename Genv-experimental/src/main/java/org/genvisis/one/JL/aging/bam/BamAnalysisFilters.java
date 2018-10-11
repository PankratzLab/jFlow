package org.genvisis.one.JL.aging.bam;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.genvisis.seq.manage.SamRecordOps;
import org.genvisis.seq.manage.SeqOps;
import org.genvisis.seq.manage.SeqOps.GC_COMP_METHOD;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Segment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Basically to store some appropriate?? filters
 */
public class BamAnalysisFilters {

  private BamAnalysisFilters() {

  }

  /**
   * @param minMapQ minimum mapping quality to use
   * @param log
   * @return a list of filters that will capture appropriate mapping and higher quality reads for
   *         normalization
   */
  public static List<SamRecordFilter> getDefualtNormReadFilters(double minMapQ, Logger log) {
    List<SamRecordFilter> filters = new FilterNGS().getStandardSAMRecordFilters(SAM_FILTER_TYPE.GENOTYPE,
                                                                                log);
    filters.add(new FilterNGS().getSamRecordMapQFilter(minMapQ));
    return filters;
  }

  /**
   * sequence content based filter
   */
  public static class PatternFilter implements SamRecordFilter {

    private Set<String> patterns;
    private int minCount;

    /**
     * @param patterns patterns to count
     * @param minCount minimum times pattern is seen (i.e. TelSeq requires 7 telomeric repeats).
     */
    public PatternFilter(Set<String> patterns, int minCount) {
      super();
      this.patterns = patterns;
      this.minCount = minCount;
    }

    /**
     * @param record
     * @return a map from pattern to number of occurrences
     */
    public Map<String, Integer> getAllCounts(SAMRecord record) {
      String read = record.getReadString();
      HashMap<String, Integer> counts = new HashMap<>();
      for (String pattern : patterns) {
        int count = SeqOps.countMotif(read, pattern, false);
        counts.put(pattern, count);
      }
      return counts;
    }

    /**
     * @param record {@link SAMRecord} to scan
     * @return the maximum count of all patterns
     */
    public int getMaxCount(SAMRecord record) {
      int maxPatternCount = 0;
      Map<String, Integer> counts = getAllCounts(record);
      for (String pattern : patterns) {
        int count = counts.get(pattern);
        if (count > maxPatternCount) {
          maxPatternCount = count;
        }
      }
      return maxPatternCount;
    }

    @Override
    public boolean filterOut(SAMRecord record) {
      return getMaxCount(record) < minCount;
    }

    @Override
    public boolean filterOut(SAMRecord first, SAMRecord second) {
      throw new UnsupportedOperationException();
    }

  }

  /**
   * Class to store min and max GC content filter for {@link SAMRecord}
   */
  public static class GCRangeFilter implements SamRecordFilter {

    private double minGC;
    private double maxGC;
    private GC_COMP_METHOD gMethod;

    /**
     * @param minGC minimum GC content threshold (proportion)
     * @param maxGC maximum GC content threshold (proportion)
     * @param gMethod {@link GC_COMP_METHOD}
     */
    public GCRangeFilter(double minGC, double maxGC, GC_COMP_METHOD gMethod) {
      super();
      this.minGC = minGC;
      this.maxGC = maxGC;
      this.gMethod = gMethod;
    }

    public double getMinGC() {
      return minGC;
    }

    public double getMaxGC() {
      return maxGC;
    }

    /**
     * @return
     */
    public GC_COMP_METHOD getgMethod() {
      return gMethod;
    }

    /**
     * @param filter
     * @return whether the other filter is contained within this range
     */
    public boolean passesOther(GCRangeFilter filter) {
      if (filter.minGC >= minGC && filter.maxGC <= maxGC) {
        return true;
      }
      return false;
    }

    @Override
    public boolean filterOut(SAMRecord record) {
      double gc = SeqOps.getProportionGC(SamRecordOps.getReadBases(record), gMethod);
      if (gc >= minGC && gc <= maxGC) {
        return false;
      }
      return true;
    }

    @Override
    public boolean filterOut(SAMRecord first, SAMRecord second) {
      throw new UnsupportedOperationException();
    }

  }

  /**
   * Filters out reads in a {@link LocusSet}
   *
   * @param <T> extends {@link Segment}
   */
  public static class RegionFilter<T extends Segment> implements SamRecordFilter {

    private LocusSet<T> exclusionSet;
    private Logger log;

    /**
     * @param exclusionSet what to exclude
     * @param log
     */
    public RegionFilter(LocusSet<T> exclusionSet, Logger log) {
      super();
      this.exclusionSet = exclusionSet;
      this.log = log;
    }

    @Override
    public boolean filterOut(SAMRecord record) {
      int[] olaps = exclusionSet.getOverlappingIndices(SamRecordOps.getReferenceSegmentForRecord(record,
                                                                                                 log));
      return olaps != null && olaps.length > 0;
    }

    @Override
    public boolean filterOut(SAMRecord first, SAMRecord second) {
      return false;
    }

  }

}
