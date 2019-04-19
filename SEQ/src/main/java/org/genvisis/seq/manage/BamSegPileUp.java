package org.genvisis.seq.manage;

import java.io.IOException;
import java.lang.ref.SoftReference;
import java.util.Collection;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;

/**
 * @author lane0212 New version of the pileup, geared toward segments
 */
public class BamSegPileUp {

  private final String bam;
  private final Map<Byte, RangeMap<Integer, Set<BamPile>>> bamPileMap;
  private final BamPile[] bamPileArray;
  private final Set<Segment> segsPiled;
  private final Map<Segment, SoftReference<String[]>> refSequences;
  private final Logger log;
  private final AggregateFilter filter;
  private final ASSEMBLY_NAME aName;
  private final FilterNGS filterNGS;
  private final ReferenceGenome referenceGenome;
  private final QueryInterval[] queryIntervals;
  private final int numThreads;

  /**
   * @param bam the bam file to pile
   * @param referenceGenomeFasta corresponding reference genome
   * @param intervals the intervals to pile on
   * @param filterNGS any filters to apply using {@link FilterNGS}
   * @param aName corresponding {@link ASSEMBLY_NAME} to determine proper contig searching
   * @param log
   */
  public BamSegPileUp(String bam, ReferenceGenome referenceGenome, Segment[] intervals,
                      QueryInterval[] queryIntervals, FilterNGS filterNGS, ASSEMBLY_NAME aName,
                      Logger log) {
    this(bam, referenceGenome, intervals, queryIntervals, filterNGS, aName,
         Runtime.getRuntime().availableProcessors(), log);
  }

  public BamSegPileUp(String bam, ReferenceGenome referenceGenome, Segment[] intervals,
                      QueryInterval[] queryIntervals, FilterNGS filterNGS, ASSEMBLY_NAME aName,
                      int numThreads, Logger log) {
    super();
    this.bam = bam;
    this.aName = aName;
    this.log = log;
    this.referenceGenome = referenceGenome;
    bamPileMap = Maps.newConcurrentMap();
    bamPileArray = new BamPile[intervals.length];

    for (Segment i : intervals) {
      if (!bamPileMap.containsKey(i.getChr())) bamPileMap.put(i.getChr(), TreeRangeMap.create());
    }

    for (int i = 0; i < intervals.length; i++) {
      Segment interval = intervals[i];
      BamPile bamPile = new BamPile(interval);
      bamPileArray[i] = bamPile;
    }
    createPileMap();
    if (queryIntervals == null) {
      SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT,
                                                 Sets.immutableEnumSet(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES));
      this.queryIntervals = BamOps.convertSegsToQI(intervals, reader.getFileHeader(), 0, true,
                                                   aName.addChr(), log);
      try {
        reader.close();
      } catch (IOException e) {}
    } else {
      this.queryIntervals = queryIntervals;
    }
    this.segsPiled = Sets.newConcurrentHashSet();
    this.refSequences = Maps.newConcurrentMap();
    this.filterNGS = filterNGS;
    filter = FilterNGS.initializeFilters(filterNGS, SAM_FILTER_TYPE.COPY_NUMBER, log);
    this.numThreads = numThreads;
  }

  private void createPileMap() {
    for (BamPile p : bamPileArray) {
      byte chr = p.getChr();
      int start = p.getStart();
      int stop = p.getStop();
      RangeMap<Integer, Set<BamPile>> chrRangeMap = bamPileMap.get(chr);
      Range<Integer> curSegmentRange = Range.closed(start, stop);
      RangeMap<Integer, Set<BamPile>> existingSubRangeMap = chrRangeMap.subRangeMap(curSegmentRange);
      Map<Range<Integer>, Set<BamPile>> existingSubRangeMapOfRanges = existingSubRangeMap.asMapOfRanges();
      if (!existingSubRangeMapOfRanges.isEmpty()) {
        RangeMap<Integer, Set<BamPile>> updatedMappings = TreeRangeMap.create();
        for (Map.Entry<Range<Integer>, Set<BamPile>> existingEntry : existingSubRangeMapOfRanges.entrySet()) {
          ImmutableSet.Builder<BamPile> newBamPileSet = ImmutableSet.builder();
          newBamPileSet.addAll(existingEntry.getValue());
          newBamPileSet.add(p);
          updatedMappings.put(existingEntry.getKey(), newBamPileSet.build());
        }
        existingSubRangeMap.putAll(updatedMappings);
      }
      RangeSet<Integer> unmappedRange = TreeRangeSet.create();
      unmappedRange.add(curSegmentRange);
      unmappedRange.removeAll(existingSubRangeMapOfRanges.keySet());
      for (Range<Integer> range : unmappedRange.asRanges()) {
        chrRangeMap.put(range, ImmutableSet.of(p));
      }
    }
  }

  /**
   * Perform the pileup of bam reads
   * 
   * @return BamPile[] of BamPile for every requested interval
   * @throws IOException when thrown by {@link SamReader}
   */
  public BamPile[] pileup() throws IOException {
    log.reportTime("Processing " + bam + " with " + numThreads + " threads.");
    try (SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT,
                                                    Sets.immutableEnumSet(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES));
         SAMRecordIterator iter = reader.query(queryIntervals, false)) {

      long t = System.nanoTime();
      AtomicLong count = new AtomicLong(0);

      ExecutorService exec = Executors.newFixedThreadPool(numThreads);
      for (int i = 0; i < numThreads; i++) {
        exec.submit(new Runnable() {

          @Override
          public void run() {
            SAMRecord record = null;
            try {
              while ((record = iter.next()) != null) {
                boolean added = process(record);
                if (added) {
                  count.incrementAndGet();
                }
              }
            } catch (NoSuchElementException e) {
              // done!
            }
          }
        });
      }
      exec.shutdown();
      try {
        exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
      } catch (InterruptedException e) {
        log.reportTimeWarning("Possible problem: " + e.getMessage());
      }

      summarizeAndFinalize();
      log.reportTime("Finalized bam piles (" + count.get() + " reads total) for " + bam + " after "
                     + ext.getTimeElapsedNanos(t));
      return bamPileArray;
    }
  }

  private boolean process(SAMRecord samRecord) {
    if (!filter.filterOut(samRecord)) {
      Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
      if (samRecordSegment.getStart() == 0 || samRecordSegment.getStop() == 0) return false;
      RangeMap<Integer, Set<BamPile>> chrRangeMap = bamPileMap.get(samRecordSegment.getChr());
      if (chrRangeMap == null) {
        log.reportTimeWarning("Read with chromosome " + samRecordSegment.getChr()
                              + " returned from BAM query but not found in requested intervals - skipping read");
        return false;
        // chrRangeMap = ImmutableRangeMap.of();
        // bamPileMap.put(samRecordSegment.getChr(), chrRangeMap);
      }
      // AtomicLong processTime = new AtomicLong(0);
      // AtomicInteger count = new AtomicInteger(0);

      // long t1 = System.nanoTime();
      chrRangeMap.subRangeMap(Range.closed(samRecordSegment.getStart(), samRecordSegment.getStop()))
                 .asMapOfRanges().values().stream().flatMap(Collection::stream) //
                 .parallel() //
                 .forEach((bamPile) -> {
                   // count.incrementAndGet();
                   // long t11 = System.nanoTime();
                   addRecordToPile(bamPile, samRecordSegment, samRecord);
                   // long t21 = System.nanoTime();
                   // processTime.addAndGet(t21 - t11);
                 });
      // long t2 = System.nanoTime();
      // long elaps = t2 - t1;
      // long ave = processTime.get() / count.get();
      // log.reportTime("Took " + ext.formatTimeElapsed(elaps, TimeUnit.NANOSECONDS) + " to process
      // "
      // + count.get() + " with an average processing time of "
      // + ext.formatTimeElapsed(ave, TimeUnit.NANOSECONDS));
      return true;
    }
    return false;
  }

  private void addRecordToPile(BamPile bamPile, Segment samRecordSegment, SAMRecord samRecord) {
    Segment cs = bamPile.getBin();
    boolean overlaps = samRecordSegment.overlaps(cs);
    if (!overlaps) {
      String error = "non overlapping record returned for query";
      log.reportError(error);
      throw new IllegalStateException(error);
    } else {
      bamPile.addRecordAtomic(samRecord, getReferenceSequence(cs), filterNGS.getPhreadScoreFilter(),
                              log);
    }
  }

  private void summarizeAndFinalize() {
    for (BamPile bamPile : bamPileArray) {
      bamPile.summarize();
    }
  }

  private String[] getReferenceSequence(Segment interval) {
    SoftReference<String[]> refSeqRef = refSequences.get(interval);
    String[] refSeq = refSeqRef == null ? null : refSeqRef.get();
    if (refSeq == null) {
      refSeq = referenceGenome.getSequenceFor(interval, aName, false);
      refSequences.put(interval, new SoftReference<String[]>(refSeq));
    }
    return refSeq;
  }

  /**
   * @author Kitty The result of a pileup. Storing the input bam filename and the output serialized
   *         file name
   */
  public static class BamPileResult {

    private final String ser;
    private final String bam;

    /**
     * @param bam full path to bam
     * @param ser full path to serialized file
     */
    public BamPileResult(String bam, String ser) {
      super();
      this.bam = bam;
      this.ser = ser;
    }

    public String getBam() {
      return bam;
    }

    public String getSer() {
      return ser;
    }

    /**
     * @param log
     * @return the de-serialized array of {@link BamPile}
     */
    public BamPile[] loadResults(Logger log) {
      return BamPile.readSerial(ser, log);
    }

  }

}
