package org.genvisis.seq.manage;

import java.io.File;
import java.io.IOException;
import java.lang.ref.SoftReference;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;
import com.google.common.collect.ImmutableRangeMap;
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
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.util.CloseableIterator;

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

  /**
   * @param bam the bam file to pile
   * @param referenceGenomeFasta corresponding reference genome
   * @param intervals the intervals to pile on
   * @param filterNGS any filters to apply using {@link FilterNGS}
   * @param aName corresponding {@link ASSEMBLY_NAME} to determine proper contig searching
   * @param log
   */
  public BamSegPileUp(String bam, String referenceGenomeFasta, Segment[] intervals,
                      QueryInterval[] queryIntervals, FilterNGS filterNGS, ASSEMBLY_NAME aName,
                      Logger log) {
    super();
    this.bam = bam;
    this.aName = aName;
    this.log = log;
    referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
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
    this.queryIntervals = queryIntervals;
    this.segsPiled = Sets.newConcurrentHashSet();
    this.refSequences = Maps.newConcurrentMap();
    this.filterNGS = filterNGS;
    filter = FilterNGS.initializeFilters(filterNGS, SAM_FILTER_TYPE.COPY_NUMBER, log);
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
    log.report("Initializing bam reading for " + bam);
    try (SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT,
                                                    Sets.immutableEnumSet(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES));
         CloseableIterator<SAMRecord> samRecordIterator = reader.queryOverlapping(queryIntervals)) {
      log.report("Iterating through reads for " + bam);
      AtomicInteger readsProcessed = new AtomicInteger(0);
      samRecordIterator.stream().parallel().forEach((samRecord) -> {
        process(samRecordIterator.next());
        int read;
        if ((read = readsProcessed.incrementAndGet()) % 10000000 == 0) {
          log.report("Processed " + read + " reads for " + bam);
        }
      });
      log.report("Finalizing bam piles for " + bam);
      summarizeAndFinalize();
      return bamPileArray;
    }
  }

  private void process(SAMRecord samRecord) {
    if (!filter.filterOut(samRecord)) {
      Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
      if (samRecordSegment.getStart() == 0 || samRecordSegment.getStop() == 0) return;
      RangeMap<Integer, Set<BamPile>> chrRangeMap = bamPileMap.get(samRecordSegment.getChr());
      if (chrRangeMap == null) {
        log.reportTimeWarning("Read with chromosome " + samRecordSegment.getChr()
                              + " returned from BAM query but not found in requested intervals");
        chrRangeMap = ImmutableRangeMap.of();
        bamPileMap.put(samRecordSegment.getChr(), chrRangeMap);
      }
      chrRangeMap.subRangeMap(Range.closed(samRecordSegment.getStart(), samRecordSegment.getStop()))
                 .asMapOfRanges().values().stream().flatMap(Collection::stream).parallel()
                 .forEach((bamPile) -> {
                   addRecordToPile(bamPile, samRecordSegment, samRecord);
                 });
    }
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
    if (segsPiled.add(cs) && segsPiled.size() % 10000 == 0) {
      log.reportTimeInfo(segsPiled.size() + " bam piles added to of " + bamPileArray.length
                         + " for " + bam);
      log.memoryUsed();
      log.memoryTotal();
      log.memoryMax();
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

  public static BamPileResult processBamFile(String bamFile, String serDir,
                                             String referenceGenomeFasta, Segment[] pileSegs,
                                             QueryInterval[] qi, FilterNGS filterNGS,
                                             ASSEMBLY_NAME aName, Logger log) throws Exception {
    String ser = serDir + ext.rootOf(bamFile) + ".ser";
    if (!Files.exists(ser)) {
      BamOps.verifyIndex(bamFile, log);
      BamSegPileUp bamSegPileUp = new BamSegPileUp(bamFile, referenceGenomeFasta, pileSegs, qi,
                                                   filterNGS, aName, log);
      BamPile[] bamPilesFinal = bamSegPileUp.pileup();
      int numMiss = 0;
      for (BamPile element : bamPilesFinal) {
        numMiss += element.getNumBasesWithMismatch();
      }
      log.reportTimeInfo(bamFile + " had " + numMiss + " mismatched bases");
      BamPile.writeSerial(bamPilesFinal, ser);
    }
    return new BamPileResult(bamFile, ser);
  }

  public static class PileupProducer extends AbstractProducer<BamPileResult> {

    private int index;
    private final String[] bamFiles;
    private final ASSEMBLY_NAME aName;
    private final String serDir;
    private final Logger log;
    private final Segment[] pileSegs;
    private final FilterNGS filterNGS;
    private final String referenceGenomeFasta;
    private final Map<String, QueryInterval[]> qiMap;

    public PileupProducer(String[] bamFiles, String serDir, String referenceGenomeFasta,
                          FilterNGS filterNGS, Segment[] pileSegs, ASSEMBLY_NAME aName,
                          Logger log) {
      super();
      this.bamFiles = bamFiles;
      this.aName = aName;
      this.serDir = serDir;
      this.referenceGenomeFasta = referenceGenomeFasta;
      this.log = log;
      this.pileSegs = pileSegs;
      this.filterNGS = filterNGS;
      new File(serDir).mkdirs();
      qiMap = new HashMap<String, QueryInterval[]>();
      setupQIs();
    }

    private void setupQIs() {
      long t1 = System.nanoTime();
      Map<SAMSequenceDictionary, QueryInterval[]> dicts = new HashMap<>();
      for (String bam : bamFiles) {
        try (SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT,
                                                        Sets.immutableEnumSet(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES))) {
          SAMSequenceDictionary sd = reader.getFileHeader().getSequenceDictionary();
          QueryInterval[] qi = null;
          for (SAMSequenceDictionary d : dicts.keySet()) {
            boolean same = false;
            try {
              sd.assertSameDictionary(d);
              same = true;
            } catch (AssertionError e) {
              // not the same
            }
            if (same) {
              qi = dicts.get(d);
              break;
            }
          }
          if (qi == null) {
            qi = BamOps.convertSegsToQI(pileSegs, reader.getFileHeader(), 0, true, aName.addChr(),
                                        log);
            dicts.put(sd, qi);
          }
          qiMap.put(bam, qi);
        } catch (IOException e) {
          throw new IllegalStateException("Couldn't read header of BAM file: " + bam, e);
        }
      }
      dicts = null;
      log.reportTime("Constructed query intervals in " + ext.getTimeElapsedNanos(t1));
    }

    @Override
    public boolean hasNext() {
      return index < bamFiles.length;
    }

    @Override
    public Callable<BamPileResult> next() {
      final String bamFile = bamFiles[index];
      index++;
      return new Callable<BamSegPileUp.BamPileResult>() {

        @Override
        public BamPileResult call() throws Exception {
          return BamSegPileUp.processBamFile(bamFile, serDir, referenceGenomeFasta, pileSegs,
                                             qiMap.get(bamFile), filterNGS, aName, log);
        }
      };
    }

    @Override
    public void shutdown() {

    }
  }
}
