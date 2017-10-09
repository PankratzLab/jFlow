package org.genvisis.seq.manage;

import java.io.File;
import java.lang.ref.SoftReference;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
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
import com.google.common.collect.Maps;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeRangeMap;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.util.CloseableIterator;

/**
 * @author lane0212 New version of the pileup, geared toward segments
 */
public class BamSegPileUp implements Iterator<BamPile> {

  private final String bam;
  private final Map<Byte, RangeMap<Integer, BamPile>> bamPileMap;
  private final BamPile[] bamPileArray;
  private final Set<Segment> segsPiled;
  private final Map<Segment, SoftReference<String[]>> refSequences;
  private final CloseableIterator<SAMRecord> samRecordIterator;
  private final Logger log;
  private final AggregateFilter filter;
  private final ASSEMBLY_NAME aName;
  private final FilterNGS filterNGS;
  private final ReferenceGenome referenceGenome;

  /**
   * @param bam the bam file to pile
   * @param referenceGenomeFasta corresponding reference genome
   * @param intervals the intervals to pile on
   * @param filterNGS any filters to apply using {@link FilterNGS}
   * @param aName corresponding {@link ASSEMBLY_NAME} to determine proper contig searching
   * @param log
   */
  public BamSegPileUp(String bam, String referenceGenomeFasta, Segment[] intervals,
                      FilterNGS filterNGS, ASSEMBLY_NAME aName, Logger log) {
    super();
    this.bam = bam;
    this.aName = aName;
    SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
    samRecordIterator = reader.queryOverlapping(BamOps.convertSegsToQI(intervals,
                                                                       reader.getFileHeader(), 0,
                                                                       false, aName.addChr(), log));
    this.log = log;
    referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
    bamPileMap = Maps.newHashMap();
    bamPileArray = new BamPile[intervals.length];
    for (int i = 0; i < intervals.length; i++) {
      Segment interval = intervals[i];
      byte chr = interval.getChr();
      int start = interval.getStart();
      int stop = interval.getStop();
      RangeMap<Integer, BamPile> chrRangeMap = bamPileMap.get(chr);
      if (chrRangeMap == null) {
        chrRangeMap = TreeRangeMap.create();
        bamPileMap.put(chr, chrRangeMap);
      }
      BamPile bamPile = new BamPile(interval);
      chrRangeMap.put(Range.closed(start, stop), bamPile);
      bamPileArray[i] = bamPile;
    }
    this.segsPiled = Sets.newHashSetWithExpectedSize(intervals.length);
    this.refSequences = Maps.newHashMap();
    this.filterNGS = filterNGS;
    filter = FilterNGS.initializeFilters(filterNGS, SAM_FILTER_TYPE.COPY_NUMBER, log);
  }

  @Override
  public boolean hasNext() {
    return samRecordIterator.hasNext();
  }

  @Override
  public BamPile next() {
    SAMRecord samRecord = samRecordIterator.next();
    if (!filter.filterOut(samRecord)) {
      Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
      if (samRecordSegment.getStart() == 0 || samRecordSegment.getStop() == 0) return null;
      RangeMap<Integer, BamPile> chrRangeMap = bamPileMap.get(samRecordSegment.getChr());
      if (chrRangeMap == null) {
        log.reportTimeWarning("Read with chromosome " + samRecordSegment.getChr()
                              + " returned from BAM query but not found in requested intervals");
        chrRangeMap = ImmutableRangeMap.of();
        bamPileMap.put(samRecordSegment.getChr(), chrRangeMap);
      }
      Collection<BamPile> overlappingPiles = chrRangeMap.subRangeMap(Range.closed(samRecordSegment.getStart(),
                                                                                  samRecordSegment.getStop()))
                                                        .asMapOfRanges().values();

      for (BamPile bamPile : overlappingPiles) {
        Segment cs = bamPile.getBin();
        boolean overlaps = samRecordSegment.overlaps(cs);
        if (!overlaps) {
          String error = "non overlapping record returned for query";
          log.reportError(error);
          throw new IllegalStateException(error);
        } else {
          bamPile.addRecord(samRecord, getReferenceSequence(cs), filterNGS.getPhreadScoreFilter(),
                            log);
        }
        if (segsPiled.add(cs) && segsPiled.size() % 1000 == 0) {
          log.reportTimeInfo(segsPiled.size() + " bam piles added to of " + bamPileArray.length
                             + " for " + bam);
          log.memoryUsed();
          log.memoryTotal();
          log.memoryMax();
        }
      }
    }
    return null;
  }

  public BamPile[] summarizeAndFinalize() {
    for (BamPile bamPile : bamPileArray) {
      bamPile.summarize();
    }
    return bamPileArray;
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
  /**
   * @author Kitty Callable for threading pileups across bam files
   */
  public static class PileUpWorker implements Callable<BamPileResult> {

    private final String bamFile;
    private final ASSEMBLY_NAME aName;
    private final String serDir;
    private final Logger log;
    private final Segment[] pileSegs;
    private final FilterNGS filterNGS;
    private final String referenceGenomeFasta;

    public PileUpWorker(String bamFile, String serDir, String referenceGenomeFasta,
                        Segment[] pileSegs, FilterNGS filterNGS, ASSEMBLY_NAME aName, Logger log) {
      super();
      this.bamFile = bamFile;
      this.aName = aName;
      this.serDir = serDir;
      this.referenceGenomeFasta = referenceGenomeFasta;
      this.pileSegs = pileSegs;
      this.filterNGS = filterNGS;
      this.log = log;

    }

    @Override
    public BamPileResult call() throws Exception {
      String ser = serDir + ext.rootOf(bamFile) + ".ser";
      if (!Files.exists(ser)) {
        BamOps.verifyIndex(bamFile, log);
        BamSegPileUp bamSegPileUp = new BamSegPileUp(bamFile, referenceGenomeFasta, pileSegs,
                                                     filterNGS, aName, log);
        while (bamSegPileUp.hasNext()) {
          bamSegPileUp.next();
        }
        BamPile[] bamPilesFinal = bamSegPileUp.summarizeAndFinalize();
        int numMiss = 0;
        for (BamPile element : bamPilesFinal) {
          numMiss += element.getNumBasesWithMismatch();
        }
        log.reportTimeInfo(bamFile + " had " + numMiss + " mismatched bases");
        BamPile.writeSerial(bamPilesFinal, ser);
      }
      return new BamPileResult(bamFile, ser);
    }
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
    }

    @Override
    public boolean hasNext() {
      return index < bamFiles.length;
    }

    @Override
    public Callable<BamPileResult> next() {
      PileUpWorker worker = new PileUpWorker(bamFiles[index], serDir, referenceGenomeFasta,
                                             pileSegs, filterNGS, aName, log);
      index++;
      return worker;
    }

    @Override
    public void shutdown() {

    }
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }

}
