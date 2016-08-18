package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.Callable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;

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
  private int numReturned;
  private final BamPile[] bamPiles;
  private final SamReader reader;
  private final Logger log;
  private final AggregateFilter filter;
  private final FilterNGS filterNGS;
  private int queryIndex;
  private final ReferenceGenome referenceGenome;

  public BamSegPileUp(String bam, String referenceGenomeFasta, Segment[] intervals,
                      FilterNGS filterNGS, Logger log) {
    super();
    this.bam = bam;
    numReturned = 0;
    reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
    this.log = log;
    referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
    bamPiles = new BamPile[intervals.length];
    this.filterNGS = filterNGS;
    filter = FilterNGS.initializeFilters(filterNGS, SAM_FILTER_TYPE.COPY_NUMBER, log);
    for (int i = 0; i < intervals.length; i++) {
      bamPiles[i] = new BamPile(intervals[i]);
    }
    queryIndex = 0;
  }

  @Override
  public boolean hasNext() {
    return queryIndex < bamPiles.length;
  }

  @Override
  public BamPile next() {
    BamPile currentPile = bamPiles[queryIndex];
    String[] currentRef = null;
    if (referenceGenome != null) {
      currentRef = referenceGenome.getSequenceFor(currentPile.getBin());
    }

    Segment cs = currentPile.getBin();
    CloseableIterator<SAMRecord> iterator =
        reader.queryOverlapping(Positions.getChromosomeUCSC(cs.getChr(), true), cs.getStart(),
                                cs.getStop());
    while (iterator.hasNext()) {
      SAMRecord samRecord = iterator.next();
      if (!filter.filterOut(samRecord)) {
        Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
        boolean overlaps = samRecordSegment.overlaps(cs);
        if (!overlaps) {
          String error = "non overlapping record returned for query";
          log.reportTimeError(error);
          throw new IllegalStateException(error);
        } else {
          currentPile.addRecord(samRecord, currentRef, filterNGS.getPhreadScoreFilter(), log);
        }
      }
    }
    queryIndex++;
    numReturned++;
    iterator.close();
    if (numReturned % 1000 == 0) {
      log.reportTimeInfo(numReturned + " queries found for " + bam);
    }
    return currentPile;
  }

  public static class BamPileResult {
    private final String ser;
    private final String bam;

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

    public BamPile[] loadResults(Logger log) {
      return BamPile.readSerial(ser, log);
    }

  }
  public static class PileUpWorker implements Callable<BamPileResult> {
    private final String bamFile;
    private final String serDir;
    private final Logger log;
    private final Segment[] pileSegs;
    private final FilterNGS filterNGS;
    private final String referenceGenomeFasta;

    public PileUpWorker(String bamFile, String serDir, String referenceGenomeFasta,
                        Segment[] pileSegs, FilterNGS filterNGS, Logger log) {
      super();
      this.bamFile = bamFile;
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
        BamSegPileUp bamSegPileUp =
            new BamSegPileUp(bamFile, referenceGenomeFasta, pileSegs, filterNGS, log);
        ArrayList<BamPile> bamPiles = new ArrayList<BamPile>();
        while (bamSegPileUp.hasNext()) {
          BamPile bamPile = bamSegPileUp.next();
          bamPile.summarize();
          bamPiles.add(bamPile);
        }
        BamPile[] bamPilesFinal = bamPiles.toArray(new BamPile[bamPiles.size()]);
        BamPile.writeSerial(bamPilesFinal, ser);
      }
      return new BamPileResult(bamFile, ser);
    }
  }

  public static class PileupProducer extends AbstractProducer<BamPileResult> {
    private int index;
    private final String[] bamFiles;
    private final String serDir;
    private final Logger log;
    private final Segment[] pileSegs;
    private final FilterNGS filterNGS;
    private final String referenceGenomeFasta;

    public PileupProducer(String[] bamFiles, String serDir, String referenceGenomeFasta,
                          FilterNGS filterNGS, Segment[] pileSegs, Logger log) {
      super();
      this.bamFiles = bamFiles;
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
      PileUpWorker worker =
          new PileUpWorker(bamFiles[index], serDir, referenceGenomeFasta, pileSegs, filterNGS, log);
      index++;
      return worker;
    }

    @Override
    public void shutdown() {
      // TODO Auto-generated method stub

    }
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }

}
