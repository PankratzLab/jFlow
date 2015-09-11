package seq.manage;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.Callable;

import seq.qc.FilterNGS;
import seq.qc.FilterNGS.SAM_FILTER_TYPE;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import common.WorkerTrain.Producer;
import filesys.Segment;

/**
 * @author lane0212 New version of the pileup, geared toward segments
 */
public class BamSegPileUp implements Iterator<BamPile> {

	private String bam;
	private ReferenceGenome referenceGenome;
	private Segment[] intervals;
	private QueryInterval[] queryIntervals;
	private int numReturned;
	private ArrayList<BamPile> bamPilesToReturn;
	private BamPile[] bamPiles;
	private SamReader reader;
	private Logger log;
	private SAMRecordIterator sIterator;
	private boolean[] addingMask;
	private AggregateFilter filter;
	private  FilterNGS filterNGS;

	public BamSegPileUp(String bam, String referenceGenomeFasta, Segment[] intervals, FilterNGS filterNGS, Logger log) {
		super();
		this.bam = bam;
		this.numReturned = 0;
		this.referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
		this.intervals = intervals;
		this.reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
		this.log = log;
		this.queryIntervals = BamOps.convertSegsToQI(intervals, reader.getFileHeader(), 0, true, log);
		log.reportTimeInfo("Optimizing " + intervals.length + " queries for pile up");
		this.sIterator = reader.query(queryIntervals, false);
		this.addingMask = new boolean[intervals.length];
		this.bamPiles = new BamPile[intervals.length];
		this.bamPilesToReturn = new ArrayList<BamPile>();
		this.filterNGS =filterNGS;
		this.filter = FilterNGS.initializeFilters(filterNGS, SAM_FILTER_TYPE.COPY_NUMBER, log);
		for (int i = 0; i < intervals.length; i++) {
			bamPiles[i] = new BamPile(intervals[i]);
			addingMask[i] = true;
		}
	}

	@Override
	public boolean hasNext() {
		while (sIterator.hasNext() && bamPilesToReturn.size() == 0) {
			SAMRecord samRecord = sIterator.next();
			if (!filter.filterOut(samRecord)) {
				Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
				String[] refSeq = referenceGenome.getSequenceFor(samRecordSegment);
				for (int i = 0; i < addingMask.length; i++) {

					if (addingMask[i]) {

						boolean overlaps = samRecordSegment.overlaps(bamPiles[i].getBin());

						if (overlaps) {
							bamPiles[i].addRecord(samRecord, refSeq, filterNGS.getPhreadScoreFilter(), log);
						} else {

							if (i == numReturned && samRecordSegment.getStart() > bamPiles[i].getBin().getStop() || !sIterator.hasNext()) {

								addingMask[i] = false;
								bamPilesToReturn.add(bamPiles[i]);
							}
						}
					}
				}

			}

		}
		if (!sIterator.hasNext()) {

			for (int i = 0; i < addingMask.length; i++) {
				if (addingMask[i]) {
					addingMask[i] = false;
					bamPilesToReturn.add(bamPiles[i]);
				}
			}
		}
		boolean hasNext = bamPilesToReturn.size() > 0;
		if (!hasNext) {
			boolean hasError = false;
			String error = "";
			if (Array.booleanArraySum(addingMask) != 0) {
				hasError = true;
				error += "Not all segments have been accounted for while searching " + bam;
			}
			if (numReturned != intervals.length || numReturned != bamPiles.length) {
				error += "\nNext not found, but not all intervals have been returned in  " + bam;
			}

			if (hasError) {
				log.reportTimeError(error);
				throw new IllegalStateException(error);
			}
		} else {
			numReturned++;
		}
		if (numReturned % 100 == 0) {
			log.reportTimeInfo(numReturned + " queries found for " + bam);
			log.memoryPercentFree();
		}
		return hasNext;
	}

	@Override
	public BamPile next() {
		// TODO Auto-generated method stub
		return bamPilesToReturn.remove(0);
	}

	public static class PileUpWorker implements Callable<BamPile[]> {
		private String bamFile;
		private String referenceGenomeFasta;
		private String serDir;
		private Logger log;
		private Segment[] pileSegs;
		private FilterNGS filterNGS;

		public PileUpWorker(String bamFile, String serDir, String referenceGenomeFasta, Segment[] pileSegs, FilterNGS filterNGS, Logger log) {
			super();
			this.bamFile = bamFile;
			this.serDir = serDir;
			this.referenceGenomeFasta = referenceGenomeFasta;
			this.pileSegs = pileSegs;
			this.filterNGS =filterNGS;
			this.log = log;

		}

		@Override
		public BamPile[] call() throws Exception {
			String ser = serDir + ext.rootOf(bamFile) + ".ser";
			if (!Files.exists(ser)) {
				BamSegPileUp bamSegPileUp = new BamSegPileUp(bamFile, referenceGenomeFasta, pileSegs, filterNGS, log);
				ArrayList<BamPile> bamPiles = new ArrayList<BamPile>();
				while (bamSegPileUp.hasNext()) {
					BamPile bamPile = bamSegPileUp.next();
					bamPile.summarize();
					bamPiles.add(bamPile);
				}
				BamPile[] bamPilesFinal = bamPiles.toArray(new BamPile[bamPiles.size()]);
				BamPile.writeSerial(bamPilesFinal, ser);
				return bamPilesFinal;
			} else {
				return BamPile.readSerial(ser, log);
			}
		}
	}

	public static class PileupProducer implements Producer<BamPile[]> {
		private int index;
		private String[] bamFiles;
		private String referenceGenomeFasta;
		private String serDir;
		private Logger log;
		private Segment[] pileSegs;
		private FilterNGS filterNGS;

		public PileupProducer(String[] bamFiles, String serDir, String referenceGenomeFasta, FilterNGS filterNGS, Segment[] pileSegs, Logger log) {
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
		public Callable<BamPile[]> next() {
			PileUpWorker worker = new PileUpWorker(bamFiles[index], serDir, referenceGenomeFasta, pileSegs, filterNGS, log);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

}
