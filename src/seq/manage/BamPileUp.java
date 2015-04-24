package seq.manage;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.Callable;

import seq.qc.FilterNGS;
import seq.qc.FilterNGS.SAM_FILTER_TYPE;
import stats.Histogram.DynamicHistogram;
import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import filesys.Segment;

public class BamPileUp implements Iterator<BamPile> {
	public enum PILE_TYPE {
		/**
		 * Will only report positions with alternate alleles, and demands a reference genome
		 */
		CONTAMINATION, /**
		 * Pileup to all positions passing filters supplied
		 */
		REGULAR;
	}

	private String bam;
	private FilterNGS filterNGS;
	private ReferenceGenome referenceGenome;
	private Segment[] intervals;
	private int binSize;
	private PILE_TYPE pileType;
	private SAM_FILTER_TYPE filterType;
	private SAMRecordIterator sIterator;
	private AggregateFilter filter;
	private SamReader reader;
	private ArrayList<BamPile> bamPiles;
	private ArrayList<BamPile> bamPilesToReturn;
	private Segment currentSegment;
	private Logger log;
	private BamPileUpSummary bamPileUpSummary;
	private QueryInterval[] queryIntervals;
	private WorkerTrain<TmpBamPile> train;

	public BamPileUp(String bam, ReferenceGenome referenceGenome, int refBinSize, FilterNGS filterNGS, Segment[] intervals, PILE_TYPE type, SAM_FILTER_TYPE filterType, Logger log) {
		super();
		this.bam = bam;
		this.binSize = refBinSize;
		this.log = log;
		this.filterNGS = filterNGS;
		this.intervals = intervals;
		this.referenceGenome = referenceGenome;
		this.pileType = type;
		this.filterType = filterType;
		init();
	}

	private void init() {
		this.reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
		log.reportTimeInfo("Optimizing " + intervals.length + " queries for pile up");
		this.queryIntervals = BamOps.convertSegsToQI(intervals, reader.getFileHeader(), 0, true, log);
		log.reportTimeInfo("Finished Optimizing " + intervals.length + " queries to " + queryIntervals.length + " intervals for pile up");

		this.sIterator = queryIntervals == null ? reader.iterator() : reader.query(queryIntervals, false);
		this.filter = initializeFilters(filterNGS, filterType, log);
		this.bamPiles = new ArrayList<BamPile>();
		this.bamPilesToReturn = new ArrayList<BamPile>();
		this.currentSegment = new Segment((byte) 0, 0, 0);
		this.train = new WorkerTrain<BamPileUp.TmpBamPile>(null, 2, 200, log);// an extra thread should be about a half hour speed up per sample
		train.setAutoShutDown(false);
		this.bamPileUpSummary = new BamPileUpSummary(log);
	}

	public QueryInterval[] getQueryIntervals() {
		return queryIntervals;
	}

	private boolean shutdown() {
		boolean cleanShut = true;
		if (sIterator.hasNext()) {
			log.reportTimeWarning("The bam file " + bam + " has more recoreds and the shutdown method was called");
			cleanShut = false;
		}
		sIterator.close();
		try {
			reader.close();
		} catch (IOException e) {
			cleanShut = false;
			log.reportException(e);
			e.printStackTrace();
		}
		return cleanShut;
	}

	@Override
	public boolean hasNext() {
		while (sIterator.hasNext() && bamPilesToReturn.size() == 0) {
			SAMRecord samRecord = sIterator.next();
			bamPileUpSummary.setTotalReads(bamPileUpSummary.getTotalReads() + 1);
			if (!filter.filterOut(samRecord)) {
				bamPileUpSummary.setReadsPiled(bamPileUpSummary.getReadsPiled() + 1);
				if (bamPileUpSummary.getReadsPiled() % 100000 == 0) {
					log.reportTimeInfo("~" + bamPileUpSummary.getReadsPiled() + " of " + bamPileUpSummary.getTotalReads() + " total reads piled to " + bamPileUpSummary.getPositionsPiled() + " positions (" + ext.getTimeElapsed(bamPileUpSummary.getTime()) + ") " + SamRecordOps.getDisplayLoc(samRecord));
					bamPileUpSummary.setTime(System.currentTimeMillis());
				}
				BamPileInitializer bamPileInitializer = new BamPileInitializer(binSize, currentSegment, samRecord, log);
				while (bamPileInitializer.hasNext()) {
					bamPiles.add(bamPileInitializer.next());
				}
				currentSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
				TmpBamPileProducer tmpBamPileProducer = new TmpBamPileProducer(samRecord, currentSegment, bamPiles.toArray(new BamPile[bamPiles.size()]), filterNGS, log);
				train.setProducer(tmpBamPileProducer);
				bamPiles = new ArrayList<BamPile>(bamPiles.size());
				while (train.hasNext()) {
					TmpBamPile tmpBamPile = train.next();
					if (tmpBamPile.overlapsCurrentRecord()) {
						bamPiles.add(tmpBamPile.getBamPile());// stays in this round
					} else {
						BamPile bamPile = tmpBamPile.getBamPile();
						bamPile.setReference(referenceGenome);
						if (filterNGS.getReadDepthFilter() == null || bamPile.getTotalDepth(false, false) > filterNGS.getReadDepthFilter()[0]) {
							int altAlleleDepth = filterNGS.getReadDepthFilter().length > 1 ? filterNGS.getReadDepthFilter()[1] : -1;
							if (pileType == PILE_TYPE.REGULAR || (bamPile.hasAltAllele(log) && bamPile.hasOnlyOneAlt(log) && bamPile.getNumAlt(log) > altAlleleDepth && bamPile.getNumRef(log) > altAlleleDepth)) {
								bamPileUpSummary.setPositionsPiled(bamPileUpSummary.getPositionsPiled() + 1);
								bamPileUpSummary.addToHistogram(bamPile.getPropRef(log));
								bamPilesToReturn.add(bamPile);
							}
						}
					}
				}
			}
		}

		return bamPilesToReturn.size() > 0;
	}

	@Override
	public BamPile next() {
		return bamPilesToReturn.remove(0);
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub

	}

	public BamPileUpSummary getBamPileUpSummary() {
		return bamPileUpSummary;
	}

	private static AggregateFilter initializeFilters(FilterNGS filterNGS, SAM_FILTER_TYPE filterType, Logger log) {
		ArrayList<SamRecordFilter> filters = filterNGS.getStandardSAMRecordFilters(filterType, log);
		filters.add(filterNGS.getSamRecordMapQFilter(filterNGS.getMappingQualityFilter()));
		AggregateFilter filter = new AggregateFilter(filters);
		return filter;
	}

	/**
	 * Feeds up new bamPiles for the next base pairs determined by the next alignment
	 *
	 */
	private static class BamPileInitializer implements Iterator<BamPile> {
		private int binSize;
		private Segment previousBin;
		private Segment nextBin;
		private Logger log;
		private Segment samRecordSeg;
		private SAMRecord samRecord;

		public BamPileInitializer(int binSize, Segment previousBin, SAMRecord samRecord, Logger log) {
			super();
			this.binSize = binSize;
			this.previousBin = previousBin;
			this.samRecord = samRecord;
			this.samRecordSeg = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
			this.nextBin = samRecordSeg.overlaps(previousBin) ? getNextBin(previousBin, binSize) : scanToNext();
		}

		private Segment scanToNext() {
			if (samRecordSeg.getChr() > 0) {
				if (previousBin.getChr() != samRecordSeg.getChr()) {
					previousBin = new Segment(samRecordSeg.getChr(), 1, 1 + binSize);
				}
				while (!samRecordSeg.overlaps(previousBin) && previousBin.getStop() <= samRecordSeg.getStop()) {
					previousBin = new Segment(previousBin.getChr(), previousBin.getStop() + 1, previousBin.getStop() + binSize);
				}
			} else {
				log.reportTimeError("Could not find valid segment for " + samRecord.toString());
			}
			return previousBin;
		}

		private static Segment getNextBin(Segment seg, int binSize) {
			return new Segment(seg.getChr(), seg.getStop() + 1, seg.getStop() + binSize);
		}

		@Override
		public boolean hasNext() {
			return nextBin.overlaps(samRecordSeg) && nextBin.getStop() <= samRecordSeg.getStop();
		}

		@Override
		public BamPile next() {
			BamPile bamPile = new BamPile(nextBin);
			nextBin = getNextBin(nextBin, binSize);
			return bamPile;
			// TODO Auto-generated method stub
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

	}

	private static class TmpBamPileProducer implements Producer<TmpBamPile> {
		private SAMRecord samRecord;
		private Segment samRecordSegment;
		private BamPile[] bamPiles;
		private FilterNGS filterNGS;
		private Logger log;
		private int index = 0;

		private TmpBamPileProducer(SAMRecord samRecord, Segment samRecordSegment, BamPile[] bamPiles, FilterNGS filterNGS, Logger log) {
			super();
			this.samRecord = samRecord;
			this.bamPiles = bamPiles;
			this.filterNGS = filterNGS;
			this.samRecordSegment = samRecordSegment;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < bamPiles.length;
		}

		@Override
		public Callable<TmpBamPile> next() {
			TmpBamPile tmpBamPile = new TmpBamPile(samRecord, samRecordSegment, bamPiles[index], filterNGS, log);
			index++;
			return tmpBamPile;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class TmpBamPile implements Callable<TmpBamPile> {
		private SAMRecord samRecord;
		private BamPile bamPile;
		private FilterNGS filterNGS;
		private boolean overlaps;
		private Segment samRecordSegment;
		private Logger log;

		private TmpBamPile(SAMRecord samRecord, Segment samRecordSegment, BamPile bamPile, FilterNGS filterNGS, Logger log) {
			super();
			this.samRecord = samRecord;
			this.bamPile = bamPile;
			this.filterNGS = filterNGS;
			this.samRecordSegment = samRecordSegment;
			this.log = log;
		}

		@Override
		public TmpBamPile call() throws Exception {
			this.overlaps = samRecordSegment.overlaps(bamPile.getBin());
			if (overlaps) {
				bamPile.addRecord(samRecord, filterNGS.getPhreadScoreFilter(), log);
			} else {
				bamPile.summarize();
			}
			return this;
		}

		private BamPile getBamPile() {
			return bamPile;
		}

		private boolean overlapsCurrentRecord() {
			return overlaps;
		}

	}

	public static class bamPileWorker implements Callable<DynamicHistogram> {
		private String bamFile;
		private Segment[] q;
		private FilterNGS filterNGS;
		private ReferenceGenome referenceGenome;
		private int binSize;
		private PILE_TYPE type;
		private SAM_FILTER_TYPE filterType;
		private Logger log;

		public bamPileWorker(String bamFile, Segment[] q, FilterNGS filterNGS, ReferenceGenome referenceGenome, int binSize, PILE_TYPE type, SAM_FILTER_TYPE filterType, Logger log) {
			super();
			this.bamFile = bamFile;
			this.q = q;
			this.filterNGS = filterNGS;
			this.referenceGenome = referenceGenome;
			this.binSize = binSize;
			this.log = log;
			this.type = type;
			this.filterType = filterType;
		}

		@Override
		public DynamicHistogram call() throws Exception {
			BamPileUp pileUp = new BamPileUp(bamFile, referenceGenome, binSize, filterNGS, q, type, filterType, log);
			String output = ext.rootOf(bamFile, false) + ".bamPile.txt";
			PrintWriter writer = Files.getAppropriateWriter(output);
			writer.println(BamPile.getHeader());
			while (pileUp.hasNext()) {
				writer.println(pileUp.next().getOuput(log));
			}
			writer.close();
			pileUp.shutdown();
			return pileUp.getBamPileUpSummary().getProRefHistogram();
		}
	}

	public static class BamPileUpSummary {
		private DynamicHistogram proRefHistogram;
		private int readsPiled;
		private int totalReads;
		private int positionsPiled;
		private long time;
		private long totalTime;
		private Logger log;

		public BamPileUpSummary(Logger log) {
			super();
			this.proRefHistogram = new DynamicHistogram(0, 1, 2);
			this.readsPiled = 0;
			this.totalReads = 0;
			this.positionsPiled = 0;
			this.time = System.currentTimeMillis();
			this.totalTime = System.currentTimeMillis();
			this.log = log;
		}

		public DynamicHistogram getProRefHistogram() {
			return proRefHistogram;
		}

		public void setProRefHistogram(DynamicHistogram proRefHistogram) {
			this.proRefHistogram = proRefHistogram;
		}

		public void addToHistogram(double propRef) {
			proRefHistogram.addDataPointToHistogram(propRef);
		}

		public int getReadsPiled() {
			return readsPiled;
		}

		public void setReadsPiled(int readsPiled) {
			this.readsPiled = readsPiled;
		}

		public int getTotalReads() {
			return totalReads;
		}

		public void setTotalReads(int totalReads) {
			this.totalReads = totalReads;
		}

		public int getPositionsPiled() {
			return positionsPiled;
		}

		public void setPositionsPiled(int positionsPiled) {
			this.positionsPiled = positionsPiled;
		}

		public long getTime() {
			return time;
		}

		public void setTime(long time) {
			this.time = time;
		}

		public long getTotalTime() {
			return totalTime;
		}

		public void setTotalTime(long totalTime) {
			this.totalTime = totalTime;
		}

		public Logger getLog() {
			return log;
		}

		public void setLog(Logger log) {
			this.log = log;
		}

	}

}

// public DynamicHistogram pileUp() {
// if (pileType == PILE_TYPE.CONTAMINATION && referenceGenome == null) {
// String error = "A reference genome must be provided to detect reference alleles, cannot detect contamination in bam file " + bam;
// log.reportTimeError(error);
// return null;
// } else if (referenceGenome == null) {
// log.reportTimeWarning("A reference genome was not provided, reference and alternate alleles will not be recorded");
// }
//
// if (pileType == PILE_TYPE.CONTAMINATION && binSize != 1) {
// String error = "A bin size of one must be used for contamination detection on  " + bam;
// log.reportTimeError(error);
// return null;
// } else {
// log.reportTimeInfo("Using a bin size of " + binSize);
// }
// WorkerTrain<TmpBamPile> train = new WorkerTrain<BamPileUp.TmpBamPile>(null, 2, 200, log);// an extra thread should be about a half hour speed up per sample
// train.setAutoShutDown(false);
// TmpBamPileProducer tmpBamPileProducer = null;
// SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
// DynamicHistogram histogram = new DynamicHistogram(0, 1, 2);
// QueryInterval[] queryIntervals = null;
// if (intervals != null) {
// log.reportTimeInfo("Optimizing " + intervals.length + " queries for pile up");
// queryIntervals = BamOps.convertSegsToQI(intervals, reader.getFileHeader(), 0, true, log);
// log.reportTimeInfo("Finished Optimizing " + intervals.length + " queries to " + queryIntervals.length + " intervals for pile up");
// }
//
// SAMRecordIterator sIterator = queryIntervals == null ? reader.iterator() : reader.query(queryIntervals, false);
// AggregateFilter filter = initializeFilters(filterNGS, filterType, log);
// Segment currentSegment = new Segment((byte) 0, 0, 0);
// ArrayList<BamPile> bamPiles = new ArrayList<BamPile>();
//
// int totalReads = 0;
// int readsPiled = 0;
// int positionsPiled = 0;
// long time = System.currentTimeMillis();
// long totalTime = System.currentTimeMillis();
//
// String output = ext.rootOf(bam, false) + ".bamPile.txt";
// PrintWriter writer = Files.getAppropriateWriter(output);
// writer.println(BamPile.getHeader());
//
// while (sIterator.hasNext()) {
// totalReads++;
// SAMRecord samRecord = sIterator.next();
// if (!filter.filterOut(samRecord)) {
// readsPiled++;
// if (readsPiled % 100000 == 0) {
// log.reportTimeInfo("~" + readsPiled + " of " + totalReads + " total reads piled to " + positionsPiled + " positions (" + ext.getTimeElapsed(time) + ") " + SamRecordOps.getDisplayLoc(samRecord));
// time = System.currentTimeMillis();
// }
// BamPileInitializer bamPileInitializer = new BamPileInitializer(binSize, currentSegment, samRecord, log);
// while (bamPileInitializer.hasNext()) {
// bamPiles.add(bamPileInitializer.next());
// }
// currentSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
// tmpBamPileProducer = new TmpBamPileProducer(samRecord, currentSegment, bamPiles.toArray(new BamPile[bamPiles.size()]), filterNGS, log);
// train.setProducer(tmpBamPileProducer);
// bamPiles = new ArrayList<BamPile>(bamPiles.size());
// while (train.hasNext()) {
// TmpBamPile tmpBamPile = train.next();
// if (tmpBamPile.overlapsCurrentRecord()) {
// bamPiles.add(tmpBamPile.getBamPile());// stays in this round
// } else {
// BamPile bamPile = tmpBamPile.getBamPile();
// if (referenceGenome != null || pileType == PILE_TYPE.CONTAMINATION) {
// bamPile.setReference(referenceGenome);
// }
// if (filterNGS.getReadDepthFilter() == null || bamPile.getTotalDepth(false, false) > filterNGS.getReadDepthFilter()[0]) {
// int altAlleleDepth = filterNGS.getReadDepthFilter().length > 1 ? filterNGS.getReadDepthFilter()[1] : -1;
// if (pileType == PILE_TYPE.REGULAR || (bamPile.hasAltAllele(log) && bamPile.hasOnlyOneAlt(log) && bamPile.getNumAlt(log) > altAlleleDepth && bamPile.getNumRef(log) > altAlleleDepth)) {
// positionsPiled++;
// histogram.addDataPointToHistogram(bamPile.getPropRef(log));
// writer.println(bamPile.getOuput(log));
// }
// }
// }
// }
// }
// }
// sIterator.close();
// train.shutdown();
// log.reportTimeInfo(readsPiled + " reads piled in total (" + ext.getTimeElapsed(totalTime) + ")");
//
// try {
// reader.close();
// } catch (IOException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// writer.close();
// return histogram;
// }
