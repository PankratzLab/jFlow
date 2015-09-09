package seq.manage;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;

import java.util.ArrayList;
import java.util.Iterator;

import seq.qc.FilterNGS;
import seq.qc.FilterNGS.SAM_FILTER_TYPE;
import common.Array;
import common.Logger;
import filesys.Segment;

/**
 * @author lane0212 New version of the pileup
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

	public BamSegPileUp(String bam, ReferenceGenome referenceGenome, Segment[] intervals, Logger log) {
		super();
		this.bam = bam;
		this.numReturned = 0;
		this.referenceGenome = referenceGenome;
		this.intervals = intervals;
		this.reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
		this.log = log;
		this.queryIntervals = BamOps.convertSegsToQI(intervals, reader.getFileHeader(), 0, true, log);
		log.reportTimeInfo("Optimizing " + intervals.length + " queries for pile up");
		this.sIterator = reader.query(queryIntervals, false);
		this.addingMask = new boolean[intervals.length];
		this.bamPiles = new BamPile[intervals.length];
		this.bamPilesToReturn = new ArrayList<BamPile>();
		this.filter = new AggregateFilter(new FilterNGS().getStandardSAMRecordFilters(SAM_FILTER_TYPE.COPY_NUMBER, log));
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

				for (int i = 0; i < addingMask.length; i++) {

					if (addingMask[i]) {

						boolean overlaps = samRecordSegment.overlaps(bamPiles[i].getBin());

						if (overlaps) {
							bamPiles[i].addRecord(samRecord, 0, log);
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

}
