package seq.manage;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

import java.util.ArrayList;
import java.util.Iterator;

import cnv.var.LocusSet;
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
	private int currentIndex;
	private ArrayList<BamPile> bamPilesToReturn;
	private BamPile[] bamPiles;
	private SamReader reader;
	private Logger log;
	private SAMRecordIterator sIterator;
	private boolean[] addingMask;

	public BamSegPileUp(String bam, ReferenceGenome referenceGenome, Segment[] intervals, Logger log) {
		super();
		this.bam = bam;
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
		for (int i = 0; i < intervals.length; i++) {
			bamPiles[i] = new BamPile(intervals[i]);
			addingMask[i] = true;
		}

	}

	@Override
	public boolean hasNext() {
		while (sIterator.hasNext() && bamPilesToReturn.size() == 0) {
			SAMRecord samRecord = sIterator.next();
			System.out.println(samRecord.toString());
			Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);

			for (int i = 0; i < addingMask.length; i++) {

				if (addingMask[i]) {

					boolean overlaps = samRecordSegment.overlaps(bamPiles[i].getBin());
					
					if (overlaps) {
						bamPiles[i].addRecord(samRecord, 0, log);
					} else {
						if (samRecordSegment.getStart() > bamPiles[i].getBin().getStop()||!sIterator.hasNext()) {
							System.out.println("DHFDSF");
							System.exit(1);
							addingMask[i] = false;
							bamPilesToReturn.add(bamPiles[i]);
						}
					}
				}

			}

		}
		return bamPilesToReturn.size() > 0;
	}

	@Override
	public BamPile next() {
		// TODO Auto-generated method stub
		return bamPilesToReturn.get(0);
	}

}
