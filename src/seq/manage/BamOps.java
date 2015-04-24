package seq.manage;

import java.io.File;
import java.io.IOException;

import common.Logger;
import common.Positions;
import filesys.Segment;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Class for common bamFile manips
 *
 */
public class BamOps {

	/**
	 * @param bamOrSam
	 *            .bam or .sam file
	 * @param stringency
	 *            Stringency validation for the records
	 * @return new reader
	 */
	public static SamReader getDefaultReader(String bamOrSam, ValidationStringency stringency) {
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		samReaderFactory.validationStringency(stringency);
		SamReader reader = samReaderFactory.open(new File(bamOrSam));
		return reader;
	}

	/**
	 * @param segs
	 *            Genvisis type segments to search
	 * @param sFileHeader
	 *            an {@link SAMFileHeader}
	 * @param bpBuffer
	 *            bp buffer to be added to the segments
	 * @param optimize
	 *            perform an extra step to make a more efficient query
	 * @param log
	 * @return array of {@link QueryInterval} that can be queried by a bamfile reader
	 */
	public static QueryInterval[] convertSegsToQI(Segment[] segs, SAMFileHeader sFileHeader, int bpBuffer, boolean optimize, Logger log) {
		QueryInterval[] qIntervals = new QueryInterval[segs.length];
		segs = Segment.sortSegments(segs);
		for (int i = 0; i < qIntervals.length; i++) {
			String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), true);
			int referenceIndex = sFileHeader.getSequenceIndex(sequenceName);
			if (referenceIndex < 0) {
				log.reportError("Error - could not find " + sequenceName + " in the sequence dictionary, halting");
				return null;
			}
			qIntervals[i] = new QueryInterval(referenceIndex, segs[i].getStart() - bpBuffer, segs[i].getStop() + bpBuffer);
		}
		if (optimize) {
			qIntervals = QueryInterval.optimizeIntervals(qIntervals);
		}
		return qIntervals;
	}

	public static Segment[] converQItoSegs(QueryInterval[] qIntervals, SAMFileHeader sFileHeader, Logger log) {
		Segment[] segs = new Segment[qIntervals.length];
		for (int i = 0; i < segs.length; i++) {
			segs[i] = new Segment(Positions.chromosomeNumber(sFileHeader.getSequence(qIntervals[i].referenceIndex).getSequenceName()), qIntervals[i].start, qIntervals[i].end);
		}
		return segs;
	}

	public static SAMFileHeader getHeader(String bamfile) {
		SamReader reader = getDefaultReader(bamfile, ValidationStringency.STRICT);
		SAMFileHeader samFileHeader = reader.getFileHeader();
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return samFileHeader;

	}

	public static String getSampleName(String bamFile) {
		SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);
		String sample = reader.getFileHeader().getReadGroups().get(0).getSample();
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return sample;
	}
}
