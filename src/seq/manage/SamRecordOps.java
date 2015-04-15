package seq.manage;

import common.Logger;
import common.Positions;
import filesys.Segment;
import htsjdk.samtools.SAMRecord;

/**
 * @author lane0212 Class for common manips of {@link SAMRecord}
 */

public class SamRecordOps {

	public static Segment getReferenceSegmentForRecord(SAMRecord samRecord, Logger log) {
		byte chr = getChromosome(samRecord, log);
		int start = samRecord.getAlignmentStart();
		int stop = samRecord.getAlignmentEnd();
		if (start == 0 || stop == 0) { 
			log.reportTimeError("Could not determing start and stop for " + samRecord.toString());
		}

		return new Segment(chr, start, stop);
	}

	public static boolean overlaps(Segment seg, SAMRecord samRecord, Logger log) {
		return seg.overlaps(getReferenceSegmentForRecord(samRecord, log));
	}

	public static byte getChromosome(SAMRecord samRecord, Logger log) {
		byte chr = Positions.chromosomeNumber(samRecord.getReferenceName());
		if (chr < 1) {
			log.reportTimeError("Could not determine chromosome for " + samRecord.toString());
		}
		return chr;
	}

	public static String getDisplayLoc(SAMRecord samRecord) {
		return "CHR: " + samRecord.getReferenceName() + " Position: " + samRecord.getAlignmentStart();
	}

	public static String[] getReadBases(SAMRecord samRecord) {
		String s = samRecord.getReadString();
		String[] sa = new String[s.length()];

		for (int i = 0; i < s.length(); i++) {
			sa[i] = s.charAt(i) + "";
		}
		return sa;
	}

	public static double[] getReadPhred(SAMRecord samRecord) {
		byte[] p = samRecord.getBaseQualities();
		double[] d = new double[p.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = p[i];
		}
		return d;
	}
}
