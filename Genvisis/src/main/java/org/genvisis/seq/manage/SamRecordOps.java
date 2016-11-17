package org.genvisis.seq.manage;

import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.Segment;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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
			log.reportError("Could not determing start and stop for " + samRecord.toString());
		}

		return new Segment(chr, start, stop);
	}

	public static boolean overlaps(Segment seg, SAMRecord samRecord, Logger log) {
		return seg.overlaps(getReferenceSegmentForRecord(samRecord, log));
	}

	public static byte getChromosome(SAMRecord samRecord, Logger log) {
		byte chr = Positions.chromosomeNumber(samRecord.getReferenceName());
		if (chr < 1) {
			log.reportError("Could not determine chromosome for " + samRecord.toString());
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

	public static class SoftClipped {
		private final String bases;
		private final Segment refSeg;// for soft clipped, this should always be a 1 bp seg

		public SoftClipped(String bases, Segment refSeg) {
			super();
			this.bases = bases;
			this.refSeg = refSeg;
		}

		public String getBases() {
			return bases;
		}

		public Segment getRefSeg() {
			return refSeg;
		}

	}

	/**
	 * @param samRecord
	 * @param log
	 * @return String[] containing any bp sequences that were soft clipped
	 */
	public static SoftClipped[] getSoftClippedBases(SAMRecord samRecord, Logger log) {
		ArrayList<SoftClipped> softies = new ArrayList<SoftClipped>();
		Cigar cigar = samRecord.getCigar();
		String[] bases = Array.decodeByteArray(samRecord.getReadBases(), log);
		int curStart = 0;
		int readIndex = 0;

		byte chr = getChromosome(samRecord, log);
		for (CigarElement cigarElement : cigar.getCigarElements()) {
			if (cigarElement.getOperator().consumesReadBases()) {
				readIndex += cigarElement.getLength();
			}
			if (cigarElement.getOperator() == CigarOperator.S) {
				String softy = Array.toStr(Array.subArray(bases, curStart, readIndex), "");
				// soft clips are at begining and ends of reads....
				int refPos = curStart == 0	? samRecord.getReferencePositionAtReadPosition(readIndex)
																		: samRecord.getReferencePositionAtReadPosition(curStart);
				softies.add(new SoftClipped(softy, new Segment(chr, refPos, refPos)));
			}
			if (cigarElement.getOperator().consumesReadBases()) {
				curStart += cigarElement.getLength();
			}
		}
		return softies.toArray(new SoftClipped[softies.size()]);
	}

}
