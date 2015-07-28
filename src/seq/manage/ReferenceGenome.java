package seq.manage;

import filesys.Segment;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import common.Array;
import common.Array.BYTE_DECODE_FORMAT;
import common.Logger;
import common.Positions;

public class ReferenceGenome {
	public static final String DEFAULT_REFERENCE = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
	private String referenceFasta;
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private byte[] currentSeq;
	private int defaultBuffer;
	private ReferenceSequence referenceSequence;
	private Logger log;

	public ReferenceGenome(String referenceFasta, Logger log) {
		super();
		this.referenceFasta = referenceFasta;
		this.log = log;
		this.indexedFastaSequenceFile = null;
		this.referenceSequence = null;
		this.defaultBuffer = 0;
		try {
			indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(referenceFasta));
		} catch (FileNotFoundException e) {
			log.reportFileNotFound(referenceFasta);
			e.printStackTrace();
		}
	}

	public IndexedFastaSequenceFile getIndexedFastaSequenceFile() {
		return indexedFastaSequenceFile;
	}

	public String getReferenceFasta() {
		return referenceFasta;
	}

	public ReferenceSequence getReferenceSequence() {
		return referenceSequence;
	}

	public int getDefaultBuffer() {
		return defaultBuffer;
	}

	public void setDefaultBuffer(int defaultBuffer) {
		this.defaultBuffer = defaultBuffer;
	}

	/**
	 * @param segs
	 *            these should be sorted in chromosomal order if speed is desired
	 * @return
	 */
	public String[][] getSequencesFor(Segment[] segs) {
		String[][] seqs = new String[segs.length][];
		for (int i = 0; i < seqs.length; i++) {
			seqs[i] = getSequenceFor(segs[i]);
		}
		return seqs;
	}

	public boolean hasContig(String contig) {
		return indexedFastaSequenceFile.getSequenceDictionary().getSequence(contig) != null;
	}

	/**
	 * @param segment
	 * @return will return null if {@link Positions#getChromosomeUCSC(int, boolean, boolean)} returns a contig not in the files {@link SAMSequenceDictionary}
	 */
	public String[] getSequenceFor(Segment segment) {
		String requestedContig = Positions.getChromosomeUCSC(segment.getChr(), true, true);
		if (hasContig(requestedContig)) {
			ReferenceSequence referenceSequence = indexedFastaSequenceFile.getSubsequenceAt(requestedContig, segment.getStart() - defaultBuffer, segment.getStop() + defaultBuffer);
			String[] requestedSeq = Array.decodeByteArray(referenceSequence.getBases(), BYTE_DECODE_FORMAT.UPPER_CASE, log);
			return requestedSeq;

		} else {
			// log.reportTimeError("Requested contig " + requestedContig + " was not in the sequence dictionary for " + referenceFasta);
			return null;
		}

	}

	/**
	 * This method is slow, it loads the complete chromosome into memory and then subsets.
	 */
	@Deprecated
	public String[] getSequenceForOld(Segment segment) {
		String requestedContig = Positions.getChromosomeUCSC(segment.getChr(), true);
		if (referenceSequence == null || currentSeq == null || !referenceSequence.getName().equals(requestedContig)) {
			referenceSequence = indexedFastaSequenceFile.getSequence(requestedContig);
			currentSeq = referenceSequence.getBases();
		}
		int start = segment.getStart() - 1 - defaultBuffer;
		int stop = segment.getStop() + defaultBuffer;
		if (start < 0) {
			log.reportTimeWarning("Buffer of " + defaultBuffer + " adjusts base pair extraction to index less than 0, adjusting start to index 0 (bp 1)");
			start = 0;
		}
		if (stop >= currentSeq.length) {
			stop = currentSeq.length - 1;
			log.reportTimeWarning("Buffer of " + defaultBuffer + " ,adjusts base pair extraction to index greater than sequence length,  adjusting stop to index " + (currentSeq.length - 1));
		}
		byte[] subTmp = null;
		try {
			subTmp = Array.subArray(currentSeq, start, stop);
		} catch (Exception e) {
			log.reportTimeError("Could not extract bases:");
			log.reportTimeError("Segment: " + segment.getUCSClocation());
			log.reportTimeError("Ref: " + referenceSequence.getName());
			log.reportTimeError("Start : " + start);
			log.reportTimeError("Stop : " + stop);
			return new String[] {};
		}
		if (!referenceSequence.getName().equals(requestedContig)) {
			log.reportTimeError("Mismatched request");
			log.reportTimeError("Segment: " + segment.getUCSClocation());
			log.reportTimeError("Ref: " + referenceSequence.getName());
			log.reportTimeError("Start : " + start);
			log.reportTimeError("Stop : " + stop);
			return new String[] {};
		}

		String[] requestedSeq = new String[subTmp.length];
		try {
			for (int i = 0; i < requestedSeq.length; i++) {
				requestedSeq[i] = new String(new byte[] { subTmp[i] }, "UTF-8").toUpperCase();
			}
		} catch (UnsupportedEncodingException e) {
			log.reportTimeError("Could not convert reference byte to string");
			e.printStackTrace();
		}
		return requestedSeq;
	}

	public double getGCContentFor(Segment seg) {
		String[] seq = getSequenceFor(seg);
		int gs = Array.countIf(seq, "G");
		int cs = Array.countIf(seq, "C");
		int gsCs = gs + cs;
		return (double) gsCs / seq.length;
	}

	public double getGCContentFor(VariantContext vc) {
		// defaultBuffer=0;
		// System.out.println("REF\t"+Array.toStr(getSequenceFor(VCOps.getSegment(vc))));
		// System.out.println("VARIANT_A\t"+vc.getAlleles().toString());
		return getGCContentFor(VCOps.getSegment(vc));
	}

}
