package seq.manage;

import filesys.Segment;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import common.Array;
import common.Logger;
import common.Positions;

public class ReferenceGenome {
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
		try {
			indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(referenceFasta));
			referenceSequence = indexedFastaSequenceFile.nextSequence();
			currentSeq = referenceSequence.getBases();
		} catch (FileNotFoundException e) {
			log.reportFileNotFound(referenceFasta);
			e.printStackTrace();
		}
	}
	
	

	public IndexedFastaSequenceFile getIndexedFastaSequenceFile() {
		return indexedFastaSequenceFile;
	}



	public void setIndexedFastaSequenceFile(
			IndexedFastaSequenceFile indexedFastaSequenceFile) {
		this.indexedFastaSequenceFile = indexedFastaSequenceFile;
	}



	public ReferenceSequence getReferenceSequence() {
		return referenceSequence;
	}



	public void setReferenceSequence(ReferenceSequence referenceSequence) {
		this.referenceSequence = referenceSequence;
	}



	public int getDefaultBuffer() {
		return defaultBuffer;
	}

	public void setDefaultBuffer(int defaultBuffer) {
		this.defaultBuffer = defaultBuffer;
	}

	public String[] getSequenceFor(Segment segment) {
		String requestedContig = Positions.getChromosomeUCSC(segment.getChr(), true);
		if (!referenceSequence.getName().equals(requestedContig)) {
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
