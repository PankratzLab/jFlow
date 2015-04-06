package seq.manage;

import filesys.Segment;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

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

	public String[] getSequenceFor(Segment segment) {
		String requestedContig = Positions.getChromosomeUCSC(segment.getChr(), true);
		if (!referenceSequence.getName().equals(requestedContig)) {
			referenceSequence = indexedFastaSequenceFile.getSequence(requestedContig);
			currentSeq = referenceSequence.getBases();
		}
		byte[] subTmp = Array.subArray(currentSeq, segment.getStart() - 1, segment.getStop());
		if(subTmp.length>1){
			System.out.println(segment.getUCSClocation());
			System.exit(1);
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

}
