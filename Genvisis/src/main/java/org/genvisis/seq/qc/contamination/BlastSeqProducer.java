package org.genvisis.seq.qc.contamination;

import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastResultsSummary;

/**
 * Set up this since .bam files and .fastq files will have essentially identical functionality
 *
 */
public abstract class BlastSeqProducer implements Producer<BlastResultsSummary[]> {
	protected Blast blast;
	protected int numSeqsPerThread;

	public BlastSeqProducer(String fastaDb, int blastWordSize, int reportWordSize, int numSeqsPerThread, Logger log) {
		super();
		this.numSeqsPerThread = numSeqsPerThread;
		this.blast = new Blast(fastaDb, blastWordSize, reportWordSize, log, true, false);
		blast.setTaxonMode(true);
	}

}
