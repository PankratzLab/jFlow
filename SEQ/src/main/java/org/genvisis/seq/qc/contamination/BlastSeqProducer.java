package org.genvisis.seq.qc.contamination;

import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastResultsSummary;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;

/**
 * Set up this since .bam files and .fastq files will have essentially identical functionality
 */
public abstract class BlastSeqProducer extends AbstractProducer<BlastResultsSummary[]> {

  protected Blast blast;
  protected int numSeqsPerThread;

  public BlastSeqProducer(String fastaDb, int blastWordSize, int reportWordSize,
                          int numSeqsPerThread, Logger log) {
    super();
    this.numSeqsPerThread = numSeqsPerThread;
    blast = new Blast(fastaDb, blastWordSize, reportWordSize, log, true, false);
    blast.setTaxonMode(true);
  }

}