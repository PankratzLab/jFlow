package org.genvisis.seq.qc.contamination;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastResultsSummary;
import org.genvisis.seq.analysis.Blast.FastaEntry;

/**
 * @author lane0212 {@link BlastSeqProducer} for fasta type files
 */
public class BlastFastaq extends BlastSeqProducer {
	private FastqReader reader;
	private int numToTest, numTested;

	/**
	 * @param fastaqFile
	 *            file with sequence to blast
	 * @param numToTest
	 *            number of sequences (in order of file ) to blast
	 * @param fastaDb
	 *            the blast db file to blast against
	 * @param blastWordSize
	 *            the word size for initial matches
	 * @param reportWordSize
	 *            word size to be reported
	 * @param numSeqsPerThread
	 *            number of sequences given to each thread
	 * @param log
	 */
	public BlastFastaq(String fastaqFile, int numToTest, String fastaDb, int blastWordSize, int reportWordSize, int numSeqsPerThread, Logger log) {
		super(fastaDb, blastWordSize, reportWordSize, numSeqsPerThread, log);
		this.numTested = 0;
		this.numToTest = numToTest;
		if (!Files.exists(fastaqFile)) {
			log.reportFileNotFound(fastaqFile);
		}
		this.reader = new FastqReader(new File(fastaqFile), true);
	}

	@Override
	public void shutdown() {
		reader.close();

	}

	@Override
	public boolean hasNext() {
		return reader.hasNext() && numTested < numToTest;
	}

	@Override
	public Callable<BlastResultsSummary[]> next() {
		int numAdded = 0;
		ArrayList<FastaEntry> curEntries = new ArrayList<FastaEntry>();
		while (reader.hasNext() && numAdded < numSeqsPerThread) {
			numTested++;
			FastqRecord record = reader.next();
			String safeRecord = ext.replaceWithLinuxSafeCharacters(record.getReadHeader(), true);// everything after a space is truncated otherwise
			FastaEntry entry = new FastaEntry(safeRecord, record.getReadString());
			curEntries.add(entry);
			numAdded++;
		}
		return new Blast.BlastWorker(blast, curEntries.toArray(new FastaEntry[curEntries.size()]), null);
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub

	}

}
