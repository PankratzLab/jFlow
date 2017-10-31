package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.List;

import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastResults;
import org.genvisis.seq.analysis.Blast.BlastResultsSummary;
import org.genvisis.seq.analysis.Blast.FastaEntry;
import org.genvisis.seq.manage.ReferenceGenome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Testing trimmed down blast
 *
 */
class BlastRR {

	public static void main(String[] args) {

		// where tmp file will end up
		String outDir = args[0];

		int segmentSize = 100;
		int stepSize = 100;

		new File(outDir).mkdirs();
		Logger logger = new Logger();
		ReferenceGenome rg = new ReferenceGenome(GENOME_BUILD.HG19, logger);
		SAMSequenceDictionary samSequenceDictionary = rg.getDictionary();

		for (SAMSequenceRecord record : samSequenceDictionary.getSequences()) {
			int chrLength = record.getSequenceLength();
			List<FastaEntry> fastaList = new LinkedList<>();
			for (int i = 0; i + segmentSize - 1 < chrLength; i += stepSize) {
				int start = i;
				int stop = i + segmentSize - 1;
				String sequence = ArrayUtils
						.toStr(rg.getSequenceFor(new Segment(record.getSequenceName(), start, stop)), "");
				if (sequence != null) {
					fastaList.add(new FastaEntry(record.getSequenceName() + "_" + start + "_" + stop, sequence));
				}
			}
			Blast blast = new Blast(rg.getReferenceFasta(), 50, 50, logger, false, false);
			blast.setEvalue(10000);

			// Previously, without blast.setEvalue(XXXX) set, we get

			// [Genvisis 0.0.0-unknown] 12:09:10PM Error - Error: Too many
			// positional arguments (1), the offending value:
			// [Genvisis 0.0.0-unknown] 12:09:10PM Error - Error:
			// (CArgException::eSynopsis) Too many positional arguments (1), the
			// offending value:
			// [Genvisis 0.0.0-unknown] 12:09:10PM Error - Irregular termination
			// in process: blastn -db
			// /Users/Kitty/.genvisis/resources/Genome/hg19/hg19.fa -outfmt 7
			// std btop -word_size 50 | java.lang.UNIXProcess@7c53a9eb
			// [Genvisis 0.0.0-unknown] 12:09:10PM Error - Unsuccessful
			// termination as indication by non-zero result code from "blastn"
			// program. Please investigate and try again. If this error
			// persists, or if you believe a non-zero response code from
			// "blastn" is not irregular, please contact the Genvisis
			// developers.

			String tmpFile = outDir + record.getSequenceName() + ".blast.tmp.gz";
			PrintWriter tmpWriter = Files.getAppropriateWriter(tmpFile);
			FastaEntry[] fastaArray = fastaList.toArray(new FastaEntry[fastaList.size()]);
			BlastResultsSummary[] blasts = blast.blastSequence(fastaArray, tmpWriter);
			tmpWriter.close();

			try {
				BufferedReader reader = Files.getAppropriateReader(tmpFile);
				while (reader.ready()) {
					BlastResults blastResults = new BlastResults(reader.readLine().trim().split("\t"), logger);
					System.out.println(blastResults.getBtop());
					// if(!blastResults.)
				}

				reader.close();
			} catch (FileNotFoundException e) {
				logger.reportException(e);
			} catch (IOException e) {
				logger.reportException(e);

			}

		}
	}

}
