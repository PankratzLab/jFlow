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
import org.genvisis.seq.analysis.Blast.FastaEntry;
import org.genvisis.seq.manage.ReferenceGenome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Testing trimmed down blast
 *
 */
public class BlastRR {

	public static void main(String[] args) {

		// where tmp file will end up
		String outDir = args[0];

		int segmentSize = 1000;
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

			String tmpFile = outDir + record.getSequenceName() + ".blast.tmp.gz";
			PrintWriter tmpWriter = Files.getAppropriateWriter(tmpFile);
			FastaEntry[] fastaArray = fastaList.toArray(new FastaEntry[fastaList.size()]);

			blast.blastSequence(fastaArray, tmpWriter);
			tmpWriter.close();

			try {
				BufferedReader reader = Files.getAppropriateReader(tmpFile);
				while (reader.ready()) {
					BlastResults blastResults = new BlastResults(reader.readLine().trim().split("\t"), logger);
					System.out.println(blastResults.getQueryID() + " matched to " + blastResults.getSubjectID()
							+ " start=" + blastResults.getSstart() + " stop=" + blastResults.getSstop());

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
