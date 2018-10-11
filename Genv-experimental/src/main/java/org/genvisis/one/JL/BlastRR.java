package org.genvisis.one.JL;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastResults;
import org.genvisis.seq.analysis.Blast.FastaEntry;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLineProcess;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.filesys.Segment;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Testing trimmed down blast - that does not rely on intermediate files
 */
public class BlastRR {

  public static void main(String[] args) {

    // where tmp file will end up
    String outDir = args[0];

    int segmentSize = 10000;
    int stepSize = 10000;

    new File(outDir).mkdirs();
    Logger logger = new Logger();
    ReferenceGenome rg = new ReferenceGenome(GenomeBuild.HG19, logger);
    SAMSequenceDictionary samSequenceDictionary = rg.getDictionary();
    for (SAMSequenceRecord record : samSequenceDictionary.getSequences()) {
      int chrLength = record.getSequenceLength();
      List<FastaEntry> fastaList = new LinkedList<>();
      for (int i = 0; i + segmentSize - 1 < chrLength; i += stepSize) {
        int start = i;
        int stop = i + segmentSize - 1;
        String sequence = ArrayUtils.toStr(rg.getSequenceFor(new Segment(record.getSequenceName(),
                                                                         start, stop)),
                                           "");
        if (sequence != null) {
          fastaList.add(new FastaEntry(record.getSequenceName() + "_" + start + "_" + stop,
                                       sequence));
        }
      }
      logger.reportTimeInfo("Blasting " + fastaList.size() + " entries from contig "
                            + record.getSequenceName());
      Blast blast = new Blast(rg.getReferenceFasta(), 500, 500, logger, false, false);
      blast.setEvalue(10000);

      FastaEntry[] fastaArray = fastaList.toArray(new FastaEntry[fastaList.size()]);

      CmdLineProcess cmdLineProcess = blast.getCmdLineProcess(fastaArray);
      int numRecords = 0;
      while (cmdLineProcess.hasNext()) {
        String line = cmdLineProcess.next();

        if (!line.startsWith("#")) {
          numRecords++;
          logger.reportTimeInfo("blast has returned " + numRecords + " matches for queries from "
                                + record.getSequenceName());
          String[] result = line.trim().split(PSF.Regex.GREEDY_WHITESPACE);
          BlastResults blastResults = new BlastResults(result, logger);
          logger.reportTimeInfo(blastResults.getQueryID() + " matched to "
                                + blastResults.getSubjectID() + " start=" + blastResults.getSstart()
                                + " stop=" + blastResults.getSstop());
        }
      }
      boolean error = cmdLineProcess.waitFor();
      if (error) {
        logger.reportError("Unsuccessful termination as indication by non-zero result code from \"blastn\" program.  Please investigate and try again.  "
                           + "If this error persists, or if you believe a non-zero response code from \"blastn\" is not irregular, please contact the Genvisis developers.");
      }

    }
  }

}
