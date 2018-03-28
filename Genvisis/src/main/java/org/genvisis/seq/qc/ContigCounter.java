/**
 * 
 */
package org.genvisis.seq.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.seq.manage.BamOps;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Class to get basic stats on each contig of a bam file
 */
public class ContigCounter {

  /**
   * @param bamFile generate contig stats for this bam
   * @param log
   * @return map of contig -> contig's stats
   */
  private static ContigStatResults statBamContigs(String bamFile, Logger log) {
    SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    samReaderFactory.validationStringency(ValidationStringency.LENIENT);
    SamReader reader = samReaderFactory.open(new File(bamFile));
    SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();

    Map<String, FlagStats> stats = new LinkedHashMap<>();
    for (SAMSequenceRecord record : dictionary.getSequences()) {
      SAMRecordIterator iter = reader.query(record.getSequenceName(), 0, record.getSequenceLength(),
                                            true);
      log.reportTimeInfo("Gathering stats for " + bamFile + " on contig "
                         + record.getSequenceName());
      FlagStats flagStats = new FlagStats();
      while (iter.hasNext()) {
        SAMRecord samRecord = iter.next();
        flagStats.stat(samRecord);
      }
      iter.close();
      stats.put(record.getSequenceName(), flagStats);
    }
    try {
      reader.close();
    } catch (IOException e) {
      log.reportException(e);
    }
    return new ContigStatResults(stats, BamOps.getSampleName(bamFile, log));
  }

  private static class ContigStatResults {

    private final Map<String, FlagStats> stats;
    private final String sample;

    /**
     * @param stats stats per contig
     * @param sample sample from bam
     */
    private ContigStatResults(Map<String, FlagStats> stats, String sample) {
      super();
      this.stats = stats;
      this.sample = sample;
    }

  }
  /**
   * thread the contig stat gathering
   */
  private static class ContigStatProducer implements Producer<ContigStatResults> {

    private int index;
    private final String[] bamFiles;
    private final Logger log;

    /**
     * @param bamFiles
     * @param log
     */
    public ContigStatProducer(String[] bamFiles, Logger log) {
      super();
      this.index = 0;
      this.bamFiles = bamFiles;
      this.log = log;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext() {
      return index < bamFiles.length;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    @Override
    public Callable<ContigStatResults> next() {
      final String bam = bamFiles[index];
      index++;
      return () -> statBamContigs(bam, log);
    }

    /*
     * (non-Javadoc)
     * @see org.genvisis.common.WorkerTrain.Producer#shutdown()
     */
    @Override
    public void shutdown() {

    }
  }

  /**
   * @param bamDir run contig stats on all bams in this directory
   * @param outDir output directory
   * @param threads number of threads, one per bam
   */
  private static void run(String bamDir, String outDir, int threads) {
    String rootOut = outDir + "ContigStats";
    new File(outDir).mkdirs();
    Logger log = new Logger(rootOut + ".log");
    String outputFile = rootOut + ".txt";
    String[] bamFiles = Files.listFullPaths(bamDir, ".bam");
    log.reportTimeInfo("Found " + bamFiles.length + " bams in " + bamDir);
    ContigStatProducer producer = new ContigStatProducer(bamFiles, log);

    log.reportTimeInfo("Output: " + outputFile);
    try (WorkerTrain<ContigStatResults> train = new WorkerTrain<>(producer, threads, 2, log);
         PrintWriter writer = Files.getAppropriateWriter(outputFile)) {

      boolean writeHeader = true;

      while (train.hasNext()) {
        ContigStatResults stats = train.next();
        for (String contig : stats.stats.keySet()) {
          if (writeHeader) {
            writer.println("SAMPLE\tCONTIG\t" + ArrayUtils.toStr(FlagStats.getHeader()));
            writeHeader = false;
          }
          writer.println(stats.sample + "\t" + contig + "\t"
                         + ArrayUtils.toStr(stats.stats.get(contig).getData()));
        }
      }
    }
  }

  public static void main(String[] args) {
    CLI c = new CLI(ContigCounter.class);
    c.addArg("bams", "directory of bams");
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "24");
    c.parseWithExit(args);
    run(c.get("bams"), c.get(CLI.ARG_OUTDIR), c.getI(CLI.ARG_THREADS));
  }
}
