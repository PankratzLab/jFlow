package org.genvisis.one.JL.aging.bam;

import java.io.IOException;
import java.util.List;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.SamRecordOps;
import org.genvisis.seq.qc.FlagStats;
import org.pankratzlab.common.Logger;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

/**
 * Class for performing analyses (implementing {@link BamAnalysis} of bams.
 */
public class BamAnalyzer {

  private List<BamAnalysis> analyses;
  private String bamFile;
  private FlagStats fullStats;
  private Logger log;

  /**
   * @param analyses {@link BamAnalysis} to perform
   * @param bamFile bam file to analyze
   * @param log
   */
  public BamAnalyzer(List<BamAnalysis> analyses, String bamFile, Logger log) {
    super();
    this.analyses = analyses;
    this.bamFile = bamFile;
    this.fullStats = new FlagStats();
    this.log = log;
  }

  /**
   * Initialize the analyses
   */
  public void init() {
    for (BamAnalysis bamAnalysis : analyses) {
      bamAnalysis.init(bamFile, log);
    }
  }

  /**
   * Calls {@link BamAnalysis#summarize()} for each analyiss
   */
  public void summmarize() {
    for (BamAnalysis bamAnalysis : analyses) {
      bamAnalysis.summarize();
    }
  }

  /**
   * Analyze the current .bam file. Calls {@link BamAnalysis#shouldAnalyze(SAMRecord, Logger)} and
   * {@link BamAnalysis#analyze(SAMRecord, Logger)} on each {@link SAMRecord} in the file
   * 
   * @throws IOException
   */
  public void analyze() throws IOException {
    SamReader reader = BamOps.getDefaultReader(bamFile, ValidationStringency.LENIENT);
    String sample = BamOps.getSampleName(bamFile, log);
    int num = 0;
    for (SAMRecord samRecord : reader) {
      fullStats.stat(samRecord);
      if (num % 1000000 == 0) {
        log.reportTimeInfo("Processing reads for sample " + sample + " , scanned " + num
                           + ", currently on " + SamRecordOps.getDisplayLoc(samRecord));
      }
      num++;
      for (BamAnalysis bamAnalysis : analyses) {
        if (bamAnalysis.shouldAnalyze(samRecord, log)) {
          bamAnalysis.analyze(samRecord, log);
        }
      }
    }
    reader.close();

  }
}
