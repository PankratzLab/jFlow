package org.genvisis.one.JL.aging.bam;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.telomere.TelSeq;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

import htsjdk.samtools.filter.SamRecordFilter;

/**
 * For comparing to TelSeq
 */
public class TestAging {

  private TestAging() {

  }
  // samtools view -s 0.1 -b SRR1738843.bam > Ten_Percent_SRR1738843.bam

  public static void main(String[] args) {
    String bam = "/Volumes/Beta/data/aric_sra/SRAPipeline/private/testBams/SRR1738843.bam";
    String outDir = ext.parseDirectoryOfFile(bam);
    Logger log = new Logger(outDir + "TL.test.log");
    String output = outDir + "tel.results.txt";
    TelSeq.telSeqIt(bam, output, (int) BamOps.estimateReadSize(bam, log).mean(),
                    new ArrayList<String>(), log);
    List<BamAnalysis> analyses = new ArrayList<>();

    PatternCounterTelomere pNgs = new PatternCounterTelomere(PatternCounterTelomere.getTelomericPattern(7),
                                                             PatternCounterTelomere.getDefaultGCTelRange(),
                                                             new ArrayList<SamRecordFilter>(),
                                                             new ArrayList<SamRecordFilter>());
    analyses.add(pNgs);

    BamAnalyzer bamAnalyzer = new BamAnalyzer(analyses, bam, log);
    bamAnalyzer.init();
    try {
      bamAnalyzer.analyze();
      bamAnalyzer.summmarize();
      for (BamAnalysis analysis : analyses) {
        Files.writeIterable(analysis.getOutput(), outDir + analysis.getRootOutputFile());
      }
    } catch (IOException e) {
      log.reportException(e);
    }

    // bamAnalyzer.

  }

}
