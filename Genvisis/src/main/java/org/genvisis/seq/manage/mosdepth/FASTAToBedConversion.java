package org.genvisis.seq.manage.mosdepth;

import java.io.PrintWriter;
import org.genvisis.CLI;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;

class FASTAToBedConversion {

  private int binSize;
  private PrintWriter bedWriter;
  private PrintWriter dropWriter;
  private Logger log;

  public FASTAToBedConversion(int binSize, String outputBed, Logger log) {
    this.binSize = binSize;
    this.bedWriter = Files.getAppropriateWriter(outputBed);
    this.dropWriter = Files.getAppropriateWriter(ext.rootOf(outputBed, false) + ".drops.bed");
    this.log = log;
  }

  void run() {
    ReferenceGenome refGen = new ReferenceGenome(GENOME_BUILD.HG19, new Logger());
    Segment[] bins = refGen.getBins(binSize).getStrictSegments();
    for (Segment seg : bins) {
      String seq = ArrayUtils.toStr(refGen.getSequenceFor(seg));
      int n = 0;
      for (int i = 0; i < seq.length(); i++) {
        if (seq.charAt(i) == 'N') {
          n++;
        }
      }
      (n == seq.length() ? dropWriter : bedWriter).println(seg.getChr() + "\t"
                                                           + (seg.getStart() - 1) + "\t"
                                                           + seg.getStop() + "\t"
                                                           + seg.getUCSClocation());

    }
    bedWriter.close();
    log.reportTime("Done writing reference genome bins to BED format.");
  }

  public static void main(String[] args) {
    CLI cli = new CLI(FASTAToBedConversion.class);
    cli.addArg("bin", "Bin Size, default 1000.");
    cli.addArg(CLI.ARG_OUTFILE, CLI.DESC_OUTFILE);

    cli.parseWithExit(args);

    //    "G:\\bamTesting\\snpSelection\\ReferenceGenomeBins.bed"
    new FASTAToBedConversion(cli.getI("bin"), cli.get(CLI.ARG_OUTFILE), new Logger()).run();
  }

}
