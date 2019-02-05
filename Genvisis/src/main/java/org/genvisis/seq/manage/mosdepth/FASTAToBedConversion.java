package org.genvisis.seq.manage.mosdepth;

import java.io.PrintWriter;
import org.genvisis.cnv.Resources;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.ReferenceGenome;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

public class FASTAToBedConversion {

  private int binSize;
  private PrintWriter bedWriter;
  private PrintWriter dropWriter;
  private Logger log;
  private GenomeBuild build;

  public FASTAToBedConversion(int binSize, String outputBed, GenomeBuild build, Logger log) {
    this.binSize = binSize;
    this.build = build;
    this.bedWriter = Files.getAppropriateWriter(outputBed);
    this.dropWriter = Files.getAppropriateWriter(ext.rootOf(outputBed, false) + ".drops.bed");
    this.log = log;
  }

  void run() {
    ReferenceGenome refGen = new ReferenceGenome(Resources.genome(build, log).getFASTA()
                                                          .getAbsolute(),
                                                 log);
    Segment[] bins = refGen.getBins(binSize).getStrictSegments();
    for (Segment seg : bins) {
      String[] bases = refGen.getSequenceFor(seg);
      if (bases == null) {
        dropWriter.println(seg.getChr() + "\t" + (seg.getStart() - 1) + "\t" + seg.getStop() + "\t"
                           + seg.getUCSClocation());
        continue;
      }
      String seq = ArrayUtils.toStr(bases);
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
    dropWriter.close();
    log.reportTime("Done writing reference genome bins to BED format.");
  }

  public static void main(String[] args) {
    CLI cli = new CLI(FASTAToBedConversion.class);
    cli.addArg("bin", "Bin Size, default 1000.");
    cli.addArg(CLI.ARG_OUTFILE, CLI.DESC_OUTFILE);
    cli.addArg("build",
               "GenomeBuild, one of " + ArrayUtils.toStr(GenomeBuild.values(), ", ") + ".");

    cli.parseWithExit(args);
    new FASTAToBedConversion(cli.getI("bin"), cli.get(CLI.ARG_OUTFILE),
                             GenomeBuild.valueOf(cli.get("build").toUpperCase()), new Logger())
                                                                                               .run();
  }

}
