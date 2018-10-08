package org.genvisis.one.JL.ngs.pca;

import java.io.File;
import java.io.IOException;
import java.util.StringJoiner;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.manage.ReferenceGenome;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CLI.Arg;
import org.pankratzlab.shared.filesys.Segment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.util.CloseableIterator;

public class NGSCount {

  private static final String EXTENSION = "extension";
  private static final String BIN_SIZE = "binSize";

  private static class BinCount extends Segment {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    private int readCount;

    /**
     * @param chr
     * @param start
     * @param stop
     */
    public BinCount(byte chr, int start, int stop) {
      super(chr, start, stop);
      // TODO Auto-generated constructor stub
    }

  }

  private static boolean statFile(String file, LocusSet<Segment> bins, FilterNGS filterNGS,
                                  Logger log) {

    SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    samReaderFactory.validationStringency(ValidationStringency.STRICT);

    try (SamReader reader = samReaderFactory.open(new File(file))) {
      AggregateFilter filter = FilterNGS.initializeFilters(filterNGS, SAM_FILTER_TYPE.COPY_NUMBER,
                                                           log);

      int counted = 0;
      for (int i = 0; i < bins.getLoci().length; i++) {

        Segment bin = bins.getLoci()[i];
        if (bin.getChr() > 0) {
          long time = System.currentTimeMillis();
          log.reportTimeInfo("Opening" + bin.getUCSClocation() + "\t" + counted + "\t" + i);

          CloseableIterator<SAMRecord> iterator = reader.queryOverlapping(bin.getChromosomeUCSC(),
                                                                          bin.getStart(),
                                                                          bin.getStop());

          log.reportTimeElapsed(time);

          log.reportTimeInfo("Scanning" + bin.getUCSClocation() + "\t" + counted + "\t" + i);
          time = System.currentTimeMillis();
          while (iterator.hasNext()) {
            SAMRecord samRecord = iterator.next();
            if (!filter.filterOut(samRecord)) {
              counted++;
              //            currentPile.addRecord(samRecord, currentRef, filterNGS.getPhreadScoreFilter(), log);
            }
          }
          log.reportTimeElapsed(time);
          log.reportTimeInfo("Finished Scanning" + bin.getUCSClocation() + "\t" + counted + "\t"
                             + i);

          iterator.close();
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    return true;
  }

  private static void run(CLI c) {
    System.setProperty("samjdk.reference_fasta", c.get(CLI.ARG_REFERENCE_GENOME));
    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "counts.log");
    ReferenceGenome referenceGenome = new ReferenceGenome(c.get(CLI.ARG_REFERENCE_GENOME), log);
    LocusSet<Segment> bins = referenceGenome.getBins(c.getI(BIN_SIZE));
    String[] files = Files.listFullPaths(c.get(CLI.ARG_INDIR), c.get(EXTENSION));
    FilterNGS filterNGS = new FilterNGS(20, 20, null);

    String outBed = outDir + ext.rootOf(c.get(CLI.ARG_REFERENCE_GENOME)) + "." + c.getI(BIN_SIZE)
                    + ".bed";
    StringJoiner joiner = new StringJoiner("\n");
    for (Segment seg : bins.getLoci()) {
      if (seg.getChr() > 0) {
        joiner.add(seg.getChromosomeUCSC() + "\t" + seg.getStart() + "\t" + seg.getStop());
      }
    }
    Files.write(joiner.toString(), outBed);
    System.exit(1);

    for (String file : files) {
      statFile(file, bins, filterNGS, log);
    }

  }

  public static void main(String[] args) {
    CLI c = new CLI(NGSCount.class);
    c.addArg(CLI.ARG_INDIR, "input directory containing bam or cram files", "dir/", true,
             Arg.STRING);

    c.addArg(CLI.ARG_REFERENCE_GENOME, CLI.DESC_REFERENCE_GENOME, "my.ref.fa", true, Arg.FILE);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "out/", true, Arg.FILE);
    c.addArgWithDefault(BIN_SIZE, "Size of bins in bp", "1000");
    c.addArgWithDefault(EXTENSION, EXTENSION, ".cram");

    c.parseWithExit(args);

    run(c);

  }

}
