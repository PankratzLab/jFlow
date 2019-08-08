package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.seq.manage.PileupProducer;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BamPile;
import org.genvisis.seq.manage.BamPile.Base;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.seq.qc.FilterNGS;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomeBuild;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

import htsjdk.variant.vcf.VCFFileReader;

public class CRAMSnpReader {

  private static final int MAX_CRAMS_PROCESSING = 3;// number of crams to read at a time
  public static final String CRAM_READS_EXT = ".reads.bed";
  public static final String HEADER = "#CHR\tSTART\tSTOP\tUCSC\tNAME\tREF\tALT\tNREF\tNALT\tNOTH\tGQ";

  private FilterNGS filterNGS = new FilterNGS(20, 20, null);
  private Logger log = new Logger();
  private ReferenceGenome refGen;
  private ASSEMBLY_NAME aName;

  private String snpVCF;
  private String outDir;

  private List<Segment> segList;
  private Map<Segment, NGSBin> segMap;
  private String[] cramFiles;

  private class NGSBin {

    int binStart;
    int binStop;
    String snpID;
    Segment snpSeg;
    String ref;
    String alt;
  }

  public static void main(String[] args) {
    CLI cli = new CLI(CRAMSnpReader.class);

    cli.addArg("vcf", "Selected-SNP VCF File", true);
    cli.addArg("cram", "CRAM file directory");
    cli.addArg("list", "File with list of cram files to process");
    cli.addGroup("cram", "list");
    cli.addArg("out", "Output file directory", true);
    cli.addArg("build", "Genome build, one of " + ArrayUtils.toStr(GenomeBuild.values(), ", "),
               true);

    if (args.length == 0) {
      System.setProperty("samjdk.reference_fasta",
                         Resources.genome(GenomeBuild.HG38, new Logger()).getFASTA().get());
      System.setProperty("reference_fasta",
                         Resources.genome(GenomeBuild.HG38, new Logger()).getFASTA().get());
      new CRAMSnpReader("G:\\bamTesting\\snpSelection\\selected_topmed.vcf",
                        "G:\\bamTesting\\01cram\\", "G:\\bamTesting\\01cram\\", GenomeBuild.HG38,
                        true).run();
    } else {
      cli.parseWithExit(args);
      System.setProperty("samjdk.reference_fasta",
                         Resources.genome(GenomeBuild.valueOf(cli.get("build").toUpperCase()),
                                          new Logger())
                                  .getFASTA().get());
      System.setProperty("reference_fasta",
                         Resources.genome(GenomeBuild.valueOf(cli.get("build").toUpperCase()),
                                          new Logger())
                                  .getFASTA().get());
      if (cli.has("list")) {
        new CRAMSnpReader(cli.get("vcf"),
                          HashVec.loadFileToStringArray(cli.get("list"), false, null,
                                                        false),
                          cli.get("out"), GenomeBuild.valueOf(cli.get("build").toUpperCase()), true)
                                                                                                    .run();
      } else {
        new CRAMSnpReader(cli.get("vcf"), cli.get("cram"), cli.get("out"),
                          GenomeBuild.valueOf(cli.get("build").toUpperCase()), true).run();
      }
    }

  }

  public CRAMSnpReader(String snpVCF, String cramDir, String outDir, GenomeBuild build,
                       boolean addChr) {
    this(snpVCF, listCRAMs(cramDir), outDir, build, addChr);
  }

  public CRAMSnpReader(String snpVCF, String[] cramFiles, String outDir, GenomeBuild build,
                       boolean addChr) {
    this.snpVCF = snpVCF;
    this.cramFiles = cramFiles;
    this.outDir = outDir;
    refGen = new ReferenceGenome(Resources.genome(build, log).getFASTA().getAbsolute(), log);
    aName = addChr ? ASSEMBLY_NAME.HG19 : ASSEMBLY_NAME.GRCH37;
  }

  public void run() {
    readVCF();
    processCRAMs();
  }

  private void processCRAMs() {
    Segment[] pileSegs = segList.toArray(new Segment[segList.size()]);
    int execNum = Math.min(cramFiles.length, MAX_CRAMS_PROCESSING);
    int numThreads = Runtime.getRuntime().availableProcessors() / execNum;
    ExecutorService exec = Executors.newFixedThreadPool(execNum);

    for (String c : cramFiles) {
      exec.submit(new Runnable() {

        @Override
        public void run() {

          long t1 = System.nanoTime();
          log.reportTime("Processing .cram file: " + c);
          try {
            String outFile = outDir + ext.rootOf(c) + CRAM_READS_EXT;
            // PrintWriter writer = new PrintWriter(new BlockCompressedOutputStream(outFile));
            PrintWriter writer = Files.getAppropriateWriter(outFile);
            writer.println(HEADER);
            BamPile[] bamPiles = PileupProducer.processBamFile(c, refGen, pileSegs, filterNGS,
                                                               aName, numThreads, log);
            for (int s = 0; s < pileSegs.length; s++) {
              NGSBin bin = segMap.get(pileSegs[s]);
              writer.print(bin.snpSeg.getChromosomeUCSC());
              writer.print("\t");
              writer.print(bin.binStart - 1);
              writer.print("\t");
              writer.print(bin.binStop);
              writer.print("\t");
              writer.print(bin.snpSeg.getUCSClocation());
              writer.print("\t");
              writer.print(bin.snpID);
              writer.print("\t");
              writer.print(bin.ref);
              writer.print("\t");
              writer.print(bin.alt);
              writer.print("\t");
              int ref = bamPiles[s].getBaseCount(Base.valueOf(bin.ref));
              writer.print(ref);
              writer.print("\t");
              int alt = bamPiles[s].getBaseCount(Base.valueOf(bin.alt));
              writer.print(alt);
              writer.print("\t");

              int rmn = 0;
              for (Base b : Base.values()) {
                if (Base.valueOf(bin.ref) != b && Base.valueOf(bin.alt) != b) {
                  rmn += bamPiles[s].getBaseCount(b);
                }
              }
              writer.print(Integer.toString(rmn));

              writer.print("\t");

              float gq = (ref + alt) / ((float) (ref + alt + rmn));
              writer.println(gq);
            }
            writer.close();
            BedOps.verifyBedIndex(outFile, log);
          } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
          log.reportTime("Processed " + c + " + in " + ext.getTimeElapsedNanos(t1));
        }
      });
    }
    exec.shutdown();
    try {
      exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
  }

  private static String[] listCRAMs(String cramDir) {
    String[] crams = new File(cramDir).list(new FilenameFilter() {

      @Override
      public boolean accept(File arg0, String arg1) {
        return arg1.endsWith(".cram");
      }
    });
    for (int i = 0; i < crams.length; i++) {
      crams[i] = ext.verifyDirFormat(cramDir) + crams[i];
    }
    return crams;
  }

  private void readVCF() {
    segList = new ArrayList<>();
    segMap = new HashMap<>();
    VCFFileReader reader = new VCFFileReader(new File(snpVCF), true);
    reader.iterator().stream().forEach(vc -> {
      int binStart = vc.getAttributeAsInt("BINSTART", -1);
      int binStop = vc.getAttributeAsInt("BINSTOP", -1);
      String name = vc.getID();
      String contig = vc.getContig();
      int start = vc.getStart();
      Segment seg = new Segment(contig, start, start);
      NGSBin bin = new NGSBin();
      bin.binStart = binStart;
      bin.binStop = binStop;
      bin.snpID = name;
      bin.snpSeg = seg;
      bin.ref = vc.getReference().getBaseString();
      bin.alt = vc.getAltAlleleWithHighestAlleleCount().getBaseString();
      segMap.put(seg, bin);
      segList.add(seg);
    });
    reader.close();
    log.reportTime("Loaded " + segList.size() + " snps to look for.");
  }

}
