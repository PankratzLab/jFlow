package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.seq.manage.PileupProducer;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BamPile;
import org.genvisis.seq.qc.FilterNGS;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;
import htsjdk.variant.vcf.VCFFileReader;

public class CRAMSnpReader {

  private FilterNGS filterNGS = new FilterNGS(20, 20, null);
  private Logger log = new Logger();
  private ReferenceGenome refGen;
  private ASSEMBLY_NAME aName;

  private String snpVCF;
  private String cramDir;
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
    System.setProperty("samjdk.reference_fasta",
                       Resources.genome(GenomeBuild.HG38, new Logger()).getFASTA().get());
    System.setProperty("reference_fasta",
                       Resources.genome(GenomeBuild.HG38, new Logger()).getFASTA().get());
    new CRAMSnpReader("G:\\bamTesting\\snpSelection\\selected_topmed.vcf",
                      "G:\\bamTesting\\00cram\\", "G:\\bamTesting\\00cram\\", GenomeBuild.HG38,
                      true).run();
  }

  public CRAMSnpReader(String snpVCF, String cramDir, String outDir, GenomeBuild build,
                       boolean addChr) {
    this.snpVCF = snpVCF;
    this.cramDir = cramDir;
    this.outDir = outDir;
    refGen = new ReferenceGenome(Resources.genome(build, log).getFASTA().getAbsolute(), log);
    aName = addChr ? ASSEMBLY_NAME.HG19 : ASSEMBLY_NAME.GRCH37;
  }

  public void run() {
    readVCF();
    listCRAMs();
    processCRAMs();
  }

  private void processCRAMs() {
    Segment[] pileSegs = segList.toArray(new Segment[segList.size()]);

    for (String c : cramFiles) {
      long t1 = System.nanoTime();
      log.reportTime("Processing .cram file: " + c);
      try {
        PrintWriter writer = Files.getAppropriateWriter(outDir + ext.rootOf(c) + ".reads.bed.gz");
        writer.println("#CHR\tSTART\tSTOP\tUCSC\tSNP\tREF\tALT\tA\tG\tC\tT\tN");
        BamPile[] bamPiles = PileupProducer.processBamFile(cramDir + c, refGen, pileSegs, filterNGS,
                                                           aName, log);
        for (int s = 0; s < pileSegs.length; s++) {
          NGSBin bin = segMap.get(pileSegs[s]);
          writer.print(bin.snpSeg.getChromosomeUCSC());
          writer.print("\t");
          writer.print(bin.binStart);
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
          writer.println(ArrayUtils.toStr(bamPiles[s].getAlleleCounts(log), "\t"));
        }
        writer.close();
      } catch (Exception e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      log.reportTime("Processed " + c + " + in " + ext.getTimeElapsedNanos(t1));
    }
  }

  private void listCRAMs() {
    cramFiles = new File(cramDir).list(new FilenameFilter() {

      @Override
      public boolean accept(File arg0, String arg1) {
        return arg1.endsWith(".cram");
      }
    });
    log.reportTime("Found " + cramFiles.length + " .cram files to process.");
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
