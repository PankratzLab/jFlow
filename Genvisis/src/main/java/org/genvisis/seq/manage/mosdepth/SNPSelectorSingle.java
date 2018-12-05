package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.util.HashSet;
import java.util.Set;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.genvisis.cnv.Resources;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.manage.BEDFileReader;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class SNPSelectorSingle extends NGSBinSNPSelector {

  public SNPSelectorSingle(String in) {
    this.inputFile = in;
  }

  String inputFile;

  public void run() {
    setup();

    VariantContextWriter vcfWriter = openWriter(inputFile);

    long startN = System.nanoTime();
    VCFHeader header;

    try (VCFFileReader reader = new VCFFileReader(new File(inputFile),
                                                  Files.exists(inputFile + ".tbi"))) {
      header = reader.getFileHeader();

      // discover which chrs are in this vcf file
      BidiMap<String, Integer> contigMap = new DualHashBidiMap<>();
      boolean useChrPrepend = header.getContigLines().get(0).getID().startsWith("chr");
      header.getContigLines().forEach(vch -> {
        // ensure parsability
        contigMap.put(vch.getID(), (int) Positions.chromosomeNumber(vch.getID()));
      });

      for (Segment bin : bins.getLoci()) {
        // no variants in file for this bin's chromosome  
        if (!contigMap.containsValue((int) bin.getChr()) || !chrs.contains((int) bin.getChr())) {
          continue;
        }

        writeSelectedForBin(bin, useChrPrepend, contigMap, reader, vcfWriter);
      }

    }

    log.reportTime("Finished processing  " + inputFile + " in " + ext.getTimeElapsedNanos(startN)
                   + "!");
    vcfWriter.close();
  }

  public static void main(String[] args) {
    CLI cli = new CLI(SNPSelectorSingle.class);

    cli.addArg("out", "Output directory");
    cli.addArg("in", "Input VCF file");
    cli.addArg("mapFile", "Mapping quality file", false);
    cli.addArg("mq", "Mapping quality tag in the filters definitions file.", "MQ", false);
    cli.addArg("af", "Allele frequency tag in the VCF and filters definitions file", "AF", false);
    cli.addArg("chrs", "Comma-delimited list of chromosomes to process.", false);
    cli.addArg("bed", "BED file with predefined regions to use.", false);
    cli.addArg("bin", "Bin size to use when splitting up reference genome, default value 1000.",
               false);
    cli.addArg("filters", "File containing filter definitions.", false);
    cli.addGroup("bed", "bin");

    cli.parseWithExit(args);

    Logger log = new Logger();
    NGSBinSNPSelector selector = new SNPSelectorSingle(cli.get("in"));

    selector.outputFile = cli.get("out");
    selector.mappingQualityFile = cli.has("mapFile") ? cli.get("mapFile") : null;
    selector.afTag = cli.has("af") ? cli.get("af") : "AF";
    selector.mqTag = cli.has("mq") ? cli.get("mq") : "MQ";
    int bin = 1000;
    String bedFile = null;
    if (cli.has("bin")) {
      bin = cli.getI("bin");
    } else if (cli.has("bed")) {
      bedFile = cli.get("bed");
    } else {
      log.reportTime("No BED file specified, nor was a bin size specified, so the default bin size of "
                     + bin + " will be used.");
    }

    selector.filterDefsFile = cli.has("filters") ? cli.get("filters") : null;

    Set<Integer> chrs = null;
    if (cli.has("chrs")) {
      String[] v = cli.get("chrs").split(",");
      chrs = new HashSet<>();
      for (String vi : v) {
        chrs.add(Integer.parseInt(vi));
      }
    } else {
      chrs = Sets.newHashSet(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                             21, 22, 23, 24, 25);
    }
    selector.chrs = chrs;
    if (bedFile != null) {
      selector.bins = new BEDFileReader(bedFile, false).loadAll(log).getStrictSegmentSet();
    } else {
      selector.bins = new ReferenceGenome(Resources.genome(GenomeBuild.HG19, log).getFASTA()
                                                   .getAbsolute(),
                                          log).getBins(bin, chrs);
    }
    selector.run();
  }

}
