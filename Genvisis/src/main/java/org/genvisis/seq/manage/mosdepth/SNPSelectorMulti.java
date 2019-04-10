package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
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

public class SNPSelectorMulti extends NGSBinSNPSelector {

  public SNPSelectorMulti(String inDir, String format, String token) {
    this.inputDir = inDir;
    this.format = format;
    this.token = token;
  }

  String inputDir;
  String format;
  String token;

  public void run() {
    setup();

    Map<Integer, String> allVCFs = new HashMap<>();

    for (int c : chrs) {
      String vcf = inputDir + format.replace(token, Integer.toString(c));
      if (!Files.exists(vcf)) {
        if (Files.exists(inputDir + format.replace(token, "0" + Integer.toString(c)))) {
          vcf = inputDir + format.replace(token, "0" + Integer.toString(c));
        } else if (Files.exists(inputDir
                                + format.replace(token,
                                                 "0" + Positions.chromosomeNumberInverse(c)))) {
          vcf = inputDir + format.replace(token, "0" + Positions.chromosomeNumberInverse(c));
        }
      }
      if (!Files.exists(vcf)) continue;
      allVCFs.put(c, vcf);
    }

    VariantContextWriter vcfWriter = openWriter(allVCFs.values().iterator().next());

    for (Entry<Integer, String> vcf : allVCFs.entrySet()) {
      long startN = System.nanoTime();
      VCFHeader header;

      try (VCFFileReader reader = new VCFFileReader(new File(vcf.getValue()),
                                                    Files.exists(vcf + ".tbi"))) {
        header = reader.getFileHeader();

        // discover which chrs are in this vcf file
        BidiMap<String, Integer> contigMap = new DualHashBidiMap<>();
        boolean useChrPrepend = header.getContigLines().get(0).getID().startsWith("chr");
        header.getContigLines().forEach(vch -> {
          // ensure parsability
          contigMap.put(vch.getID(), (int) Positions.chromosomeNumber(vch.getID()));
        });

        Set<Integer> found = new HashSet<>();
        for (Segment bin : bins.getLoci()) {
          // no variants in file for this bin's chromosome
          if (bin.getChr() != vcf.getKey().intValue()
              || !contigMap.containsValue((int) bin.getChr())
              || !chrs.contains((int) bin.getChr())) {
            // System.err.println("Error - bin contig " + bin.getChr()
            // + " was not present in the VCF!");
            continue;
          }
          found.add((int) bin.getChr());
          writeSelectedForBin(bin, useChrPrepend, contigMap, reader, vcfWriter);
        }
        log.reportTimeWarning("Found " + found.size() + " chrs (" + found.toString() + ") in vcf "
                              + vcf.getValue());
      }

      log.reportTime("Finished processing chr " + vcf.getKey() + " in "
                     + ext.getTimeElapsedNanos(startN) + "!");
      writeStats();
    }
    vcfWriter.close();
  }

  public static void main(String[] args) {
    CLI cli = new CLI(SNPSelectorMulti.class);

    cli.addArg("out", "Output directory");
    cli.addArg("in", "Input directory");
    cli.addArg("format", "File format with special token where the chromosome number goes");
    cli.addArg("token", "Special filename token to replace with chromosome number", "##", false);
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
    NGSBinSNPSelector selector = new SNPSelectorMulti(cli.get("in"), cli.get("format"),
                                                      cli.get("token"));

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
      selector.bins = new ReferenceGenome(Resources.genome(GenomeBuild.HG38, log).getFASTA()
                                                   .getAbsolute(),
                                          log).getBins(bin, chrs);
    }
    selector.run();
  }

}
