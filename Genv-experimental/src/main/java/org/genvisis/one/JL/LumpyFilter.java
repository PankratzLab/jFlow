package org.genvisis.one.JL;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.Segment;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Test filtering of SV vcf using svtyper GQ of lumpy calls
 */
public class LumpyFilter {

  private LumpyFilter() {

  }

  /**
   * @param args
   */
  public static void main(String[] args) {

    String vcf = "/Volumes/Beta/data/aric_sra/SRAPipeline/private/testBams/lumpyTest/H_UM-Schiffman-129-SS-129lumpy.gt.sort.vcf";
    Logger log = new Logger(ext.parseDirectoryOfFile(vcf) + "filt.log");
    SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome("/Volumes/Beta/ref/all_sequences.fa",
                                                                      log).getIndexedFastaSequenceFile()
                                                                          .getSequenceDictionary();
    VCFOps.verifyIndex(vcf, new Logger());
    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    VariantContextWriter writer = VCFOps.initBuilder(VCFOps.getAppropriateRoot(vcf, false)
                                                     + ".filt.vcf", VCFOps.DEFUALT_WRITER_OPTIONS,
                                                     samSequenceDictionary)
                                        .build();
    // VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY,
    // log);
    VCFHeader newVCFHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(),
                                           reader.getFileHeader().getGenotypeSamples());
    newVCFHeader.setSequenceDictionary(samSequenceDictionary);

    HashMap<String, List<String>> typeLocs = new HashMap<>();
    writer.writeHeader(newVCFHeader);
    for (VariantContext vc : reader) {
      if (vc.getGenotype(0).getGQ() > 50) {
        writer.add(vc);
        Segment seg = VCOps.getSegment(vc);
        String type = VCOps.getAnnotationsFor(new String[] {"SVTYPE"}, vc, ".")[0];
        if (!typeLocs.containsKey(type)) {
          typeLocs.put(type, new ArrayList<>());
        }
        typeLocs.get(type).add(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop() + "\t"
                               + vc.getGenotype(0).getGQ());
      }
    }
    reader.close();
    writer.close();
    for (String type : typeLocs.keySet()) {
      Files.writeIterable(typeLocs.get(type),
                          VCFOps.getAppropriateRoot(vcf, false) + "_" + type + ".filt.bed");
    }
  }

}
