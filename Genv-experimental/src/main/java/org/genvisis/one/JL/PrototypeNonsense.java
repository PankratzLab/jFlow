package org.genvisis.one.JL;

import org.genvisis.seq.manage.ReferenceGenome;
import org.pankratzlab.common.Logger;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

/**
 * 
 *
 */
public class PrototypeNonsense {

  public static void main(String[] args) {
    VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile("output.vcf");
    builder.clearOptions();
    builder.setOption(Options.INDEX_ON_THE_FLY);
    builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);
    VCFHeader vcfHeader = new VCFHeader();
    SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome("refGenome.fasta",
                                                                      new Logger()).getIndexedFastaSequenceFile()
                                                                                   .getSequenceDictionary();

    builder.setReferenceDictionary(samSequenceDictionary);
    vcfHeader.setSequenceDictionary(samSequenceDictionary);
    VariantContextWriter writer = builder.build();
    writer.writeHeader(vcfHeader);

    for (int i = 0; i < args.length; i++) {/// your data
      VariantContextBuilder builderVc = new VariantContextBuilder();
      builderVc.chr("chr");
      builderVc.start(1);
      builderVc.start(3);
      // builder.alleles(null);
      // builder...args. ...

      // add your info to builder

      writer.add(builderVc.make());
    }
    writer.close();
  }

}
