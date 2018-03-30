package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class AddAD {

  public static void main(String[] args) {

    VCFFormatHeaderLine formatAD = new VCFFormatHeaderLine("AD", 1, VCFHeaderLineType.String,
                                                           "MOCK AD");
    String vcf = "a.vcf.gz";
    String outVCF = VCFOps.getAppropriateRoot(vcf, false) + ".AD.vcf.gz";
    Logger log = new Logger(outVCF + ".log");
    VCFFileReader reader = new VCFFileReader(new File(vcf));
    VCFHeader header = reader.getFileHeader();
    VariantContextWriter writer = VCFOps.initWriter(outVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                    header.getSequenceDictionary());
    header.addMetaDataLine(formatAD);
    writer.writeHeader(header);
    int num = 0;
    for (VariantContext vc : reader) {
      num++;
      if (num % 1000 == 0) {
        log.reportTimeInfo(Integer.toString(num));
      }
      VariantContextBuilder builder = new VariantContextBuilder(vc);
      GenotypesContext gc = vc.getGenotypes();
      List<Genotype> newGenos = new ArrayList<>();
      for (Genotype g : gc) {
        newGenos.add(new GenotypeBuilder(g).AD(new int[] {50, 50}).make());
      }
      builder.genotypes(newGenos);
      writer.add(builder.make());
    }
    reader.close();
    writer.close();

  }

}
