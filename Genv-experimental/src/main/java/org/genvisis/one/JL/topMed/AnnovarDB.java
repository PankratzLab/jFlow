package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.StringJoiner;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class AnnovarDB {

  public static void main(String[] args) {
    //    https://github.com/kkshaxqd/myperlscript/blob/master/compileAnnnovarIndex.pl
    //    perl compileAnnnovarIndex.pl /scratch.global/topmed/popAFs/hg38_topmed.txt 1000 >/scratch.global/topmed/popAFs/hg38_topmed.txt.idx

    String vcf = "/scratch.global/topmed/popAFs/freeze.5b.pass_and_fail.gtonly.minDP10_pcs.withPops.v_Apr_13_18.TOPMed.commit.51f93c82bddd9a61a73865433b96bb103c78f1db.conf.0.95.Mahal.3.AF.concat.sort.vcf.gz";

    VCFFileReader reader = new VCFFileReader(new File(vcf));
    String out = "/scratch.global/topmed/popAFs/hg38_topmed.txt";
    StringJoiner header = new StringJoiner("\t");
    header.add("#Chr");
    header.add("Start");
    header.add("End");
    header.add("Alt");

    Collection<VCFInfoHeaderLine> infos = reader.getFileHeader().getInfoHeaderLines();
    for (VCFInfoHeaderLine info : infos) {
      header.add("TOPMed_freeze_5b." + info.getID());
    }

    PrintWriter writer = Files.getAppropriateWriter(out);
    writer.println(header.toString());
    int num = 0;
    for (VariantContext vc : reader) {
      num++;
      if (num % 100000 == 0) {
        System.out.println(num + " variants scanned");
      }
      StringJoiner variant = new StringJoiner("\t");
      variant.add(vc.getContig());
      variant.add(Integer.toString(vc.getStart()));
      variant.add(Integer.toString(vc.getEnd()));
      if (vc.getAlternateAlleles().size() != 1) {
        throw new IllegalArgumentException("invalid alt alleles "
                                           + ArrayUtils.toStr(vc.getAlternateAlleles()));
      }
      variant.add(Integer.toString(vc.getEnd()));
      variant.add(vc.getAlternateAllele(0).getBaseString());
      for (VCFInfoHeaderLine info : infos) {
        if (vc.hasAttribute(info.getID())) {
          variant.add(vc.getAttributeAsString(info.getID(), "."));
        } else {
          variant.add(".");
        }
      }
      writer.println(variant.toString());
    }
    writer.close();
    reader.close();
  }

}
