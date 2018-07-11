package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.StringJoiner;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Positions;
import org.genvisis.filesys.Segment;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class AnnovarDB {

  public static void main(String[] args) {
    //    https://github.com/kkshaxqd/myperlscript/blob/master/compileAnnnovarIndex.pl
    //    perl compileAnnnovarIndex.pl /scratch.global/topmed/popAFs/hg38_topmed.txt 1000 >/scratch.global/topmed/popAFs/hg38_topmed.txt.idx
    HashSet<String> use = new HashSet<>();
    //    --expression TOPMed_freeze_5b.AC \
    //    --expression TOPMed_freeze_5b.AF \
    //    --expression TOPMed_freeze_5b.AF_EM_POP_DEF_African_Americans \
    //    --expression TOPMed_freeze_5b.AF_EM_POP_DEF_EAS \
    //    --expression TOPMed_freeze_5b.AF_EM_POP_DEF_Hispanics \
    //    --expression TOPMed_freeze_5b.AF_EM_POP_DEF_SAS \
    //    --expression TOPMed_freeze_5b.AF_EM_POP_DEF_Whites \
    //    
    //    --expression TOPMed_freeze_5b.N_Alleles_EM_POP_DEF_African_Americans \
    //    --expression TOPMed_freeze_5b.N_Alleles_EM_POP_DEF_EAS \
    //    --expression TOPMed_freeze_5b.N_Alleles_EM_POP_DEF_Hispanics \
    //    --expression TOPMed_freeze_5b.N_Alleles_EM_POP_DEF_SAS \
    //    --expression TOPMed_freeze_5b.N_Alleles_EM_POP_DEF_Whites \
    use.add("AC");
    use.add("AF");
    use.add("AF_EM_POP_DEF_African_Americans");
    use.add("AF_EM_POP_DEF_EAS");
    use.add("AF_EM_POP_DEF_Hispanics");
    use.add("AF_EM_POP_DEF_Whites");
    use.add("AF_EM_POP_DEF_SAS");
    use.add("N_Alleles_EM_POP_DEF_African_Americans");
    use.add("N_Alleles_EM_POP_DEF_EAS");
    use.add("N_Alleles_EM_POP_DEF_Hispanics");
    use.add("N_Alleles_EM_POP_DEF_Whites");
    use.add("N_Alleles_EM_POP_DEF_SAS");

    String vcf = "/scratch.global/topmed/popAFs/freeze.5b.pass_and_fail.gtonly.minDP10_pcs.withPops.v_Apr_13_18.TOPMed.commit.51f93c82bddd9a61a73865433b96bb103c78f1db.conf.0.95.Mahal.3.AF.concat.sort.vcf.gz";

    VCFFileReader reader = new VCFFileReader(new File(vcf));
    String out = "/scratch.global/topmed/popAFs/hg38_topmed.txt";
    StringJoiner header = new StringJoiner("\t");
    header.add("#Chr");
    header.add("Start");
    header.add("End");
    header.add("Ref");
    header.add("Alt");

    Collection<VCFInfoHeaderLine> infos = reader.getFileHeader().getInfoHeaderLines();
    List<String> ids = new ArrayList<>();
    for (VCFInfoHeaderLine info : infos) {
      if (use.contains(info.getID())) {
        ids.add(info.getID());
        header.add("TOPMed_freeze_5b." + info.getID());
      }
    }

    PrintWriter writer = Files.getAppropriateWriter(out);
    writer.println(header.toString());
    int num = 0;
    //    chr1  964576  964612

    //    Segment[] tests = new Segment[] {new Segment("chr1", 964570, 964615),
    //                                     new Segment("chr17", 48037710, 48037718),
    //                                     new Segment("chr18", 54224213, 54224215),};
    //    for (Segment seq : tests) {
    //      CloseableIterator<VariantContext> iter = reader.query(Positions.getChromosomeUCSC(seq.getChr(),
    //                                                                                        true),
    //                                                            seq.getStart(), seq.getStop());
    CloseableIterator<VariantContext> iter = reader.iterator();
    while (iter.hasNext()) {
      VariantContext vc = iter.next();
      num++;
      if (num % 1000000 == 0) {
        System.out.println(num + " variants scanned");
      }
      StringJoiner variant = new StringJoiner("\t");
      variant.add(Integer.toString(Positions.chromosomeNumber(vc.getContig())));
      int start = vc.getStart();
      int end = vc.getEnd();
      String ref = "";
      String alt = "";
      if (vc.isSimpleInsertion()) {
        ref = "-";
        String sub = vc.getAlternateAllele(0).getBaseString();
        alt = sub.substring(1, sub.length());
      } else if (vc.isSimpleDeletion()) {
        ref = vc.getReference().getBaseString();
        alt = "-";
        start = start + 1;
      } else {
        ref = vc.getReference().getBaseString();
        alt = vc.getAlternateAllele(0).getBaseString();
      }
      variant.add(Integer.toString(start));
      variant.add(Integer.toString(end));
      variant.add(ref);
      variant.add(alt);

      if (vc.getAlternateAlleles().size() != 1) {
        reader.close();
        throw new IllegalArgumentException("invalid alt alleles "
                                           + ArrayUtils.toStr(vc.getAlternateAlleles()));

      }

      for (String info : ids) {
        if (use.contains(info)) {
          if (vc.hasAttribute(info)) {
            variant.add(vc.getAttributeAsString(info, "."));
          } else {
            variant.add(".");
          }
        }
      }
      writer.println(variant.toString());
      //      }
    }
    writer.close();
    reader.close();
  }

}
