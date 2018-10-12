package org.genvisis.one.JL.mtDNA;

import java.io.File;
import org.genvisis.seq.manage.VCOps;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.core.CLI;
import org.pankratzlab.shared.filesys.Segment;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class ReverseAnnotate {

  private static void run(String vcf, String regions) {
    VCFFileReader reader = new VCFFileReader(new File(vcf), false);
    String out = ext.parseDirectoryOfFile(vcf) + ext.rootOf(regions) + "annot.txt";
    Segment[] segs = Segment.loadRegions(regions, 0, 1, 2, true);
    StringBuilder builder = new StringBuilder("#CHR\tBP1\tBP2\tID");
    for (VariantContext vc : reader) {
      for (Segment seg : segs) {
        if (VCOps.getSegment(vc).overlaps(seg)) {
          String[] toGet = new String[] {"Uniprot_name", "OXPHOS_complex"};
          String[] annot = VCOps.getAnnotationsFor(toGet, vc, "NOTHING");
          if (!annot[0].equals(".")) {
            builder.append("\n" + seg.getChromosomeUCSC() + "\t" + seg.getStart() + "\t"
                           + seg.getStop() + "\t" + annot[0]);
            builder.append("\n" + seg.getChromosomeUCSC() + "\t" + seg.getStart() + "\t"
                           + seg.getStop() + "\t" + annot[1]);
          }
        }
      }
    }
    reader.close();
    Files.write(builder.toString(), out);
    // cat mtGenesannot.txt |sort|uniq|grep -v 16569 |grep -v "9207 ATP8" |grep -v "8572
    // ATP6">mtUniport.reg
  }

  public static void main(String[] args) {
    CLI c = new CLI(ReverseAnnotate.class);
    String vcf = "a.vcf";
    String regions = "a.regions";

    c.addArgWithDefault("vcf", "vcf", vcf);
    c.addArgWithDefault("regions", "regions", regions);

    c.parseWithExit(args);
    run(c.get("vcf"), c.get("regions"));
  }

}
