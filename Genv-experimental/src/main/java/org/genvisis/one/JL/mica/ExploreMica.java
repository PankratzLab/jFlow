package org.genvisis.one.JL.mica;

import java.io.File;
import org.genvisis.seq.manage.VCOps;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.Segment;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author Kitty
 */
public class ExploreMica {

  private ExploreMica() {

  }

  public static void main(String[] args) {
    String[] vcfFiles = new String[] {"/Volumes/Beta/data/MICA/data/for_pankratz/allLocal_noWater_noHap.test.norm.decomp.vcf"};

    for (String vcfFile : vcfFiles) {
      Logger log = new Logger(ext.parseDirectoryOfFile(vcfFile) + "log.log");
      Segment mSeg = new Segment("chr6:31,378,315-31,380,254");
      VCFFileReader reader = new VCFFileReader(new File(vcfFile), false);
      for (VariantContext vc : reader) {
        if (VCOps.getSegment(vc).overlaps(mSeg)) {
          if (vc.isIndel()) {
            log.report(VCOps.getSegment(vc).getUCSClocation() + "\t" + vc.getID() + "\t"
                       + vc.getReference().getBaseString() + "\t"
                       + vc.getAlternateAlleles().toString());
            log.report(vc.toStringWithoutGenotypes());
          }
        }
      }
      reader.close();
    }

  }
}
