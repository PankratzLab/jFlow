package org.genvisis.one.JL;

import java.io.File;
import java.util.Map;
import org.genvisis.common.Logger;
import org.genvisis.seq.analysis.mutMap.VCMut;
import org.genvisis.seq.analysis.mutMap.VCMut.PARSE_METHOD;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SNPEffMap {

  public static void main(String[] args) {
    VCFFileReader reader = new VCFFileReader(new File("/Volumes/Beta/data/Cushings/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz"));

    // HashSet<String> snpeffs = new HashSet<>();
    // int i = 0;

    Map<String, String> effmap = VCMut.getSNPEFFOntologyMap(new Logger());
    System.out.println(effmap.keySet());
    // System.exit(1);
    for (VariantContext vc : reader) {
      VCMut vcMut = VCMut.fromVariantContext(vc, PARSE_METHOD.SNP_EFF, effmap);
      if (vcMut.isHasFullInfo()) {
        System.out.println(vcMut.toString());
      }
      // i++;
      // String eff = VCOps.getAnnotationsFor(new String[] {
      // "SNPEFF_EFFECT" }, vc, ".")[0];
      // // SNPEFF_EFFECT
      // snpeffs.add(eff);
      // if (i % 100000 == 0) {
      // System.out.println(snpeffs.toString());
      // }

    }
    // System.out.println(ArrayUtils.toStr(snpeffs, "\n"));

  }

}
