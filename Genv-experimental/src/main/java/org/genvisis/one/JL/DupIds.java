package org.genvisis.one.JL;

import java.io.File;
import java.util.HashSet;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;
import org.genvisis.seq.manage.VCFOps.PLINK_SET_MODE;
import org.genvisis.seq.manage.VCOps;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Would not recommend even looking at this code
 */
public class DupIds {

  private DupIds() {

  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    String vcf = "/Volumes/Beta/data/Cushings/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g_ids.vcf.gz";
    String varSet = ".variant";
    String outDir = "/Volumes/Beta/data/Cushings/dupDetect/";
    Logger log = new Logger(outDir + "log.log");
    double maf = 0.05;
    String outVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + varSet + "." + maf + ".vcf.gz";

    new File(outDir).mkdirs();

    if (!Files.exists(outVCF)) {
      String[] samples = VCFOps.getSamplesInFile(vcf);
      HashSet<String> sampToUse = new HashSet<>();
      for (String samp : samples) {
        if (samp.endsWith(varSet)) {
          sampToUse.add(samp);
        }
      }
      VCFFileReader reader = new VCFFileReader(new File(vcf), true);
      VariantContextWriter writer = VCFOps.initWriter(outVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                      reader.getFileHeader()
                                                            .getSequenceDictionary());

      VCFOps.copyHeader(reader, writer, sampToUse, HEADER_COPY_TYPE.SUBSET_STRICT, log);

      int num = 0;
      int written = 0;
      for (VariantContext vc : reader) {
        num++;
        try {
          double max = Double.parseDouble(VCOps.getAnnotationsFor(new String[] {"PopFreqMax"}, vc,
                                                                  ".")[0]);
          if (!vc.isFiltered() && !vc.isIndel() && vc.isBiallelic() && max > maf) {
            written++;
            writer.add(VCOps.getSubset(vc, sampToUse));
          }
          if (num % 10000 == 0) {
            log.reportTimeInfo("num=" + num + " written = " + written);
          }
        } catch (NumberFormatException e) {}
      }
      reader.close();
      writer.close();
    }
    VCFOps.convertToPlinkSet("", outVCF, "plink", PLINK_SET_MODE.GWAS_QC, log);
  }
}
