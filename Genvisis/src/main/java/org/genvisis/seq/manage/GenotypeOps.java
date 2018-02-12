package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.common.Logger;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;

public class GenotypeOps {

  public static List<Allele> getNoCall() {
    ArrayList<Allele> noCall = new ArrayList<Allele>();
    noCall.add(Allele.NO_CALL);
    noCall.add(Allele.NO_CALL);
    return noCall;
  }

  /**
   * @param annosToGet genotype annotations to extract
   * @param g {@link Genotype}
   * @param defaultValue missing/failed extractions will be set to this
   * @param log
   * @return
   */
  public static String[] getGenoAnnotationsFor(String[] annosToGet, Genotype g, String defaultValue,
                                               Logger log) {
    String[] annos = new String[annosToGet.length];
    for (int i = 0; i < annos.length; i++) {
      String tmp = defaultValue;
      if (g.hasAnyAttribute(annosToGet[i])) {
        try {
          tmp = g.getAnyAttribute(annosToGet[i]).toString();
        } catch (NullPointerException npe) {
          log.reportTimeWarning("Could not extract annotation " + annosToGet[i] + " for genotype "
                                + g.toString());
          log.reportException(npe);
        }
      }
      annos[i] = tmp;
    }
    return annos;
  }

  public static String[][] getGenoFormatKeys(String vcf, Logger log) {
    VCFFileReader reader = new VCFFileReader(new File(vcf), false);
    String[] annotationKeys = new String[reader.getFileHeader().getFormatHeaderLines().size()];
    String[] descriptions = new String[reader.getFileHeader().getFormatHeaderLines().size()];
    int index = 0;
    for (VCFFormatHeaderLine vcfHeaderLine : reader.getFileHeader().getFormatHeaderLines()) {
      annotationKeys[index] = vcfHeaderLine.getID();
      descriptions[index] = vcfHeaderLine.getDescription();
      index++;
    }
    reader.close();
    return new String[][] {annotationKeys, descriptions};

  }

}
