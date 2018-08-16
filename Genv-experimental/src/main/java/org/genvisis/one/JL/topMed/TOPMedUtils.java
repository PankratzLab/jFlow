package org.genvisis.one.JL.topMed;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.one.JL.topMed.TOPMedUtils.GeneImpact;
import org.genvisis.one.JL.topMed.TOPMedUtils.IMPACT;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * 
 *
 */
public class TOPMedUtils {

  static final String MAFS = "mafs";
  static final String ANN = "ANN";
  static final String ANN_DELIM = ",";
  static final String ANN_DELIM_SUB = "\\|";
  static final String ANN_BLANK = ".";

  static final int IMPACT_INDEX = 2;
  static final int GENE_INDEX = 3;

  enum IMPACT {
    MODIFIER, LOW, MODERATE, HIGH;
  }

  static class GeneImpact {

    String gene;
    IMPACT impact;

    /**
     * @param gene
     * @param impact
     */
    GeneImpact(String gene, IMPACT impact) {
      super();
      this.gene = gene;
      this.impact = impact;
    }
  }

  static List<GeneImpact> getAllGeneImpacts(VariantContext vc) {
    List<GeneImpact> geneImpacts = new ArrayList<>();
    String ann = vc.getAttributeAsString(TOPMedUtils.ANN, TOPMedUtils.ANN_BLANK);
    if (!TOPMedUtils.ANN_BLANK.equals(ann)) {
      String[] anns = ann.split(TOPMedUtils.ANN_DELIM);
      for (String a : anns) {
        String[] subAnn = a.split(TOPMedUtils.ANN_DELIM_SUB);
        IMPACT impact = IMPACT.valueOf(subAnn[TOPMedUtils.IMPACT_INDEX]);
        geneImpacts.add(new GeneImpact(subAnn[TOPMedUtils.GENE_INDEX], impact));
      }
    }
    return geneImpacts;
  }
}
