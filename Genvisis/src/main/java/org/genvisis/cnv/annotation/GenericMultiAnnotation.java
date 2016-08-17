package org.genvisis.cnv.annotation;

import java.util.Hashtable;

import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author lane0212 loading for existing annotations (SNP_EFF,etc) that do not require object
 *         parsing
 */
public class GenericMultiAnnotation implements AnnotationParser {

  /**
   * Note that these hard-coded BASICs may or may not be present in our anno file and will become a
   * "."
   */
  public static final String[] SNP_EFF_BASIC =
      new String[] {"SNPEFF_GENE_NAME", "SNPEFF_IMPACT", "SNPEFF_TRANSCRIPT_ID"};
  public static final String[] FREQ_BASIC =
      new String[] {"popfreq_max_20150413", "popfreq_all_20150413", "esp6500siv2_all",
                    "esp6500siv2_aa", "esp6500siv2_ea"};

  private final Hashtable<String, Integer> index;
  private final String[] locusNames;
  private final String[] keysToLoad;
  private String[][] data;
  private boolean found;// Warning, is true if any of the locusNames are found

  /**
   * @param locusNames will typically be the marker names of interest
   * @param keysToLoad annotation keys to populate
   */
  public GenericMultiAnnotation(String[] locusNames, String[] keysToLoad) {
    super();
    this.locusNames = locusNames;
    this.keysToLoad = keysToLoad;
    found = false;
    index = new Hashtable<String, Integer>();
    for (int i = 0; i < locusNames.length; i++) {
      index.put(locusNames[i], i);
    }
  }

  public String[] getLocusNames() {
    return locusNames;
  }

  @Override
  public boolean isFound() {
    return found;
  }

  @Override
  public void parseAnnotation(VariantContext vc, Logger log) {
    data[index.get(vc.getID())] = VCOps.getAnnotationsFor(keysToLoad, vc, ".");
  }

  @Override
  public void setFound(boolean found) {
    this.found = found;
  }

  @Override
  public boolean shouldAnnotateWith(VariantContext vc, Logger log) {
    return index.containsKey(vc.getID());
  }

}
