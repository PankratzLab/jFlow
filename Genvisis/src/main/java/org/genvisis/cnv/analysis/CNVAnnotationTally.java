/**
 * 
 */
package org.genvisis.cnv.analysis;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.annotation.segments.SegmentAnnotationKeys;
import org.genvisis.filesys.CNVariantAnnotated;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ListMultimap;

/**
 * Class for counting annotations from {@link CNVariantAnnotated}, hopefully
 */
public class CNVAnnotationTally {

  private ListMultimap<Integer, CNVariantAnnotated> locs;

  private HashMultimap<Integer, String> inds;

  /**
   * j
   */
  public CNVAnnotationTally() {
    super();
    this.locs = ArrayListMultimap.create();
    this.inds = HashMultimap.create();
  }

  /**
   * This method assumes that IIDs are unique
   */
  private void add(CNVariantAnnotated cnVariantAnnotated) {
    locs.put(cnVariantAnnotated.getCN(), cnVariantAnnotated);
    inds.put(cnVariantAnnotated.getCN(), cnVariantAnnotated.getIndividualID());

  }

  public ListMultimap<Integer, CNVariantAnnotated> getLocs() {
    return locs;
  }

  public HashMultimap<Integer, String> getInds() {
    return inds;
  }

  /**
   * @param cList list of {@link CNVariantAnnotated} to tally
   * @param key the annotation key to tally
   * @return count map holding map per CN seen
   */
  public static Map<String, CNVAnnotationTally> tallyAnnotation(List<CNVariantAnnotated> cList,
                                                                SegmentAnnotationKeys key) {
    Map<String, CNVAnnotationTally> tally = new HashMap<>();
    for (CNVariantAnnotated cnVarAn : cList) {
      if (cnVarAn.getSegmentAnotation().getAttributes().containsKey(key.toString())) {
        List<String> anns = cnVarAn.getSegmentAnotation().getAttributes().get(key.toString());
        for (String ann : anns) {
          prepTally(tally, ann);
          tally.get(ann).add(cnVarAn);
        }
      }
    }
    return tally;
  }

  private static void prepTally(Map<String, CNVAnnotationTally> tally, String ann) {
    if (!tally.containsKey(ann)) {
      tally.put(ann, new CNVAnnotationTally());
    }
  }

}
