/**
 * 
 */
package org.genvisis.seq.manage;

/**
 * utility class for genomic sequence operations
 */
public class SeqOps {

  /**
   * Methods for how to compute GC content
   */
  public enum GC_COMP_METHOD {
    /**
     * GC content only computed using non-missing bases
     */
    GCTA_ONLY,
    /**
     * This is for compatibility for TelSeq, which includes "N"s in the total for GC content calcs
     */
    N_COUNTS_FOR_TOTAL
  }

  /**
   * @param seq sequence
   * @param gMethod {@link GC_COMP_METHOD} to use
   * @return proportion of sequence with G or C
   */
  public static double getProportionGC(String[] seq, GC_COMP_METHOD gMethod) {
    int gs = 0;
    int cs = 0;
    int as = 0;
    int ts = 0;
    int ns = 0;
    for (int i = 0; i < seq.length; i++) {
      if ("G".equalsIgnoreCase(seq[i])) {
        gs++;
      } else if ("C".equalsIgnoreCase(seq[i])) {
        cs++;
      } else if ("A".equalsIgnoreCase(seq[i])) {
        as++;
      } else if ("T".equalsIgnoreCase(seq[i])) {
        ts++;
      } else if ("N".equalsIgnoreCase(seq[i]) && gMethod == GC_COMP_METHOD.N_COUNTS_FOR_TOTAL) {
        ns++;
      } else if (!"N".equalsIgnoreCase(seq[i]) && !"".equalsIgnoreCase(seq[i].trim())) {
        throw new IllegalArgumentException("Invalid base " + seq[i]);
      }
    }
    int gsCs = gs + cs;
    int asTs = as + ts;
    int total = asTs + gsCs + ns;
    return (double) gsCs / total;
  }

  /**
   * @param seq count the number of occurrences of motif in this sequence
   * @param motif
   * @param caseSensitive If false, motif and sequence will be converted to lower case
   * @return number of times this motif is seen
   */
  public static int countMotif(String seq, String motif, boolean caseSensitive) {
    int lastIndex = 0;
    int count = 0;
    String seql = caseSensitive ? seq : seq.toLowerCase();
    String motifl = caseSensitive ? motif : motif.toLowerCase();
    while (lastIndex != -1) {
      lastIndex = seql.indexOf(motifl, lastIndex);
      if (lastIndex != -1) {
        count++;
        lastIndex += motifl.length();
      }
    }
    return count;
  }

}
