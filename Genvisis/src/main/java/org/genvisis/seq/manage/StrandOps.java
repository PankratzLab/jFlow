package org.genvisis.seq.manage;

import java.util.List;

import org.genvisis.bioinformatics.Sequence;
import org.genvisis.common.Array;
import org.genvisis.common.ext;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class StrandOps {

  public enum CONFIG {
                      STRAND_CONFIG_SAME_ORDER_SAME_STRAND, STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND, STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND, STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND, STRAND_CONFIG_DIFFERENT_ALLELES, STRAND_CONFIG_BOTH_NULL, STRAND_CONFIG_SPECIAL_CASE, STRAND_CONFIG_AMBIGOUS, STRAND_CONFIG_UNKNOWN;
  }

  public static final String[] VALID_ALLELES = {"A", "C", "G", "T", "I", "D"};

  public static final String[] NULL_ALLELES = {".", "-", "N", "NA", "0"};

  // confusing terminology, here flipped means opposite strand, and opposite means flipped allele
  public static CONFIG determineStrandConfig(String[] alleles, String[] referenceAlleles) {
    String[] flipped;
    boolean[] nullChecks;
    int index;

    if (ext.indexOfStr(alleles[0], VALID_ALLELES) >= 0
        && ext.indexOfStr(alleles[1], VALID_ALLELES) >= 0) {
      if (isAmbiguous(alleles)) {
        return CONFIG.STRAND_CONFIG_AMBIGOUS;
      }

      if (referenceAlleles[0] == null) {
        referenceAlleles[0] = alleles[0];
        referenceAlleles[1] = alleles[1];
        return CONFIG.STRAND_CONFIG_SAME_ORDER_SAME_STRAND;
      } else if (referenceAlleles[1] == null) {
        if (alleles[0].equals(referenceAlleles[0])) {
          referenceAlleles[1] = alleles[1];
          // return STRAND_CONFIG_SAME;
          return CONFIG.STRAND_CONFIG_SPECIAL_CASE;
        } else if (alleles[1].equals(referenceAlleles[0])) {
          referenceAlleles[1] = alleles[0];
          // return STRAND_CONFIG_OPPOSITE;
          return CONFIG.STRAND_CONFIG_SPECIAL_CASE;
        } else {
          flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
          if (flipped[0].equals(referenceAlleles[0])) {
            referenceAlleles[1] = flipped[1];
            // return STRAND_CONFIG_SAME_FLIPPED;
            return CONFIG.STRAND_CONFIG_SPECIAL_CASE;
          } else if (flipped[1].equals(referenceAlleles[0])) {
            referenceAlleles[1] = flipped[0];
            // return STRAND_CONFIG_OPPOSITE_FLIPPED;
            return CONFIG.STRAND_CONFIG_SPECIAL_CASE;
          } else {
            return CONFIG.STRAND_CONFIG_DIFFERENT_ALLELES;
          }
        }
      } else {
        if (alleles[0].equals(referenceAlleles[0]) && alleles[1].equals(referenceAlleles[1])) {
          return CONFIG.STRAND_CONFIG_SAME_ORDER_SAME_STRAND;
        } else if (alleles[0].equals(referenceAlleles[1])
                   && alleles[1].equals(referenceAlleles[0])) {
          return CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
        } else {
          flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
          if (flipped[0].equals(referenceAlleles[0]) && flipped[1].equals(referenceAlleles[1])) {
            return CONFIG.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND;
          } else if (flipped[0].equals(referenceAlleles[1])
                     && flipped[1].equals(referenceAlleles[0])) {
            return CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
          } else {
            return CONFIG.STRAND_CONFIG_DIFFERENT_ALLELES;
          }
        }
      }
    } else {
      nullChecks = new boolean[] {false, false};
      for (int i = 0; i < nullChecks.length; i++) {
        if (ext.indexOfStr(alleles[i], NULL_ALLELES) >= 0) {
          nullChecks[i] = true;
        } else if (ext.indexOfStr(alleles[i], VALID_ALLELES) == -1) {
          return CONFIG.STRAND_CONFIG_SPECIAL_CASE;
        }
      }
      if (Array.booleanArraySum(nullChecks) == 1) {
        index = nullChecks[0] ? 1 : 0;
        if (referenceAlleles[0] == null) {
          referenceAlleles[0] = alleles[index];
          return index == 0 ? CONFIG.STRAND_CONFIG_SAME_ORDER_SAME_STRAND
                            : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
        } else if (referenceAlleles[1] == null) {
          if (alleles[index].equals(referenceAlleles[0])) {
            return index == 0 ? CONFIG.STRAND_CONFIG_SAME_ORDER_SAME_STRAND
                              : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
          } else {
            flipped = new String[] {Sequence.flip(alleles[index])};
            if (flipped[0].equals(referenceAlleles[0])) {
              return index == 0 ? CONFIG.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND
                                : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
            } else {
              return CONFIG.STRAND_CONFIG_DIFFERENT_ALLELES;
            }
          }
        } else {
          if (alleles[index].equals(referenceAlleles[0])) {
            return index == 0 ? CONFIG.STRAND_CONFIG_SAME_ORDER_SAME_STRAND
                              : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
          } else if (alleles[index].equals(referenceAlleles[1])) {
            return index == 1 ? CONFIG.STRAND_CONFIG_SAME_ORDER_SAME_STRAND
                              : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
          } else {
            flipped = new String[] {Sequence.flip(alleles[index])};
            if (flipped[0].equals(referenceAlleles[0])) {
              return index == 0 ? CONFIG.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND
                                : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
            } else if (flipped[0].equals(referenceAlleles[1])) {
              return index == 1 ? CONFIG.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND
                                : CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
            } else {
              return CONFIG.STRAND_CONFIG_DIFFERENT_ALLELES;
            }
          }
        }
      } else if (Array.booleanArraySum(nullChecks) == 2) {
        return CONFIG.STRAND_CONFIG_BOTH_NULL;
      } else {
        return CONFIG.STRAND_CONFIG_DIFFERENT_ALLELES;
      }
    }
  }

  public static String flipIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles) {
    if (strand == Strand.NEGATIVE) {
      if (b.equals("A")) {
        return "T";
      } else if (b.equals("G")) {
        return "C";
      } else if (b.equals("C")) {
        return "G";
      } else if (b.equals("T")) {
        return "A";
      } else {
        if (!ignoreInvalidAlleles) {
          throw new IllegalArgumentException("Invalid base for strand flip " + b);
        }
        return b;
      }
    } else {
      return b;
    }
  }


  public static String[] flipIfNeeded(String[] b, Strand strand, boolean ignoreInvalidAlleles) {
    String[] flipped = new String[b.length];
    for (int i = 0; i < flipped.length; i++) {
      flipped[i] = flipIfNeeded(b[i], strand, ignoreInvalidAlleles);
    }
    return flipped;
  }

  public static String flipsIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles) {
    return flipsIfNeeded(b, strand, ignoreInvalidAlleles, false);
  }

  public static String flipsIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles,
                                     boolean reverse) {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < b.length(); i++) {
      sb.append(flipIfNeeded(b.charAt(i) + "", strand, ignoreInvalidAlleles));
    }
    return reverse ? sb.reverse().toString() : sb.toString();
  }

  public static boolean isAmbiguous(String[] alleles) {
    boolean ambiguous = false;
    String a1 = alleles[0];
    String a2 = alleles[1];
    String[][] ambiguousDefs = new String[][] {{"A", "T"}, {"T", "A"}, {"C", "G"}, {"G", "C"}};
    for (String[] ambiguousDef : ambiguousDefs) {
      if (a1.equals(ambiguousDef[0]) && a2.equals(ambiguousDef[1])) {
        ambiguous = true;
        break;
      }

    }
    return ambiguous;
  }

  public static boolean isAmbiguous(VariantContext vc) {
    boolean ambiguous = false;
    String ref = vc.getReference().getDisplayString();
    List<Allele> alt = vc.getAlternateAlleles();
    String[][] ambiguousDefs = new String[][] {{"A", "T"}, {"T", "A"}, {"C", "G"}, {"G", "C"}};
    for (Allele b : alt) {
      for (String[] ambiguousDef : ambiguousDefs) {
        if (ref.equals(ambiguousDef[0]) && b.getDisplayString().equals(ambiguousDef[1])) {
          ambiguous = true;
          break;
        }
      }
    }

    return ambiguous;
  }
}
