package org.genvisis.seq.manage;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.bioinformatics.Sequence;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class StrandOps {

  private static final Map<String, String> STRAND_FLIPS;

  static {
    ImmutableMap.Builder<String, String> strandFlipsBuilder = ImmutableMap.builder();
    strandFlipsBuilder.put("A", "T");
    strandFlipsBuilder.put("T", "A");
    strandFlipsBuilder.put("C", "G");
    strandFlipsBuilder.put("G", "C");
    STRAND_FLIPS = strandFlipsBuilder.build();
  }

  public static Allele flipIfNeeded(Allele allele, Strand strand) {
    if (strand == Strand.NEGATIVE && !allele.isSymbolic()) {
      return Allele.create(flipEach(allele.getDisplayString(), true).toString(),
                           allele.isReference());
    }
    return allele;
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
    StringBuilder sb;

    if (strand == Strand.NEGATIVE) {
      sb = flipEach(b, ignoreInvalidAlleles);
    } else {
      sb = new StringBuilder(b);
    }
    return reverse ? sb.reverse().toString() : sb.toString();
  }

  private static StringBuilder flipEach(String b, boolean ignoreInvalidAlleles) {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < b.length(); i++) {
      sb.append(flip(String.valueOf(b.charAt(i)), ignoreInvalidAlleles));
    }
    return sb;
  }

  public static String flipIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles) {
    if (strand == Strand.NEGATIVE) {
      return flip(b, ignoreInvalidAlleles);
    } else {
      return b;
    }
  }

  private static String flip(String b, boolean ignoreInvalidAlleles) {
    if (b.length() == 1) {
      return STRAND_FLIPS.get(b);
    }
    if (!ignoreInvalidAlleles) {
      throw new IllegalArgumentException("Invalid base for strand flip: " + b);
    }
    return b;
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

  public enum AlleleOrder {
    SAME, OPPOSITE, UNKNOWN, AMBIGUOUS;
  }

  public enum AlleleStrand {
    SAME, FLIPPED, UNKNOWN, AMBIGUOUS;
  }

  public enum AlleleMatch {
    REF1,
    REF1_FLIP,
    REF2,
    REF2_FLIP,
    NO_MATCH,
    AMBIGUOUS,
    UNKNOWN,
    UNKNOWN_REF1,
    UNKNOWN_REF2,
    INVALID,
    INVALID_REF;
  }

  public enum AlleleStatus {
    VALID, NULL, INVALID;
  }

  public static class AlleleConfig {

    public AlleleConfig(AlleleStatus a1, AlleleStatus a2, AlleleMatch m1, AlleleMatch m2) {
      this.allele1 = a1;
      this.allele2 = a2;
      this.match1 = m1;
      this.match2 = m2;

      valid = (m1 == AlleleMatch.REF1 && m2 == AlleleMatch.REF2)
              || (m1 == AlleleMatch.REF2 && m2 == AlleleMatch.REF1)
              || (m1 == AlleleMatch.REF1_FLIP && m2 == AlleleMatch.REF2_FLIP)
              || (m1 == AlleleMatch.REF2_FLIP && m2 == AlleleMatch.REF1_FLIP)
              || (m1 == AlleleMatch.REF1 && m2 == AlleleMatch.UNKNOWN)
              || (m1 == AlleleMatch.REF2 && m2 == AlleleMatch.UNKNOWN)
              || (m1 == AlleleMatch.REF1_FLIP && m2 == AlleleMatch.UNKNOWN)
              || (m1 == AlleleMatch.REF2_FLIP && m2 == AlleleMatch.UNKNOWN)
              || (m2 == AlleleMatch.REF1 && m1 == AlleleMatch.UNKNOWN)
              || (m2 == AlleleMatch.REF2 && m1 == AlleleMatch.UNKNOWN)
              || (m2 == AlleleMatch.REF1_FLIP && m1 == AlleleMatch.UNKNOWN)
              || (m2 == AlleleMatch.REF2_FLIP && m1 == AlleleMatch.UNKNOWN)
              || (m1 == AlleleMatch.AMBIGUOUS && m2 == AlleleMatch.AMBIGUOUS);
    }

    final AlleleStatus allele1;
    final AlleleStatus allele2;

    final AlleleMatch match1;
    final AlleleMatch match2;

    final boolean valid;

  }

  public enum CONFIG {
    SAME_ORDER_SAME_STRAND("Match", AlleleOrder.SAME, AlleleStrand.SAME),
    SAME_ORDER_FLIPPED_STRAND("Flipped Strand", AlleleOrder.SAME, AlleleStrand.FLIPPED),
    OPPOSITE_ORDER_SAME_STRAND("Opposite Order", AlleleOrder.OPPOSITE, AlleleStrand.SAME),
    OPPOSITE_ORDER_FLIPPED_STRAND("Opposite Order / Flipped Strand", AlleleOrder.OPPOSITE,
                                  AlleleStrand.FLIPPED),
    DIFFERENT_ALLELES("Different Alleles", AlleleOrder.UNKNOWN, AlleleStrand.UNKNOWN),
    BOTH_NULL("Both NULL", AlleleOrder.UNKNOWN, AlleleStrand.UNKNOWN),
    SPECIAL_CASE("Special Case", AlleleOrder.UNKNOWN, AlleleStrand.UNKNOWN),
    AMBIGUOUS("Ambiguous", AlleleOrder.UNKNOWN, AlleleStrand.UNKNOWN),
    UNKNOWN("Unknown", AlleleOrder.UNKNOWN, AlleleStrand.UNKNOWN),
    INVALID("Invalid", AlleleOrder.UNKNOWN, AlleleStrand.UNKNOWN);

    private final String desc;
    private final AlleleOrder alleleOrder;
    private final AlleleStrand alleleStrand;

    /**
     * @param desc
     * @param alleleOrder
     * @param alleleStrand
     */
    private CONFIG(String desc, AlleleOrder alleleOrder, AlleleStrand alleleStrand) {
      this.desc = desc;
      this.alleleOrder = alleleOrder;
      this.alleleStrand = alleleStrand;
    }

    public String getDescription() {
      return desc;
    }

    public AlleleOrder getAlleleOrder() {
      return alleleOrder;
    }

    public AlleleStrand getAlleleStrand() {
      return alleleStrand;
    }

  }

  public static boolean isAmbiguous(String[] alleles) {
    String a1 = alleles[0];
    String a2 = alleles[1];
    if (AMBIG_MAP.containsKey(a1) && AMBIG_MAP.get(a1).equals(a2)) {
      return true;
    }
    return false;
  }

  private static final Map<String, String> AMBIG_MAP = ImmutableMap.of("A", "T", "T", "A", "C", "G",
                                                                       "G", "C");
  private static final Set<String> NULL_ALLELE_SET = ImmutableSet.of(".", "-", "N", "NA", "0");

  private static AlleleStatus getAlleleStatus(String allele) {
    if (NULL_ALLELE_SET.contains(allele)) {
      return AlleleStatus.NULL;
    } else if (Sequence.validAllele(allele)) {
      return AlleleStatus.VALID;
    } else {
      return AlleleStatus.INVALID;
    }
  }

  public static AlleleConfig strandConfig(String[] alleles, String[] referenceAlleles) {
    if (alleles.length != 2 || referenceAlleles.length != 2) {
      throw new IllegalArgumentException("Expected an array of length 2, got {alleles: "
                                         + alleles.length + "; ref: " + referenceAlleles.length
                                         + "}");
    }

    AlleleInternal a1 = new AlleleInternal(alleles[0]);
    AlleleInternal a2 = new AlleleInternal(alleles[1]);
    AlleleInternal r1 = new AlleleInternal(referenceAlleles[0]);
    AlleleInternal r2 = new AlleleInternal(referenceAlleles[1]);

    if (r1.status == AlleleStatus.INVALID || r2.status == AlleleStatus.INVALID
        || r1.allele.equals(r2.allele)) {
      return new AlleleConfig(a1.status, a2.status, AlleleMatch.INVALID_REF,
                              AlleleMatch.INVALID_REF);
    }
    AlleleMatch aM1, aM2;

    CHECK c11 = check(a1, r1);
    CHECK c12 = check(a1, r2);
    CHECK c21 = check(a2, r1);
    CHECK c22 = check(a2, r2);

    aM1 = match(c11, c12);
    aM2 = match(c21, c22);

    return new AlleleConfig(a1.status, a2.status, aM1, aM2);
  }

  private static enum CHECK {
    MATCH, FLIP, NO_MATCH, UNKNOWN
  }

  private static class AlleleInternal {

    private final String allele;
    private final AlleleStatus status;

    public AlleleInternal(String a) {
      this.status = getAlleleStatus(a);
      if (status == AlleleStatus.NULL) {
        this.allele = "-";
      } else {
        this.allele = a.toUpperCase();
      }
    }
  }

  private static AlleleMatch match(CHECK c1, CHECK c2) {
    if ((c1 == CHECK.MATCH && c2 == CHECK.FLIP) || (c1 == CHECK.FLIP && c2 == CHECK.MATCH)) {
      return AlleleMatch.AMBIGUOUS;
    }
    if (c1 == CHECK.MATCH && c2 != CHECK.MATCH) {
      return AlleleMatch.REF1;
    }
    if (c1 != CHECK.MATCH && c2 == CHECK.MATCH) {
      return AlleleMatch.REF2;
    }
    if (c1 == CHECK.FLIP && c2 != CHECK.FLIP) {
      return AlleleMatch.REF1_FLIP;
    }
    if (c1 != CHECK.FLIP && c2 == CHECK.FLIP) {
      return AlleleMatch.REF2_FLIP;
    }
    if ((c1 == CHECK.MATCH && c2 == CHECK.MATCH) || (c1 == CHECK.FLIP && c2 == CHECK.FLIP)) {
      return AlleleMatch.INVALID_REF;
    }
    if (c1 == CHECK.NO_MATCH && c2 == CHECK.NO_MATCH) {
      return AlleleMatch.NO_MATCH;
    }
    if (c1 == CHECK.UNKNOWN && c2 != CHECK.UNKNOWN) {
      return AlleleMatch.UNKNOWN_REF1;
    }
    if (c1 != CHECK.UNKNOWN && c2 == CHECK.UNKNOWN) {
      return AlleleMatch.UNKNOWN_REF2;
    }
    if (c1 == CHECK.UNKNOWN && c2 == CHECK.UNKNOWN) {
      return AlleleMatch.UNKNOWN;
    }
    throw new IllegalStateException("Invalid check conditions");
  }

  private static CHECK check(AlleleInternal all, AlleleInternal ref) {
    if (all.allele.length() != ref.allele.length()) {
      return CHECK.NO_MATCH;
    }
    if (all.status == AlleleStatus.VALID) {
      if (all.allele.length() == 1) {
        if (all.allele.equals(ref.allele)) {
          return CHECK.MATCH;
        } else if (Sequence.flip(all.allele).equals(ref.allele)) {
          return CHECK.FLIP;
        } else {
          if (ref.status == AlleleStatus.NULL) {
            return CHECK.UNKNOWN;
          } else {
            return CHECK.NO_MATCH;
          }
        }
      }
      CHECK stat = check(new AlleleInternal(all.allele.substring(0, 1)),
                         new AlleleInternal(ref.allele.substring(0, 1)));
      if (stat == CHECK.NO_MATCH) {
        return stat;
      }
      CHECK type = stat;
      for (int i = 1; i < all.allele.length(); i++) {
        CHECK newStat = check(new AlleleInternal(all.allele.substring(i, i + 1)),
                              new AlleleInternal(ref.allele.substring(i, i + 1)));
        if (stat == CHECK.MATCH && newStat == CHECK.FLIP) {
          return CHECK.NO_MATCH;
        } else if (stat == CHECK.FLIP && newStat == CHECK.MATCH) {
          return CHECK.NO_MATCH;
        } else if (newStat == CHECK.NO_MATCH) {
          return CHECK.NO_MATCH;
        } else if (newStat == CHECK.UNKNOWN) {
          stat = newStat;
        } else if (stat == CHECK.UNKNOWN) {
          if (type == CHECK.UNKNOWN) {
            type = newStat;
          } else if (type != newStat) {
            return CHECK.NO_MATCH;
          }
        }
      }
      return stat;
    }
    return CHECK.UNKNOWN;
  }

  // confusing terminology, here flipped means opposite strand, and opposite means flipped allele
  public static CONFIG determineStrandConfig(String[] alleles, String[] referenceAlleles) {
    String[] flipped;
    boolean[] nullChecks;
    int index;

    if (alleles.length != 2 || referenceAlleles.length != 2) {
      throw new IllegalArgumentException("Expected an array of length 2, got {alleles: "
                                         + alleles.length + "; ref: " + referenceAlleles.length
                                         + "}");
    }

    boolean[] validAlleles = {Sequence.validAllele(alleles[0]), Sequence.validAllele(alleles[1])};

    if (validAlleles[0] && validAlleles[1]) {
      if (isAmbiguous(alleles)) {
        return CONFIG.AMBIGUOUS;
      }

      if (referenceAlleles[0] == null) {
        referenceAlleles[0] = alleles[0];
        referenceAlleles[1] = alleles[1];
        return CONFIG.SAME_ORDER_SAME_STRAND;
      } else if (referenceAlleles[1] == null) {
        if (alleles[0].equals(referenceAlleles[0])) {
          referenceAlleles[1] = alleles[1];
          // return STRAND_CONFIG_SAME;
          return CONFIG.SPECIAL_CASE;
        } else if (alleles[1].equals(referenceAlleles[0])) {
          referenceAlleles[1] = alleles[0];
          // return STRAND_CONFIG_OPPOSITE;
          return CONFIG.SPECIAL_CASE;
        } else {
          flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
          if (flipped[0].equals(referenceAlleles[0])) {
            referenceAlleles[1] = flipped[1];
            // return STRAND_CONFIG_SAME_FLIPPED;
            return CONFIG.SPECIAL_CASE;
          } else if (flipped[1].equals(referenceAlleles[0])) {
            referenceAlleles[1] = flipped[0];
            // return STRAND_CONFIG_OPPOSITE_FLIPPED;
            return CONFIG.SPECIAL_CASE;
          } else {
            return CONFIG.DIFFERENT_ALLELES;
          }
        }
      } else {
        if (alleles[0].equals(referenceAlleles[0]) && alleles[1].equals(referenceAlleles[1])) {
          return CONFIG.SAME_ORDER_SAME_STRAND;
        } else if (alleles[0].equals(referenceAlleles[1])
                   && alleles[1].equals(referenceAlleles[0])) {
          return CONFIG.OPPOSITE_ORDER_SAME_STRAND;
        } else {
          flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
          if (flipped[0].equals(referenceAlleles[0]) && flipped[1].equals(referenceAlleles[1])) {
            return CONFIG.SAME_ORDER_FLIPPED_STRAND;
          } else if (flipped[0].equals(referenceAlleles[1])
                     && flipped[1].equals(referenceAlleles[0])) {
            return CONFIG.OPPOSITE_ORDER_FLIPPED_STRAND;
          } else {
            return CONFIG.DIFFERENT_ALLELES;
          }
        }
      }
    } else {
      nullChecks = new boolean[] {false, false};
      for (int i = 0; i < nullChecks.length; i++) {
        if (alleles[i] == null || NULL_ALLELE_SET.contains(alleles[i])) {
          nullChecks[i] = true;
        } else if (!validAlleles[i]) {
          return CONFIG.SPECIAL_CASE;
        }
      }
      if (nullChecks[0] && nullChecks[1]) {
        return CONFIG.BOTH_NULL;
      }
      if (ArrayUtils.booleanArraySum(nullChecks) == 1) {
        index = nullChecks[0] ? 1 : 0;
        if (referenceAlleles[0] == null) {
          referenceAlleles[0] = alleles[index];
          return index == 0 ? CONFIG.SAME_ORDER_SAME_STRAND : CONFIG.OPPOSITE_ORDER_SAME_STRAND;
        } else if (referenceAlleles[1] == null) {
          if (alleles[index].equals(referenceAlleles[0])) {
            return index == 0 ? CONFIG.SAME_ORDER_SAME_STRAND : CONFIG.OPPOSITE_ORDER_SAME_STRAND;
          } else {
            flipped = new String[] {Sequence.flip(alleles[index])};
            if (flipped[0].equals(referenceAlleles[0])) {
              return index == 0 ? CONFIG.SAME_ORDER_FLIPPED_STRAND
                                : CONFIG.OPPOSITE_ORDER_FLIPPED_STRAND;
            } else {
              return CONFIG.DIFFERENT_ALLELES;
            }
          }
        } else {
          if (alleles[index].equals(referenceAlleles[0])) {
            return index == 0 ? CONFIG.SAME_ORDER_SAME_STRAND : CONFIG.OPPOSITE_ORDER_SAME_STRAND;
          } else if (alleles[index].equals(referenceAlleles[1])) {
            return index == 1 ? CONFIG.SAME_ORDER_SAME_STRAND : CONFIG.OPPOSITE_ORDER_SAME_STRAND;
          } else {
            flipped = new String[] {Sequence.flip(alleles[index])};
            if (flipped[0].equals(referenceAlleles[0])) {
              return index == 0 ? CONFIG.SAME_ORDER_FLIPPED_STRAND
                                : CONFIG.OPPOSITE_ORDER_FLIPPED_STRAND;
            } else if (flipped[0].equals(referenceAlleles[1])) {
              return index == 1 ? CONFIG.SAME_ORDER_FLIPPED_STRAND
                                : CONFIG.OPPOSITE_ORDER_FLIPPED_STRAND;
            } else {
              return CONFIG.DIFFERENT_ALLELES;
            }
          }
        }
      } else {
        return CONFIG.DIFFERENT_ALLELES;
      }
    }
  }
}
