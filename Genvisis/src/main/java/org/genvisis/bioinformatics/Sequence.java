package org.genvisis.bioinformatics;

import java.util.Map;
import java.util.Set;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;

public class Sequence {

  /*
   * @format:off
   */
  public static final char[] ALLELES = {'A', 'C', 'T', 'G'};
  private static final Set<Character> VALID_NUCLES = ImmutableSet.of('A', 'C', 'T', 'G', 'N', 'a', 'c','t', 'g', 'n');
  private static final Map<Character, Character> BASE_MAP = new ImmutableMap.Builder<Character, Character>()
                                                                            .put('A','T')
                                                                            .put('C', 'G')
                                                                            .put('G', 'C')
                                                                            .put('T', 'A')
                                                                            .put('I','I')
                                                                            .put('D','D')
                                                                            .put('a','t')
                                                                            .put('c','g')
                                                                            .put('g','c')
                                                                            .put('t','a')
                                                                            .put('i','i')
                                                                            .put('d','d')
                                                                            .build();
  /*
   * @format:on
   */

  public static char flip(char allele) {
    Character c = BASE_MAP.get(allele);
    return c == null ? allele : c;
  }

  public static String flip(String allele) {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < allele.length(); i++) {
      sb.append(flip(allele.charAt(i)));
    }
    return sb.toString();
  }

  public static boolean validBase(String allele) {
    return allele.length() == 1 && validBase(allele.charAt(0));
  }

  public static boolean validAllele(String allele) {
    if (allele == null) return false;
    if (allele.length() == 1) {
      return validBase(allele);
    }

    for (int i = 0; i < allele.length(); i++) {
      if (!VALID_NUCLES.contains(allele.charAt(i))) {
        return false;
      }
    }

    return true;
  }

  public static boolean validBase(char allele) {
    return BASE_MAP.keySet().contains(allele);
  }
}
