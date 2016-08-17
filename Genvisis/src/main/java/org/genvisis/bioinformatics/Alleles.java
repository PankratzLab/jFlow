package org.genvisis.bioinformatics;

import org.genvisis.common.ext;

public class Alleles {
  public static String[] getAlleleFreqForA1(String a1, String[] bothAllelesAndFreqForFirst) {
    int decimalPos, sigFigs;
    double freq;

    if (a1.charAt(0) == bothAllelesAndFreqForFirst[0].charAt(0)) {
      return bothAllelesAndFreqForFirst;
    }

    if (a1.charAt(0) == bothAllelesAndFreqForFirst[1].charAt(0)) {
      if (bothAllelesAndFreqForFirst[2].equals("NA")) {
        return new String[] {a1, bothAllelesAndFreqForFirst[0], "NA"};
      } else {
        decimalPos = bothAllelesAndFreqForFirst[2].indexOf(".");
        if (decimalPos == -1) {
          if (bothAllelesAndFreqForFirst[2].equals("1")
              || bothAllelesAndFreqForFirst[2].equals("2")) {
            sigFigs = 0;
          } else {
            System.err.println("Error - invalid frequency: " + bothAllelesAndFreqForFirst[2]);
            return null;
          }
        } else {
          sigFigs = bothAllelesAndFreqForFirst[2].length() - decimalPos - 1;
        }
        try {
          freq = 1 - Double.parseDouble(bothAllelesAndFreqForFirst[2]);
          if (freq < 0 || freq > 1) {
            throw new NumberFormatException();
          }
        } catch (NumberFormatException nfe) {
          System.err.println("Error - invalid frequency: " + bothAllelesAndFreqForFirst[2]);
          return null;
        }
        return new String[] {a1, bothAllelesAndFreqForFirst[0], ext.formDeci(freq, sigFigs, true)};
      }
    }

    System.err.println("Error - Allele 1 (" + a1
                       + ") does not match with either of the alleles provided ("
                       + bothAllelesAndFreqForFirst[0] + "/" + bothAllelesAndFreqForFirst[1] + ")");

    return null;
  }

  public static char[] getAllelesInOrder(char a1, char[] bothAlleles) {
    if (a1 == bothAlleles[0]) {
      return new char[] {a1, bothAlleles[1]};
    }

    if (a1 == bothAlleles[1]) {
      return new char[] {a1, bothAlleles[0]};
    }

    System.err.println("Error - Allele 1 (" + a1
                       + ") does not match with either of the alleles provided (" + bothAlleles[0]
                       + "/" + bothAlleles[1] + ")");

    return null;
  }
}
