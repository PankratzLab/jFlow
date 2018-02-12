/**
 * 
 */
package org.genvisis.cnv.analysis.pod;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import org.genvisis.cnv.analysis.pod.PODAnalysis.GenoCompResult;

/**
 * Implementing logic described in table 1 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3680018/
 */

class PODGenotype {

  enum Genotype {

    AA(2, 0), AB(1, 1), BB(0, 2), NC(0, 0);

    private int numA;
    private int numB;

    private Genotype(int numA, int numB) {
      this.numA = numA;
      this.numB = numB;
    }

    int getNumA() {
      return numA;
    }

    int getNumB() {
      return numB;
    }

    static GenoCompResult getSharedAlleleCount(Genotype off, Genotype p) {
      int equals = off == p ? 1 : 0;
      int bothCalled = off != Genotype.NC && p != Genotype.NC ? 1 : 0;
      switch (p) {
        case AA:
          return new GenoCompResult(equals, off.getNumA(), bothCalled);
        case AB:
          return new GenoCompResult(equals, off.getNumA() + off.getNumB(), bothCalled);
        case BB:
          return new GenoCompResult(equals, off.getNumB(), bothCalled);
        case NC:
          return new GenoCompResult(0, 0, 0);
        default:
          throw new IllegalArgumentException("Invalid genotype");

      }
    }

    static Genotype fromByte(byte geno) {
      switch (geno) {
        case 0:
          return AA;
        case 1:
          return AB;
        case 2:
          return BB;
        case -1:
          return NC;

        default:
          throw new IllegalArgumentException("Invalid genotype " + geno);
      }
    }
  }

  enum BAF_EFFECT {
    INCREASE, DECREASE;
  }

  enum POD {
    MATERNAL, PATERNAL, NONE;
  }

  /**
   * Informative genotypes for parent 2 when parent 1 is {@link Genotype#BB}
   */
  private static final Set<Genotype> BBInformative = new HashSet<>(Arrays.asList(Genotype.AA,
                                                                                 Genotype.AB,
                                                                                 Genotype.NC));

  /**
   * Informative genotypes for parent 2 when parent 1 is {@link Genotype#AA}
   */
  private static final Set<Genotype> AAInformative = new HashSet<>(Arrays.asList(Genotype.BB,
                                                                                 Genotype.AB,
                                                                                 Genotype.NC));

  static POD getPodEffect(Genotype mother, Genotype father, BAF_EFFECT bEffect) {
    if ((mother == Genotype.NC && father == Genotype.NC) || mother == father) {
      // not informative
      return POD.NONE;
    }
    switch (bEffect) {
      case DECREASE:
        return getDecreasedBAFPODEffect(mother, father);
      case INCREASE:
        return getIncreasedBAFPODEffect(mother, father);
      default:
        throw new IllegalArgumentException("Invalid baf effect " + bEffect);

    }

  }

  private static POD getIncreasedBAFPODEffect(Genotype mother, Genotype father) {

    if (mother == Genotype.AA && AAInformative.contains(father)) {
      return POD.PATERNAL;
    } else if (father == Genotype.AA && AAInformative.contains(mother)) {
      return POD.MATERNAL;

    } else {
      return POD.NONE;
    }
  }

  private static POD getDecreasedBAFPODEffect(Genotype mother, Genotype father) {

    if (mother == Genotype.BB && BBInformative.contains(father)) {
      return POD.PATERNAL;
    } else if (father == Genotype.BB && BBInformative.contains(mother)) {
      return POD.MATERNAL;

    } else {
      return POD.NONE;
    }

  }

}
