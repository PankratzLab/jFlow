/**
 * 
 */
package org.genvisis.cnv.analysis.pod;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.cnv.analysis.pod.PODGenotype.BAF_EFFECT;
import org.genvisis.cnv.analysis.pod.PODGenotype.Genotype;
import org.genvisis.cnv.analysis.pod.PODGenotype.POD;
import org.pankratzlab.common.ArrayUtils;

/**
 * Class to determine informative bafs
 */
class InformativeBAF {

  /**
   * Default sd cut
   */
  static final double CHEBYSHEV = Math.sqrt(2);

  enum BAF_STATUS {
    INFORMATIVE, NON_INFORMATIVE;
  }

  enum BAF_STRATEGY {
    HET_ONLY;
  }

  static class InformativeResult {

    private List<Integer> informatives;
    private List<BAF_EFFECT> effects;
    private NormalDistribution normalDistribution;
    private List<POD> pods;

    InformativeResult(List<Integer> informatives, List<BAF_EFFECT> effects,
                      NormalDistribution normalDistribution) {
      super();
      this.informatives = informatives;
      this.effects = effects;
      this.normalDistribution = normalDistribution;
      this.pods = new ArrayList<>();
    }

    List<Integer> getInformatives() {
      return informatives;
    }

    List<POD> getPods() {
      return pods;
    }

    void setPods(List<POD> pods) {
      this.pods = pods;
    }

    List<BAF_EFFECT> getEffects() {
      return effects;
    }

    NormalDistribution getNormalDistribution() {
      return normalDistribution;
    }
  }

  static InformativeResult getInformativeIndices(double[] bafs, byte[] genotypes,
                                                 BAF_STRATEGY bStrategy, double sdCut) {

    ArrayList<Integer> tmp = new ArrayList<>();

    for (int i = 0; i < genotypes.length; i++) {
      switch (bStrategy) {
        case HET_ONLY:
          if (Genotype.fromByte(genotypes[i]) == Genotype.AB) {
            tmp.add(i);
          }
          break;
        default:
          throw new IllegalArgumentException("Invalid BAF strategy " + bStrategy);
      }
    }

    NormalDistribution normalDistribution = getDistribution(ArrayUtils.subArray(bafs,
                                                                                tmp.stream()
                                                                                   .mapToInt(i -> i)
                                                                                   .toArray()));
    ArrayList<Integer> informatives = new ArrayList<>();
    ArrayList<BAF_EFFECT> effects = new ArrayList<>();

    double comp = sdCut * normalDistribution.getStandardDeviation();
    for (Integer i : tmp) {
      if (Genotype.fromByte(genotypes[i]) != Genotype.AB && bStrategy == BAF_STRATEGY.HET_ONLY) {
        throw new IllegalArgumentException("Invalid baf outlier test for non het genotype");
      }
      if (bafs[i] < normalDistribution.getMean() - comp) {
        informatives.add(i);
        effects.add(BAF_EFFECT.DECREASE);
      } else if (bafs[i] > normalDistribution.getMean() + comp) {
        informatives.add(i);
        effects.add(BAF_EFFECT.INCREASE);
      }
    }

    return new InformativeResult(informatives, effects, normalDistribution);

  }

  private static NormalDistribution getDistribution(double[] bafs) {
    return new NormalDistribution(ArrayUtils.mean(bafs, true), ArrayUtils.stdev(bafs, true));
  }

}
