package org.genvisis.cnv.stats.cdf;

import org.apache.commons.math3.distribution.NormalDistribution;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;

public class Cdf {

  public static double cdf(OpdfGaussian dist, ObservationReal n) {
    return new NormalDistribution(dist.mean(),
                                  Math.sqrt(dist.variance())).cumulativeProbability(n.value);
  }

}
