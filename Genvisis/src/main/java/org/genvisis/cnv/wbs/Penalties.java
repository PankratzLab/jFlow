/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.genvisis.cnv.wbs.ChangePoint.ChangePointComparable;

import com.google.common.primitives.Doubles;

/**
 * Information criterion penalties
 */
class Penalties {

  static final double DEFAULT_ALPHA = 1.01;

  private Penalties() {

  }

  enum Penalty_type {
    SSIC, BIC, MBIC;
  }

  static double getPenalty(double[] x, List<ChangePoint> wlist, Penalty_type pType) {
    switch (pType) {
      case BIC:
        return bicPenalty(x, wlist);
      case MBIC:
        return mbicPenalty(x, wlist);
      case SSIC:
        return ssicPenalty(x, DEFAULT_ALPHA, wlist);
      default:
        throw new IllegalArgumentException("Invalid penalty type " + pType);

    }

  }

  /**
   * https://github.com/cran/wbs/blob/8919e7c35389b92e27b6948572271e0843b5f808/R/penalties.R#L17-L28
   * We only support "log"
   */
  private static double ssicPenalty(double[] x, double alpha, List<ChangePoint> wlist) {
    double pen = Math.pow(Math.log(x.length), alpha);
    return (double) wlist.size() * pen;
  }

  private static double bicPenalty(double[] x, List<ChangePoint> wlist) {
    return wlist.size() * Math.log(x.length);
  }

  private static double mbicPenalty(double[] x, List<ChangePoint> wlist) {
    int k = wlist.size();
    if (k > 0) {
      return ((double) 3 / 2) * k * Math.log(x.length) + ((double) 1 / 2) * sumLogDiff(x, wlist);
    } else {
      return ((double) 1 / 2) * sumLogDiff(x, wlist);
    }
  }

  private static double sumLogDiff(double[] x, List<ChangePoint> wlist) {
    List<Double> cps = new ArrayList<>();
    cps.add(new Double(0.0));
    Collections.sort(wlist, new ChangePointComparable());
    for (ChangePoint changePoint : wlist) {
      cps.add((double) changePoint.getChangePoint());
    }
    cps.add(new Double(x.length));

    double[] diff = WBSUtilities.computeLagOneDifference(Doubles.toArray(cps));
    double sum = 0;
    for (int i = 0; i < diff.length; i++) {
      sum += Math.log(diff[i] / x.length);
    }
    Collections.sort(wlist);
    return sum;

  }
}
