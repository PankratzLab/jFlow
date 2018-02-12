/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import org.genvisis.cnv.wbs.ChangePoint.ChangePointComparable;
import org.genvisis.cnv.wbs.ChangePoint.MinthComparable;
import org.genvisis.cnv.wbs.Penalties.Penalty_type;
import org.genvisis.common.ArrayUtils;
import com.google.common.primitives.Doubles;

/**
 * Port of https://github.com/cran/wbs/blob/master/R/changepoints.R
 */
class ChangePoints {

  // https://www.math.ucla.edu/~anderson/rw1001/library/base/html/mad.html
  // https://en.wikipedia.org/wiki/Median_absolute_deviation
  private static final double MAD_FACTOR = 1.4826;
  static final double DEFAULT_TH_CONST = 1.3;

  private double sigma;
  private double th;
  private double thConst;
  private int kMax;
  private List<ChangePoint> standardTHChangePoints;
  private List<PenalizedChangePoints> penalizedChangePoints;

  private ChangePoints(double sigma, double th, double thConst, int kMax,
                       List<ChangePoint> standardTHChangePoints) {
    super();
    this.sigma = sigma;
    this.th = th;
    this.thConst = thConst;
    this.kMax = kMax;
    this.standardTHChangePoints = standardTHChangePoints;
    this.penalizedChangePoints = new ArrayList<>();
  }

  static class PenalizedChangePoints {

    private List<ChangePoint> penalizedChangePoints;
    private Penalty_type penaltyType;
    private List<Double> icCurve;

    private PenalizedChangePoints(Penalty_type penaltyType) {
      super();
      this.penalizedChangePoints = new ArrayList<>();
      this.penaltyType = penaltyType;
      this.icCurve = new ArrayList<>();
    }

    List<ChangePoint> getPenalizedChangePoints() {
      return penalizedChangePoints;
    }

    Penalty_type getPenaltyType() {
      return penaltyType;
    }

    List<Double> getIcCurve() {
      return icCurve;
    }

    // pretty sure we just find min index

    // https://github.com/cran/wbs/blob/8919e7c35389b92e27b6948572271e0843b5f808/R/changepoints.R#L209
    private void parseChangePoints(List<ChangePoint> wlist) {
      double minVal = Double.MAX_VALUE;
      int minIndex = 0;
      for (int i = 0; i < icCurve.size(); i++) {
        if (icCurve.get(i).doubleValue() <= minVal) {
          minVal = icCurve.get(i).doubleValue();
          minIndex = i;
        }
      }
      if (minIndex > 0) {
        penalizedChangePoints = wlist.subList(0, minIndex);
      }
    }
  }

  List<PenalizedChangePoints> getPenalizedChangePoints() {
    return penalizedChangePoints;
  }

  List<ChangePoint> getStandardTHChangePoints() {
    return standardTHChangePoints;
  }

  double getSigma() {
    return sigma;
  }

  double getTh() {
    return th;
  }

  double getThConst() {
    return thConst;
  }

  int getkMax() {
    return kMax;
  }

  static double computeSigma(double[] x) {
    return ArrayUtils.mad(WBSUtilities.rootTwo(WBSUtilities.computeLagOneDifference(x)),
                          MAD_FACTOR);

  }

  static double computTh(double thConst, double sigma, int n, int kmax, List<ChangePoint> wlist) {
    if (kmax > 0) {
      if (kmax < wlist.size()) {
        Collections.sort(wlist, new MinthComparable());
        ChangePoint kmaxCP = wlist.get(kmax - 1);
        double min = getMinLagOneCPDiff(wlist);
        Collections.sort(wlist); // Sort by cusum

        return kmaxCP.getMinth() - min / 2;

      } else {
        return 0;
      }
    } else {
      return thConst * sigma * Math.sqrt(2 * Math.log(n));

    }
  }

  private static double getMinLagOneCPDiff(List<ChangePoint> wlist) {
    List<Double> minths = new ArrayList<>();
    HashSet<Double> have = new HashSet<>();
    for (int i = 0; i < wlist.size(); i++) {
      if (!have.contains(wlist.get(i).getMinth())) {
        Double m = new Double(wlist.get(i).getMinth());
        have.add(m);
        minths.add(m);
      }
    }
    double[] aMinths = WBSUtilities.computeLagOneDifference(ArrayUtils.reverse(Doubles.toArray(minths)));
    return ArrayUtils.min(aMinths);
  }

  /**
   * Standard Binary Segmentation changepoint algorithm
   */
  static ChangePoints changepointsSbs(List<ChangePoint> wlist, double[] x, double thConst,
                                      int kmax) {
    double sigma = computeSigma(x);

    double th = computTh(thConst, sigma, x.length, kmax, wlist);

    List<ChangePoint> filteredChangePoints = new ArrayList<>();
    for (ChangePoint changePoint : wlist) {
      if (changePoint.getMinth() > th) {
        filteredChangePoints.add(changePoint);
      }
    }

    Collections.sort(filteredChangePoints); // Sort by cusum

    return new ChangePoints(sigma, th, thConst, filteredChangePoints.size(), filteredChangePoints);
  }

  static ChangePoints changepointsWbs(List<ChangePoint> wlist, double[] x, double thConst,
                                      List<Penalty_type> penaltyTypes, int kmax) {

    ChangePoints changePoints = changepointsSbs(wlist, x, thConst, -1);

    if (!penaltyTypes.isEmpty() && !changePoints.standardTHChangePoints.isEmpty()) {

      ChangePoints cptCand = changepointsSbs(wlist, x, thConst, kmax);

      List<PenalizedChangePoints> penalizedChangePoints = new ArrayList<>();
      for (Penalty_type type : penaltyTypes) {
        PenalizedChangePoints p = new PenalizedChangePoints(type);
        p.icCurve.add((double) (x.length / 2) * Math.log(ArrayUtils.variance(x)));
        penalizedChangePoints.add(p);
      }

      for (int i = 1; i <= cptCand.getStandardTHChangePoints().size(); i++) {
        List<ChangePoint> sub = cptCand.getStandardTHChangePoints().subList(0, i);
        double minLogLik = computeMinLogLik(x, sub);
        for (int j = 0; j < penaltyTypes.size(); j++) {
          penalizedChangePoints.get(j).icCurve.add(new Double(minLogLik
                                                              + Penalties.getPenalty(x, sub,
                                                                                     penaltyTypes.get(j))));
        }
      }
      for (PenalizedChangePoints pcp : penalizedChangePoints) {
        pcp.parseChangePoints(wlist);
      }
      changePoints.penalizedChangePoints = penalizedChangePoints;
      return changePoints;

    } else {
      return changePoints;
    }
  }

  static double computeMinLogLik(double[] x, List<ChangePoint> wlist) {

    double sumSquared = 0;
    double[] means = computeMeansBetweenChangepoints(x, wlist);
    for (int i = 0; i < means.length; i++) {
      sumSquared += Math.pow(x[i] - means[i], 2);
    }
    return (double) x.length / 2 * Math.log(sumSquared / x.length);
  }

  static double[] computeMeansBetweenChangepoints(double[] x, List<ChangePoint> wlist) {
    double[] means = new double[x.length];
    // sort by chpt index
    Collections.sort(wlist, new ChangePointComparable());

    assignMean(means, x, wlist.get(0)

                              .getChangePoint(),
               0);

    assignMean(means, x, x.length, wlist.get(wlist.size() - 1).getChangePoint());

    for (int i = 0; i < wlist.size() - 1; i++) {
      assignMean(means, x, wlist.get(i + 1).getChangePoint(), wlist.get(i).getChangePoint());
    }
    // put back in order
    Collections.sort(wlist);
    return means;
  }

  private static void assignMean(double[] means, double[] x, int stop, int start) {
    double mean = ArrayUtils.mean(ArrayUtils.subArray(x, start, stop));
    for (int j = start; j < stop; j++) {
      means[j] = mean;
    }
  }

}
