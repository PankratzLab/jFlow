/**
 * 
 */
package org.genvisis.cnv.hmm;

import org.pankratzlab.common.Logger;

/**
 * Essentially a port of
 * https://github.com/WGLab/PennCNV/blob/122691c8178b1ac803c60c802dacec83dc2593a4/kext/khmm.c#L650
 */
public class BaumWelchNP {

  private BaumWelchNP() {

  }

  private static final boolean testing = true;
  private static final double DELTA = 1;
  private static final double PARAM_CHANGE = 1e-9;

  /**
   * @param pennHmm an {@link PennHmm}
   * @param o1 lrrs for training
   * @param o2 bafs for training
   * @param pfb pfbs for training
   * @param snpdist snp distances for training
   * @param copyNumberOnlyDef copy number only probe definitions
   * @param log {@link Logger}
   */
  static void BaumWelchNP_CHMM(PennHmm pennHmm, double[] o1, double[] o2, double[] pfb,
                               final int[] snpdist, final boolean[] copyNumberOnlyDef, Logger log) {
    PennHmm.validateLengths(pennHmm, o1, o2, pfb, snpdist, copyNumberOnlyDef);
    double[][] biot = new double[pennHmm.getN()][o1.length];
    for (int t = 0; t < o1.length; t++) {
      double biosumtemp = 0;
      double[] biotemp = new double[pennHmm.getN()];
      for (int i = 0; i < pennHmm.getN(); i++) {
        if (copyNumberOnlyDef[t] || o2[t] > 1 || pfb[t] > 1) {
          biotemp[i] = Math.exp(PennHmm.b1iot(i, pennHmm.getB3(), o1[t]));
        } else {
          biotemp[i] = Math.exp(PennHmm.b1iot(i, pennHmm.getB1(), o1[t]));
        }
        biosumtemp += biotemp[i];
      }
      for (int i = 0; i < pennHmm.getN(); i++) {
        biot[i][t] = (biotemp[i] / biosumtemp);

      }
      biosumtemp = 0;
      for (int i = 0; i < pennHmm.getN(); i++) {
        if (copyNumberOnlyDef[t] || o2[t] > 1 || pfb[t] > 1) {
          biotemp[i] = 1;
        } else {
          biotemp[i] = Math.exp(PennHmm.b2iot(i, pennHmm.getB2(), pfb[t], o2[t]));
        }
        biosumtemp += biotemp[i];

      }
      for (int i = 0; i < pennHmm.getN(); i++) {
        biot[i][t] *= (biotemp[i] / biosumtemp);
      }
    }

    double[][] alpha = new double[o1.length][pennHmm.getN()];
    double[] scale = new double[o1.length];

    double logprobinit = ForwardWithScale_CHMM(pennHmm, o1, o2, pfb, snpdist, copyNumberOnlyDef,
                                               biot, alpha, scale, log);

    double logprobcurrent = logprobinit;
    double logprobLast = logprobinit;

    double[][] beta = new double[o1.length][pennHmm.getN()];
    BackwardWithScale_CHMM(pennHmm, o1, o2, pfb, snpdist, copyNumberOnlyDef, biot, beta, scale,
                           log);

    double[][] gamma = new double[o1.length][pennHmm.getN()];
    ComputeGamma_CHMM(pennHmm, alpha, beta, gamma, log);

    log.reportTimeInfo("Running ComputeXi_CHMM ");

    double[][][] xi = new double[o1.length][pennHmm.getN()][pennHmm.getN()];
    ComputeXi_CHMM(pennHmm, o1, o2, snpdist, copyNumberOnlyDef, biot, alpha, beta, xi, log);
    log.reportTimeInfo("Finished running ComputeXi_CHMM ");
    double delta = 0;
    int iter = 0;
    double adjust = 0.0;
    do {
      for (int i = 0; i < pennHmm.getN(); i++) {
        double denominatorA = 0.0;
        for (int t = 0; t < o1.length - 1; t++) {
          denominatorA += gamma[t][i];
        }
        for (int j = 0; j < pennHmm.getN(); j++) {
          double numeratorA = 0.0;
          for (int t = 0; t < o1.length - 1; t++) {
            numeratorA += xi[t][i][j];
          }
          pennHmm.getA()[i][j] = PARAM_CHANGE + (1 - PARAM_CHANGE) * numeratorA / denominatorA;
        }
      }
      double inflation = Math.floor(1e9 * (1 - Math.exp(-5000 / PennHmm.STATE_CHANGE))) / 1e9;
      double inflation_s4 = Math.floor(1e9 * (1 - Math.exp(-5000 / PennHmm.STATE_CHANGE / 1000)))
                            / 1e9;

      for (int i = 0; i < pennHmm.getN(); i++) {
        for (int j = 0; j < pennHmm.getN(); j++) {
          if (i != j && ((pennHmm.getA()[i][j] >= inflation && i != 3)
                         || (pennHmm.getA()[i][j] >= inflation_s4 && i == 3))) {
            if (i != 3) {
              adjust = pennHmm.getA()[i][j] - inflation;
            }
            if (i == 3) {
              adjust = pennHmm.getA()[i][j] - inflation_s4;
            }
            pennHmm.getA()[i][j] -= adjust;
            pennHmm.getA()[i][i] += adjust;

          }
        }
      }
      log.reportTimeInfo("Iteration: " + iter);

      logprobcurrent = ForwardWithScale_CHMM(pennHmm, o1, o2, pfb, snpdist, copyNumberOnlyDef, biot,
                                             alpha, scale, log);
      BackwardWithScale_CHMM(pennHmm, o1, o2, pfb, snpdist, copyNumberOnlyDef, biot, beta, scale,
                             log);

      ComputeGamma_CHMM(pennHmm, alpha, beta, gamma, log);
      ComputeXi_CHMM(pennHmm, o1, o2, snpdist, copyNumberOnlyDef, biot, alpha, beta, xi, log);
      delta = logprobcurrent - logprobLast;

      iter++;

      log.reportTimeInfo("Finished Baum-Welch iteration=" + iter + " delta=" + delta
                         + " current_logprobf=" + logprobcurrent + "last_logprobf=" + logprobLast);

      log.reportTimeInfo("\n" + PennHmm.PrintCHMM(pennHmm, log));

      logprobLast = logprobcurrent;

    } while (delta > DELTA);

  }

  private static void ComputeXi_CHMM(PennHmm pennHmm, double[] o1, double[] o2, final int[] snpdist,
                                     final boolean[] copyNumberOnlyDef, double[][] biot,
                                     double[][] alpha, double[][] beta, double[][][] xi,
                                     Logger log) {

    for (int t = 0; t < o1.length - 1; t++) {
      double sum = 0.0;
      PennHmm converted = PennHmm.convertHMMTransition(pennHmm, snpdist[t]);
      for (int i = 0; i < pennHmm.getN(); i++) {
        for (int j = 0; j < pennHmm.getN(); j++) {
          xi[t][i][j] = alpha[t][i] * beta[t + 1][j] * converted.getA()[i][j] * biot[j][t + 1];
          sum += xi[t][i][j];
        }
      }
      for (int i = 0; i < pennHmm.getN(); i++) {
        for (int j = 0; j < pennHmm.getN(); j++) {
          xi[t][i][j] /= sum;
        }
      }
    }
  }

  //  

  private static void ComputeGamma_CHMM(PennHmm pennHmm, double[][] alpha, double[][] beta,
                                        double[][] gamma, Logger log) {

    for (int t = 0; t < alpha.length; t++) {
      double denominator = 0.0;
      for (int j = 0; j < pennHmm.getN(); j++) {
        gamma[t][j] = alpha[t][j] * beta[t][j];
        denominator += gamma[t][j];
      }
      for (int i = 0; i < pennHmm.getN(); i++) {
        gamma[t][i] = gamma[t][i] / denominator;
      }
    }

  }

  private static void BackwardWithScale_CHMM(PennHmm pennHmm, double[] o1, double[] o2,
                                             double[] pfb, final int[] snpdist,
                                             final boolean[] copyNumberOnlyDef, double[][] biot,
                                             double[][] beta, double[] scale, Logger log) {
    for (int i = 0; i < pennHmm.getN(); i++) {
      beta[o1.length - 1][i] = 1.0 / scale[o1.length - 1];
    }

    for (int t = o1.length - 2; t >= 0; t--) {
      PennHmm converted = PennHmm.convertHMMTransition(pennHmm, snpdist[t]);
      for (int i = 0; i < pennHmm.getN(); i++) {
        double sum = 0.0;
        for (int j = 0; j < pennHmm.getN(); j++) {
          sum += converted.getA()[i][j] * biot[j][t + 1] * beta[t + 1][j];
        }
        beta[t][i] = sum / scale[t];
      }
    }
  }

  /**
   * Something like
   * https://github.com/WGLab/PennCNV/blob/08e6c03a1af9797544147727d1bf1a8d2f7e6452/kext/khmm.c#L575
   */
  private static double ForwardWithScale_CHMM(PennHmm pennHmm, double[] o1, double[] o2,
                                              double[] pfb, final int[] snpdist,
                                              final boolean[] copyNumberOnlyDef, double[][] biot,
                                              double[][] alpha, double[] scale, Logger log) {

    //    1. Initialization
    log.reportTimeInfo("1. Initialization");
    scale[0] = 0;
    for (int i = 0; i < pennHmm.getN(); i++) {
      alpha[0][i] = pennHmm.getPi()[i] * biot[i][0];

      scale[0] += alpha[0][i];
    }
    for (int i = 0; i < pennHmm.getN(); i++) {
      alpha[0][i] /= scale[0];

    }
    //    2. Induction 
    log.reportTimeInfo("2. Induction ");

    for (int t = 0; t < o1.length - 1; t++) {
      scale[t + 1] = 0.0;
      PennHmm converted = PennHmm.convertHMMTransition(pennHmm, snpdist[t]);
      for (int j = 0; j < pennHmm.getN(); j++) {
        double sum = 0;
        for (int i = 0; i < pennHmm.getN(); i++) {
          if (converted.getA()[i][j] < 0) {
            throw new IllegalStateException("found value less than 0 in A marix");
          }
          sum += alpha[t][i] * converted.getA()[i][j];
        }
        alpha[t + 1][j] = sum * biot[j][t + 1];
        scale[t + 1] += alpha[t + 1][j];
      }

      for (int j = 0; j < pennHmm.getN(); j++) {
        alpha[t + 1][j] /= scale[t + 1];
      }
    }
    log.reportTimeInfo("3. Termination ");
    double pprob = 0.0;
    for (int t = 0; t < o1.length; t++) {
      pprob += Math.log(scale[t]);
    }
    log.reportTimeInfo("Log prob init=" + pprob);
    if (testing) {
      double penn = -4228044.64024354;
      log.reportTimeInfo("PennCNV prob init=" + penn);
      log.reportTimeInfo("diff=" + Math.abs(pprob - penn));

    }
    return pprob;

  }

}
