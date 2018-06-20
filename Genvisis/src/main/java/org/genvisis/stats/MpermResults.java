package org.genvisis.stats;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.common.ArrayUtils;

public class MpermResults {

  public enum TestType {
    PLINK {

      public double permute(MpermResults res, List<Double> newPheno) {
        List<Integer> aIndices = res.getAIndices();
        List<Integer> bIndices = res.getBIndices();
        boolean twoSided = res.getTwoSided();
        double stat = res.getStat();
        double m0 = res.getM0();
        double m1 = res.getM1();
        double stat_j = -1;
        if (aIndices.size() > 0) {
          double aSum_j = 0.0;
          // calculate our "aSum" based on the random selection
          for (Integer k : aIndices) {
            aSum_j += newPheno.get(k);
          }
          double bSum_j = (m1 * aIndices.size() + m0 * bIndices.size()) - aSum_j;
          stat_j = aSum_j / aIndices.size() - bSum_j / bIndices.size();
          if (stat < 0 && !twoSided) stat_j = 0;
        }
        // increment r1, the number of greater differences in this region in the same direction as the original change
        res.updateR1(stat_j);
        return stat_j;
      }

      public double computeP(MpermResults res) {
        int r1 = res.getR1();
        int mperm = res.getMperm();
        return (double) (r1 + 1) / (double) (mperm + 1);
      }

      public double computeStat(MpermResults res) {
        double stat = res.getM1() - res.getM0();
        if (stat < 0 && !res.getTwoSided()) stat = 0;
        return stat;
      }
    },
    MW {

      public double worst() {
        return Double.MAX_VALUE;
      }

      public double betterStat(double a, double b) {
        return Math.min(a, b);
      }

      public double permute(MpermResults res, List<Double> newPheno) {
        List<Integer> aIndices = res.getAIndices();
        List<Integer> bIndices = res.getBIndices();
        double stat_j = 1;
        if (aIndices.size() > 0) {
          double u = MannWhitneyUTest.mannWhitneyU(newPheno, aIndices, bIndices);
          stat_j = MannWhitneyUTest.calculatePValue(u, aIndices.size(), bIndices.size());
        } else {
          stat_j = 1;
        }
        return stat_j;
      }

      public double computeP(MpermResults res) {
        return res.getStat();
      }

      public double computeStat(MpermResults res) {
        double stat = 1;
        double[] aArray = res.aArray;
        double[] bArray = res.bArray;
        if (aArray.length > 0) {
          double u = MannWhitneyUTest.mannWhitneyU(aArray, bArray);
          stat = MannWhitneyUTest.calculatePValue(u, aArray.length, bArray.length);
        }
        return stat;
      }
    },
    LOGISTIC {

      public double permute(MpermResults res, List<Double> newPheno) {
        List<Integer> aIndices = res.getAIndices();
        int[] cnvs = res.getCnvs();
        double stat_j = 0;
        if (aIndices.size() > 0) {
          double[] phenoArray = new double[newPheno.size()];
          for (int i = 0; i < newPheno.size(); i++) {
            phenoArray[i] = newPheno.get(i);
          }
          LogisticRegression lr = new LogisticRegression(phenoArray,
                                                         ArrayUtils.toDoubleArray(cnvs));
          stat_j = lr.getBetas()[0];
        } else {
          stat_j = 0;
        }
        return stat_j;
      }

      public double computeP(MpermResults res) {
        double stat = res.getStat();
        double se = res.getSE();
        return ProbDist.NormDist(Math.abs(stat / se));
      }

      public double computeStat(MpermResults res) {
        double stat = 0;
        if (res.aIndices.size() > 0) {
          List<Double> pheno = res.pheno;
          double[] phenoArray = new double[pheno.size()];
          for (int i = 0; i < pheno.size(); i++) {
            phenoArray[i] = pheno.get(i);
          }
          LogisticRegression lr = new LogisticRegression(phenoArray,
                                                         ArrayUtils.toDoubleArray(res.cnvs));
          stat = lr.getBetas()[0];
          res.se = lr.getSEofBs()[0];
          res.p1 = lr.getSigs()[0];
        } else {
          stat = 0;
          res.se = 0;
          res.p1 = Double.NaN;
        }
        return stat;
      }
    };

    public double worst() {
      return Double.MIN_VALUE;
    }

    public double betterStat(double a, double b) {
      return Math.max(a, b);
    }

    public double permute(MpermResults res, List<Double> newPheno) {
      return worst();
    }

    public double computeP(MpermResults res) {
      return Double.NaN;
    }

    public double computeStat(MpermResults res) {
      return worst();
    }
  }

  private byte chr;
  private int pos, r1, r2, mperm;
  private double m1, m0, stat, p1, p2, se;
  private String name;
  private double[] aArray, bArray;
  private boolean twoSided;
  private List<Integer> aIndices, bIndices;
  private List<Double> pheno;
  private TestType type;
  private int[] cnvs;

  public MpermResults(byte chr, int pos, String name, List<Double> pheno, int[] cnvs, int mperm,
                      boolean twoSided, TestType type) {

    this.chr = chr;
    this.pos = pos;
    this.name = name;
    this.mperm = mperm;
    this.r1 = 0;
    this.r2 = 0;
    this.twoSided = twoSided;
    this.type = type;
    this.cnvs = cnvs;
    this.pheno = pheno;

    this.aIndices = new ArrayList<Integer>();
    this.bIndices = new ArrayList<Integer>();

    for (int i = 0; i < cnvs.length; i++) {
      if (cnvs[i] == 1) {
        aIndices.add(i);
      } else {
        bIndices.add(i);
      }
    }

    this.aArray = new double[aIndices.size()];
    this.bArray = new double[bIndices.size()];

    for (int i = 0; i < aIndices.size(); i++) {
      aArray[i] = pheno.get(aIndices.get(i));
    }

    for (int i = 0; i < bIndices.size(); i++) {
      bArray[i] = pheno.get(bIndices.get(i));

    }

    computeMeans();
  }

  protected double getSE() {
    return se;
  }

  protected int getMperm() {
    return mperm;
  }

  protected int getR1() {
    return r1;
  }

  protected int[] getCnvs() {
    return cnvs;
  }

  protected double getM0() {
    return m0;
  }

  protected double getM1() {
    return m1;
  }

  protected boolean getTwoSided() {
    return twoSided;
  }

  protected List<Integer> getBIndices() {
    return bIndices;
  }

  protected List<Integer> getAIndices() {
    return aIndices;
  }

  public void computeMeans() {
    m1 = ArrayUtils.mean(aArray);
    m0 = ArrayUtils.mean(bArray);

    stat = type.computeStat(this);
  }

  public double getP() {
    if (p1 == 0) computeP();
    return p1;
  }

  public double getP2() {
    if (p2 == 0) computeP2();
    return p2;
  }

  public double getStat() {
    return stat;
  }

  public void computeP() {
    if (aArray.length == 0) {
      p1 = Double.NaN;
    } else {
      p1 = type.computeP(this);
    }
  }

  public void computeP2() {
    if (aArray.length == 0) {
      p2 = Double.NaN;
    } else {
      p2 = (double) (r2 + 1) / (double) (mperm + 1);
    }
  }

  public double permute(List<Double> newPheno) {
    return type.permute(this, newPheno);
  }

  public void updateR1(double bestStat) {
    if (Double.isNaN(stat) || type.betterStat(stat, bestStat) == bestStat) r1++;
  }

  public void updateR2(double bestStat) {
    if (Double.isNaN(stat) || type.betterStat(stat, bestStat) == bestStat) r2++;
  }

  public String toString() {
    return chr + "\t" + pos + "\t" + name + "\t" + aIndices.size() + "\t" + m1 + "\t" + m0 + "\t"
           + getP() + "\t" + getP2();
  }
}
