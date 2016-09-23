package org.genvisis.cnv.hmm;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariant.CNVBuilder;
import org.genvisis.filesys.LocusSet;
import org.genvisis.stats.Stats;

import com.google.common.primitives.Ints;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;

/*
 * This file was adapted from several subroutines from the UMDHMM package by Tapas Kanungo (Date: 15
 * December 1997) The original UMDHMM package was downloaded from
 * http://www.kanungo.com/software/software.html. The citation for the UMDHMM program is
 * "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai
 * (editor), Cambridge University Press, 1999." The basic framework (including variable name,
 * subroutine name) is highly similar to the original UMDHMM package, but the actual implementation
 * is completely different as no "discrete symbol emission" is used in PennCNV.
 */
/**
 * @author lane0212 <br>
 *         Mimics the hmm functionality used in in PennCNV (http://penncnv.openbioinformatics.org/)
 *         <br>
 *         Really, we mimic the kext C package here <br>
 *         See below for an example .hmm file <br>
 *         NOTE<br>
 *         the model has 6 states;<br>
 *         Index 0: CN 0<br>
 *         Index 1: CN 1<br>
 *         Index 2: CN 2<br>
 *         Index 3: LOH!<br>
 *         Index 4: CN 3<br>
 *         Index 5: CN 4<br>
 *
 *         NOTE: <br>
 *         kext starts at index 1
 *
 *
 *
 */

// M=6
// N=6
// A:
// 0.936719716 0.006332139 0.048770575 0.000000001 0.008177573 0.000000001
// 0.000801036 0.949230924 0.048770575 0.000000001 0.001168245 0.000029225
// 0.000004595 0.000047431 0.999912387 0.000000001 0.000034971 0.000000621
// 0.000049998 0.000049998 0.000049998 0.999750015 0.000049998 0.000049998
// 0.000916738 0.001359036 0.048770575 0.000000001 0.948953653 0.000000002
// 0.000000001 0.000000001 0.027257213 0.000000001 0.000000004 0.972742785
// B:
// 0.950000 0.000001 0.050000 0.000001 0.000001 0.000001
// 0.000001 0.950000 0.050000 0.000001 0.000001 0.000001
// 0.000001 0.000001 0.999995 0.000001 0.000001 0.000001
// 0.000001 0.000001 0.050000 0.950000 0.000001 0.000001
// 0.000001 0.000001 0.050000 0.000001 0.950000 0.000001
// 0.000001 0.000001 0.050000 0.000001 0.000001 0.950000
// pi:
// 0.000001 0.000500 0.999000 0.000001 0.000500 0.000001
// B1_mean:
// -3.527211 -0.664184 0.000000 100.000000 0.395621 0.678345
// B1_sd:
// 1.329152 0.284338 0.159645 0.211396 0.209089 0.191579
// B1_uf:
// 0.010000
// B2_mean:
// 0.000000 0.250000 0.333333 0.500000 0.500000
// B2_sd:
// 0.016372 0.042099 0.045126 0.034982 0.304243
// B2_uf:
// 0.010000
// B3_mean:
// -2.051407 -0.572210 0.000000 0.000000 0.361669 0.626711
// B3_sd:
// 2.132843 0.382025 0.184001 0.200297 0.253551 0.353183
// B3_uf:
// 0.010000

public class PennHmm {
  public static final int LOH_FLAG = 99;
  private static final double NOT_ZERO_PI = 0.000000001; // 1e-9
  private static final double STATE_CHANGE = 100000.0;
  private static final double VITHUGE = 100000000000.0;
  private static final double FLOAT_MINIMUM = 1.175494351e-38;
  private final int M;
  private final int N;// number of states
  private double[] pi;
  private double[][] a;
  private final BStatus B1;// for LRR measure from SNP markers
  private final BStatus B2;// for BAF measure from SNP markers
  private final BStatus B3;// for LRR measure from CN only markers

  private final Logger log;

  public PennHmm(int m, int n, double[] pi, double[][] a, BStatus b1, BStatus b2, BStatus b3,
                 Logger log) {
    super();
    M = m;
    N = n;
    this.pi = pi;
    this.a = a;
    this.B1 = b1;
    this.B2 = b2;
    this.B3 = b3;
    this.log = log;
  }

  /**
   * @param pennHmm make a clone essentially
   */
  public PennHmm(PennHmm pennHmm) {
    M = pennHmm.M;
    N = pennHmm.N;
    pi = pennHmm.pi;
    a = pennHmm.a;
    B1 = new BStatus(pennHmm.B1.b_mean, pennHmm.B1.b_sd, pennHmm.B1.b_uf);
    B2 = new BStatus(pennHmm.B2.b_mean, pennHmm.B2.b_sd, pennHmm.B2.b_uf);
    B3 = new BStatus(pennHmm.B3.b_mean, pennHmm.B3.b_sd, pennHmm.B3.b_uf);
    log = pennHmm.log;
  }

  @Override
  public String toString() {
    return "PennHmm [M=" + M + ", N=" + N + ", pi=" + Arrays.toString(pi) + ", a="
           + Arrays.toString(a) + ", B1=" + B1.toString() + ", B2=" + B2.toString() + ", B3="
           + B3.toString() + "]";
  }

  public int getN() {
    return N;
  }

  public double[][] getA() {
    return a;
  }

  public void setA(double[][] a) {
    this.a = a;
  }

  public Logger getLog() {
    return log;
  }

  public BStatus getB1() {
    return B1;
  }

  public double[] getPi() {
    return pi;
  }

  public BStatus getB2() {
    return B2;
  }

  public BStatus getB3() {
    return B3;
  }

  /**
   * take log of pi values<br>
   * WARNING : modifies internal pi array;
   */
  private void logPi() {
    double[] loggedPi = new double[pi.length];
    for (int i = 0; i < loggedPi.length; i++) {
      double pt = pi[i];
      if (pt == 0) {/* eliminate problems with zero probability */
        pt = NOT_ZERO_PI;
      }
      loggedPi[i] = Math.log(pt);
    }
    pi = loggedPi;
  }

  /**
   * @param state hmm state
   * @param bStatus B1, B2...
   * @param o the intensity observation (LRR value)
   * @return
   */
  private static double b1iot(int state, BStatus bStatus, double o) {
    double p = 0;
    p = bStatus.getB_uf();
    p +=
      (1 - bStatus.getB_uf()) * bStatus.getGaussians()[state].probability(new ObservationReal(o));

    if (p == 0) {
      p = FLOAT_MINIMUM;

    }
    return Math.log(p);
  }

  /**
   * @param state hmm state
   * @param bStatus
   * @param pfb
   * @param b baf value to test
   * @return log odds, I think, for the current data and state
   */
  private static double b2iot(int state, BStatus bStatus, double pfb, double b) {

    double uf = bStatus.getB_uf();
    double p = uf;
    // from kext....
    // double mean0 = bStatus.getB_mean()[1];//note the index starts at 1
    // double mean25 = bStatus.getB_mean()[2];
    // double mean33 = bStatus.getB_mean()[3];
    // double mean50 = bStatus.getB_mean()[4];
    // double mean50_state1 = bStatus.getB_mean()[5];
    // double sd0 = bStatus.getB_sd()[1];
    // double sd25 = bStatus.getB_sd()[2];
    // double sd33 = bStatus.getB_sd()[3];
    // double sd50 = bStatus.getB_sd()[4];
    // double sd50_state1 = bStatus.getB_sd()[5];

    // if (state == 1) {
    // if (b==0) {
    // p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
    // } else if (b==1) {
    // p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
    // } else {
    // p+= (1-uf) * pdf_normal (b, mean50_state1, sd50_state1);
    // }
    // }
    ObservationReal o = new ObservationReal(b);
    if (state == 0) {
      OpdfGaussian opdfGaussian = bStatus.getGaussians()[4];

      if (b == 0) {
        p += (1 - uf) * Stats.cdf(opdfGaussian, new ObservationReal(0));
      } else if (b == 1) {
        p += (1 - uf) * Stats.cdf(opdfGaussian, new ObservationReal(0));
      } else {
        p += (1 - uf) * opdfGaussian.probability(o);
      }

    } else if (state == 1) {
      if (b == 0) {
        p += (1 - uf) * (1 - pfb) / 2;
      } else if (b == 1) {
        p += (1 - uf) * pfb / 2;
      } else {
        OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
        OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0],
                                                          Math.pow(bStatus.getB_sd()[0], 2));

        p += (1 - uf) * (1 - pfb) * opdfGaussian.probability(o);
        p += (1 - uf) * pfb * opdfGaussianMinus.probability(o);
      }
    } else if (state == 2) {
      if (b == 0) {
        p += (1 - uf) * (1 - pfb) * (1 - pfb) / 2;
      } else if (b == 1) {
        p += (1 - uf) * pfb * pfb / 2;
      } else {
        OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
        OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0],
                                                          Math.pow(bStatus.getB_sd()[0], 2));
        OpdfGaussian opdfGaussian5 = bStatus.getGaussians()[3];
        p += (1 - uf) * (1 - pfb) * (1 - pfb) * opdfGaussian.probability(o);
        p += (1 - uf) * 2 * pfb * (1 - pfb) * opdfGaussian5.probability(o);
        p += (1 - uf) * pfb * pfb * opdfGaussianMinus.probability(o);
      }
    } else if (state == 3) {
      if (b == 0) {
        p += (1 - uf) * (1 - pfb) / 2;
      } else if (b == 1) {
        p += (1 - uf) * pfb / 2;
      } else {
        OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
        OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0],
                                                          Math.pow(bStatus.getB_sd()[0], 2));

        p += (1 - uf) * (1 - pfb) * opdfGaussian.probability(o);
        p += (1 - uf) * pfb * opdfGaussianMinus.probability(o);
      }
    } else if (state == 4) {
      if (b == 0) {
        p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) / 2;
      } else if (b == 1) {
        p += (1 - uf) * pfb * pfb * pfb / 2;
      } else {
        OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
        OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0],
                                                          Math.pow(bStatus.getB_sd()[0], 2));
        OpdfGaussian opdfGaussian33 = bStatus.getGaussians()[2];
        OpdfGaussian opdfGaussianMinus33 = new OpdfGaussian(1 - bStatus.getB_mean()[2],
                                                            Math.pow(bStatus.getB_sd()[2], 2));

        p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * opdfGaussian.probability(o);
        p += (1 - uf) * 3 * (1 - pfb) * (1 - pfb) * pfb * opdfGaussian33.probability(o);
        p += (1 - uf) * 3 * (1 - pfb) * pfb * pfb * opdfGaussianMinus33.probability(o);
        p += (1 - uf) * pfb * pfb * pfb * opdfGaussianMinus.probability(o);
      }
    } else if (state == 5) {
      if (b == 0) {
        p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) / 2;
      } else if (b == 1) {
        p += (1 - uf) * pfb * pfb * pfb * pfb / 2;
      } else {
        OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
        OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0],
                                                          Math.pow(bStatus.getB_sd()[0], 2));
        OpdfGaussian opdfGaussian25 = bStatus.getGaussians()[1];
        OpdfGaussian opdfGaussianMinus25 = new OpdfGaussian(1 - bStatus.getB_mean()[1],
                                                            Math.pow(bStatus.getB_sd()[1], 2));
        OpdfGaussian opdfGaussian5 = bStatus.getGaussians()[3];

        p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) * opdfGaussian.probability(o);
        p += (1 - uf) * 4 * (1 - pfb) * (1 - pfb) * (1 - pfb) * pfb * opdfGaussian25.probability(o);
        p += (1 - uf) * 6 * (1 - pfb) * (1 - pfb) * pfb * pfb * opdfGaussian5.probability(o);
        p += (1 - uf) * 4 * (1 - pfb) * pfb * pfb * pfb * opdfGaussianMinus25.probability(o);
        p += (1 - uf) * pfb * pfb * pfb * pfb * opdfGaussianMinus.probability(o);
      }
    }
    if (!Numbers.isFinite(p)) {
      String error = "Non-finite p, " + p;
      throw new IllegalStateException(error);

    }
    if (p == 0) {
      p = FLOAT_MINIMUM;
    }
    return Math.log(p);
  }

  /**
   * Load the hmm file, must be in exact PennCNV/Kext format
   */
  public static PennHmm loadPennHmm(String hmmFile, Logger log) {
    if (!Files.exists(hmmFile)) {
      log.reportFileNotFound(hmmFile);
      return null;
    } else {
      try {
        BufferedReader reader = Files.getAppropriateReader(hmmFile);
        String M = reader.readLine().trim();
        int m = -1;
        int n = -1;
        double[][] a = null;
        if (!M.startsWith("M=")) {
          log.reportTimeError("cannot read M annotation from HMM file");
          return null;
        } else {
          m = Integer.parseInt(M.replaceAll("M=", ""));
        }
        String N = reader.readLine().trim();
        if (!N.startsWith("N=")) {
          log.reportTimeError("cannot read N annotation from HMM file");
          return null;
        } else {
          n = Integer.parseInt(N.replaceAll("N=", ""));
        }
        String aTag = reader.readLine().trim();
        if (!aTag.startsWith("A:")) {
          log.reportTimeError("cannot read A: annotation from HMM file");
          return null;
        } else {
          a = loadMatrix(reader, m, n);
        }
        String bTag = reader.readLine().trim();
        if (!bTag.startsWith("B:")) {
          log.reportTimeError("cannot read B: annotation from HMM file");
          return null;
        } else {
          loadMatrix(reader, m, n);
        }
        double[] pi = loadTwoLineDoubleArray("pi", reader, log);
        BStatus B1 = loadBstatus("B1", reader, log);
        BStatus B2 = loadBstatus("B2", reader, log);
        BStatus B3 = loadBstatus("B3", reader, log);
        reader.close();
        return new PennHmm(m, n, pi, a, B1, B2, B3, log);
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + hmmFile + "\" not found in current directory");
        return null;
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + hmmFile + "\"");
        return null;
      }
    }
  }

  public static ViterbiResult ViterbiLogNP_CHMM(PennHmm pennHmm, double[] o1, double[] o2,
                                                double[] pfb, final int[] snpdist,
                                                final boolean[] copyNumberOnlyDef) {
    if (o1.length != o2.length || o1.length != pfb.length || o1.length != snpdist.length
        || o1.length != copyNumberOnlyDef.length) {
      String error = "BUG: mismatched array lengths";
      pennHmm.getLog().reportTimeError(error);
      error += "\nO1 Length: " + o1.length + "\n";
      error += "O2 Length: " + o2.length + "\n";
      error += "pfb Length: " + pfb.length + "\n";

      error += "SnpDist Length: " + snpdist.length + "\n";
      error += "CN Length: " + copyNumberOnlyDef.length + "\n";

      throw new IllegalArgumentException(error);
    }

    double[][] biot = new double[pennHmm.getN()][o1.length];
    int[] q = new int[o1.length];
    double[][] delta = new double[o1.length][pennHmm.getN()];
    int[][] psi = new int[o1.length][pennHmm.getN()];

    PennHmm pennHmmLog = new PennHmm(pennHmm);
    pennHmmLog.logPi();
    for (int i = 0; i < pennHmmLog.getN(); i++) {
      for (int t = 0; t < o1.length; t++) {
        if (copyNumberOnlyDef[t] || o2[t] > 1 || pfb[t] > 1) {
          biot[i][t] = b1iot(i, pennHmmLog.getB3(), o1[t]);
        } else {
          double bioTmp = b1iot(i, pennHmmLog.getB1(), o1[t]);
          bioTmp += b2iot(i, pennHmmLog.getB2(), pfb[t], o2[t]);
          biot[i][t] = bioTmp;
        }
      }
    }

    for (int i = 0; i < pennHmmLog.getN(); i++) {
      delta[0][i] = pennHmmLog.getPi()[i] + biot[i][0];
      psi[0][i] = 0;
    }
    for (int t = 1; t < o1.length; t++) {
      PennHmm converted = convertHMMTransition(pennHmmLog, snpdist[t - 1]);// account for physical
                                                                           // distance between
                                                                           // markers
      for (int j = 0; j < converted.getN(); j++) {
        double maxval = -1 * VITHUGE;
        int maxvalind = 1;
        for (int i = 0; i < converted.getN(); i++) {
          double val = delta[t - 1][i] + Math.log(converted.getA()[i][j]);
          if (val > maxval) {
            maxval = val;
            maxvalind = i;
          }
        }
        delta[t][j] = maxval + biot[j][t];
        psi[t][j] = maxvalind;
      }
    }

    double pprob = -1 * VITHUGE;
    q[o1.length - 1] = 1;
    for (int i = 0; i < pennHmmLog.getN(); i++) {
      if (delta[o1.length - 1][i] > pprob) {
        pprob = delta[o1.length - 1][i];
        q[o1.length - 1] = i;
      }
    }
    for (int t = o1.length - 2; t >= 0; t--) {
      q[t] = psi[t + 1][q[t + 1]];
    }
    return new ViterbiResult(q, delta);
  }

  /**
   * Stores the q state array and delta matrix
   *
   */
  public static class ViterbiResult {
    private final int[] q;
    private final double[][] delta;
    private ArrayList<int[]> indexStateChange;

    /**
     * @param q states
     * @param delta transitions
     */
    public ViterbiResult(int[] q, double[][] delta) {
      super();
      this.q = q;
      this.delta = delta;
    }

    public double[][] getDelta() {
      return delta;
    }

    public ArrayList<int[]> getIndexStateChange() {
      return indexStateChange;
    }

    public int[] getQ() {
      return q;
    }

    /**
     * @param proj
     * @param fid assign all cnvs consolidated from the state sequence to this family id
     * @param iid assign all cnvs consolidated from the state sequence to this individual id
     * @param currentChr assign all cnvs consolidated from the state sequence to this chromosome
     * @param positions physical postions corresponding to the states
     * @param normalState Copy number normal
     * @return
     */
    public LocusSet<CNVariant> analyzeStateSequence(Project proj, String fid, String iid,
                                                    byte currentChr, int[] positions,
                                                    String[] names, int normalState,
                                                    boolean reverse, boolean verbose) {
      CNVariant.CNVBuilder builder = new CNVBuilder();
      builder.familyID(fid);
      builder.individualID(iid);
      indexStateChange = new ArrayList<int[]>();
      ArrayList<CNVariant> tmp = new ArrayList<CNVariant>();
      if (positions.length != q.length) {
        String error = "Have " + q.length + " state sequences, but " + positions.length
                       + " positions";
        proj.getLog().reportTimeError(error);
        throw new IllegalArgumentException(error);
      } else {
        boolean foundSignal = false;
        int currentFind = 2;
        int[] states = q;
        if (reverse) {
          states = Array.reverse(states);
        }
        for (int i = 0; i < states.length; i++) {

          int currentCN = getStateCN(states[i]);

          if (currentCN != normalState) {// CN 3 denotes LOH.

            if (foundSignal && currentCN != currentFind) {// new, adjacent cnv
              builder.stop(positions[i - 1]);
              indexStateChange.add(new int[] {builder.getStartIndex(), i - 1});
              tmp.add(builder.build());
              builder = new CNVBuilder();
              builder.familyID(fid);
              builder.individualID(iid);
              builder.chr(currentChr);
              builder.start(positions[i]);
              builder.startIndex(i);
              builder.numMarkers(1);
              builder.cn(currentCN);
              currentFind = currentCN;
            } else if (foundSignal) {// continue with previous cnv
              builder.numMarkers(builder.getNumMarkers() + 1);
            } else {// new cnv
              builder = new CNVBuilder();
              builder.familyID(fid);
              builder.individualID(iid);
              builder.chr(currentChr);
              builder.start(positions[i]);
              builder.startIndex(i);

              builder.numMarkers(1);
              builder.cn(currentCN);
              foundSignal = true;
              currentFind = currentCN;
            }
          } else if (foundSignal) {
            currentFind = normalState;
            builder.stop(positions[i - 1]);
            indexStateChange.add(new int[] {builder.getStartIndex(), i - 1});
            tmp.add(builder.build());
            builder = new CNVBuilder();
            builder.familyID(fid);
            builder.individualID(iid);
            foundSignal = false;
          } else {
            currentFind = normalState;
            builder = new CNVBuilder();
            builder.familyID(fid);
            builder.individualID(iid);
            foundSignal = false;
          }
        }
        if (foundSignal) {
          builder.stop(positions[positions.length - 1]);
          indexStateChange.add(new int[] {builder.getStartIndex(), q.length - 1});
          tmp.add(builder.build());
        }
      }

      LocusSet<CNVariant> cnvs = new LocusSet<CNVariant>(tmp.toArray(new CNVariant[tmp.size()]),
                                                         true, proj.getLog()) {

        /**
         * 
         */
        private static final long serialVersionUID = 1L;

      };

      int numTotalMarkers = 0;
      for (int i = 0; i < cnvs.getLoci().length; i++) {
        numTotalMarkers += cnvs.getLoci()[i].getNumMarkers();
      }
      int numNonNormalStates = 0;
      for (int element : q) {
        if (element != normalState) {
          // && q[i] != 3
          numNonNormalStates++;
        }
      }
      if (numNonNormalStates != numTotalMarkers) {
        String error = "BUG: detected " + numNonNormalStates
                       + " non-normal states, but collapsed to " + numTotalMarkers + " markers";
        error += "\nSample: " + fid + "\t" + fid;
        proj.getLog().reportTimeError(error);
        for (int i = 0; i < tmp.size(); i++) {
          proj.getLog().reportTimeError(tmp.get(i).toPlinkFormat());

        }
        throw new IllegalStateException(error);
      } else if (verbose) {
        proj.getLog()
            .reportTimeInfo("Found " + cnvs.getLoci().length + " cnvs over " + numTotalMarkers
                            + " total markers covering " + cnvs.getBpCovered()
                            + " bp on chromosome " + currentChr);
      }
      return cnvs;
    }
  }

  private static BStatus loadBstatus(String bTag, BufferedReader reader,
                                     Logger log) throws IOException {
    double[] bmean = loadTwoLineDoubleArray(bTag + "_mean:", reader, log);

    if (bmean != null) {
      double[] bSd = loadTwoLineDoubleArray(bTag + "_sd:", reader, log);
      if (bSd != null) {
        double[] b_uf = loadTwoLineDoubleArray(bTag + "_uf:", reader, log);
        if (b_uf != null) {
          BStatus bStatus = new BStatus(bmean, bSd, b_uf[0]);
          return bStatus;
        }
      }
    }

    return null;
  }

  private static double[] loadTwoLineDoubleArray(String tag, BufferedReader reader,
                                                 Logger log) throws IOException {
    String atag = reader.readLine();
    if (!atag.startsWith(tag)) {
      String error = "cannot read " + tag + " annotation from HMM file, invalid file format";
      log.reportTimeError(error);
      throw new IllegalStateException(error);
    } else {
      String[] line = reader.readLine().trim().split("[\\s]+");
      return Array.toDoubleArray(line);
    }
  }

  private static double[][] loadMatrix(BufferedReader reader, int m, int n) throws IOException {
    double[][] a;
    a = new double[n][m];
    for (int i = 0; i < n; i++) {
      String[] line = reader.readLine().trim().split("[\\s]+");
      double[] tmp = Array.toDoubleArray(line);
      a[i] = tmp;
    }
    return a;
  }

  /**
   * Stores the B1* etc entries for the data distributions of each state
   *
   */
  private static class BStatus {
    private final double[] b_mean;
    private double[] b_sd;
    private final OpdfGaussian[] gaussians;
    private final double b_uf;

    private BStatus(double[] b_mean, double[] b_sd, double b_uf) {
      super();
      this.b_mean = b_mean;
      this.b_sd = b_sd;
      this.b_uf = b_uf;
      gaussians = new OpdfGaussian[b_mean.length];
      for (int i = 0; i < b_mean.length; i++) {
        if (!Numbers.isFinite(b_sd[i]) || !Numbers.isFinite(b_mean[i]) || b_sd[i] <= 0) {
          throw new IllegalArgumentException("Invalid b status mean:" + b_mean[i] + " sd: "
                                             + b_sd[i]);
        }
        gaussians[i] = new OpdfGaussian(b_mean[i], Math.pow(b_sd[i], 2));
      }
    }

    private OpdfGaussian[] getGaussians() {
      return gaussians;
    }

    private double[] getB_mean() {
      return b_mean;
    }

    private double[] getB_sd() {
      return b_sd;
    }

    private double getB_uf() {
      return b_uf;
    }



    @Override
    public String toString() {
      return "BStatus [b_mean=" + Arrays.toString(b_mean) + ", b_sd=" + Arrays.toString(b_sd)
             + ", gaussians=" + Arrays.toString(gaussians) + ", b_uf=" + b_uf + "]";
    }

    /**
     * and create new distributions from the new sds
     */
    private void setB_sd(double[] b_sdToSet) {
      this.b_sd = b_sdToSet;
      for (int i = 0; i < b_mean.length; i++) {
        if (!Numbers.isFinite(b_sd[i]) || !Numbers.isFinite(b_mean[i]) || b_sd[i] <= 0) {
          throw new IllegalArgumentException("Invalid b status mean:" + b_mean[i] + " sd: "
                                             + b_sd[i]);
        }
        gaussians[i] = new OpdfGaussian(b_mean[i], Math.pow(b_sd[i], 2));
      }
    }

  }

  private static PennHmm convertHMMTransition(PennHmm pennHmm, int dist) {
    PennHmm converted = new PennHmm(pennHmm);
    double D = STATE_CHANGE;
    double offdiagonal_sum = 0;
    double[][] tmpA = new double[pennHmm.getA().length][pennHmm.getA()[0].length];
    for (int i = 0; i < converted.getN(); i++) {
      offdiagonal_sum = 0;
      for (int j = 0; j < converted.getN(); j++) {
        if (i != j) {
          if (i == 3) {
            tmpA[i][j] = pennHmm.getA()[i][j] * (1 - Math.exp(-dist / D / 1000))
                         / (1 - Math.exp(-5000 / D / 1000));
          } else {
            tmpA[i][j] =
                       pennHmm.getA()[i][j] * (1 - Math.exp(-dist / D)) / (1 - Math.exp(-5000 / D));
          }
          if (tmpA[i][j] > 1) {
            pennHmm.getLog()
                   .reportTimeWarning("Off-diagonal cell A[%i][%i] (%f to %f by %i) in transition matrix is over boundary of 1 (HMM model is not optimized). Assign 0.999 as the value instead.\n"
                                      + i + "\t" + j + "\t" + pennHmm.getA()[i][j] + "\t"
                                      + tmpA[i][j] + "\t" + dist);
            tmpA[i][j] = 0.999; /*
                                 * maximum possible off-diagonal value (since state3 frequency is
                                 * 0.999)
                                 */
          }
          offdiagonal_sum += tmpA[i][j];
        }
      }
      if (offdiagonal_sum >= 1) {
        for (int j = 0; j < converted.getN(); j++) {
          tmpA[i][j] /= (offdiagonal_sum / 0.999);
        }
        offdiagonal_sum = 0.999;
      }
      tmpA[i][i] = 1 - offdiagonal_sum;
    }
    converted.setA(tmpA);
    return converted;
  }

  public static LocusSet<CNVariant> scoreCNVsSameChr(PennHmm pennHmm, LocusSet<CNVariant> cLocusSet,
                                                     int[] posChr, double[] lrrChr,
                                                     double[] bafsChr, double[] pfbsChr,
                                                     boolean[] copyNumberOnlyDef, int[] q,
                                                     int normalState, Logger log) {
    ArrayList<CNVariant> scored = new ArrayList<CNVariant>();

    for (int i = 0; i < cLocusSet.getLoci().length; i++) {
      CNVariant current = cLocusSet.getLoci()[i];
      ArrayList<Integer> indicestmp = new ArrayList<Integer>();
      for (int j = 0; j < posChr.length; j++) {// TODO, speed up this search
        int pos = posChr[j];
        if (pos >= current.getStart() && pos <= current.getStop()) {
          int stateCN = getStateCN(q[j]);
          if (stateCN == current.getCN()) {// can be different when two markers have identical
                                           // positions and only one of the two supports a cnv start
                                           // or end
            indicestmp.add(j);
          } else {
            if (current.getStart() != pos && current.getStop() != pos) {
              throw new IllegalStateException("Found mismatched states within a cnv");
            }
          }
        }
        if (pos > current.getStop()) {
          break;
        }
      }
      int hmmState = current.getCN();
      if (hmmState > 2) {
        hmmState++;// LOH stored as state 3
      }

      int[] indices = Ints.toArray(indicestmp);
      if (indices.length != current.getNumMarkers()) {
        String error = "BUG: could not reconstruct original markers, found " + indices.length
                       + " and should have found " + current.getNumMarkers() + "Sample FID: "
                       + current.getFamilyID() + " IID: " + current.getIndividualID();
        log.reportTimeError(error);
        throw new IllegalStateException(error);
      }

      double score = getLocScore(pennHmm, Array.subArray(lrrChr, indices),
                                 Array.subArray(bafsChr, indices), Array.subArray(pfbsChr, indices),
                                 Array.subArray(copyNumberOnlyDef, indices), hmmState);
      CNVBuilder builder = new CNVBuilder(current);
      builder.score(score);
      scored.add(builder.build());
    }
    return new LocusSet<CNVariant>(scored.toArray(new CNVariant[scored.size()]), true, log) {

      /**
       * 
       */
      private static final long serialVersionUID = 1L;

    };
  }

  private static int getStateCN(int q) {
    int stateCN = q;
    if (stateCN == 3) {
      stateCN = LOH_FLAG;
    } else if (stateCN > 3) {
      stateCN--;
    }
    return stateCN;
  }

  private static double getLocScore(PennHmm pennHmm, double[] o1, double[] o2, double[] pfb,
                                    boolean[] copyNumberOnlyDef, int actualStateIndex) {
    int[] tests = new int[] {0, 1, 2, 4, 5};
    double confState = Double.NaN;
    double maxOther = -1 * Double.MAX_VALUE;

    for (int test : tests) {
      double tmp = GetStateProb_CHMM(pennHmm, o1, o2, pfb, copyNumberOnlyDef, test);
      if (Numbers.isFinite(tmp)) {
        if (test == actualStateIndex) {
          confState = tmp;
        } else if (tmp > maxOther) {
          maxOther = tmp;
        }
      } else {
        throw new IllegalStateException("Invalid tmp score " + tmp);
      }
    }
    return confState - maxOther;
  }

  private static double GetStateProb_CHMM(PennHmm pennHmm, double[] o1, double[] o2, double[] pfb,
                                          boolean[] copyNumberOnlyDef, int state) {
    double logProb = 0;
    for (int i = 0; i < o1.length; i++) {
      if (copyNumberOnlyDef[i] || o2[i] > 1 || pfb[i] > 1) {
        // TODO PennCNV uses B1 here, why not B3?
        logProb += b1iot(state, pennHmm.getB1(), o1[i]);
      } else {
        logProb += b1iot(state, pennHmm.getB1(), o1[i]);
        logProb += b2iot(state, pennHmm.getB2(), pfb[i], o2[i]);
      }
    }
    return logProb;
  }

  /**
   * Parameterize the {@link PennHmm} by an empirical sd measure
   * 
   * @param pennHmm
   * @param sdo must be strictly positive value
   * @param log
   * @return
   */
  public static PennHmm adjustBSD(PennHmm pennHmm, double sdo, Logger log)
  /*
   * adjust the CHMM model so that the standard deviation of B1 and B2 and B3 match the observed
   * data (by an empirical method)
   */
  {
    if (Double.isNaN(sdo) || sdo <= 0) {
      throw new IllegalArgumentException("Invalid standard deviation " + sdo);
    }
    PennHmm pennAdjusted = new PennHmm(pennHmm);
    double[] b1AdjustedSd = new double[pennHmm.getB1().getB_sd().length];
    double[] b2AdjustedSd = new double[pennHmm.getB2().getB_sd().length];
    double[] b3AdjustedSd = new double[pennHmm.getB3().getB_sd().length];

    double ratio = sdo / pennHmm.getB1().getB_sd()[2];
    log.reportTimeInfo("Adjusting hmm sd measures by " + sdo + " and ratio " + ratio);

    for (int i = 0; i < pennHmm.getN(); i++) {

      b1AdjustedSd[i] = pennHmm.getB1().getB_sd()[i] * ratio;
    }
    for (int i = 0; i < pennHmm.getN() - 1; i++) {
      b2AdjustedSd[i] = pennHmm.getB2().getB_sd()[i] * ratio;// fewer states for baf
    }
    for (int i = 0; i < pennHmm.getN(); i++) {
      b3AdjustedSd[i] = pennHmm.getB3().getB_sd()[i] * ratio;
    }
    pennAdjusted.getB1().setB_sd(b1AdjustedSd);
    pennAdjusted.getB2().setB_sd(b2AdjustedSd);
    pennAdjusted.getB3().setB_sd(b3AdjustedSd);
    return pennAdjusted;
  }
}
