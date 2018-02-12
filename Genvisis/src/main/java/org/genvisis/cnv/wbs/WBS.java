/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;

/**
 * Wild Binary Segmentation http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf
 * https://cran.r-project.org/web/packages/wbs/wbs.pdf Port of
 * https://github.com/cran/wbs/blob/8919e7c35389b92e27b6948572271e0843b5f808/src/wbs.c
 */
public class WBS {

  static final int DEFAULT_M = 5000;
  static final int DEFAULT_SEED = 42;

  private WBS() {

  }

  private static void ipiArgMax(double[] res, int n, IPArg iArg) {

    int i;
    int k;
    int maxCount = 0;
    double tmp, maxFabs;

    iArg.ipargmax = 0;
    maxFabs = -1;

    for (i = 0; i < n - 1; i++) {
      tmp = Math.abs(res[i]);

      if (tmp > maxFabs) {
        iArg.ipargmax = i;
        maxFabs = tmp;
        maxCount = 1;
      } else if (tmp == maxFabs) maxCount++;
    }

    /* if there are multiple points maximizing cusums, we take the median */

    if (maxCount > 1) {
      maxCount = maxCount / 2 + (maxCount % 2);

      k = 0;
      i = 0;

      while ((i < (n - 1)) && (k < maxCount)) {
        i++;
        if (Math.abs(res[i]) == maxFabs) k++;
      }

      iArg.ipargmax = i;
    }

    iArg.ipmax = res[iArg.ipargmax];

  }

  private static void binarySegmentationRecursive(double[] x, int n, int s, int e, double[] res,
                                                  double[] iplus, double[] iminus, double[] ipres,
                                                  double minth, int scale) {

    int len = e - s + 1;
    int cptcand;

    IPArg ipArg = new IPArg();

    if (len > 1) {

      wbsIpi(ArrayUtils.subArray(x, s - 1), len, ipres, iplus, iminus, ipArg);

      cptcand = ipArg.ipargmax + s;

      res[idx(cptcand, 1, n - 1)] = (double) s;
      res[idx(cptcand, 2, n - 1)] = (double) e;
      res[idx(cptcand, 3, n - 1)] = (double) cptcand;
      res[idx(cptcand, 4, n - 1)] = ipArg.ipmax;

      if (minth > Math.abs(ipArg.ipmax) || minth < 0) {
        minth = Math.abs(ipArg.ipmax);
      }

      res[idx(cptcand, 5, n - 1)] = minth;
      res[idx(cptcand, 6, n - 1)] = (double) scale;

      binarySegmentationRecursive(x, n, s, cptcand, res, iplus, iminus, ipres, minth, scale + 1);
      binarySegmentationRecursive(x, n, cptcand + 1, e, res, iplus, iminus, ipres, minth,
                                  scale + 1);

    }

  }

  static void bsRecWrapper(double[] x, int n, double[] res) {

    double[] iplus = new double[n - 1];
    double[] iminus = new double[n - 1];
    double[] ipres = new double[n - 1];

    /* negative value of minth serves as infinity */
    binarySegmentationRecursive(x, n, 1, n, res, iplus, iminus, ipres, -1.0, 1);

  }

  private static void wbsIntegratedRecursive(double[] x, int n, int s, int e, double[] res,
                                             double[] iplus, double[] iminus, double[] ipres,
                                             double[] wbsres, int[] index, int indexn, int M,
                                             double minth, int scale) {
    int len = e - s + 1;
    IPArg ipArg = new IPArg();

    if (len > 1) {
      if (indexn > 0) {
        int cptcand;

        wbsIpi(ArrayUtils.subArray(x, s - 1), len, ipres, iplus, iminus, ipArg);

        if (Math.abs(ipArg.ipmax) < (wbsres[idx(index[0], 5, M)])) {
          cptcand = (int) wbsres[idx(index[0], 3, M)];
          res[idx(cptcand, 1, n - 1)] = (double) s;
          res[idx(cptcand, 2, n - 1)] = (double) e;
          res[idx(cptcand, 3, n - 1)] = (double) cptcand;
          res[idx(cptcand, 4, n - 1)] = wbsres[idx(index[0], 4, M)];

          if (minth > wbsres[idx(index[0], 5, M)] || minth < 0) {
            minth = wbsres[idx(index[0], 5, M)];
          }

        } else {

          cptcand = ipArg.ipargmax + s;
          res[idx(cptcand, 1, n - 1)] = (double) s;
          res[idx(cptcand, 2, n - 1)] = (double) e;
          res[idx(cptcand, 3, n - 1)] = (double) cptcand;
          res[idx(cptcand, 4, n - 1)] = (double) (ipArg.ipmax);

          if (minth > Math.abs(ipArg.ipmax) || minth < 0) {
            minth = Math.abs(ipArg.ipmax);
          }
        }

        res[idx(cptcand, 5, n - 1)] = minth;
        res[idx(cptcand, 6, n - 1)] = (double) scale;

        /* left */
        int[] indexl = new int[indexn];
        int[] indexr = new int[indexn];
        int indexnl = 0;
        int indexnr = 0;
        int i;

        for (i = 1; i <= indexn; i++) {
          if ((wbsres[idx(index[i - 1], 1, M)] >= s)
              && (wbsres[idx(index[i - 1], 2, M)] <= cptcand)) {
            indexl[indexnl] = index[i - 1];
            indexnl++;
          } else if ((wbsres[idx(index[i - 1], 1, M)] >= (cptcand + 1))
                     && (wbsres[idx(index[i - 1], 2, M)] <= e)) {
            indexr[indexnr] = index[i - 1];
            indexnr++;
          }
        }

        if (indexnl != 0) {
          indexl = Arrays.copyOf(indexl, indexnl);
          wbsIntegratedRecursive(x, n, s, cptcand, res, iplus, iminus, ipres, wbsres, indexl,
                                 indexnl, M, minth, scale + 1);
        } else {

          binarySegmentationRecursive(x, n, s, cptcand, res, iplus, iminus, ipres, minth,
                                      scale + 1);
        }

        if (indexnr != 0) {
          indexr = Arrays.copyOf(indexr, indexnr);

          wbsIntegratedRecursive(x, n, cptcand + 1, e, res, iplus, iminus, ipres, wbsres, indexr,
                                 indexnr, M, minth, scale + 1);

        } else {

          binarySegmentationRecursive(x, n, cptcand + 1, e, res, iplus, iminus, ipres, minth,
                                      scale + 1);

        }

      } else {

        binarySegmentationRecursive(x, n, s, e, res, iplus, iminus, ipres, minth, scale);
      }
    }
  }

  /**
   * ipargmax and ipmax are passed by reference in c, so this is a ref wrapper
   */
  private static class IPArg {

    int ipargmax;
    double ipmax;

  }

  private static void wbsIpi(double[] x, int n, double[] res, double[] iplus, double[] iminus,
                             IPArg iArg) {

    double sumx;
    double factor;
    int i;
    double nDbl = (double) (n);
    double oneOverN = 1.0 / nDbl;
    double nSquared = nDbl * nDbl;
    double iDbl;
    double iplusoneInv;

    sumx = 0;
    for (i = 1; i < n; i++) {
      sumx += x[i];
    }

    iminus[0] = 1.0 / Math.sqrt(nSquared - nDbl) * sumx;
    iplus[0] = Math.sqrt(1.0 - oneOverN) * x[0];
    res[0] = iplus[0] - iminus[0];

    for (i = 1; i < n - 1; i++) {
      iDbl = (double) i;
      iplusoneInv = 1.0 / (iDbl + 1.0);

      factor = Math.sqrt((nDbl - iDbl - 1.0) * iDbl * iplusoneInv / (nDbl - iDbl));
      iplus[i] = iplus[i - 1] * factor + x[i] * Math.sqrt(iplusoneInv - oneOverN);
      iminus[i] = iminus[i - 1] / factor - x[i] / Math.sqrt(nSquared * iplusoneInv - nDbl);
      res[i] = iplus[i] - iminus[i];
    }

    ipiArgMax(res, n, iArg);

  }

  private static int idx(int i, int j, int ld) {
    return (((j) - 1) * (ld)) + ((i) - 1);
  }

  /**
   * WBS integrated wrapper - called in trials at
   * https://github.com/PankratzLab/Analysis/blob/master/WBS/WBS_TestClip.md types taken from R code
   * Note for why the int/augmented version The optional augmentation of Ms,e by {0} is done to
   * ensure that the algorithm also exam- ines the entire current interval [s, e], and not only its
   * randomly drawn subintervals, in case [s,e] only contains one change-point and hence it is
   * optimal to examine [s,e] in its entirety. We note that unlike the BS procedure, the WBS
   * algorithm (in the case without the optional augmentation) returns estimated change-points in
   * the order corresponding to decreasing maxima of |X ̃b |, which is due to the maximisation over
   * m. There is no corresponding sm ,em maximisation in the BS procedure, which means that the
   * maxima of the CUSUM statis- tics corresponding to estimated change-points in the latter
   * procedure are not necessarily arranged in decreasing order. <br>
   * http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf
   */
  static List<ChangePoint> wbsIntegratedRecursiveWrapper(double[] x, int[][] intervals,
                                                         Logger log) {
    int M = intervals[0].length;
    if (ArrayUtils.containsMissingValue(x)) {
      throw new IllegalArgumentException("input array cannot contain missing values");
    }
    List<ChangePoint> wList = new ArrayList<>();

    int n = x.length;
    double[] res = new double[6 * x.length];
    double[] iplus = new double[n - 1];
    double[] iminus = new double[n - 1];
    double[] ipres = new double[n - 1];
    double[] wbsres = new double[M * 5];
    int[] index = new int[M];
    IPArg ipArg = new IPArg();
    /* find cpt candidates on given intervals */

    for (int i = 1; i <= M; i++) {
      int start = intervals[0][i - 1];
      int end = intervals[1][i - 1];
      wbsIpi(ArrayUtils.subArray(x, start - 1), end - start + 1, ipres, iplus, iminus, ipArg);
      int cptcand = ipArg.ipargmax + start;
      wbsres[idx(i, 1, M)] = (double) start;
      wbsres[idx(i, 2, M)] = (double) end;
      wbsres[idx(i, 3, M)] = (double) cptcand;
      wbsres[idx(i, 4, M)] = ipArg.ipmax;
      wbsres[idx(i, 5, M)] = Math.abs(ipArg.ipmax);
      index[i - 1] = i;
    }
    /* sort elementd from the one with the largest abs(cusum) */

    index = getSortedIndex(M, wbsres, index);

    /* standard wbs part */

    wbsIntegratedRecursive(x, n, 1, n, res, iplus, iminus, ipres, wbsres, index, M, M, -1.0, 1);

    int space = res.length / 6;
    for (int j = 0; j < space - 1; j++) {// Last entry is inverse repeat of first
      wList.add(new ChangePoint((int) res[j], (int) res[j + space - 1],
                                (int) res[j + 2 * (space - 1)], res[j + 3 * (space - 1)],
                                res[j + 4 * (space - 1)], (int) res[j + 5 * (space - 1)]));

    }
    log.reportTimeInfo(M + " Random intervals giving " + space + " possible change point(s)");
    return wList;

  }

  private static void wbsRecursive(double[] x, int n, int s, int e, double[] res, double[] wbsres,
                                   int[] index, int indexn, int M, int scale) {

    int len = e - s + 1;

    if (len > 1) {

      if (indexn > 0) {

        int cptcand;

        cptcand = (int) wbsres[idx(index[0], 3, M)];

        res[idx(cptcand, 1, n - 1)] = wbsres[idx(index[0], 1, M)];
        res[idx(cptcand, 2, n - 1)] = wbsres[idx(index[0], 2, M)];
        res[idx(cptcand, 3, n - 1)] = cptcand;
        res[idx(cptcand, 4, n - 1)] = wbsres[idx(index[0], 4, M)];
        res[idx(cptcand, 5, n - 1)] = wbsres[idx(index[0], 5, M)];
        res[idx(cptcand, 6, n - 1)] = (double) scale;

        /* left */
        int[] indexl = new int[indexn];
        int[] indexr = new int[indexn];
        int indexnl = 0;
        int indexnr = 0;
        int i;

        for (i = 1; i <= indexn; i++) {
          if ((wbsres[idx(index[i - 1], 1, M)] >= s)
              && (wbsres[idx(index[i - 1], 2, M)] <= cptcand)) {
            indexl[indexnl] = index[i - 1];
            indexnl++;
          } else if ((wbsres[idx(index[i - 1], 1, M)] >= (cptcand + 1))
                     && (wbsres[idx(index[i - 1], 2, M)] <= e)) {
            indexr[indexnr] = index[i - 1];
            indexnr++;
          }
        }

        if (indexnl != 0) {
          indexl = Arrays.copyOf(indexl, indexnl);
          wbsRecursive(x, n, s, cptcand, res, wbsres, indexl, indexnl, M, scale + 1);
        }

        if (indexnr != 0) {
          indexr = Arrays.copyOf(indexr, indexnr);
          wbsRecursive(x, n, cptcand + 1, e, res, wbsres, indexr, indexnr, M, scale + 1);
        }

      }
    }
  }

  static void wbsRecursiveWrapper(double[] x, int n, double[] res, int[] intervals, int M) {
    if (x == null || x.length >= 0) {
      throw new UnsupportedOperationException();
    }
    double[] iplus = new double[n - 1];
    double[] iminus = new double[n - 1];
    double[] ipres = new double[n - 1];
    double[] wbsres = new double[M * 5];
    int[] index = new int[M];
    int i;
    int s;
    int e;
    int cptcand;

    IPArg ipArg = new IPArg();
    /* find cpt candidates on given intervals */

    for (i = 1; i <= M; i++) {

      s = intervals[idx(i, 1, M)];
      e = intervals[idx(i, 2, M)];

      wbsIpi(ArrayUtils.subArray(x, s - 1), e - s + 1, ipres, iplus, iminus, ipArg);
      cptcand = ipArg.ipargmax + s;

      wbsres[idx(i, 1, M)] = (double) s;
      wbsres[idx(i, 2, M)] = (double) e;
      wbsres[idx(i, 3, M)] = (double) cptcand;
      wbsres[idx(i, 4, M)] = ipArg.ipmax;
      wbsres[idx(i, 5, M)] = Math.abs(ipArg.ipmax);
      index[i - 1] = i;

    }
    /* sort elementd from the one with the largest abs(cusum) */
    index = getSortedIndex(M, wbsres, index);

    /* standard wbs part */

    wbsRecursive(x, n, 1, n, res, wbsres, index, M, M, 1);

  }

  /**
   * Trying to accomplish this
   * https://github.com/cran/wbs/blob/8919e7c35389b92e27b6948572271e0843b5f808/src/wbs.c#L248-L254
   * 
   * @param M
   * @param wbsres
   * @param index
   * @return
   */
  private static int[] getSortedIndex(int M, double[] wbsres, int[] index) {
    double[] tmp = new double[M];

    System.arraycopy(wbsres, idx(1, 5, M), tmp, 0, M);

    int[] indices = Sort.getSortedIndices(tmp);
    return ArrayUtils.reverse(Sort.getOrdered(index, indices));
  }
}
