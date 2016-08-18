package org.genvisis.cnv.hmm;

// package cnv.hmm;
//
// import java.util.ArrayList;
// import java.util.Arrays;
// import java.util.Hashtable;
// import java.util.Random;
//
// import org.apache.commons.math3.distribution.HypergeometricDistribution;
// import org.apache.commons.math3.distribution.NormalDistribution;
// import org.apache.commons.math3.random.MersenneTwister;
// import org.apache.commons.math3.util.CombinatoricsUtils;
//
// import common.Array;
// import common.Sort;
//
/// **
// * @author lane0212 Utility class for the circular binary segmentation. Ported from C ->
// https://github.com/Illumina/canvas/ which was in turn ported from R/fortran ala
// https://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html
// *
// */
// public class CBSUtils {
// private static final double M_LN2 = 0.693147180559945309417;
// private static final double DBL_EPSILON = 2.2204460492503131E-16; // dangerous...the smallest x
// such that 1 + x != 1
// private static final double M_LN_SQRT_2PI = 0.918938533204672741780329736406;
// private static final double M_LOG10_2 = 0.301029995663981195213738894724; // log10(2)
//
// private static double PExceed(int nPerm, int n1s, int[] sbdry, int sbdryOffset, double pExcd) {
// double dlcnk;
// int i, n, k, n1, k1, n2, k2, n3, k3;
//
// n = nPerm;
// k = n1s;
// n1 = nPerm - sbdry[sbdryOffset];
// dlcnk = CombinatoricsUtils.binomialCoefficientLog(n, k);
// pExcd = Math.exp(CombinatoricsUtils.binomialCoefficientLog(n1, k) - dlcnk);
// if (n1s >= 2) {
// n1 = sbdry[sbdryOffset];
// n = nPerm - sbdry[sbdryOffset + 1];
// k = n1s - 1;
// pExcd += Math.exp(Math.log(n1) + CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk);
// }
// if (n1s >= 3) {
// n1 = sbdry[sbdryOffset];
// n2 = sbdry[sbdryOffset + 1];
// n = nPerm - sbdry[sbdryOffset + 2];
// k = n1s - 2;
// pExcd += Math.exp(Math.log(n1) + Math.log(n1 - 1.0) - Math.log(2.0) +
// CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk) + Math.exp(Math.log(n1) + Math.log(n2 -
// n1) + CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk);
// }
// if (n1s > 3) {
// for (i = 4; i <= n1s; i++) {
// n1 = sbdry[sbdryOffset + i - 4];
// k1 = i - 1;
// k2 = i - 2;
// k3 = i - 3;
// n2 = sbdry[sbdryOffset + i - 3];
// n3 = sbdry[sbdryOffset + i - 2];
// n = nPerm - sbdry[sbdryOffset + i - 1];
// k = n1s - i + 1;
// pExcd += Math.exp(CombinatoricsUtils.binomialCoefficientLog(n1, k1) +
// CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk) +
// Math.exp(CombinatoricsUtils.binomialCoefficientLog(n1, k2) + Math.log(n3 - n1) +
// CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk) +
// Math.exp(CombinatoricsUtils.binomialCoefficientLog(n1, k3) + Math.log(n2 - n1) + Math.log(n3 -
// n2) + CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk) +
// Math.exp(CombinatoricsUtils.binomialCoefficientLog(n1, k3) + Math.log(n2 - n1) - Math.log(2.0) +
// Math.log(n2 - n1 - 1.0) + CombinatoricsUtils.binomialCoefficientLog(n, k) - dlcnk);
// }
// }
// return pExcd;
// }
//
// public enum SegmentSplitUndo {
// None, Prune, SDUndo
// }
//
// public enum SegmentationMethod {
// Wavelets, CBS
// }
//
// // / <summary>
// // / Contains static funtions found in R/changepoints.R and src/changepoints.f
// // / </summary>
//
// // / <summary>
// // / Outputs:
// // / lengthSeg
// // / segmentMeans
// // / </summary>
// // / <param name="genomeData"></param>
// // / <param name="sbdry"></param>
// // / <param name="lengthSeg">segment lengths</param>
// // / <param name="segmentMeans">segment means</param>
// // / <param name="dataType">"logratio" or "binary"</param>
// // / <param name="alpha"></param>
// // / <param name="nPerm"></param>
// // / <param name="pMethod"></param>
// // / <param name="minWidth"></param>
// // / <param name="kMax"></param>
// // / <param name="nMin"></param>
// // / <param name="trimmedSD"></param>
// // / <param name="undoSplits">"none" or "prune" or "sdundo"</param>
// // / <param name="undoPrune"></param>
// // / <param name="undoSD"></param>
// // / <param name="verbose"></param>
// // / <param name="nGrid"></param>
// // / <param name="tol"></param>
// public static int[] ChangePoints(double[] genomeData, int[] sbdry, int[] lengthSeg, double[]
// segmentMeans, MersenneTwister rnd, String dataType, double alpha, int nPerm, String pMethod, int
// minWidth, int kMax, int nMin, double trimmedSD, SegmentSplitUndo undoSplits, double undoPrune,
// double undoSD, int verbose, int nGrid, double tol) {
// // int n = genomeData.Length;
// if (trimmedSD <= 0) {
//
// trimmedSD = Array.mad(Array.Diff(genomeData, 1), 1.4826) / Math.sqrt(2);
// }
// // start with the whole
// ArrayList<Integer> segEnd = new ArrayList<Integer>();
// segEnd.add(0); // inclusive
// segEnd.add(genomeData.length); // exclusive
// int k = segEnd.size();
// ArrayList<Integer> changeLocations = new ArrayList<Integer>();
// int nChangePoints = 0;
// int[] iChangePoint = null;
// while (k > 1) {
// int currentN = segEnd.get(k - 1) - segEnd.get(k - 2);
// if (verbose >= 3) {
// // Console.Write(".... current segment: {0} - {1} \n", segEnd[k - 2] + 1, segEnd[k - 1]);
// }
// if (currentN >= 2 * minWidth) {
// double[] currentGenomeData = new double[currentN];
// java.lang.System.arraycopy(genomeData, segEnd.get(k - 2), currentGenomeData, 0, currentN);
// // check whether hybrid method needs to be used
// boolean hybrid = false;
// double delta = 0.0;
// if (pMethod.equals("hybrid") && nMin < currentN) {
// hybrid = true;
// delta = (kMax + 1.0) / currentN;
// }
//
// // if all values of current.genomdat are the same don't segment
// if (Array.max(currentGenomeData) == Array.min(currentGenomeData)) {
// nChangePoints = 0;
// } else {
// // centering the current data will save a lot of computations later
// double currentAverage = Array.mean(currentGenomeData);
// currentGenomeData = Array.InplaceSub(currentGenomeData, currentAverage);
// // need total sum of squares too
// double currentTSS = Array.WeightedSumOfSquares(currentGenomeData, null);
// iChangePoint = FindChangePoints(currentGenomeData, currentTSS, nPerm, alpha, nChangePoints,
// iChangePoint, dataType.equals("binary"), hybrid, minWidth, kMax, delta, nGrid, sbdry, tol, rnd);
// }
// } else {
// nChangePoints = 0;
// }
// // Save the change location
// // segEnd[k - 1] will be removed when nChangePoints == 0
// if (nChangePoints == 0) {
// changeLocations.add(segEnd.get(k - 1));
// }
// // Offset iChangePoint by segEnd[k - 2]
// for (int i = 0; i < nChangePoints; i++) {
// iChangePoint[i] += segEnd.get(k - 2);
// }
// switch (nChangePoints) // switch by the number of change points
// {
// case 0: // no change point
// segEnd.remove(k - 1); // Remove the last element
// break;
// case 1: // one change point
// segEnd.add(k - 1, iChangePoint[0]);
// break;
// case 2: // two change points
// // TODO altered
// ArrayList<Integer> tmp = new ArrayList<Integer>();
// for (int i = 0; i < iChangePoint.length; i++) {
// tmp.add(iChangePoint[i]);
// }
// segEnd.addAll(k - 1, tmp);
// break;
// default:
// System.err.println("There should be 0, 1, or 2 change points");
// break;
// }
// k = segEnd.size();
// if (verbose >= 3) {
// System.err.println(".... segments to go: {0} \n");
// }
// }
// ArrayList<Integer> tmpR = new ArrayList<Integer>();
// int rIndex = changeLocations.size() - 1;
// for (int i = 0; i < changeLocations.size(); i++) {
// tmpR.add(changeLocations.get(rIndex));
// rIndex--;
// }
// // TODO altered
// // changeLocations.Reverse(); // changeLocations is no longer needed
// ArrayList<Integer> segEnds = changeLocations;
// int nSeg = segEnds.size();
// segEnds.add(0, 0);
// lengthSeg = Array.Diff(Array.toIntArray(segEnds), 1);
//
// if (nSeg > 1) {
// if (undoSplits == SegmentSplitUndo.Prune) {
// lengthSeg = ChangePointsPrune(genomeData, lengthSeg, undoPrune);
// }
// if (undoSplits == SegmentSplitUndo.SDUndo) {
// lengthSeg = ChangePointsSDUndo(genomeData, lengthSeg, trimmedSD, undoSD);
// }
// }
// segmentMeans = new double[lengthSeg.length];
// int ll = 0, uu = 0;
// for (int i = 0; i < lengthSeg.length; i++) {
// uu += lengthSeg[i];
// // Works even if weights == null
// segmentMeans[i] = WeightedAverage(genomeData, null, ll, uu);
// ll = uu;
// }
// return lengthSeg;
// }
//
// private static int[] ChangePointsSDUndo(double[] genomeData, int[] lengthSeg, double trimmedSD,
// double changeSD) {
// changeSD *= trimmedSD;
//
// int[] changePointLocationstmp = Array.CumulativeSum(lengthSeg);
// ArrayList<Integer> changePointLocations = new ArrayList<Integer>();
// for (int i = 0; i < changePointLocationstmp.length; i++) {
// changePointLocations.add(changePointLocationstmp[i]);
// }
// boolean sdUndo = true;
// while (sdUndo) {
// int k = changePointLocations.size();
// if (k > 1) {
// ArrayList<Integer> starts = new ArrayList<Integer>(); // make a copy of changePointLocations
// for (int i = 0; i < changePointLocations.size(); i++) {
// starts.add(changePointLocations.get(i));
// }
// starts.remove(k - 1); // Remove the element at k - 1
// starts.add(0, 0); // Insert 0 to the front of tmp
// // starts will be used as the start indices (0-based, inclusive) ==> no need to add 1
// // changePointLocations will be used as the end indices (0-based, exclusive)
// double[] segmentMedians = new double[k];
// for (int i = 0; i < k; i++) // for each segment
// {
// segmentMedians[i] = Array.median(Array.subArray(genomeData, starts.get(i),
// changePointLocations.get(i)));
//
// }
// double[] absDiffSegmentMedians = new double[segmentMedians.length];
// for (int i = 0; i < absDiffSegmentMedians.length; i++) {
// absDiffSegmentMedians[i] = Math.abs(segmentMedians[i]);
// }
//
// double min = 0;
// int iMin = ArgMin(absDiffSegmentMedians, min);
// if (min < changeSD) {
// changePointLocations.remove(iMin);
// } else {
// sdUndo = false;
// }
// } else {
// sdUndo = false;
// }
// }
// changePointLocations.add(0, 0);
// return Array.Diff(Array.toIntArray(changePointLocations), 1); // segment lengths
// }
//
// // / <returns>argmin_i x[i]</returns>
// public static int ArgMin(double[] x, double min) {
// if (x == null) {
// min = Double.NaN;
// return -1;
// }
// min = x[0];
// int iMin = 0;
// for (int i = 1; i < x.length; i++) {
// if (x[i] < min) {
// min = x[i];
// iMin = i;
// }
// }
// return iMin;
// }
//
// // / <summary>
// // / R function changepoints.prune, which calls Fortran subroutine prune
// // / </summary>
// // / <param name="genomeData"></param>
// // / <param name="lengthSeg"></param>
// // / <param name="changeCutoff"></param>
// // / <returns></returns>
// private static int[] ChangePointsPrune(double[] genomeData, int[] lengthSeg, double changeCutoff)
// {
// double[] sx = new double[lengthSeg.length]; // segment sums
// int[] loc = new int[lengthSeg.length - 1]; // one element for each change point
// int[][] loc1 = new int[2][lengthSeg.length - 1]; // one element for each change point
// int prunedNChangePoints = 0; // number of change points after pruning
//
// int i, j, k, kmj; // kmj: k minus j
// double ssq, wssqk, wssq1, wssqj;
// boolean jleft;
//
// ssq = Array.PartialSumOfPowers(genomeData, 2, 0, genomeData.length); // sum of squares
// k = 0; // segment start index
// for (i = 0; i < lengthSeg.length; i++) // for each segment
// {
// sx[i] = Array.PartialSumOfPowers(genomeData, 1, k, lengthSeg[i]);
// k += lengthSeg[i];
// }
// // k = nseg - 1 == number of change points
// for (i = 0; i < loc.length; i++) {
// loc[i] = i + 1;
// loc1[1][i] = i + 1;
// }
// wssqk = ssq - ErrorSumOfSquares(lengthSeg, sx, loc.length, loc);
// // j: number of change points
// for (j = loc.length - 1; j > 0; j--) // j (= loc.Length - 1, ..., 1) is not an index
// {
// kmj = loc.length - j;
// jleft = true;
// for (i = 0; i < j; i++) {
// loc[i] = i + 1;
// loc1[0][i] = i + 1;
// }
// wssqj = ssq - ErrorSumOfSquares(lengthSeg, sx, j, loc);
// while (jleft) {
// Combination(j, kmj, loc, jleft);
// // TODO modified
// wssq1 = ssq - ErrorSumOfSquares(lengthSeg, sx, j, loc);
// if (wssq1 <= wssqj) {
// wssqj = wssq1;
// for (i = 0; i < j; i++) {
// loc1[0][i] = loc[i];
// }
// }
// }
// if (wssqj / wssqk > 1 + changeCutoff) {
// prunedNChangePoints = j + 1; // go back to the previous j
// for (i = 0; i < prunedNChangePoints; i++) {
// loc[i] = loc1[1][i];
// }
// break;
// } else {
// for (i = 0; i < j; i++) {
// loc1[1][i] = loc1[0][i];
// }
// }
// }
// int[] cumSumLengthSeg = Array.CumulativeSum(lengthSeg);
// int[] prunedChangePoints = new int[prunedNChangePoints + 2];
// for (i = 0; i < prunedNChangePoints; i++) { // loc[i] is an 1-based index
// prunedChangePoints[i + 1] = cumSumLengthSeg[loc[i] - 1];
// }
// prunedChangePoints[0] = 0;
// prunedChangePoints[prunedChangePoints.length - 1] = genomeData.length;
// return Array.Diff(prunedChangePoints, 1);
// }
//
// // / <summary>
// // / Fortran subroutine fndcpt.
// // / Ternary segmentation with permutation reference distribution
// // / </summary>
// // / <param name="genomeData"></param>
// // / <param name="totalSumOfSquares"></param>
// // / <param name="nPerm"></param>
// // / <param name="cutoffPValue"></param>
// // / <param name="nChangePoints"></param>
// // / <param name="iChangePoint"></param>
// // / <param name="isBinary"></param>
// // / <param name="hybrid"></param>
// // / <param name="al0"></param>
// // / <param name="hk"></param>
// // / <param name="delta"></param>
// // / <param name="nGrid"></param>
// // / <param name="sbdry"></param>
// // / <param name="tol"></param>
//
// // private static class ChangePoints {
// // private int nChangePoints;
// // private int[] iChangePoint;
// // }
// private static int[] FindChangePoints(double[] genomeData, double totalSumOfSquares, int nPerm,
// double cutoffPValue, int nChangePoints, int[] iChangePoint, boolean isBinary, boolean hybrid, int
// al0, int hk, double delta, int nGrid, int[] sbdry, double tol, MersenneTwister rnd) {
// double[] px = new double[genomeData.length]; // permuted genomeData
// double[] sx = new double[genomeData.length];
// iChangePoint = new int[2]; // up to 2 change points
//
// // nrej: # of non-rejected tests
// // nrejc: # of non-rejected tests cutoff
// int np, nrej, nrejc, n1, n2, n12, l, k;
// int[] iseg = new int[2];
// // segment lengths: iseg[0], iseg[1] - iseg[0], genomeData.Length - iseg[1]
// double ostat = 0;
// double ostat1, pstat, tPValue, pValue1, pValue2;
//
// nrej = 0;
// nChangePoints = 0;
//
// TMaxO(genomeData, totalSumOfSquares, sx, iseg, ostat, al0, isBinary);
// ostat1 = Math.sqrt(ostat);
// ostat *= 0.99999;
// // if maximal t-statistic is too small (for now use 0.1) don't split
// if (ostat1 <= 0.1) {
// return iChangePoint;
// } // call rndend() before return?
// // if maximal t-statistic is too large (for now use 7.0) split
// // also make sure it's not affected by outliers i.e. small seglength
// l = Math.min(iseg[1] - iseg[0], genomeData.length - iseg[1] + iseg[0]);
// if (!((ostat1 >= 7.0) && (l >= 10))) {
// // o.w calculate p-value and decide if & how data are segmented
// if (hybrid) {
// pValue1 = TailP(ostat1, delta, genomeData.length, nGrid, tol);
// if (pValue1 > cutoffPValue) {
// return iChangePoint;
// } // pValue1 is the lower bound
// pValue2 = cutoffPValue - pValue1;
// nrejc = (int) (pValue2 * nPerm);
// k = nrejc * (nrejc + 1) / 2 + 1;
// for (np = 1; np <= nPerm; np++) {
// XPerm(genomeData, px, rnd);
// pstat = HTMaxP(hk, totalSumOfSquares, px, sx, al0, isBinary);
// if (ostat <= pstat) {
// nrej++;
// k++;
// }
// if (nrej > nrejc) {
// return iChangePoint;
// }
// if (np >= sbdry[k - 1]) {
// break;
// }
// }
// } else {
// nrejc = (int) (cutoffPValue * nPerm);
// k = nrejc * (nrejc + 1) / 2 + 1;
// for (np = 1; np <= nPerm; np++) {
// XPerm(genomeData, px, rnd);
// pstat = TMaxP(totalSumOfSquares, px, sx, al0, isBinary);
// if (ostat <= pstat) {
// nrej++;
// k++;
// }
// if (nrej > nrejc) {
// return iChangePoint;
// }
// if (np >= sbdry[k - 1]) {
// break;
// }
// }
// }
// }
// // 200
// if (iseg[1] == genomeData.length) // The second change point is the right boundary
// {
// nChangePoints = 1;
// iChangePoint[0] = iseg[0];
// } else {
// if (iseg[0] == 0) // The first change point is the left boundary
// {
// nChangePoints = 1;
// iChangePoint[0] = iseg[1];
// } else {
// l = 0;
// n1 = iseg[0];
// n12 = iseg[1];
// n2 = n12 - n1;
// // |-- n1 = iseg[0] --|-- n2 = n12 - n1 --|
// // |-- n12 = iseg[1] --|
// tPValue = TPermP(n1, n2, n12, genomeData, l, px, nPerm, rnd);
// if (tPValue <= cutoffPValue) {
// nChangePoints = 1;
// iChangePoint[0] = iseg[0];
// }
// l = iseg[0];
// n12 = genomeData.length - iseg[0];
// n2 = genomeData.length - iseg[1];
// n1 = n12 - n2;
// // |-- n1 = n12 - n2 --|-- n2 = n - iseg[1] --|
// // |-- n12 = n - iseg[0] --|
// tPValue = TPermP(n1, n2, n12, genomeData, l, px, nPerm, rnd);
// if (tPValue <= cutoffPValue) {
// nChangePoints++;
// iChangePoint[nChangePoints - 1] = iseg[1];
// }
// }
// }
// return iChangePoint;
// // 500
// }
//
// // / <summary>
// // / Permute x and store the results in px
// // / </summary>
// // / <param name="x">the array to be permuted</param>
// // / <param name="px">the permuted array</param>
// private static void XPerm(double[] x, double[] px, MersenneTwister rnd) {
// int i, j;
// double cc;
//
// for (i = 0; i < x.length; i++) {
// px[i] = x[i];
// }
//
// for (i = x.length - 1; i >= 0; i--) {
// cc = rnd.nextDouble();
// j = (int) (cc * (i + 1));
// j = (j > i) ? i : j;
// Swap(px[i], px[j]);
// }
// }
//
// public static double TrimmedVariance(Hashtable<String, double[]> scoresByChr, double trim) {
// int n = 0;
// for (String chr : scoresByChr.keySet()) {
// n += scoresByChr.get(chr).length;
// }
// double[] diff = new double[n - 1];
// int i = 0;
// double last = Double.NaN;
// for (String chr : scoresByChr.keySet()) {
// if (scoresByChr.get(chr).length <= 0) {
// continue;
// }
// if (i > 0) {
// diff[i] = scoresByChr.get(chr)[0] - last;
// i++;
// }
//
// java.lang.System.arraycopy(Array.Diff(scoresByChr.get(chr), 1), 0, diff, i,
// scoresByChr.get(chr).length - 1);
// i += (scoresByChr.get(chr).length - 1);
// last = scoresByChr.get(chr)[scoresByChr.get(chr).length - 1];
// }
// int nKeep = (int) Math.round((1 - 2 * trim) * (n - 1));
// // R code: inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))
//
// Array.InplaceAbs(diff);
// Arrays.sort(diff);
// return InflationFactor(trim) * Array.PartialSumOfPowers(diff, 2, 0, nKeep) / (2 * nKeep);
// }
//
// // / <summary>
// // / Truncated N(0, 1) to (x1, x2), where P(X <= x1) = trim and P(X <= x2) = 1 - trim
// // / </summary>
// // / <param name="trim">tail probability</param>
// // / <returns>approximately (E[X^2] = 1) / Et[X^2]</returns>
// private static double InflationFactor(double trim) {
// NormalDistribution norm = new NormalDistribution(); // N(0, 1)
// double a = norm.inverseCumulativeProbability(1 - trim);
// double step = 2 * a / 10000;
// double[] x1s = Seq(-a + step / 2, a - step / 2, 10000);
// // Truncated N(0, 1) to P(X <= x) = trim and P(X <= x) = 1 - trim
// double eX2 = 0.0;
// for (int i = 0; i < x1s.length; i++) {
// eX2 += (x1s[i] * x1s[i]) * norm.density(x1s[i]);
// }
// eX2 = eX2 * step / (1 - 2 * trim);
// // eX2 now approximates Et[X^2]: E[X^2] of the truncated N(0, 1)
// return 1 / eX2; // approx (E[X^2] = 1) / Et[X^2]
// // == 1 / (1 + (-a * dnorm(-a) - a * dnorm(a)) / (1 - 2*trim) - ((dnorm(-1) - dnorm(a)) / (1 - 2*
// trim))^2)
// // According to http://en.wikipedia.org/wiki/Truncated_normal_distribution
// }
//
// public static void TMaxO(double[] genomeData, double totalSumOfSquares, double[] sx, int[] iseg,
// double ostat, int al0, boolean isBinary) {
// // look at the partial sums in blocks of size sqrt(n)
// int ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1, nb2, bi, bj, ilo, ihi, jlo, jhi,
// alenmax, i2j, sxmxi, alenlo, alenhi, tmaxi, tmaxj, ixlo, ixhi, nal0;
// double psum, psmin, psmax, psmin0, psmax0, bssmax, bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1,
// sij2, sijmx0, absx, sxmx, bijbss, rnov2, psdiff;
//
// // use local arrays for working within blocks
// // block partial sum max and min
// double[] bpsmax, bpsmin;
// // location of the max and min
// // bb: block boundary
// int[] bb, ibmin, ibmax;
//
// // t statistic corresponding to max for block i,j (and max possible)
// // double precision, allocatable :: bssbij(:), bssijmax(:)
// double[] bssbij, bssijmax;
// // row, column and order vector for reordering bssbij
// int[] bloci, blocj, loc, alen;
//
// // calculate number of blocks (nb) and block boundaries (vector bb)
// rn = genomeData.length;
// if (genomeData.length >= 50) {
// nb = (int) Math.round(Math.sqrt(genomeData.length));
// } else {
// nb = 1;
// }
// // the number of pairwise block comparison
// nb2 = nb * (nb + 1) / 2;
// // allocate memory
// bpsmax = new double[nb];
// bpsmin = new double[nb];
// bb = new int[nb];
// ibmin = new int[nb];
// ibmax = new int[nb];
// bssbij = new double[nb2];
// bssijmax = new double[nb2];
// bloci = new int[nb2];
// blocj = new int[nb2];
// loc = new int[nb2];
// alen = new int[nb2];
//
// // block boundaries
// for (i = 0; i < nb; i++) {
// bb[i] = (int) Math.round(rn * ((i + 1.0) / nb));
// }
// ilo = 1; // not just an index
// psum = 0;
// psmin0 = 0;
// psmax0 = 0;
// ipsmin0 = genomeData.length;
// ipsmax0 = genomeData.length;
// for (j = 0; j < nb; j++) // j is an index and only an index
// {
// sx[ilo - 1] = psum + genomeData[ilo - 1];
// psmin = sx[ilo - 1];
// ipsmin = ilo;
// psmax = sx[ilo - 1];
// ipsmax = ilo;
// for (i = ilo + 1; i <= bb[j]; i++) {
// sx[i - 1] = sx[i - 2] + genomeData[i - 1];
// if (sx[i - 1] < psmin) {
// psmin = sx[i - 1];
// ipsmin = i;
// }
// if (sx[i - 1] > psmax) {
// psmax = sx[i - 1];
// ipsmax = i;
// }
// }
// // store the block min, max and locations
// ibmin[j] = ipsmin;
// ibmax[j] = ipsmax;
// bpsmin[j] = psmin;
// bpsmax[j] = psmax;
// // adjust global min, max and locations
// if (psmin < psmin0) {
// psmin0 = psmin;
// ipsmin0 = ipsmin;
// }
// if (psmax > psmax0) {
// psmax0 = psmax;
// ipsmax0 = ipsmax;
// }
// // reset ilo to be the block boundary + 1
// psum = sx[bb[j] - 1];
// ilo = bb[j] + 1;
// }
// // calculate bss for max s_i - min s_i
// psdiff = psmax0 - psmin0;
// rj = Math.abs(ipsmax0 - ipsmin0);
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssmax = rnjov1 * Math.pow((psdiff - 0.5), 2);
// } else {
// bssmax = rnjov1 * Math.pow(psdiff, 2);
// }
// tmaxi = Math.min(ipsmax0, ipsmin0);
// tmaxj = Math.max(ipsmax0, ipsmin0);
// // if the segment is all constant then psdiff = 0 and so bssmax = 0
// if (psdiff <= 0) {
// bssmax = 0;
// // go to 120
// } else {
// // for a pair of blocks (i,j) calculate the max absolute t-statistic
// // at the (min_i, max_j) and (max_i, min_j) locations
// // for other indices the t-statistic can be bounded using this
// //
// // if a block doesn't have the potential to exceed bssmax ignore it
// // calculate the bsslim for each block and include ones >= bssmax
// rnov2 = rn / 2;
// l = 0; // l is an index and counter in the following nested for loop
// nal0 = genomeData.length - al0;
// for (i = 1; i <= nb; i++) // i is not just an 1-based index
// {
// for (j = i; j <= nb; j++) // j is not just an 1-based index
// {
// // calculate bsslim
// ilo = (i == 1) ? 1 : bb[i - 2] + 1;
// ihi = bb[i - 1];
// jlo = (j == 1) ? 1 : bb[j - 2] + 1;
// jhi = bb[j - 1];
// alenhi = jhi - ilo;
// if (alenhi > nal0) {
// alenhi = nal0;
// }
// rjhi = alenhi;
// alenlo = (i == j) ? 1 : jlo - ihi;
// if (alenlo < al0) {
// alenlo = al0;
// }
// // max S_k over block j - min S_k over block i
// sij1 = Math.abs(bpsmax[j - 1] - bpsmin[i - 1]);
// // max S_k over block i - min S_k over block j
// sij2 = Math.abs(bpsmax[i - 1] - bpsmin[j - 1]);
// // if i = j then sij1 and sij2 are the same
// sijmx0 = Math.max(sij1, sij2);
// rjlo = alenlo;
// rnjov1 = rn / Math.min(rjlo * (rn - rjlo), rjhi * (rn - rjhi));
// if (isBinary) {
// bsslim = rnjov1 * Math.pow(sijmx0 - 0.5, 2);
// } else {
// bsslim = rnjov1 * Math.pow(sijmx0, 2);
// }
// // if its as large as bssmax add block
// if (bssmax <= bsslim) {
// loc[l] = l + 1;
// bloci[l] = i;
// blocj[l] = j;
// bssijmax[l] = bsslim;
// // max sij in the (i,j) block, t-statistic etc
// if (sij1 > sij2) {
// alen[l] = Math.abs(ibmax[j - 1] - ibmin[i - 1]);
// rj = alen[l];
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssbij[l] = rnjov1 * Math.pow(sij1 - 0.5, 2);
// } else {
// bssbij[l] = rnjov1 * Math.pow(sij1, 2);
// }
// } else {
// alen[l] = Math.abs(ibmin[j - 1] - ibmax[i - 1]);
// rj = alen[l];
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssbij[l] = rnjov1 * Math.pow(sij2 - 0.5, 2);
// } else {
// bssbij[l] = rnjov1 * Math.pow(sij2, 2);
// }
// }
// l++;
// }
// }
// }
// nb1 = l;
// // Now sort the t-statistics by their magnitude
// // TODO, same functionality?
//
// int[] indices = Sort.quicksort(bssbij);
// int[] tLoc = new int[loc.length];
// for (int m = 0; m < tLoc.length; m++) {
// tLoc[m] = loc[indices[m]];
// }
// loc = tLoc;
// Arrays.sort(bssbij);
// // now go through the blocks in reverse order (largest down)
// for (l = nb1 - 1; l >= 0; l--) // l is an index
// {
// if (loc[l] - 1 >= 0) {
// k = loc[l] - 1; // k is an index in the for loop
// // need to check a block only if it has potential to increase bss
// // rjlo is the smalllest (j-i) in the block and rjhi is the largest
// bsslim = bssijmax[k];
// if (bssmax <= bsslim) {
// // bi, bj give the block location
// bi = bloci[k];
// bj = blocj[k];
// // max arc length of interest in block
// alenmax = alen[k];
// ilo = (bi == 1) ? 1 : bb[bi - 2] + 1;
// ihi = bb[bi - 1];
// jlo = (bj == 1) ? 1 : bb[bj - 2] + 1;
// jhi = bb[bj - 1];
// alenhi = jhi - ilo;
// if (alenhi > nal0) {
// alenhi = nal0;
// }
// rjhi = alenhi;
// alenlo = (bi == bj) ? 1 : (jlo - ihi);
// if (alenlo < al0) {
// alenlo = al0;
// }
// rjlo = alenlo;
// // if arc length is larger than n/2 make it n - arc length
// if (alenmax > genomeData.length - alenmax) {
// alenmax = genomeData.length - alenmax;
// }
// // if alenlo <= n/2 start from (ihi, jlo) and go up
// // if alenhi >= n/2 start from (ilo, jhi) and go down
// if ((rjlo <= rnov2) && (alenlo <= alenmax)) {
// for (i2j = alenlo; i2j <= alenmax; i2j++) {
// // excess calculations to set range of i
// ixlo = Math.max(0, jlo - ilo - i2j);
// ixhi = Math.max(0, ihi + i2j - jhi);
// sxmx = 0;
// sxmxi = ilo + ixlo - 1;
// for (i = ilo + ixlo; i <= ihi - ixhi; i++) {
// j = i + i2j;
// absx = Math.abs(sx[j - 1] - sx[i - 1]);
// if (sxmx < absx) {
// sxmx = absx;
// sxmxi = i;
// }
// }
// rj = i2j;
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bijbss = rnjov1 * Math.pow(sxmx - 0.5, 2);
// } else {
// bijbss = rnjov1 * Math.pow(sxmx, 2);
// }
// if (bijbss > bssmax) {
// bssmax = bijbss;
// tmaxi = sxmxi;
// tmaxj = sxmxi + i2j;
// }
// }
// }
// // make arclength n - arc length
// alenmax = genomeData.length - alenmax;
// if ((rjhi >= rnov2) && (alenhi >= alenmax)) {
// for (i2j = alenhi; i2j >= alenmax; i2j--) {
// // excess calcultaions to set range of i
// ixlo = Math.max(0, jlo - ilo - i2j);
// ixhi = Math.max(0, ihi + i2j - jhi);
// sxmx = 0;
// sxmxi = ilo + ixlo - 1;
// for (i = ilo + ixlo; i <= ihi - ixhi; i++) {
// j = i + i2j;
// absx = Math.abs(sx[j - 1] - sx[i - 1]);
// if (sxmx < absx) {
// sxmx = absx;
// sxmxi = i;
// }
// }
// rj = i2j;
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bijbss = rnjov1 * Math.pow(sxmx - 0.5, 2);
// } else {
// bijbss = rnjov1 * Math.pow(sxmx, 2);
// }
// if (bijbss > bssmax) {
// bssmax = bijbss;
// tmaxi = sxmxi;
// tmaxj = sxmxi + i2j;
// }
// }
// }
// }
// }
// }
// }
// // 120
// if (isBinary) {
// if (totalSumOfSquares <= 0.0001) {
// totalSumOfSquares = 1.0;
// }
// bssmax = bssmax / (totalSumOfSquares / rn);
// } else {
// if (totalSumOfSquares <= bssmax + 0.0001) {
// totalSumOfSquares = bssmax + 1.0;
// }
// bssmax = bssmax / ((totalSumOfSquares - bssmax) / (rn - 2.0));
// }
// ostat = bssmax;
// iseg[0] = tmaxi;
// iseg[1] = tmaxj;
// }
//
// public static double ErrorSumOfSquares(int[] lengthSeg, double[] segmentSums, int k, int[]
// locations) {
// // double[] sx = segmentSums; // TODO: remove sx
// double segsx;
// int segnx, i, j;
//
// double errorSumOfSquares = 0.0;
// segsx = 0.0; // sum of signals across some segments
// segnx = 0; // sum of segment lengths across some segments
// for (i = 0; i < locations[0]; i++) {
// segsx += segmentSums[i];
// segnx += lengthSeg[i];
// }
// // (sum of signals across some segments)^2 / (sum of segment lengths across some segments)
// errorSumOfSquares += Math.pow(segsx, 2) / segnx;
// for (j = 1; j < k; j++) {
// segsx = 0.0;
// segnx = 0;
// for (i = locations[j - 1]; i < locations[j]; i++) {
// segsx += segmentSums[i];
// segnx += lengthSeg[i];
// }
// errorSumOfSquares += Math.pow(segsx, 2) / segnx;
// }
// segsx = 0.0;
// segnx = 0;
// for (i = locations[k - 1]; i < lengthSeg.length; i++) {
// segsx += segmentSums[i];
// segnx += lengthSeg[i];
// }
// errorSumOfSquares += Math.pow(segsx, 2) / segnx;
// return errorSumOfSquares;
// }
//
// public static void Combination(int r, int nmr, int[] locations, boolean rleft) {
// int i, j;
//
// i = r - 1; // i is an index
// while (locations[i] == nmr + i + 1) {
// i--;
// }
// locations[i]++;
// for (j = i + 1; j < r; j++) {
// locations[j] = locations[j - 1] + 1;
// }
// if (locations[0] == nmr + 1) {
// rleft = false;
// }
// }
//
// public static double WeightedAverage(double[] x, double[] w, int iStart, int iEnd) {
// if (iEnd == -1) {
// iEnd = x.length;
// }
// double sumWeight = 0.0;
// double weightedSum = 0.0;
// for (int i = iStart; i < iEnd; i++) {
// double wi = (w == null) ? 1.0 : w[i];
// sumWeight += wi;
// weightedSum += x[i] * wi;
// }
// return weightedSum / sumWeight;
// }
//
// public static double TailP(double b, double delta, int m, int nGrid, double tol) {
// double t, tl, dincr, bsqrtm, x, nux;
//
// dincr = (0.5 - delta) / nGrid;
// bsqrtm = b / Math.sqrt(m);
//
// tl = 0.5 - dincr;
// t = 0.5 - 0.5 * dincr;
// double tailP = 0.0;
// for (int i = 0; i < nGrid; i++) {
// tl = tl + dincr;
// t = t + dincr;
// x = bsqrtm / Math.sqrt(t * (1 - t));
// nux = Nu(x, tol);
// tailP = tailP + Math.pow(nux, 2) * IntegralInvT1tSq(tl, dincr);
// }
// tailP = 9.973557E-2 * Math.pow(b, 3) * Math.exp(-Math.pow(b, 2) / 2) * tailP;
// // since test is two-sided need to multiply tailp by 2
// tailP = 2.0 * tailP;
//
// return tailP;
// }
//
// private static double Nu(double x, double tol) {
// double nu;
// double lnu0, lnu1, dk, xk;
// int i, k;
//
// // fpnorm(x): pnorm(x, 0, 1, lower.tail=TRUE, log.p=FALSE)
// // calculates P(X <= x)
// NormalDistribution norm = new NormalDistribution(); // N(0, 1)
// if (x > 0.01) {
// lnu1 = Math.log(2.0) - 2 * Math.log(x);
// lnu0 = lnu1;
// k = 2;
// dk = 0;
// for (i = 0; i < k; i++) {
// dk = dk + 1;
// xk = -x * Math.sqrt(dk) / 2.0;
// lnu1 = lnu1 - 2.0 * norm.cumulativeProbability(xk) / dk;
// }
// while (Math.abs((lnu1 - lnu0) / lnu1) > tol) {
// lnu0 = lnu1;
// for (i = 0; i < k; i++) {
// dk = dk + 1;
// xk = -x * Math.sqrt(dk) / 2.0;
// lnu1 = lnu1 - 2.0 * norm.cumulativeProbability(xk) / dk;
// }
// k *= 2;
// }
// } else {
// lnu1 = -0.583 * x;
// }
// nu = Math.exp(lnu1);
// return nu;
// }
//
// private static double IntegralInvT1tSq(double x, double a) {
// double y;
// double integral;
//
// y = x + a - 0.5;
// integral = (8.0 * y) / (1.0 - 4.0 * Math.pow(y, 2)) + 2.0 * Math.log((1.0 + 2.0 * y) / (1.0 - 2.0
// * y));
// y = x - 0.5;
// integral = integral - (8.0 * y) / (1.0 - 4.0 * Math.pow(y, 2)) - 2.0 * Math.log((1.0 + 2.0 * y) /
// (1.0 - 2.0 * y));
// return integral;
// }
//
// public static double HTMaxP(int k, double totalSumOfSquares, double[] px, double[] sx, int al0,
// boolean isBinary) {
// int i, j, nmj;
// double rn, rj, absx, sxmx, bssmx, psmin, psmax, psdiff, bsslim, rnjov1;
// int ipsmin, ipsmax;
//
// // create blocks of size k (or k+1) to span 1 thru n
// // block partial sum max and min
// double[] bpsmax, bpsmin;
// // location of the max and min
// int[] bb;
// // variables to work on block specific data
// int nb, ilo, ihi, l;
// double psum, psdiffsq;
//
// rn = px.length;
// // number of blocks of size k (plus fraction since n/k may not be integer)
// nb = (int) (rn / k);
// // allocate memory
// bpsmax = new double[nb];
// bpsmin = new double[nb];
// bb = new int[nb];
// // block boundaries
// for (i = 0; i < nb; i++) {
// // TODO alter
// bb[i] = (int) Math.round(rn * ((double) i + 1) / nb);
// }
//
// // don't need global min and max
// // find the max, min of partial sums and their locations within blocks
// ilo = 1;
// psum = 0;
// double hTMaxP = 0.0;
// for (j = 0; j < nb; j++) // j is just an index in this for loop
// {
// sx[ilo - 1] = psum + px[ilo - 1];
// psmin = sx[ilo - 1];
// ipsmin = ilo;
// psmax = sx[ilo - 1];
// ipsmax = ilo;
// for (i = ilo; i < bb[j]; i++) {
// sx[i] = sx[i - 1] + px[i];
// if (sx[i] < psmin) {
// psmin = sx[i];
// ipsmin = i + 1;
// }
// if (sx[i] > psmax) {
// psmax = sx[i];
// ipsmax = i + 1;
// }
// }
// // store the block min, max and locations
// bpsmin[j] = psmin;
// bpsmax[j] = psmax;
// // reset ilo to be the block boundary + 1
// psum = sx[bb[j] - 1];
// ilo = bb[j] + 1;
// // calculate the bss at the block max & min pr
// i = Math.abs(ipsmin - ipsmax);
// if ((i <= k) && (i >= al0)) {
// rj = i;
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssmx = rnjov1 * Math.pow(bpsmax[j] - bpsmin[j] - 0.5, 2);
// } else {
// bssmx = rnjov1 * Math.pow(bpsmax[j] - bpsmin[j], 2);
// }
// if (hTMaxP < bssmx) {
// hTMaxP = bssmx;
// }
// }
// }
// // check the first block
// ilo = 1;
// ihi = bb[0];
// psdiff = bpsmax[0] - bpsmin[0];
// if (isBinary) {
// psdiffsq = Math.pow(psdiff - 0.5, 2);
// } else {
// psdiffsq = Math.pow(psdiff, 2);
// }
// for (j = al0; j <= k; j++) {
// rj = j;
// rnjov1 = rn / (rj * (rn - rj));
// bsslim = rnjov1 * psdiffsq;
// if (bsslim < hTMaxP) {
// break;
// } // go to 50
// sxmx = 0.0;
// for (i = ilo; i <= ihi - j; i++) {
// absx = Math.abs(sx[i + j - 1] - sx[i - 1]);
// if (sxmx < absx) {
// sxmx = absx;
// }
// }
// if (isBinary) {
// bssmx = rnjov1 * Math.pow(Math.abs(sxmx) - 0.5, 2);
// } else {
// bssmx = rnjov1 * Math.pow(sxmx, 2);
// }
// if (hTMaxP < bssmx) {
// hTMaxP = bssmx;
// }
// }
// // 50 now the minor arcs spanning the end (n)
// psdiff = Math.max(Math.abs(bpsmax[0] - bpsmin[nb - 1]), Math.abs(bpsmax[nb - 1] - bpsmin[0]));
// if (isBinary) {
// psdiffsq = Math.pow(psdiff - 0.5, 2);
// } else {
// psdiffsq = Math.pow(psdiff, 2);
// }
// for (j = al0; j <= k; j++) {
// rj = j;
// rnjov1 = rn / (rj * (rn - rj));
// bsslim = rnjov1 * psdiffsq;
// if (bsslim < hTMaxP) {
// break;
// } // go to 100
// sxmx = 0.0;
// nmj = px.length - j;
// for (i = 0; i < j; i++) // i is just an index in this for loop
// {
// absx = Math.abs(sx[i + nmj] - sx[i]);
// if (sxmx < absx) {
// sxmx = absx;
// }
// }
// if (isBinary) {
// bssmx = rnjov1 * Math.pow(Math.abs(sxmx) - 0.5, 2);
// } else {
// bssmx = rnjov1 * Math.pow(sxmx, 2);
// }
// if (hTMaxP < bssmx) {
// hTMaxP = bssmx;
// }
// }
// // 100 now the other blocks
// for (l = 1; l < nb; l++) // l is just an index in this for loop
// {
// ilo = bb[l - 1] + 1;
// ihi = bb[l];
// psdiff = bpsmax[l] - bpsmin[l];
// if (isBinary) {
// psdiffsq = Math.pow(psdiff - 0.5, 2);
// } else {
// psdiffsq = Math.pow(psdiff, 2);
// }
// for (j = al0; j <= k; j++) {
// rj = j;
// rnjov1 = rn / (rj * (rn - rj));
// bsslim = rnjov1 * psdiffsq;
// if (bsslim < hTMaxP) {
// break;
// } // go to 150
// sxmx = 0.0;
// for (i = ilo; i <= ihi - j; i++) {
// absx = Math.abs(sx[i + j - 1] - sx[i - 1]);
// if (sxmx < absx) {
// sxmx = absx;
// }
// }
// if (isBinary) {
// bssmx = rnjov1 * Math.pow(Math.abs(sxmx) - 0.5, 2);
// } else {
// bssmx = rnjov1 * Math.pow(sxmx, 2);
// }
// if (hTMaxP < bssmx) {
// hTMaxP = bssmx;
// }
// }
// // 150
// psdiff = Math.max(Math.abs(bpsmax[l] - bpsmin[l - 1]), Math.abs(bpsmax[l - 1] - bpsmin[l]));
// if (isBinary) {
// psdiffsq = Math.pow(psdiff - 0.5, 2);
// } else {
// psdiffsq = Math.pow(psdiff, 2);
// }
// for (j = al0; j <= k; j++) {
// rj = j;
// rnjov1 = rn / (rj * (rn - rj));
// bsslim = rnjov1 * psdiffsq;
// if (bsslim < hTMaxP) {
// break;
// } // go to 200
// sxmx = 0.0;
// nmj = px.length - j;
// for (i = ilo - j; i <= ilo - 1; i++) {
// absx = Math.abs(sx[i + j - 1] - sx[i - 1]);
// if (sxmx < absx) {
// sxmx = absx;
// }
// }
// if (isBinary) {
// bssmx = rnjov1 * Math.pow(Math.abs(sxmx) - 0.5, 2);
// } else {
// bssmx = rnjov1 * Math.pow(sxmx, 2);
// }
// if (hTMaxP < bssmx) {
// hTMaxP = bssmx;
// }
// }
// }// 200
// if (isBinary) {
// if (totalSumOfSquares <= 0.0001) {
// totalSumOfSquares = 1.0;
// }
// hTMaxP = hTMaxP / (totalSumOfSquares / rn);
// } else {
// if (totalSumOfSquares <= hTMaxP + 0.0001) {
// totalSumOfSquares = hTMaxP + 1.0;
// }
// hTMaxP = hTMaxP / ((totalSumOfSquares - hTMaxP) / (rn - 2.0));
// }
//
// return hTMaxP;
// }
//
// public static double TMaxP(double totalSumOfSquares, double[] px, double[] sx, int al0, boolean
// isBinary) {
// // look at the partial sums in blocks of size sqrt(n)
//
// int ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1, nb2, bi, bj, ilo, ihi, jlo, jhi,
// alenmax, i2j, alenlo, alenhi, ixlo, ixhi, nal0;
// double psum, psmin, psmax, psmin0, psmax0, bssmax, bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1,
// sij2, sijmx0, absx, sxmx, bijbss, rnov2, psdiff;
//
// // use local arrays for working within blocks
// // block partial sum max and min
// double[] bpsmax, bpsmin;
// // location of the max and min
// int[] bb, ibmin, ibmax;
//
// // t statistic corresponding to max for block i,j (and max possible)
// double[] bssbij, bssijmax;
// // row, column and order vector for reordering bssbij
// int[] bloci, blocj, loc, alen;
//
// // calculate number of blocks (nb) and block boundaries (vector bb)
// rn = px.length;
// if (px.length >= 50) {
// nb = (int) Math.round(Math.sqrt(px.length));
// } else {
// nb = 1;
// }
//
// // the number of paiwise block comparison
// nb2 = nb * (nb + 1) / 2;
// // allocate memory
// bpsmax = new double[nb];
// bpsmin = new double[nb];
// bb = new int[nb];
// ibmin = new int[nb];
// ibmax = new int[nb];
// bssbij = new double[nb2];
// bssijmax = new double[nb2];
// bloci = new int[nb2];
// blocj = new int[nb2];
// loc = new int[nb2];
// alen = new int[nb2];
//
// // block boundaries
// for (i = 0; i < nb; i++) {
// // TODO check
// bb[i] = (int) Math.round(rn * i + 1 / nb);
// }
//
// // find the max, min of partial sums and their locations within blocks
// ilo = 1;
// psum = 0;
// psmin0 = 0;
// psmax0 = 0;
// ipsmin0 = px.length;
// ipsmax0 = px.length;
// for (j = 0; j < nb; j++) // j is just an index
// {
// sx[ilo - 1] = psum + px[ilo - 1];
// psmin = sx[ilo - 1];
// ipsmin = ilo;
// psmax = sx[ilo - 1];
// ipsmax = ilo;
// for (i = ilo; i < bb[j]; i++) {
// sx[i] = sx[i - 1] + px[i];
// if (sx[i] < psmin) {
// psmin = sx[i];
// ipsmin = i + 1;
// }
// if (sx[i] > psmax) {
// psmax = sx[i];
// ipsmax = i + 1;
// }
// }
// // store the block min, max and locations
// ibmin[j] = ipsmin;
// ibmax[j] = ipsmax;
// bpsmin[j] = psmin;
// bpsmax[j] = psmax;
// // adjust global min, max and locations
// if (psmin < psmin0) {
// psmin0 = psmin;
// ipsmin0 = ipsmin;
// }
// if (psmax > psmax0) {
// psmax0 = psmax;
// ipsmax0 = ipsmax;
// }
// // reset ilo to be the block boundary + 1
// psum = sx[bb[j] - 1];
// ilo = bb[j] + 1;
// }
//
// // calculate bss for max s_i - min s_i
// psdiff = psmax0 - psmin0;
// rj = Math.abs(ipsmax0 - ipsmin0);
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssmax = rnjov1 * Math.pow(psdiff - 0.5, 2);
// } else {
// bssmax = rnjov1 * Math.pow(psdiff, 2);
// }
//
// // for a pair of blocks (i,j) calculate the max absolute t-statistic
// // at the (min_i, max_j) and (max_i, min_j) locations
// // for other indices the t-statistic can be bounded using this
//
// // if a block doesn't have the potential to exceed bssmax ignore it
// // calculate the bsslim for each block and include ones >= bssmax
//
// rnov2 = rn / 2;
// l = 0;
// nal0 = px.length - al0;
// for (i = 1; i <= nb; i++) {
// for (j = i; j <= nb; j++) {
// // calculate bsslim
// if (i == 1) {
// ilo = 1;
// } else {
// ilo = bb[i - 2] + 1;
// }
// ihi = bb[i - 1];
// if (j == 1) {
// jlo = 1;
// } else {
// jlo = bb[j - 2] + 1;
// }
// jhi = bb[j - 1];
// alenhi = jhi - ilo;
// if (alenhi > nal0) {
// alenhi = nal0;
// }
// rjhi = alenhi;
// if (i == j) {
// alenlo = 1;
// } else {
// alenlo = jlo - ihi;
// }
// if (alenlo < al0) {
// alenlo = al0;
// }
// // max S_k over block j - min S_k over block i
// sij1 = Math.abs(bpsmax[j - 1] - bpsmin[i - 1]);
// // max S_k over block i - min S_k over block j
// sij2 = Math.abs(bpsmax[i - 1] - bpsmin[j - 1]);
// // if i = j then sij1 and sij2 are the same
// sijmx0 = Math.max(sij1, sij2);
// rjlo = alenlo;
// rnjov1 = rn / Math.min(rjlo * (rn - rjlo), rjhi * (rn - rjhi));
// if (isBinary) {
// bsslim = rnjov1 * Math.pow(sijmx0 - 0.5, 2);
// } else {
// bsslim = rnjov1 * Math.pow(sijmx0, 2);
// }
// // if its as large as bssmax add block
// if (bssmax <= bsslim) {
// loc[l] = l;
// bloci[l] = i;
// blocj[l] = j;
// bssijmax[l] = bsslim;
// // max sij in the (i,j) block, t-statistic etc
// if (sij1 > sij2) {
// alen[l] = Math.abs(ibmax[j - 1] - ibmin[i - 1]);
// rj = alen[l];
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssbij[l] = rnjov1 * Math.pow(sij1 - 0.5, 2);
// } else {
// bssbij[l] = rnjov1 * Math.pow(sij1, 2);
// }
// } else {
// alen[l] = Math.abs(ibmin[j - 1] - ibmax[i - 1]);
// rj = alen[l];
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bssbij[l] = rnjov1 * Math.pow(sij2 - 0.5, 2);
// } else {
// bssbij[l] = rnjov1 * Math.pow(sij2, 2);
// }
// }
// l = l + 1;
// }
// }
// }
// nb1 = l;
// // Now sort the t-statistics by their magnitude
// int[] indices = Sort.quicksort(bssbij);
// int[] tLoc = new int[loc.length];
// for (int m = 0; m < tLoc.length; m++) {
// tLoc[m] = loc[indices[m]];
// }
// loc = tLoc;
// Arrays.sort(bssbij);
//
// // now go through the blocks in reverse order (largest down)
// for (l = nb1 - 1; l >= 0; l--) // l is an index in this for loop
// {
// k = loc[l] - 1; // k is an index in this for loop
// // need to check a block only if it has potential to increase bss
// // rjlo is the smalllest (j-i) in the block and rjhi is the largest
// bsslim = bssijmax[k];
// if (bssmax <= bsslim) {
// // bi, bj give the block location
// bi = bloci[k];
// bj = blocj[k];
// // max arc length of interest in block
// alenmax = alen[k];
// if (bi == 1) {
// ilo = 1;
// } else {
// ilo = bb[bi - 2] + 1;
// }
// ihi = bb[bi - 1];
// if (bj == 1) {
// jlo = 1;
// } else {
// jlo = bb[bj - 2] + 1;
// }
// jhi = bb[bj - 1];
// alenhi = jhi - ilo;
// if (alenhi > nal0) {
// alenhi = nal0;
// }
// rjhi = alenhi;
// if (bi == bj) {
// alenlo = 1;
// } else {
// alenlo = jlo - ihi;
// }
// if (alenlo < al0) {
// alenlo = al0;
// }
// rjlo = alenlo;
// // if arc length is larger than n/2 make is n - arc length
// if (alenmax > px.length - alenmax) {
// alenmax = px.length - alenmax;
// }
// // if alenlo <= n/2 start from (ihi, jlo) and go up
// // if alenhi >= n/2 start from (ilo, jhi) and go down
// if ((rjlo <= rnov2) && (alenlo <= alenmax)) {
// for (i2j = alenlo; i2j <= alenmax; i2j++) {
// // excess calcultaions to set range of i
// ixlo = Math.max(0, jlo - ilo - i2j);
// ixhi = Math.max(0, ihi + i2j - jhi);
// sxmx = 0;
// for (i = ilo + ixlo - 1; i < ihi - ixhi; i++) // i: 0-based index
// {
// j = i + i2j; // j: 0-based index
// absx = Math.abs(sx[j] - sx[i]);
// if (sxmx < absx) {
// sxmx = absx;
// }
// }
// rj = i2j;
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bijbss = rnjov1 * Math.pow(sxmx - 0.5, 2);
// } else {
// bijbss = rnjov1 * Math.pow(sxmx, 2);
// }
// if (bijbss > bssmax) {
// bssmax = bijbss;
// }
// }
// }
// // make arclength n - arc length
// alenmax = px.length - alenmax;
// if ((rjhi >= rnov2) && (alenhi >= alenmax)) {
// for (i2j = alenhi; i2j >= alenmax; i2j--) {
// // excess calcultaions to set range of i
// ixlo = Math.max(0, jlo - ilo - i2j);
// ixhi = Math.max(0, ihi + i2j - jhi);
// sxmx = 0;
// for (i = ilo + ixlo - 1; i < ihi - ixhi; i++) // i: 0-based index
// {
// j = i + i2j; // j: 0-based index
// absx = Math.abs(sx[j] - sx[i]);
// if (sxmx < absx) {
// sxmx = absx;
// }
// }
// rj = i2j;
// rnjov1 = rn / (rj * (rn - rj));
// if (isBinary) {
// bijbss = rnjov1 * Math.pow(sxmx - 0.5, 2);
// } else {
// bijbss = rnjov1 * Math.pow(sxmx, 2);
// }
// if (bijbss > bssmax) {
// bssmax = bijbss;
// }
// }
// }
// }
// }
// double tMaxP = 0.0;
// if (isBinary) {
// if (totalSumOfSquares <= 0.0001) {
// totalSumOfSquares = 1.0;
// }
// tMaxP = bssmax / (totalSumOfSquares / rn);
// } else {
// if (totalSumOfSquares <= bssmax + 0.0001) {
// totalSumOfSquares = bssmax + 1.0;
// }
// tMaxP = bssmax / ((totalSumOfSquares - bssmax) / (rn - 2.0));
// }
//
// return tMaxP;
// }
//
// public static double TPermP(int n1, int n2, int n, double[] genomeData, int gDOffset, double[]
// px, int nPerm, MersenneTwister rnd) {
// int np, i, m1, j, nrej;
// double xsum1, xsum2, xbar, ostat, pstat, rn1, rn2, rm1, tstat, tss, rn, cc;
//
// rn1 = n1;
// rn2 = n2;
// rn = rn1 + rn2;
// if (n1 == 1 || n2 == 1) {
// nrej = nPerm;
// // go to 110
// } else {
// xsum1 = 0.0;
// tss = 0.0;
// for (i = 0; i < n1; i++) // i: 0-based index in the for loop
// {
// px[i] = genomeData[gDOffset + i];
// xsum1 = xsum1 + genomeData[gDOffset + i];
// tss = tss + Math.pow(genomeData[gDOffset + i], 2);
// }
// xsum2 = 0.0;
// for (i = n1; i < n; i++) {
// px[i] = genomeData[gDOffset + i];
// xsum2 = xsum2 + genomeData[gDOffset + i];
// tss = tss + Math.pow(genomeData[gDOffset + i], 2);
// }
// xbar = (xsum1 + xsum2) / rn;
// tss = tss - rn * Math.pow(xbar, 2);
// if (n1 <= n2) {
// m1 = n1;
// rm1 = rn1;
// ostat = 0.99999 * Math.abs(xsum1 / rn1 - xbar);
// tstat = Math.pow(ostat, 2) * rn1 * rn / rn2;
// } else {
// m1 = n2;
// rm1 = rn2;
// ostat = 0.99999 * Math.abs(xsum2 / rn2 - xbar);
// tstat = Math.pow(ostat, 2) * rn2 * rn / rn1;
// }
// nrej = 0;
// tstat = tstat / ((tss - tstat) / (rn - 2.0));
// // if observed t is large (> 5) don't bother with permutation p-value
// // also make sure there are enough observations i.e. m1 >= 10
// if ((tstat > 25) && (m1 >= 10)) { // go to 110
// } else {
// for (np = 0; np < nPerm; np++) {
// xsum1 = 0;
// for (i = n - 1; i >= n - m1; i--) {
// cc = rnd.nextDouble();
// j = (int) (cc * (i + 1));
// j = (j > i) ? i : j;
// Swap(px[i], px[j]);
// xsum1 = xsum1 + px[i];
// }
// pstat = Math.abs(xsum1 / rm1 - xbar);
// if (ostat <= pstat) {
// nrej = nrej + 1;
// }
// }
// }
// }
// // 110
// double tPermP = (double) nrej / nPerm;
//
// return tPermP;
// }
//
// public static void Swap(double lhs, double rhs) {
// double temp;
// temp = lhs;
// lhs = rhs;
// rhs = temp;
// }
//
// public static double[] Seq(double from, double to, int length) {
// if (length == 0) {
// return null;
// }
// double[] seq = new double[length];
// double step = (to - from) / (length - 1);
// seq[0] = from;
// seq[length - 1] = to;
// for (int i = 1; i < length - 1; i++) {
// seq[i] = seq[i - 1] + step;
// }
// return seq;
// }
//
// public static void ComputeBoundary(int nPerm, double alpha, double eta, int[] sbdry) {
// System.out.println(((int) Math.floor(nPerm * alpha) + 1));
// System.out.println(sbdry.length);
// System.out.println(eta);
//
// ComputeBoundary(eta, nPerm, (int) Math.floor(nPerm * alpha) + 1, sbdry, 1E-2);
// }
//
// public static int[] GetFiniteIndices(double[] scores) {
// ArrayList<Integer> indexList = new ArrayList<Integer>();
// for (int i = 0; i < scores.length; i++) {
// if (!Numbers.isFinite(scores[i])) {
// continue;
// }
// indexList.add(i);
// }
// return Array.toIntArray(indexList);
// }
//
// // / <summary>
// // / Ported from Fortran
// // / </summary>
// // / <param name="eta"></param>
// // / <param name="nPerm"></param>
// // / <param name="maxOnes"></param>
// // / <param name="sbdry"></param>
// // / <param name="tol"></param>
// public static void ComputeBoundary(double eta, int nPerm, int maxOnes, int[] sbdry, double tol) {
// sbdry = new int[(maxOnes * (maxOnes + 1) / 2)];
// double[] etaStar = new double[maxOnes];
// double eta0 = 0;
// double etaLo = 0;
// double etaHi = 0;
// double pLo = 0;
// double pHi = 0;
// double pExcd = 0;
// int j, l;
//
// l = 0;
// sbdry[0] = nPerm - (int) (nPerm * eta);
//
// etaStar[0] = eta;
// eta0 = eta;
// for (j = 2; j <= maxOnes; j++) {
// etaHi = eta0 * 1.1;
//
// sbdry =EtaBoundary(nPerm, etaHi, j, sbdry, l + 1);
// pHi = PExceed(nPerm, j, sbdry, l + 1, pHi);
// etaLo = eta0 * 0.25;
// sbdry= EtaBoundary(nPerm, etaLo, j, sbdry, l + 1);
// pLo = PExceed(nPerm, j, sbdry, l + 1, pLo);
// while ((etaHi - etaLo) / etaLo > tol) {
// eta0 = etaLo + (etaHi - etaLo) * (eta - pLo) / (pHi - pLo);
// sbdry= EtaBoundary(nPerm, eta0, j, sbdry, l + 1);
// pExcd = PExceed(nPerm, j, sbdry, l + 1, pExcd);
// if (pExcd > eta) {
// etaHi = eta0;
// pHi = pExcd;
// } else {
// etaLo = eta0;
// pLo = pExcd;
// }
// }
// etaStar[j - 1] = eta0;
// l += j;
// }
// }
//
// private static int[] EtaBoundary(int nPerm, double eta0, int n1s, int[] sbdry, int sbdryOffset) {
// int i, k;
// // double dn;
// double tProb;
//
// double dn = nPerm - n1s; // Hypergeometric parameterizes differently
// k = 0;
// // ArrayList<Integer> sbdrytmp = new ArrayList<Integer>();
// // for (int j = 0; j < sbdry.length; j++) {
// // sbdrytmp.add(sbdry[j]);
// // System.out.println("DFS\t" + j + "\t" + sbdry[j]);
// // }
//
// for (i = 1; i <= nPerm; i++) {
// // Hypergeometric(total # of balls, # white balls, # balls drawn)
// // var hyper = new Hypergeometric(Convert.ToInt32(nPerm), Convert.ToInt32(n1s),
// // Convert.ToInt32(i)); // Hypergeometric distribution
// // P(X < k + 0.1) == P(X <= k): k = # of white balls drawn
// // tProb = hyper.CumulativeDistribution(k + 0.1);
// // Console.WriteLine("Hypergeometric({0}, {1}, {2}, {3}) = {4}", nPerm, n1s, i, k, tProb);
//
// // tProb = R.phyper(k, n1s, dn, i, true, false);
// // k = vector of quantiles representing the number of white balls drawn without replacement from
// an urn which contains both black and white balls.
// // n1s = the number of white balls in the urn.
// // dn = the number of black balls in the urn.
// // i =the number of balls drawn from the urn.
// // lower tail, not log
// // TODO????
//
// HypergeometricDistribution hypergeometricDistribution = new HypergeometricDistribution((int) (n1s
// + dn), k, (int) (n1s + dn));
// tProb = hypergeometricDistribution.upperCumulativeProbability(i);
// //tProb =1-tProb;
// System.out.println(tProb+"\t"+i+"\t"+hypergeometricDistribution.getNumberOfSuccesses()+"\t"+hypergeometricDistribution.getPopulationSize());
// // tProb = R.phyper(k, n1s, dn, i, true, false);
// // Console.WriteLine("phyper({0}, {1}, {2}, {3}) = {4}", k, n1s, dn, i, tProb);
// // var hyper = new HypergeometricDistribution(Convert.ToInt32(n1s),
// // Convert.ToInt32(dn), Convert.ToInt32(i));
// // tProb = hyper.DistributionFunction(Convert.ToInt32(k));
// if (tProb <= eta0) {
// System.out.println((sbdryOffset + k) + "\t" + i);
// sbdry[sbdryOffset + k] = i;
// k += 1;
// }
// }
// return sbdry;
// }
// //
// // // Not good! Depends on Math.NET Numerics
// // public static double phyper(double x, double NR, double NB, double n, boolean lowerTail,
// boolean logP) {
// // /* Sample of n balls from NR red and NB black ones; x are red */
// // double d, pd;
// // if (Double.isNaN(x) || Double.isNaN(NR) || Double.isNaN(NB) || Double.isNaN(n)) {
// // return x + NR + NB + n;
// // }
// // x = Math.floor(x + 1e-7);
// // NR = Math.floor(NR + 0.5);
// // NB = Math.floor(NB + 0.5);
// // n = Math.floor(n + 0.5);
// //
// // if (NR < 0 || NB < 0 || !Double.isFinite(NR + NB) || n < 0 || n > NR + NB) {
// // return Double.NaN;
// // }
// //
// // if (x * (NR + NB) > n * NR) {
// // /* Swap tails. */
// // double oldNB = NB;
// // NB = NR;
// // NR = oldNB;
// // x = n - x - 1;
// // lowerTail = !lowerTail;
// // }
// //
// // if (x < 0) {
// // return R_DT_0(lowerTail, logP);
// // }
// // if (x >= NR || x >= n) {
// // return R_DT_1(lowerTail, logP);
// // }
// //
// // d = R.dhyper(x, NR, NB, n, logP);
// // pd = R.pdhyper(x, NR, NB, n, logP);
// //
// // return logP ? R_DT_Log(d + pd, lowerTail) : R_D_Lval(d * pd, lowerTail);
// // }
// //
// // public static double dhyper(double x, double r, double b, double n, boolean giveLog) {
// // if (Double.isNaN(x) || Double.isNaN(r) || Double.isNaN(b) || Double.isNaN(n)) {
// // return x + r + b + n;
// // }
// // if (R_D_negInonint(x)) {
// // return R_D__0(giveLog);
// // }
// // x = R_D_forceint(x);
// // r = R_D_forceint(r);
// // b = R_D_forceint(b);
// // n = R_D_forceint(n);
// // if (n < x || r < x || (n - x) > b) {
// // return R.R_D__0(giveLog);
// // }
// // if (n == 0) {
// // return (x == 0) ? R.R_D__1(giveLog) : R.R_D__0(giveLog);
// // }
// // double p = n / (r + b);
// // double q = (r + b - n) / (r + b);
// //
// // double p1 = R.dbinom_raw(x, r, p, q, giveLog);
// // double p2 = R.dbinom_raw(n - x, b, p, q, giveLog);
// // double p3 = R.dbinom_raw(n, r + b, p, q, giveLog);
// //
// // return ((giveLog) ? p1 + p2 - p3 : p1 * p2 / p3);
// // }
// //
// // public static double dbinom_raw(double x, double n, double p, double q, bool giveLog) {
// // double lf, lc;
// // if (p == 0) {
// // return (x == 0) ? R_D__1(giveLog) : R_D__0(giveLog);
// // }
// // if (q == 0) {
// // return (x == n) ? R_D__1(giveLog) : R_D__0(giveLog);
// // }
// // if (x == 0) {
// // if (n == 0) {
// // return R_D__1(giveLog);
// // }
// // lc = (p < 0.1) ? (-bd0(n, n * q) - n * p) : (n * Math.log(q));
// // return R_D_exp(lc, giveLog);
// // }
// // if (x == n) {
// // lc = (q < 0.1) ? -bd0(n, n * p) - n * q : n * Math.log(p);
// // return R_D_exp(lc, giveLog);
// // }
// // if (x < 0 || x > n) {
// // return R_D__0(giveLog);
// // }
// // // n*p or n*q can underflow to zero if n and p or q are small. This
// // // used to occur in dbeta, and gives NaN as from R 2.3.0.
// // lc = stirlerr(n) - R.stirlerr(x) - R.stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
// //
// // // f = (M_2PI*x*(n-x))/n; could overflow or underflow
// // // Upto R 2.7.1:
// // // lf = log(M_2PI) + log(x) + log(n-x) - log(n);
// // // -- following is much better for x << n :
// // lf = Math.Log(R.M_2PI) + Math.Log(x) + log1p(-x / n);
// //
// // return R.R_D_exp(lc - 0.5 * lf, giveLog);
// // }
// //
// // public static double stirlerr(double n) {
// // double S0 = 0.083333333333333333333; // 1/12
// // double S1 = 0.00277777777777777777778; // 1/360
// // double S2 = 0.00079365079365079365079365; // 1/1260
// // double S3 = 0.000595238095238095238095238; // 1/1680
// // double S4 = 0.0008417508417508417508417508; // 1/1188
// //
// // // error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
// //
// // double[] sferr_halves = { 0.0, // n=0 - wrong, place holder only
// // 0.1534264097200273452913848, // 0.5
// // 0.0810614667953272582196702, // 1.0
// // 0.0548141210519176538961390, // 1.5
// // 0.0413406959554092940938221, // 2.0
// // 0.03316287351993628748511048, // 2.5
// // 0.02767792568499833914878929, // 3.0
// // 0.02374616365629749597132920, // 3.5
// // 0.02079067210376509311152277, // 4.0
// // 0.01848845053267318523077934, // 4.5
// // 0.01664469118982119216319487, // 5.0
// // 0.01513497322191737887351255, // 5.5
// // 0.01387612882307074799874573, // 6.0
// // 0.01281046524292022692424986, // 6.5
// // 0.01189670994589177009505572, // 7.0
// // 0.01110455975820691732662991, // 7.5
// // 0.010411265261972096497478567, // 8.0
// // 0.009799416126158803298389475, // 8.5
// // 0.009255462182712732917728637, // 9.0
// // 0.008768700134139385462952823, // 9.5
// // 0.008330563433362871256469318, // 10.0
// // 0.007934114564314020547248100, // 10.5
// // 0.007573675487951840794972024, // 11.0
// // 0.007244554301320383179543912, // 11.5
// // 0.006942840107209529865664152, // 12.0
// // 0.006665247032707682442354394, // 12.5
// // 0.006408994188004207068439631, // 13.0
// // 0.006171712263039457647532867, // 13.5
// // 0.005951370112758847735624416, // 14.0
// // 0.005746216513010115682023589, // 14.5
// // 0.005554733551962801371038690 // 15.0
// // };
// // double nn;
// //
// // if (n <= 15.0) {
// // nn = n + n;
// // if (nn == (int) nn)
// // return (sferr_halves[(int) nn]);
// // return (lgammafn(n + 1.0) - (n + 0.5) * Math.log(n) + n - M_LN_SQRT_2PI);
// // }
// //
// // nn = n * n;
// // if (n > 500)
// // return ((S0 - S1 / nn) / n);
// // if (n > 80)
// // return ((S0 - (S1 - S2 / nn) / nn) / n);
// // if (n > 35)
// // return ((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
// // // 15 < n <= 35 :
// // return ((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
// // }
// //
// // public static double bd0(double x, double np) {
// // if (!R_finite(x) || !R_finite(np) || np == 0.0) {
// // return Double.NaN;
// // }
// // if (Math.abs(x - np) < 0.1 * (x + np)) {
// // double v = (x - np) / (x + np);
// // double s = (x - np) * v;/* s using v -- change by MM */
// // double ej = 2 * x * v;
// // v = v * v;
// // for (int j = 1;; j++) { /* Taylor series */
// // ej *= v;
// // double s1 = s + ej / ((j << 1) + 1);
// // if (s1 == s) /* last term was effectively 0 */
// // return (s1);
// // s = s1;
// // }
// // }
// // /* else: | x - np | is not too small */
// // return (x * Math.log(x / np) + np - x);
// // }
// //
// // public static boolean R_finite(double x) {
// // return !(!Numbers.isFinite(x));
// // }
// //
// // public static double R_D_forceint(double x) {
// // return Math.floor(x + 0.5);
// // }
// //
// // public static boolean R_D_negInonint(double x) {
// // return (x < 0.0 || R_D_nonint(x));
// // }
// //
// // public static boolean R_D_nonint(double x) {
// // return (Math.abs(x - Math.floor(x + 0.5)) > 1E-7);
// // }
// //
// // public static double R_D_Lval(double p, boolean lowerTail) {
// // return (lowerTail ? (p) : (0.5 - (p) + 0.5));
// // }
// //
// // public static double R_DT_0(boolean lowerTail, boolean log) {
// // return (lowerTail ? R_D__0(log) : R_D__1(log));
// // }
// //
// // public static double R_DT_1(boolean lowerTail, boolean log) {
// // return (lowerTail ? R_D__1(log) : R_D__0(log));
// // }
// //
// // public static double R_D__0(boolean log) {
// // return log ? Double.NEGATIVE_INFINITY : 0.0;
// // }
// //
// // public static double R_D__1(boolean log) {
// // return log ? 0.0 : 1.0;
// // }
// //
// // public static double R_DT_Log(double p, boolean lowerTail) {
// // return (lowerTail ? (p) : R_Log1_Exp(p));
// // }
// //
// // /* log(1 - exp(x)) in more stable form than log1p(- R_D_qIv(x))) : */
// // public static double R_Log1_Exp(double x) {
// // return Math.log(-Math.exp(x) - 1);
// // }
// //
// // // Computes exp(x) - 1
// // public static double expm1(double x) {
// // double y, a = Math.abs(x);
// //
// // if (a < DBL_EPSILON)
// // return x;
// // if (a > 0.697)
// // return Math.exp(x) - 1; /* negligible cancellation */
// //
// // if (a > 1E-8)
// // y = Math.exp(x) - 1;
// // else
// // /* Taylor expansion, more accurate in this range */
// // y = (x / 2 + 1) * x;
// //
// // /* Newton step for solving log(1 + y) = x for y : */
// // /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
// // y -= (1 + y) * (Math.log(1 + y) - x);
// // return y;
// // }
// //
// // public static double lgammafn_sign(double x, int[] sgn) {
// // double ans, y, sinpiy;
// // double xmax = 0.0;
// // double dxrel = 0.0;
// //
// // if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
// // xmax = d1mach(2) / Math.log(d1mach(2));/* = 2.533 e305 for IEEE double */
// // dxrel = Math.sqrt(d1mach(4));/* sqrt(Eps) ~ 1.49 e-8 for IEEE double */
// // }
// //
// // if (sgn != null)
// // sgn[0] = 1;
// // if (Double.isNaN(x))
// // return x;
// //
// // if (x < 0 && (Math.floor(-x) % 2.0) == 0) {
// // if (sgn != null)
// // sgn[0] = -1;
// // }
// // if (x <= 0 && x == (int) (x)) { // Negative integer argument
// // return Double.POSITIVE_INFINITY; // +Inf, since lgamma(x) = log|gamma(x)|
// // }
// //
// // y = Math.abs(x);
// //
// // if (y < 1e-306)
// // return -Math.log(x); // denormalized range, R change
// // if (y <= 10)
// // return Math.Log(Math.Abs(gammafn(x)));
// // // ELSE y = |x| > 10 ----------------------
// //
// // if (y > xmax) {
// // return Double.POSITIVE_INFINITY;
// // }
// //
// // if (x > 0) { // i.e. y = x > 10
// // if (x > 1e17)
// // return (x * (Math.Log(x) - 1.0));
// // else if (x > 4934720.0)
// // return (M_LN_SQRT_2PI + (x - 0.5) * Math.Log(x) - x);
// // else
// // return M_LN_SQRT_2PI + (x - 0.5) * Math.Log(x) - x + lgammacor(x);
// // }
// // // else: x < -10; y = -x
// // sinpiy = Math.Abs(Math.Sin(Math.PI * y));
// // if (sinpiy == 0) { /*
// // * Negative integer argument === Now UNNECESSARY: caught above
// // */
// // MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y={0}]\n", y.ToString());
// // return ML_ERR_return_NAN();
// // }
// //
// // ans = M_LN_SQRT_PId2 + (x - 0.5) * Math.Log(y) - x - Math.Log(sinpiy) - lgammacor(y);
// //
// // if (Math.Abs((x - Math.Truncate(x - 0.5)) * ans / x) < dxrel) {
// //
// // /*
// // * The answer is less than half precision because the argument is too near a negative integer.
// // */
// // ML_ERROR(ME_PRECISION, "lgamma");
// // }
// //
// // return ans;
// // }
// //
// // private static double d1mach(int i) {
// // switch (i) {
// // case 1:
// // return Double.MIN_VALUE;
// // case 2:
// // return Double.MAX_VALUE;
// //
// // case 3: /*
// // * = FLT_RADIX ^ - DBL_MANT_DIG for IEEE: = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON
// // */
// // return 0.5 * DBL_EPSILON;
// //
// // case 4: /*
// // * = FLT_RADIX ^ (1- DBL_MANT_DIG) = for IEEE: = 2^-52 = DBL_EPSILON
// // */
// // return DBL_EPSILON;
// //
// // case 5:
// // return M_LOG10_2;
// //
// // default:
// // return 0.0;
// // }
// // }
// // \ public static double gammafn(double x)
// // {
// // double[] gamcs = new double[42]{
// // +.8571195590989331421920062399942e-2,
// // +.4415381324841006757191315771652e-2,
// // +.5685043681599363378632664588789e-1,
// // -.4219835396418560501012500186624e-2,
// // +.1326808181212460220584006796352e-2,
// // -.1893024529798880432523947023886e-3,
// // +.3606925327441245256578082217225e-4,
// // -.6056761904460864218485548290365e-5,
// // +.1055829546302283344731823509093e-5,
// // -.1811967365542384048291855891166e-6,
// // +.3117724964715322277790254593169e-7,
// // -.5354219639019687140874081024347e-8,
// // +.9193275519859588946887786825940e-9,
// // -.1577941280288339761767423273953e-9,
// // +.2707980622934954543266540433089e-10,
// // -.4646818653825730144081661058933e-11,
// // +.7973350192007419656460767175359e-12,
// // -.1368078209830916025799499172309e-12,
// // +.2347319486563800657233471771688e-13,
// // -.4027432614949066932766570534699e-14,
// // +.6910051747372100912138336975257e-15,
// // -.1185584500221992907052387126192e-15,
// // +.2034148542496373955201026051932e-16,
// // -.3490054341717405849274012949108e-17,
// // +.5987993856485305567135051066026e-18,
// // -.1027378057872228074490069778431e-18,
// // +.1762702816060529824942759660748e-19,
// // -.3024320653735306260958772112042e-20,
// // +.5188914660218397839717833550506e-21,
// // -.8902770842456576692449251601066e-22,
// // +.1527474068493342602274596891306e-22,
// // -.2620731256187362900257328332799e-23,
// // +.4496464047830538670331046570666e-24,
// // -.7714712731336877911703901525333e-25,
// // +.1323635453126044036486572714666e-25,
// // -.2270999412942928816702313813333e-26,
// // +.3896418998003991449320816639999e-27,
// // -.6685198115125953327792127999999e-28,
// // +.1146998663140024384347613866666e-28,
// // -.1967938586345134677295103999999e-29,
// // +.3376448816585338090334890666666e-30,
// // -.5793070335782135784625493333333e-31
// // };
// //
// // int i, n;
// // double y;
// // double sinpiy, value;
// //
// // int ngam = 0;
// // double xmin = 0, xmax = 0.0, xsml = 0.0, dxrel = 0.0;
// //
// // // Initialize machine dependent constants, the first time gamma() is called.
// // // FIXME for threads !
// // if (ngam == 0)
// // {
// // ngam = chebyshev_init(gamcs, 42, DBL_EPSILON / 20); //was .1*d1mach(3)
// // gammalims(ref xmin, ref xmax); //-> ./gammalims.c
// // xsml = Math.Exp(fmax2(Math.Log(DBL_MIN), -Math.Log(DBL_MAX)) + 0.01);
// // // = exp(.01)*DBL_MIN = 2.247e-308 for IEEE
// // dxrel = Math.Sqrt(DBL_EPSILON); //was sqrt(d1mach(4))
// // }
// //
// // if (Double.IsNaN(x)) return x;
// //
// // // If the argument is exactly zero or a negative integer
// // // then return NaN.
// // if (x == 0 || (x < 0 && x == (long)x))
// // {
// // ML_ERROR(ME_DOMAIN, "gammafn");
// // return Double.NaN;
// // }
// //
// // y = Math.Abs(x);
// //
// // if (y <= 10)
// // {
// // // Compute gamma(x) for -10 <= x <= 10
// // // Reduce the interval and find gamma(1 + y) for 0 <= y < 1
// // // first of all. */
// //
// // n = (int)x;
// // if (x < 0) --n;
// // y = x - n;/* n = floor(x) ==> y in [ 0, 1 ) */
// // --n;
// // value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
// // if (n == 0)
// // return value;/* x = 1.dddd = 1+y */
// //
// // if (n < 0)
// // {
// // // compute gamma(x) for -10 <= x < 1
// //
// // // exact 0 or "-n" checked already above
// //
// // // The answer is less than half precision
// // // because x too near a negative integer.
// // if (x < -0.5 && Math.Abs(x - (int)(x - 0.5) / x) < dxrel)
// // {
// // ML_ERROR(ME_PRECISION, "gammafn");
// // }
// //
// // // The argument is so close to 0 that the result would overflow.
// // if (y < xsml)
// // {
// // ML_ERROR(ME_RANGE, "gammafn");
// // if (x > 0) return Double.PositiveInfinity;
// // else return Double.NegativeInfinity;
// // }
// //
// // n = -n;
// //
// // for (i = 0; i < n; i++)
// // {
// // value /= (x + i);
// // }
// // return value;
// // }
// // else
// // {
// // // gamma(x) for 2 <= x <= 10
// //
// // for (i = 1; i <= n; i++)
// // {
// // value *= (y + i);
// // }
// // return value;
// // }
// // }
// // else
// // {
// // // gamma(x) for y = |x| > 10.
// //
// // if (x > xmax)
// // { // Overflow
// // ML_ERROR(ME_RANGE, "gammafn");
// // return Double.PositiveInfinity;
// // }
// //
// // if (x < xmin)
// // { // Underflow
// // ML_ERROR(ME_UNDERFLOW, "gammafn");
// // return 0.0;
// // }
// //
// // if (y <= 50 && y == (int)y)
// // { // compute (n - 1)!
// // value = 1.0;
// // for (i = 2; i < y; i++) value *= i;
// // }
// // else
// // { // normal case
// // value = Math.Exp((y - 0.5) * Math.Log(y) - y + M_LN_SQRT_2PI +
// // ((2 * y == (int)2 * y) ? stirlerr(y) : lgammacor(y)));
// // }
// // if (x > 0)
// // return value;
// //
// // if (Math.Abs((x - (int)(x - 0.5)) / x) < dxrel)
// // {
// //
// // /* The answer is less than half precision because */
// // /* the argument is too near a negative integer. */
// //
// // ML_ERROR(ME_PRECISION, "gammafn");
// // }
// //
// // sinpiy = Math.Sin(Math.PI * y);
// // if (sinpiy == 0)
// // { /* Negative integer arg - overflow */
// // ML_ERROR(ME_RANGE, "gammafn");
// // return Double.PositiveInfinity;
// // }
// //
// // return -Math.PI / (y * sinpiy * value);
// // }
// // }
// }
