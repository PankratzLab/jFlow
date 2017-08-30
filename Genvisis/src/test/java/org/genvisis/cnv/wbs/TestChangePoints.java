/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.List;

import org.genvisis.common.Logger;
import org.junit.Test;

import junit.framework.Assert;

/**
 * 
 */
public class TestChangePoints {
	/**
	 * 
	 */
	static final double[] MINTHS = new double[] {6.817240610450304, 5.237307724561478,
																							 4.366603647474328};
	/**
	 * 
	 */
	static final int[] CPTS = new int[] {150, 128, 168};
	private static final double[] expectedMeans = new double[] {0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															0.0164107194566035,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															-1.15620205572082,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															1.07842263907085,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591,
																															-0.013626005675591};

	@Test
	public void testComputeSigma() {
		double sigma = ChangePoints.computeSigma(TestWBS.TEST_DATA);
		// 0.9803384 in R with mad(diff(x)/sqrt(2))
		Assert.assertEquals(0.9803383614231946, sigma, TestWBS.DELTA);
	}

	@Test
	public void testComputTh() {
		int[][] randomIntervals = WBSUtilities.randomIntervals(TestWBS.TEST_DATA.length, WBS.DEFAULT_M,
																													 WBS.DEFAULT_SEED);
		List<ChangePoint> wList = WBS.wbsIntegratedRecursiveWrapper(TestWBS.TEST_DATA,
																											randomIntervals,
																											new Logger());
		double sigma = ChangePoints.computeSigma(TestWBS.TEST_DATA);
		double th = ChangePoints.computTh(ChangePoints.DEFAULT_TH_CONST, sigma,
																			TestWBS.TEST_DATA.length, -1, wList);
		// 4.304432 in R with
		// https://github.com/cran/wbs/blob/8919e7c35389b92e27b6948572271e0843b5f808/R/changepoints.R#L96
		Assert.assertEquals(4.304431734978156, th, TestWBS.DELTA);

		double tho = ChangePoints.computTh(ChangePoints.DEFAULT_TH_CONST, sigma,
																			 TestWBS.TEST_DATA.length, 1000, wList);

		Assert.assertEquals(0, tho, TestWBS.DELTA);

		double thk = ChangePoints.computTh(ChangePoints.DEFAULT_TH_CONST, sigma,
																			 TestWBS.TEST_DATA.length, 50, wList);
		Assert.assertEquals(1.6022492545358957, thk, TestWBS.DELTA);
	}

	@Test
	public void testChangepointsSbs() {
		int[][] randomIntervals = WBSUtilities.randomIntervals(TestWBS.TEST_DATA.length, WBS.DEFAULT_M,
																													 WBS.DEFAULT_SEED);
		List<ChangePoint> wList = WBS.wbsIntegratedRecursiveWrapper(TestWBS.TEST_DATA,
																											randomIntervals,
																											new Logger());
		// res.tmp From R
		// [,1] [,2]
		// [1,] 150 6.817241
		// [2,] 128 5.237308
		// [3,] 168 4.366604

		ChangePoints sbsChangePoints = ChangePoints.changepointsSbs(wList, TestWBS.TEST_DATA,

																																ChangePoints.DEFAULT_TH_CONST, -1);
		Assert.assertEquals(3, sbsChangePoints.getStandardTHChangePoints().size());

		int i = 0;
		for (ChangePoint wbsChangePoint : sbsChangePoints.getStandardTHChangePoints()) {
			Assert.assertEquals(CPTS[i], wbsChangePoint.getChangePoint());
			Assert.assertEquals(MINTHS[i], wbsChangePoint.getMinth(), TestWBS.DELTA);
			i++;
		}
	}

	@Test
	public void testComputeMeansBetweenChangepoints() {
		ChangePoints sbsChangePoints = ChangePoints.changepointsSbs(WBS.wbsIntegratedRecursiveWrapper(TestWBS.TEST_DATA,
																																												WBSUtilities.randomIntervals(TestWBS.TEST_DATA.length,
																																																										 WBS.DEFAULT_M,
																																																										 WBS.DEFAULT_SEED),
																																												new Logger()),
																																TestWBS.TEST_DATA,

																																ChangePoints.DEFAULT_TH_CONST, -1);

		double[] actual = ChangePoints.computeMeansBetweenChangepoints(TestWBS.TEST_DATA,
																																	 sbsChangePoints.getStandardTHChangePoints());
		Assert.assertEquals(expectedMeans.length, actual.length);

		for (int i = 0; i < actual.length; i++) {
			Assert.assertEquals(expectedMeans[i], actual[i], 1e-13);

		}


	}

	@Test
	public void TestComputeMinLogLik() {
		ChangePoints sbsChangePoints = ChangePoints.changepointsSbs(WBS.wbsIntegratedRecursiveWrapper(TestWBS.TEST_DATA,
																																												WBSUtilities.randomIntervals(TestWBS.TEST_DATA.length,
																																																										 WBS.DEFAULT_M,
																																																										 WBS.DEFAULT_SEED),
																																												new Logger()),
																																TestWBS.TEST_DATA,

																																ChangePoints.DEFAULT_TH_CONST, -1);

		Assert.assertEquals(-11.758062621950765,
												ChangePoints.computeMinLogLik(TestWBS.TEST_DATA,
																											sbsChangePoints.getStandardTHChangePoints()),
												TestWBS.DELTA);

	}


	// 5.803962 ssic
	// 5.703782 bic
}