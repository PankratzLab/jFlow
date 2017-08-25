/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.Arrays;

import org.genvisis.cnv.wbs.ChangePoints.PenalizedChangePoints;
import org.genvisis.cnv.wbs.Penalties.Penalty_type;
import org.genvisis.common.Logger;
import org.junit.Test;

import junit.framework.Assert;

/**
 * 
 *
 */
public class TestPenalties {


	private static final double[] ssicPenalty10 = new double[] {13.7369557145235,
																															16.4593666708363,
																															9.73237253432228,
																															5.65382432668327,
																															10.9068240999732,
																															8.67204118914339,
																															12.0315623218365,
																															15.199128698087,
																															20.4425461247536,
																															26.2228964160272,
																															30.9001999949233};
	private static final double[] bicPenalty10 = new double[] {13.7369557145235,
																														 16.3591868292811,
																														 9.53201285121201,
																														 5.35328480201787,
																														 10.5061047337527,
																														 8.17114198136773,
																														 11.4304832725057,
																														 14.4978698072011,
																														 19.6411073923125,
																														 25.321277842031,
																														 29.8984015793719

	};

	private static final double[] mbicPenalty10 = new double[] {13.7369557145235,
																															18.5179308860493,
																															13.156975619571,
																															10.35951676357,
																															16.9252533136841,
																															15.171865966202,
																															19.4396678496892,
																															23.6986163621237,
																															30.4756869417539,
																															37.0899236154787,
																															43.094460763058

	};


	@Test
	public void TestIcCurves() {
		ChangePoints wbsChangePoints = ChangePoints.changepointsWbs(WBS.wbsIntegratedRecursiveWrapper(TestWBS.TEST_DATA,
																																												WBSUtilities.randomIntervals(TestWBS.TEST_DATA.length,
																																																										 WBS.DEFAULT_M,
																																																										 WBS.DEFAULT_SEED),
																																												new Logger()),
																																TestWBS.TEST_DATA,

																																ChangePoints.DEFAULT_TH_CONST,
																																Arrays.asList(Penalty_type.values()),
																																10);
		Assert.assertEquals(3, wbsChangePoints.getPenalizedChangePoints().size());

		for (PenalizedChangePoints penalizedChangePoints : wbsChangePoints.getPenalizedChangePoints()) {
			double[] test = null;
			switch (penalizedChangePoints.getPenaltyType()) {
				case BIC:
					test = bicPenalty10;
					break;
				case MBIC:
					test = mbicPenalty10;
					break;
				case SSIC:
					test = ssicPenalty10;
					break;
				default:
					break;
			}
			Assert.assertEquals(test.length, penalizedChangePoints.getIcCurve().size());
			for (int i = 0; i < test.length; i++) {
				// System.out.println(i + "\t" + penalizedChangePoints.getPenaltyType() + "\t"
				// + penalizedChangePoints.getIcCurve().get(i));
				Assert.assertEquals(test[i],
														penalizedChangePoints.getIcCurve().get(i),
														1e-12);
			}
			int i = 0;
			Assert.assertEquals(3, penalizedChangePoints.getPenalizedChangePoints().size());

			for (ChangePoint wbsChangePoint : penalizedChangePoints.getPenalizedChangePoints()) {
				Assert.assertEquals(TestChangePoints.CPTS[i], wbsChangePoint.getChangePoint());
				Assert.assertEquals(TestChangePoints.MINTHS[i], wbsChangePoint.getMinth(), TestWBS.DELTA);
				i++;
			}



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
}
