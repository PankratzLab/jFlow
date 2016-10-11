/**
 *
 */
package org.genvisis.stats;

import org.junit.Test;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import junit.framework.Assert;

/**
 * Test for the {@link Stats} class
 *
 */
public class TestStats {
	private static final double DELTA = 1e-15;// Epsilon of error

	/**
	 * Ensure that cdf function is operating as expected <br>
	 * Can also check http://www.danielsoper.com/statcalc/calculator.aspx?id=53
	 * 
	 */
	@Test
	public void cdfTest() {
		double c = Stats.cdf(new OpdfGaussian(0, 1.0), new ObservationReal(0));
		Assert.assertEquals(.5, c, DELTA);
		double c2 = Stats.cdf(new OpdfGaussian(1.5, Math.pow(2, 2)), new ObservationReal(0));
		Assert.assertEquals(0.22662735237686826, c2, DELTA);
	}
}
