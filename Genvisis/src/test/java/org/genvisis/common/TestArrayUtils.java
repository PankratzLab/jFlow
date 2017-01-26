package org.genvisis.common;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit tests for {@link ArrayUtils};
 */
public class TestArrayUtils {
	@Test
	public void stdevTest() {
		double stdev = ArrayUtils.stdev(new double[]{2, 4, 4, 4, 5, 5, 7, 9}, true);
		Assert.assertEquals(2.0, stdev, 0.0000001);
	}
}
