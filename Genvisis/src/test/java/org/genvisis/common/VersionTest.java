package org.genvisis.common;

import org.genvisis.common.HttpUpdate.Version;
import org.junit.Test;

import junit.framework.Assert;

/**
 * Version test for {@link HttpUpdate}.
 */
public class VersionTest {

	@Test
	public void testVersions() {
		Version v1 = new Version("25.1.0");
		Version v2 = new Version("1.1.1");
		Version v3 = new Version("1.2.0");
		Version v4 = new Version("1.2.1");

		testPair(v2, v1);
		testPair(v2, v3);
		testPair(v3, v4);
	}

	private void testPair(Version small, Version big) {
		Assert.assertTrue(small.isLessThan(big));
		Assert.assertTrue(big.isGreaterThan(small));
		Assert.assertFalse(small.isGreaterThan(big));
		Assert.assertFalse(big.isLessThan(small));
	}
}
