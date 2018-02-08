package org.genvisis.common;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import org.junit.Test;

/**
 * Tests for the {@link Bundled} class
 */
public class TestBundled {
	private static final String VALID_RSRC = "testPed.dat";
	private static final String INVALID_RSRC = "this_file_does_not_exist.txt";

	@Test
	public void testURLFound() {
		assertNotNull(Bundled.get(VALID_RSRC));
	}

	@Test
	public void testStreamFound() {
		assertNotNull(Bundled.getStream(VALID_RSRC));
	}

	@Test
	public void testRsrcFound() {
		assertNotNull(Bundled.getFile(VALID_RSRC));
	}

	@Test
	public void testRsrcNotFound() {
		assertNull(Bundled.getFile(INVALID_RSRC));
	}
}
