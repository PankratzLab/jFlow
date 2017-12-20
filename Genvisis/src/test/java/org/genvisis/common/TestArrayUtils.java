package org.genvisis.common;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

/**
 * Tests for methods in {@link ArrayUtils}
 */
public class TestArrayUtils {

	@Test
	public void testAppendToArray() {
		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {"1", "2"}, new String[] {"3", "4"}));

		assertArrayEquals(new Object[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new Object[] {"1", "2"}, new Object[] {"3", "4"}));

		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {"1", "2"}, "3", "4"));

		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {"1", "2", "3"}, new String[] {"4"}));

		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {"1", "2", "3", "4"}, new String[] {}));

		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {"1", "2", "3"}, "4"));

		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {}, "1", "2", "3", "4"));

		assertArrayEquals(new String[] {"1", "2", "3", "4"},
											ArrayUtils.appendToArray(new String[] {}, "1", "2", "3", "4"));
	}
}
