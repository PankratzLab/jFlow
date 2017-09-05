package org.genvisis.cnv.filesys;

import org.genvisis.common.Elision;
import org.junit.Test;

import junit.framework.Assert;

/**
 * Unit tests for the {@link Compression} class
 */
public class TestCompression {

	/**
	 * @see Compression#lrrCompress
	 */
	@Test
	public void testLrrCompress() {
		// Test compressions within range
		for (float expected : new float[] {Float.NaN, 0, 1, -1,
																			 (float) (Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT
																								+ 0.005),
																			 (float) (Math.abs(Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT)
																								- 0.005)}) {
			Assert.assertEquals("Failed to compress value: " + expected + ";", 0,
													Compression.lrrCompress(expected, new byte[3], 0));
		}

		// Test compressions outside of accepted range
		for (float expected : new float[] {(float) (Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT
																								- 0.005),
																			 (float) Math.abs(Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT
																												- 0.005)}) {
			Assert.assertEquals(-1, Compression.lrrCompress(expected, new byte[3], 0));
		}
	}

	/**
	 * @see Compression#lrrDecompress
	 */
	@Test
	public void testLrrDecompress() throws Elision {
		for (float expected : new float[] {Float.NaN, 0, 1, -1, 12, -12,
																			 (float) (Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT
																								+ 0.0006f)}) {
			testDecompress(expected, Compression.lrrCompress(expected), Compression::lrrDecompress);
		}
	}

	/**
	 * Helper method to compare an expected value with its decompressed value
	 */
	private void testDecompress(float expected, byte[] compressed,
															java.util.function.Function<byte[], Float> f) throws Elision {
		Assert.assertEquals(expected, f.apply(compressed));
	}

}
