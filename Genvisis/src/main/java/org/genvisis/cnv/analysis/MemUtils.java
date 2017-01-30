package org.genvisis.cnv.analysis;

/**
 * Helper class for computing memory stats. See http://stackoverflow.com/a/12807848/1027800
 */
public final class MemUtils {

	private MemUtils() {
		// prevent instantiation of utility class
	}

	/**
	 * @return How much memory is actually in use right now
	 */
	public static long allocatedMem() {
		Runtime r = Runtime.getRuntime();
		return r.totalMemory() - r.freeMemory();
	}

	/**
	 * @return An estimate of how much memory can be allocated before hitting {@link OutOfMemoryError}.
	 */
	public static long availableMem() {
		return Runtime.getRuntime().maxMemory() - allocatedMem();
	}
}
