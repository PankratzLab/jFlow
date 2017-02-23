package org.genvisis.cnv.util;

/**
 * For system utility methods that are not available in J6
 */
@Deprecated
public final class Java6Helper {

	private Java6Helper() {
		// prevent instantiation
	}

	/**
	 * Requires J8
	 * http://grepcode.com/file/repository.grepcode.com/java/root/jdk/openjdk/8u40-b25/java/lang/Float.java#Float.isFinite%28float%29
	 */
	public static boolean isFinite(Float f) {
		return !f.isInfinite() && !f.isNaN();
	}

	/**
	 * Requires J8
	 * http://grepcode.com/file/repository.grepcode.com/java/root/jdk/openjdk/8u40-b25/java/lang/Double.java#Double.isFinite%28double%29
	 */
	public static boolean isFinite(Double d) {
		return !d.isInfinite() && !d.isNaN();
	}

	/**
	 * Requires J7
	 * http://grepcode.com/file/repository.grepcode.com/java/root/jdk/openjdk/7-b147/java/lang/Integer.java#Integer.compare%28int%2Cint%29
	 */
	public static int compare(int x, int y) {
		return (x < y) ? -1 : ((x == y) ? 0 : 1);
	}
}
