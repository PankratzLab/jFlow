package org.genvisis.common;

/**
 * Static utility class for working with numbers.
 *
 */
public final class Numbers {

	private Numbers() {
		// prevent instantiation of static utility class
	}

	public static boolean isFinite(Double d) {
		return !d.isInfinite() && !d.isNaN();
	}

	public static boolean isFinite(Float f) {
		return !f.isInfinite() && !f.isNaN();
	}
}
