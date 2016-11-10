package org.genvisis.cnv.util;

/**
 * For system utility methods that will not be available until J7.
 */
@Deprecated
public final class Java6Helper {

	private Java6Helper() { 
		// prevent instantiation
	}

	/**
	 * As Integer.compare in J7+
	 */
	public static int compare(int x, int y) {
    return (x < y) ? -1 : ((x == y) ? 0 : 1);
	}
}
