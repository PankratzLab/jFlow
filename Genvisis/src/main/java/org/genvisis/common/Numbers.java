package org.genvisis.common;

import java.text.NumberFormat;
import java.text.ParseException;
import java.util.Locale;

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

	public static int compare(int i1, int i2) {
		return new Integer(i1).compareTo(new Integer(i2));
	}

	public static int parseWithLocale(String n) {
		try {
			return NumberFormat.getInstance(Locale.US).parse(n).intValue();
		} catch (ParseException e) {
			throw new NumberFormatException(e.getMessage());
		}
	}
}
