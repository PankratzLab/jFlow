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

	/**
	 * See {@link #parseWithLocale(String, Locale)}. Uses {@link Locale#US}.
	 */
	public static int parseWithLocale(String n) {
		return parseWithLocale(n, Locale.US);
	}

	/**
	 * @param n Number string
	 * @param l {@link Locale} to use when parsing
	 * @return Integer value of the parsed number
	 */
	public static int parseWithLocale(String n, Locale l) {
		try {
			return NumberFormat.getInstance(l).parse(n).intValue();
		} catch (ParseException e) {
			throw new NumberFormatException(e.getMessage());
		}
	}
}
