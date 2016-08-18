package org.genvisis.cnv.util;

import htsjdk.tribble.annotation.Strand;

/**
 * Static utility helper methods for working with CNVs
 */
public final class CNVHelper {

	/**
	 * @param strand
	 * @return Single-character string equivalent of the given Strand enum.
	 */
	public static String decode(Strand strand) {
		switch (strand) {
			case NEGATIVE: return "-";
			case POSITIVE: return "+";
			case NONE:
			default: return "!";
		}

	}
}
