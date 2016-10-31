package org.genvisis.common;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

/**
 * String comparator that uses the double value of both strings for comparison. Will only parse a
 * given string once, and if a string is not number-compatible, will use standard string comparison.
 *
 * @author hinerm
 */
public class SciStringComparator implements Comparator<String> {
	private Map<String, Double> doubleVals = new HashMap<String, Double>();
	private boolean isNumeric = true;

	@Override
	public int compare(String o1, String o2) {
		if (!isNumeric) {
			return o1.compareTo(o2);
		}

		Double d1 = get(o1);
		Double d2 = get(o2);

		if (!isNumeric) {
			return o1.compareTo(o2);
		}

		return d1.compareTo(d2);
	}

	private Double get(String o2) {
		Double d = doubleVals.get(o2);
		if (d == null) {
			try {
				d = Double.parseDouble(o2);
				if (Double.isNaN(d)) {
					throw new IllegalArgumentException();
				}
			} catch (Exception e) {
				isNumeric = false;
				return null;
			}
			doubleVals.put(o2, d);
		}
		return d;
	}


}
