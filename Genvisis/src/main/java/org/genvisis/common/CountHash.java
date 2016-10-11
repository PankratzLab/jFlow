// The speed of CountHash is constant for a fixed number of entries. CountVector has a comparable
// speed up to 50 unique values and then becomes much slower: 50% slower at 100 values, twice as
// slow at 200 values, 8x at 1000 values, nearly linear
package org.genvisis.common;

import java.util.Enumeration;
import java.util.Hashtable;

public class CountHash {
	private final Hashtable<String, String> hash;
	private final Hashtable<String, String> orderOfEntry;
	private String[] values;
	private int totalCount;

	public CountHash() {
		hash = new Hashtable<String, String>();
		orderOfEntry = new Hashtable<String, String>();
		totalCount = 0;
	}

	public void add(String value) {
		int count;

		if (hash.containsKey(value)) {
			count = Integer.parseInt(hash.get(value));
		} else {
			count = 0;
			orderOfEntry.put(hash.size() + "", value);
		}
		count++;
		hash.put(value, count + "");
		totalCount++;
	}

	public Hashtable<String, String> getHash() {
		return hash;
	}

	public void remove(String value, boolean reportIfAbsent) {
		if (hash.containsKey(value)) {
			totalCount -= Integer.parseInt(hash.get(value));
			hash.remove(value);
		} else if (reportIfAbsent) {
			System.err.println("Error - '"	+ value
													+ "' was not seen in this instance of CountHash and therefore could not be removed");
		}
	}

	public void clear() {
		hash.clear();
		totalCount = 0;
	}

	public String[] getValues() {
		if (values == null) {
			values = new String[hash.size()];
			for (int i = 0; i < values.length; i++) {
				values[i] = orderOfEntry.get(i + "");
			}
		}
		return values;
	}

	public void sortValuesAlphanumerically() {
		values = HashVec.getKeys(hash);
		values = Sort.putInOrder(values, Array.isAllNumbers(values));
	}

	public int getSize() {
		return hash.size();
	}

	public int getSizeOfCountGreaterThan(int minCount) {
		int size;
		Enumeration<String> keys;

		size = 0;
		keys = hash.keys();
		while (keys.hasMoreElements()) {
			if (Integer.parseInt(hash.get(keys.nextElement())) >= minCount) {
				size++;
			}
		}
		return size;
	}

	public int getSizeOfCountEquals(int count) {
		int size;
		Enumeration<String> keys;

		size = 0;
		keys = hash.keys();
		while (keys.hasMoreElements()) {
			if (Integer.parseInt(hash.get(keys.nextElement())) == count) {
				size++;
			}
		}
		return size;
	}

	public int getTotalCount() {
		return totalCount;
	}

	public int[] getCounts() {
		int[] counts;

		getValues();
		counts = new int[values.length];
		for (int i = 0; i < counts.length; i++) {
			counts[i] = Integer.parseInt(hash.get(values[i]));
		}
		return counts;
	}

	public int getCount(String value) {
		if (hash.containsKey(value)) {
			return Integer.parseInt(hash.get(value));
		} else {
			return 0;
		}
	}
}
