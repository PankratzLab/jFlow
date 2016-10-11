// The speed of CountHash is constant for a fixed number of entries. CountVector has a comparable
// speed up to 50 unique values and then becomes much slower: 50% slower at 100 values, twice as
// slow at 200 values, 8x at 1000 values, nearly linear
package org.genvisis.common;

import java.util.Hashtable;
import java.util.Vector;

import com.google.common.primitives.Ints;

public class CountVector {
	private final Vector<String> v;
	private final IntVector iv;
	private int[] order;

	public CountVector() {
		v = new Vector<String>();
		iv = new IntVector();
		order = null;
	}

	public CountVector(Hashtable<String, String> hash) {
		this();

		String[] keys;

		keys = HashVec.getKeys(hash);
		for (String key : keys) {
			add(hash.get(key));
		}
	}

	public void add(String str) {
		int index = v.indexOf(str);

		if (index == -1) {
			v.add(str);
			iv.add(1);
		} else {
			iv.set(index, iv.get(index) + 1);
		}
		order = null;
	}

	public void clear() {
		v.clear();
		iv.clear();
		order = null;
	}

	public String[] getValues() {
		if (order == null) {
			return Array.toStringArray(v);
		} else {
			return Array.toStringArray(v, order);
		}
	}

	public int[] getCounts() {
		if (order == null) {
			return Ints.toArray(iv);
		} else {
			return Vectors.orderedArray(iv, order);
		}
	}

	public int getSize() {
		return v.size();
	}

	public String[] list() {
		String[] results = new String[v.size()];

		for (int i = 0; i < v.size(); i++) {
			results[i] = v.elementAt(i) + " (n=" + iv.elementAt(i) + ")";
		}

		return results;
	}

	public Hashtable<String, String> convertToHash() {
		Hashtable<String, String> hash;

		hash = new Hashtable<String, String>();
		for (int i = 0; i < v.size(); i++) {
			hash.put(v.elementAt(i), iv.elementAt(i) + "");
		}

		return hash;
	}

	public void sort(boolean ascending) {
		order = Sort.quicksort(Ints.toArray(iv), ascending ? Sort.ASCENDING : Sort.DESCENDING);
	}

	public static CountVector[] initArray(int size) {
		CountVector[] array = new CountVector[size];

		for (int i = 0; i < array.length; i++) {
			array[i] = new CountVector();
		}

		return array;
	}

	public int getCount(String value) {
		for (int i = 0; i < v.size(); i++) {
			if (value.equals(v.elementAt(i))) {
				return iv.elementAt(i);
			}
		}

		return 0;
	}

}
