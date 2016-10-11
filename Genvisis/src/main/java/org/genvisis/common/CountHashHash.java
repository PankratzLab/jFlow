package org.genvisis.common;

import java.util.Hashtable;

public class CountHashHash {
	private final Hashtable<String, CountHash> hashes;

	public CountHashHash() {
		hashes = new Hashtable<String, CountHash>();
	}

	public void add(String key, String value) {
		if (!hashes.containsKey(key)) {
			hashes.put(key, new CountHash());
		}
		hashes.get(key).add(value);
	}

	public void clear() {
		hashes.clear();
	}

	public String[] getKeys() {
		return HashVec.getKeys(hashes);
	}

	public String[] getValues(String key) {
		if (hashes.containsKey(key)) {
			return hashes.get(key).getValues();
		} else {
			System.err.println("Error - This CountHash does not contain a key called '" + key + "'");
			return null;
		}
	}

	public int[] getCounts(String key) {
		if (hashes.containsKey(key)) {
			return hashes.get(key).getCounts();
		} else {
			System.err.println("Error - This CountHash does not contain a key called '" + key + "'");
			return null;
		}
	}
}
