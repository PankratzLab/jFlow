// Shouldn't this be a single Hashtable for faster random access?
package common;

import java.util.*;

public class CountVector {
	private Vector<String> v;
	private IntVector iv;
	private int[] order;

	public CountVector() {
		v = new Vector<String>();
		iv = new IntVector();
		order = null;
	}

	public CountVector(Hashtable<String,String> hash) {
		this();

		String[] keys;
		
		keys = HashVec.getKeys(hash);
		for (int i = 0; i < keys.length; i++) {
			add(hash.get(keys[i]));
		}
	}

	public void add(String str) {
		int index = v.indexOf(str);

		if (index==-1) {
			v.add(str);
			iv.add(1);
		} else {
			iv.incrementAt(index);
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
			return iv.toArray();
		} else {
			return iv.toArray(order);
		}
	}

	public int getSize() {
		return v.size();
	}

	public String[] list() {
		String[] results = new String[v.size()];

		for (int i = 0; i<v.size(); i++) {
			results[i] = v.elementAt(i)+" (n="+iv.elementAt(i)+")";
		}

		return results;
	}

	public Hashtable<String,String> convertToHash() {
		Hashtable<String,String> hash;
		
		hash = new Hashtable<String,String>();
		for (int i = 0; i<v.size(); i++) {
			hash.put(v.elementAt(i), iv.elementAt(i)+"");
        }
		
		return hash;
	}
	
	public void sort(boolean ascending) {
		order = Sort.quicksort(iv.toArray(), ascending?Sort.ASCENDING:Sort.DESCENDING);
	}

	public static CountVector[] initArray(int size) {
		CountVector[] array = new CountVector[size];

		for (int i = 0; i<array.length; i++) {
			array[i] = new CountVector();
		}

		return array;
	}
}
