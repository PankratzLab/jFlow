package org.genvisis.common;

import java.util.*;

public class StringVector {
	private Vector<String> v;

	public StringVector() {
		v = new Vector<String>();
	}

	public StringVector(int initialSize) {
		v = new Vector<String>(initialSize);
	}

	public StringVector(String[] initialArray) {
		v = new Vector<String>();
		for (int i = 0; i<initialArray.length; i++) {
			add(initialArray[i]);
		}
	}

	public void add(String value) {
		v.add(value);
	}

	public void addIfAbsent(String value) {
		if (!v.contains(value)) {
			v.add(value);
		}
	}

	public boolean contains(String s) {
		return v.contains(s);
	}

	public String elementAt(int i) {
		return v.elementAt(i);
	}

	public int indexOf(String elem) {
		return v.indexOf(elem);
	}

	public int indexOf(String elem, int index) {
		return v.indexOf(elem, index);
	}

	public void insertElementAt(String value, int at) {
		v.insertElementAt(value, at);
	}

	public int lastIndexOf(String elem) {
		return v.lastIndexOf(elem);
	}

	public void clear() {
		v.removeAllElements();
	}

	public void removeAllElements() {
		v.removeAllElements();
	}

	public boolean removeElement(String s) {
		return v.removeElement(s);
	}

	public void removeElementAt(int i) {
		v.removeElementAt(i);
	}

	public void setElementAt(String value, int index) {
		v.setElementAt(value, index);
	}

	public void addToElementAt(String value, int index) {
		String val = elementAt(index);
		val += value;
		v.setElementAt(val, index);
	}

	public String popFirst() {
		String str = v.firstElement();
		v.removeElementAt(0);
		return str;
	}

	public String popLast() {
		String str = v.lastElement();
		v.removeElementAt(v.size()-1);
		return str;
	}

	public String popAt(int at) {
		String str = v.elementAt(at);
		v.removeElementAt(at);
		return str;
	}

	public int size() {
		return v.size();
	}

	public Object clone() {
		StringVector newSV = new StringVector();
		for (int i = 0; i<v.size(); i++) {
			newSV.add(v.elementAt(i));
		}
		return newSV;
	}

	public String[] toArray() {
		String[] array = new String[v.size()];
		for (int i = 0; i<v.size(); i++) {
			array[i] = v.elementAt(i);
		}
		return array;
	}

	public static String[][] toMatrix(StringVector[] svs) {
    	String[][] matrix = new String[svs.length][];
    
    	for (int i = 0; i<svs.length; i++) {
    		matrix[i] = svs[i].toArray();
    	}
    
    	return matrix;
    }

	public static StringVector[] newStringVectors(int size) {
		StringVector[] svs;
		
		svs = new StringVector[size];
		for (int i = 0; i<size; i++) {
			svs[i] = new StringVector();
		}
		
		return svs;
	}
}
