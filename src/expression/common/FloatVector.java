package common;

import java.util.*;

public class FloatVector {
	private Vector<String> v;

	public FloatVector() {
		v = new Vector<String>();
	}

	public FloatVector(int initialSize) {
		v = new Vector<String>(initialSize);
	}

	public FloatVector(float[] initialArray) {
		v = new Vector<String>();
		for (int i = 0; i<initialArray.length; i++) {
			add(initialArray[i]);
		}
	}

	public void add(float value) {
		add(value, false);
	}

	public int add(float value, boolean inAlphabeticalOrder) {
		if (inAlphabeticalOrder) {
			for (int i = 0; i<=v.size(); i++) {
				if (i==v.size()||Float.parseFloat(v.elementAt(i))>value) {
					v.insertElementAt(value+"", i);
					return i;
				}
			}
		} else {
			v.add(value+"");
			return v.size()-1;
		}

		return -1;
	}

	public void addIfAbsent(float value) {
		if (!v.contains(value+"")) {
			v.add(value+"");
		}
	}

	public boolean contains(float s) {
		return v.contains(s+"");
	}

	public float elementAt(int i) {
		return Float.parseFloat(v.elementAt(i));
	}

	public int indexOf(float elem) {
		return v.indexOf(elem+"");
	}

	public int indexOf(float elem, int index) {
		return v.indexOf(elem+"", index);
	}

	public void insertElementAt(float value, int at) {
		v.insertElementAt(value+"", at);
	}

	public int lastIndexOf(float elem) {
		return v.lastIndexOf(elem+"");
	}

	public void clear() {
		v.removeAllElements();
	}

	public void removeAllElements() {
		v.removeAllElements();
	}

	public boolean removeElement(float s) {
		return v.removeElement(s+"");
	}

	public void removeElementAt(int i) {
		v.removeElementAt(i);
	}

	public void setElementAt(float value, int index) {
		v.setElementAt(value+"", index);
	}

	public void addToElementAt(float value, int index) {
		float val = elementAt(index);
		val += value;
		v.setElementAt(val+"", index);
	}

	public float popFirst() {
		float d = Float.parseFloat(v.firstElement());
		v.removeElementAt(0);
		return d;
	}

	public float popLast() {
		float d = Float.parseFloat(v.lastElement());
		v.removeElementAt(v.size()-1);
		return d;
	}

	public float popAt(int at) {
		float d = Float.parseFloat(v.elementAt(at));
		v.removeElementAt(at);
		return d;
	}

	public int size() {
		return v.size();
	}

	public Object clone() {
		FloatVector newDV = new FloatVector();
		for (int i = 0; i<v.size(); i++) {
			newDV.add(Float.parseFloat(v.elementAt(i)));
		}
		return newDV;
	}

	public float[] toArray() {
		float[] array = new float[v.size()];
		for (int i = 0; i<v.size(); i++) {
			array[i] = Float.parseFloat(v.elementAt(i));
		}
		return array;
	}
	
	public static double[] toDoubleArray(float[] array) {
		double[] newArray = new double[array.length];
		
		for (int i = 0; i<array.length; i++) {
			newArray[i] = array[i];
        }
		
		return newArray;
	}
}
