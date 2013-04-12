package common;

import java.io.Serializable;
import java.util.*;

public class DoubleVector implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private Vector<String> v;

	public DoubleVector() {
		v = new Vector<String>();
	}

	public DoubleVector(int initialSize) {
		v = new Vector<String>(initialSize);
	}

	public DoubleVector(double[] initialArray) {
		v = new Vector<String>();
		for (int i = 0; i<initialArray.length; i++) {
			add(initialArray[i]);
		}
	}

	public void add(double value) {
		add(value, false);
	}

	public void add(String value) {
		add(Double.parseDouble(value), false);
	}

	public int add(double value, boolean inAlphabeticalOrder) {
		if (inAlphabeticalOrder) {
			for (int i = 0; i<=v.size(); i++) {
				if (i==v.size()||Double.parseDouble(v.elementAt(i))>value) {
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

	public void addIfAbsent(double value) {
		if (!v.contains(value+"")) {
			v.add(value+"");
		}
	}

	public boolean contains(double s) {
		return v.contains(s+"");
	}

	public double elementAt(int i) {
		return Double.parseDouble(v.elementAt(i));
	}

	public int indexOf(double elem) {
		return v.indexOf(elem+"");
	}

	public int indexOf(double elem, int index) {
		return v.indexOf(elem+"", index);
	}

	public void insertElementAt(double value, int at) {
		v.insertElementAt(value+"", at);
	}

	public int lastIndexOf(double elem) {
		return v.lastIndexOf(elem+"");
	}

	public void clear() {
		v.removeAllElements();
	}

	public void removeAllElements() {
		v.removeAllElements();
	}

	public boolean removeElement(double s) {
		return v.removeElement(s+"");
	}

	public void removeElementAt(int i) {
		v.removeElementAt(i);
	}

	public void setElementAt(double value, int index) {
		v.setElementAt(value+"", index);
	}

	public void addToElementAt(double value, int index) {
		double val = elementAt(index);
		val += value;
		v.setElementAt(val+"", index);
	}

	public double popFirst() {
		double d = Double.parseDouble(v.firstElement());
		v.removeElementAt(0);
		return d;
	}

	public double popLast() {
		double d = Double.parseDouble(v.lastElement());
		v.removeElementAt(v.size()-1);
		return d;
	}

	public double popAt(int at) {
		double d = Double.parseDouble(v.elementAt(at));
		v.removeElementAt(at);
		return d;
	}

	public int size() {
		return v.size();
	}

	public Object clone() {
		DoubleVector newDV = new DoubleVector();
		for (int i = 0; i<v.size(); i++) {
			newDV.add(Double.parseDouble(v.elementAt(i)));
		}
		return newDV;
	}

	public double[] toArray() {
		double[] array = new double[v.size()];
		for (int i = 0; i<v.size(); i++) {
			array[i] = Double.parseDouble(v.elementAt(i));
		}
		return array;
	}
	
	public static double[][] toDoubleMatrix(DoubleVector[] dvs) {
		double[][] matrix = new double[dvs.length][];
    
    	for (int i = 0; i<dvs.length; i++) {
    		matrix[i] = dvs[i].toArray();
    	}
    
    	return matrix;
    }

	public static DoubleVector[] newDoubleVectors(int size) {
		DoubleVector[] dvs;
		
		dvs = new DoubleVector[size];
		for (int i = 0; i<size; i++) {
			dvs[i] = new DoubleVector();
		}
		
		return dvs;
	}
	
}
