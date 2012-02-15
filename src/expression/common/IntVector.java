package common;

import java.util.*;

public class IntVector {
	private Vector<String> v;

	public IntVector() {
		v = new Vector<String>();
	}

	public IntVector(int initialSize) {
		v = new Vector<String>(initialSize);
	}

	public IntVector(int[] initialArray) {
		v = new Vector<String>();
		for (int i = 0; i<initialArray.length; i++) {
			add(initialArray[i]);
		}
	}

	public void add(int value) {
		add(value, false);
	}

	public int add(int value, boolean inAlphabeticalOrder) {
		if (inAlphabeticalOrder) {
			for (int i = 0; i<=v.size(); i++) {
				if (i==v.size()||Integer.parseInt(v.elementAt(i))>value) {
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

	public void addIfAbsent(int value) {
		if (!v.contains(value+"")) {
			v.add(value+"");
		}
	}

	public boolean contains(int s) {
		return v.contains(s+"");
	}

	public int elementAt(int i) {
		return Integer.parseInt(v.elementAt(i));
	}

	public int indexOf(int elem) {
		return v.indexOf(elem+"");
	}

	public int indexOf(int elem, int index) {
		return v.indexOf(elem+"", index);
	}

	public void insertElementAt(int value, int at) {
		v.insertElementAt(value+"", at);
	}

	public int lastIndexOf(int elem) {
		return v.lastIndexOf(elem+"");
	}

	public void clear() {
		v.removeAllElements();
	}

	public void removeAllElements() {
		v.removeAllElements();
	}

	public boolean removeElement(int s) {
		return v.removeElement(s+"");
	}

	public void removeElementAt(int i) {
		v.removeElementAt(i);
	}

	public void setElementAt(int value, int index) {
		v.setElementAt(value+"", index);
	}

	public void incrementAt(int index) {
		int val = elementAt(index);
		val++;
		v.setElementAt(val+"", index);
	}

	public int popFirst() {
		int i = Integer.parseInt(v.firstElement());
		v.removeElementAt(0);
		return i;
	}

	public int popLast() {
		int i = Integer.parseInt(v.lastElement());
		v.removeElementAt(v.size()-1);
		return i;
	}

	public int popAt(int at) {
		int i = Integer.parseInt(v.elementAt(at));
		v.removeElementAt(at);
		return i;
	}

	public int size() {
		return v.size();
	}

	public IntVector clone() {
		IntVector newIV = new IntVector();
		for (int i = 0; i<v.size(); i++) {
			newIV.add(Integer.parseInt(v.elementAt(i)));
		}
		return newIV;
	}

	public int[] toArray() {
		int[] array = new int[v.size()];
		for (int i = 0; i<v.size(); i++) {
			array[i] = Integer.parseInt(v.elementAt(i));
		}
		return array;
	}
	
	public int[] toArray(int[] order) {
		int[] array;
		
		array = new int[v.size()];
		if (order.length != array.length) {
			System.err.println("Error - order does not have the same number of elements (n="+order.length+") as the Vector (n="+array.length+")");
			return null;
		}
		for (int i = 0; i<array.length; i++) {
			array[i] = Integer.parseInt(v.elementAt(order[i]));
        }

		return array;
	}

	public static int[][] toIntMatrix(IntVector[] ivs) {
    	int[][] matrix = new int[ivs.length][];
    
    	for (int i = 0; i<ivs.length; i++) {
    		matrix[i] = ivs[i].toArray();
    	}
    
    	return matrix;
    }

	public static IntVector[] newIntVectors(int size) {
		IntVector[] ivs;
		
		ivs = new IntVector[size];
		for (int i = 0; i<size; i++) {
			ivs[i] = new IntVector();
		}

		return ivs;
	}
}
