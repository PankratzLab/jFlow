package common;

import java.util.*;

public class ByteVector {
	private Vector<String> v;

	public ByteVector() {
		v = new Vector<String>();
	}

	public ByteVector(byte initialSize) {
		v = new Vector<String>(initialSize);
	}

	public ByteVector(byte[] initialArray) {
		v = new Vector<String>();
		for (byte i = 0; i<initialArray.length; i++) {
			add(initialArray[i]);
		}
	}

	public void add(byte value) {
		add(value, false);
	}

	public int add(byte value, boolean inAlphabeticalOrder) {
		if (inAlphabeticalOrder) {
			for (byte i = 0; i<=v.size(); i++) {
				if (i==v.size()||Byte.parseByte(v.elementAt(i))>value) {
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

	public void addIfAbsent(byte value) {
		if (!v.contains(value+"")) {
			v.add(value+"");
		}
	}

	public boolean contains(byte s) {
		return v.contains(s+"");
	}

	public byte elementAt(byte i) {
		return Byte.parseByte(v.elementAt(i));
	}

	public int indexOf(byte elem) {
		return v.indexOf(elem+"");
	}

	public int indexOf(byte elem, byte index) {
		return v.indexOf(elem+"", index);
	}

	public void insertElementAt(byte value, byte at) {
		v.insertElementAt(value+"", at);
	}

	public int lastIndexOf(byte elem) {
		return v.lastIndexOf(elem+"");
	}

	public void clear() {
		v.removeAllElements();
	}

	public void removeAllElements() {
		v.removeAllElements();
	}

	public boolean removeElement(byte s) {
		return v.removeElement(s+"");
	}

	public void removeElementAt(byte i) {
		v.removeElementAt(i);
	}

	public void setElementAt(byte value, byte index) {
		v.setElementAt(value+"", index);
	}

	public void incrementAt(byte index) {
		byte val = elementAt(index);
		val++;
		v.setElementAt(val+"", index);
	}

	public byte popFirst() {
		byte i = Byte.parseByte(v.firstElement());
		v.removeElementAt(0);
		return i;
	}

	public byte popLast() {
		byte i = Byte.parseByte(v.lastElement());
		v.removeElementAt(v.size()-1);
		return i;
	}

	public byte popAt(byte at) {
		byte i = Byte.parseByte(v.elementAt(at));
		v.removeElementAt(at);
		return i;
	}

	public int size() {
		return v.size();
	}

	public ByteVector clone() {
		ByteVector newIV = new ByteVector();
		for (byte i = 0; i<v.size(); i++) {
			newIV.add(Byte.parseByte(v.elementAt(i)));
		}
		return newIV;
	}

	public byte[] toArray() {
		byte[] array = new byte[v.size()];
		for (byte i = 0; i<v.size(); i++) {
			array[i] = Byte.parseByte(v.elementAt(i));
		}
		return array;
	}

}
