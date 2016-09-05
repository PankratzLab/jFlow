package org.genvisis.one.jl;

import java.util.ArrayList;
import java.util.concurrent.Callable;


public abstract class CB<E, T> implements Callable<T> {
	private int bufferSize;
	private ArrayList<E> buffer;

	public CB(int bufferSize) {
		this.bufferSize = bufferSize;
		this.buffer = new ArrayList<E>(bufferSize);
	}

	public void addToBuffer(final E e) {
		buffer.add(e);
	}

	public void clearBuffer() {
		buffer.clear();
	}

	public int getBufferSize() {
		return bufferSize;
	}

	public boolean isFull() {
		return bufferSize == buffer.size();
	}

	public ArrayList<E> getBuffer() {
		return buffer;
	}

}
