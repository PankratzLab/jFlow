package org.genvisis.one.jl;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;

public class CallableBufferSet<E, T> {

	private int numThreads, timeOutDay;
	private Logger log;
	private CallableBuffer<E, T>[] buffers;
	private boolean full;

	public CallableBufferSet(int numThreads, int timeOutDay, CallableBuffer<E, T>[] buffers, Logger log) {
		super();
		this.numThreads = numThreads;
		this.timeOutDay = timeOutDay;
		this.buffers = buffers;
		this.full = false;
		this.log = log;
	}

	public void addToBuffers(final E e) throws IllegalStateException {
		boolean added = false;
		for (int i = 0; i < buffers.length; i++) {
			if (!added && !buffers[i].isFull()) {
				buffers[i].addToBuffer(e);
				added = true;
				if (i == buffers.length - 1 && buffers[i].isFull()) {
					full = true;
				}
				break;
			}
		}
		if (!added && full) {
			throw new IllegalStateException("All buffers were full, the isFull method must be checked after each addition, and then excecuted to clear");
		}
	}

	public boolean isFull() {
		return full;
	}

	public void clearBuffers() {
		for (int i = 0; i < buffers.length; i++) {
			buffers[i].clearBuffer();

		}
	}

	public ArrayList<T> execute() {
		WorkerHive<T> hive = new WorkerHive<T>(numThreads, timeOutDay, log);
		hive.addCallables(buffers);
		hive.execute(true);
		ArrayList<T> results = hive.getResults();
		full = false;
		return results;

	}
	
	

	public static abstract class CallableBuffer<E, T> implements Callable<T> {
		private int bufferSize;
		private ArrayList<E> buffer;

		public CallableBuffer(int bufferSize) {
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

}
