package org.genvisis.one.JL;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerHive;

public class CallableBufferSet<E, T> {

  private final int numThreads, timeOutDay;
  private final Logger log;
  private final CallableBuffer<E, T>[] buffers;
  private boolean full;

  public CallableBufferSet(int numThreads, int timeOutDay, CallableBuffer<E, T>[] buffers,
                           Logger log) {
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
    for (CallableBuffer<E, T> buffer : buffers) {
      buffer.clearBuffer();

    }
  }

  public ArrayList<T> execute() {
    WorkerHive<T> hive = new WorkerHive<>(numThreads, timeOutDay, log);
    hive.addCallables(buffers);
    hive.execute(true);
    ArrayList<T> results = hive.getResults();
    full = false;
    return results;

  }

  public static abstract class CallableBuffer<E, T> implements Callable<T> {

    private final int bufferSize;
    private final ArrayList<E> buffer;

    public CallableBuffer(int bufferSize) {
      this.bufferSize = bufferSize;
      this.buffer = new ArrayList<>(bufferSize);
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
