package org.genvisis.one.JL;

import java.util.ArrayList;
import java.util.concurrent.Callable;


public abstract class CB<E, T> implements Callable<T> {
  private final int bufferSize;
  private final ArrayList<E> buffer;

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

  public ArrayList<E> getBuffer() {
    return buffer;
  }

  public int getBufferSize() {
    return bufferSize;
  }

  public boolean isFull() {
    return bufferSize == buffer.size();
  }

}
