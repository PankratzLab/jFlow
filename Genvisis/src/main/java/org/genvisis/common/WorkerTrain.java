package org.genvisis.common;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Class used to add callables to a thread pool and retrieve in order
 *
 * @param <E> the object type returned by {@link WorkerTrain#next()} WARNING - make sure to call
 *          shutdown after processing, I think
 */
public class WorkerTrain<E> implements Iterator<E> {

  private final ExecutorService executor;
  private final Producer<E> producer;
  private final BlockingQueue<Future<E>> bq;
  private final int numThreads;
  private final int qBuffer;
  private final Thread trainLoader;
  private final Logger log;
  private boolean autoShutDown;

  /**
   * @param producer dishes up {@link Callable}s to run on the thread pool
   * @param nThreads thread pool size
   * @param buffer the amount of data that will be added to the processing queue, if less than
   *          double the number of threads, double the number of threads will be used as the buffer
   *          size
   * @param log
   */
  public WorkerTrain(Producer<E> producer, int nThreads, int buffer, Logger log) {
    this.numThreads = nThreads;
    this.qBuffer = Math.max(numThreads * 2, buffer);
    this.executor = Executors.newFixedThreadPool(numThreads);
    this.bq = new ArrayBlockingQueue<>(qBuffer, true);
    this.producer = producer;
    this.log = log;
    this.trainLoader = new Thread(this::loadTrain);
    this.autoShutDown = true;

    start();
  }

  /**
   * @param producer dishes up {@link Callable}s to run on the thread pool
   * @param nThreads thread pool size
   * @param log
   */
  public WorkerTrain(Producer<E> producer, int nThreads, Logger log) {
    this(producer, nThreads, 0, log);
  }

  private void start() {
    trainLoader.start();
  }

  private void loadTrain() {
    while (producer.hasNext()) {
      if (Thread.interrupted()) return;
      try {
        bq.put(executor.submit(producer.next()));
      } catch (InterruptedException e) {
        Thread.currentThread().interrupt();
      }
    }
  }

  public void shutdown() {
    trainLoader.interrupt();
    executor.shutdown();
    producer.shutdown();
  }

  public void setAutoShutDown(boolean autoShutDown) {
    this.autoShutDown = autoShutDown;
  }

  public Logger getLog() {
    return log;
  }

  @Override
  public boolean hasNext() {
    boolean hasNext = producer.hasNext() || !bq.isEmpty();
    if (!hasNext && autoShutDown) shutdown();
    return hasNext;
  }

  @Override
  public E next() {
    if (hasNext()) {
      try {
        return bq.take().get();
      } catch (InterruptedException e) {
        log.reportException(e);
        shutdown();
        Thread.currentThread().interrupt();
        return null;
      } catch (ExecutionException e) {
        log.reportException(e);
        shutdown();
        return null;
      }
    } else {
      shutdown();
      throw new NoSuchElementException();
    }
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }

  /**
   * Abstract {@link WorkerTrain.Producer} class. Provides a no-op {@link #shutdown()} and a
   * {@link #remove()} that throws {@link UnsupportedOperationException}.
   *
   * @param <E> - Generic type of the underlying {@link Callable}.
   */
  public abstract static class AbstractProducer<E> implements Producer<E> {

    @Override
    public void shutdown() {
      // No-op
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }

  /**
   * @param <E> iterator returns a {@link Callable} of this type
   */
  public interface Producer<E> extends Iterator<Callable<E>> {

    /**
     * Will be called by the {@link WorkerTrain} when there are no longer callables available. Can
     * be used to close readers etc...
     */
    public void shutdown();
  }
}
