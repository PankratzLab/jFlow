package org.pankratzlab.common;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Class used to add callables to a thread pool and retrieve in order
 *
 * @param <E> the object type returned by {@link WorkerTrain#next()}
 */
public class WorkerTrain<E> implements Iterator<E>, AutoCloseable {

  private final ExecutorService executor;
  private final Producer<E> producer;
  private final BlockingQueue<Future<E>> bq;
  private final int numThreads;
  private final int qBuffer;
  private final Thread trainLoader;
  private final Logger log;
  private final AtomicBoolean producerHasNext;
  private final ReentrantLock lock;

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
    producerHasNext = new AtomicBoolean(producer.hasNext());
    lock = new ReentrantLock(true);

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
      // We know there's something in the producer so hasNext will be true until the next item is
      // enqueued
      producerHasNext.set(true);
      Future<E> next = executor.submit(producer.next());
      boolean done = false;
      while (!done) {
        if (Thread.interrupted()) return;
        // The queue may be full; we do not want to maintain a lock on the queue, which would
        // deadlock with next(), so we wait until there's room and our entry is accepted
        try {
          lock.lock();
          // Updating the queue and updating producerHasNext must be atomic to avoid race conditions
          if (bq.offer(next)) {
            done = true;
            producerHasNext.set(producer.hasNext());
          }
        } finally {
          lock.unlock();
        }
      }
    }
  }

  @Override
  public void close() {
    trainLoader.interrupt();
    executor.shutdown();
    producer.shutdown();
  }

  public Logger getLog() {
    return log;
  }

  @Override
  public boolean hasNext() {
    // If the producer has unqueued tasks or the queue has remaining values, then next() will
    // eventually return something. But we cannot check the producer itself directly, as the act of
    // retrieving from the producer and enqueuing cannot be atomic without introducing potential
    // deadlock. So we use a boolean surrogate whose update is synchronized with enqueuing and 
    // dequeuing actions.
    boolean hasNext;
    try {
      lock.lock();
      hasNext = producerHasNext.get() || !bq.isEmpty();
    } finally {
      lock.unlock();
    }

    if (!hasNext) close();
    return hasNext;
  }

  @Override
  public E next() {
    if (hasNext()) {
      try {

        Future<E> next = null;

        while (next == null) {
          // We know the queue has something to return, but it may not have received it yet if the
          // producer is slow. So we wait for it to arrive.
          next = bq.poll();
        }

        return next.get();
      } catch (InterruptedException e) {
        log.reportException(e);
        close();
        Thread.currentThread().interrupt();
        return null;
      } catch (ExecutionException e) {
        log.reportException(e);
        close();
        return null;
      }
    } else {
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
