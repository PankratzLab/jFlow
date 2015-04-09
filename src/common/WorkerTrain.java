package common;

import java.util.Iterator;
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
 * @param <E>
 *            the object type returned by {@link WorkerTrain#next()} WARNING - make sure to call shutdown after processing, I think
 */
public class WorkerTrain<E> implements Iterator<E> {
	private ExecutorService executor;
	private Producer<E> producer;
	private BlockingQueue<Future<E>> bq;
	private int currentQSize, numThreads, qBuffer;
	private Logger log;
	private boolean autoShutDown;

	/**
	 * @param producer
	 *            dishes up {@link Callable}s to run on the thread pool
	 * @param nThreads
	 *            thread pool size
	 * @param buffer
	 *            the amount of data that will be added to the processing que, if <= 0 or less than the number of threads, the number of threads will be used as the buffer size
	 * @param log
	 */
	public WorkerTrain(Producer<E> producer, int nThreads, int buffer, Logger log) {
		this.numThreads = nThreads;
		this.qBuffer = buffer <= 0 ? nThreads : Math.max(numThreads, buffer);// always utilize all threads given
		this.executor = Executors.newFixedThreadPool(numThreads);
		this.bq = new ArrayBlockingQueue<Future<E>>(qBuffer, true);
		this.producer = producer;
		this.log = log;
		this.currentQSize = 0;
		this.autoShutDown = true;
	}

	public void shutdown() {
		executor.shutdown();
		producer.shutdown();
	}

	public void setAutoShutDown(boolean autoShutDown) {
		this.autoShutDown = autoShutDown;
	}

	public void setProducer(Producer<E> producer) {
		this.producer = producer;
	}

	public Logger getLog() {
		return log;
	}

	@Override
	public boolean hasNext() {
		while (producer.hasNext() && currentQSize < qBuffer) {// add if possible
			bq.add(executor.submit(producer.next()));
			currentQSize++;
		}
		boolean hasNext = !bq.isEmpty();
		if (!hasNext && autoShutDown) {
			shutdown();
		}
		return hasNext;
	}

	@Override
	public E next() {
		try {
			E e = bq.take().get();
			currentQSize--;
			return e;
		} catch (ExecutionException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		shutdown();
		return null;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub

	}

	/**
	 *
	 * @param <E>
	 *            iterator returns a {@link Callable} of this type
	 */
	public interface Producer<E> extends Iterator<Callable<E>> {
		/**
		 * Will be called by the {@link WorkerTrain} when there are no longer callables available. Can be used to close readers etc...
		 */
		public void shutdown();
	}
}
