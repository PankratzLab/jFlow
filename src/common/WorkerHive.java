package common;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Trying to consolidate some common multi threading tasks
 */
public class WorkerHive<T> {
	private int timeOutDays;
	private ExecutorService executor;
	private List<Future<T>> futures;
	private List<Callable<T>> bees;
	private ArrayList<T> results;
	private int reportEvery = 0;
	private Logger log;

	/**
	 * @param nThreads
	 *            constructs an executor service with this many threads
	 * @param timeOutDay
	 *            time the job will run for before quitting
	 * @param log
	 */
	public WorkerHive(int nThreads, int timeOutDay, Logger log) {
		super();
		this.log = log;
		this.executor = Executors.newFixedThreadPool(nThreads);
		this.bees = new ArrayList<Callable<T>>();
		this.futures = new ArrayList<Future<T>>();
		this.results = new ArrayList<T>();
	}

	/**
	 * @param bee
	 *            a single job for the que
	 */
	public void addCallable(Callable<T> bee) {
		bees.add(bee);
	}

	/**
	 * @param beesToAdd
	 *            multiple jobs for the que
	 */
	public void addCallables(List<Callable<T>> beesToAdd) {
		bees.addAll(beesToAdd);
	}

	/**
	 * @param reportEvery
	 *            a message will be display everytime this many jobs are completed
	 */
	public void setReportEvery(int reportEvery) {
		this.reportEvery = reportEvery;
	}

	/**
	 * @param beesToAdd
	 *            multiple jobs for the que
	 */
	public void addCallables(Callable<T>[] beesToAdd) {
		for (int i = 0; i < beesToAdd.length; i++) {
			addCallable(beesToAdd[i]);
		}
	}

	/**
	 * @return an array list of objects representing the called method of the completed jobs
	 */
	public ArrayList<T> getResults() {
		return results;
	}

	// public void run(Runnable[] rs){
	// executor.invokeAll(tasks)
	// }

	/**
	 * @param awaitTermination call the shutdown method
	 */
	public void execute(boolean awaitTermination) {

		if (bees.size() > 0) {
			for (int i = 0; i < bees.size(); i++) {
				futures.add(executor.submit(bees.get(i)));
			}

			for (int i = 0; i < futures.size(); i++) {
				try {
					results.add(futures.get(i).get());
					if (reportEvery > 0 && i % reportEvery == 0 && i != 0) {
						log.reportTimeInfo("Finished " + (i) + " of " + futures.size());
					}
				} catch (InterruptedException e) {
					log.reportTimeError("Could not complete job on internal index " + i);
					log.reportException(e);
				} catch (ExecutionException e) {
					log.reportTimeError("Could not complete job on internal index " + i);
					log.reportException(e);
				}
			}
			if (awaitTermination) {
				executor.shutdown();
				try {
					executor.awaitTermination(timeOutDays, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					log.reportException(e);
				}
			}
		} else {
			log.reportTimeError("No jobs were submitted to run");
		}
	}
}
