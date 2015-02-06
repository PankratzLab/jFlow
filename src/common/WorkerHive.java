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
 * 
 * @param <T>
 *            parameterize and run callables returning this type
 */
public class WorkerHive<T> {
	private int timeOutDays;
	private ExecutorService executor;
	private List<Future<T>> futures;
	private List<Callable<T>> bees;
	private ArrayList<T> results;
	private int reportEvery = 0;
	private Logger log;

	public WorkerHive(int nThreads, int timeOutDay, Logger log) {
		super();
		this.log = log;
		this.executor = Executors.newFixedThreadPool(nThreads);
		this.bees = new ArrayList<Callable<T>>();
		this.futures = new ArrayList<Future<T>>();
		this.results = new ArrayList<T>();
	}

	public void addCallable(Callable<T> bee) {
		bees.add(bee);
	}

	public void addCallables(List<Callable<T>> beesToAdd) {
		bees.addAll(beesToAdd);
	}

	public void setReportEvery(int reportEvery) {
		this.reportEvery = reportEvery;
	}

	public void addCallables(Callable<T>[] beesToAdd) {
		for (int i = 0; i < beesToAdd.length; i++) {
			addCallable(beesToAdd[i]);
		}
	}

	public ArrayList<T> getResults() {
		return results;
	}

	// public void run(Runnable[] rs){
	// executor.invokeAll(tasks)
	// }

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
