package org.genvisis.common;

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
  private final ExecutorService executor;
  private List<Future<T>> futures;
  private List<Callable<T>> bees;
  private ArrayList<T> results;
  private int reportEvery = 0;
  private final Logger log;

  /**
   * @param nThreads constructs an executor service with this many threads
   * @param timeOutDay time the job will run for before quitting
   * @param log
   */
  public WorkerHive(int nThreads, int timeOutDay, Logger log) {
    super();
    this.log = log;
    this.executor = Executors.newFixedThreadPool(nThreads);
    this.bees = new ArrayList<>();
    this.futures = new ArrayList<>();
    this.results = new ArrayList<>();
  }

  /**
   * @param bee a single job for the que
   */
  public void addCallable(Callable<T> bee) {
    bees.add(bee);
  }

  public void clear() {
    this.bees = new ArrayList<>();
    this.futures = new ArrayList<>();
    this.results = new ArrayList<>();
  }

  /**
   * @param beesToAdd multiple jobs for the que
   */
  public void addCallables(List<Callable<T>> beesToAdd) {
    bees.addAll(beesToAdd);
  }

  /**
   * @param reportEvery a message will be display everytime this many jobs are completed
   */
  public void setReportEvery(int reportEvery) {
    this.reportEvery = reportEvery;
  }

  /**
   * @param beesToAdd multiple jobs for the que
   */
  public void addCallables(Callable<T>[] beesToAdd) {
    for (Callable<T> element : beesToAdd) {
      addCallable(element);
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
          if (reportEvery > 0 && i % reportEvery == 0) {
            log.reportTimeInfo("Finished " + (i + 1) + " of " + futures.size());
          }
        } catch (InterruptedException e) {
          log.reportError("Could not complete job on internal index " + i);
          log.reportException(e);
        } catch (ExecutionException e) {
          log.reportError("Could not complete job on internal index " + i);
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
      log.reportError("No jobs were submitted to run");
    }
  }
}
