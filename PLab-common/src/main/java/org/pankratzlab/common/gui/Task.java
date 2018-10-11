package org.pankratzlab.common.gui;

import javax.swing.SwingWorker;

/**
 * Abstract {@link SwingWorker} extension for connecting to {@link TaskManager}. Allows incremental
 * stepping and dynamic adjustments to task scope.
 *
 * @param <T> See {@link SwingWorker}
 * @param <V> See {@link SwingWorker}
 */
public abstract class Task<T, V> extends SwingWorker<T, V> {

  private int maxSteps;
  private int currentStep = 0;
  private int lastProgress = -1;
  private final String channel;
  private final String id;

  /**
   * Create an anonymous task with no unique identifier
   *
   * @see {@link #Task(String, String, int)}
   */
  public Task(String channel, int numSteps) {
    this("ANON", channel, numSteps);
  }

  /**
   * @param id Optional identifier for this {@code Task}
   * @param channel String identifier for where to broadcast changes
   * @param numSteps Max steps in this {@code Task}
   */
  public Task(String id, String channel, int numSteps) {
    this.id = id;
    this.channel = channel;
    maxSteps = numSteps;
    TaskManager.registerTask(this);
  }

  /**
   * Signifies one step (of the max for this task) work. Use this instead of
   * {@link #setProgress(int)}
   */
  public void doStep() {
    currentStep++;
    // Compute the relative [0, 100] progress value
    int progress = (int) (((double) currentStep / maxSteps) * 100);
    if (progress != lastProgress) {
      // If the value has changed, notify
      lastProgress = progress;
      setProgress(lastProgress);
    }
  }

  /**
   * @param numSteps Number of max steps to add to this {@code Task}
   */
  public void addSteps(int numSteps) {
    maxSteps += numSteps;
  }

  /**
   * @return The unique identifier for this particular {@code Task}
   */
  public String id() {
    return id;
  }

  /**
   * @return The {@link TaskManager} channel this {@code Task} reports on.
   */
  public String channel() {
    return channel;
  }
}
