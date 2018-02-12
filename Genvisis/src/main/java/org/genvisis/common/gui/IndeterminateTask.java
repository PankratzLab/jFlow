package org.genvisis.common.gui;

/**
 * Specialization of {@link Task} with an unknown length, in that it will never call
 * {@link #doStep()}.
 *
 * @see {@link Task}
 */
public abstract class IndeterminateTask<T, V> extends Task<T, V> {

  /**
   * Create a task with no steps.
   *
   * @see {@link Task#Task(String, int)}
   */
  public IndeterminateTask(String channel) {
    super(channel, 0);
  }

  @Override
  public void doStep() {
    // No-op
  }

  @Override
  public void addSteps(int numSteps) {
    // No-op
  }
}
