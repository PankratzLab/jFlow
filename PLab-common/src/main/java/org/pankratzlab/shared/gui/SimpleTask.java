package org.pankratzlab.shared.gui;

import javax.swing.SwingWorker;

/**
 * Trivial {@link Task} type-erasure extension that uses {@link String} for both of the
 * {@link SwingWorker} type parameters.
 */
public abstract class SimpleTask extends Task<String, String> {

  /**
   * @see {@link Task#Task(String, int)}
   */
  public SimpleTask(String channel, int numSteps) {
    super(channel, numSteps);
  }
}
