package org.pankratzlab.common.gui;

import org.pankratzlab.common.Logger;

/**
 * {@link TaskListener} wrapping a {@link Logger}. Logs progress and status notes with configuration
 * to report to the console or not.
 * <p>
 * NB: to report to a text UI component, it should be linked to the {@code Logger} itself with
 * {@link Logger#linkTextArea(javax.swing.JTextArea)}
 * </p>
 */
public class LoggerListener extends AbstractTaskListener {

  private final Logger log;
  private final boolean logText;

  /**
   * @see {@link #LoggerListener(Logger, boolean, String...)}
   */
  public LoggerListener(Logger log, String... ids) {
    this(log, false, ids);
  }

  /**
   * @param log {@link Logger} instance to report to
   * @param logText If true, attached {@code Logger} will also report to console. Default:
   *          {@code false}
   * @param channels See {@link AbstractTaskListener#AbstractTaskListener(String...)}
   */
  public LoggerListener(Logger log, boolean logText, String... channels) {
    super(channels);
    this.log = log;
    this.logText = logText;
  }

  @Override
  public void listen(Task<?, ?> task) {
    super.listen(task);
  }

  @Override
  public void showProgress(int progress) {
    report("Progress: " + progress);
  }

  @Override
  public void cancelled() {
    report("Task cancelled.");
  }

  @Override
  public void start() {
    report("Task started.");
  }

  @Override
  public void done() {
    report("Complete.", false);
  }

  /**
   * Report a message for active tasks.
   *
   * @see {@link #report(String, boolean)}
   */
  private void report(String msg) {
    report(msg, true);
  }

  /**
   * Helper method to report the {@link Task#id()}s and {@link Task#channel()}s, with a status
   * message
   *
   * @param msg Message to append
   * @param active If true, only report active tasks.
   */
  private void report(String msg, boolean active) {
    String ids = active ? activeIDs() : ids();
    String channels = active ? activeChannels() : channels();
    log.report("ID(s): " + ids + "; CH(s): " + channels + "; " + msg, true, logText);
  }
}
