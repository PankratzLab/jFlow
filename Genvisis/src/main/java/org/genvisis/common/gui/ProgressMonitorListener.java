package org.genvisis.common.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ProgressMonitor;
import javax.swing.Timer;

/**
 * {@link TaskListener} implementation wrapping a {@link ProgressMonitor}.
 * <p>
 * NB: {@code ProgressMonitor} provides a UI, potentially with a "Cancel" button. Consider what your
 * {@link Task} needs to clean up if canceled.
 * </p>
 */
public class ProgressMonitorListener extends AbstractTaskListener {

  // Min delay before popping a monitor
  private final static int MIN_WAIT = 500;

  // Smallest task time that will pop a monitor
  private final static int MIN_TIME = 2000;

  // Time between checks for user-cancellation
  private static final int CANCEL_DELAY = 500;

  // Timer to check if the progress monitor has been canceled by the user
  private final Timer cancelCheck = new Timer(CANCEL_DELAY, new ActionListener() {

    @Override
    public void actionPerformed(ActionEvent e) {
      if (monitor.isCanceled()) {
        doCancel();
      }
    }
  });

  // The current monitor
  private ProgressMonitor monitor = null;
  private final Object message;
  private String note;
  private final int minWait;
  private final int minTime;

  /**
   * @see #ProgressMonitorListener(Object, String, int, int, String...)
   */
  public ProgressMonitorListener(String... ids) {
    this(null, ids);
  }

  /**
   * @see #ProgressMonitorListener(Object, String, int, int, String...)
   */
  public ProgressMonitorListener(Object message, String... ids) {
    this(message, "", ids);
  }

  /**
   * @see #ProgressMonitorListener(Object, String, int, int, String...)
   */
  public ProgressMonitorListener(Object message, String note, String... ids) {
    this(message, note, MIN_WAIT, ids);
  }

  /**
   * @see #ProgressMonitorListener(Object, String, int, int, String...)
   */
  public ProgressMonitorListener(int minWait, String... ids) {
    this(minWait, MIN_TIME, ids);
  }

  /**
   * @see #ProgressMonitorListener(Object, String, int, int, String...)
   */
  public ProgressMonitorListener(int minWait, int minTime, String... ids) {
    this(null, "", minWait, minTime, ids);
  }

  /**
   * @see #ProgressMonitorListener(Object, String, int, int, String...)
   */
  public ProgressMonitorListener(Object message, String note, int minWait, String... ids) {
    this(message, note, minWait, MIN_TIME, ids);
  }

  /**
   * @param message Default: {@code null}
   * @param note Default: {@code ""}
   * @param minWait Guaranteed time, in milliseconds, before a monitor will pop up
   * @param minTime Tasks less than this time, in milliseconds, will not pop a monitor
   * @param ids See {@link AbstractTaskListener#AbstractTaskListener(String...)}
   * @see ProgressMonitor#ProgressMonitor(java.awt.Component, Object, String, int, int)
   */
  public ProgressMonitorListener(Object message, String note, int minWait, int minTime,
                                 String... ids) {
    super(ids);
    this.message = message;
    this.note = note;
    this.minWait = minWait;
    this.minTime = minTime;
  }

  @Override
  public void showProgress(int progress) {
    monitor.setProgress(progress);
  }

  @Override
  public void cancelled() {
    // Nothing special to do if canceled
  }

  @Override
  public void start() {
    // If we aren't currently monitoring anything, create a new ProgressMonitor and start the
    // appropriate timers.
    if (monitor == null) {
      // ProgressMonitor will not show with a progress value less than the minimum.
      monitor = new ProgressMonitor(null, message, note, MIN_PROGRESS, MAX_PROGRESS);
      monitor.setMillisToDecideToPopup(minWait);
      monitor.setMillisToPopup(minTime);
      monitor.setProgress(0);
      cancelCheck.restart();
    }
  }

  @Override
  public void done() {
    monitor.close();
    cancelCheck.stop();
    // ProgressMonitor instances can not be reused
    monitor = null;
  }

}
