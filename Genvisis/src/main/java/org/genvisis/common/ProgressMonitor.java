package org.genvisis.common;

import java.util.EmptyStackException;
import java.util.HashMap;
import java.util.Stack;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

import org.genvisis.common.ProgressMonitor.DISPLAY_MODE;

class Task {

  public Task(String name) {
    this.name = name;
    creationTime = System.nanoTime();
    lastUpdate = System.nanoTime();
    lastDisplay = System.nanoTime();
  }

  private final long creationTime;

  public long getCreationTime() {
    return creationTime;
  }

  private final String name;

  public String getName() {
    return name;
  }

  private String label;

  public void setLabel(String lbl) {
    label = lbl;
  }

  public String getLabel() {
    return label;
  }

  private boolean indeterm;

  public void setIndeterminate(boolean indet) {
    indeterm = indet;
  }

  public boolean getIndeterminate() {
    return indeterm;
  }

  private int numExp; // TODO should be final?

  public void setExpectedUpdateCount(int cnt) {
    numExp = cnt;
  }

  public int getExpectedUpdateCount() {
    return numExp;
  }

  private volatile int updateCount;

  public int getUpdateCount() {
    return updateCount;
  }

  private int timeoutMins = ProgressMonitor.DEFAULT_TIMEOUT_MINS;

  public void setTimeout(int timeout) {
    timeoutMins = timeout;
  }

  public int getTimeout() {
    return timeoutMins;
  }

  private volatile long lastUpdate;

  protected void updateTime() {
    lastUpdate = System.nanoTime();
  }

  public long getLastUpdate() {
    return lastUpdate;
  }

  private volatile long lastDisplay;

  protected synchronized void updateDisplayTime() {
    lastDisplay = System.nanoTime();
  }

  public long getLastDisplayUpdateTime() {
    return lastDisplay;
  }

  private volatile boolean hasWarned = false;

  public void warned() {
    hasWarned = true;
  }

  public boolean alreadyWarned() {
    return hasWarned;
  }

  private volatile DISPLAY_MODE displayMode;

  public void setDisplayMode(DISPLAY_MODE mode) {
    displayMode = mode;
  }

  public DISPLAY_MODE getDisplayMode() {
    return displayMode;
  }

  // private final Object CANCEL_LOCK = new Object();
  // private volatile boolean cancelRequested = false;
  // public void requestCancel() { synchronized(CANCEL_LOCK) { cancelRequested = true; }}
  // public boolean cancelRequested() { synchronized(CANCEL_LOCK) { return cancelRequested; }}

  protected void updateOnce() {
    update(1);
  }

  protected void update(int numUpdates) {
    updateCount += numUpdates;
    if (updateCount > numExp) {
      numExp = updateCount;
    }
    updateTime();
  }

}


public class ProgressMonitor {

  public static enum DISPLAY_MODE {
                                   GUI_AND_CONSOLE, GUI_ONLY, CONSOLE_ONY;
  }

  public static final int DEFAULT_TIMEOUT_MINS = 10;
  private static final int INDET_ELAPSED_LOG_SECONDS = 180; // 3 minutes between indeterminate task
                                                            // updates
  private static final int DET_ELAPSED_LOG_SECONDS = 10; // seconds between determinate task updates

  HashMap<String, Task> taskMap = new HashMap<String, Task>();
  // subclass for no-duplicates behavior ("HashStack")
  Stack<String> taskUpdateStack = new Stack<String>() {
    ;
    private static final long serialVersionUID = 1L;

    @Override
    public String push(String item) {
      this.remove(item);
      return super.push(item);
    };
  };

  private final JProgressBar internalProgBar;
  private final Logger internalLogger;

  public ProgressMonitor(JProgressBar progBar, Logger log) {
    internalProgBar = progBar;
    internalLogger = log;
    Timer timer = new Timer("TaskMonitorTimer", true);
    timer.schedule(monitorTask, 60 * 1000, 60 * 1000); // wait a minute to begin monitoring, and
                                                       // wait a minute between each check
  }

  /**
   * Start monitoring a new indeterminate task with a default timeout
   *
   * @param taskName Unique String ID for task
   * @param label Initial display label for task (can be changed later)
   */
  public synchronized void beginIndeterminateTask(String taskName, String label,
                                                  DISPLAY_MODE mode) {
    beginTask(taskName, label, true, 0, DEFAULT_TIMEOUT_MINS, mode);
  }

  /**
   * Begin monitoring a new task, with a default timeout
   *
   * @param taskName Unique String ID for task
   * @param label Initial display label for task (can be changed later)
   * @param indeterminate Is the task indeterminate or not?
   * @param expectedUpdateCount Number of expected update calls (used to display progress and
   *        calculate percentage complete)
   */
  public synchronized void beginDeterminateTask(String taskName, String label,
                                                int expectedUpdateCount, DISPLAY_MODE mode) {
    beginTask(taskName, label, false, expectedUpdateCount, DEFAULT_TIMEOUT_MINS, mode);
  }

  /**
   * Begin tracking a new task
   *
   * @param taskName Unique String ID for task
   * @param label Initial display label for task (can be changed later)
   * @param indeterminate Is the task indeterminate or not?
   * @param expectedUpdateCount Number of expected update calls (used to display progress and
   *        calculate percentage complete)
   * @param timeout Time (in minutes) before a notice is displayed to the user
   */
  private synchronized void beginTask(String taskName, String label, boolean indeterminate,
                                      int expectedUpdateCount, int timeout, DISPLAY_MODE mode) {
    if (taskMap.containsKey(taskName)) {
      // likely programmer error, or trying to run a task twice at the same time
      // TODO, sometimes causes conflict with opening multiple ScatterPlots at the same time...
      throw new IllegalArgumentException("Error - task with key [" + taskName
                                         + "] already exists!");
    }

    Task newTask = new Task(taskName);
    newTask.setLabel(label);
    newTask.setIndeterminate(indeterminate);
    newTask.setExpectedUpdateCount(expectedUpdateCount);
    newTask.setTimeout(timeout);
    newTask.setDisplayMode(mode);
    taskMap.put(taskName, newTask);

    updateDisplay(newTask);
  }

  public synchronized void endTask(String taskName) {
    Task myTask = taskMap.get(taskName);
    if (myTask == null) {
      throw new IllegalArgumentException("No available task named [" + taskName + "]");
    }
    taskMap.remove(taskName);
    taskUpdateStack.remove(taskName);

    String nextTask = null;
    try {
      nextTask = taskUpdateStack.peek();
    } catch (EmptyStackException e) {
    }

    String msg = ext.getTime() + "]\tTask '" + myTask.getName() + "' with status '"
                 + myTask.getLabel() + "' is complete";
    if (internalLogger != null) {
      internalLogger.report(msg);
    } else {
      System.out.println(msg);
    }

    if (nextTask == null) {
      clearDisplay();
    } else {
      Task next = taskMap.get(nextTask);
      if (next == null) {
        throw new IllegalArgumentException("No available task named [" + nextTask + "]");
      }
      updateDisplay(next);
    }
  }

  public void updateTask(String taskName) {
    updateTask(taskName, 1);
  }

  public synchronized void updateTask(String taskName, int numUpdates) {
    Task myTask = taskMap.get(taskName);
    if (myTask == null) {
      throw new IllegalArgumentException("No available task named [" + taskName + "]");
    }

    if (!myTask.getIndeterminate()) {
      myTask.update(numUpdates);
    } else {
      myTask.updateTime();
    }

    updateDisplay(myTask);
  }

  public synchronized String getCurrentTaskName() {
    return taskUpdateStack.peek();
  }

  // public synchronized void cancelTask(String taskName) {
  // Task myTask = taskMap.get(taskName);
  // if (myTask == null) {
  // throw new IllegalArgumentException("No available task named [" + taskName + "]");
  // }
  // myTask.requestCancel();
  // }

  public synchronized void changeTaskLabelNoUpdate(String taskName, String newLabel) {
    Task myTask = taskMap.get(taskName);
    if (myTask == null) {
      throw new IllegalArgumentException("No available task named [" + taskName + "]");
    }

    myTask.setLabel(newLabel);

    updateDisplay(myTask);
  }

  public synchronized void changeTaskLabelWithUpdate(String taskName, String newLabel) {
    Task myTask = taskMap.get(taskName);
    if (myTask == null) {
      throw new IllegalArgumentException("No available task named [" + taskName + "]");
    }

    myTask.setLabel(newLabel);

    updateTask(taskName);
  }

  private void clearDisplay() {
    if (internalProgBar != null) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          internalProgBar.setStringPainted(false);
          internalProgBar.setString(null);
          internalProgBar.setMinimum(0);
          internalProgBar.setMaximum(0);
          internalProgBar.setValue(0);
          internalProgBar.setIndeterminate(false);
          internalProgBar.setVisible(false);
        }
      });
    } else {
      // TODO console logging?
    }
  }

  private synchronized void updateDisplay(final Task task) {
    taskUpdateStack.push(task.getName());
    if (internalProgBar != null && (task.getDisplayMode() == DISPLAY_MODE.GUI_AND_CONSOLE
                                    || task.getDisplayMode() == DISPLAY_MODE.GUI_ONLY)) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          internalProgBar.setVisible(true);
          internalProgBar.setStringPainted(task.getLabel() != null);

          double rawPct = 100d * (task.getUpdateCount()) / (task.getExpectedUpdateCount());

          String pct = (task.getIndeterminate() ? "" : " (" + ext.formDeci(rawPct, 0) + "%)");

          internalProgBar.setString(task.getLabel() + pct);
          internalProgBar.setMinimum(0);
          internalProgBar.setMaximum(task.getExpectedUpdateCount());
          internalProgBar.setValue(task.getUpdateCount());
          internalProgBar.setIndeterminate(task.getIndeterminate());
          internalProgBar.repaint();

          task.updateDisplayTime();
        }
      });
    }
    if (task.getDisplayMode() == DISPLAY_MODE.GUI_AND_CONSOLE
        || task.getDisplayMode() == DISPLAY_MODE.CONSOLE_ONY) {
      long currTime = System.nanoTime();
      if (task.getIndeterminate()) {
        long elapsedLong = currTime - task.getLastDisplayUpdateTime();// task.getCreationTime();
        elapsedLong /= 60;
        elapsedLong /= 1000;
        elapsedLong /= 1000;
        int elapsed = (int) elapsedLong;
        if (elapsed > INDET_ELAPSED_LOG_SECONDS) {
          String msg = ext.getTime() + "]\tTask '" + task.getName() + "' with status '"
                       + task.getLabel() + "' has been updated";
          if (internalLogger != null) {
            internalLogger.report(msg);
          } else {
            System.out.println(msg);
          }
        }
      } else {
        long elapsedLong = currTime - task.getLastDisplayUpdateTime();// task.getCreationTime();
        elapsedLong /= 60d;
        elapsedLong /= 1000d;
        elapsedLong /= 1000d;
        int elapsed = (int) elapsedLong;
        if (elapsed > DET_ELAPSED_LOG_SECONDS) {
          double rawPct = 100d * (task.getUpdateCount()) / (task.getExpectedUpdateCount());
          String pct = (task.getIndeterminate() ? "" : " (" + ext.formDeci(rawPct, 2) + "%)");
          String msg = ext.getTime() + "]\tTask '" + task.getName() + "' with status '"
                       + task.getLabel() + "' is " + pct + " complete";
          if (task.getExpectedUpdateCount() > 100 && task.getUpdateCount() > 0) {
            if (internalLogger != null) {
              internalLogger.report(msg);
            } else {
              System.out.println(msg);
            }
          } else {
            if (internalLogger != null) {
              internalLogger.report(msg);
            } else {
              System.out.println(msg);
            }
          }
        }
      }

      task.updateDisplayTime();
    }
  }

  TimerTask monitorTask = new TimerTask() {
    @Override
    public void run() {
      if (taskUpdateStack.size() == 0) {
        return;
      }
      String[] keys = new String[taskUpdateStack.size()];
      taskUpdateStack.copyInto(keys);
      for (String str : keys) {
        Task task = taskMap.get(str);
        if (task == null) {
          continue; // probably removed before we got this far
        }
        synchronized (this) {
          int elapsed = (int) ((System.nanoTime() - task.getLastUpdate()) / (60 * 1000000000));
          if (elapsed > task.getTimeout() && !task.alreadyWarned()) {
            task.warned();
            String msg = "Alert! A long-running task has failed to progress in the past "
                         + task.getTimeout() + " minutes.";
            String[] opts = new String[] {"Continue", "Exit Genvisis"};
            int opt = JOptionPane.showOptionDialog(internalProgBar, msg, "Warning - stalled task!",
                                                   JOptionPane.YES_NO_OPTION,
                                                   JOptionPane.WARNING_MESSAGE, null, opts, 0);
            if (opt == 0 || opt == JOptionPane.CLOSED_OPTION) {
              continue;
            } else {
              System.exit(1);
            }
          }
        }
      }
      // TODO monitor jobs, warning for long running (could be programmer forgetting to endTask()!)
      // ProgressMonitor can't control running processes (pause or stop), but could work as a
      // System.exit()
    }
  };
}

