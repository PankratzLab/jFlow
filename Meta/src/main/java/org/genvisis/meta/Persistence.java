package org.genvisis.meta;

import java.util.HashSet;
import java.util.Set;
import javax.swing.SwingWorker;

/**
 * Helper class for keeping an application (the JVM) alive while processing things on what would
 * otherwise be daemon threads (e.g. via {@link SwingWorker})
 */
public final class Persistence {

  private static final Set<String> activeIDs = new HashSet<>();

  private Persistence() {
    // prevent creation of utility class
  }

  /**
   * @see #start(String)
   */
  public static void start(final Class<?> id) {
    start(id.getName());
  }

  /**
   * Add a block on JVM termination with the given ID. This must be removed with
   * {@link #stop(String)}.
   *
   * @param id ID to use when keeping the JVM alive
   */
  public static void start(final String id) {
    activeIDs.add(id);
    Thread daemonThread = new Thread(() -> wait(id));
    daemonThread.setDaemon(false);
    daemonThread.start();
  }

  /**
   * @see #stop(String)
   */
  public static void stop(final Class<?> id) {
    stop(id.getName());
  }

  /**
   * Remove the block on JVM termination with the given ID.
   *
   * @param id ID used in {@link #start(String)} method.
   */
  public static void stop(final String id) {
    activeIDs.remove(id);
  }

  private static void wait(final String id) {
    while (activeIDs.contains(id)) {
      try {
        Thread.sleep(50);
      } catch (InterruptedException e) {
        // Wait for signal
      }
    }
  }
}
