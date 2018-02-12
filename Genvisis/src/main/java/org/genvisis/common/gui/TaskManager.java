package org.genvisis.common.gui;

import java.lang.ref.WeakReference;
import java.util.Iterator;
import javax.swing.SwingWorker;
import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;

/**
 * Middle-layer utility method for managing computation-intensive {@link Task}s in a GUI
 * environment. This utility class is purely to facilitate {@code Task} and {@code TaskListener}s
 * linking up without direct knowledge of each other.
 * <p>
 * The general workflow is:
 * <ol>
 * <li>A {@link TaskListener} is created and registered to listen to {@code String}-identified
 * channel(s)</li>
 * <li>One or more {@link Task}s are created and each is registered to broadcast on particular
 * channel</li>
 * <li>Each {@code Task} is used as {@link SwingWorker}, performing work off-EDT.</li>
 * <li>When the {@link Task#doStep()} method is called, {@code TaskListener}s are updated as
 * appropriate.</li>
 * </ol>
 * </p>
 */
public final class TaskManager {

  /**
   * Map of listeners to keys. {@link WeakReference}s are used to avoid unintentionally preserving
   * listeners.
   */
  private static final Multimap<String, WeakReference<TaskListener>> listeners = MultimapBuilder.hashKeys()
                                                                                                .hashSetValues()
                                                                                                .build();

  private TaskManager() {
    // Prevent instantiation of static utility class
  }

  /**
   * Register a {@link TaskListener} to respond to property changes on the given channel
   */
  public static void registerListener(String channel, TaskListener abstractTaskListener) {
    listeners.put(channel, new WeakReference<TaskListener>(abstractTaskListener));
  }

  /**
   * Register a {@link Task} with any {@link TaskListener}s listening to the same channel.
   */
  public static <T, V> void registerTask(Task<T, V> task) {
    // Filter listeners, removing already collected weak references
    Iterator<WeakReference<TaskListener>> it = listeners.get(task.channel()).iterator();
    while (it.hasNext()) {
      WeakReference<TaskListener> next = it.next();
      TaskListener listener = next.get();
      if (listener == null) {
        it.remove();
      } else {
        // Register the task with each valid listener
        listener.listen(task);
      }
    }
  }
}
