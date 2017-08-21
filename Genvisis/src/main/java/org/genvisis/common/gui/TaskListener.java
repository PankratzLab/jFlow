package org.genvisis.common.gui;

import java.beans.PropertyChangeListener;

/**
 * Interface for listening to one or more {@link Task}s managed by the {@link TaskManager}
 */
public interface TaskListener extends PropertyChangeListener {

	/**
	 * Time (in milliseconds) for a listener to wait after completion of all tasks. If a new task is
	 * started before this delay it will be merged with the current step counts.
	 */
	public static final int HIDE_DELAY = 2000;

	/**
	 * All {@link TaskListeners} should use this min progress
	 */
	public static final int MIN_PROGRESS = 0;

	/**
	 * All {@link TaskListeners} should use this max progress
	 */
	public static final int MAX_PROGRESS = 100;

	/**
	 * @param task Start listening to the given {@link Task}
	 */
	void listen(Task<?, ?> task);

	/**
	 * Update the UI with the current progress (also implicitly ends any indeterminate rendering)
	 *
	 * @param Progress value (from 0 to 100)
	 */
	void showProgress(int progress);

	/**
	 * Forcibly cancel <b>all</b> {@link Task}s associated with this listener.
	 */
	void doCancel();

	/**
	 * Hook that is called when a {@link Task} is cancelled.
	 */
	void cancelled();

	/**
	 * Hook for when a task is starting. NB: Should start in indeterminate mode (see
	 * {@link IndeterminateTask}) if supported.
	 */
	void start();

	/**
	 * Hook for when all tasks are complete, whether canceled or finished naturally. The listener
	 * should reset after {@link #HIDE_DELAY}.
	 */
	void done();

	/**
	 * @return Combined string of {@link Task#channel()}s for all registered tasks
	 */
	String channels();

	/**
	 * @return Combined string of {@link Task#channel()}s for all registered tasks that are currently running (not finished/cancelled)
	 */
	String activeChannels();

	/**
	 * @return Combined string of {@link Task#id()}s for all registered tasks
	 */
	String ids();

	/**
	 * @return Combined string of {@link Task#id()}s for all registered tasks that are currently running (not finished/cancelled)
	 */
	String activeIDs();
}
