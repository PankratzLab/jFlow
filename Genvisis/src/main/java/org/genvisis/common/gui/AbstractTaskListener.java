package org.genvisis.common.gui;

import java.beans.PropertyChangeEvent;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.swing.SwingWorker.StateValue;

/**
 * Abstract superclass for {@link TaskListener} implementations. Provides boilerplate
 * implementations for common operations.
 */
public abstract class AbstractTaskListener implements TaskListener {

	private Map<Task, Integer> taskProgress = new HashMap<>();
	private Set<Task> running = new HashSet<>();

	/**
	 * @param channels This {@link TaskListener} will be
	 *        {@link TaskManager#registerListener(String, TaskListener)}ed on each of these.
	 */
	public AbstractTaskListener(String... channels) {
		for (String ch : channels) {
			TaskManager.registerListener(ch, this);
		}
	}

	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		// Ensure we have a SwingWorker source
		Object o = evt.getSource();
		if (!(o instanceof Task)) {
			// Not a SwingWorker
			return;
		}
		Task source = (Task) o;

		// Process event
		if ("progress".equals(evt.getPropertyName())) {
			taskProgress.put(source, (Integer) evt.getNewValue());
			int progress = 0;
			for (Integer p : taskProgress.values()) {
				progress += p;
			}
			// Compute overall effective progress for all tasks
			showProgress(progress / taskProgress.size());
		} else if (StateValue.STARTED.equals(evt.getNewValue())) {
			start();
		} else if (source.isCancelled() || StateValue.DONE.equals(evt.getNewValue())) {
			// Whether canceled or complete, this task is no longer running
			running.remove(source);
			source.removePropertyChangeListener(this);

			// Trigger the cancelled hook if appropriate
			// We want to call cancel-specific behavior before re-triggering start() below
			if (source.isCancelled()) {
				cancelled();
			}

			if (running.isEmpty()) {
				// If so, reset this listener and perform any post-task actions
				done();
				taskProgress.clear();
			} else {
				// Re-trigger start state, entering indeterminate mode until a progress event is received.
				start();
			}
		}
	}

	@Override
	public void doCancel() {
		for (Task t : running) {
			t.cancel(true);
		}
	}

	@Override
	public void listen(Task<?, ?> task) {
		task.addPropertyChangeListener(this);
		taskProgress.put(task, 0);
		running.add(task);
	}

	@Override
	public String channels() {
		return makeString(taskProgress.keySet(), Task::channel);
	}

	@Override
	public String activeChannels() {
		return makeString(running, Task::channel);
	}

	@Override
	public String ids() {
		return makeString(taskProgress.keySet(), Task::id);
	}

	@Override
	public String activeIDs() {
		return makeString(running, Task::id);
	}

	/**
	 * Helper method to generate a combined string from a {@link Task} method.
	 */
	private String makeString(Set<Task> tasks, Function<Task, String> toString) {
		return tasks.stream().map(t -> toString.apply(t)).distinct().sorted()
								.collect(Collectors.joining(", "));
	}
}
