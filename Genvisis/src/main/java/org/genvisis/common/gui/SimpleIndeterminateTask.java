package org.genvisis.common.gui;

/**
 * Trivial {@link IndeterminateTask} type-erasure extension that uses {@link String} for both of the
 * {@link SwingWorker} type parameters.
 *
 * @see {@link IndeterminateTask}
 * @see {@link SimpleTask}
 */
public abstract class SimpleIndeterminateTask extends IndeterminateTask<String, String> {

	/**
	 * @see {@link IndeterminateTask#IndeterminateTask(String)}
	 */
	public SimpleIndeterminateTask(String channel) {
		super(channel);
	}
}
