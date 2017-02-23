package org.genvisis.common;

import java.util.ArrayList;
import java.util.List;

/**
 * Abstract superclass for {@link StartupCheck} implementations. Use {@link #addMessage(String)} to
 * mitigate the need for explicit management of error messages.
 */
public abstract class AbstractStartupCheck implements StartupCheck {

	private List<String> messages = new ArrayList<String>();

	@Override
	public List<String> check() {
		doCheck();
		return messages;
	}

	/**
	 * Build the warning messages for this startup check. Ensures {@link #warningHeader()} is the
	 * first message, and all subsequent items are indented.
	 *
	 * @param msg String to append to this {@code StartupCheck}'s warning message.
	 */
	protected void addMessage(String msg) {
		if (messages.isEmpty()) {
			messages.add(warningHeader());
		}
		messages.add("\t" + msg);
	}

	protected abstract String warningHeader();

	protected abstract void doCheck();
}
