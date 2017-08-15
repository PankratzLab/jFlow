package org.genvisis.common.gui;

import javax.swing.JTextArea;

/**
 * {@link TaskListener} that prints a simple string markup of a progress bar. Can write to either a
 * {@link JTextArea} or system out. Will attempt to update a single line as progress is reported -
 * by caret tracking, a {@code JTextArea} can continue to be updated by other sources. Consoles
 * should expect to be less robust.
 * <p>
 * NB: Console output uses a carriage return, which may not be respected in some consoles (e.g.
 * Eclipse)
 * </p>
 * <p>
 * NB2: this is not a very elegant way to manage both swing and console output. Not clear how much
 * of this should be the role of a logger.
 * </p>
 */
public class TextListener extends AbstractTaskListener {

	// Fields for TextArea management
	private final JTextArea textArea;
	private int caret = -1;
	private int endLine = -1;

	/**
	 * Create a {@code TaskProgressListener} that reports to {@link System#out}
	 *
	 * @see #TextListener(JTextArea, String...)
	 */
	public TextListener(String... channels) {
		this(null, channels);
	}

	/**
	 * Create a {@code TextListener} that reports to {@link JTextArea}
	 *
	 * @param textArea Text area to report progress
	 * @param channels see {@link AbstractTaskListener#AbstractTaskListener(String...)}
	 */
	public TextListener(JTextArea textArea, String... channels) {
		super(channels);
		this.textArea = textArea;
	}

	@Override
	public void showProgress(int progress) {
		writeProgressString(progress);
	}

	/**
	 * Write an intermediate progress message (not the end of reporting)
	 *
	 * @see #writeProgressString(int, boolean)
	 */
	private void writeProgressString(int progress) {
		writeProgressString(progress, false);
	}

	/**
	 * Write a progress message
	 *
	 * @param progress Progress value
	 * @param end If true, add appropriate end-of-message text
	 */
	private void writeProgressString(int progress, boolean end) {
		StringBuilder sb = new StringBuilder("ID(s): " + ids() + "; CH(s): " + channels()
																				 + "; ");

		// Fill out completed progress units
		for (int i = 0; i < progress / 10; i++) {
			sb.append("=");
		}

		// Add a half unit
		if (progress % 10 >= 5) {
			sb.append("-");
		}

		// Add end message
		if (end) {
			sb.append(" done!");
		}

		if (textArea != null) {
			// Update text area
			sb.append("\n");

			String line = sb.toString();

			if (endLine < 0) {
				textArea.insert(line, caret);
			} else {
				textArea.replaceRange(line, caret, endLine);
			}
			textArea.setCaretPosition(textArea.getDocument().getLength());

			endLine = caret + line.length();
		} else {
			// Report to sysout
			sb.append(end ? '\n' : '\r');
			System.out.print(sb.toString());
		}
	}

	@Override
	public void cancelled() {
		// Nothing to do here
	}

	@Override
	public void start() {
		if (textArea != null && caret == -1) {
			caret = textArea.getCaretPosition();
		}
	}

	@Override
	public void done() {
		// Write a complete line, terminated with a \n
		writeProgressString(100, true);
	}
}
