package org.genvisis.common.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.swing.JProgressBar;
import javax.swing.Timer;

/**
 * {@link TaskListener} that wraps a {@link JProgressBar}
 */
public class JProgressBarListener extends AbstractTaskListener {

	private final boolean hideWhenDone;
	private final JProgressBar progressBar;

	// Timer to hide and reset the progress bar after a delay
	private final Timer timer = new Timer(HIDE_DELAY, new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			progressBar.setValue(0);
			progressBar.setVisible(false);
		}
	});

	public JProgressBarListener(String... channels) {
		this(true, channels);
	}

	/**
	 * Usage note: do not call <code>setVisible(false)</code> on the internal JProgressBar within a
	 * few lines of instantiation, as the internal component listener will not fire, and the first
	 * time the component is set to visible it will be rehidden automatically.
	 * 
	 * @see {@link AbstractTaskListener#AbstractTaskListener(String...)}
	 */
	public JProgressBarListener(boolean hideWhenDone, String... channels) {
		super(channels);
		this.hideWhenDone = hideWhenDone;
		// #showProgress is scaled from 0 to 100 by SwingWorker.
		progressBar = new JProgressBar(0, 100);
		// We make the progress bar visible at the beginning to be counted for real estate when
		progressBar.setVisible(true);

		// The first time the UI is actually shown, we hide it.
		progressBar.addComponentListener(new ComponentListener() {

			@Override
			public void componentShown(ComponentEvent e) {
				hide();
			}

			@Override
			public void componentResized(ComponentEvent e) {
				hide();
			}

			@Override
			public void componentMoved(ComponentEvent e) {
				// No-op
			}

			@Override
			public void componentHidden(ComponentEvent e) {
				// No-op
			}

			/**
			 * Helper method to do a one-time hide after packing
			 */
			private void hide() {
				progressBar.setVisible(false);
				progressBar.removeComponentListener(this);
			}
		});
	}

	/**
	 * Get the progress bar. NB: should only be used for UI sizing and placement.
	 *
	 * @return The backing {@link JProgressBar}
	 */
	public JProgressBar getBar() {
		return progressBar;
	}

	@Override
	public void showProgress(int progress) {
		// Ensure indeterminate state is cleared
		progressBar.setIndeterminate(false);
		progressBar.setValue(progress);
	}

	@Override
	public void cancelled() {}

	@Override
	public void start() {
		// Start in indeterminate mode
		progressBar.setIndeterminate(true);
		// Ensure the bar is visible
		progressBar.setVisible(true);
		// Stop any existing hide timer
		timer.stop();
	}

	@Override
	public void done() {
		// Set progress bar to full
		progressBar.setValue(progressBar.getMaximum());
		// Ensure indeterminate state is cleared
		progressBar.setIndeterminate(false);
		// Start the hide timer
		if (hideWhenDone) {
			timer.restart();
		}
	}
}
