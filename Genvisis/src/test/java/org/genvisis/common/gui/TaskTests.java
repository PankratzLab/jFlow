package org.genvisis.common.gui;

import java.awt.Dimension;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextArea;

import org.genvisis.common.Logger;

/**
 * Interactive tests for {@link TaskManager}.
 * <p>
 * NB: not written as unit tests since these are just to demonstrate API use and confirm GUI
 * behavior.
 * </p>
 */
public class TaskTests {

	/**
	 * Ensure an empty progress bar reserves UI space appropriately
	 */
	private static void testProgressBar() {
		JProgressBarListener listener = new JProgressBarListener();
		makeFrame("testProgBar", listener);
	}

	/**
	 * Test the {@link TextListener} reporting to standard out. This test can be used to validate
	 * headless performance of the task framework.
	 */
	private static void testConsole() {
		GraphicsEnvironment.getLocalGraphicsEnvironment();
		// Note if this is running headlessly or not
		System.out.println("Headless? " + String.valueOf(GraphicsEnvironment.isHeadless()));
		String channel = "Text console";
		int count = 20;
		new TextListener(channel);
		Task task = new SimpleTask(channel, count) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count);
				return null;
			}
		};
		task.execute();
		try {
			// No UI so we have to wait for the task thread to terminate to avoid JVM termination
			task.get();
		} catch (InterruptedException | ExecutionException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Test the {@link TextListener} reporting to a {@link JTextArea}
	 */
	private static void testTextArea() {
		String channel = "Text testArea";
		int count = 20;
		JTextArea text = new JTextArea(25, 50);
		makeFrame("testTextArea", text);
		new TextListener(text, channel);
		Task task = new SimpleTask(channel, count) {
			@Override
			protected String doInBackground() throws Exception {
				for (int i = 0; i < count; i++) {
					doStep();
					text.append("stepped!\n");
					// Sleeping allows a smaller number of steps to show progress
					Thread.sleep(150);
				}
				return null;
			}
		};
		task.execute();
	}

	/**
	 * Test the {@link LoggerListener}.
	 */
	private static void testLogger() {
		String channel = "Log Test";
		int count = 20;
		JTextArea text = new JTextArea(25, 35);
		makeFrame("testLogger", text);
		Logger log = new Logger();
		log.linkTextArea(text);
		// NB: by default, the listener only reports to linked text areas. To report to console, pass
		// "true" flag to constructor.
		new LoggerListener(log, channel);
		Task task = new SimpleTask(channel, count) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count);
				return null;
			}
		};
		task.execute();
	}

	/**
	 * Test the {@link ProgressMonitorListener}. Optionally can test canceling tasks as well
	 */
	private static void testProgressMonitor() {
		String channel = "id";
		int count = 10;
		// Vacuous frame to keep the JVM alive
		makeFrame("testMonitor");
		// Create a progress monitor with no delay/min time
		// We don't need to add ProgressMonitorListeners to the UI so we don't need to
		// keep a reference
		// However, remember that the TaskManager keeps weak references to listeners.
		new ProgressMonitorListener(0, 0, channel);
		Task task = new SimpleTask(channel, count) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count);
				return null;
			}

			@Override
			public void done() {
				// The done method will run back on the EDT when the task is finished.
				if (isCancelled()) {
					// If cancelled (i.e. via the progressmonitor) show a message
					JOptionPane.showMessageDialog(null, "Cancelled!");
				}
			}
		};
		task.execute();
	}

	/**
	 * Test running a determinate and indeterminate task simultaneously on one Id
	 */
	private static void testMixedIndeterminate() {
		String channel = "id";
		int count1 = 10;
		int count2 = 8;
		JProgressBarListener listener = new JProgressBarListener(channel);
		makeFrame("testMixedInd", listener);
		Task task1 = new SimpleIndeterminateTask(channel) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count1);
				return null;
			}
		};
		Task task2 = new SimpleTask(channel, count2) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count2);
				return null;
			}
		};
		task1.execute();
		task2.execute();
	}

	/**
	 * Test running a progress bar task with no set length
	 */
	private static void testIndeterminate() {
		JProgressBarListener l1 = new JProgressBarListener();
		makeFrame("testInd", l1);
		String channel = "id";
		int count = 10;
		TaskManager.registerListener(channel, l1);
		Task task = new SimpleIndeterminateTask(channel) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count);
				return null;
			}
		};
		task.execute();
	}

	/**
	 * Test the {@link Task#addSteps(int)} functionality, where the max step count is increased while
	 * the task is running (kind of like a Windows file transfer time estimate).
	 */
	private static void testExpanding() {
		String channel = "id1";
		int count = 10;
		JProgressBarListener listener = new JProgressBarListener(channel);
		makeFrame("testExpand", listener);
		SimpleTask task = new SimpleTask(channel, 0) {
			@Override
			protected String doInBackground() throws Exception {
				for (int i = 0; i < count; i++) {
					addSteps(count);
					doSteps(this, count);
				}
				return null;
			}
		};
		task.execute();
	}

	/**
	 * Test creating a series of tasks up front and running each sequentially
	 */
	private static void testSequential() {
		final JProgressBarListener listener = new JProgressBarListener();
		final int count = 10;

		makeFrame("testSequential", listener);

		Thread t = new Thread(new Runnable() {
			@Override
			public void run() {
				List<Task> tasks = new ArrayList<>();

				// Create and register the tasks
				for (int i = 0; i < 10; i++) {
					String ch = String.valueOf(i);
					TaskManager.registerListener(ch, listener);
					tasks.add(new SimpleTask(ch, count) {
						@Override
						protected String doInBackground() throws Exception {
							doSteps(this, count);
							return null;
						}
					});
				}

				// Run each task sequentially
				for (Task t : tasks) {
					try {
						t.execute();
						t.get();
					} catch (InterruptedException | ExecutionException e) {
						e.printStackTrace();
					}
				}
			}
		});
		t.setDaemon(true);
		t.start();
	}

	/**
	 * Test creating a series of discrete tasks that are created after the previous task ends.
	 */
	private static void testDynamic() {
		final JProgressBarListener listener = new JProgressBarListener();
		final int count = 5;
		final List<String> channels = new ArrayList<>();

		// We define a set of tasks that our listener will list to
		for (int i = 0; i < 5; i++) {
			String ch = String.valueOf(i);
			channels.add(ch);
			TaskManager.registerListener(ch, listener);
		}

		// Create the window
		makeFrame("testDynamic", listener);

		// Because we want to create tasks after the previous task finishes, we have to
		// block until a
		// task ends. Thus we need to do task creation/running off the EDT
		Thread t = new Thread(new Runnable() {
			@Override
			public void run() {
				for (String ch : channels) {
					try {
						SimpleTask task = new SimpleTask(ch, count) {
							@Override
							protected String doInBackground() throws Exception {
								doSteps(this, count);
								return null;
							}

						};
						task.execute();
						task.get();
						// Wait for the progress bar to clear
						synchronized (task) {
							task.wait(JProgressBarListener.HIDE_DELAY + 500);
						}
					} catch (InterruptedException | ExecutionException e) {
						e.printStackTrace();
					}
				}
			}
		});
		// Start the thread as a daemon to avoid keeping the JVM alive if the window is
		// closed early.
		t.setDaemon(true);
		t.start();
	}

	/**
	 * Test running tasks in parallel with multiple listeners per task
	 */
	private static void testMulti() {
		String ch1 = "id1";
		String ch2 = "id2";
		int count1 = 40;
		int count2 = 10;

		// Listening to task 1 and 2
		JProgressBarListener listener12 = new JProgressBarListener(ch1, ch2);

		// Listening to task 1 only
		JProgressBarListener listener1 = new JProgressBarListener(ch1);

		// Listening to task 2 only
		JProgressBarListener listener2 = new JProgressBarListener(ch2);

		// Create a window to display the progress bars
		makeFrame("testMulti", listener12, listener1, listener2);

		// Create a task for each ID
		SimpleTask task1 = new SimpleTask(ch1, count1) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count1);
				return null;
			}
		};

		SimpleTask task2 = new SimpleTask(ch2, count2) {
			@Override
			protected String doInBackground() throws Exception {
				doSteps(this, count2);
				return null;
			}
		};

		// Run the tasks
		task1.execute();
		task2.execute();
	}

	/**
	 * Helper method to call the {@link Task#doStep()} method a number of times.
	 *
	 * @param t Task to perform steps
	 * @param count Number of steps to perform
	 */
	private static void doSteps(Task t, int count) throws InterruptedException {
		for (int i = 0; i < count; i++) {
			if (t.isCancelled()) {
				return;
			}
			t.doStep();
			// Sleeping allows a smaller number of steps to show progress
			Thread.sleep(150);
		}
	}

	/**
	 * Helper method to create and show a {@link JFrame} to display a number of components (e.g.
	 * {@link JProgressBar}s).
	 */
	private static void makeFrame(String title) {
		makeFrame(title, new JComponent[0]);
	}

	private static void makeFrame(String title, JProgressBarListener... listeners) {
		makeFrame(title, Arrays.stream(listeners).map(JProgressBarListener::getBar)
													 .toArray(JComponent[]::new));
	}

	private static void makeFrame(String title, JComponent... components) {
		JFrame window = new JFrame();
		window.setMinimumSize(new Dimension(200, 0));
		window.setTitle(title);

		JPanel pane = new JPanel();
		pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));
		for (JComponent c : components) {
			pane.add(c);
		}
		window.add(pane);
		window.pack();
		window.setLocationRelativeTo(null);
		window.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		window.setVisible(true);
	}

	/**
	 * Uncomment a test to run it
	 */
	public static void main(String... args) {
		// testMulti();
		// testDynamic();
		// testProgressBar();
		// testSequential();
		// testExpanding();
		// testIndeterminate();
		// testMixedIndeterminate();
		// testProgressMonitor();
		 testLogger();
//		 testTextArea();
//		testConsole();
	}
}
