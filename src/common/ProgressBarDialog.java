package common;

import java.awt.*;
import javax.swing.*;

//TODO change to a SwingWorker so that this updates in realtime (especially for QQPlot) http://www.javacreed.com/swing-worker-example/
public class ProgressBarDialog extends JDialog {
	public static final long serialVersionUID = 1L;

	private JProgressBar pb;
	private long timeStarted;
	private long timeDelay;

	public ProgressBarDialog(String frameText, int min, int max, int width, int height) {
		this(frameText, min, max, width, height, 0);
	}
	
	public ProgressBarDialog(String frameText, int min, int max, int width, int height, int timeDelay) {
		super((JFrame)null, frameText);
		
		this.timeDelay = timeDelay;
		timeStarted = System.currentTimeMillis();

		pb = new JProgressBar(min, max);
		pb.setPreferredSize(new Dimension(175, 20));
		pb.setStringPainted(true);
		pb.setValue(0);

		JPanel center_panel = new JPanel();
		center_panel.add(pb);

		getContentPane().add(center_panel, BorderLayout.CENTER);
		pack();
		setVisible(true);
		setLocation(width/2-87, height/2-10);
		// setLocationRelativeTo(null); // center on screen

		toFront(); // raise above other java windows
	}
	
	public void setTimeDelay(long delay) {
		this.timeDelay = delay;
	}

	public void setProgress(int value) {
		if (System.currentTimeMillis() - timeStarted > timeDelay) {
			pb.setValue(value);
			pb.setString((int)((double)value/(double)(pb.getMaximum()-pb.getMinimum())*100)+"%");
		}
	}

	public void close() {
		dispose();
	}
}
