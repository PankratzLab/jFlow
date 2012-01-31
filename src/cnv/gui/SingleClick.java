package cnv.gui;

import java.awt.event.MouseEvent;

public class SingleClick implements Runnable {
	private boolean alive;
	private ClickListener listener;
	private int side;
	private int ms;

	public SingleClick(ClickListener listener, int side, int ms) {
		this.listener = listener;
		this.side = side;
		this.ms = ms;
		this.alive = true;
	}

	public void cancel() {
		alive = false;
	}

	public boolean isAlive() {
		return alive;
	}

	public void run() {
		try {
			Thread.sleep(ms);
			if (alive) {
				if (side==MouseEvent.BUTTON1) {
					listener.singleLeftClick();
				}
				if (side==MouseEvent.BUTTON3) {
					listener.singleRightClick();
				}
			}
			alive = false;
			
		} catch (InterruptedException ie) {
			ie.printStackTrace();
		}

	}
}
