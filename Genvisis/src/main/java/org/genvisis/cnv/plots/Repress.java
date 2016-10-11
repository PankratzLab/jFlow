package org.genvisis.cnv.plots;

public class Repress implements Runnable {
	private final MosaicPanel panel;
	private boolean alive;
	private final int millis;

	public Repress(MosaicPanel panel, int millis) {
		this.panel = panel;
		this.millis = millis;
		alive = true;
	}

	public void cancel() {
		alive = false;
	}

	@Override
	public void run() {
		try {
			System.out.println("Repressed");
			panel.setFlow(false);
			Thread.sleep(millis);
			if (alive) {
				panel.createImage();
			}
		} catch (InterruptedException ie) {
			ie.printStackTrace();
		}

	}
}
