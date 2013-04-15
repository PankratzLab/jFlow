package cnv.gui;

import java.io.File;

import javax.swing.JOptionPane;

import common.Files;

import cnv.filesys.ClusterFilterCollection;

public class AutoSaveClusterFilterCollection implements Runnable {

	private ClusterFilterCollection collection;
	private String filename;
	private int period;
	private boolean kill;
	
	public AutoSaveClusterFilterCollection(ClusterFilterCollection collection, String filename, int periodInSeconds) {
		this.collection = collection;
		this.filename = filename;
		this.period = periodInSeconds;
		this.kill = false;
	}

	public void saveNow() {
		collection.serialize(filename);
	}
	
	public void kill() {
		new File(filename).delete();
		kill = true;
	}
	
	public void run() {
		do {
			collection.serialize(filename);
			try {
				Thread.sleep(period*1000);
			} catch (InterruptedException ie) {}
		} while (!kill && Files.exists(filename, false));
		
		if (!kill) {
			JOptionPane.showMessageDialog(null, "Temporary file containing new ClusterFilters has been deleted; Auto-save has been turned off", "Error", JOptionPane.ERROR_MESSAGE);
		}
	}
}
