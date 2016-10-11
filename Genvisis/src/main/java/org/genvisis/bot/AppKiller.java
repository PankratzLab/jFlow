package org.genvisis.bot;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Date;

import org.genvisis.common.ext;

public class AppKiller implements Runnable {
	private final String plugFile;

	public AppKiller(String plugFile) {
		this.plugFile = plugFile;

		try {
			new PrintWriter(new FileWriter(plugFile)).close();
		} catch (Exception e) {
		}
	}

	@Override
	public void run() {
		long time;

		time = new Date().getTime();
		while (new File(plugFile).exists()) {
			try {
				Thread.sleep(1000);
			} catch (InterruptedException ie) {
			}
		}
		System.err.println("Plug was pulled after " + ext.getTimeElapsed(time));
		System.exit(1);
	}
}
