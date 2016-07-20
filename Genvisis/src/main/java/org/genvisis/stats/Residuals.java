package org.genvisis.stats;

import org.genvisis.common.*;

public class Residuals {
	public static void createAndWinsorize(String filename, double sdThresholdForWinsorization, Logger log) {
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		Logger log;
		double sdThresholdForWinsorization = -1;

		String usage = "\n" +
		"stats.Residuals requires 1-2 arguments\n" +
		"   (1) filename of (i.e. file=pheno.csv (not the default))\n" +
		"   (2) (optional) threshold for Winsorization in standard deviation units (i.e. sdThresh=4.0 (not the default))\n" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			createAndWinsorize(filename, sdThresholdForWinsorization, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}