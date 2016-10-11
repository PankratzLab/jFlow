package org.genvisis.link.bat;

import java.io.File;
import java.io.IOException;
import java.util.Date;

import org.genvisis.common.ext;
import org.genvisis.link.LinkageFormat;
import org.genvisis.link.Recode;

public class SystemLinkage {
	public static void createAll() {
		long time;

		time = new Date().getTime();

		System.out.println(ext.getTime() + "\tstart");
		// new park.init.sordid();
		for (int i = 1; i <= 23; i++) {
			create("chromosome" + i + ".dat");
		}
		// new alleleDistribution("plateList.dat", true, 1, 23);
		Recode.recode("", 1, 23, 1);

		System.out.println(ext.getTime() + "\tstop");
		System.out.println("Total time: " + ext.getTimeElapsed(time));
	}

	public static void create(String filename) {
		int chr;

		filename = new File(filename).getName();

		if (filename.startsWith("chromosome")) {
			chr = -1;
			try {
				if (filename.length() == 15) {
					chr = Integer.valueOf(filename.substring(10, 11)).intValue();
				} else if (filename.length() == 16) {
					chr = Integer.valueOf(filename.substring(10, 12)).intValue();
				}
			} catch (NumberFormatException nfe) {
			}
			if (chr < 1 || chr > 23) {
				System.err.println("Error - There should only be 1 or 2 digits between chromosome and .dat, and it should be between 1 and 23");
			} else if (!new File("marker.database").exists()) {
				System.err.println("Error - could not find " + "marker.database" + " in current directory");
				System.err.println("        this is required to put markers in order for the map file");
			} else if (!new File("struct.dat").exists()) {
				System.err.println("Error - could not find " + "struct.dat" + " in current directory");
				System.err.println("        this is required to create a .pre file");
			} else {
				try {
					Filesystem.create(chr, false);
					LinkageFormat.create(chr, false);
				} catch (IOException ioe) {
					System.err.println("Error - these really should be dealt with locally...");
					ioe.printStackTrace();
				}
				return;
			}
		} else {
			System.err.println("Error - this only works on a chromosomeX.dat file");
			System.err.println("      - you have '"	+ filename + "' with a length of "
													+ filename.length());
		}

		ext.waitForResponse();
	}

	public static void main(String[] args) throws IOException {
		try {
			if (args.length == 0) {
				createAll();
			} else {
				create(args[0]);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
