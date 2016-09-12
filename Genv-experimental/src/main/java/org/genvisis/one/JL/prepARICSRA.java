package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;

import org.genvisis.common.Files;

public class prepARICSRA {

	private static void ff() {

	}

	private enum PLATFORM {
		/**
		 * That's whole exome
		 */
		WXS,
		/**
		 * Whole genome
		 */
		WGS
	}

	private static void prep(String sraRunTable) {
//		String[] extract = new String[]{}
		try {
			BufferedReader reader = Files.getAppropriateReader(sraRunTable);
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;

		String sraRunTable = "/Volumes/Work/data/aric_sra/prep/SraRunTable.txt";

		String usage = "\n" + "this requires 0-1 arguments\n" + "   (1) SRA data table (i.e. sraRunTable="
				+ sraRunTable + " (default))\n";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("sraRunTable=")) {
				sraRunTable = args[i].split("=")[1];
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
			prep(sraRunTable);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
