package org.genvisis.one;

public class CARDIA_Analyses {
	
	private static void makePhenosFromClipboard() {
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "CARDIA_Analyses.dat";

		String usage = "\n" + "one.CARDIA_Analyses requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			makePhenosFromClipboard();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
