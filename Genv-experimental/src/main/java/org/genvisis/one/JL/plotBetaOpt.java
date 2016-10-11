package org.genvisis.one.JL;

/**
 * @author lane0212 See if I can get away from R...
 *
 */
public class plotBetaOpt {
	// https://docs.oracle.com/javafx/2/charts/bar-chart.htm

	public static void plotBetas(String file) {


	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename =
										"C:/data/misc/betaPlots/Whites_WBC_TOTAL_SingleSNPmatched.final_beta_summary.txt";

		String usage = "\n"	+ "one.JL.plotBetaOpt requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
