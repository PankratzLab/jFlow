package one.JL;

/**
 * @author lane0212 See if I can get away from R...
 *
 */
public class plotBetaOpt {
	// https://docs.oracle.com/javafx/2/charts/bar-chart.htm
	
	public static void plotBetas(String file){
		
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/data/misc/betaPlots/Whites_WBC_TOTAL_SingleSNPmatched.final_beta_summary.txt";

		String usage = "\n" +
				"one.JL.plotBetaOpt requires 0-1 arguments\n" +
				"   (1) filename (i.e. file=" + filename + " (default))\n" +
				"";

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
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
