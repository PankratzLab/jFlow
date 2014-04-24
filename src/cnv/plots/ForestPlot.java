package cnv.plots;

import java.util.*;

import common.*;

public class ForestPlot {

	private static void createForest(String plotLabel, Vector<ForestTree> trees, Logger log) {
		// TODO Auto-generated method stub
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "ForestPlot.dat";
		String logfile = null;
		Logger log;

		String usage = "\n" + "cnv.plots.ForestPlot requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

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
			
			Vector<ForestTree> trees;
			
			trees = new Vector<ForestTree>();
			trees.add(new ForestTree("PROGENI/GenePD", 0.296154, 0.0834038, 0, 0));
			trees.add(new ForestTree("NGRC", 0.105856, 0.0559677, 0, 0));
			trees.add(new ForestTree("23andMe", 0.213202, 0.027064, 0, 0));
			trees.add(new ForestTree("Summary", 0.156191625, 0.022027313, 0, 1));
			
			log = new Logger(logfile);
			createForest("rs12338", trees, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class ForestTree {
	public String label;
	public double beta;
	public double stderr;
	public int color;
	public int shape;
	
	public ForestTree(String label, double beta, double stderr, int color, int shape) {
		this.label = label;
		this.beta = beta;
		this.stderr = stderr;
		this.color = color;
		this.shape = shape;
	}

	public String getLabel() {
		return label;
	}

	public double getBeta() {
		return beta;
	}

	public double getStderr() {
		return stderr;
	}

	public int getColor() {
		return color;
	}

	public int getShape() {
		return shape;
	}
}
