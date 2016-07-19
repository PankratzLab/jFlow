// -Xms1024M -Xmx1024M
package cnv.plots;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import common.*;

public class QQPlotFrame extends JFrame implements ActionListener {
	public static final long serialVersionUID = 1L;
	
	/**
	 * 
	 * @param plotLabel label for the plot
	 * @param labels labels for each set of pvals
	 * @param pvals array of pvals for each label in labels
	 * @param log10 use -log10 of p values
	 * @param rotated display QQ rotated
	 * @param symmetric force symmetric axes
	 * @param maxValue maximum -log10 pval to display
	 * @param log
	 */
	public QQPlotFrame(String plotLabel, QQPlot qqPlot) {
		super(plotLabel);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		getContentPane().add(qqPlot.getQqPanel(), BorderLayout.CENTER);

		setBounds(20, 20, 1000, 720);
		setVisible(true);
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();

		System.err.println("Error - unknown command '"+command+"'");
	}


	/**
	 * 
	 * @param filenames files to load pvals from
	 * @param plotLabel label for the plot
	 * @param displayQuantiles display quantiles instead of -log10 pvals
	 * @param displayStandardQQ display standard QQ plot with -log10 pvals
	 * @param displayRotatedQQ display rotated QQ plot with -log10 pvals
	 * @param maxToPlot -log10(p) at which to start truncating
	 * @param symmetric force symmetric axes
	 * @param maxValue maximum -log10 p-value to plot
	 * @param log
	 */
	public static void loadPvals(String[] filenames, String plotLabel, boolean displayQuantiles, boolean displayStandardQQ, boolean displayRotatedQQ, double maxToPlot, boolean symmetric, float maxValue, Logger log) {
		QQPlot qqPlot = QQPlot.loadPvals(filenames, plotLabel, displayQuantiles, displayStandardQQ, displayRotatedQQ, maxToPlot, symmetric, maxValue, log);
		if (qqPlot == null) return;
		new QQPlotFrame(plotLabel, qqPlot);
	}
	
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String[] filenames = QQPlot.DEFAULT_FILES;
		String computeDir = "";
		String computePrefix = null;
		int max = -1;
		boolean displayQuantiles = false;
		boolean displayStandardQQ = true;
		boolean displayRotatedQQ = false;
		double maxToPlot = -1;
		boolean symmetric = false;
		String plotLabel = "Q-Q Plot";
		float maxValue = Float.MAX_VALUE;
		String logfile = null;
		String outFile = null;

		String usage = "\n"+
		"plot.QQPlot requires 0-1 arguments\n"+
		"   (1) name of files with p-values (i.e. files="+Array.toStr(filenames, ";")+" (default))\n"+
		"   (2) -log10(p) at which to start truncating (i.e. maxToPlot=10 (default: -1))\n"+
		"   (3) make symmetric (i.e. -symmetric (not the default))\n"+
		"   (4) name of plot, for frame (i.e. plotLabel="+plotLabel+" (default))\n"+
		"   (5) maximum -log10 p-value to plot (i.e. maxValue=Infinity (default))\n"+
		"   (6) (optional) log file (i.e. log="+logfile+" (default))\n"+
		"   (7) (optional) file to export image to instead of displaying (i.e. outFile="+outFile+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("files=")) {
				filenames = args[i].substring(6).split(";");
				numArgs--;
			} else if (args[i].startsWith("prefix=")) {
				computePrefix = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("max=")) {
				max = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("maxToPlot=")) {
				maxToPlot = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-quantiles")) {
				displayQuantiles = true;
				numArgs--;
			} else if (args[i].startsWith("-standardQQ")) {
				displayStandardQQ = true;
				numArgs--;
			} else if (args[i].startsWith("-rotatedQQ")) {
				displayRotatedQQ = true;
				numArgs--;
			} else if (args[i].startsWith("-symmetric")) {
				symmetric = true;
				numArgs--;
			} else if (args[i].startsWith("plotLabel=")) {
				plotLabel = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("maxValue=")) {
				maxValue = ext.parseFloatArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("outFile=")) {
				outFile = ext.parseStringArg(args[i], null);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		System.out.println("Found "+filenames.length+" files to parse: \n"+Array.toStr(filenames, "\n"));
		
		try {
			if (computePrefix != null) {
				QQPlot.computeCI(computeDir, computePrefix, max, new Logger(logfile));
			} else if (outFile != null) {
				QQPlot qqPlot = QQPlot.loadPvals(filenames, plotLabel, displayQuantiles, displayStandardQQ, displayRotatedQQ, maxToPlot, symmetric, maxValue, new Logger(logfile));
	    		qqPlot.screenCap(outFile);
			} else {
				loadPvals(filenames, plotLabel, displayQuantiles, displayStandardQQ, displayRotatedQQ, maxToPlot, symmetric, maxValue, new Logger(logfile));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
