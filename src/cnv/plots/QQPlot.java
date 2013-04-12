// -Xms1024M -Xmx1024M
package cnv.plots;

import java.io.*;
import common.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import filesys.SerialFloatArray;

import stats.ProbDist;

public class QQPlot extends JFrame implements ActionListener {
	public static final long serialVersionUID = 1L;
	
	public static final String[] DEFAULT_FILES = {
//		"C:\\Documents and Settings\\npankrat\\My Documents\\Downloads\\allelic.txt"
//		"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\dominant.txt"
		"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\no_covariates.txt,0",
//		"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\4PCs.txt",
		"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\E2_E4.txt,0",
//		"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\4PCs_E2_E4.txt",
//		"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\E4_binary.txt"
	};

	public static final boolean JAR = false;
	public static final Color[] COLOR_SCHEME = {
		Color.BLACK, 
		Color.GRAY,
		new Color(55, 129, 252), // dark blue
		new Color(140, 20, 180), // deep purple
		new Color(201, 30, 10), new Color(217, 109, 194), // deep red/pink
		new Color(33, 31, 53), // new Color(255, 255, 255), // dark dark / light light
		new Color(94, 88, 214), // light purple
		new Color(189, 243, 61), // light green
		new Color(33, 87, 0),  // dark green

		new Color(45, 120, 150),
		new Color(240, 150, 100),
		new Color(100, 240, 0),

		};
	

	public QQPlot(String plotLabel, String[] labels, double[][] pvals, boolean log10, boolean rotated, boolean symmetric) {
		super(plotLabel);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		System.out.println("Loading data for "+ext.listWithCommas(labels));

//		JPanel plotsPanel = new JPanel();
//		plotsPanel.setLayout(new GridLayout(1, 2));

		QQPanel panelA = new QQPanel(pvals, log10, rotated);
		panelA.setSymmetricAxes(symmetric);
		panelA.setColorScheme(COLOR_SCHEME);
		if (!log10) {
			panelA.setForcePlotXmax(1);
			panelA.setForcePlotYmax(1);
		}
////		panelA.setPreferredSize(new Dimension(500,500));
//		plotsPanel.add(panelA);
		getContentPane().add(panelA, BorderLayout.CENTER);

//		QQPanel panelB = new QQPanel(pvals, true);
////		panelB.setPreferredSize(new Dimension(500,500));
//		plotsPanel.add(panelB);

//		getContentPane().add(plotsPanel, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(pvals.length, 1));

		JLabel label;
//		label = new JLabel("Q-Q plot", JLabel.CENTER);
//		label.setFont(new Font("Arial", 0, 20));
//		descrPanel.add(label);
		
		for (int i = 0; i<pvals.length; i++) {
			label = new JLabel("lambda = "+ext.formDeci(ProbDist.ChiDistReverse(Array.median(pvals[i]), 1)/ProbDist.ChiDistReverse(0.50, 1), 4)+" ("+labels[i]+")", JLabel.CENTER);
			label.setForeground(pvals.length==1?COLOR_SCHEME[0]:COLOR_SCHEME[i+2]);
			label.setFont(new Font("Arial", 0, 20));
			descrPanel.add(label);
        }

		descrPanel.setBackground(Color.WHITE);
		getContentPane().add(descrPanel, BorderLayout.NORTH);

		repaint();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
		panelA.createImage();
		panelA.updateUI();
//		panelB.createImage();
//		panelB.updateUI();
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();

		System.err.println("Error - unknown command '"+command+"'");
	}

//	public static String formatFilenames(String[] filenames, int[] indices, String[] labels) {
//		String str;
//		
//		if (filenames.length != indices.length || filenames.length != labels.length) {
//			System.err.println("Error - mismatched array lengths ("+filenames.length+"/"+indices.length+"/"+labels.length+")");
//			return null;
//		}
//		
//		str="";
//		for (int i = 0; i < filenames.length; i++) {
//			str += (i==0?"":";")+filenames[i]+","+indices[i]+"="+labels[i];
//		}
//		
//		return str;		
//	}
//	
	public static void loadPvals(String[] filenames, String plotLabel, boolean displayQuantiles, boolean displayStandardQQ, boolean displayRotatedQQ, double maxToPlot, boolean symmetric) {
		BufferedReader reader;
		String[] labels;
		double[][] pvals;
		int count;
		String trav;
		boolean error, header;
		int[] cols;
		double minPval;
		String delimiter;
		String temp;
		
		if (maxToPlot < 0) {
			minPval = -1;
		} else {
			minPval = 1/Math.pow(10, maxToPlot);
		}
		
		error = false;
		labels = new String[filenames.length];
		pvals = new double[filenames.length][];
		cols = Array.intArray(filenames.length, 0);
		for (int i = 0; i<filenames.length; i++) {
			if (filenames[i].indexOf("=") > 0) {
				labels[i] = filenames[i].substring(filenames[i].lastIndexOf("=")+1);
				filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf("="));
			} else {
				labels[i] = ext.rootOf(filenames[i]);
			}
			if (filenames[i].indexOf(",")  > 0) {
				try {
					cols[i] = Integer.parseInt(filenames[i].substring(filenames[i].lastIndexOf(",")+1));
					filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf(","));
				} catch (Exception e) {}
			}
			System.out.println("Loading "+labels[i]);
			delimiter = Files.determineDelimiter(filenames[i], new Logger());
			try {
				reader = Files.getReader(filenames[i], JAR, true, true);
				count = 0;
				temp = reader.readLine();
				try {
					trav = temp.trim().split(delimiter, -1)[cols[i]];
				} catch (Exception e) {
					System.err.println("Error - could not parse "+filenames[i]+" completely:");
					System.err.println(temp);
					e.printStackTrace();
					return;
				}
				try {
					if (!ext.isMissingValue(trav)) {
						Double.parseDouble(trav);
						header = false;
						count++;
					} else {
						header = false;
					}
				} catch (NumberFormatException nfe) {
					header = true;
				}
				while (reader.ready()) {
					trav = reader.readLine().trim().split(delimiter, -1)[cols[i]];
					if (!ext.isMissingValue(trav)) {
						count++;
					}
				}
				reader.close();

				reader = Files.getReader(filenames[i], JAR, true, true);
				pvals[i] = new double[count];
				count = 0;
				if (header) {
					reader.readLine();					
				}
				while (reader.ready()) {
					trav = reader.readLine().trim().split(delimiter, -1)[cols[i]];
					if (!ext.isMissingValue(trav)) {
						pvals[i][count] = Double.parseDouble(trav);
						if (pvals[i][count] < minPval) {
							pvals[i][count] = minPval;
						}
							
						count++;
					}
				}
				reader.close();
				
				pvals[i] = Sort.putInOrder(pvals[i]);
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error - missing file: \""+filenames[i]+"\"");
				error = false;
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+filenames[i]+"\"");
				error = false;
			}
        }

		if (error) {
			return;
		}
		
		if (displayQuantiles) {
			new QQPlot(plotLabel, labels, pvals, false, false, false);
		}
		if (displayStandardQQ) {
			new QQPlot(plotLabel, labels, pvals, true, false, symmetric);
		}
		if (displayRotatedQQ) {
			new QQPlot(plotLabel, labels, pvals, true, true, false);
		}
	}

	public static void computeCI(String dir, String prefix, int max) {
		PrintWriter writer;
		int count, length;
		float[] array;
		int[] order;
		float[][] dists;		
		
		if (max == -1) {
			max = Integer.MAX_VALUE;
		}
		count = 1;
		while (new File(dir+prefix+"."+count+".results").exists() && count < max) {
			count++;
		}
		count--;
		System.out.println("Found "+(count == 0?"no":count+"")+" replicates to process");
		dists = new float[0][0];
		length = -1;
		for (int i = 1; i <= count; i++) {
			if (i%10 == 0) {
				System.out.println(i);
			}
			array = SerialFloatArray.load(dir+prefix+"."+i+".results", false).getArray();
			order = Sort.quicksort(array);
			if (i == 1) {
				length = array.length;
				dists = new float[length][count];
			} else if (array.length != length) {
				System.err.println("Error - mismatched number of p-values at rep "+i);
			}
			for (int j = 0; j < array.length; j++) {
				dists[j][i-1] = array[order[j]];
				if (j == array.length-1 && (dists[j][i-1]+"").equalsIgnoreCase("NaN")) {
					System.err.println("Error - rep "+i+" has NaNs");
				}
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(dir+"CIs.xln"));
			writer.println("2.5%ile\t5%ile\t95%ile\t97.5%ile");
			for (int i = 0; i < dists.length; i++) {
				writer.println(Array.toStr(Array.quants(dists[i], new double[] {0.025, 0.05, 0.95, 0.975})));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "CIs.xln");
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String[] filenames = DEFAULT_FILES;
		String computeDir = "";
		String computePrefix = null;
		int max = -1;
		boolean displayQuantiles = false;
		boolean displayStandardQQ = true;
		boolean displayRotatedQQ = true;
		double maxToPlot = -1;
		boolean symmetric = false;
		String plotLabel = "Q-Q Plot";

		String usage = "\n"+
		"plot.QQPlot requires 0-1 arguments\n"+
		"   (1) name of files with p-values (i.e. files="+Array.toStr(filenames, ";")+" (default))\n"+
		"   (2) -log10(p) at which to start truncating (i.e. maxToPlot=10 (default: -1))\n"+
		"   (3) make symmetric (i.e. -symmetric (not the default))\n"+
		"   (4) name of plot, for frame (i.e. plotLabel="+plotLabel+" (default))\n"+
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
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		System.out.println("Found "+filenames.length+" files to parse: \n"+Array.toStr(filenames, "\n"));
		
		try {
			if (computePrefix != null) {
				computeCI(computeDir, computePrefix, max);
			} else {
				loadPvals(filenames, plotLabel, displayQuantiles, displayStandardQQ, displayRotatedQQ, maxToPlot, symmetric);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
