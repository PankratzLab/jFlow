// -Xms1024M -Xmx1024M
package cnv.plots;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import common.*;
import filesys.SerialFloatArray;

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
	

	public QQPlot(String plotLabel, String[] labels, double[][] pvals, boolean log10, boolean rotated, boolean symmetric, float maxValue, Logger log) {
		super(plotLabel);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		log.report("Loading data for "+ext.listWithCommas(labels));

//		JPanel plotsPanel = new JPanel();
//		plotsPanel.setLayout(new GridLayout(1, 2));

		QQPanel panelA = new QQPanel(pvals, log10, rotated, maxValue);
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
			label = new JLabel("lambda = "+ext.formDeci(Array.lambda(pvals[i]), 4)+" ("+labels[i]+")", JLabel.CENTER);
			label.setForeground(pvals.length==1?COLOR_SCHEME[0]:COLOR_SCHEME[i+2]);
			label.setFont(new Font("Arial", 0, 20));
			descrPanel.add(label);
        }

		descrPanel.setBackground(Color.WHITE);
		getContentPane().add(descrPanel, BorderLayout.NORTH);

		repaint();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
//		I thought this createImage() led to an unnecessary double rendering, but then didn't get anything to view
//		However it is double rendering (at least in QQPlot) and removing createImage() removes both renderings
		panelA.createImage();
		panelA.updateUI(); // this adds the additional rendering, but if this isn't here then a partial image or no image is shown until the pane is moused over
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
	public static void loadPvals(String[] filenames, String plotLabel, boolean displayQuantiles, boolean displayStandardQQ, boolean displayRotatedQQ, double maxToPlot, boolean symmetric, float maxValue, Logger log) {
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
		int invalids;
		
		if (filenames == null || filenames.length == 0) {
			JOptionPane.showMessageDialog(null, "There are no files selected for viewing in QQ plot; please add the names of the files you want to view after QQ_FILENAMES= in the project's .properties file (also be sure to uncomment the property by removing the has symbol (\"#\"))", "Error", JOptionPane.ERROR_MESSAGE);
			log.reportError("There are no files selected for viewing in QQ plot; please add the names of the files you want to view after QQ_FILENAMES= in the project's .properties file (also be sure to uncomment the property by removing the has symbol (\"#\"))");
			return;
		}
		
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
					labels[i] += " '"+Files.getHeaderOfFile(filenames[i], log)[cols[i]]+"'";
				} catch (Exception e) {}
			}
			log.report("Loading "+labels[i]);
			if (!Files.exists(filenames[i])) {
				JOptionPane.showMessageDialog(null, "Error - file "+filenames[i]+" does not exist", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			delimiter = Files.determineDelimiter(filenames[i], log);
			try {
				reader = Files.getReader(filenames[i], JAR, true, true);
				count = 0;
				temp = reader.readLine();
				try {
					if (delimiter.equals(",")) {
						trav = ext.splitCommasIntelligently(temp, true, log)[cols[i]];
					} else {
						trav = temp.trim().split(delimiter, -1)[cols[i]];
					}
				} catch (Exception e) {
					log.reportError("Error - could not parse "+filenames[i]+" completely:");
					log.reportError(temp);
					log.reportException(e);
					return;
				}
				invalids = 0;
				try {
					if (!ext.isMissingValue(trav)) {
						if (Double.parseDouble(trav) <= 0) {
							JOptionPane.showMessageDialog(null, "Error - one of the p-values in file "+filenames[i]+" is near zero ("+trav+") for line:\n"+ext.replaceAllWith(temp, delimiter, "  "), "Error", JOptionPane.ERROR_MESSAGE);
							invalids++;
						}
						header = false;
						count++;
					} else {
						header = false;
					}
				} catch (NumberFormatException nfe) {
					header = true;
				}
				while (reader.ready()) {
					temp = reader.readLine();
					if (delimiter.equals(",")) {
						trav = ext.splitCommasIntelligently(temp, true, log)[cols[i]];
					} else {
						trav = temp.trim().split(delimiter, -1)[cols[i]];
					}
					if (!ext.isMissingValue(trav)) {
						if (Double.parseDouble(trav) <= 0) {
							if (invalids < 3) {
								JOptionPane.showMessageDialog(null, "Error - one of the p-values in file "+filenames[i]+" is near zero ("+trav+") for line:\n"+ext.replaceAllWith(temp, delimiter, "  "), "Error", JOptionPane.ERROR_MESSAGE);
							}
							invalids++;
						} else {
							count++;
						}
					}
				}
				reader.close();
				if (invalids > 2) {
					JOptionPane.showMessageDialog(null, "There were "+invalids+" total markers that had a p-value at or near zero for file "+filenames[i], "Error", JOptionPane.ERROR_MESSAGE);
				}

				reader = Files.getReader(filenames[i], JAR, true, true);
				pvals[i] = new double[count];
				count = 0;
				if (header) {
					reader.readLine();					
				}
				while (reader.ready()) {
					if (delimiter.equals(",")) {
						trav = ext.splitCommasIntelligently(reader.readLine(), true, log)[cols[i]];
					} else {
						trav = reader.readLine().trim().split(delimiter, -1)[cols[i]];
					}
					if (!ext.isMissingValue(trav)) {
						pvals[i][count] = Double.parseDouble(trav);
						if (pvals[i][count] > 0) {
							if (pvals[i][count] < minPval) {
								pvals[i][count] = minPval;
							}
							count++;
						}
					}
				}
				reader.close();
				if (count != pvals[i].length) {
					JOptionPane.showMessageDialog(null, "Error - mismatched number of values: "+count+" of "+pvals[i].length+" were valid", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
				
				pvals[i] = Sort.putInOrder(pvals[i]);
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error - missing file: \""+filenames[i]+"\"");
				error = false;
			} catch (IOException ioe) {
				log.reportError("Error reading file \""+filenames[i]+"\"");
				error = false;
			}
        }

		if (error) {
			return;
		}
		
		if (displayQuantiles) {
			new QQPlot(plotLabel, labels, pvals, false, false, false, maxValue, log);
		}
		if (displayStandardQQ) {
			new QQPlot(plotLabel, labels, pvals, true, false, symmetric, maxValue, log);
		}
		if (displayRotatedQQ) {
			new QQPlot(plotLabel, labels, pvals, true, true, false, maxValue, log);
		}
	}

	public static void computeCI(String dir, String prefix, int max, Logger log) {
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
		log.report("Found "+(count == 0?"no":count+"")+" replicates to process");
		dists = new float[0][0];
		length = -1;
		for (int i = 1; i <= count; i++) {
			if (i%10 == 0) {
				log.report(i+"");
			}
			array = SerialFloatArray.load(dir+prefix+"."+i+".results", false).getArray();
			order = Sort.quicksort(array);
			if (i == 1) {
				length = array.length;
				dists = new float[length][count];
			} else if (array.length != length) {
				log.reportError("Error - mismatched number of p-values at rep "+i);
			}
			for (int j = 0; j < array.length; j++) {
				dists[j][i-1] = array[order[j]];
				if (j == array.length-1 && (dists[j][i-1]+"").equalsIgnoreCase("NaN")) {
					log.reportError("Error - rep "+i+" has NaNs");
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
			log.reportError("Error writing to " + "CIs.xln");
			log.reportException(e);
		}
	}
	
	public static void main(String[] args) {
//		args = new String[] {"files=perfectPolymorphs.xln,5;perfectPolymorphs.xln,6;perfectPolymorphs.xln,12;perfectPolymorphs.xln,13;perfectPolymorphs.xln,19;perfectPolymorphs.xln,20;perfectPolymorphs.xln,26;perfectPolymorphs.xln,27"};
		args = new String[] {"files=perfectPolymorphs.xln,5;perfectPolymorphs.xln,6;perfectPolymorphs.xln,12"};

		int numArgs = args.length;
		String[] filenames = DEFAULT_FILES;
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

		String usage = "\n"+
		"plot.QQPlot requires 0-1 arguments\n"+
		"   (1) name of files with p-values (i.e. files="+Array.toStr(filenames, ";")+" (default))\n"+
		"   (2) -log10(p) at which to start truncating (i.e. maxToPlot=10 (default: -1))\n"+
		"   (3) make symmetric (i.e. -symmetric (not the default))\n"+
		"   (4) name of plot, for frame (i.e. plotLabel="+plotLabel+" (default))\n"+
		"   (5) maximum -log10 p-value to plot (i.e. maxValue=Infinity (default))\n"+
		"   (6) (optinal) log file (i.e. log="+logfile+" (default))\n"+
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
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		System.out.println("Found "+filenames.length+" files to parse: \n"+Array.toStr(filenames, "\n"));
		
		try {
			if (computePrefix != null) {
				computeCI(computeDir, computePrefix, max, new Logger(logfile));
			} else {
				loadPvals(filenames, plotLabel, displayQuantiles, displayStandardQQ, displayRotatedQQ, maxToPlot, symmetric, maxValue, new Logger(logfile));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
