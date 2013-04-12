// -Xms1024M -Xmx1024M
package cnv.plots;

import java.io.*;
import java.util.Vector;

import common.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class GenericPlot extends JFrame implements ActionListener {
	public static final long serialVersionUID = 1L;
	
	public static final String[] TRANSFORMS = {"none", "-log10", "ln"};
	public static final int N0_TRANSFORM = 0;
	public static final int NEG_LOG10_TRANSFORM = 1;
	public static final int NAT_LOG_TRANSFORM = 2;
//	public static final String DEFAULT_PARAMETERS = "xLabel=hg18_Pval,yLabel=hg19_Pval";
	
//	public static final String[] DEFAULT_FILES = {
//		"D:\\tWork\\Consortium\\Megas\\impRareComparisons.xln=CIDR,4:-log10,7:-log10",
//		"D:\\tWork\\Consortium\\Megas\\impRareComparisons.xln=Miami,10:-log10,13:-log10",
//		"D:\\tWork\\Consortium\\Megas\\impRareComparisons.xln=NGRC,16:-log10,19:-log10"
//	};

//	public static final String[] DEFAULT_FILES = {
//		"D:\\tWork\\Consortium\\Megas\\impCommonComparisons.xln=CIDR,4:-log10,7:-log10",
//		"D:\\tWork\\Consortium\\Megas\\impCommonComparisons.xln=Miami,10:-log10,13:-log10",
//		"D:\\tWork\\Consortium\\Megas\\impCommonComparisons.xln=NGRC,16:-log10,19:-log10"
//	};

	public static final String DEFAULT_PARAMETERS = "xLabel=rawX,yLabel=rawY";
	public static final String[] DEFAULT_FILES = {
		"D:\\tWork\\Consortium\\Megas\\impCommonComparisons.xln=CIDR,4,7",
//		"D:\\tWork\\Consortium\\Megas\\impCommonComparisons.xln=Miami,10:-log10,13:-log10",
//		"D:\\tWork\\Consortium\\Megas\\impCommonComparisons.xln=NGRC,16:-log10,19:-log10"
	};
	
	public static final boolean JAR = false;
	public static final Color[] COLOR_SCHEME = {Color.BLACK, Color.GRAY,
			new Color(55, 129, 252), // dark blue
			new Color(140, 20, 180), // deep purple
			new Color(201, 30, 10), new Color(217, 109, 194), // deep red/pink
			new Color(33, 31, 53), new Color(255, 255, 255), // dark dark / light light
			new Color(94, 88, 214), // light purple
			new Color(189, 243, 61), // light green
			new Color(33, 87, 0)}; // dark green
	

	public GenericPlot(String[] labels, float[][][] coords, String xLabel, String yLabel, float xMax, float yMax) {
		super("Plot");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		System.out.println("Loading data for "+ext.listWithCommas(labels));

		GenericPanel panelA = new GenericPanel(coords, xLabel, yLabel);
		panelA.setColorScheme(COLOR_SCHEME);
		if (xMax != Float.MIN_VALUE) {
			panelA.setForcePlotXmax(xMax);
		}
		if (yMax != Float.MIN_VALUE) {
			panelA.setForcePlotYmax(yMax);
		}
		getContentPane().add(panelA, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(labels.length, 1));

		JLabel label;
		
		for (int i = 0; i<labels.length; i++) {
			label = new JLabel(labels[i], JLabel.CENTER);
			label.setForeground(labels.length==1?COLOR_SCHEME[0]:COLOR_SCHEME[i+2]);
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
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();

		System.err.println("Error - unknown command '"+command+"'");
	}

	public static void loadPvals(String[] filenames, String parameters) {
		BufferedReader reader;
		String[] labels, line;
		float[][][] coords;
		float[] values;
		float xMax, yMax;
		String temp, filename, xLabel, yLabel;
		int count;
		String trav;
		boolean error;
		int[] cols, transforms;
		Vector<float[]> v;
		int estimatedNumberOfRecords;

		xLabel = "";
		yLabel = "";
		xMax = Float.MIN_VALUE;
		yMax = Float.MIN_VALUE;
		line = parameters.trim().split(",");
		for (int i = 0; i < line.length; i++) {
			if (line[i].startsWith("xMax=")) {
				xMax = ext.parseFloatArg(line[i]);
			}
			if (line[i].startsWith("yMax=")) {
				yMax = ext.parseFloatArg(line[i]);
			}
			if (line[i].startsWith("xLabel=")) {
				xLabel = ext.parseStringArg(line[i], "");
			}
			if (line[i].startsWith("yLabel=")) {
				yLabel = ext.parseStringArg(line[i], "");
			}			
		}
		
		error = false;
		labels = new String[filenames.length];
		coords = new float[filenames.length][][];
		
		for (int i = 0; i<filenames.length; i++) {
			line = filenames[i].trim().split(",");
			if (line.length != 1 && line.length != 3) {
				System.err.println("Error - each file set must either contain the x and y coordinates in the first two columns, or column indices for the x and y axes must be provided and must be separated by commas; file set '"+filenames[i]+"' contained "+line.length+" tokens");
				return;
			} else {
				if (line[0].indexOf("=") > 0) {
					labels[i] = line[0].substring(line[0].indexOf("=")+1);
					filename = line[0].substring(0, line[0].indexOf("="));
				} else {
					labels[i] = ext.rootOf(line[0]);
					filename = line[0];
				}				
				
				cols = new int[] {0,1};
				transforms = new int[] {0,0};
				for (int j = 0; j < 2; j++) {
					if (line[1+j].indexOf(":") > 0) {
						trav = line[1+j].substring(line[1+j].lastIndexOf(":")+1);
						line[1+j] = line[1+j].substring(0, line[1+j].lastIndexOf(":"));
						transforms[j] = ext.indexOfStr(trav, TRANSFORMS);
						if (transforms[j] == -1) {
							System.err.println("Error - invalid transform option ('"+trav+"'); valid options are "+ext.listWithCommas(TRANSFORMS));
							transforms[j] = 0;
						}
					}
					try {
						cols[j] = Integer.parseInt(line[1+j]);
					} catch (NumberFormatException nfe) {
						System.err.println("Error - invalid column number '"+line[1+j]+"' listed for the "+(j==0?"x":"y")+"-axis for file set '"+filenames[i]+"'");
						return;
					}						
				}
			}

			System.out.println("Loading "+labels[i]);
			try {
				reader = Files.getReader(filename, JAR, true, true);
				count = 0;
				v = null;
				while (reader.ready()) {
					temp = reader.readLine();
					line = temp.trim().split("[\\s]+");
					
					if (!ext.isMissingValue(line[cols[0]]) && !ext.isMissingValue(line[cols[1]])) {
						try {
							values = new float[] {Float.parseFloat(line[cols[0]]), Float.parseFloat(line[cols[1]])};
							for (int j = 0; j < transforms.length; j++) {
								if (transforms[j] == NEG_LOG10_TRANSFORM || transforms[j] == NAT_LOG_TRANSFORM) {
									if (values[j] <= 0) {
										System.err.println("Error - cannot perform a log transform on zero or negative numbers; pair ("+line[cols[0]]+", "+line[cols[0]]+" will be not be plotted");
									}
									if (transforms[j] == NEG_LOG10_TRANSFORM) {
										values[j] = (float)(-1*Math.log10(values[j]));
									} else if (transforms[j] == NAT_LOG_TRANSFORM) {
										values[j] = (float)Math.log(values[j]);
									}
								}
							}
							if (v == null) {
								estimatedNumberOfRecords = (int)((double)new File(filename).length() / (double)temp.length());
								System.out.println("Estimated number of records is "+estimatedNumberOfRecords+"; budgeting for "+(int)(estimatedNumberOfRecords*1.1));
								v = new Vector<float[]>((int)(estimatedNumberOfRecords*1.2));
							}
							v.add(values);
						} catch (NumberFormatException nfe) {}
					}
					count++;
				}
				reader.close();
				
				coords[i] = Matrix.toFloatArrays(v);

				System.out.println("Actual number of records was "+count+" of which "+coords[i].length+" were valid");
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
		new GenericPlot(labels, coords, xLabel, yLabel, xMax, yMax);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String[] filenames = DEFAULT_FILES;
		String parameters = DEFAULT_PARAMETERS;

		String usage = "\n"+
		"plots.GenericPlot requires 0-1 arguments\n"+
		"   (1) name of files with pairs of values (i.e. file="+Array.toStr(filenames, ",")+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			loadPvals(filenames, parameters);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
