// filter based on chromosome
package cnv.plots;

import java.io.*;
import java.util.*;
import java.awt.*;

import javax.swing.*;

import common.*;
import cnv.filesys.Project;
import cnv.gui.ColorIcon;
import cnv.gui.WrapLayout;
import cnv.qc.SexChecks;

public class SexPlot extends JFrame {
	public static final long serialVersionUID = 1L;

	SexPanel sexPanel;

	public SexPlot(Project proj, String[][] samples, double[][] data) {
		this(proj, samples, data, null, null);
	}

	public SexPlot(Project proj, String[][] samples, double[][] data, byte[] sexes, byte[] estimatedSexes) {
		super("Genvisis - Sex Plot - " + proj.getNameOfProject());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		SexPanel sexPanel = new SexPanel(proj, samples, data, sexes, estimatedSexes);
		// panel.setToolTipText("");
		getContentPane().add(sexPanel, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(2, 1));

		JLabel label = new JLabel("Sex Plot", JLabel.CENTER);
		label.setFont(new Font("Arial", 0, 20));
		descrPanel.add(label);
		descrPanel.setBackground(Color.WHITE);

		getContentPane().add(descrPanel, BorderLayout.NORTH);

		JPanel legendPanel = colorLegendPanel();
		getContentPane().add(legendPanel, BorderLayout.SOUTH);

		// TODO extra paint appears to be unnecessary
//		repaint();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
//		unnecessary leads to a double rendering
//		sexPanel.createImage();
	}

	private JPanel colorLegendPanel() {
		JPanel colorLegendPanel = new JPanel(new WrapLayout(FlowLayout.CENTER, 0, 0));
		colorLegendPanel.setBackground(Color.WHITE);//BACKGROUND_COLOR);
	
	    JLabel legend = new JLabel("Color Key: ");
		legend.setFont(new Font("Arial", 0, 14));
		colorLegendPanel.add(legend);
		
		for (int i=0; i<SexPanel.COLOR_SCHEME.length; i++){
			JPanel enclosure = new JPanel();
			enclosure.setBackground(Color.WHITE);
			enclosure.add(new JLabel(new ColorIcon(12,12,SexPanel.COLOR_SCHEME[i])));
			enclosure.add(new JLabel(SexPanel.COLOR_SCHEME_MEANING[i]));
			colorLegendPanel.add(enclosure);
		}
		return colorLegendPanel;
	}

	public static void loadGenderResults(Project proj) {
		BufferedReader reader;
		String[] line;
		Vector<String[]> samples;
		Vector<double[]> datapoints;
		Vector<Byte> sexes;
		Vector<Byte> estimatedSexes;
	
		samples = new Vector<String[]>();
		datapoints = new Vector<double[]>();
		sexes = new Vector<Byte>();
		estimatedSexes = new Vector<Byte>();
		try {
			reader = Files.getReader(proj.PROJECT_DIRECTORY.getValue()+"sexCheck.xln", proj.JAR_STATUS.getValue(), true, false);
			if (reader == null) {
				return;
			}
			ext.checkHeader(reader.readLine().trim().split("\t"), SexChecks.SEX_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t", -1);
				if (ext.isMissingValue(line[4]) || ext.isMissingValue(line[8]) ) {
					System.err.println("Error - sample '"+line[0]+"' does not have a valid meanLRR for X or Y");
				} else {
					samples.add(new String[] {line[0], "chr23"});
					datapoints.add(new double[] {Double.parseDouble(line[5]), Double.parseDouble(line[9])});
					sexes.add(Byte.parseByte(line[3]));
					estimatedSexes.add(Byte.parseByte(line[4]));
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.PROJECT_DIRECTORY.getValue()+"sexCheck.xln"+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.PROJECT_DIRECTORY.getValue()+"sexCheck.xln"+"\"");
			return;
		}
	
		new SexPlot(proj, Matrix.toStringArrays(samples), Matrix.toDoubleArrays(datapoints), Array.toByteArray(sexes), Array.toByteArray(estimatedSexes));
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;

		String usage = "\n"+
		"cnv.plot.SexPlot requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			loadGenderResults(new Project(filename, false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
