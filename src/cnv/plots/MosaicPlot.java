// -Xms1024M -Xmx1024M
// filter based on chromosome
package cnv.plots;

import java.io.*;
import java.util.*;
import common.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import cnv.filesys.Project;

public class MosaicPlot extends JFrame implements ActionListener {
	public static final long serialVersionUID = 1L;
	public static final String[] MOSAICISM_HEADER = {"Sample", "Band", "LRR N", "mean LRR", "BAF N", "SD of BAF (0.15-0.85)", "IQR of BAF (0.15-0.85)", "%Homo"};

	public MosaicPlot(Project proj, String[][] samples, double[][] data) {
		super("Mosaicism Plot");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		MosaicPanel panel = new MosaicPanel(proj, samples, data);
		// panel.setToolTipText("");
		getContentPane().add(panel, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(2, 1));

		JLabel label = new JLabel("Mosaicism Plot", JLabel.CENTER);
		label.setFont(new Font("Arial", 0, 20));
		descrPanel.add(label);

		label = new JLabel("Only those B Allele Frequency (BAF) values between 0.15 and 0.85 are used in these calculations", JLabel.CENTER);
		label.setFont(new Font("Arial", 0, 14));
		descrPanel.add(label);
		descrPanel.setBackground(Color.WHITE);

		getContentPane().add(descrPanel, BorderLayout.NORTH);

		repaint();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
		panel.createImage();
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();

		System.err.println("Error - unknown command '"+command+"'");
	}

	public static void loadMosaicismResults(Project proj) {
		BufferedReader reader;
		String[] line;
		Vector<String[]> samples;
		Vector<double[]> datapoints;

		samples = new Vector<String[]>();
		datapoints = new Vector<double[]>();
		try {
			reader = Files.getReader(proj.getFilename(Project.MOSAIC_RESULTS_FILENAME), proj.getJarStatus(), true, true);
			ext.checkHeader(reader.readLine().trim().split("\t"), MOSAICISM_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!line[5].equals(".")&&!line[6].equals(".")&&Integer.parseInt(line[1].substring(3, line[1].length()-1))<23) {
					samples.add(new String[] {line[0], line[1]});
					datapoints.add(new double[] {Double.parseDouble(line[5]), Double.parseDouble(line[6])});
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getFilename(Project.MOSAIC_RESULTS_FILENAME)+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getFilename(Project.MOSAIC_RESULTS_FILENAME)+"\"");
			System.exit(2);
		}

		new MosaicPlot(proj, Matrix.toStringArrays(samples), Matrix.toDoubleArrays(datapoints));
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;

		String usage = "\n"+
		"plot.MosaicPlot requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
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
			loadMosaicismResults(new Project(filename, false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
