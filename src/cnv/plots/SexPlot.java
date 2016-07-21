// filter based on chromosome
package cnv.plots;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import javax.swing.*;

import common.*;
import cnv.filesys.Project;
import cnv.gui.ColorIcon;
import cnv.gui.WrapLayout;
import cnv.qc.SexChecks;

public class SexPlot extends JFrame {
	public static final long serialVersionUID = 1L;
	
	private static final String[] SEX_CHECKS_REQUIREMENTS = {"Sample", "Sex", SexChecks.EST_SEX_HEADER, "Check", "Median X LRR", "Median Y LRR", "Excluded"};
	
	SexPanel sexPanel;

	public SexPlot(Project proj, String[] samples, double[][] data, byte[] sexes, byte[] estimatedSexes, boolean[] excluded) {
		super("Genvisis - Sex Plot - " + proj.PROJECT_NAME.getValue());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		sexPanel = new SexPanel(proj, samples, data, sexes, estimatedSexes, excluded);
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
		
		setJMenuBar(createJMenuBar());

		// TODO extra paint appears to be unnecessary
//		repaint();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
//		unnecessary leads to a double rendering
//		sexPanel.createImage();
	}
	
	private JMenuBar createJMenuBar() {

		JMenuBar menuBar = new JMenuBar();
	
		JMenu viewMenu = new JMenu("View");
		viewMenu.setMnemonic(KeyEvent.VK_V);
		
		final JCheckBoxMenuItem showExcludedSwitch = new JCheckBoxMenuItem();
		
		AbstractAction showExcludedAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
            	sexPanel.setShowExcluded(showExcludedSwitch.isSelected());
                sexPanel.paintAgain();          
            }
        };
        showExcludedSwitch.setAction(showExcludedAction);
		showExcludedSwitch.setText("Show Excluded Samples");
		showExcludedSwitch.setMnemonic(KeyEvent.VK_X);
		viewMenu.add(showExcludedSwitch);
		
		menuBar.add(viewMenu);
		
		return menuBar;
	
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

	public static void loadSexCheckResults(Project proj) {
		Vector<String> samples = new Vector<String>();
		Vector<double[]> datapoints = new Vector<double[]>();
		Vector<Byte> sexes = new Vector<Byte>();
		Vector<Byte> estimatedSexes = new Vector<Byte>();
		Vector<Boolean> excluded = new Vector<Boolean>();
		try {
			BufferedReader reader = Files.getReader(proj.SEXCHECK_RESULTS_FILENAME.getValue(), proj.JAR_STATUS.getValue(), true, false);
			if (reader == null) {
				return;
			}
			String[] line = reader.readLine().trim().split("\t");
			if (!ext.checkHeader(line, SexChecks.SEX_HEADER, false, proj.getLog(), false)) {
				proj.message("The header in file '"+proj.SEXCHECK_RESULTS_FILENAME.getValue()+"' is not as expected and may cause problems; see log for more detail");
			}
			int[] indices = ext.indexFactors(SEX_CHECKS_REQUIREMENTS, SexChecks.SEX_HEADER, false, false);
			for (int index : indices) {
				if (index == -1) {
					return;
				}
			}
			
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t", -1);
				if (ext.isMissingValue(line[indices[4]]) || ext.isMissingValue(line[indices[5]]) ) {
					System.err.println("Error - sample '"+line[indices[0]]+"' does not have a valid medianLRR for X or Y");
				} else {
					samples.add(line[indices[0]]);
					datapoints.add(new double[] {Double.parseDouble(line[indices[4]]), Double.parseDouble(line[indices[5]])});
					sexes.add(Byte.parseByte(line[indices[1]]));
					estimatedSexes.add(Byte.parseByte(line[indices[2]]));
					excluded.add(line[indices[6]].equals("1"));
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.SEXCHECK_RESULTS_FILENAME.getValue()+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.SEXCHECK_RESULTS_FILENAME.getValue()+"\"");
			return;
		}
	
		new SexPlot(proj, samples.toArray(new String[samples.size()]), Matrix.toDoubleArrays(datapoints), Array.toByteArray(sexes), Array.toByteArray(estimatedSexes), Array.toBooleanArray(excluded));
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
			loadSexCheckResults(new Project(filename, false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
