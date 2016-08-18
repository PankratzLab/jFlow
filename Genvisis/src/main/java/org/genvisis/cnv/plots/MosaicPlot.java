// -Xms1024M -Xmx1024M
// filter based on chromosome
package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Vector;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Files;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class MosaicPlot extends JFrame implements ActionListener {
  public static final long serialVersionUID = 1L;
  public static final String[] MOSAICISM_HEADER =
      {"Sample", "Arm", "LRR N", "mean LRR", "BAF N", "SD of BAF (0.15-0.85)",
       "IQR of BAF (0.15-0.85)", "%Homo", "ForcedCallArmPercentMosaicism", "BpWeightedAverageArm",
       "BpWeightedAverageCalled", "NumberRegionsDetected", "BpCalledMosaic", "BpInArm",
       "ProportionArmCalledMosaic"};

  MosaicPanel panel;

  public MosaicPlot(Project proj, String[][] samples, double[][] data) {
    super("Genvisis - Mosaicism Plot - " + proj.PROJECT_NAME.getValue());
    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

    panel = new MosaicPanel(proj, samples, data);
    // panel.setToolTipText("");
    getContentPane().add(panel, BorderLayout.CENTER);

    JPanel descrPanel = new JPanel();
    descrPanel.setLayout(new GridLayout(2, 1));

    JLabel label = new JLabel("Mosaicism Plot", JLabel.CENTER);
    label.setFont(new Font("Arial", 0, 20));
    descrPanel.add(label);

    label =
        new JLabel("Only those B Allele Frequency (BAF) values between 0.15 and 0.85 are used in these calculations",
                   JLabel.CENTER);
    label.setFont(new Font("Arial", 0, 14));
    descrPanel.add(label);
    descrPanel.setBackground(Color.WHITE);

    getContentPane().add(descrPanel, BorderLayout.NORTH);

    setJMenuBar(createMenuBar());

    setBounds(20, 20, 1000, 720);
    setVisible(true);

    // unnecessary leads to a double rendering
    // panel.createImage();
  }

  private JMenuBar createMenuBar() {
    JMenuBar menuBar = new JMenuBar();
    JMenu menu = new JMenu("File");

    JCheckBoxMenuItem hideExcluded = new JCheckBoxMenuItem();
    hideExcluded.setMnemonic(KeyEvent.VK_E);
    hideExcluded.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        panel.hideExcluded = ((JCheckBoxMenuItem) e.getSource()).isSelected();
        panel.paintAgain();
      }
    });
    hideExcluded.setText("Hide Excluded");
    menu.add(hideExcluded);

    menuBar.add(menu);

    return menuBar;
  }

  @Override
  public void actionPerformed(ActionEvent ae) {
    String command = ae.getActionCommand();

    System.err.println("Error - unknown command '" + command + "'");
  }

  public static void loadMosaicismResults(Project proj) {
    BufferedReader reader;
    String[] line;
    Vector<String[]> samples;
    Vector<double[]> datapoints;
    SampleData sampleData;
    String[] classes;

    String mosaicFile = proj.MOSAIC_RESULTS_FILENAME.getValue();

    if (!Files.exists(mosaicFile, proj.JAR_STATUS.getValue())) {
      JOptionPane.showMessageDialog(null, "Could not find file: " + mosaicFile, "Error",
                                    JOptionPane.ERROR_MESSAGE);
      return;
    }

    sampleData = proj.getSampleData(2, false);
    classes = sampleData.getClasses();
    if (ext.indexOfStr("mask", classes) >= 0) {
      // TODO left incomplete, what was the goal of this added code??
    }

    samples = new Vector<String[]>();
    datapoints = new Vector<double[]>();
    try {
      reader = Files.getReader(mosaicFile, proj.JAR_STATUS.getValue(), true, true);
      if (!ext.checkHeader(reader.readLine().trim().split("\t"), MOSAICISM_HEADER, false)) {
        proj.message("Different Mosaicism header than expected (see log); this could blow up");
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        if (!line[5].equals(".") && !line[6].equals(".")
            && Integer.parseInt(line[1].substring(3, line[1].length() - 1)) < 23) {
          samples.add(new String[] {line[0], line[1]});
          datapoints.add(new double[] {Double.parseDouble(line[5]), Double.parseDouble(line[6])});
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + mosaicFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + mosaicFile + "\"");
      System.exit(2);
    }

    new MosaicPlot(proj, Matrix.toStringArrays(samples), Matrix.toDoubleArrays(datapoints));
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;

    String usage = "\n" + "plot.MosaicPlot requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
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
