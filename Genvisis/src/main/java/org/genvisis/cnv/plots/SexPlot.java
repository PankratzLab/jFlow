// filter based on chromosome
package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Vector;
import javax.swing.AbstractAction;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.ColorIcon;
import org.genvisis.cnv.gui.WrapLayout;
import org.genvisis.cnv.qc.SexChecks;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.gui.UITools;
import com.google.common.primitives.Booleans;
import com.google.common.primitives.Bytes;

public class SexPlot extends JFrame {

  public static final long serialVersionUID = 1L;

  private static final String[] SEX_CHECKS_REQUIREMENTS = {"Sample", "Sex",
                                                           SexChecks.EST_SEX_HEADER, "Check",
                                                           "Median X LRR", "Median Y LRR",
                                                           "Excluded", "Check", "Note"};

  SexPanel sexPanel;

  public SexPlot(Project proj, String[] samples, double[][] data, byte[] sexes,
                 byte[] estimatedSexes, boolean[] excluded, boolean[] uncertains, String[] notes) {
    super("Genvisis - Sex Plot - " + proj.PROJECT_NAME.getValue());
    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

    sexPanel = new SexPanel(proj, samples, data, sexes, estimatedSexes, excluded, uncertains,
                            notes);
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
    // repaint();

    setMinimumSize(new Dimension(20, 20));
    UITools.setSize(this, new Dimension(1000, 720));
    pack();
    setVisible(true);
    // unnecessary leads to a double rendering
    // sexPanel.createImage();
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

    JMenu actionMenu = new JMenu("Actions");
    actionMenu.setMnemonic(KeyEvent.VK_A);

    final JMenuItem uncertainsToTrailer = new JMenuItem();

    AbstractAction uncertainsToTrailerAction = new AbstractAction() {

      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        sexPanel.uncertainsToTrailer();
      }
    };
    uncertainsToTrailer.setAction(uncertainsToTrailerAction);
    uncertainsToTrailer.setText("Uncertains To Trailer Plot");
    uncertainsToTrailer.setMnemonic(KeyEvent.VK_U);
    actionMenu.add(uncertainsToTrailer);

    menuBar.add(actionMenu);

    return menuBar;

  }

  private JPanel colorLegendPanel() {
    JPanel colorLegendPanel = new JPanel(new WrapLayout(FlowLayout.CENTER, 0, 0));
    colorLegendPanel.setBackground(Color.WHITE);// BACKGROUND_COLOR);

    JLabel legend = new JLabel("Color Key: ");
    legend.setFont(new Font("Arial", 0, 14));
    colorLegendPanel.add(legend);

    for (int i = 0; i < SexPanel.COLOR_SCHEME.length; i++) {
      JPanel enclosure = new JPanel();
      enclosure.setBackground(Color.WHITE);
      enclosure.add(new JLabel(new ColorIcon(12, 12, SexPanel.COLOR_SCHEME[i])));
      enclosure.add(new JLabel(SexPanel.COLOR_SCHEME_MEANING[i]));
      colorLegendPanel.add(enclosure);
    }
    return colorLegendPanel;
  }

  public static void loadSexCheckResults(Project proj) {
    Vector<String> samples = new Vector<>();
    Vector<double[]> datapoints = new Vector<>();
    Vector<Byte> sexes = new Vector<>();
    Vector<Byte> estimatedSexes = new Vector<>();
    Vector<Boolean> excluded = new Vector<>();
    Vector<Boolean> uncertains = new Vector<>();
    Vector<String> notes = new Vector<>();
    try {
      BufferedReader reader = Files.getReader(proj.SEXCHECK_RESULTS_FILENAME.getValue(), true,
                                              false);
      if (reader == null) {
        return;
      }
      String[] line = reader.readLine().trim().split("\t");
      if (!ext.checkHeader(line, SexChecks.SEX_HEADER, false, proj.getLog(), false)) {
        proj.message("The header in file '" + proj.SEXCHECK_RESULTS_FILENAME.getValue()
                     + "' is not as expected and may cause problems; see log for more detail");
      }
      int[] indices = ext.indexFactors(SEX_CHECKS_REQUIREMENTS, line, false);
      for (int index : indices) {
        if (index == -1) {
          return;
        }
      }

      String temp;
      while ((temp = reader.readLine()) != null) {
        line = temp.trim().split("\t", -1);
        if (ext.isMissingValue(line[indices[4]]) || ext.isMissingValue(line[indices[5]])) {
          System.err.println("Error - sample '" + line[indices[0]]
                             + "' does not have a valid medianLRR for X or Y");
        } else {
          samples.add(line[indices[0]]);
          datapoints.add(new double[] {Double.parseDouble(line[indices[4]]),
                                       Double.parseDouble(line[indices[5]])});
          sexes.add(Byte.parseByte(line[indices[1]]));
          estimatedSexes.add(Byte.parseByte(line[indices[2]]));
          excluded.add(line[indices[6]].equals("1"));
          uncertains.add(line[indices[7]].equals("1"));
          notes.add(line[indices[8]]);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + proj.SEXCHECK_RESULTS_FILENAME.getValue()
                         + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + proj.SEXCHECK_RESULTS_FILENAME.getValue()
                         + "\"");
      return;
    }

    new SexPlot(proj, samples.toArray(new String[samples.size()]),
                Matrix.toDoubleArrays(datapoints), Bytes.toArray(sexes),
                Bytes.toArray(estimatedSexes), Booleans.toArray(excluded),
                Booleans.toArray(uncertains), ArrayUtils.toStringArray(notes));
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;

    String usage = "\n" + "cnv.plot.SexPlot requires 0-1 arguments\n"
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
      loadSexCheckResults(new Project(filename));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
