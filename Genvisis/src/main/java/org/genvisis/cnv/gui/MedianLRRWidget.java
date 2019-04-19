package org.genvisis.cnv.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.Collections;
import java.util.concurrent.ExecutionException;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToggleButton;

import org.genvisis.cnv.analysis.MedianLRRWorker;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.PropertyEditorButton;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.cnv.prop.FileProperty;
import org.genvisis.cnv.prop.PropertyKeys;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.gui.UITools;

public class MedianLRRWidget extends JFrame implements Runnable {

  private static final long serialVersionUID = 1L;
  // once a job has been started, used to track completion
  private volatile int computeComplete = 42;
  public static final String FILENAME = "Enter Analysis Name";

  public static final String[] CLASSES_TO_DUMP = {"IID"};
  private static final int PREFERRED_WIDTH = 500;

  private MedianLRRWorker medianLRRWorker;

  // Median LRR worker params

  private final String initRegion;
  private final Project proj;

  public MedianLRRWidget(Project proj, String initRegion) {
    this.proj = proj;
    this.initRegion = initRegion;

  }

  @Override
  public void run() {
    createAndShowGUI();
  }

  public void createAndShowGUI() {
    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    setTitle("Median Log R Ratio Settings");
    WindowListener exitListener = new WindowAdapter() {

      // check for running job before exit
      @Override
      public void windowClosing(WindowEvent e) {
        String btns[] = {"Exit Anyway", "Cancel"};
        if (computeComplete != 42 && !medianLRRWorker.isDone()) {
          int promptResult = JOptionPane.showOptionDialog(null,
                                                          "A Thread is computing Median Log R Ratios\nThread will continue if you exit",
                                                          "Warning - Running thread",
                                                          JOptionPane.DEFAULT_OPTION,
                                                          JOptionPane.WARNING_MESSAGE, null, btns,
                                                          btns[1]);
          System.out.println(promptResult);
          if (promptResult == 0) {
            dispose();
          }
        } else {
          dispose();
        }
      }
    };
    addWindowListener(exitListener);
    MainPanel mainPanel = new MainPanel();
    add(mainPanel);

    pack();
    // fix the width and height expand the height to allow the progress bar
    UITools.setSize(this, new Dimension(PREFERRED_WIDTH, getHeight() + 40));
    pack();
    UITools.centerComponent(this);
    setVisible(true);
  }

  private class MainPanel extends JPanel {

    private static final long serialVersionUID = 2L;
    private static final int CORRECTION_INDEX = 0;
    private static final int TRANSFORMATION_INDEX = 1;
    private static final int FONT_SIZE = 14;

    private JTabbedPane tabbedPane;
    private JCheckBox recomputeLRR = new JCheckBox("Recompute LRR");
    private JCheckBox correctLRR = new JCheckBox("Correct LRR");
    private JCheckBox correctXY = new JCheckBox("Correct XY");
    private ButtonGroup chrXStrategy = new ButtonGroup();

    private JRadioButton rawValuesButton = null;
    private ButtonGroup transformations = new ButtonGroup();
    private ButtonGroup transformScopes = new ButtonGroup();

    private JCheckBox homozygousCheckBox = new JCheckBox("Process homozygous markers only");
    private final JTextField regionTextField;
    private final JTextArea fileInputArea;
    private final JProgressBar progressBar;

    private MainPanel() {
      setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
      tabbedPane = new JTabbedPane();
      add(tabbedPane);

      JPanel correctionPane = new JPanel();
      JPanel transformPane = new JPanel();
      correctionPane.setLayout(new BoxLayout(correctionPane, BoxLayout.Y_AXIS));
      transformPane.setLayout(new BoxLayout(transformPane, BoxLayout.Y_AXIS));
      tabbedPane.add("Correct values", correctionPane);
      tabbedPane.add("Transform values", transformPane);

      progressBar = new JProgressBar(0, 100);

      // Build Correction panel
      for (CHROMOSOME_X_STRATEGY x : CHROMOSOME_X_STRATEGY.values()) {
        JRadioButton btn = makeButton(x.name());
        btn.setToolTipText(x.getToolTip());
        btn.setEnabled(false);
        chrXStrategy.add(btn);
        if (x.equals(CHROMOSOME_X_STRATEGY.ARTIFICIAL)) {
          btn.setSelected(true);
        }
      }

      createGroup(chrXStrategy, true, recomputeLRR, correctLRR, correctXY);
      requiresProperty(proj.INTENSITY_PC_FILENAME, correctLRR, correctXY);

      addLabel(correctionPane, "Recompute LRRs?");
      correctionPane.add(recomputeLRR);
      addLabel(correctionPane, "Correct LRRs?");
      correctionPane.add(correctLRR);
      correctionPane.add(correctXY);
      addLabel(correctionPane, "Sex-specific correction strategy for chrX (if present)");
      for (AbstractButton b : Collections.list(chrXStrategy.getElements())) {
        correctionPane.add(b);
      }

      homozygousCheckBox.setSelected(false);
      homozygousCheckBox.setToolTipText("Only Homozygous calls will be included in the median value computation (does not change other metrics)");
      correctionPane.add(homozygousCheckBox);

      // Build Transformation panel
      for (String transform : Transforms.TRANFORMATIONS) {
        JRadioButton btn = makeButton(transform);
        if ("Raw Values".equals(transform)) {
          rawValuesButton = btn;
          btn.setSelected(true);
        }
        transformations.add(btn);
      }

      ActionListener scopesListener = createGroup(transformScopes, false, rawValuesButton);
      addLabel(transformPane, "Select a transformation method");
      for (AbstractButton b : Collections.list(transformations.getElements())) {
        b.addActionListener(scopesListener);
        transformPane.add(b);
      }
      addLabel(transformPane, "Transform across genome, or chromosome?");
      for (String scope : Transforms.SCOPES) {
        JRadioButton btn = makeButton(scope);
        if (Transforms.SCOPES[0].equals(scope)) {
          btn.setSelected(true);
        }
        btn.setEnabled(false);
        transformPane.add(btn);
        transformScopes.add(btn);
      }

      // Build the common panel area
      // FIXME can't get these last labels to left-justify..
      JPanel bottomPane = new JPanel();
      bottomPane.setLayout(new BoxLayout(bottomPane, BoxLayout.Y_AXIS));
      addLabel(bottomPane, "Output File Prefix");
      fileInputArea = new JTextArea(initRegion);
      fileInputArea.setBorder(BorderFactory.createLineBorder(Color.BLACK));
      fileInputArea.setMaximumSize(new Dimension(Integer.MAX_VALUE, getMinimumSize().height));
      bottomPane.add(fileInputArea);

      addLabel(bottomPane, "Input Regions of Interest");
      regionTextField = new JTextField(initRegion);
      regionTextField.setBorder(BorderFactory.createLineBorder(Color.BLACK));
      regionTextField.setMaximumSize(regionTextField.getPreferredSize());
      regionTextField.setToolTipText("Enter UCSC region or probeset_id;(Your Region Name);marker names(;)");
      JScrollPane scroll = new JScrollPane(regionTextField);
      bottomPane.add(scroll);

      JButton computeButton = new JButton("Compute");
      computeButton.setToolTipText("Compute Median Log R Ratios");
      computeButton.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          validateAndRun();
        }
      });

      JButton plotButton = new JButton("2D Plot");
      plotButton.setToolTipText("Launch 2D plot to visualize results");
      plotButton.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          open2DPlot();
        }
      });

      JPanel buttons = new JPanel();
      buttons.add(computeButton);
      buttons.add(plotButton);
      // TODO should really re-validate buttons after changing properties
      buttons.add(new PropertyEditorButton(proj, PropertyKeys.KEY_INTENSITY_PC_NUM_COMPONENTS,
                                           PropertyKeys.KEY_INTENSITY_PC_FILENAME));
      buttons.setAlignmentX(CENTER_ALIGNMENT);
      bottomPane.add(buttons);

      add(bottomPane);
      pack();
    }

    /**
     * Open the last MedianLRRWorker's output in {@link TwoDPlot}
     */
    private void open2DPlot() {
      if (computeComplete == 42) {
        JOptionPane.showMessageDialog(this, "Please compute Median values before visualizing");
      } else {
        String fileNameToVisualize;
        revalidate();
        try {
          fileNameToVisualize = medianLRRWorker.get();
          TwoDPlot twoDplot = TwoDPlot.createGUI(proj, true, true);
          twoDplot.showSpecificFile(fileNameToVisualize, // 2, 3);
                                    MedianLRRWorker.CLASSES_TO_DUMP.length + 1,
                                    MedianLRRWorker.CLASSES_TO_DUMP.length + 2);
          twoDplot.updateGUI();

        } catch (InterruptedException e) {
          JOptionPane.showMessageDialog(this,
                                        "Thread was interupted when computing Median Log R Ratios");
          e.printStackTrace();
        } catch (ExecutionException e) {
          JOptionPane.showMessageDialog(this, "There was an error computing Median Log R Ratios");
          e.printStackTrace();
        }
      }
    }

    /**
     * Validate the GUI settings and if OK, start a new MedianLRRWorker
     */
    private void validateAndRun() {
      if (medianLRRWorker != null && !medianLRRWorker.isDone()) {
        JOptionPane.showMessageDialog(this, "Thread is busy computing median Log R Ratios");
        return;
      }

      int transformType = 0;
      int transformScope = 0;
      CHROMOSOME_X_STRATEGY xStrategy = CHROMOSOME_X_STRATEGY.ARTIFICIAL;
      boolean doRecompute = false;
      boolean doLRR = false;
      boolean doXY = false;

      // Update the output file name based on whether this is a correciton or transformation
      String customName = fileInputArea.getText();
      if (customName == null) {
        customName = "";
      }

      if (tabbedPane.getSelectedIndex() == CORRECTION_INDEX) {
        if (recomputeLRR.isSelected()) {
          doRecompute = true;
          customName += "_RECOMPUTE_LRR";
        }

        if (correctLRR.isSelected()) {
          customName += "_CORRECT_LRR";
          doLRR = true;
        }

        if (correctXY.isSelected()) {
          customName += "_CORRECT_XY";
          doXY = true;
        }

        if (doRecompute || doLRR || doXY) {
          xStrategy = CHROMOSOME_X_STRATEGY.valueOf(chrXStrategy.getSelection().getActionCommand());
          customName += "_" + chrXStrategy.getSelection().getActionCommand();
        }
      } else if (tabbedPane.getSelectedIndex() == TRANSFORMATION_INDEX) {
        String transformTypeName = transformations.getSelection().getActionCommand();
        String transformScopeName = transformScopes.getSelection().getActionCommand();
        customName += "_" + transformTypeName;
        customName += "_" + transformScopeName;
        transformType = ext.indexOfStr(transformTypeName, Transforms.TRANFORMATIONS);
        transformScope = ext.indexOfStr(transformScopeName, Transforms.SCOPES);
      }

      customName = ext.replaceWithLinuxSafeCharacters(customName, true);

      // Check if output file would overwrite
      boolean valid;

      if (customName.isEmpty()) {
        valid = false;
        JOptionPane.showMessageDialog(null,
                                      "Please specify an output name or select correction/transformation options");
      } else if (!MedianLRRWorker.checkExists(proj, customName)) {
        valid = true;
      } else {
        valid = false;
        String overwrite[] = {"Overwrite", "Cancel"};
        int promptResult = JOptionPane.showOptionDialog(this,
                                                        "The Files for Analysis " + customName
                                                              + " Exist",
                                                        "Warning - Analysis Files Exist",
                                                        JOptionPane.DEFAULT_OPTION,
                                                        JOptionPane.WARNING_MESSAGE, null,
                                                        overwrite, overwrite[1]);
        if (promptResult == 0) {
          valid = true;
        }
      }
      if (valid) {
        startJob(customName, transformType, transformScope, doRecompute, doLRR, doXY, xStrategy);
      }
    }

    /**
     * Build and run the {@link MedianLRRWorker}
     */
    private void startJob(String outputName, int transformType, int transformScope,
                          boolean recomputeLRR, boolean correctLRR, boolean correctXY,
                          CHROMOSOME_X_STRATEGY xStrategy) {
      add(progressBar);
      progressBar.setVisible(true);
      progressBar.setStringPainted(true);
      computeComplete = 0;

      medianLRRWorker = new MedianLRRWorker(proj, regionTextField.getText().split("\n"),
                                            transformType, transformScope, outputName, progressBar,
                                            recomputeLRR, correctLRR, correctXY,
                                            homozygousCheckBox.isSelected(), xStrategy,
                                            proj.getLog());
      medianLRRWorker.execute();

      revalidate();
    }

    /**
     * Helper method to create and add a JLabel with particular formatting to a panel
     */
    private JLabel addLabel(JPanel pane, String text) {
      Font f = new Font("Arial", 0, FONT_SIZE);
      JLabel label = new JLabel(text);
      Dimension d = new Dimension(PREFERRED_WIDTH - 10, FONT_SIZE + 2);
      label.setPreferredSize(d);
      label.setFont(f);
      pane.add(label);
      return label;
    }

    private JRadioButton makeButton(String text) {
      JRadioButton b = new JRadioButton(text);
      b.setActionCommand(text);
      return b;
    }

    /**
     * @return A {@link DisableGroup} action listener constructed with the given parameters
     */
    private ActionListener createGroup(ButtonGroup toggleIf, boolean checked,
                                       JToggleButton... boxes) {
      if (boxes == null || boxes.length == 0) {
        return null;
      }

      return new DisableGroup(toggleIf, checked, boxes);
    }

    /**
     * @return A {@link RequiresFileProperty} action listener constructed with the given parameters
     */
    private ActionListener requiresProperty(FileProperty property, JToggleButton... toggles) {
      if (toggles == null || toggles.length == 0) {
        return null;
      }

      return new RequiresFileProperty(property, toggles);
    }
  }

  /**
   * ActionListener that listens to an arbitrary number of JToggleButtons and keeps them disabled if
   * the given project property is unset/invalid
   */
  private static class RequiresFileProperty extends ToggleListener {

    private final FileProperty prop;

    public RequiresFileProperty(FileProperty property, JToggleButton... toggles) {
      super(toggles);

      this.prop = property;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
      if (!Files.exists(prop.getValue())) {
        JOptionPane.showMessageDialog(null,
                                      "This option requires a valid setting for the project property: "
                                            + prop.getName());
        for (JToggleButton btn : toggles()) {
          btn.setSelected(false);
        }
      }
    }
  }

  /**
   * ActionListener that listens to an arbitrary number of JToggleButtons and disables a
   * corresponding ButtonGroup if none of them are in the specified state
   */
  private static class DisableGroup extends ToggleListener {

    private final ButtonGroup group;
    private final boolean state;

    public DisableGroup(ButtonGroup toggleIfChecked, boolean state, JToggleButton... toggles) {
      super(toggles);
      this.group = toggleIfChecked;
      this.state = state;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
      boolean disable = false;
      for (int i = 0; !disable && i < toggles().length; i++) {
        disable = toggles()[i].isSelected() == state;
      }

      for (AbstractButton b : Collections.list(group.getElements())) {
        b.setEnabled(disable);
      }
    }
  }

  /**
   * Abstract superclass for action listeners that list to a number of JToggleButtons
   */
  private abstract static class ToggleListener implements ActionListener {

    private final JToggleButton[] toggles;

    public ToggleListener(JToggleButton... toggles) {
      this.toggles = toggles;
      for (JToggleButton toggle : toggles) {
        toggle.addActionListener(this);
      }
    }

    public JToggleButton[] toggles() {
      return toggles;
    }
  }
}
