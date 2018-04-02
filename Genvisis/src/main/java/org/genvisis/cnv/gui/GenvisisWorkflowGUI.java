package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedSet;
import java.util.Vector;
import java.util.concurrent.ConcurrentMap;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import net.miginfocom.swing.MigLayout;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.AndRequirementSet;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Step.FINAL_CODE;
import org.genvisis.cnv.workflow.StepTask;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.ext;
import org.genvisis.common.gui.JProgressBarListener;
import org.genvisis.common.gui.Task;
import org.genvisis.qsub.Qsub;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class GenvisisWorkflowGUI extends JDialog {

  private static final long serialVersionUID = 1L;

  private final JPanel contentPanel = new JPanel();

  private final ConcurrentMap<Step, JCheckBox> checkBoxes = Maps.newConcurrentMap();
  private final ConcurrentMap<Step, JLabel> descLabels = Maps.newConcurrentMap();
  private final ConcurrentMap<Step, Map<Requirement, JLabel>> requirementsLabels = Maps.newConcurrentMap();
  private final ConcurrentMap<Step, JAccordionPanel> panels = Maps.newConcurrentMap();
  public ConcurrentMap<Step, Map<Requirement, JComponent>> varFields = Maps.newConcurrentMap();
  public ConcurrentMap<Step, JProgressBar> progBars = Maps.newConcurrentMap();
  public ConcurrentMap<Step, JProgressBarListener> progListeners = Maps.newConcurrentMap();
  public ConcurrentMap<Step, Task<Void, Void>> progTasks = Maps.newConcurrentMap();
  public ConcurrentMap<Step, ArrayList<JButton>> fileBtns = Maps.newConcurrentMap();
  public ConcurrentMap<Step, JLabel> alreadyRunLbls = Maps.newConcurrentMap();
  public ConcurrentMap<Step, JButton> cancelStepBtns = Maps.newConcurrentMap();

  Project proj;

  private static final String TOP_LABEL = "Genvisis Project Workflow:";

  volatile boolean cancelled = false;
  volatile Set<Step> selected;
  SortedSet<Step> steps;
  int DEFAULT_SCROLL_SPEED = 16;

  /**
   * Create the dialog.
   * 
   * @param steps TODO
   */
  public GenvisisWorkflowGUI(Project proj2, final Launch launch, final SortedSet<Step> steps) {
    proj = proj2;
    if (proj == null) {
      doClose();
      return;
    } else {
      launch.loadProjects();
      launch.setIndexOfCurrentProject(ext.removeDirectoryInfo(proj.getPropertyFilename()));
      proj = launch.loadProject();
    }
    proj.getLog().report("Launching Genvisis Project Pipeline");
    this.steps = steps;
    selected = Sets.newTreeSet(steps);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    JPanel optionPanel = new JPanel();
    JScrollPane scrollPane = new JScrollPane(optionPanel);
    scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
    scrollPane.getVerticalScrollBar().setUnitIncrement(DEFAULT_SCROLL_SPEED);
    contentPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
    contentPanel.add(scrollPane, "cell 0 2, grow");
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    optionPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
    {
      JLabel lblKitAndKaboodle = new JLabel(TOP_LABEL);
      lblKitAndKaboodle.setFont(new Font("Arial", Font.BOLD, 16));
      contentPanel.add(lblKitAndKaboodle, "cell 0 0,alignx center");
    }
    {
      int i = 1;
      for (Step step : steps) {
        final int index = i++;
        JAccordionPanel panel = createPanel(step, index);
        optionPanel.add(panel, "cell 0 " + (i) + ", alignx left, growx");
        panels.put(step, panel);
      }
    }
    Insets btnInsets = new Insets(0, 1, 0, 1);
    {
      JPanel buttonPane = new JPanel();
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      buttonPane.setLayout(new MigLayout("hidemode 2", "[][][][][][][][grow][][][][]", "[]"));

      JLabel lblSelect = new JLabel("Select:");
      buttonPane.add(lblSelect, "flowx,cell 0 0");

      btnSelectAll = new JButton("All");
      btnSelectAll.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent arg0) {
          for (Entry<Step, JCheckBox> entry : checkBoxes.entrySet()) {
            entry.getValue().setSelected(true);
          }
          selected.addAll(steps);
          refreshLabels(GenvisisWorkflowGUI.this, steps);
        }
      });
      btnSelectAll.setMargin(btnInsets);
      buttonPane.add(btnSelectAll, "cell 0 0");

      btnDeselectAll = new JButton("None");
      btnDeselectAll.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          for (Entry<Step, JCheckBox> entry : checkBoxes.entrySet()) {
            entry.getValue().setSelected(false);
          }
          selected.clear();
          refreshLabels(GenvisisWorkflowGUI.this, steps);
        }
      });
      btnDeselectAll.setMargin(btnInsets);
      buttonPane.add(btnDeselectAll, "cell 1 0");

      btnSelectValid = new JButton("Valid");
      btnSelectValid.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          Set<Step> selectedSteps = Sets.newHashSet();
          for (Entry<Step, JCheckBox> entry : checkBoxes.entrySet()) {
            selectedSteps.add(entry.getKey()); // pretend everything is selected
          }
          Map<Step, Map<Requirement, String>> variables = getVariables();
          for (final Step step : steps) {
            if (step == null || checkBoxes.get(step) == null || varFields.get(step) == null) {
              continue;
            }
            if (!step.checkIfOutputExists(variables)) {
              boolean check = step.hasRequirements(proj, selectedSteps, variables);
              checkBoxes.get(step).setSelected(check);
              selected.add(step);
              if (check) {
                selectedSteps.add(step);
              } else {
                selected.remove(step);
                selectedSteps.remove(step);
              }
            } else {
              selectedSteps.remove(step);
              selected.remove(step);
              checkBoxes.get(step).setSelected(false);
            }
          }
          refreshLabels(GenvisisWorkflowGUI.this, steps);
        }
      });
      btnSelectValid.setMargin(btnInsets);
      buttonPane.add(btnSelectValid, "cell 2 0");

      JSeparator separator = new JSeparator();
      separator.setOrientation(SwingConstants.VERTICAL);
      buttonPane.add(separator, "cell 3 0,growy");

      JLabel lblCollapse = new JLabel("Collapse:");
      buttonPane.add(lblCollapse, "cell 4 0");

      JButton btnAll = new JButton("All");
      btnAll.setMargin(btnInsets);
      btnAll.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          for (Entry<Step, JAccordionPanel> entry : panels.entrySet()) {
            entry.getValue().shrink();
          }
        }
      });
      buttonPane.add(btnAll, "cell 5 0");

      JButton btnNone = new JButton("None");
      btnNone.setMargin(btnInsets);
      btnNone.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          for (Entry<Step, JAccordionPanel> entry : panels.entrySet()) {
            entry.getValue().expand();
          }
        }
      });
      buttonPane.add(btnNone, "cell 6 0");

      progVal = new JProgressBar(0, steps.size());
      progVal.setValue(0);
      progVal.setStringPainted(true);
      progVal.setString("Validating Steps...");
      progVal.setVisible(false);
      buttonPane.add(progVal, "growx, cell 7 0");

      AbstractAction listener = new AbstractAction() {

        private static final long serialVersionUID = 1L;

        @Override
        public void actionPerformed(ActionEvent e) {
          if (running) {
            return;
          }
          if (e.getActionCommand().equals("Close")) {
            cancelled = true;
            doClose();
          } else if (e.getActionCommand().equals("Export")) {
            exportToText(false);
          } else if (e.getActionCommand().equals("ExportDefaults")) {
            exportToText(true);
          } else {
            run();
          }
        }
      };

      btnExportToText = new JButton("Export To Text");
      btnExportToText.setActionCommand("Export");
      btnExportToText.addActionListener(listener);
      btnExportToText.setMargin(btnInsets);

      btnExportToTextDefaults = new JButton("Export Defaults");
      btnExportToTextDefaults.setActionCommand("ExportDefaults");
      btnExportToTextDefaults.addActionListener(listener);
      btnExportToTextDefaults.setMargin(btnInsets);
      if (Files.programExists("qsub")) {
        buttonPane.add(btnExportToText, "cell 8 0, alignx right");
        buttonPane.add(btnExportToTextDefaults, "cell 9 0, alignx right");
      } else {
        buttonPane.add(btnExportToText, "cell 9 0, alignx right");
      }
      btnOk = new JButton("Run");
      btnOk.setMargin(btnInsets);
      btnOk.setActionCommand("Run");
      btnOk.setMnemonic(KeyEvent.VK_O);
      buttonPane.add(btnOk, "cell 10 0,alignx right");
      getRootPane().setDefaultButton(btnOk);
      btnCancel = new JButton("Close");
      btnCancel.setMargin(btnInsets);
      btnCancel.setActionCommand("Close");
      btnCancel.setMnemonic(KeyEvent.VK_C);
      buttonPane.add(btnCancel, "cell 11 0,alignx left");

      btnOk.addActionListener(listener);
      btnCancel.addActionListener(listener);

    }
    pack();
    setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    addWindowListener(new WindowAdapter() {

      @Override
      public void windowClosing(WindowEvent e) {
        if (running) {
          return;
        }
        doClose();
      }
    });
    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {
        Map<Step, Map<Requirement, String>> variables = getVariables();
        for (Step step : steps) {
          if (step.checkIfOutputExists(variables)) {
            checkBoxes.get(step).setSelected(false);
            alreadyRunLbls.get(step).setVisible(true);
            panels.get(step).shrink();
            selected.remove(step);
          }
        }
      }
    });
    refreshLabels(this, steps);
    setMinimumSize(new Dimension(100, 100));
    UITools.setSize(this, new Dimension(750, 850));
    setTitle(TOP_LABEL);
    addWindowFocusListener(new WindowFocusListener() {

      @Override
      public void windowLostFocus(WindowEvent e) {}

      @Override
      public void windowGainedFocus(WindowEvent e) {
        launch.toFront();
        GenvisisWorkflowGUI.this.toFront();
      }
    });

    InputMap inMap = getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    ActionMap actMap = getRootPane().getActionMap();

    KeyStroke escKey = KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0);
    AbstractAction escapeAction = new AbstractAction() {

      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        doClose();
      }
    };
    inMap.put(escKey, "Action.escape");
    actMap.put("Action.escape", escapeAction);

    pack();
  }

  protected void doClose() {
    // TODO check if running, may not be able to stop current step. Ask if should interrupt. Maybe
    // do cleanup if desired?
    cancelled = true;
    setVisible(false);
    dispose();
  }

  private JAccordionPanel createPanel(final Step step, final int index) {
    final JAccordionPanel panel = new JAccordionPanel();

    final JCheckBox chckbx = new JCheckBox();
    chckbx.setAction(new StepRefresher(GenvisisWorkflowGUI.this, step) {

      /**
      *
      */
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        if (chckbx.isSelected()) {
          selected.add(step);
        } else {
          selected.remove(step);
        }
        super.actionPerformed(e);
      }
    });
    chckbx.setFont(chckbx.getFont().deriveFont(Font.PLAIN, 14));
    Grafik.scaleCheckBoxIcon(chckbx);
    chckbx.setSelected(selected.contains(step));
    chckbx.setToolTipText(step.getDescription());
    chckbx.setText(index + ": " + step.getName());

    checkBoxes.put(step, chckbx);
    panel.topPanel.add(chckbx, "cell 0 0");

    JProgressBarListener progListener = new JProgressBarListener(false, step.getName());
    JProgressBar stepProgBar = progListener.getBar();
    progBars.put(step, stepProgBar);
    progListeners.put(step, progListener);
    panel.topPanel.add(stepProgBar, "cell 1 0, alignx right, hidemode 3, split 2");
    JLabel alreadyRanLbl = new JLabel("Output Already Exists!");
    alreadyRanLbl.setVisible(false);
    alreadyRunLbls.put(step, alreadyRanLbl);
    panel.topPanel.add(alreadyRanLbl, "cell 1 0, alignx right, hidemode 3");
    JButton cancelStepButton = new JButton();
    cancelStepButton.addActionListener(e -> {
      Task<Void, Void> st = progTasks.get(step);
      if (st != null) {
        st.cancel(true);
      }
    });
    cancelStepButton.setIcon(Grafik.getImageIcon("images/redx.png"));
    cancelStepButton.setVisible(false);
    cancelStepButton.setIconTextGap(0);
    cancelStepButton.setMargin(new Insets(1, 1, 1, 1));
    cancelStepBtns.put(step, cancelStepButton);
    panel.topPanel.add(cancelStepButton, "cell 1 0, alignx right, hidemode 3");

    JLabel descLbl = new JLabel("<html><center><p>" + step.getDescription()
                                + "</p></center></html>");
    descLbl.setVerticalAlignment(SwingConstants.TOP);
    descLbl.setHorizontalAlignment(SwingConstants.LEFT);
    descLbl.setFont(descLbl.getFont().deriveFont(Font.PLAIN));
    descLabels.put(step, descLbl);

    panel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
    String rows = "[][]";
    for (int i = 0; i < step.getRequirements().getFlatRequirementsList().size(); i++) {
      rows = rows + "[]";
    }
    panel.contentPanel.setLayout(new MigLayout("", "[200px,grow]push[200px,grow]", rows));
    panel.contentPanel.add(descLbl, "cell 0 0");

    RequirementSet reqs = step.getRequirements();

    if (reqs.size() > 0) {
      JLabel reqLbl = new JLabel("Requires:");
      panel.contentPanel.add(reqLbl, "cell 0 1");

      Map<Requirement, JLabel> reqLbls = new java.util.HashMap<>();
      requirementsLabels.put(step, reqLbls);
      int rowIndex = 2;
      addLabels(step, reqs, panel.contentPanel, rowIndex, new ArrayList<Integer>());

      Map<Requirement, JComponent> reqInputFields = Maps.newLinkedHashMap();
      varFields.put(step, reqInputFields);
      rowIndex = 2;
      addReqFields(step, reqs, panel.contentPanel, rowIndex);
    }

    return panel;

  }

  private Map<Requirement, JComponent> addToPanel(Step step, Requirement req, JPanel panel,
                                                  int rowIndex) {
    Map<Requirement, JComponent> reqInputFields = Maps.newLinkedHashMap();
    if (req.getType() == Requirement.RequirementInputType.BOOL) {
      JCheckBox checkBox = new JCheckBox();
      checkBox.setAction(new StepRefresher(GenvisisWorkflowGUI.this, step));
      checkBox.setFont(checkBox.getFont().deriveFont(14));
      Grafik.scaleCheckBoxIcon(checkBox);
      checkBox.setVerticalAlignment(SwingConstants.TOP);
      checkBox.setHorizontalAlignment(SwingConstants.RIGHT);
      boolean sel = Boolean.parseBoolean(req.getDefaultValue().toString());
      checkBox.setSelected(sel);
      reqInputFields.put(req, checkBox);
      panel.add(checkBox, "alignx right, aligny center, growx, gapleft 20, cell 1 " + rowIndex);
    } else if (req.getType() == Requirement.RequirementInputType.ENUM) {
      Object o = req.getDefaultValue();
      Enum<?>[] vals = ((Enum<?>) o).getClass().getEnumConstants();
      JComboBox<Enum<?>> combo = new JComboBox<>(vals);
      combo.setAction(new StepRefresher(GenvisisWorkflowGUI.this, step));
      combo.setFont(combo.getFont().deriveFont(14));
      combo.setSelectedItem(o);
      reqInputFields.put(req, combo);
      panel.add(combo, "alignx right, aligny center, growx, gapleft 20, cell 1 " + rowIndex);
    } else if (req instanceof Requirement.ListSelectionRequirement) {
      Requirement.ListSelectionRequirement listSelectionReq = (Requirement.ListSelectionRequirement) req;
      JList<String> jList = new JList<>(new Vector<>(listSelectionReq.getOptions()));
      jList.setFont(jList.getFont().deriveFont(14));
      for (String defaultOption : listSelectionReq.getDefaultOptions()) {
        jList.setSelectedValue(defaultOption, false);
      }
      reqInputFields.put(req, jList);
      panel.add(new JScrollPane(jList),
                "alignx right, aligny center, growx, gapleft 20, cell 1 " + rowIndex);
    } else if (req.getType() != Requirement.RequirementInputType.NONE) {
      JTextField textField = new JTextField();
      textField.getDocument().addDocumentListener(new TextChangedListener() {

        @Override
        public void changedUpdate(DocumentEvent e) {
          refreshLabels(GenvisisWorkflowGUI.this, step.getRelatedSteps());
        }
      });
      textField.setText(req.getDefaultValue().toString());
      reqInputFields.put(req, textField);
      panel.add(textField,
                "alignx right, aligny center, growx, gapleft 20, split 1, cell 1 " + rowIndex);
      if (req.getType() == Requirement.RequirementInputType.FILE
          || req.getType() == Requirement.RequirementInputType.DIR) {
        JButton fileBtn = new JButton();
        fileBtn.setAction(new AbstractAction() {

          private static final long serialVersionUID = 1L;

          @Override
          public void actionPerformed(ActionEvent e) {
            JTextField fileField = (JTextField) varFields.get(step).get(req);

            String current = fileField.getText();

            String dir = "".equals(current) ? proj.PROJECT_DIRECTORY.getValue(false, false)
                                            : ext.parseDirectoryOfFile(current);
            JFileChooser chooser = new JFileChooser(dir);
            chooser.setMultiSelectionEnabled(false);
            Requirement.RequirementInputType type = req.getType();

            if (type == Requirement.RequirementInputType.FILE) {
              chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
              chooser.setDialogTitle("Select File");
            } else if (type == Requirement.RequirementInputType.DIR) {
              chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
              chooser.setDialogTitle("Select Directory");
            }

            int retValue = chooser.showDialog(GenvisisWorkflowGUI.this, "Select");

            if (retValue == JFileChooser.CANCEL_OPTION) {
              refreshLabels(GenvisisWorkflowGUI.this, step.getRelatedSteps());
              return;
            } else {
              File newFile = chooser.getSelectedFile();
              String txt = ext.verifyDirFormat(newFile.getAbsolutePath());
              if (type == Requirement.RequirementInputType.FILE) {
                txt = txt.substring(0, txt.length() - 1);
              }
              fileField.setText(txt);
            }
            refreshLabels(GenvisisWorkflowGUI.this, step.getRelatedSteps());
          }
        });
        fileBtn.setText("...");
        fileBtn.setMargin(new Insets(1, 4, 0, 3));
        panel.add(fileBtn, "cell 1 " + rowIndex);
        ArrayList<JButton> list = fileBtns.get(step);
        if (list == null) {
          list = new ArrayList<>();
          fileBtns.put(step, list);
        }
        list.add(fileBtn);
      }
    }
    return reqInputFields;
  }

  /**
   * @param step Step
   * @param rs RequirementSet
   * @param panel Panel to which to add components
   * @param rowIndex row index at which to start adding components
   * @return last used row index
   */
  private int addReqFields(Step step, RequirementSet rs, JPanel panel, int rowIndex) {
    List<Requirement> levelReqs = rs.getRequirements();
    List<RequirementSet> levelReqSets = rs.getRequirementSets();

    for (int i = 0; i < levelReqs.size(); i++) {
      varFields.get(step).putAll(addToPanel(step, levelReqs.get(i), panel, rowIndex));
      if (i < levelReqs.size() - 1 || levelReqSets.size() > 0) {
        rowIndex += 2; // skip past join labels (e.g. "AND" / "OR")
      }
    }
    for (int i = 0; i < levelReqSets.size(); i++) {
      rowIndex = addReqFields(step, levelReqSets.get(i), panel, rowIndex);
      rowIndex++;
      if (i < levelReqSets.size() - 1) {
        rowIndex += 1; // skip past join labels
      }
    }
    return rowIndex;
  }

  /**
   * Add requirement labels (and recursively add sub-labels) to a panel;
   * 
   * @param rs RequirementSet
   * @param panel JPanel
   * @param rowIndex rowIndex on which to start adding labels
   * @param argSetNums Argument index prefixes for sub-requirement indices
   * @return last used row index
   */
  private int addLabels(Step step, RequirementSet rs, JPanel panel, int rowIndex,
                        List<Integer> argSetNums) {
    List<Requirement> levelReqs = rs.getRequirements();
    List<RequirementSet> levelReqSets = rs.getRequirementSets();
    String join = rs.getJoinString();

    String prefStr = ArrayUtils.toStr(argSetNums, ".");
    for (int i = 0; i < levelReqs.size(); i++) {
      JLabel indLbl = new JLabel(prefStr + (levelReqs.size() > 1 ? ("." + (i + 1)) : "") + ": ");
      indLbl.setFont(indLbl.getFont().deriveFont(Font.PLAIN, 9));
      panel.add(indLbl, "gapleft 25, aligny top, split 1, cell 0 " + rowIndex);
      JLabel requirementLbl = new JLabel("<html><p>" + sanitize(levelReqs.get(i).getDescription())
                                         + "</p></html>");
      requirementLbl.setFont(requirementLbl.getFont().deriveFont(Font.PLAIN, 9));
      panel.add(requirementLbl, "aligny top, cell 0 " + rowIndex);
      requirementsLabels.get(step).put(levelReqs.get(i), requirementLbl);
      if (i < levelReqs.size() - 1 || levelReqSets.size() > 0) {
        JLabel joinLbl = new JLabel(join);
        joinLbl.setFont(joinLbl.getFont().deriveFont(Font.PLAIN, 10));
        panel.add(joinLbl, "gapleft 18, cell 0 " + rowIndex);
        rowIndex++;
      }
    }
    // for each subset of requirements
    for (int i = 0; i < levelReqSets.size(); i++) {
      Integer add = Integer.valueOf(i + 1);
      // track argument index to which we're adding sub-requirements (i.e. {1, 1} becomes "1.1.N"
      // for each of N sub-requirements)
      argSetNums.add(add);
      // recursively add sub-requirements, recording the rowIndex of the last subitem added to the
      // panel
      rowIndex = addLabels(step, levelReqSets.get(i), panel, rowIndex, argSetNums);
      // remove last argument index but keep higher ones (i.e. {1, 1} becomes {1} and next iteration
      // becomes {1, 2}
      argSetNums.remove(argSetNums.size() - 1);
      rowIndex++;
      if (i < levelReqSets.size() - 1) {
        // if we have more requirements, add a label for the logical operator of this RequirementSet
        JLabel joinLbl = new JLabel(join);
        joinLbl.setFont(joinLbl.getFont().deriveFont(Font.PLAIN, 10));
        panel.add(joinLbl, "gapleft 18, cell 0 " + rowIndex + " 2 1");

        if (rs instanceof AndRequirementSet) {
          JSeparator jSep = new JSeparator(JSeparator.HORIZONTAL);
          jSep.setForeground(jSep.getForeground().brighter());
          panel.add(jSep, "cell 0 " + rowIndex + " 2 1, pad 1 0 0 0, growx");
        }

        rowIndex++;
      }
    }
    return rowIndex;
  }

  private String sanitize(String str) {
    return str.replaceAll("<", "&lt;");
  }

  public Set<Step> getSelectedOptions() {
    return selected;
  }

  public boolean getCancelled() {
    return cancelled;
  }

  public void startStep(Step step) {
    if (alreadyRunLbls.get(step).isVisible()) {
      alreadyRunLbls.get(step).setVisible(false);
    }
    progBars.get(step).setString("Running...");
    progBars.get(step).setStringPainted(true);
    cancelStepBtns.get(step).setVisible(true);
    boolean enableCancel = !step.getFlags()
                                .contains(org.genvisis.cnv.workflow.Requirement.Flag.MULTITHREADED);
    cancelStepBtns.get(step).setEnabled(enableCancel);
    String cancel1 = "Cancel";
    // TODO change when able to cancel:
    String cancel2 = "Unable to cancel multi-threaded operation";
    cancelStepBtns.get(step).setToolTipText(enableCancel ? cancel1 : cancel2);
  }

  public void endStep(Step step, FINAL_CODE code) {
    String resp = code.getMessage();
    progBars.get(step).setString(resp);
    cancelStepBtns.get(step).setVisible(false);
    cancelStepBtns.get(step).setEnabled(false);
  }

  /**
   * Helper {@link DocumentListener} that redirects both {@link #insertUpdate(DocumentEvent)} and
   * {@link #removeUpdate(DocumentEvent)} to {@link #changedUpdate(DocumentEvent)}.
   */
  public abstract static class TextChangedListener implements DocumentListener {

    @Override
    public void insertUpdate(DocumentEvent e) {
      changedUpdate(e);
    }

    @Override
    public void removeUpdate(DocumentEvent e) {
      changedUpdate(e);
    }
  }
  /**
   * Helper {@link ActionListener} to refresh one or more UI steps.
   */
  public static class StepRefresher extends AbstractAction {

    /**
    *
    */
    private static final long serialVersionUID = 1L;
    private final transient Set<Step> stepsToRefresh;
    private final GenvisisWorkflowGUI refrenceGUI;

    /**
     * @param gui UI which is displaying the given STEPs.
     * @param steps STEPs to validate by this {@link ActionListener}.
     */
    public StepRefresher(final GenvisisWorkflowGUI gui, final Step... steps) {
      refrenceGUI = gui;
      stepsToRefresh = new HashSet<>();
      for (final Step s : steps) {
        stepsToRefresh.addAll(s.getRelatedSteps());
      }

    }

    @Override
    public void actionPerformed(ActionEvent e) {
      refreshLabels(refrenceGUI, stepsToRefresh);
    }
  }

  private Set<Step> getAllRelatedSteps(final Collection<Step> refreshSteps) {
    HashSet<Step> allSteps = new HashSet<>();
    allSteps.addAll(refreshSteps);
    boolean go = true;
    while (go) {
      int stepCnt = allSteps.size();
      for (Step s : this.steps) {
        if (!Collections.disjoint(refreshSteps, s.getRelatedSteps())) {
          allSteps.add(s);
        }
      }
      if (allSteps.size() == stepCnt) {
        go = false;
      }
    }
    return Collections.unmodifiableSet(allSteps);
  }

  /**
   * Validate all elements of the given {@link Step}s and refresh the specified UI.
   *
   * @param stepsToRefresh
   */
  public static void refreshLabels(final GenvisisWorkflowGUI gui, Collection<Step> steps) {
    final Collection<Step> stepsToRefresh = gui.getAllRelatedSteps(steps);
    new Thread(new Runnable() {

      @Override
      public void run() {
        try {
          SwingUtilities.invokeAndWait(new Runnable() {

            @Override
            public void run() {
              gui.progVal.setValue(0);
              gui.progVal.setVisible(true);
            }
          });
        } catch (InvocationTargetException e) {} catch (InterruptedException e) {}
        final Color greenDark = Color.GREEN.darker();
        final Color dark = Color.GRAY;
        final Set<Step> selectedSteps = Sets.newHashSet();
        for (Entry<Step, JCheckBox> entry : gui.checkBoxes.entrySet()) {
          if (entry.getValue().isSelected()) {
            selectedSteps.add(entry.getKey());
          }
        }
        int i = 0;
        for (final Step step : stepsToRefresh) {
          if (step == null || gui.checkBoxes.get(step) == null || gui.varFields.get(step) == null) {
            continue;
          }
          final Map<Step, Map<Requirement, String>> variables = gui.getVariables();
          final int update = ++i;
          if (!step.checkIfOutputExists(variables) || gui.checkBoxes.get(step).isSelected()) {
            boolean check = step.hasRequirements(gui.proj, selectedSteps, variables);
            gui.descLabels.get(step).setForeground(check ? greenDark : Color.RED);
            gui.checkBoxes.get(step).setForeground(check ? greenDark : Color.RED);
            final Map<Requirement, JLabel> reqLbls = gui.requirementsLabels.get(step);
            try {
              SwingUtilities.invokeAndWait(new Runnable() {

                @Override
                public void run() {
                  for (Requirement req : step.getRequirements().getFlatRequirementsList()) {
                    boolean met = req.checkRequirement(variables.get(step).get(req), selectedSteps,
                                                       variables);
                    reqLbls.get(req).setForeground(met ? greenDark : Color.RED);
                  }
                  gui.progVal.setValue(update);
                }
              });
            } catch (InvocationTargetException e) {} catch (InterruptedException e) {}
          } else {
            final Map<Requirement, JLabel> reqLbls = gui.requirementsLabels.get(step);
            try {
              SwingUtilities.invokeAndWait(new Runnable() {

                @Override
                public void run() {
                  gui.checkBoxes.get(step).setSelected(false);
                  gui.alreadyRunLbls.get(step).setVisible(true);
                  gui.descLabels.get(step).setForeground(dark);
                  gui.checkBoxes.get(step).setForeground(dark);
                  for (Requirement req : step.getRequirements().getFlatRequirementsList()) {
                    reqLbls.get(req).setForeground(dark);
                  }
                  gui.progVal.setValue(update);
                }
              });
            } catch (InvocationTargetException e) {} catch (InterruptedException e) {}
          }
        }
        try {
          SwingUtilities.invokeAndWait(new Runnable() {

            @Override
            public void run() {
              gui.progVal.setValue(0);
              gui.progVal.setVisible(false);
            }
          });
        } catch (

        InvocationTargetException e) {} catch (InterruptedException e) {}

      }
    }).start();

  }

  private volatile boolean running = false;

  private JButton btnSelectAll;
  private JButton btnDeselectAll;
  private JButton btnSelectValid;
  private JProgressBar progVal;

  private Thread runThread;

  private JButton btnExportToText;
  private JButton btnExportToTextDefaults;
  private JButton btnOk;
  private JButton btnCancel;

  private void lockup(final boolean lock) {
    SwingUtilities.invokeLater(() -> {
      for (JCheckBox box : checkBoxes.values()) {
        box.setEnabled(!lock);
      }
      for (Map<Requirement, ? extends JComponent> flds : varFields.values()) {
        for (JComponent fld : flds.values()) {
          fld.setEnabled(!lock);
        }
      }
      for (ArrayList<JButton> btns : fileBtns.values()) {
        for (JButton btn : btns) {
          btn.setEnabled(!lock);
        }
      }
      btnSelectAll.setEnabled(!lock);
      btnDeselectAll.setEnabled(!lock);
      btnSelectValid.setEnabled(!lock);
      btnCancel.setEnabled(!lock);
      btnOk.setEnabled(!lock);
      btnExportToText.setEnabled(!lock);
      btnExportToTextDefaults.setEnabled(!lock);
    });
  }

  private void exportToText(final boolean useDefaults) {
    running = true;
    new Thread(new Runnable() {

      @Override
      public void run() {
        lockup(true);

        try {
          Set<Step> selectedSteps = getSelectedOptions();
          Map<Step, Map<Requirement, String>> variables = getVariables();
          if (checkRequirementsAndNotify(variables)) {
            StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");
            Set<Requirement.Flag> flags = EnumSet.noneOf(Requirement.Flag.class);
            for (Step step : selectedSteps) {
              flags.addAll(step.getFlags());
              String cmd = step.getCommandLine(proj, variables);
              GenvisisWorkflow.addStepInfo(output, step, cmd);
            }
            boolean hasQsub = Files.programExists("qsub");
            String file = proj.PROJECT_DIRECTORY.getValue() + "GenvisisPipeline";
            String suggFile = file + ext.getTimestampForFilename() + ".pbs";
            String command = output.toString();
            if (useDefaults || !hasQsub) {
              Qsub.qsubDefaults(suggFile, command);
            } else {
              if (hasQsub) {
                file = Qsub.qsubGUI(suggFile, command);
                if (file != null) {
                  if (!file.endsWith(".qsub") && !file.endsWith(".pbs")) {
                    file = file + ".pbs";
                  }
                }
              }
            }
            if (file != null) {
              proj.message("GenvisisPipeline commands written to " + file, "Command File Written",
                           JOptionPane.INFORMATION_MESSAGE);
            }
          }
        } catch (Exception e) {
          proj.getLog().reportException(e);
          ;
          proj.message("WARNING: Failed to create Genvisis workflow script. See log.");
        }

        end();
      }
    }).start();
  }

  public void nextStep(StepTask currentTask, FINAL_CODE returnCode, Set<Step> selectedSteps,
                       Map<Step, Map<Requirement, String>> variables) {
    Throwable e;
    Step currentStep = currentTask.getStep();
    progTasks.remove(currentStep);

    Step nextStep = null;
    if (returnCode == FINAL_CODE.FAILED) {
      e = currentTask.getFailureException();
      StringBuilder failureMessage = new StringBuilder("Error Occurred on Step ").append(currentStep.getName())
                                                                                 .append(":");
      if (e != null) {
        proj.getLog().reportException(e, 0);
        failureMessage.append("\n").append(e.getMessage());
      } else if (currentTask.getFailed()) {
        for (String msg : currentTask.getFailureMessages()) {
          failureMessage.append("\n").append(msg);
        }
      } else if (!currentStep.checkIfOutputExists(variables)) {
        failureMessage.append("\nUnknown error occurred.");
      }
      failureMessage.append("\nPlease check project log for more details.");
      String[] opts = {"Continue", "Retry", "Cancel"};
      int opt = JOptionPane.showOptionDialog(GenvisisWorkflowGUI.this, failureMessage.toString(),
                                             "Error!", JOptionPane.YES_NO_OPTION,
                                             JOptionPane.ERROR_MESSAGE, null, opts,
                                             opts[opts.length - 1]);
      if (opt == JOptionPane.CLOSED_OPTION || opt == opts.length - 1) { // closed or cancel
        nextStep = null;
      } else if (opt == 1) { // retry
        // TODO implement Step.cleanupAfterFailure()
        nextStep = currentStep;
      } else {
        nextStep = findNext(currentStep, selectedSteps);
      }
    } else if (returnCode == FINAL_CODE.CANCELLED) {
      boolean continueExec = death(currentStep, selectedSteps);
      if (!continueExec) {
        // TODO implement Step.cleanupAfterFailure()
        nextStep = null;
      } else {
        nextStep = findNext(currentStep, selectedSteps);
      }
    } else {
      nextStep = findNext(currentStep, selectedSteps);
    }
    if (nextStep != null) {
      runStep(nextStep, selectedSteps, variables);
    } else {
      end();
    }
  }

  private Step findNext(Step step, Set<Step> selectedSteps) {
    Iterator<Step> iter = selectedSteps.iterator();
    while (iter.hasNext()) {
      if (step.equals(iter.next()) && iter.hasNext()) {
        return iter.next();
      }
    }
    return null;
  }

  private void end() {
    lockup(false);
    running = false;
  }

  private void runStep(Step step, Set<Step> options,
                       Map<Step, Map<Requirement, String>> variables) {
    Task<Void, Void> stepTask = step.createTask(GenvisisWorkflowGUI.this, proj, variables, options);
    progTasks.put(step, stepTask);
    stepTask.execute();
  }

  private void run() {
    running = true;
    runThread = new Thread(new Runnable() {

      @Override
      public void run() {
        lockup(true);

        Set<Step> options = getSelectedOptions();
        Map<Step, Map<Requirement, String>> variables = getVariables();
        if (checkRequirementsAndNotify(variables)) {
          Step first = options.iterator().next();
          runStep(first, options, variables);
        } else {
          end();
        }
      }
    });
    runThread.start();
  }

  protected boolean death(Step step, Set<Step> options) {
    step.gracefulDeath(proj);
    // TODO remove message when gracefulDeath is implemented for each step
    JOptionPane.showMessageDialog(GenvisisWorkflowGUI.this,
                                  "Error - cleanup of cancelled steps is not implemented.  Please clean or remove any generated files and try again.",
                                  "Error", JOptionPane.ERROR_MESSAGE);
    boolean foundMore = false;
    for (Step step2 : steps.tailSet(step)) {
      if (options.contains(step2)) {
        foundMore = true;
        break;
      }
    }
    boolean continueExec = false;
    if (foundMore) {
      int opt = JOptionPane.showConfirmDialog(GenvisisWorkflowGUI.this,
                                              "A step was cancelled.  Do you wish to continue?",
                                              "Step Cancelled", JOptionPane.YES_NO_OPTION);
      if (opt == JOptionPane.YES_OPTION) {
        continueExec = true;
      }
    }
    return continueExec;
  }

  private boolean checkRequirementsAndNotify(Map<Step, Map<Requirement, String>> variables) {
    Set<Step> options = getSelectedOptions();

    boolean passesChecks = true;
    List<String> reqMsgs = Lists.newArrayList();
    for (Step step : options) {
      if (!getResources(step.getRequirements())
          || !step.hasRequirements(proj, options, variables)) {
        reqMsgs.add(checkBoxes.get(step).getText());
        passesChecks = false;
      }
    }
    if (!reqMsgs.isEmpty()) {
      StringBuilder msg = new StringBuilder("Prerequisite requirements are unmet for the following steps:");
      for (String str : reqMsgs) {
        msg.append("\n").append(str);
      }
      proj.message(msg.toString());
    }
    return passesChecks;
  }

  /**
   * @param requirementEntries entry for {@link Step} from variables.
   */
  private boolean getResources(RequirementSet requirements) {
    boolean allAvailable = true;
    for (Requirement requirement : requirements.getFlatRequirementsList()) {
      if (requirement instanceof Requirement.ResourceRequirement) {
        Requirement.ResourceRequirement resReq = (Requirement.ResourceRequirement) requirement;
        Resource resource = resReq.getResource();
        String reqPath = resource.get();
        if (reqPath == null) {
          allAvailable = false;
        }
      }
    }

    return allAvailable;
  }

  private Map<Step, Map<Requirement, String>> getVariables() {

    Map<Step, Map<Requirement, String>> returnVars = Maps.newHashMap();
    for (Entry<Step, Map<Requirement, JComponent>> entry : varFields.entrySet()) {
      Map<Requirement, String> values = Maps.newHashMap();
      returnVars.put(entry.getKey(), values);
      for (Entry<Requirement, JComponent> reqComp : entry.getValue().entrySet()) {
        String val = "";
        Requirement req = reqComp.getKey();
        JComponent j = reqComp.getValue();
        if (j instanceof JTextField) {
          val = ((JTextField) j).getText().trim();
        } else if (j instanceof JCheckBox) {
          val = Boolean.toString(((JCheckBox) j).isSelected());
        } else if (j instanceof JSpinner) {
          val = ((JSpinner) j).getValue().toString();
        } else if (j instanceof JComboBox) {
          val = ((JComboBox<?>) j).getSelectedItem().toString();
        } else if (j instanceof JList<?>) {
          val = Requirement.ListSelectionRequirement.createArgValString(((JList<?>) j).getSelectedValuesList());
        }
        values.put(req, val);
      }
    }
    return returnVars;
  }
}
