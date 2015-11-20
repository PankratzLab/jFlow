package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import net.miginfocom.swing.MigLayout;
import cnv.Launch;
import cnv.filesys.Project;
import cnv.manage.GenvisisPipeline;
import cnv.manage.GenvisisPipeline.RequirementInputType;
import cnv.manage.GenvisisPipeline.STEP;
import common.Array;
import common.Grafik;
import common.ext;

import javax.swing.JSeparator;

import java.awt.event.ActionListener;

public class GenvisisPipelineGUI extends JDialog {

    private static final long serialVersionUID = 1L;

    private final JPanel contentPanel = new JPanel();
    

    private HashMap<STEP, JCheckBox> checkBoxes = new HashMap<STEP, JCheckBox>();
    private HashMap<STEP, JLabel> descLabels = new HashMap<STEP, JLabel>();
    private HashMap<STEP, ArrayList<JLabel>> requirementsLabels = new HashMap<STEP, ArrayList<JLabel>>();
    private HashMap<STEP, JAccordionPanel> panels = new HashMap<GenvisisPipeline.STEP, JAccordionPanel>();
    public ConcurrentHashMap<STEP, ArrayList<? extends JComponent>> varFields = new ConcurrentHashMap<GenvisisPipeline.STEP, ArrayList<? extends JComponent>>();
    public HashMap<STEP, JProgressBar> progBars = new HashMap<GenvisisPipeline.STEP, JProgressBar>();
    public HashMap<STEP, ArrayList<JButton>> fileBtns = new HashMap<GenvisisPipeline.STEP, ArrayList<JButton>>();
    public HashMap<STEP, JLabel> alreadyRunLbls = new HashMap<GenvisisPipeline.STEP, JLabel>();
    
    Project proj;
    
    private static final String TOP_LABEL = "Genvisis Project Pipeline:";
    
    volatile boolean cancelled = false;
    volatile boolean[] selected;
    STEP[] steps;
    
    /**
     * Create the dialog.
     */
    public GenvisisPipelineGUI(Project proj2, final Launch launch) {
        if (proj2 == null) {
            this.proj = createNewProject();
        } else {
            this.proj = proj2;
        }
        if (this.proj == null) {
            doClose();
            return;
        } else {
            launch.loadProjects();
            launch.setIndexOfCurrentProject(ext.removeDirectoryInfo(this.proj.getPropertyFilename()));
            this.proj = launch.loadProject();
        }
        this.proj.getLog().report("Launching Genvisis Project Pipeline");
        this.steps = GenvisisPipeline.getStepsForProject(this.proj);
        selected = Array.booleanArray(this.steps.length, true);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        JPanel optionPanel = new JPanel();
        JScrollPane scrollPane = new JScrollPane(optionPanel);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
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
            for (int i = 0; i < this.steps.length; i++) {
                final int index = i;
                JAccordionPanel panel = createPanel(index);
                optionPanel.add(panel, "cell 0 " + (i) + ", alignx left, growx");
                panels.put(this.steps[i], panel);
            }
        }
        {
            JPanel buttonPane = new JPanel();
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            buttonPane.setLayout(new MigLayout("", "[][][][][][][grow][47px][59px]", "[23px]"));
            
            JLabel lblSelect = new JLabel("Select:");
            buttonPane.add(lblSelect, "flowx,cell 0 0");
            
            btnSelectAll = new JButton("All");
            btnSelectAll.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent arg0) {
                    for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
                        entry.getValue().setSelected(true);
                    }
                    for (int i = 0; i < selected.length; i++) {
                        selected[i] = true;
                    }
                }
            });
            buttonPane.add(btnSelectAll, "cell 0 0");
            
            btnDeselectAll = new JButton("None");
            btnDeselectAll.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
                        entry.getValue().setSelected(false);
                    }
                    for (int i = 0; i < selected.length; i++) {
                        selected[i] = false;
                    }
                }
            });
            buttonPane.add(btnDeselectAll, "cell 1 0");
            
            JSeparator separator = new JSeparator();
            separator.setOrientation(SwingConstants.VERTICAL);
            buttonPane.add(separator, "cell 2 0,growy");
            
            JLabel lblCollapse = new JLabel("Collapse:");
            buttonPane.add(lblCollapse, "cell 3 0");
            
            JButton btnAll = new JButton("All");
            btnAll.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    for (Entry<STEP, JAccordionPanel> entry : panels.entrySet()) {
                        entry.getValue().shrink();
                    }
                }
            });
            buttonPane.add(btnAll, "cell 4 0");
            
            JButton btnNone = new JButton("None");
            btnNone.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    for (Entry<STEP, JAccordionPanel> entry : panels.entrySet()) {
                        entry.getValue().expand();
                    }
                }
            });
            buttonPane.add(btnNone, "cell 5 0");
            JButton okButton = new JButton("OK");
            okButton.setActionCommand("OK");
            okButton.setMnemonic(KeyEvent.VK_O);
            buttonPane.add(okButton, "cell 7 0,alignx left,aligny top");
            getRootPane().setDefaultButton(okButton);
            JButton cancelButton = new JButton("Close");
            cancelButton.setActionCommand("Close");
            cancelButton.setMnemonic(KeyEvent.VK_C);
            buttonPane.add(cancelButton, "cell 8 0,alignx left,aligny top");
            
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
                    } else {
                        run();
                    }
                }
            };
            okButton.addActionListener(listener);
            cancelButton.addActionListener(listener);
            
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
                HashMap<STEP, ArrayList<String>> variables = getVariables();
                for (int i = 0; i < steps.length; i++) {
                    STEP step = steps[i];
                    if (step.checkIfOutputExists(proj, variables)) {
                        checkBoxes.get(step).setSelected(false);
                        alreadyRunLbls.get(step).setVisible(true);
                        panels.get(step).shrink();
                        selected[i] = false;
                    }
                }
            }
        });
        refreshLabels();
        setBounds(100, 100, 700, 850);
        setTitle(TOP_LABEL);
        addWindowFocusListener(new WindowFocusListener() {
            @Override
            public void windowLostFocus(WindowEvent e) {}
            @Override
            public void windowGainedFocus(WindowEvent e) {
                launch.toFront();
                GenvisisPipelineGUI.this.toFront();
            }
        });
        // panels start out visible to help with spacing (otherwise the containing jscrollpane is too small)
//        for (JPanel panel : panels.values()) {
//            ((JAccordionPanel) panel).shrink();
//        }
    }
    
    protected void doClose() {
        // TODO check if running, may not be able to stop current step.  Ask if should interrupt.  Maybe do cleanup if desired?
        cancelled = true;
        this.setVisible(false);
        this.dispose();
    }

    AbstractAction fileSelectAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            String[] parts = e.getActionCommand().split(":");
            int stepIndex = Integer.parseInt(parts[0]);
            int fieldIndex = Integer.parseInt(parts[1]);
            JTextField fileField = (JTextField) varFields.get(GenvisisPipelineGUI.this.steps[stepIndex]).get(fieldIndex);
            
            String current = fileField.getText();
            
            String dir = current.equals("") ? proj.PROJECT_DIRECTORY.getValue(false, false) : ext.parseDirectoryOfFile(current); 
            JFileChooser chooser = new JFileChooser(dir);
            chooser.setMultiSelectionEnabled(false);
            RequirementInputType[][] reqs = GenvisisPipelineGUI.this.steps[stepIndex].reqTypes;
            int accum = 0;
            RequirementInputType type = null;
            outer: for (int i = 0; i < reqs.length; i++) {
                for (int j = 0; j < reqs[i].length; j++) {
                    if (reqs[i][j] == RequirementInputType.NONE) {
                        continue;
                    }
                    if (accum == fieldIndex) {
                        type = reqs[i][j];
                        break outer;
                    }
                    accum++;
                }
            }

            if (type == RequirementInputType.FILE) {
                chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
                chooser.setDialogTitle("Select File");
            } else if (type == RequirementInputType.DIR) {
                chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                chooser.setDialogTitle("Select Directory");
            }
            
            int retValue = chooser.showDialog(GenvisisPipelineGUI.this, "Select");
            
            if (retValue == JFileChooser.CANCEL_OPTION) {
                refreshLabels();
                return;
            } else {
                File newFile = chooser.getSelectedFile();
                String txt = ext.verifyDirFormat(newFile.getAbsolutePath());
                if (type == RequirementInputType.FILE) {
                    txt = txt.substring(0, txt.length() - 1);
                }
                fileField.setText(txt);
            }
            refreshLabels();
        }
    };
    
    private JAccordionPanel createPanel(final int index) {
        STEP step = this.steps[index];
        final JAccordionPanel panel = new JAccordionPanel();
        
        final JCheckBox chckbx = new JCheckBox();
        chckbx.setAction(new AbstractAction() {
            static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                selected[index] = chckbx.isSelected();
                refreshLabels();
            }
        });
        chckbx.setFont(chckbx.getFont().deriveFont(Font.PLAIN, 14));
        Grafik.scaleCheckBoxIcon(chckbx);
        chckbx.setSelected(selected[index]);
        chckbx.setToolTipText(step.stepDesc);
        chckbx.setText((index + 1) + ": " + step.stepName);
        
        checkBoxes.put(step, chckbx);
        panel.topPanel.add(chckbx, "cell 0 0");
        
        JProgressBar stepProgBar = new JProgressBar();
        progBars.put(step, stepProgBar);
        stepProgBar.setVisible(false);
        panel.topPanel.add(stepProgBar, "cell 1 0, alignx right, hidemode 3, split 1");
        JLabel alreadyRanLbl = new JLabel("Output Already Exists!");
        alreadyRanLbl.setVisible(false);
        alreadyRunLbls.put(step, alreadyRanLbl);
        panel.topPanel.add(alreadyRanLbl, "cell 1 0, alignx right, hidemode 3");

        JLabel descLbl = new JLabel("<html><center><p>" + step.stepDesc + "</p></center></html>");
        descLbl.setVerticalAlignment(SwingConstants.TOP);
        descLbl.setHorizontalAlignment(SwingConstants.LEFT);
        descLbl.setFont(descLbl.getFont().deriveFont(Font.PLAIN));
        descLabels.put(step, descLbl);
        
        panel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
        String rows = "[][]";
        for (int i = 0; i < step.getRequirements().length; i++) {
            rows = rows + "[]";
        }
        panel.contentPanel.setLayout(new MigLayout("", "[200px,grow]push[200px,grow]", rows));
        panel.contentPanel.add(descLbl, "cell 0 0");
        
        String[][] reqs = step.getRequirements();
        
        String reqSteps = "abcdefghijklmnop"; // probably won't need that many
        if (reqs.length > 0) {
            JLabel reqLbl = new JLabel("Requires:");
            panel.contentPanel.add(reqLbl, "cell 0 1");
            
            ArrayList<JLabel> reqLbls = new ArrayList<JLabel>();
            int rowIndex = 2;
            for (int i = 0; i < reqs.length; i++) {
                // AND
                
                for (int j = 0; j < reqs[i].length; j++) {
                    // OR
                    JLabel indLbl = new JLabel("" + (i + 1) + (reqs[i].length > 1 ? reqSteps.charAt(j) : "") + ". ");
                    indLbl.setFont(indLbl.getFont().deriveFont(Font.PLAIN, 9));
                    panel.contentPanel.add(indLbl, "gapleft 25, aligny top, split 1, cell 0 " + rowIndex);
                    JLabel requirementLbl = new JLabel("<html><p>" + reqs[i][j] + "</p></html>");
                    requirementLbl.setFont(requirementLbl.getFont().deriveFont(Font.PLAIN, 9));
                    panel.contentPanel.add(requirementLbl, "aligny top, cell 0 " + rowIndex);
                    reqLbls.add(requirementLbl);
                    rowIndex++;
                    if (j < reqs[i].length - 1) {
                        JLabel orLbl = new JLabel("OR");
                        orLbl.setFont(orLbl.getFont().deriveFont(Font.PLAIN, 10));
                        panel.contentPanel.add(orLbl, "gapleft 18, cell 0 " + rowIndex);
                        rowIndex++;
                    }
                    
                }
                requirementsLabels.put(step, reqLbls);
                
                if (i < reqs.length - 1) {
                    JLabel andLbl = new JLabel("AND");
                    andLbl.setFont(andLbl.getFont().deriveFont(Font.PLAIN, 10));
                    panel.contentPanel.add(andLbl, "gapleft 7, cell 0 " + rowIndex);
                    rowIndex++;
                }
                
            }

            RequirementInputType[][] inputTypes = step.getRequirementInputTypes();
            ArrayList<JComponent> reqInputFields = new ArrayList<JComponent>();
            int reqIndex = 0;
            rowIndex = 2;
            for (int i = 0; i < inputTypes.length; i++) {
                for (int j = 0; j < inputTypes[i].length; j++) {
                    
                    if (inputTypes[i][j] == RequirementInputType.BOOL) {
                        JCheckBox checkBox = new JCheckBox();
                        checkBox.setFont(checkBox.getFont().deriveFont(14));
                        Grafik.scaleCheckBoxIcon(checkBox);
                        checkBox.setVerticalAlignment(SwingConstants.TOP);
                        checkBox.setHorizontalAlignment(SwingConstants.RIGHT);
                        reqIndex++;
                        reqInputFields.add(checkBox);
                        panel.contentPanel.add(checkBox, "alignx right, aligny center, growx, gapleft 20, cell 1 " + rowIndex);
                    } else if (inputTypes[i][j] != RequirementInputType.NONE) {
                        JTextField textField = new JTextField();
//                        textField.setHorizontalAlignment(JTextField.RIGHT);
//                        textField.addKeyListener(new KeyAdapter() {
//                            @Override
//                            public void keyTyped(KeyEvent e) {
//                                super.keyTyped(e);
//                                refreshLabels();
//                            }
//                        });
                        textField.addCaretListener(new CaretListener() {
                            @Override
                            public void caretUpdate(CaretEvent e) {
                                refreshLabels();
                            }
                        });
                        textField.setText(step.getRequirementDefaults(proj)[reqIndex].toString());
                        reqIndex++;
                        reqInputFields.add(textField);
                        panel.contentPanel.add(textField, "alignx right, aligny center, growx, gapleft 20, split 1, cell 1 " + rowIndex);
                        if(inputTypes[i][j] == RequirementInputType.FILE || inputTypes[i][j] == RequirementInputType.DIR) {
                            JButton fileBtn = new JButton();
                            fileBtn.setAction(fileSelectAction);
                            fileBtn.setText("...");
                            fileBtn.setActionCommand((index) + ":" + (reqInputFields.size() - 1));
                            fileBtn.setMargin(new Insets(1, 4, 0, 3));
                            panel.contentPanel.add(fileBtn, "cell 1 " + rowIndex);
                            ArrayList<JButton> list = fileBtns.get(step);
                            if (list == null) {
                                list = new ArrayList<JButton>();
                                fileBtns.put(step, list);
                            }
                            list.add(fileBtn);
                        }
                    }
                    rowIndex++;
                    if (j < reqs[i].length - 1) {
                        rowIndex++;
                    }
                }
                rowIndex++;
            }
            varFields.put(step, reqInputFields);
            
        }
        
        return panel;
        
    }
    
    public boolean[] getSelectedOptions() {
        return selected;
    }

    public boolean getCancelled() {
        return cancelled;
    }
    
    public void startStep(STEP step) {
        if (alreadyRunLbls.get(step).isVisible()) {
            alreadyRunLbls.get(step).setVisible(false);
        }
        progBars.get(step).setString("Running...");
        progBars.get(step).setStringPainted(true);
        progBars.get(step).setIndeterminate(true);
        progBars.get(step).setVisible(true);
    }
    
    public void endStep(STEP step, boolean failed) {
        progBars.get(step).setString(failed ? "Failed!" : "Complete!");
        progBars.get(step).setIndeterminate(false);
    }
    
    
    public void refreshLabels() {
        new Thread(new Runnable() {
            @Override
            public void run() {
                final Color greenDark = Color.GREEN.darker();
                final Color dark = Color.GRAY;
                for (final STEP step : GenvisisPipelineGUI.this.steps) {
                    if (step == null || checkBoxes.get(step) == null || varFields.get(step) == null) {
                        continue;
                    }
                    HashMap<STEP, Boolean> selectedSteps = new HashMap<GenvisisPipeline.STEP, Boolean>();
                    for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
                        selectedSteps.put(entry.getKey(), entry.getValue().isSelected());
                    }
                    HashMap<STEP, ArrayList<String>> variables = getVariables();
                    if (!step.checkIfOutputExists(proj, variables) || checkBoxes.get(step).isSelected()) {
                        boolean check = step.hasRequirements(proj, selectedSteps, variables);
                        descLabels.get(step).setForeground(check ? greenDark : Color.RED);
                        checkBoxes.get(step).setForeground(check ? greenDark : Color.RED);
                        final ArrayList<JLabel> reqLbls = requirementsLabels.get(step);
                        final boolean[][] reqVals = step.checkRequirements(proj, selectedSteps, variables);
                        SwingUtilities.invokeLater(new Runnable() {
                            @Override
                            public void run() {
                                int lblIndex = 0;
                                for (int i = 0; i < reqVals.length; i++) {
                                    for (int j = 0; j < reqVals[i].length; j++) {
                                        reqLbls.get(lblIndex).setForeground(reqVals[i][j] ? greenDark : Color.RED);
                                        lblIndex++;
                                    }
                                }
                            }
                        });
                    } else {
                        final boolean[][] reqVals = step.checkRequirements(proj, selectedSteps, variables);
                        final ArrayList<JLabel> reqLbls = requirementsLabels.get(step);
                        SwingUtilities.invokeLater(new Runnable() {
                            @Override
                            public void run() {
                                int lblIndex = 0;
                                checkBoxes.get(step).setSelected(false);
                                alreadyRunLbls.get(step).setVisible(true);
                                descLabels.get(step).setForeground(dark);
                                checkBoxes.get(step).setForeground(dark);
                                for (int i = 0; i < reqVals.length; i++) {
                                    for (int j = 0; j < reqVals[i].length; j++) {
                                        reqLbls.get(lblIndex).setForeground(dark);
                                        lblIndex++;
                                    }
                                }
                            }
                        });
                    }
                }
                
            }
        }).start();
        
    }
    
    private volatile boolean running = false;

    private JButton btnSelectAll;

    private JButton btnDeselectAll;
    
    private void lockup(final boolean lock) {
        try {
            SwingUtilities.invokeAndWait(new Runnable() {
                @Override
                public void run() {
                    for (JCheckBox box : checkBoxes.values()) {
                        box.setEnabled(!lock);
                    }
                    for (ArrayList<? extends JComponent> flds : varFields.values()) {
                        for (JComponent fld : flds) {
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
                }
            });
        } catch (InvocationTargetException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
    
    private void run() {
        running = true;
        new Thread(new Runnable() {
            @Override
            public void run() {
                lockup(true);
                
                boolean[] options = getSelectedOptions();

                HashMap<STEP, Boolean> selectedSteps = new HashMap<GenvisisPipeline.STEP, Boolean>();
                for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
                    selectedSteps.put(entry.getKey(), entry.getValue().isSelected());
                }
                HashMap<STEP, ArrayList<String>> variables = getVariables();
                if (checkRequirementsAndNotify(selectedSteps, variables)) {
                    for (int i = 0; i < options.length; i++) {
                        if (options[i]) {
                            startStep(GenvisisPipelineGUI.this.steps[i]);
                            Exception e = null;
                            try {
                                GenvisisPipelineGUI.this.steps[i].run(proj, variables);
                            } catch (Exception e1) {
                                e = e1;
                            }
                            boolean failed = false;
                            if (e != null || GenvisisPipelineGUI.this.steps[i].getFailed() || !GenvisisPipelineGUI.this.steps[i].checkIfOutputExists(proj, variables)) {
                                failed = true;
                            }
                            endStep(GenvisisPipelineGUI.this.steps[i], failed);
                            StringBuilder failureMessage = new StringBuilder("Error Occurred on Step ").append(i + 1);
                            if (e != null) {
                                proj.getLog().reportException(e);
                                failureMessage.append("\n").append(e.getMessage());
                            } else if (GenvisisPipelineGUI.this.steps[i].getFailed()) {
                                for (String msg : GenvisisPipelineGUI.this.steps[i].getFailureMessages()) {
                                    failureMessage.append("\n").append(msg);
                                }
                            } else if (!GenvisisPipelineGUI.this.steps[i].checkIfOutputExists(proj, variables)) {
                                failureMessage.append("\nUnknown error occurred.");
                            }
                            if (failed) {
                                failureMessage.append("\nPlease check project log for more details.");
                                String[] opts = {"Continue", "Retry", "Cancel"};
                                int opt = JOptionPane.showOptionDialog(GenvisisPipelineGUI.this, failureMessage.toString(), "Error!", JOptionPane.YES_NO_OPTION, JOptionPane.ERROR_MESSAGE, null, opts, opts[2]);
                                if (opt == JOptionPane.CLOSED_OPTION || opt == 2) { // closed or cancel
                                    for (STEP step : GenvisisPipelineGUI.this.steps) {
                                        step.resetRun();
                                    }
                                    break;
                                } else if (opt == 1) { // retry
                                    GenvisisPipelineGUI.this.steps[i].resetRun();
                                    i--;
                                } else {
                                    // continue
                                }
                            }
                            
                        }
                    }
                }
                
                lockup(false);
                running = false;
            }
        }).start();
        
    }
    
    private boolean checkRequirementsAndNotify(HashMap<STEP, Boolean> selectedSteps, HashMap<STEP, ArrayList<String>> variables) {
        boolean[] options = getSelectedOptions();
        
        ArrayList<String> reqMsgs = new ArrayList<String>();
        for (int i = 0; i < options.length; i++) {
            if (options[i]) {
                STEP step = this.steps[i];
                if (!step.hasRequirements(proj, selectedSteps, variables)) {
                    reqMsgs.add((i + 1) + ". " + step.stepName);
                }
            }
        }
        boolean retVal = true;
        if (reqMsgs.size() > 0) {
            StringBuilder msg = new StringBuilder("Prerequisite requirements are unmet for the following steps:");
            for (String str : reqMsgs) {
                msg.append("\n").append(str);
            }
            proj.message(msg.toString());
            retVal = false;
        }
        return retVal;
    }
    
    private HashMap<STEP, ArrayList<String>> getVariables() {
        HashMap<STEP, ArrayList<String>> returnVars = new HashMap<GenvisisPipeline.STEP, ArrayList<String>>();
        for (Entry<STEP, ArrayList<? extends JComponent>> entry : varFields.entrySet()) {
            ArrayList<String> values = new ArrayList<String>();
            returnVars.put(entry.getKey(), values);
            for (JComponent j : entry.getValue()) {
                String val = "";
                if (j instanceof JTextField) {
                    val = ((JTextField) j).getText().trim();
                } else if (j instanceof JCheckBox) {
                    val = ((JCheckBox) j).isSelected() + "";
                } else if (j instanceof JSpinner) {
                    val = ((JSpinner) j).getValue().toString();
                }
                values.add(val);
            }
        }
        return returnVars;
    }
    
    private Project createNewProject() {
        ProjectCreationGUI createGUI = new ProjectCreationGUI();
        createGUI.setModal(true);
        createGUI.setVisible(true);
        
        if (createGUI.wasCancelled()) {
            return null;
        } else {
            return createGUI.getCreatedProject();
        }
    }
    

}
