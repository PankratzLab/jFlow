package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;

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
import cnv.manage.KitAndKaboodle;
import cnv.manage.KitAndKaboodle.RequirementInputType;
import cnv.manage.KitAndKaboodle.STEP;
import common.Array;
import common.Grafik;
import common.ext;

public class KitAndKaboodleGUI extends JDialog {

    private static final long serialVersionUID = 1L;

    private final JPanel contentPanel = new JPanel();
    

    private HashMap<STEP, JCheckBox> checkBoxes = new HashMap<STEP, JCheckBox>();
    private HashMap<STEP, JLabel> descLabels = new HashMap<STEP, JLabel>();
    private HashMap<STEP, ArrayList<JLabel>> requirementsLabels = new HashMap<STEP, ArrayList<JLabel>>();
    private HashMap<STEP, JAccordionPanel> panels = new HashMap<KitAndKaboodle.STEP, JAccordionPanel>();
    public HashMap<STEP, ArrayList<? extends JComponent>> varFields = new HashMap<KitAndKaboodle.STEP, ArrayList<? extends JComponent>>();
    public HashMap<STEP, JProgressBar> progBars = new HashMap<KitAndKaboodle.STEP, JProgressBar>();
    public HashMap<STEP, ArrayList<JButton>> fileBtns = new HashMap<KitAndKaboodle.STEP, ArrayList<JButton>>();
    public HashMap<STEP, JLabel> alreadyRunLbls = new HashMap<KitAndKaboodle.STEP, JLabel>();
    
    Project proj;
    
    private static final String TOP_LABEL = "Genvisis Project Pipeline:";
    
    volatile boolean cancelled = false;
    volatile boolean[] selected;
    STEP[] steps;
    
    /**
     * Create the dialog.
     */
    public KitAndKaboodleGUI(Project proj2, final Launch launch) {
        if (proj2 == null) {
            this.proj = createNewProject();
        } else {
            this.proj = proj2;
        }
        if (this.proj == null) {
            doClose();
            return;
        } else {
            launch.setIndexOfCurrentProject(this.proj.getNameOfProject());
        }
        this.steps = KitAndKaboodle.getStepsForProject(this.proj);
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
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            JButton okButton = new JButton("OK");
            okButton.setActionCommand("OK");
            buttonPane.add(okButton);
            getRootPane().setDefaultButton(okButton);
            JButton cancelButton = new JButton("Close");
            cancelButton.setActionCommand("Close");
            buttonPane.add(cancelButton);
            
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
        for (int i = 0; i < this.steps.length; i++) {
            STEP step = this.steps[i];
            if (step.checkIfOutputExists(this.proj, varFields)) {
                checkBoxes.get(step).setSelected(false);
                alreadyRunLbls.get(step).setVisible(true);
                panels.get(step).shrink();
                selected[i] = false;
            }
        }
        refreshLabels();
        setBounds(100, 100, 660, 850);
        setTitle(TOP_LABEL);
        addWindowFocusListener(new WindowFocusListener() {
            @Override
            public void windowLostFocus(WindowEvent e) {}
            @Override
            public void windowGainedFocus(WindowEvent e) {
                launch.toFront();
                KitAndKaboodleGUI.this.toFront();
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
            JTextField fileField = (JTextField) varFields.get(KitAndKaboodleGUI.this.steps[stepIndex]).get(fieldIndex);
            
            String current = fileField.getText();
            
            String dir = current.equals("") ? proj.PROJECT_DIRECTORY.getValue(false, false) : ext.parseDirectoryOfFile(current); 
            JFileChooser chooser = new JFileChooser(dir);
            chooser.setMultiSelectionEnabled(false);
            RequirementInputType[][] reqs = KitAndKaboodleGUI.this.steps[stepIndex].reqTypes;
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
            
            int retValue = chooser.showDialog(KitAndKaboodleGUI.this, "Select");
            
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
        Color greenDark = Color.GREEN.darker();
        Color dark = Color.GRAY;
        for (STEP step : this.steps) {
            if (step == null || checkBoxes.get(step) == null || varFields.get(step) == null) {
                continue;
            }
            if (!step.checkIfOutputExists(proj, varFields) || checkBoxes.get(step).isSelected()) {
                if (step.hasRequirements(proj, checkBoxes, varFields)) {
                    descLabels.get(step).setForeground(greenDark);
                    checkBoxes.get(step).setForeground(greenDark);
                } else {
                    descLabels.get(step).setForeground(Color.RED);
                    checkBoxes.get(step).setForeground(Color.RED);
                }
                ArrayList<JLabel> reqLbls = requirementsLabels.get(step);
                boolean[][] reqVals = step.checkRequirements(proj, checkBoxes, varFields);
                int lblIndex = 0;
                for (int i = 0; i < reqVals.length; i++) {
                    for (int j = 0; j < reqVals[i].length; j++) {
                        reqLbls.get(lblIndex).setForeground(reqVals[i][j] ? greenDark : Color.RED);
                        lblIndex++;
                    }
                }
            } else {
                checkBoxes.get(step).setSelected(false);
                alreadyRunLbls.get(step).setVisible(true);
                descLabels.get(step).setForeground(dark);
                checkBoxes.get(step).setForeground(dark);
                boolean[][] reqVals = step.checkRequirements(proj, checkBoxes, varFields);
                int lblIndex = 0;
                ArrayList<JLabel> reqLbls = requirementsLabels.get(step);
                for (int i = 0; i < reqVals.length; i++) {
                    for (int j = 0; j < reqVals[i].length; j++) {
                        reqLbls.get(lblIndex).setForeground(dark);
                        lblIndex++;
                    }
                }
            }
        }
        
    }
    
    private volatile boolean running = false;
    
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
                
                if (checkRequirementsAndNotify()) {
                    for (int i = 0; i < options.length; i++) {
                        if (options[i]) {
                            startStep(KitAndKaboodleGUI.this.steps[i]);
                            Exception e = null;
                            try {
                                KitAndKaboodleGUI.this.steps[i].run(proj, varFields);
                            } catch (Exception e1) {
                                e = e1;
                            }
                            boolean failed = false;
                            if (e != null || KitAndKaboodleGUI.this.steps[i].getFailed() || !KitAndKaboodleGUI.this.steps[i].checkIfOutputExists(proj, varFields)) {
                                failed = true;
                            }
                            endStep(KitAndKaboodleGUI.this.steps[i], failed);
                            StringBuilder failureMessage = new StringBuilder("Error Occurred on Step ").append(i + 1);
                            if (e != null) {
                                proj.getLog().reportException(e);
                                failureMessage.append("\n").append(e.getMessage());
                            } else if (KitAndKaboodleGUI.this.steps[i].getFailed()) {
                                for (String msg : KitAndKaboodleGUI.this.steps[i].getFailureMessages()) {
                                    failureMessage.append("\n").append(msg);
                                }
                            } else if (!KitAndKaboodleGUI.this.steps[i].checkIfOutputExists(proj, varFields)) {
                                failureMessage.append("\nUnknown error occurred.");
                            }
                            if (failed) {
                                failureMessage.append("\nPlease check project log for more details.");
                                String[] opts = {"Continue", "Retry", "Cancel"};
                                int opt = JOptionPane.showOptionDialog(KitAndKaboodleGUI.this, failureMessage.toString(), "Error!", JOptionPane.YES_NO_OPTION, JOptionPane.ERROR_MESSAGE, null, opts, opts[2]);
                                if (opt == JOptionPane.CLOSED_OPTION || opt == 2) { // closed or cancel
                                    for (STEP step : KitAndKaboodleGUI.this.steps) {
                                        step.resetRun();
                                    }
                                    break;
                                } else if (opt == 1) { // retry
                                    KitAndKaboodleGUI.this.steps[i].resetRun();
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
    
    private boolean checkRequirementsAndNotify() {
        boolean[] options = getSelectedOptions();
        
        ArrayList<String> reqMsgs = new ArrayList<String>();
        for (int i = 0; i < options.length; i++) {
            if (options[i]) {
                STEP step = this.steps[i];
                if (!step.hasRequirements(proj, checkBoxes, varFields)) {
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
