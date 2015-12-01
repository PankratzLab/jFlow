package cnv.gui;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.Insets;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import common.Array;
import common.Files;
import common.Grafik;
import common.ext;
import cnv.filesys.SourceFileHeaderData;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import cnv.manage.MitoPipeline;
import cnv.manage.SourceFileParser;
import net.miginfocom.swing.MigLayout;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import javax.swing.JComboBox;

public class ProjectCreationGUI extends JDialog {

    private static final long serialVersionUID = 1L;
    private static final String XY_TOOLTIP = "<html>Suggested values for X/Y correction:<br />Illumina: use the default of 1.<br />Affymetrix: use 100.<br />DBGAP: use 2000.</html>";
    private JPanel contentPane;
    private JTextField txtFldProjName;
    private JTextField txtFldProjDir;
    private JTextField txtFldSrcDir;
    private JTextField txtFldSrcExt;
    private JTextField txtFldTgtMkrs;
//    private JSpinner spinnerLrrSd;
    private Project proj;
    private volatile boolean cancelled = false;
    private Action fileSelectAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;

        @Override
        public void actionPerformed(ActionEvent e) {
            JFileChooser jfc = new JFileChooser();
            String curr = ext.pwd();
            JTextField fld = null;
            if (e.getActionCommand().equals("PROJECT")) {
                String txt = txtFldProjDir.getText().trim();
                if (!txt.equals("") && !txt.equals("./")) {
                    curr = txt;
                }
                jfc.setCurrentDirectory(new File(curr));
                jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                jfc.setDialogTitle("Select Project Directory");
                jfc.setMultiSelectionEnabled(false);
                fld = txtFldProjDir;
            } else if (e.getActionCommand().equals("SOURCE")) {
                String txt = txtFldSrcDir.getText().trim();
                if (!txt.equals("") && !txt.equals("./")) {
                    curr = txt;
                } else {
                    txt = txtFldProjDir.getText().trim();
                    if (!txt.equals("") && !txt.equals("./")) {
                        curr = txt;
                    }
                }
                jfc.setCurrentDirectory(new File(curr));
                jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                jfc.setDialogTitle("Select Source Directory");
                jfc.setMultiSelectionEnabled(false);
                fld = txtFldSrcDir;
            } else if (e.getActionCommand().equals("TARGET")) {
                String txt = txtFldProjDir.getText().trim();
                if (!txt.equals("") && !txt.equals("./")) {
                    curr = txt;
                }
                jfc.setCurrentDirectory(new File(curr));
                jfc.setDialogTitle("Select Target Markers File");
                jfc.setMultiSelectionEnabled(false);
                fld = txtFldTgtMkrs;
            } 
            
            int val = jfc.showDialog(ProjectCreationGUI.this, "Select");
            if (val == JFileChooser.APPROVE_OPTION) {
                File f = jfc.getSelectedFile();
                if (fld != null) {
                    String txt = ext.verifyDirFormat(f.getAbsolutePath());
                    try {
                        txt = ext.verifyDirFormat(f.getCanonicalPath());
                    } catch (IOException e1) {}
                    if (jfc.getFileSelectionMode() == FileChooser.FILES_ONLY && txt.length() > 1 && txt.endsWith("/")) {
                        txt = txt.substring(0, txt.length() - 1);
                    }
                    fld.setText(txt);
                }
                
            }
            
        }
    };
    
    private JComboBox<Project.ARRAY> comboBoxArrayType;
//    private JComboBox<String> comboBoxArrayType;
    private JLabel lblSrcFileStatus;
    private JSpinner spinnerXY;

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    ProjectCreationGUI frame = new ProjectCreationGUI();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }

    /**
     * Create the frame.
     */
    public ProjectCreationGUI() {
        setTitle("Genvisis: Create New Project");
        Project proj = new Project();
        
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setBounds(100, 100, 550, 500);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("", "[grow][10px:10px:10px][grow][grow]", "[grow][][grow][][][grow][][][][][grow][][][][grow][]"));
        
        JLabel lblGenvisisProjectCreation = new JLabel("Genvisis Project Creation");
        lblGenvisisProjectCreation.setFont(new Font("Arial", Font.BOLD, 16));
        contentPane.add(lblGenvisisProjectCreation, "cell 0 0 4 1,alignx center");
        
        JSeparator separator = new JSeparator();
        contentPane.add(separator, "cell 0 1 4 1,growx");
        
        JLabel lblProjectName = new JLabel("Project Name:");
        contentPane.add(lblProjectName, "cell 0 3,alignx trailing");
        
        txtFldProjName = new JTextField(proj.PROJECT_NAME.getDefaultValueString());
        contentPane.add(txtFldProjName, "cell 2 3,growx");
        txtFldProjName.setColumns(10);
        
        JLabel lblProjectDirectory = new JLabel("Project Directory:");
        contentPane.add(lblProjectDirectory, "cell 0 4,alignx trailing");
        
        txtFldProjDir = new JTextField(proj.PROJECT_DIRECTORY.getDefaultValueString());
        contentPane.add(txtFldProjDir, "flowx,cell 2 4,growx");
        txtFldProjDir.setColumns(10);
        
        JLabel lblSourceFileDirectory = new JLabel("Source File Directory:");
        contentPane.add(lblSourceFileDirectory, "cell 0 6,alignx trailing");
        
        CaretListener checkSource = new CaretListener() {
            @Override
            public void caretUpdate(CaretEvent e) {
                updateSourceFileNotice();
            }
        };
        
        txtFldSrcDir = new JTextField(proj.SOURCE_DIRECTORY.getDefaultValueString());
        txtFldSrcDir.addCaretListener(checkSource);
        contentPane.add(txtFldSrcDir, "flowx,cell 2 6,growx");
        txtFldSrcDir.setColumns(10);
        
        JLabel lblSourceFileExtension = new JLabel("Source File Extension:");
        contentPane.add(lblSourceFileExtension, "cell 0 7,alignx trailing");
        
        txtFldSrcExt = new JTextField(proj.SOURCE_FILENAME_EXTENSION.getDefaultValueString());
        txtFldSrcExt.addCaretListener(checkSource);
        contentPane.add(txtFldSrcExt, "cell 2 7,growx");
        txtFldSrcExt.setColumns(10);
        
        JLabel lblArrayType = new JLabel("Array Type:");
        contentPane.add(lblArrayType, "cell 0 8,alignx trailing");
        
        comboBoxArrayType = new JComboBox<Project.ARRAY>(Project.ARRAY.values());
//        comboBoxArrayType = new JComboBox<String>();
        comboBoxArrayType.setFont(comboBoxArrayType.getFont().deriveFont(Font.PLAIN));
        contentPane.add(comboBoxArrayType, "cell 2 8,growx");
        
        lblSrcFileStatus = new JLabel("");
        contentPane.add(lblSrcFileStatus, "cell 2 10,alignx right,aligny top");
        
        JLabel lblXyCorrectionRatio = new JLabel("X/Y Correction Ratio:");
        contentPane.add(lblXyCorrectionRatio, "cell 0 11,split 2, alignx right");
        contentPane.add(Grafik.getToolTipIconLabel(XY_TOOLTIP), "cell 0 11, alignx right");
        
        spinnerXY = new JSpinner();
        spinnerXY.setModel(new SpinnerNumberModel(proj.XY_SCALE_FACTOR.getDefaultValue().doubleValue(), 0.001, 10000000000d, 1.0));
        contentPane.add(spinnerXY, "cell 2 11,growx");
        
//        JLabel lblLogrRatioStddev = new JLabel("Log-R Ratio Std.Dev. Cut-off Threshold:");
//        contentPane.add(lblLogrRatioStddev, "cell 0 12,alignx trailing");
        
//        spinnerLrrSd = new JSpinner();
//        spinnerLrrSd.setModel(new SpinnerNumberModel(proj.LRRSD_CUTOFF.getDefaultValue().doubleValue(), 0.0, 3.0, 0.1));
//        contentPane.add(spinnerLrrSd, "cell 2 12,growx");
        
        JLabel lblTargetMarkersFile = new JLabel("[Optional] Target Markers File:");
        contentPane.add(lblTargetMarkersFile, "cell 0 12,alignx trailing");
        
        txtFldTgtMkrs = new JTextField(proj.TARGET_MARKERS_FILENAMES.getDefaultValueString());
        contentPane.add(txtFldTgtMkrs, "flowx,cell 2 12,growx");
        txtFldTgtMkrs.setColumns(10);
        
        Insets fileBtnInsets = new Insets(0, 3, 0, 3);
        JButton btnProjDir = new JButton(fileSelectAction);
        btnProjDir.setMargin(fileBtnInsets);
        btnProjDir.setText("...");
        btnProjDir.setActionCommand("PROJECT");
        contentPane.add(btnProjDir, "cell 2 4");
        
        JButton btnSrcDir = new JButton(fileSelectAction);
        btnSrcDir.setMargin(fileBtnInsets);
        btnSrcDir.setText("...");
        btnSrcDir.setActionCommand("SOURCE");
        contentPane.add(btnSrcDir, "cell 2 6");
        
        JButton btnTgtMkrs = new JButton(fileSelectAction);
        btnTgtMkrs.setMargin(fileBtnInsets);
        btnTgtMkrs.setText("...");
        btnTgtMkrs.setActionCommand("TARGET");
        contentPane.add(btnTgtMkrs, "cell 2 12");
        
        JSeparator separator_1 = new JSeparator();
        contentPane.add(separator_1, "cell 0 15 4 1,growx");
        
        JPanel panel = new JPanel();
        contentPane.add(panel, "south");
        panel.setLayout(new MigLayout("", "[grow][]", "[]"));
        
        JButton btnCreate = new JButton("Validate and Create");
        btnCreate.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (checkValues()) {
                    if (createProject()) {
                        doClose(false);
                    }
                }
            }
        });
        btnCreate.setMnemonic(KeyEvent.VK_V);
        panel.add(btnCreate, "flowx,cell 1 0");
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                doClose(true);
            }
        });
        btnCancel.setMnemonic(KeyEvent.VK_C);
        panel.add(btnCancel, "cell 1 0");
        
        updateSourceFileNotice();
    }
    
    private boolean checkValues() {
        String name = txtFldProjName.getText().trim();
        String projDir = txtFldProjDir.getText().trim();
        String srcDir = txtFldSrcDir.getText().trim();
        String srcExt = txtFldSrcExt.getText().trim();
//        String idHdr = txtFldIDHdr.getText().trim();
//        double lrrSd = ((Double)spinnerLrrSd.getValue()).doubleValue();
//        double xy = ((Double)spinnerXY.getValue()).doubleValue();
        String tgtMkrs = txtFldTgtMkrs.getText().trim();

        boolean validProjDir = false;
        try {
            File f = (new File(projDir));
            f.getCanonicalPath();
            validProjDir = true; // f.exists();
        } catch (IOException e) {}
        boolean validSrcDir = false;
        try {
            File f = (new File(srcDir));
            f.getCanonicalPath();
            validSrcDir = f.exists();
        } catch (IOException e) {}
        boolean validTgtMkrs = false;
        try {
            File f = (new File(tgtMkrs));
            f.getCanonicalPath();
            validTgtMkrs = true;
        } catch (IOException e) {}
        
        boolean[] checks = {!name.equals("") && !name.equals((new Project()).PROJECT_NAME.getDefaultValueString()), 
                validProjDir,
                validSrcDir,
                !srcExt.equals(""),
//                !idHdr.equals("") ,
                validTgtMkrs};//Files.exists(tgtMkrs)};
        
        if (Array.booleanArraySum(checks) < checks.length) {
            StringBuilder errorMsg = new StringBuilder();
            for (int i = 0; i < checks.length; i++) {
                if (!checks[i]) {
                    errorMsg.append("\n").append(getErrorFor(i));
                }
            }
            
            JOptionPane.showMessageDialog(null, errorMsg.toString(), "Error", JOptionPane.ERROR_MESSAGE);
            return false;
        }
        return true;
        
    }
    
    private void updateSourceFileNotice() {
        int cnt = getValidSourceFileCount();
        if (cnt <= 0) {
            lblSrcFileStatus.setText("No source files found");
            lblSrcFileStatus.setForeground(Color.RED);
        } else {
            lblSrcFileStatus.setText(cnt + " source file(s) found");
            lblSrcFileStatus.setForeground(Color.GREEN.darker());
        }
        lblSrcFileStatus.setVisible(true);
    }
    
    private int getValidSourceFileCount() {
        String srcDir = txtFldSrcDir.getText().trim();
        final String srcExt = txtFldSrcExt.getText().trim();

        boolean validSrcDir = false;
        try {
            File f = (new File(srcDir));
            f.getCanonicalPath();
            validSrcDir = f.exists();
        } catch (IOException e) {}
        if (!validSrcDir || "".equals(srcExt)) return -1;
        
        String[] files = (new File(srcDir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File arg0, String arg1) {
                return arg1.endsWith(srcExt);
            }
        });
        
        return files.length;
        // TODO check if files are valid, check if file headers are valid, check if file headers contain ID field
        // TODO use FinalReport header class (to be constructed)
    }
    
    private String getErrorFor(int index) {
        switch (index) {
            case 0: // Project Name
                return "Project name must have a value, and must not be 'New Project'";
            case 1: // Project dir
                return "Project directory must be a valid directory path";
            case 2: // src dir
                return "Source directory must be a valid and existing directory path";
            case 3: // src ext
                return "Source file extension cannot be empty";
            case 4: // id hdr
                return "ID header cannot be empty";
            case 5: // tgt mkrs
                return "Target markers must be a valid file";
            default:
                return "Unknown error - please check all options and try again";
        }
    }
    
    private boolean createProject() {
        String name = txtFldProjName.getText().trim();
        String projDir = txtFldProjDir.getText().trim();
        String srcDir = txtFldSrcDir.getText().trim();
        String srcExt = txtFldSrcExt.getText().trim();
//        String idHdr = txtFldIDHdr.getText().trim();
//        double lrrSd = ((Double)spinnerLrrSd.getValue()).doubleValue();
        double xy = ((Double)spinnerXY.getValue()).doubleValue();
        String tgtMkrs = txtFldTgtMkrs.getText().trim();
        
        HashMap<String, SourceFileHeaderData> headers = SourceFileHeaderData.validate(srcDir, srcExt, true, new common.Logger());
        if (headers == null) {
            // errors found in headers - check output and retry?
            return false;
        }
        SourceFileHeaderData reportHdr = null;
        for (SourceFileHeaderData d : headers.values()) {
            if (reportHdr == null) {
                reportHdr = d;
                break;
            }
        }
        // TODO do column assignment
        SourceFileHeaderGUI gui = new SourceFileHeaderGUI(reportHdr);
        gui.setModal(true);
        gui.setVisible(true);
        if (gui.wasCancelled()) {
            return false;
        }
        String sourceDelim = null;
        int sampCol = -100;
        String[] cols = null;
        for (SourceFileHeaderData d : headers.values()) {
            if (sourceDelim == null) {
                sourceDelim = d.getSourceFileDelimiter();
            }
            if (cols == null) {
                cols = d.cols;
            }
            d.colSnpIdent = gui.getSelectedSNPIndex();
            d.colSampleIdent = gui.getSelectedSampleID();
            if (sampCol == -100) {
                sampCol = d.colSampleIdent;
            }
            d.colGeno1 = gui.getSelectedGeno1();
            d.colGeno2 = gui.getSelectedGeno2();
            d.colGenoAB1 = gui.getSelectedAB1();
            d.colGenoAB2 = gui.getSelectedAB2();
            d.colBAF = gui.getSelectedBAF();
            d.colLRR = gui.getSelectedLRR();
            d.colGC = gui.getSelectedGC();
            d.colR = gui.getSelectedR();
            d.colTheta = gui.getSelectedTheta();
            d.colX = gui.getSelectedX();
            d.colY = gui.getSelectedY();
            d.colXRaw = gui.getSelectedXRaw();
            d.colYRaw = gui.getSelectedYRaw();
        }
        
        String path = MitoPipeline.initGenvisisProject();
        File file = new File(projDir);
        
        Project dummy = new Project();
        dummy.setGuiState(true);
        if (!file.exists()) {
            if (file.mkdirs()) {
                dummy.getLog().report("\n" + ext.getTime() + " Created directory " + projDir);
            } else {
                dummy.getLog().reportError("Error - failed to create  " + projDir + ", please manually create it unless it already exists");
                dummy.message("Error - failed to create  " + projDir + ", please manually create it unless it already exists");
                return false;
            }
        }
        String filename = path + name + MitoPipeline.PROJECT_EXT;
        if (Files.exists(filename)) {
            dummy.getLog().reportError("Project " + name + " already exists");
            dummy.message("Error - Project " + name + " already exists");
            return false;
        } else {
            Files.write((new Project()).PROJECT_NAME.getName() + "=" + name, filename);
        }
        Project actualProj = new Project(filename, false);
        actualProj.PROJECT_NAME.setValue(name);
        actualProj.PROJECT_DIRECTORY.setValue(projDir);
        actualProj.SOURCE_DIRECTORY.setValue(srcDir);
        actualProj.SOURCE_FILENAME_EXTENSION.setValue(srcExt);
//        actualProj.LRRSD_CUTOFF.setValue(lrrSd);
        actualProj.XY_SCALE_FACTOR.setValue(xy);
        actualProj.TARGET_MARKERS_FILENAMES.setValue(new String[]{ext.removeDirectoryInfo(tgtMkrs)});
        actualProj.ARRAY_TYPE.setValue((ARRAY) comboBoxArrayType.getSelectedItem());
        // if (abLookup != null && Files.exists(projectDirectory + abLookup)) {
        // proj.setProperty(proj.AB_LOOKUP_FILENAME, ext.removeDirectoryInfo(abLookup));
        // }
        actualProj.ID_HEADER.setValue(sampCol == SourceFileHeaderGUI.FILENAME_IND ? SourceFileParser.FILENAME_AS_ID_OPTION : cols[sampCol]);
        actualProj.SOURCE_FILE_DELIMITER.setValue(SOURCE_FILE_DELIMITERS.getDelimiter(sourceDelim));
        actualProj.saveProperties();
        actualProj.setSourceFileHeaders(headers);
        this.proj = actualProj;
        return true;
    }
    
    public Project getCreatedProject() {
        return this.proj;
    }
    
    private void doClose(boolean cancel) {
        this.cancelled = cancel;
        setVisible(false);
    }

    public boolean wasCancelled() {
        return cancelled;
    }

}
