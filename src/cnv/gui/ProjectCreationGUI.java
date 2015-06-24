package cnv.gui;

import java.awt.EventQueue;
import java.awt.Font;
import java.awt.Insets;
import java.io.File;
import java.io.IOException;

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

import common.Array;
import common.Files;
import common.ext;
import cnv.filesys.Project;
import cnv.manage.MitoPipeline;
import net.miginfocom.swing.MigLayout;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class ProjectCreationGUI extends JDialog {

    private static final long serialVersionUID = 1L;
    private JPanel contentPane;
    private JTextField txtFldProjName;
    private JTextField txtFldProjDir;
    private JTextField txtFldSrcDir;
    private JTextField txtFldSrcExt;
    private JTextField txtFldIDHdr;
    private JTextField txtFldTgtMkrs;
    private JSpinner spinnerLrrSd;
    private Project proj;
    volatile boolean cancelled = false;
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
                String txt = txtFldProjDir.getText().trim();
                if (!txt.equals("") && !txt.equals("./")) {
                    curr = txt;
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
                    if (e.getActionCommand().equals("TARGET") && txt.length() > 1 && txt.endsWith("/")) {
                        txt = txt.substring(0, txt.length() - 1);
                    }
                    fld.setText(txt);
                }
                
            }
            
        }
    };

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
        Project proj = new Project();
        
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setBounds(100, 100, 544, 332);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("", "[grow][10px:10px:10px][grow][grow]", "[][][grow][][][][][][][][][grow][][]"));
        
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
        contentPane.add(lblSourceFileDirectory, "cell 0 5,alignx trailing");
        
        txtFldSrcDir = new JTextField(proj.SOURCE_DIRECTORY.getDefaultValueString());
        contentPane.add(txtFldSrcDir, "flowx,cell 2 5,growx");
        txtFldSrcDir.setColumns(10);
        
        JLabel lblSourceFileExtension = new JLabel("Source File Extension:");
        contentPane.add(lblSourceFileExtension, "cell 0 6,alignx trailing");
        
        txtFldSrcExt = new JTextField(proj.SOURCE_FILENAME_EXTENSION.getDefaultValueString());
        contentPane.add(txtFldSrcExt, "cell 2 6,growx");
        txtFldSrcExt.setColumns(10);
        
        JLabel lblIdHeader = new JLabel("ID Header:");
        contentPane.add(lblIdHeader, "cell 0 7,alignx trailing");
        
        txtFldIDHdr = new JTextField(proj.ID_HEADER.getDefaultValueString());
        contentPane.add(txtFldIDHdr, "cell 2 7,growx");
        txtFldIDHdr.setColumns(10);
        
        JLabel lblLogrRatioStddev = new JLabel("Log-R Ratio Std.Dev. Cut-off Threshold:");
        contentPane.add(lblLogrRatioStddev, "cell 0 8,alignx trailing");
        
        spinnerLrrSd = new JSpinner();
        spinnerLrrSd.setModel(new SpinnerNumberModel(proj.LRRSD_CUTOFF.getDefaultValue().doubleValue(), 0.0, 3.0, 0.0));
        contentPane.add(spinnerLrrSd, "cell 2 8,growx");
        
        JLabel lblTargetMarkersFile = new JLabel("Target Markers File:");
        contentPane.add(lblTargetMarkersFile, "cell 0 9,alignx trailing");
        
        txtFldTgtMkrs = new JTextField(proj.TARGET_MARKERS_FILENAME.getDefaultValueString());
        contentPane.add(txtFldTgtMkrs, "flowx,cell 2 9,growx");
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
        contentPane.add(btnSrcDir, "cell 2 5");
        
        JButton btnTgtMkrs = new JButton(fileSelectAction);
        btnTgtMkrs.setMargin(fileBtnInsets);
        btnTgtMkrs.setText("...");
        btnTgtMkrs.setActionCommand("TARGET");
        contentPane.add(btnTgtMkrs, "cell 2 9");
        
        JSeparator separator_1 = new JSeparator();
        contentPane.add(separator_1, "cell 0 12 4 1,growx");
        
        JPanel panel = new JPanel();
        contentPane.add(panel, "south");
        panel.setLayout(new MigLayout("", "[grow][]", "[]"));
        
        JButton btnCreate = new JButton("Create");
        btnCreate.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (checkValues()) {
                    if (createProject()) {
                        doClose(false);
                    }
                }
            }
        });
        panel.add(btnCreate, "flowx,cell 1 0");
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                doClose(true);
            }
        });
        panel.add(btnCancel, "cell 1 0");
    }
    
    
    private boolean checkValues() {
        String name = txtFldProjName.getText().trim();
        String projDir = txtFldProjDir.getText().trim();
        String srcDir = txtFldSrcDir.getText().trim();
        String srcExt = txtFldSrcExt.getText().trim();
        String idHdr = txtFldIDHdr.getText().trim();
//        double lrrSd = ((Double)spinnerLrrSd.getValue()).doubleValue();
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
                !idHdr.equals("") ,
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
        String idHdr = txtFldIDHdr.getText().trim();
        double lrrSd = ((Double)spinnerLrrSd.getValue()).doubleValue();
        String tgtMkrs = txtFldTgtMkrs.getText().trim();

        String path = MitoPipeline.initGenvisisProject();
        File file = new File(projDir);
        
        Project dummy = new Project();
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
            // TODO Error - project already exists
            // Files.backup(ext.removeDirectoryInfo(filename), path, path + "backup/", false);
            dummy.getLog().reportError("Project " + name + " already exists");
            dummy.message("Error - Project " + name + " already exists");
//            dummy.getLog().reportError("Using project file " + filename + ", you may also specify project filename using the command line argument \"proj=\"");
            return false;
        } else {
            // log.report("Project properties file can be found at " + filename);
            // if (proj != null) {
            Files.write((new Project()).PROJECT_NAME.getName() + "=" + name, filename);
            // }
        }
        Project actualProj = new Project(filename, false);
        actualProj.PROJECT_NAME.setValue(name);
        actualProj.PROJECT_DIRECTORY.setValue(projDir);
        actualProj.SOURCE_DIRECTORY.setValue(srcDir);
        actualProj.SOURCE_FILENAME_EXTENSION.setValue(srcExt);
        actualProj.ID_HEADER.setValue(idHdr);
        actualProj.LRRSD_CUTOFF.setValue(lrrSd);
        actualProj.TARGET_MARKERS_FILENAME.setValue(ext.removeDirectoryInfo(tgtMkrs));
        // if (abLookup != null && Files.exists(projectDirectory + abLookup)) {
        // proj.setProperty(proj.AB_LOOKUP_FILENAME, ext.removeDirectoryInfo(abLookup));
        // }
        actualProj.saveProperties();
        this.proj = actualProj;
        return true;
    }
    
    public Project getCreatedProject() {
        return this.proj;
    }
    
    private void doClose(boolean cancel) {
        //
        cancelled = cancel;
        setVisible(false);
    }

}
