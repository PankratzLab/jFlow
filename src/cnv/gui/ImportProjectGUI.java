package cnv.gui;

import java.awt.BorderLayout;

import javax.swing.AbstractAction;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;

import java.awt.Font;

import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.Action;

import common.Files;
import common.Grafik;
import common.ext;

import java.awt.Insets;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;

public class ImportProjectGUI extends JDialog {

    private static final long serialVersionUID = 1L;
    
    private final JPanel contentPanel = new JPanel();
    private JTextField txtFldProjName;
    private JTextField txtFldProjDir;
    private JTextField txtFldSampleDataFile;
    private JLabel lblFoundProjectStatus;
    private JLabel lblFoundSampleDataStatus;
    private JLabel lblFoundSamplesStatus;
    private JLabel lblFoundSampleListStatus;
    private JLabel lblFoundMarkerLookupStatus;
    private JLabel lblFoundDataStatus;
    private JLabel lblFoundMarkerListStatus;
    private JLabel lblFoundTransposedStatus;

    private ImageIcon redX = Grafik.getImageIcon("images/redx.png", true);
    private ImageIcon tick = Grafik.getImageIcon("images/tick.png", true);
    
    
    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        try {
            ImportProjectGUI dialog = new ImportProjectGUI();
            dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
            dialog.setVisible(true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    

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
            } else if (e.getActionCommand().equals("SAMPLEDATA")) {
                String txt = txtFldSampleDataFile.getText().trim();
                if (!txt.equals("") && !txt.equals("./")) {
                    curr = txt;
                } else {
                    txt = txtFldProjDir.getText().trim();
                    if (!txt.equals("") && !txt.equals("./")) {
                        curr = txt;
                    }
                }
                jfc.setCurrentDirectory(new File(curr));
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setDialogTitle("Select SampleData File");
                jfc.setMultiSelectionEnabled(false);
                fld = txtFldSampleDataFile;
            }            
            int val = jfc.showDialog(ImportProjectGUI.this, "Select");
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

    /**
     * Create the dialog.
     */
    public ImportProjectGUI() {
        setTitle("Genvisis: Import Existing Project Files");
        
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setBounds(100, 100, 550, 500);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new MigLayout("", "[grow][10px:10px:10px][grow 70][grow 40]", "[grow][][grow][][][][][][grow][grow][][][][grow][]"));
        CaretListener caretListener = new CaretListener() {
            @Override
            public void caretUpdate(CaretEvent arg0) {
                updateFound(getStatuses());
            }
        };
        Insets btnInsets = new Insets(-1, 3, -2, 3);
        {
            JLabel lblGenvisisProjectImport = new JLabel("Genvisis Project Import");
            lblGenvisisProjectImport.setFont(new Font("Arial", Font.BOLD, 16));
            contentPanel.add(lblGenvisisProjectImport, "cell 0 0 5 1,alignx center");
        }
        {
            JSeparator separator = new JSeparator();
            contentPanel.add(separator, "cell 0 1 5 1,growx");
        }
        {
            JLabel lblProjectName = new JLabel("Project Name:");
            contentPanel.add(lblProjectName, "cel   l 0 3,alignx right");
        }
        {
            txtFldProjName = new JTextField("New Project");
            txtFldProjName.setColumns(10);
            contentPanel.add(txtFldProjName, "cell 2 3,growx");
        }
        {
            JLabel lblProjectDirectory = new JLabel("Project Directory:");
            contentPanel.add(lblProjectDirectory, "cell 0 4,alignx right");
        }
        {
            txtFldProjDir = new JTextField("./");
            txtFldProjDir.setColumns(10);
            txtFldProjDir.addCaretListener(caretListener);
            contentPanel.add(txtFldProjDir, "flowx,cell 2 4,growx");
        }
        {
            JButton btnSelectProjDir = new JButton(fileSelectAction);
            btnSelectProjDir.setText("...");
            btnSelectProjDir.setMargin(btnInsets);
            btnSelectProjDir.setActionCommand("PROJECT");
            contentPanel.add(btnSelectProjDir, "cell 2 4");
        }
        {
            JLabel lblSampleDataFile = new JLabel("<html><code>SampleData</code> File:</html>");
            contentPanel.add(lblSampleDataFile, "cell 0 5,alignx right");
        }
        {
            txtFldSampleDataFile = new JTextField("./");
            txtFldSampleDataFile.setColumns(10);
            txtFldSampleDataFile.addCaretListener(caretListener);
            contentPanel.add(txtFldSampleDataFile, "flowx,cell 2 5,growx");
        }
        {
            JButton btnSelectSampleDataFile = new JButton(fileSelectAction);
            btnSelectSampleDataFile.setText("...");
            btnSelectSampleDataFile.setMargin(btnInsets);
            btnSelectSampleDataFile.setActionCommand("SAMPLEDATA");
            contentPanel.add(btnSelectSampleDataFile, "cell 2 5");
        }
        {
            JSeparator separator = new JSeparator();
            contentPanel.add(separator, "cell 0 7 4 1,growx");
        }
        {
            JPanel panel = new JPanel();
            contentPanel.add(panel, "cell 0 8 4 7,grow");
            panel.setLayout(new MigLayout("", "[grow][][10px:10px:10px][][grow]", "[grow][][][][][grow]"));
            {
                JLabel lblValidProjectDirectory = new JLabel("Valid Project Directory:");
                panel.add(lblValidProjectDirectory, "cell 0 1,alignx right");
                lblFoundProjectStatus = new JLabel(redX);
                panel.add(lblFoundProjectStatus, "cell 1 1,alignx left");
            }
            {
                JLabel lblValidSampledataFile = new JLabel("<html>Valid <code>SampleData</code> File:</html>");
                panel.add(lblValidSampledataFile, "cell 3 1,alignx right");
                lblFoundSampleDataStatus = new JLabel(redX);
                panel.add(lblFoundSampleDataStatus, "cell 4 1,alignx left");
            }
            {
                JLabel lblFoundSamplesDirectory = new JLabel("<html>Found <code>samples/</code> Directory:</html>");
                panel.add(lblFoundSamplesDirectory, "cell 0 2,alignx right");
                lblFoundSamplesStatus = new JLabel(redX);
                panel.add(lblFoundSamplesStatus, "cell 1 2,alignx left");
            }
            {
                JLabel lblFoundSampleList = new JLabel("<html>Found <code>SampleList</code> File:</html>");
                panel.add(lblFoundSampleList, "cell 3 2,alignx right");
                lblFoundSampleListStatus = new JLabel(redX);
                panel.add(lblFoundSampleListStatus, "cell 4 2,alignx left");
            }
            {
                JLabel lblFoundTransposedDirectory = new JLabel("<html>Found <code>transposed/</code> Directory:</html>");
                panel.add(lblFoundTransposedDirectory, "cell 0 3,alignx right");
                lblFoundTransposedStatus = new JLabel(redX);
                panel.add(lblFoundTransposedStatus, "cell 1 3,alignx left");
            }
            {
                JLabel lblFoundMarkerList = new JLabel("<html>Found <code>MarkerSet</code> File:</html>");
                panel.add(lblFoundMarkerList, "cell 3 3,alignx right");
                lblFoundMarkerListStatus = new JLabel(redX);
                panel.add(lblFoundMarkerListStatus, "cell 4 3,alignx left");
            }
            {
                JLabel lblFoundDataDirectory = new JLabel("<html>Found <code>data/</code> Directory:</html>");
                panel.add(lblFoundDataDirectory, "cell 0 4,alignx right");
                lblFoundDataStatus = new JLabel(redX);
                panel.add(lblFoundDataStatus, "cell 1 4,alignx left");
            }
            {
                JLabel lblFoundMarkerlookupFile = new JLabel("<html>Found <code>MarkerLookup</code> File:</html>");
                panel.add(lblFoundMarkerlookupFile, "cell 3 4,alignx right");
                lblFoundMarkerLookupStatus = new JLabel(redX);
                panel.add(lblFoundMarkerLookupStatus, "cell 4 4,alignx left");
            }
        }
        {
            JPanel buttonPane = new JPanel();
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            buttonPane.setLayout(new MigLayout("", "[grow][47px][65px]", "[][23px]"));
            {
                JSeparator separator = new JSeparator();
                buttonPane.add(separator, "cell 0 0 3 1,growx");
            }
            {
                JButton okButton = new JButton("OK");
                okButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent arg0) {
                    }
                });
                okButton.setActionCommand("OK");
                buttonPane.add(okButton, "cell 1 1,alignx left,aligny top");
                getRootPane().setDefaultButton(okButton);
            }
            {
                JButton cancelButton = new JButton("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent arg0) {
                    }
                });
                cancelButton.setActionCommand("Cancel");
                buttonPane.add(cancelButton, "cell 2 1,alignx left,aligny top");
            }
        }
    }
    
    private boolean checkProjectName() {
        String name = txtFldProjName.getText().trim();
        return !name.equals("") && !"New Project".equals(name);
    }
    
    private boolean[] getStatuses() {
        String baseDir = ext.verifyDirFormat(txtFldProjDir.getText().trim());
        String sampleData = txtFldSampleDataFile.getText().trim();
        boolean foundProject = Files.exists(baseDir);
        File sampDFile = new File(sampleData);
        boolean foundSampleData = sampDFile.exists() && sampDFile.isFile() && sampDFile.canRead();
        boolean foundSamples = Files.exists(baseDir + "samples/");
        boolean foundData = Files.exists(baseDir + "data/");
        boolean foundTransposed = Files.exists(baseDir + "transposed/");
        boolean foundSampleList = Files.exists(baseDir + "data/samples.bis");
        boolean foundMarkerList = Files.exists(baseDir + "data/markers.bim");
        boolean foundMarkerLookup = Files.exists(baseDir + "data/markerLookup.bml");
        return new boolean[]{
                /*0*/ foundProject,
                /*1*/ foundSampleData,
                /*2*/ foundSamples,
                /*3*/ foundData,
                /*4*/ foundTransposed,
                /*5*/ foundSampleList,
                /*6*/ foundMarkerList,
                /*7*/ foundMarkerLookup
        };
    }
    
    private void updateFound(final boolean[] statuses) {
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                lblFoundProjectStatus.setIcon(statuses[0] ? tick : redX);
                lblFoundSampleDataStatus.setIcon(statuses[1] ? tick : redX);
                lblFoundSamplesStatus.setIcon(statuses[2] ? tick : redX);
                lblFoundDataStatus.setIcon(statuses[3] ? tick : redX);
                lblFoundTransposedStatus.setIcon(statuses[4] ? tick : redX);
                lblFoundSampleListStatus.setIcon(statuses[5] ? tick : redX);
                lblFoundMarkerListStatus.setIcon(statuses[6] ? tick : redX);
                lblFoundMarkerLookupStatus.setIcon(statuses[7] ? tick : redX);
            }
        });
    }
    
    
    
}
