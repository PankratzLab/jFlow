package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.File;

import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;
import cnv.filesys.Project;

import common.Array;
import common.Files;
import common.Grafik;
import common.ext;

public class PlinkExportOptions extends JDialog {

    private static final long serialVersionUID = 1L;

    private Project proj;
    
    private final JPanel contentPanel = new JPanel();
    private JTextField textFieldPlinkFileroot;
    private JLabel lblNameConflict;
    private JComboBox<String> comboBoxTargetMarkers;
    private JComboBox<String> comboBoxClusterFilters;
    private ButtonGroup exportTypeGroup;
    private JCheckBox chckbxOverwrite;

    private JButton okButton;
    private JButton cancelButton;

    private volatile boolean cancelled = false;

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        try {
            PlinkExportOptions dialog = new PlinkExportOptions(new Project());
            dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
            dialog.setVisible(true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static final String NO_CLUSTER_FILTERS = "(--Do not apply any cluster filter--)";
    private static final String ALL_MARKERS = "(--Export all markers--)";
    private static final String NEW_MARKERS_LIST = "(--Create new targetMarkers file--)";
    
    /**
     * Create the dialog.
     */
    public PlinkExportOptions(final Project proj) {
        this.proj = proj;
        setTitle("PLINK Export Options");
        setBounds(100, 100, 300, 250);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new MigLayout("", "[grow][grow]", "[][][][][][][][][][][][][]"));
        {
            lblPedigreeFile = new JLabel("Pedigree File:");
            contentPanel.add(lblPedigreeFile, "flowx,cell 0 0");
        }
        {
            textFieldPedigreeFile = new JTextField(proj.PEDIGREE_FILENAME.getValue());
            textFieldPedigreeFile.addKeyListener(new KeyAdapter() {
                @Override
                public void keyReleased(KeyEvent e) {
                    super.keyReleased(e);
                    refreshPedigreeStatus();
                }
                
            });
            contentPanel.add(textFieldPedigreeFile, "cell 0 1 2 1,growx");
            textFieldPedigreeFile.setColumns(10);
            refreshPedigreeStatus();
        }
        {
            lblExportFiletype = new JLabel("Export Filetype:");
            contentPanel.add(lblExportFiletype, "flowx,cell 0 6");
        }
        {
            exportTypeGroup = new ButtonGroup();
        }
        {
            rdbtnText = new JRadioButton("Text");
            exportTypeGroup.add(rdbtnText);
            contentPanel.add(rdbtnText, "flowx,cell 0 7 2 1,alignx center");
        }
        {
            JLabel lblClusterFilterCollection = new JLabel("Cluster Filter File:");
            contentPanel.add(lblClusterFilterCollection, "flowx,cell 0 4");
        }
        {
            comboBoxClusterFilters = new JComboBox<String>(getClusterFiltersOptions());
            comboBoxClusterFilters.setFont(comboBoxClusterFilters.getFont().deriveFont(Font.PLAIN));
            comboBoxClusterFilters.setSelectedItem(NO_CLUSTER_FILTERS);
            contentPanel.add(comboBoxClusterFilters, "cell 0 5 2 1,growx");
        }
        {
            JLabel lblTargetMarkersFile = new JLabel("Target Markers File:");
            contentPanel.add(lblTargetMarkersFile, "flowx,cell 0 2");
        }
        {
            comboBoxTargetMarkers = new JComboBox<String>(getTargetMarkersOptions());
            comboBoxTargetMarkers.addItemListener(new ItemListener() {
                @Override
                public void itemStateChanged(ItemEvent arg0) {
                    if (arg0.getStateChange() == ItemEvent.SELECTED) {
                        String val = (String) comboBoxTargetMarkers.getSelectedItem();
                        if (NEW_MARKERS_LIST.equals(val)) {
                            NewMarkerListDialog nmld = new NewMarkerListDialog(proj.getMarkerNames(), proj.PROJECT_DIRECTORY.getValue());
                            nmld.setModal(true);
                            nmld.setVisible(true);
                            if (nmld.getReturnCode() == JOptionPane.YES_OPTION) {
                                String mkrFile = nmld.getFileName();
                                proj.TARGET_MARKERS_FILENAMES.addValue(mkrFile);
                                comboBoxTargetMarkers.setModel(new DefaultComboBoxModel<String>(proj.TARGET_MARKERS_FILENAMES.getValue()));
                                comboBoxTargetMarkers.setSelectedItem(mkrFile);
                            } else {
                                comboBoxTargetMarkers.setSelectedItem(ALL_MARKERS);
                            }
                        }
                    }
                }
            });
            comboBoxTargetMarkers.setFont(comboBoxTargetMarkers.getFont().deriveFont(Font.PLAIN));
            comboBoxTargetMarkers.setSelectedItem(ALL_MARKERS);
            contentPanel.add(comboBoxTargetMarkers, "cell 0 3 2 1,growx");
        }
        {
            JLabel lblPlinkOutputFileroot = new JLabel("PLINK Output Fileroot:");
            contentPanel.add(lblPlinkOutputFileroot, "flowx,cell 0 8");
        }
        {
            textFieldPlinkFileroot = new JTextField();
            textFieldPlinkFileroot.setFont(textFieldPlinkFileroot.getFont().deriveFont(Font.PLAIN));
            textFieldPlinkFileroot.setText("plink");
            textFieldPlinkFileroot.addKeyListener(new KeyAdapter() {
                @Override
                public void keyReleased(KeyEvent e) {
                    super.keyReleased(e);
                    updatePlinkStatus();
                }
            });
            contentPanel.add(textFieldPlinkFileroot, "cell 0 9 2 1,growx");
        }
        {
            lblNameConflict = new JLabel("Error - files already exist!");
            lblNameConflict.setForeground(Color.RED.darker());
            lblNameConflict.setVisible(false);
            contentPanel.add(lblNameConflict, "flowx,cell 0 11 2 1,alignx right");
        }
        {
            chckbxOverwrite = new JCheckBox();
            chckbxOverwrite.setActionCommand("overwrite");
            chckbxOverwrite.addActionListener(buttonListener);
            chckbxOverwrite.setText("Overwrite");
            chckbxOverwrite.setEnabled(false);
            contentPanel.add(chckbxOverwrite, "cell 0 11 2 1,alignx right");
        }
        {
            rdbtnBinary = new JRadioButton("Binary");
            exportTypeGroup.add(rdbtnBinary);
            rdbtnBinary.setSelected(true);
            contentPanel.add(rdbtnBinary, "cell 0 7 2 1");
        }
        {
            tooltipPedigree = Grafik.getToolTipIconLabel("The pedigree file has the standard 6 columns (Family ID, Individual ID, Father's ID, Mother's ID, Sex, Affection/Phenotype) as well as a 7th column with the DNA/Sample ID from the raw data to match it up to. There is no header row. See manual for more detail.");
            contentPanel.add(tooltipPedigree, "cell 0 0");
        }
        {
            tooltipTgtMkrs = Grafik.getToolTipIconLabel("The target marker file is optional and has one marker name per row. If there are any additional columns, then they are ignored. See manual for more detail.");
            contentPanel.add(tooltipTgtMkrs, "cell 0 2");
        }
        {
            tooltipClusterFilters = Grafik.getToolTipIconLabel("The list of cluster filter files is generated from any file ending with \"*clusterFilters.ser\" in the project's data/ directory. The clusters can be manually added from within the ScatterPlot module. See manual for more detail.");
            contentPanel.add(tooltipClusterFilters, "cell 0 4");
        }
        {
            tooltipExportType = Grafik.getToolTipIconLabel("PLINK format can either be compressed in a binary format or in full text. See the PLINK website for more detail.");
            contentPanel.add(tooltipExportType, "cell 0 6");
        }
        {
            tooltipFileroot = Grafik.getToolTipIconLabel("The root of the filenames to be generated. If these already exist, you must click the checkbox to overwrite them.");
            contentPanel.add(tooltipFileroot, "cell 0 8");
        }
        {
            JPanel buttonPane = new JPanel();
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            {
                okButton = new JButton();
                okButton.addActionListener(buttonListener);
                okButton.setText("OK");
                okButton.setActionCommand("OK");
                buttonPane.add(okButton);
                getRootPane().setDefaultButton(okButton);
            }
            {
                cancelButton = new JButton();
                cancelButton.addActionListener(buttonListener);
                cancelButton.setText("Cancel");
                cancelButton.setActionCommand("Cancel");
                buttonPane.add(cancelButton);
            }
        }
        
        pack();
        updatePlinkStatus();
    }
    
    private void refreshPedigreeStatus() {
        String text = textFieldPedigreeFile.getText();
        if (new File(text).exists()) {
            textFieldPedigreeFile.setForeground(Color.GREEN.darker());
        } else {
            textFieldPedigreeFile.setForeground(Color.RED.darker());
        }
    }

    private ActionListener buttonListener = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent arg0) {
            if (arg0.getActionCommand().equals("Cancel")) {
                close(true);
            } else if (arg0.getActionCommand().equals("overwrite")) {
                if (!isValidPlinkRoot()) {
                    okButton.setEnabled(chckbxOverwrite.isSelected());
                }
            } else {
                if (checkPedigree()) {
                    close(false);
                } else {
                    JOptionPane.showMessageDialog(PlinkExportOptions.this, "Error - Pedigree file doesn't exist.", "Missing Pedigree File", JOptionPane.ERROR_MESSAGE, null);
                }
            }
        }
    };
    private JLabel lblExportFiletype;
    private JRadioButton rdbtnText;
    private JRadioButton rdbtnBinary;
    private JLabel lblPedigreeFile;
    private JTextField textFieldPedigreeFile;
    private JLabel tooltipPedigree;
    private JLabel tooltipTgtMkrs;
    private JLabel tooltipClusterFilters;
    private JLabel tooltipExportType;
    private JLabel tooltipFileroot;

    
    private void close(boolean cancelled) {
        this.cancelled = cancelled;
        this.setVisible(false);
    }
    
    protected boolean checkPedigree() {        
        String text = textFieldPedigreeFile.getText();
        return new File(text).exists();
    }

    private void updatePlinkStatus() {
        if (isValidPlinkRoot()) {
            lblNameConflict.setVisible(false);
            chckbxOverwrite.setEnabled(false);
            chckbxOverwrite.setSelected(false);
            okButton.setEnabled(true);
        } else {
            lblNameConflict.setVisible(true);
            chckbxOverwrite.setEnabled(true);
            okButton.setEnabled(false);
        }
    }
    
    private boolean isValidPlinkRoot() {
        String plinkRoot = textFieldPlinkFileroot.getText().trim();
        String plinkRootDir = proj.PROJECT_DIRECTORY.getValue() + plinkRoot;
        String[] exts = getFileExtensions();
        boolean exists = false;
        for (String ext : exts) {
            if (new File(plinkRootDir + ext).exists()) {
                exists = true;
                break;
            }
        }
        return !exists;
    }
    
    private String[] getClusterFiltersOptions() {
        return Array.addStrToArray(NO_CLUSTER_FILTERS, Files.list(proj.DATA_DIRECTORY.getValue(false, true), null, ext.removeDirectoryInfo(proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME)), false, proj.JAR_STATUS.getValue()));
    }
    
    private String[] getTargetMarkersOptions() {
        String[] values = proj.TARGET_MARKERS_FILENAMES.getValue();
        String[] retVals = new String[values.length + 2];
        for (int i = 0; i < values.length; i++) {
            retVals[i] = values[i];
        }
        retVals[retVals.length - 1] = NEW_MARKERS_LIST;
        retVals[retVals.length - 2] = ALL_MARKERS;
        return retVals;
    }
    
    public String getClusterFilterSelection() {
        String value = (String) comboBoxClusterFilters.getSelectedItem();
        if (NO_CLUSTER_FILTERS.equals(value)) {
            return null;
        } 
        return value;
    }

    public boolean getShouldWrite() {
        return isValidPlinkRoot() || chckbxOverwrite.isSelected();
    }
    
    public String getPlinkRoot() {
        return getShouldWrite() ? textFieldPlinkFileroot.getText().trim() : null;
    }
    
    /**
     * Checks the pedigree file text field value against the current value of the project's PEDIGREE_FILENAME property.  If they are different and if the new value is a valid file, set the PEDIGREE_FILENAME value to the new file path.
     * @return the value of proj.PEDIGREE_FILENAME
     */
    public String getPedigree() {
        String pedProj = proj.PEDIGREE_FILENAME.getValue();
        String currPed = textFieldPedigreeFile.getText();
        if ((new File(currPed)).exists() && !pedProj.equals(currPed)) {
            proj.PEDIGREE_FILENAME.setValue(currPed);
        }
        return proj.PEDIGREE_FILENAME.getValue();
    }
    
    public boolean exportAsBinary() {
        return rdbtnBinary.isSelected();
    }
    
    public String[] getFileExtensions() {
        return rdbtnBinary.isSelected() ? new String[]{".bim", ".bed", ".fam"} : new String[]{".map", ".ped"};
    }
    
    public boolean getCancelled() {
        return cancelled;
    }

    public String getTargetMarkersFile() {
        String val = (String) comboBoxTargetMarkers.getSelectedItem();
        if (ALL_MARKERS.equals(val)) {
            return null;
        } else if (NEW_MARKERS_LIST.equals(val)) {
            return null;
        } else {
            return val;
        }
    }
    
}
