package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import net.miginfocom.swing.MigLayout;
import cnv.filesys.Project;
import common.Array;
import common.Files;
import common.ext;

public class PlinkExportOptions extends JDialog {

    private static final long serialVersionUID = 1L;

    private Project proj;
    
    private final JPanel contentPanel = new JPanel();
    private JTextField textFieldPlinkFileroot;
    private JLabel lblNameConflict;
    private JComboBox<String> comboBoxTargetMarkers;
    private JComboBox<String> comboBoxClusterFilters;

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

    /**
     * Create the dialog.
     */
    public PlinkExportOptions(Project proj) {
        this.proj = proj;
        setTitle("PLINK Export Options");
        setBounds(100, 100, 300, 250);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new MigLayout("", "[grow][grow]", "[][][][][][][][][]"));
        {
            JLabel lblClusterFilterCollection = new JLabel("Cluster Filter File:");
            contentPanel.add(lblClusterFilterCollection, "cell 0 0");
        }
        {
            comboBoxClusterFilters = new JComboBox<String>(getClusterFiltersOptions());
            comboBoxClusterFilters.setFont(comboBoxClusterFilters.getFont().deriveFont(Font.PLAIN));
            comboBoxClusterFilters.setSelectedItem(getClusterFiltersDefaultValue());
            contentPanel.add(comboBoxClusterFilters, "cell 0 1 2 1,growx");
        }
        {
            JLabel lblTargetMarkersFile = new JLabel("Target Markers File:");
            contentPanel.add(lblTargetMarkersFile, "cell 0 2");
        }
        {
            comboBoxTargetMarkers = new JComboBox<String>(getTargetMarkersOptions());
            comboBoxTargetMarkers.setFont(comboBoxTargetMarkers.getFont().deriveFont(Font.PLAIN));
            comboBoxTargetMarkers.setSelectedItem(getTargetMarkersDefaultValue());
            contentPanel.add(comboBoxTargetMarkers, "cell 0 3 2 1,growx");
        }
        {
            JLabel lblPlinkOutputFileroot = new JLabel("PLINK Output Fileroot:");
            contentPanel.add(lblPlinkOutputFileroot, "cell 0 4");
        }
        {
            textFieldPlinkFileroot = new JTextField();
            textFieldPlinkFileroot.setFont(textFieldPlinkFileroot.getFont().deriveFont(Font.PLAIN));
            textFieldPlinkFileroot.setText("plinkBinary");
            textFieldPlinkFileroot.addCaretListener(new CaretListener() {
                @Override
                public void caretUpdate(CaretEvent e) {
                    updatePlinkStatus();
                }
            });
            contentPanel.add(textFieldPlinkFileroot, "cell 0 5 2 1,growx");
        }
        {
            lblNameConflict = new JLabel("Error - files already exist!");
            lblNameConflict.setForeground(Color.RED.darker());
            lblNameConflict.setVisible(false);
            contentPanel.add(lblNameConflict, "flowx,cell 0 6 2 1,alignx right");
        }
        {
            chckbxOverwrite = new JCheckBox();
            chckbxOverwrite.setActionCommand("overwrite");
            chckbxOverwrite.addActionListener(buttonListener);
            chckbxOverwrite.setText("Overwrite");
            chckbxOverwrite.setEnabled(false);
            contentPanel.add(chckbxOverwrite, "cell 0 6,alignx right");
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
                close(false);
            }
        }
    };

    
    private void close(boolean cancelled) {
        this.cancelled = cancelled;
        this.setVisible(false);
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
        if (new File(plinkRootDir + ".bed").exists() || new File(plinkRootDir + ".fam").exists() || new File(plinkRootDir + ".bim").exists()) {
            return false;
        }
        return true;
    }
    
    private Object getClusterFiltersDefaultValue() {
        return proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME);
    }

    private String[] getClusterFiltersOptions() {
        return Array.addStrToArray("(--Do not apply any cluster filter--)", Files.list(proj.DATA_DIRECTORY.getValue(false, true), null, ext.removeDirectoryInfo(proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME)), false, proj.JAR_STATUS.getValue()));
    }
    
    private Object getTargetMarkersDefaultValue() {
        return proj.TARGET_MARKERS_FILENAME.getValue();
    }
    
    private String[] getTargetMarkersOptions() {
        return new String[]{proj.TARGET_MARKERS_FILENAME.getValue()};
    }
    
    public String getClusterFilterSelection() {
        String value = (String) comboBoxClusterFilters.getSelectedItem();
        if ("(--Do not apply any cluster filter--)".equals(value)) {
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

    public boolean getCancelled() {
        return cancelled;
    }
    
}
