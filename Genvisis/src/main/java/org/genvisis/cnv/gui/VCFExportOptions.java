package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.gui.UITools;

import net.miginfocom.swing.MigLayout;

public class VCFExportOptions extends JDialog {

  private static final String SAMP = "samp";
  private static final String OK = "OK";
  private static final String OVERWRITE = "overwrite";
  private static final String CANCEL = "Cancel";

  private static final long serialVersionUID = 1L;

  private final Project proj;

  private final JPanel contentPanel = new JPanel();
  private JTextField textFieldVCFFileroot;
  private JLabel lblNameConflict;
  private JComboBox<String> comboBoxTargetMarkers;
  private JComboBox<String> comboBoxClusterFilters;
  private JCheckBox chckbxOverwrite;
  private JLabel lblGcCutoff;
  private JSpinner gcSpinner;

  private JButton okButton;
  private JButton cancelButton;

  private volatile boolean cancelled = false;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      VCFExportOptions dialog = new VCFExportOptions(new Project());
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static final String NO_CLUSTER_FILTERS = "(--Do not apply any cluster filter--)";
  private static final String ALL_MARKERS = "(--Export all markers--)";
  private static final String NEW_MARKERS_LIST = "(--Create new targetMarkers file--)";
  private static final String LOAD_MARKERS_LIST = "(--Load new targetMarkers file--)";

  /**
   * Create the dialog.
   */
  public VCFExportOptions(final Project proj) {
    this.proj = proj;
    setTitle("VCF Export Options");
    setMinimumSize(new Dimension(100, 100));
    UITools.setSize(this, new Dimension(350, 450));
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(new MigLayout("", "[grow]", "[][][][][][][][][][][]"));
    {
      lblSampleFile = new JLabel("Samples to export file:");
      contentPanel.add(lblSampleFile, "flowx,cell 0 0");
    }
    {
      textFieldSampleListFile = new JTextField();
      textFieldSampleListFile.addKeyListener(new KeyAdapter() {

        @Override
        public void keyReleased(KeyEvent e) {
          super.keyReleased(e);
          refreshStatus(textFieldSampleListFile);
        }
      });
      contentPanel.add(textFieldSampleListFile, "flowx,cell 0 1,growx");
      textFieldSampleListFile.setColumns(10);
      refreshStatus(textFieldSampleListFile);
    }
    {
      btnSampleFile = new JButton();
      btnSampleFile.addActionListener(buttonListener);
      btnSampleFile.setActionCommand(SAMP);
      btnSampleFile.setText("...");
      contentPanel.add(btnSampleFile, "cell 0 1");
    }
    {
      btnSelectChromosomesTo = new JButton("Select Chromosomes to Export");
      btnSelectChromosomesTo.addActionListener(new ActionListener() {

        public void actionPerformed(ActionEvent arg0) {
          String[] chrs = new String[27];
          for (int i = 0; i < chrs.length; i++) {
            chrs[i] = Positions.CHR_CODES[i];
          }
          IncludeExcludeGUI iegui = new IncludeExcludeGUI(VCFExportOptions.this, chrs,
                                                          ArrayUtils.booleanArray(chrs.length,
                                                                                  true));
          iegui.setModal(true);
          iegui.setVisible(true);
          if (iegui.getCloseCode() == JOptionPane.OK_OPTION) {
            boolean[] sel = iegui.getSelected();
            chrsToWrite = new Byte[ArrayUtils.booleanArraySum(sel)];
            int ind = 0;
            for (int i = 0; i < sel.length; i++) {
              if (sel[i]) {
                chrsToWrite[ind++] = Positions.chromosomeNumber(chrs[i]);
              }
            }
          }
        }
      });
      contentPanel.add(btnSelectChromosomesTo, "cell 0 4,growx");
    }
    {
      JLabel lblClusterFilterCollection = new JLabel("Cluster Filter File:");
      contentPanel.add(lblClusterFilterCollection, "flowx,cell 0 5");
    }
    {
      comboBoxClusterFilters = new JComboBox<>();// (getClusterFiltersOptions());
      comboBoxClusterFilters.setFont(comboBoxClusterFilters.getFont().deriveFont(Font.PLAIN));
      comboBoxClusterFilters.setSelectedItem(NO_CLUSTER_FILTERS);
      contentPanel.add(comboBoxClusterFilters, "cell 0 6,growx");
    }
    {
      JLabel lblTargetMarkersFile = new JLabel("Markers to export file:");
      contentPanel.add(lblTargetMarkersFile, "flowx,cell 0 2");
    }
    {
      comboBoxTargetMarkers = new JComboBox<>();// (getTargetMarkersOptions());
      comboBoxTargetMarkers.addItemListener(new ItemListener() {

        @Override
        public void itemStateChanged(ItemEvent arg0) {
          if (arg0.getStateChange() == ItemEvent.SELECTED) {
            String val = (String) comboBoxTargetMarkers.getSelectedItem();
            if (NEW_MARKERS_LIST.equals(val)) {
              ListEditor nmld = ListEditor.createMarkerListCreator(proj);
              nmld.setModal(true);
              nmld.setVisible(true);
              if (nmld.getReturnCode() == JOptionPane.YES_OPTION) {
                String mkrFile = nmld.getFileName();
                proj.TARGET_MARKERS_FILENAMES.addValue(mkrFile);
                comboBoxTargetMarkers.setModel(new DefaultComboBoxModel<>(getTargetMarkersOptions()));
                comboBoxTargetMarkers.setSelectedItem(mkrFile);
                proj.saveProperties();
              } else {
                comboBoxTargetMarkers.setSelectedItem(ALL_MARKERS);
              }
            } else if (LOAD_MARKERS_LIST.equals(val)) {
              JFileChooser jfc = new JFileChooser(proj.PROJECT_DIRECTORY.getValue());
              jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
              jfc.setMultiSelectionEnabled(true);
              jfc.setDialogTitle("Import targetMarkers File");

              int selVal = jfc.showDialog(VCFExportOptions.this, "Select");
              if (selVal == JFileChooser.APPROVE_OPTION) {
                File[] newFiles = jfc.getSelectedFiles();
                for (File newFile : newFiles) {
                  proj.TARGET_MARKERS_FILENAMES.addValue(newFile.getAbsolutePath());
                }
                comboBoxTargetMarkers.setModel(new DefaultComboBoxModel<>(getTargetMarkersOptions()));
                comboBoxTargetMarkers.setSelectedItem(newFiles[0].getAbsolutePath());
                proj.saveProperties();
              } else {
                comboBoxTargetMarkers.setSelectedItem(ALL_MARKERS);
              }
            }
          }
        }
      });
      comboBoxTargetMarkers.setFont(comboBoxTargetMarkers.getFont().deriveFont(Font.PLAIN));
      comboBoxTargetMarkers.setSelectedItem(ALL_MARKERS);
      contentPanel.add(comboBoxTargetMarkers, "cell 0 3,growx");
    }
    {
      JLabel lblVCFOutputFileroot = new JLabel("VCF Output Fileroot:");
      contentPanel.add(lblVCFOutputFileroot, "flowx,cell 0 7");
    }
    {
      textFieldVCFFileroot = new JTextField();
      textFieldVCFFileroot.setFont(textFieldVCFFileroot.getFont().deriveFont(Font.PLAIN));
      textFieldVCFFileroot.setText(proj.PROJECT_DIRECTORY.getValue()
                                   + ext.replaceWithLinuxSafeCharacters(proj.PROJECT_NAME.getValue()));
      textFieldVCFFileroot.addKeyListener(new KeyAdapter() {

        @Override
        public void keyReleased(KeyEvent e) {
          super.keyReleased(e);
          updateVCFStatus();
        }
      });
      contentPanel.add(textFieldVCFFileroot, "cell 0 8,growx");
    }
    {
      lblNameConflict = new JLabel("Error - files already exist!");
      lblNameConflict.setForeground(Color.RED.darker());
      lblNameConflict.setVisible(false);
      contentPanel.add(lblNameConflict, "flowx,cell 0 9,alignx right");
    }
    {
      chckbxOverwrite = new JCheckBox();
      chckbxOverwrite.setActionCommand(OVERWRITE);
      chckbxOverwrite.addActionListener(buttonListener);
      chckbxOverwrite.setText("Overwrite");
      chckbxOverwrite.setEnabled(false);
      contentPanel.add(chckbxOverwrite, "cell 0 9,alignx right");
    }
    {
      panel = new JPanel();
      contentPanel.add(panel, "cell 0 10,grow");
      panel.setLayout(new MigLayout("ins 0", "[][][grow]", "[][]"));
      lblGcCutoff = new JLabel("GC Cutoff:");
      panel.add(lblGcCutoff, "cell 0 0");
      {
        separator = new JSeparator();
        separator.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator, "cell 1 0 1 2,growy");
      }
      {
        chckbxSplitChromosomes = new JCheckBox("Split Chromosomes?");
        panel.add(chckbxSplitChromosomes, "cell 2 0,alignx right");
      }

      SpinnerNumberModel model = new SpinnerNumberModel(proj.GC_THRESHOLD.getValue().doubleValue(),
                                                        0, 1, 0.01);
      gcSpinner = new JSpinner(model);
      panel.add(gcSpinner, "cell 0 1");
      {
        chckbxExportchrIn = new JCheckBox("Export \"chr\" in contigs?");
        panel.add(chckbxExportchrIn, "cell 2 1,alignx right");
      }
      JFormattedTextField jftf = ((JSpinner.DefaultEditor) gcSpinner.getEditor()).getTextField();
      jftf.setColumns(4);
    }
    {
      tooltipExportType = Grafik.getToolTipIconLabel("<html><p width=\"380\">The minimum GC quality score for a sample to be exported.  NOT IMPLEMENTED YET.</p></html>");
      panel.add(tooltipExportType, "cell 0 0");
    }
    {
      tooltipPedigree = Grafik.getToolTipIconLabel("<html><p width=\"380\">A list of Sample IDs for which to export data. Leave blank to export data for all samples.</p></html>");
      contentPanel.add(tooltipPedigree, "cell 0 0");
    }
    {
      tooltipTgtMkrs = Grafik.getToolTipIconLabel("<html><p width=\"380\">The target markers file is optional and has one marker name per row. If there are any additional columns, then they are ignored. See manual for more detail.</p></html>");
      contentPanel.add(tooltipTgtMkrs, "cell 0 2");
    }
    {
      tooltipClusterFilters = Grafik.getToolTipIconLabel("<html><p width=\"380\">The list of cluster filter files is generated from any file ending with \"*clusterFilters.ser\" in the project's data/ directory. The clusters can be manually added from within the ScatterPlot module. See manual for more detail.</p></html>");
      contentPanel.add(tooltipClusterFilters, "cell 0 5");
    }
    {
      tooltipFileroot = Grafik.getToolTipIconLabel("<html><p width=\"380\">The root of the VCF file to be generated. If this already exists, you must click the checkbox to overwrite it.</p></html>");
      contentPanel.add(tooltipFileroot, "cell 0 7");
    }
    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        okButton = new JButton();
        okButton.addActionListener(buttonListener);
        okButton.setText(OK);
        okButton.setActionCommand(OK);
        buttonPane.add(okButton);
        getRootPane().setDefaultButton(okButton);
      }
      {
        cancelButton = new JButton();
        cancelButton.addActionListener(buttonListener);
        cancelButton.setText(CANCEL);
        cancelButton.setActionCommand(CANCEL);
        buttonPane.add(cancelButton);
      }
    }
    addWindowListener(new WindowAdapter() {

      @Override
      public void windowClosing(WindowEvent e) {
        close(true);
      }
    });

    pack();
    updateVCFStatus();
  }

  private void refreshStatus(JTextField fileField) {
    String text = fileField.getText();
    if (new File(text).exists()) {
      fileField.setForeground(Color.GREEN.darker());
    } else {
      fileField.setForeground(Color.RED.darker());
    }
  }

  private final ActionListener buttonListener = new ActionListener() {

    @Override
    public void actionPerformed(ActionEvent arg0) {
      if (arg0.getActionCommand().equals(CANCEL)) {
        close(true);
      } else if (arg0.getActionCommand().equals(OVERWRITE)) {
        if (!isValidVCFRoot()) {
          okButton.setEnabled(chckbxOverwrite.isSelected());
        }
      } else if (arg0.getActionCommand().equals(OK)) {
        if (checkSampleList()) {
          // Note that the actual export happens in org.genvisis.cnv.Launch
          close(false);
        } else {
          JOptionPane.showMessageDialog(VCFExportOptions.this,
                                        "Error - Sample list file doesn't exist.",
                                        "Missing Sample List File", JOptionPane.ERROR_MESSAGE,
                                        null);
        }
      } else if (arg0.getActionCommand().equals(SAMP)) {
        selectFile(textFieldSampleListFile);
      }
    }
  };

  private void selectFile(JTextField field) {
    String value = field.getText();
    JFileChooser fileChooser = new JFileChooser();
    fileChooser.setMultiSelectionEnabled(false);
    if (value != null && !"".equals(value)) {
      if (Files.exists(value)) {
        fileChooser.setSelectedFile(new File(value));
      } else if (Files.exists(ext.parseDirectoryOfFile(value))) {
        fileChooser.setSelectedFile(new File(ext.parseDirectoryOfFile(value)));
      }
    }
    fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
    if (fileChooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
      field.setText(fileChooser.getSelectedFile().getAbsolutePath());
    }
  }

  private JLabel lblSampleFile;
  private JTextField textFieldSampleListFile;
  private JLabel tooltipPedigree;
  private JLabel tooltipTgtMkrs;
  private JLabel tooltipClusterFilters;
  private JLabel tooltipExportType;
  private JLabel tooltipFileroot;
  private JButton btnSampleFile;
  private JPanel panel;
  private JSeparator separator;
  private JCheckBox chckbxSplitChromosomes;
  private JCheckBox chckbxExportchrIn;
  private JButton btnSelectChromosomesTo;
  private Byte[] chrsToWrite;

  private void close(boolean cancelled) {
    this.cancelled = cancelled;
    setVisible(false);
  }

  protected boolean checkSampleList() {
    String text = textFieldSampleListFile.getText();
    return text.equals("") || new File(text).exists();
  }

  private void updateVCFStatus() {
    if (isValidVCFRoot()) {
      lblNameConflict.setVisible(false);
      chckbxOverwrite.setEnabled(false);
      chckbxOverwrite.setSelected(false);
      okButton.setEnabled(true);
    } else {
      lblNameConflict.setVisible(true);
      chckbxOverwrite.setEnabled(true);
      okButton.setEnabled(chckbxOverwrite.isSelected());
    }
  }

  private boolean isValidVCFRoot() {
    String vcfRoot = textFieldVCFFileroot.getText().trim();
    String vcfRootDir = (Files.isRelativePath(vcfRoot) ? proj.PROJECT_DIRECTORY.getValue() : "")
                        + vcfRoot;
    boolean exists = false;
    if (new File(vcfRootDir + ".vcf.gz").exists()) {
      exists = true;
    }
    return !exists;
  }

  private String[] getClusterFiltersOptions() {
    return ArrayUtils.addStrToArray(NO_CLUSTER_FILTERS,
                                    Files.list(proj.DATA_DIRECTORY.getValue(false, true), null,
                                               ext.removeDirectoryInfo(proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME)),
                                               false));
  }

  private String[] getTargetMarkersOptions() {
    String[] values = proj.TARGET_MARKERS_FILENAMES.getValue();
    String[] retVals = new String[values.length + 3];
    for (int i = 1; i <= values.length; i++) {
      retVals[i] = values[i - 1];
    }
    retVals[0] = ALL_MARKERS;
    retVals[retVals.length - 1] = LOAD_MARKERS_LIST;
    retVals[retVals.length - 2] = NEW_MARKERS_LIST;
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
    return isValidVCFRoot() || chckbxOverwrite.isSelected();
  }

  public String getVCFRoot() {
    return getShouldWrite() ? textFieldVCFFileroot.getText().trim() : null;
  }

  public double getGC() {
    return (Double) gcSpinner.getValue();
  }

  public boolean getSplitChrs() {
    return chckbxSplitChromosomes.isSelected();
  }

  public boolean getWriteChrContigs() {
    return chckbxExportchrIn.isSelected();
  }

  public Byte[] getChrsToWrite() {
    return chrsToWrite;
  }

  public boolean getCancelled() {
    return cancelled;
  }

  public String getSampleListFile() {
    return textFieldSampleListFile.getText();
  }

  public String getTargetMarkersFile() {
    String val = (String) comboBoxTargetMarkers.getSelectedItem();
    if (ALL_MARKERS.equals(val)) {
      return null;
    } else if (NEW_MARKERS_LIST.equals(val)) {
      return null;
    } else if (LOAD_MARKERS_LIST.equals(val)) {
      return null;
    } else {
      return val;
    }
  }

}
