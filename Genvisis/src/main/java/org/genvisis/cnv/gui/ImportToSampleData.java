package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import javax.swing.DefaultListCellRenderer;
import javax.swing.DefaultListModel;
import javax.swing.DefaultListSelectionModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

import com.google.common.collect.Sets;

import net.miginfocom.swing.MigLayout;

public class ImportToSampleData extends JDialog {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private Project proj;
  private Logger log;

  private final JPanel contentPanel = new JPanel();
  private JTextField fldFile;
  private JTextField fldMissingValue;
  private JCheckBox chkAllowMissing;
  private JLabel lblMissingValue;

  private ActionListener listMissing = new ActionListener() {
    public void actionPerformed(ActionEvent e) {
      lblMissingValue.setEnabled(chkAllowMissing.isSelected());
      fldMissingValue.setEnabled(chkAllowMissing.isSelected());
    }
  };

  private ActionListener listImport = new ActionListener() {
    public void actionPerformed(ActionEvent e) {
      String file = fldFile.getText();
      int[] cols = list.getSelectedIndices();
      if (cols.length == 0) return;
      String[] headers = new String[cols.length];
      for (int i = 0; i < cols.length; i++) {
        headers[i] = list.getModel().getElementAt(cols[i]);
        cols[i] = cols[i] + 1;
      }
      String missing = fldMissingValue.getText().trim();
      Hashtable<String, String> hash = HashVec.loadFileToHashString(file, 0, cols, "\t", true);
      Set<String> allSamples = Sets.newHashSet(proj.getSamples());
      Set<String> foundSamples = new HashSet<>(hash.keySet());
      foundSamples.removeAll(allSamples);
      allSamples.removeAll(hash.keySet());
      if (foundSamples.size() > 0) {
        JOptionPane.showMessageDialog(ImportToSampleData.this,
                                      "<html>Error, invalid samples found:<br />"
                                                               + ArrayUtils.toStr(foundSamples,
                                                                                  "<br />")
                                                               + "</html>",
                                      "Error - Extra Samples", JOptionPane.ERROR_MESSAGE);
        return;
      }
      if (allSamples.size() > 0) {
        if (!chkAllowMissing.isSelected()) {
          JOptionPane.showMessageDialog(ImportToSampleData.this,
                                        "<html>Error, missing samples found:<br />"
                                                                 + ArrayUtils.toStr(allSamples,
                                                                                    "<br />")
                                                                 + "</html>",
                                        "Error - Samples Missing", JOptionPane.ERROR_MESSAGE);
          return;
        }
      }
      proj.getSampleData(false).addData(hash, "DNA", headers, missing, "\t", proj.getLog());
      setVisible(false);
      dispose();
    }
  };

  private ActionListener listSelectFile = new ActionListener() {
    public void actionPerformed(ActionEvent e) {
      JFileChooser chooser = new JFileChooser(proj.PROJECT_DIRECTORY.getValue());
      chooser.setMultiSelectionEnabled(false);

      chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
      chooser.setDialogTitle("Select File");

      int retValue = chooser.showDialog(ImportToSampleData.this, "Select");

      if (retValue == JFileChooser.CANCEL_OPTION) {
        return;
      } else {
        File newFile = chooser.getSelectedFile();
        String file = newFile.getPath();
        String[] header = Files.getHeaderOfFile(file, log);
        int[] inds = ext.indexFactors(Aliases.DNA, header, false);
        inds = ArrayUtils.removeAllValues(inds, -1);
        if (inds.length == 0) {
          JOptionPane.showMessageDialog(ImportToSampleData.this,
                                        "Error - no valid column header found.  Looking for one of "
                                                                 + ArrayUtils.toStr(Aliases.DNA)
                                                                 + ".",
                                        "Error - Invalid Header", JOptionPane.ERROR_MESSAGE);
          return;
        }
        if (inds.length > 1) {
          String[] hdrs = new String[inds.length];
          for (int i = 0; i < inds.length; i++) {
            hdrs[i] = header[inds[i]];
          }
          JOptionPane.showMessageDialog(ImportToSampleData.this,
                                        "Error - multiple valid column headers found:\n "
                                                                 + ArrayUtils.toStr(hdrs) + ".",
                                        "Error - Invalid Header", JOptionPane.ERROR_MESSAGE);
          return;
        }
        int headerInd = inds[0];
        String[] otherHeaders = ArrayUtils.removeFromArray(header, headerInd);

        String[] sampleDataHeader = Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(),
                                                          log);
        preexisting = Sets.newHashSet(sampleDataHeader);

        DefaultListModel<String> newModel = new DefaultListModel<>();
        for (String o : otherHeaders) {
          newModel.addElement(o);
        }

        list.setModel(newModel);
        separator.setVisible(true);
        lblSelectCols.setVisible(true);
        scrollPane.setVisible(true);
        fldFile.setText(newFile.getPath());
        setSize(420, 400);
      }
    }
  };
  private JList<String> list;
  private JSeparator separator;
  private JScrollPane scrollPane;
  private Set<String> preexisting = new HashSet<>();
  private JLabel lblSelectCols;

  public ImportToSampleData(Project proj) {
    this.proj = proj;
    this.log = proj.getLog();
    setModal(true);
    setTitle("Import Data to SampleData");
    setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    setBounds(100, 100, 420, 200);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(new MigLayout("hidemode 3", "[99.00][grow][]", "[][][][][][][grow]"));

    JLabel lblFileToImport = new JLabel("File to Import:");
    contentPanel.add(lblFileToImport, "cell 0 0");

    JLabel lblFileInfo = Grafik.getToolTipIconLabel("File must have a header, and the first column must be one of: "
                                                    + ArrayUtils.toStr(Aliases.DNA, ", ") + ".");
    contentPanel.add(lblFileInfo, "cell 0 0");

    fldFile = new JTextField();
    fldFile.setEditable(false);
    contentPanel.add(fldFile, "cell 0 1 2 1,grow");
    fldFile.setColumns(10);

    JButton btnFile = new JButton(">");
    btnFile.addActionListener(listSelectFile);
    contentPanel.add(btnFile, "cell 2 1,grow");

    JLabel lblAllowMissing = new JLabel("Allow Missing:");
    contentPanel.add(lblAllowMissing, "cell 0 2,alignx center");

    lblMissingValue = new JLabel("Missing Value:");
    lblMissingValue.setEnabled(false);
    contentPanel.add(lblMissingValue, "cell 1 2");

    chkAllowMissing = new JCheckBox("");
    chkAllowMissing.addActionListener(listMissing);
    contentPanel.add(chkAllowMissing, "cell 0 3,alignx center");

    fldMissingValue = new JTextField(".");
    fldMissingValue.setFont(new Font("Monospaced", Font.PLAIN,
                                     fldMissingValue.getFont().getSize()));
    fldMissingValue.setEnabled(false);
    contentPanel.add(fldMissingValue, "cell 1 3 2 1,grow");
    fldMissingValue.setColumns(10);

    separator = new JSeparator();
    separator.setVisible(false);
    contentPanel.add(separator, "cell 0 4 3 1,growx");

    lblSelectCols = new JLabel("Select Columns to Add:");
    lblSelectCols.setVisible(false);
    contentPanel.add(lblSelectCols, "cell 0 5");

    scrollPane = new JScrollPane();
    scrollPane.setVisible(false);
    contentPanel.add(scrollPane, "cell 0 6 3 1,grow");

    list = new JList<>();
    list.setCellRenderer(new DefaultListCellRenderer() {
      private static final long serialVersionUID = 1L;

      @Override
      public Component getListCellRendererComponent(JList<?> list, Object value, int index,
                                                    boolean isSelected, boolean cellHasFocus) {
        Component c = super.getListCellRendererComponent(list, value, index, isSelected,
                                                         cellHasFocus);
        if (preexisting.contains(list.getModel().getElementAt(index))) {
          c.setEnabled(false);
          c.setFocusable(false);
          if (c instanceof JComponent) {
            ((JComponent) c).setToolTipText("Already present in SampleData");
          }
        } else {
          c.setEnabled(true);
          c.setFocusable(true);
          if (c instanceof JComponent) {
            ((JComponent) c).setToolTipText(null);
          }
        }
        return c;
      }
    });
    list.setSelectionModel(new DefaultListSelectionModel() {
      private static final long serialVersionUID = 1L;

      @Override
      public void setSelectionInterval(int index0, int index1) {
        setValueIsAdjusting(true);
        super.setSelectionInterval(index0, index1);
        for (int i = index0; i <= index1; i++) {
          if (preexisting.contains(list.getModel().getElementAt(i))) {
            removeSelectionInterval(i, i);
          }
        }
        setValueIsAdjusting(false);
      }

      @Override
      public void addSelectionInterval(int index0, int index1) {
        setValueIsAdjusting(true);
        super.addSelectionInterval(index0, index1);
        for (int i = index0; i <= index1; i++) {
          if (preexisting.contains(list.getModel().getElementAt(i))) {
            removeSelectionInterval(i, i);
          }
        }
        setValueIsAdjusting(false);
      }
    });
    scrollPane.setViewportView(list);

    JPanel buttonPane = new JPanel();
    buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
    getContentPane().add(buttonPane, BorderLayout.SOUTH);

    JButton okButton = new JButton("Import");
    okButton.addActionListener(listImport);
    okButton.setActionCommand("Import");
    buttonPane.add(okButton);
    getRootPane().setDefaultButton(okButton);

    JButton cancelButton = new JButton("Cancel");
    cancelButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
      }
    });
    cancelButton.setActionCommand("Cancel");
    buttonPane.add(cancelButton);
  }

}
