package org.genvisis.cnv.gui;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import org.genvisis.common.Files;
import org.genvisis.common.Positions;

import net.miginfocom.swing.MigLayout;

public class NewRegionListDialog extends JDialog implements ActionListener {

  private static final long serialVersionUID = 1L;
  private static final String DEFAULT_MISSING_UCSC_LOCATION = "chr1";

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    EventQueue.invokeLater(new Runnable() {
      @Override
      public void run() {
        try {
          NewRegionListDialog frame =
              new NewRegionListDialog(new String[0], null, true, DEFAULT_MISSING_UCSC_LOCATION);
          frame.setVisible(true);
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    });
  }

  private final JPanel contentPane;
  private final JTextField textField;
  private final JTextArea textArea;
  private final JButton btnCreate;
  private final JButton btnCancel;
  private final JButton button;
  private int returnCode = JOptionPane.DEFAULT_OPTION;
  private HashSet<String> idSet = null;
  private String missingValue = DEFAULT_MISSING_UCSC_LOCATION;
  private boolean allowMissing = true;
  private String dir = "";

  public NewRegionListDialog(String[] sampleNames, final String dir, boolean allowMissing) {
    this(sampleNames, dir, allowMissing, DEFAULT_MISSING_UCSC_LOCATION);
  }

  /**
   * Create the frame.
   */
  public NewRegionListDialog(String[] sampleNames, final String dir, boolean allowMissing,
      String missingValue) {
    this.allowMissing = allowMissing;
    this.missingValue =
        allowMissing && missingValue == null ? DEFAULT_MISSING_UCSC_LOCATION : missingValue;
    this.dir = dir;
    setTitle("Create New UCSC Regions List");
    setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
    setBounds(100, 100, 550, 300);
    contentPane = new JPanel();
    contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    setContentPane(contentPane);
    contentPane.setLayout(new MigLayout("", "[424px,grow]", "[14px][238px][][][][][grow]"));

    JScrollPane scrollPane = new JScrollPane();
    contentPane.add(scrollPane, "cell 0 1,grow");

    textArea = new JTextArea();
    scrollPane.setViewportView(textArea);

    String label = sampleNames == null
        ? "UCSC regions (one per line, with optional tab-separated comments): "
        : "Sample IDs with (optional) UCSC regions and comments (tab-separated, one set per line): ";
    JLabel lblMarkerNamesone = new JLabel(label);
    contentPane.add(lblMarkerNamesone, "cell 0 0,growx,aligny top");

    JLabel lblFileName = new JLabel("File name (relative to project directory):");
    contentPane.add(lblFileName, "cell 0 2");

    textField = new JTextField();
    contentPane.add(textField, "cell 0 3,growx");
    textField.setColumns(10);

    JSeparator separator = new JSeparator();
    contentPane.add(separator, "cell 0 5,growx");

    JPanel panel = new JPanel();
    contentPane.add(panel, "south,alignx right,growy");
    panel.setLayout(new MigLayout("alignx right, ins 5 0 3 0", "[65px][65px]", "[23px]"));

    btnCreate = new JButton();
    btnCreate.addActionListener(this);
    btnCreate.setText("Create");
    panel.add(btnCreate, "cell 0 0,alignx left,aligny top");

    btnCancel = new JButton();
    btnCancel.addActionListener(this);
    btnCancel.setText("Cancel");
    panel.add(btnCancel, "cell 1 0,alignx left,aligny top");

    button = new JButton();
    Action fileSelectAction = new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent arg0) {
        JFileChooser jfc = new JFileChooser(dir);
        jfc.setMultiSelectionEnabled(false);
        int opt = jfc.showSaveDialog(NewRegionListDialog.this);
        if (opt == JFileChooser.APPROVE_OPTION) {
          textField.setText(jfc.getSelectedFile().getAbsolutePath());
        }
      }
    };
    button.setAction(fileSelectAction);
    button.setText("...");
    contentPane.add(button, "cell 0 3");

    if (sampleNames != null) {
      idSet = new HashSet<String>();
      for (String id : sampleNames) {
        idSet.add(id);
      }
    }
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    if (e.getSource() == btnCreate) {
      boolean success = doCreate();
      if (!success) {
        return;
      } else {
        returnCode = JOptionPane.YES_OPTION;
        setVisible(false);
      }
    } else {
      returnCode = JOptionPane.NO_OPTION;
      setVisible(false);
    }
  }

  private boolean doCreate() {
    int sampInd = idSet == null ? -1 : 0;
    int rgnIndex = idSet == null ? 0 : 1;

    String[][] rgnsWithComments = getRegions();
    ArrayList<String> invalidPositions = new ArrayList<String>();
    ArrayList<String> invalidIDs = new ArrayList<String>();
    for (String[] rgn : rgnsWithComments) {
      int[] pos;
      if (rgn.length == 1 && sampInd != -1 && allowMissing) {
        pos = Positions.parseUCSClocation(missingValue);
      } else {
        pos = Positions.parseUCSClocation(rgn[rgnIndex]);
      }
      if (pos == null) {
        invalidPositions.add(rgn[rgnIndex]);
      }
      if (sampInd != -1) {
        if (!idSet.contains(rgn[sampInd])) {
          invalidIDs.add(rgn[sampInd]);
        }
      }
    }
    if (invalidPositions.size() > 0) {
      String[] options = {"Ignore and Continue", "Return"};
      StringBuilder msg = new StringBuilder("Warning - ").append(invalidPositions.size())
          .append(" regions are not valid UCSC regions:");
      for (String inv : invalidPositions) {
        msg.append("\n").append(inv);
      }
      int opt = JOptionPane.showOptionDialog(this, msg.toString(), "Warning - invalid regions!",
          JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE, null, options, options[1]);
      if (opt != 0) {
        return false;
      }
    }
    if (invalidIDs.size() > 0) {
      String[] options = {"Ignore and Continue", "Return"};
      StringBuilder msg = new StringBuilder("Warning - ").append(invalidIDs.size())
          .append(" sample IDs are not valid:");
      for (String inv : invalidIDs) {
        msg.append("\n").append(inv);
      }
      int opt = JOptionPane.showOptionDialog(this, msg.toString(), "Warning - invalid sample IDs!",
          JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE, null, options, options[1]);
      if (opt != 0) {
        return false;
      }
    }

    String filepath = textField.getText();
    if (Files.isRelativePath(filepath)) {
      filepath = dir + filepath;
    }
    File newFile = new File(filepath);
    boolean reachable = false;
    try {
      filepath = newFile.getCanonicalPath();
      reachable = true;
    } catch (IOException e) {
      reachable = false;
    }
    if (!reachable) {
      JOptionPane.showMessageDialog(null, "Error - file name [" + filepath + "] is invalid.",
          "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    } else if (Files.exists(filepath)) {
      JOptionPane.showMessageDialog(null, "Error - file [" + filepath + "] already exists.",
          "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }

    Files.writeMatrix(rgnsWithComments, filepath, "\t");
    textField.setText(filepath);
    return true;
  }

  public String getFileName() {
    return textField.getText();
  }

  public String[][] getRegions() {
    String[] rgnLines = textArea.getText().split("\n");
    String[][] rgnsWithComments = new String[rgnLines.length][];
    for (int i = 0; i < rgnLines.length; i++) {
      rgnsWithComments[i] = rgnLines[i].trim().split("\t");
    }
    return rgnsWithComments;
  }

  public int getReturnCode() {
    return returnCode;
  }



}
