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

import net.miginfocom.swing.MigLayout;

public class NewMarkerListDialog extends JDialog implements ActionListener {

  private static final long serialVersionUID = 1L;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    EventQueue.invokeLater(new Runnable() {
      @Override
      public void run() {
        try {
          NewMarkerListDialog frame = new NewMarkerListDialog(new String[0], null);
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
  private final HashSet<String> markerSet;
  private final JButton btnCreate;
  private final JButton btnCancel;
  private int returnCode = JOptionPane.DEFAULT_OPTION;


  private final JButton button;

  /**
   * Create the frame.
   */
  public NewMarkerListDialog(String[] markers, final String dir) {
    setTitle("Create New Marker List");
    setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
    setBounds(100, 100, 450, 300);
    contentPane = new JPanel();
    contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    setContentPane(contentPane);
    contentPane.setLayout(new MigLayout("", "[424px,grow]", "[14px][238px][][][][][grow]"));

    JScrollPane scrollPane = new JScrollPane();
    contentPane.add(scrollPane, "cell 0 1,grow");

    textArea = new JTextArea();
    scrollPane.setViewportView(textArea);

    JLabel lblMarkerNamesone =
        new JLabel("Marker names (one per line, with tab-separated comments): ");
    contentPane.add(lblMarkerNamesone, "cell 0 0,growx,aligny top");

    JLabel lblFileName = new JLabel("File name (full path):");
    contentPane.add(lblFileName, "cell 0 2");

    textField = new JTextField();
    contentPane.add(textField, "flowx,cell 0 3,growx");
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
        int opt = jfc.showSaveDialog(NewMarkerListDialog.this);
        if (opt == JFileChooser.APPROVE_OPTION) {
          textField.setText(jfc.getSelectedFile().getAbsolutePath());
        }
      }
    };
    button.setAction(fileSelectAction);
    button.setText("...");
    contentPane.add(button, "cell 0 3");

    markerSet = new HashSet<String>();
    for (String marker : markers) {
      markerSet.add(marker);
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
    String[][] mkrsWithComments = getMarkers();
    ArrayList<String> invalid = new ArrayList<String>();
    for (String[] mkr : mkrsWithComments) {
      if (!markerSet.contains(mkr[0])) {
        invalid.add(mkr[0]);
      }
    }
    if (invalid.size() > 0) {
      String[] options = {"Ignore and Continue", "Return"};
      StringBuilder msg =
          new StringBuilder("Warning - ").append(invalid.size())
                                         .append(" markers are not present in the current marker set:");
      for (String inv : invalid) {
        msg.append("\n").append(inv);
      }
      int opt = JOptionPane.showOptionDialog(this, msg.toString(), "Warning - invalid markers!",
                                             JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE,
                                             null, options, options[1]);
      if (opt != 0) {
        return false;
      }
    }

    String filepath = textField.getText();
    String dir = "";
    File newFile = new File(dir + filepath);
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
    } else if (Files.exists(dir + filepath)) {
      JOptionPane.showMessageDialog(null, "Error - file [" + filepath + "] already exists.",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }

    Files.writeMatrix(mkrsWithComments, filepath, "\t");
    textField.setText(filepath);
    return true;
  }

  public String getFileName() {
    return textField.getText();
  }

  public String[][] getMarkers() {
    String[] mkrLines = textArea.getText().split("\n");
    String[][] mkrsWithComments = new String[mkrLines.length][];
    for (int i = 0; i < mkrLines.length; i++) {
      mkrsWithComments[i] = mkrLines[i].trim().split("\t");
    }
    return mkrsWithComments;
  }

  public int getReturnCode() {
    return returnCode;
  }



}
