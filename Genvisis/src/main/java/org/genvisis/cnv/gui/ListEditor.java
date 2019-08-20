package org.genvisis.cnv.gui;

import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.gui.UITools;

import net.miginfocom.swing.MigLayout;

public class ListEditor extends JDialog implements ActionListener {

  private static final long serialVersionUID = 1L;
  private final JPanel contentPane;
  private final JTextField textField;
  private final JTextArea textArea;
  private final JButton btnCreate;
  private final JButton btnCancel;
  private int returnCode = JOptionPane.NO_OPTION;
  private final JButton button;

  Validator validator = new TextListValidator();

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    EventQueue.invokeLater(new Runnable() {

      @Override
      public void run() {
        try {
          String title = "Marker List";
          String desc = "Marker names (one per line, with tab-separated comments): ";
          ListEditor frame = new ListEditor(title, desc, "./", true,
                                            "D:/data/ny_choanal/omni2.5v1.2/data/test.txt",
                                            new MarkerValidator(null));
          frame.setVisible(true);
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    });
  }

  public static ListEditor createMarkerListCreator(Project proj) {
    String title = "Marker List";
    String listDesc = "Marker names (one per line, with tab-separated comments): ";
    String defaultDir = proj.PROJECT_DIRECTORY.getValue();
    String[] allowedValues = proj.getMarkerNames();
    boolean isEditor = false;
    String fileToEdit = null;
    return new ListEditor(title, listDesc, defaultDir, isEditor, fileToEdit,
                          new MarkerValidator(allowedValues));
  }

  public static ListEditor createMarkerListEditor(Project proj, String markerFile) {
    String title = "Marker List";
    String listDesc = "Marker names (one per line, with tab-separated comments): ";
    String defaultDir = proj.PROJECT_DIRECTORY.getValue();
    String[] allowedValues = proj.getMarkerNames();
    boolean isEditor = true;
    String fileToEdit = markerFile;
    return new ListEditor(title, listDesc, defaultDir, isEditor, fileToEdit,
                          new MarkerValidator(allowedValues));
  }

  public static ListEditor createRegionListCreator(String[] sampleNames, String dir,
                                                   boolean allowMissing) {
    String title = "UCSC Regions List";
    String listDesc = sampleNames == null ? "UCSC regions (one per line, with optional tab-separated comments): "
                                          : "Sample IDs with (optional) UCSC regions and comments (tab-separated, one set per line): ";
    String defaultDir = dir;
    boolean isEditor = false;
    String fileToEdit = null;
    return new ListEditor(title, listDesc, defaultDir, isEditor, fileToEdit,
                          new RegionValidator(sampleNames,
                                              RegionValidator.DEFAULT_MISSING_UCSC_LOCATION,
                                              allowMissing));
  }

  public static ListEditor createRegionListEditor(String[] sampleNames, String dir,
                                                  boolean allowMissing, String fileToEdit) {
    String title = "UCSC Regions List";
    String listDesc = sampleNames == null ? "UCSC regions (one per line, with optional tab-separated comments): "
                                          : "Sample IDs with (optional) UCSC regions and comments (tab-separated, one set per line): ";
    String defaultDir = dir;
    boolean isEditor = true;
    return new ListEditor(title, listDesc, defaultDir, isEditor, fileToEdit,
                          new RegionValidator(sampleNames,
                                              RegionValidator.DEFAULT_MISSING_UCSC_LOCATION,
                                              allowMissing));
  }

  /**
   * Create the frame.
   */
  public ListEditor(String title, String listDesc, String defaultDir, boolean isEditor,
                    String fileToEdit, Validator validator) {
    setTitle((isEditor ? "Edit " : "Create New ") + title);
    this.validator = validator == null ? new TextListValidator() : validator;

    setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
    setMinimumSize(new Dimension(100, 100));
    UITools.setSize(this, new Dimension(550, 300));
    contentPane = new JPanel();
    contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    setContentPane(contentPane);
    contentPane.setLayout(new MigLayout("", "[550px,grow]", "[14px][238px][][][][][grow]"));

    JScrollPane scrollPane = new JScrollPane();
    contentPane.add(scrollPane, "cell 0 1,grow");

    textArea = new JTextArea();
    scrollPane.setViewportView(textArea);

    JLabel lblMarkerNamesone = new JLabel(listDesc);
    contentPane.add(lblMarkerNamesone, "cell 0 0,growx,aligny top");

    // TODO change to "(relative to project directory)" if project exists:
    JLabel lblFileName = new JLabel("File name (full path):");
    contentPane.add(lblFileName, "cell 0 2");

    textField = new JTextField();
    contentPane.add(textField, "flowx,cell 0 3,growx");
    textField.setColumns(10);

    if (isEditor) {
      if (fileToEdit != null && !"".equals(fileToEdit)) {
        if (Files.exists(fileToEdit)) {
          String[] text = HashVec.loadFileToStringArray(fileToEdit, false, null, false);
          textArea.setText(ArrayUtils.toStr(text, "\n"));

          textField.setText(fileToEdit);
        }
      }
    }

    JSeparator separator = new JSeparator();
    contentPane.add(separator, "cell 0 5,growx");

    JPanel panel = new JPanel();
    contentPane.add(panel, "south,alignx right,growy");
    panel.setLayout(new MigLayout("alignx right, ins 5 0 3 0", "[65px][65px]", "[23px]"));

    btnCreate = new JButton();
    btnCreate.addActionListener(this);
    btnCreate.setText(isEditor ? "Save" : "Create");
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
        JFileChooser jfc = new JFileChooser(defaultDir);
        jfc.setMultiSelectionEnabled(false);
        int opt = jfc.showSaveDialog(ListEditor.this);
        if (opt == JFileChooser.APPROVE_OPTION) {
          textField.setText(jfc.getSelectedFile().getAbsolutePath());
        }
      }
    };
    button.setAction(fileSelectAction);
    button.setText("...");
    contentPane.add(button, "cell 0 3");

    pack();
  }

  public String getFileName() {
    return textField.getText();
  }

  /**
   * @return Either JOptionPane.YES_OPTION or JOptionPane.NO_OPTION;
   */
  public int getReturnCode() {
    return returnCode;
  }

  public String[][] getValues() {
    String[] lines = textArea.getText().split("\n");
    String[][] splitLines = new String[lines.length][];
    for (int i = 0; i < lines.length; i++) {
      splitLines[i] = lines[i].trim().split("\t", -1);
    }
    return splitLines;
  }

  static abstract class Validator {

    abstract boolean validate(String[] line);

    abstract boolean reportInvalid(List<String[]> invalidLines);

    public boolean validate(String[][] lines) {
      ArrayList<String[]> invalidLines = new ArrayList<>();
      for (String[] line : lines) {
        boolean success = validate(line);
        if (!success) {
          invalidLines.add(line);
        }
      }
      return reportInvalid(invalidLines);
    }

  }

  static class TextListValidator extends Validator {

    public TextListValidator() {}

    public boolean validate() {
      return true;
    }

    @Override
    boolean reportInvalid(List<String[]> invalidLines) {
      return true;
    }

    @Override
    boolean validate(String[] line) {
      return true;
    }
  }

  public static class MarkerValidator extends Validator {

    Set<String> allowedMkrs;

    public MarkerValidator(String[] allowedMarkers) {
      if (allowedMarkers == null) {
        allowedMkrs = null;
      } else {
        allowedMkrs = new HashSet<>();
        for (String s : allowedMarkers) {
          this.allowedMkrs.add(s);
        }
      }
    }

    @Override
    boolean validate(String[] line) {
      return allowedMkrs == null || allowedMkrs.contains(line[0]);
    }

    @Override
    boolean reportInvalid(List<String[]> invalidLines) {
      if (invalidLines.size() > 0) {
        String[] options = {"Ignore and Continue", "Return"};
        StringBuilder msg = new StringBuilder("Warning - ").append(invalidLines.size())
                                                           .append(" markers were not present in the list of approved markers:");
        for (String[] inv : invalidLines) {
          msg.append("\n").append(inv[0]);
        }
        int opt = JOptionPane.showOptionDialog(null, msg.toString(), "Warning - invalid values!",
                                               JOptionPane.YES_NO_OPTION,
                                               JOptionPane.WARNING_MESSAGE, null, options,
                                               options[1]);
        if (opt != 0) {
          return false;
        }
      }
      return true;
    }

  }

  public static class RegionValidator extends Validator {

    private static final String DEFAULT_MISSING_UCSC_LOCATION = "chr1";
    private String missingValue = DEFAULT_MISSING_UCSC_LOCATION;
    private boolean allowMissing = false;

    Set<String> samples;

    ArrayList<String> invalidPositions = new ArrayList<>();
    ArrayList<String> invalidIDs = new ArrayList<>();

    public RegionValidator(String[] sampleNames, String missingVal, boolean allowMissing) {
      if (missingVal != null) {
        this.missingValue = missingVal;
      }
      this.allowMissing = allowMissing;
      if (sampleNames == null) {
        samples = null;
      } else {
        samples = new HashSet<>();
        for (String s : sampleNames) {
          samples.add(s);
        }
      }
    }

    @Override
    boolean validate(String[] line) {
      int sampInd = samples == null ? -1 : 0;
      int rgnIndex = samples == null ? 0 : 1;

      int[] pos;
      if (line.length == 1 && sampInd != -1 && allowMissing) {
        pos = Positions.parseUCSClocation(missingValue);
      } else {
        pos = Positions.parseUCSClocation(line[rgnIndex]);
      }
      if (pos == null) {
        invalidPositions.add(line[rgnIndex]);
      }
      if (sampInd != -1) {
        if (samples != null && !samples.contains(line[sampInd])) {
          invalidIDs.add(line[sampInd]);
        }
      }
      return true;
    }

    @Override
    boolean reportInvalid(List<String[]> invalidLines) {
      boolean success = true;

      if (invalidPositions.size() > 0) {
        String[] options = {"Ignore and Continue", "Return"};
        StringBuilder msg = new StringBuilder("Warning - ").append(invalidPositions.size())
                                                           .append(" regions are not valid UCSC regions:");
        for (String inv : invalidPositions) {
          msg.append("\n").append(inv);
        }
        int opt = JOptionPane.showOptionDialog(null, msg.toString(), "Warning - invalid regions!",
                                               JOptionPane.YES_NO_OPTION,
                                               JOptionPane.WARNING_MESSAGE, null, options,
                                               options[1]);
        if (opt != 0) {
          success = false;
        }
      }
      if (success) {
        if (invalidIDs.size() > 0) {
          String[] options = {"Ignore and Continue", "Return"};
          StringBuilder msg = new StringBuilder("Warning - ").append(invalidIDs.size())
                                                             .append(" sample IDs are not valid:");
          for (String inv : invalidIDs) {
            msg.append("\n").append(inv);
          }
          int opt = JOptionPane.showOptionDialog(null, msg.toString(),
                                                 "Warning - invalid sample IDs!",
                                                 JOptionPane.YES_NO_OPTION,
                                                 JOptionPane.WARNING_MESSAGE, null, options,
                                                 options[1]);
          if (opt != 0) {
            success = false;
          }
        }
      }

      invalidPositions = new ArrayList<>();
      invalidIDs = new ArrayList<>();

      return success;
    }
  }

  private boolean doCreate() {
    String[][] values = getValues();
    if (validator.validate(values)) {
      String filepath = textField.getText();
      String dir = "";
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
        String msg = "Error - file [" + filepath
                     + "] already exists.\n\nWould you like to overwrite this file?\n";
        String ttl = "Error - File Already Exists!";
        String[] opts = {"Overwrite", "Cancel"};

        int val = JOptionPane.showOptionDialog(ListEditor.this, msg, ttl,
                                               JOptionPane.OK_CANCEL_OPTION,
                                               JOptionPane.WARNING_MESSAGE, null, opts, opts[1]);
        if (val == JOptionPane.CLOSED_OPTION || val == 1) {
          return false;
        }
      }

      Files.writeMatrix(values, filepath, "\t");
      textField.setText(filepath);
      return true;
    }
    return false;
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

}
