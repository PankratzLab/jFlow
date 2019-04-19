package org.genvisis.one.ben.imagetag;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.border.EmptyBorder;

import org.genvisis.one.ben.imagetag.AnnotatedImage.Annotation;
import org.pankratzlab.common.Files;

import net.miginfocom.swing.MigLayout;

public class AnnotationExportDialog extends JDialog {

  private final JPanel contentPanel = new JPanel();
  private JList<String> annotationList;
  private JRadioButton exportSamples;
  private JRadioButton exportImgFiles;
  private JTextField fileField;
  private JButton selectOutFile;

  private Action selectFileAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      JFileChooser jfc = new JFileChooser();
      jfc.setMultiSelectionEnabled(false);
      int code = jfc.showSaveDialog(AnnotationExportDialog.this);
      if (code == JFileChooser.APPROVE_OPTION) {
        fileField.setText(jfc.getSelectedFile().getAbsolutePath());
      }
    }
  };
  private JButton cancelButton;
  private JButton okButton;
  private IAnnotator annotator;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      AnnotationExportDialog dialog = new AnnotationExportDialog(new Annotator());
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public AnnotationExportDialog(IAnnotator annotator) {
    this.annotator = annotator;
    setBounds(100, 100, 500, 400);
    setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    setTitle("Export Annotations");
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setLayout(new MigLayout("", "[][grow][]", ""));
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);

    contentPanel.add(new JLabel("Select Annotations:"), "cell 1 1");

    DefaultListModel<String> dlm = new DefaultListModel<>();
    for (Annotation a : annotator.getAnnotations()) {
      dlm.addElement(a.annotation);
    }
    annotationList = new JList<>(dlm);
    annotationList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    contentPanel.add(new JScrollPane(annotationList), "cell 1 2, grow");

    contentPanel.add(new JLabel("Select Output Type:"), "cell 1 3");

    exportSamples = new JRadioButton("Samples (export sample if any image has the selected annotation(s))");
    exportImgFiles = new JRadioButton("Images (export specific images with selected annotation(s))");
    ButtonGroup bg1 = new ButtonGroup();
    bg1.add(exportSamples);
    bg1.add(exportImgFiles);

    contentPanel.add(exportSamples, "cell 1 4");
    contentPanel.add(exportImgFiles, "cell 1 5");

    contentPanel.add(new JLabel("Select Output File:"), "cell 1 6");
    fileField = new JTextField();
    contentPanel.add(fileField, "cell 1 7, grow, split 2");
    selectOutFile = new JButton();
    selectOutFile.setAction(selectFileAction);
    selectOutFile.setText(">");
    contentPanel.add(selectOutFile, "cell 1 7, pad 0 0 -1 0");

    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        okButton = new JButton();
        okButton.addActionListener((e) -> {
          AnnotationExportDialog.this.runExport();
        });
        okButton.setText("OK");
        getRootPane().setDefaultButton(okButton);
        buttonPane.add(okButton);
      }
      {
        cancelButton = new JButton();
        cancelButton.addActionListener((e) -> {
          close();
        });
        cancelButton.setText("Cancel");
        buttonPane.add(cancelButton);
      }
    }
  }

  private void runExport() {
    List<String> annots = annotationList.getSelectedValuesList();
    if (annots.isEmpty()) {
      JOptionPane.showMessageDialog(AnnotationExportDialog.this, "Error - no annotations selected.",
                                    "Error!", JOptionPane.ERROR_MESSAGE);
      return;
    }

    String outputFile = fileField.getText();
    if ("".equals(outputFile.trim())) {
      JOptionPane.showMessageDialog(AnnotationExportDialog.this,
                                    "Error - no output file specified.", "Error!",
                                    JOptionPane.ERROR_MESSAGE);
      return;
    }

    List<Annotation> annot = new ArrayList<>();
    for (Annotation a : annotator.getAnnotations()) {
      if (annots.contains(a.annotation)) {
        annot.add(a);
      }
    }

    boolean exportSample = exportSamples.isSelected();

    PrintWriter writer = Files.getAppropriateWriter(outputFile);
    HashMap<String, HashMap<String, AnnotatedImage>> annMap = annotator.getAnnotationMap();
    for (Entry<String, HashMap<String, AnnotatedImage>> fcsAnnot : annMap.entrySet()) {
      sample: for (AnnotatedImage ai : fcsAnnot.getValue().values()) {
        for (Annotation a : ai.getAnnotations()) {
          if (annot.contains(a)) {
            if (exportSample) {
              writer.println(fcsAnnot.getKey());
              break sample;
            } else {
              writer.println(ai.getImageFile());
              continue sample;
            }
          }
        }
      }
    }
    writer.close();

    close();
  }

  private void close() {
    AnnotationExportDialog.this.setVisible(false);
    AnnotationExportDialog.this.dispose();
  }
}
