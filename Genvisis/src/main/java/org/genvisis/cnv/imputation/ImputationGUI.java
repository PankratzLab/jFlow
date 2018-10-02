package org.genvisis.cnv.imputation;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.StringJoiner;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.FileChooser;
import org.genvisis.cnv.gui.JAccordionPanel;
import org.genvisis.cnv.imputation.ImputationPipeline.IMPUTATION_PIPELINE_PATH;
import org.genvisis.cnv.imputation.ImputationPipeline.ImputationPipeRunner;
import org.genvisis.cnv.imputation.ImputationPipeline.KeepDrops;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.utils.qsub.Qsub;
import com.google.common.collect.ImmutableSet;
import net.miginfocom.swing.MigLayout;

public class ImputationGUI extends JDialog {

  private boolean init = true;
  private static final Font COMPONENT_FONT = new Font("Arial", Font.PLAIN, 11);
  private static final Font HEADER_FONT = new Font("Arial", Font.PLAIN, 13);
  private final JPanel contentPanel = new JPanel();
  private final ButtonGroup buttonGroup = new ButtonGroup();
  private JRadioButton rdbtnVcf;
  private JRadioButton rdbtnPlink;
  private JCheckBox chckbxShape;
  private JCheckBox chckbxMini;
  private JCheckBox chckbxUseGRC;
  private JLabel lblUseGRC;
  private ItemListener fileTypeListener = new ItemListener() {

    @Override
    public void itemStateChanged(ItemEvent e) {
      if (init) return;
      if (e.getStateChange() == ItemEvent.SELECTED) {
        boolean plink = e.getItem() == rdbtnPlink;
        chckbxMini.setVisible(plink);
        chckbxShape.setVisible(plink);
        chckbxUseGRC.setVisible(!plink);
        lblUseGRC.setVisible(!plink);
      }
    }
  };
  private JTextField txtFldMkrKeep;
  private JTextField txtFldMkrDrop;
  private JTextField txtFldSampKeep;
  private JTextField txtFldSampDrop;
  private JAccordionPanel keepDropPanel;
  private JAccordionPanel fileTypePanel;
  private JTextField txtFldRefFile;
  private JTextField txtFldOutDir;
  private static final String refFileDesc = "<html><p>A Reference Panel / Site List file, containing at minimum the following columns: mkr, chr, pos, ref, and alt.</p><p>These may be available at http://www.well.ox.ac.uk/~wrayner/tools/</html>";
  private JList<Integer> chrList;
  Project proj;
  Logger log = new Logger();

  /**
   * Create the dialog.
   */
  public ImputationGUI(Project proj) {
    this.proj = proj;
    this.log = proj.getLog();
    setTitle("Genvisis - Imputation Export - " + proj.PROJECT_NAME.getValue());
    setBounds(100, 100, 400, 600);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(new MigLayout("hidemode 3", "[grow]", "[][grow]"));
    setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
    addWindowListener(new WindowAdapter() {

      @Override
      public void windowClosing(WindowEvent e) {
        close();
      }
    });

    {
      JLabel lblImputationExport = new JLabel("<html><u>Imputation Export</u></html>");
      lblImputationExport.setFont(new Font("Arial", Font.PLAIN, 18));
      contentPanel.add(lblImputationExport, "cell 0 0,alignx center");

      JScrollPane scrollPane = new JScrollPane();
      scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
      contentPanel.add(scrollPane, "cell 0 1,grow");
      JPanel panel = new JPanel();
      scrollPane.setViewportView(panel);
      panel.setLayout(new MigLayout("hidemode 3, ins 0 5 0 5", "[grow][grow]",
                                    "[][][][][][][][][grow]"));
      {
        fileTypePanel = new JAccordionPanel();
        fileTypePanel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
        panel.add(fileTypePanel, "cell 0 1 2 1,grow");
      }
      {
        rdbtnPlink = new JRadioButton("PLINK");
        rdbtnPlink.addItemListener(fileTypeListener);
        fileTypePanel.contentPanel.setLayout(new MigLayout("hidemode 3", "[grow][grow]",
                                                           "[][][][][]"));
        rdbtnPlink.setFont(COMPONENT_FONT);
        rdbtnPlink.setSelected(true);
        buttonGroup.add(rdbtnPlink);
        fileTypePanel.contentPanel.add(rdbtnPlink, "cell 0 0,alignx right,aligny top");
      }
      {
        rdbtnVcf = new JRadioButton("VCF");
        rdbtnVcf.addItemListener(fileTypeListener);
        rdbtnVcf.setFont(COMPONENT_FONT);
        buttonGroup.add(rdbtnVcf);
        fileTypePanel.contentPanel.add(rdbtnVcf, "cell 1 0,alignx left,aligny top");
      }

      JLabel lblOptions = new JLabel("<html><u>Options:</u></html>");
      lblOptions.setFont(COMPONENT_FONT);
      fileTypePanel.contentPanel.add(lblOptions, "cell 0 1 2 1,alignx center");
      {
        chckbxShape = new JCheckBox("Generate ShapeIT Script");
        chckbxShape.addItemListener(new ItemListener() {

          @Override
          public void itemStateChanged(ItemEvent e) {
            if (e.getStateChange() == ItemEvent.DESELECTED) {
              if (chckbxMini.isSelected()) {
                chckbxMini.setSelected(false);
              }
            }
          }
        });
        chckbxShape.setFont(COMPONENT_FONT);
        fileTypePanel.contentPanel.add(chckbxShape, "cell 0 2 2 1,alignx center,aligny top");
      }
      {
        chckbxMini = new JCheckBox("Generate MiniMac Script (requires ShapeIT)");
        chckbxMini.addItemListener(new ItemListener() {

          @Override
          public void itemStateChanged(ItemEvent e) {
            if (e.getStateChange() == ItemEvent.SELECTED) {
              if (!chckbxShape.isSelected()) {
                chckbxShape.setSelected(true);
              }
            }
          }
        });
        chckbxMini.setFont(COMPONENT_FONT);
        fileTypePanel.contentPanel.add(chckbxMini, "cell 0 3 2 1,alignx center,aligny top");
      }
      chckbxUseGRC = new JCheckBox("Use GRC instead of HG");
      chckbxUseGRC.setFont(COMPONENT_FONT);
      fileTypePanel.contentPanel.add(chckbxUseGRC, "cell 0 2 2 1,alignx center,aligny top");
      chckbxUseGRC.setVisible(false);
      lblUseGRC = new JLabel("(GRC contig names do not include 'chr')");
      lblUseGRC.setFont(COMPONENT_FONT.deriveFont(10f));
      fileTypePanel.contentPanel.add(lblUseGRC, "cell 0 3 2 1,alignx center,aligny center");
      lblUseGRC.setVisible(false);
      JLabel lblFileType = new JLabel("<html>&nbsp;-&nbsp;<u>File Type:</u></html>");
      fileTypePanel.topPanel.add(lblFileType, "cell 0 0");
      lblFileType.setFont(HEADER_FONT);
      keepDropPanel = new JAccordionPanel();
      panel.add(keepDropPanel, "cell 0 2 2 1,growx");
      keepDropPanel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
      keepDropPanel.contentPanel.setLayout(new MigLayout("", "[grow][grow]", "[][][][][][][][]"));
      keepDropPanel.shrink();
      JLabel lblMarkersToKeep = new JLabel("Markers to Keep File:   (optional)");
      lblMarkersToKeep.setFont(COMPONENT_FONT);
      keepDropPanel.contentPanel.add(lblMarkersToKeep, "cell 0 0 3 1");
      JLabel lblMarkersToDrop = new JLabel("Markers to Drop File:   (optional)");
      lblMarkersToDrop.setFont(COMPONENT_FONT);
      keepDropPanel.contentPanel.add(lblMarkersToDrop, "cell 0 2 3 1");
      txtFldMkrDrop = new JTextField();
      keepDropPanel.contentPanel.add(txtFldMkrDrop, "flowx,cell 0 3 3 1,growx");
      txtFldMkrDrop.setColumns(10);
      JButton btnMkrDrop = new JButton(">");
      btnMkrDrop.addActionListener(e -> {
        String[] files = new FileChooser(ImputationGUI.this, proj.PROJECT_DIRECTORY.getValue(),
                                         false, false, "Select Marker Drops File", log).getFiles();
        if (files == null) {
          return;
        } else {
          txtFldMkrDrop.setText(files[0]);
        }
      });
      keepDropPanel.contentPanel.add(btnMkrDrop, "pad 2 0 -2 0,cell 0 3 3 1,alignx right");
      JLabel lblSamplesToKeep = new JLabel("Samples to Keep File:   (optional)");
      lblSamplesToKeep.setFont(COMPONENT_FONT);
      keepDropPanel.contentPanel.add(lblSamplesToKeep, "cell 0 4 3 1");
      txtFldSampKeep = new JTextField();
      keepDropPanel.contentPanel.add(txtFldSampKeep, "flowx,cell 0 5 3 1,growx");
      txtFldSampKeep.setColumns(10);
      JButton btnSampKeep = new JButton(">");
      btnSampKeep.addActionListener(e -> {
        String[] files = new FileChooser(ImputationGUI.this, proj.PROJECT_DIRECTORY.getValue(),
                                         false, false, "Select Sample Keeps File", log).getFiles();
        if (files == null) {
          return;
        } else {
          txtFldSampKeep.setText(files[0]);
        }
      });
      keepDropPanel.contentPanel.add(btnSampKeep, "pad 2 0 -2 0,cell 0 5 3 1,alignx right");
      JLabel lblSamplesToDropfile = new JLabel("Samples to Drop File:   (optional)");
      lblSamplesToDropfile.setFont(COMPONENT_FONT);
      keepDropPanel.contentPanel.add(lblSamplesToDropfile, "cell 0 6 3 1");
      txtFldMkrKeep = new JTextField();
      keepDropPanel.contentPanel.add(txtFldMkrKeep, "flowx,cell 0 1 3 1,growx");
      txtFldMkrKeep.setColumns(10);
      JButton btnMkrKeep = new JButton(">");
      btnMkrKeep.addActionListener(e -> {
        String[] files = new FileChooser(ImputationGUI.this, proj.PROJECT_DIRECTORY.getValue(),
                                         false, false, "Select Marker Keeps File", log).getFiles();
        if (files == null) {
          return;
        } else {
          txtFldMkrKeep.setText(files[0]);
        }
      });
      keepDropPanel.contentPanel.add(btnMkrKeep, "pad 2 0 -3 0,cell 0 1 3 1,alignx right");
      txtFldSampDrop = new JTextField();
      keepDropPanel.contentPanel.add(txtFldSampDrop, "flowx,cell 0 7 3 1,growx");
      txtFldSampDrop.setColumns(10);
      JButton btnSampDrop = new JButton(">");
      btnSampDrop.addActionListener(e -> {
        String[] files = new FileChooser(ImputationGUI.this, proj.PROJECT_DIRECTORY.getValue(),
                                         false, false, "Select Sample Drops File", log).getFiles();
        if (files == null) {
          return;
        } else {
          txtFldSampDrop.setText(files[0]);
        }
      });
      keepDropPanel.contentPanel.add(btnSampDrop, "pad 2 0 -2 0,cell 0 7 3 1,alignx right");
      JLabel lblKeepDrop = new JLabel("<html>&nbsp;-&nbsp;<u>Keep/Drop Files</u></html>");
      lblKeepDrop.setFont(HEADER_FONT);
      keepDropPanel.topPanel.add(lblKeepDrop, "cell 0 0, pad 0 0 0 14");

      JAccordionPanel chrsPanel = new JAccordionPanel();
      panel.add(chrsPanel, "cell 0 3 2 1,growx");
      chrsPanel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
      chrsPanel.contentPanel.setLayout(new MigLayout("", "[grow]", "[][grow]"));

      JLabel lblChrsOptionalNote = new JLabel("(optional - leave blank for all)");
      lblChrsOptionalNote.setFont(COMPONENT_FONT.deriveFont(10f));
      chrsPanel.contentPanel.add(lblChrsOptionalNote, "cell 0 0,alignx left,aligny top");
      JScrollPane chrScroll = new JScrollPane();
      chrsPanel.contentPanel.add(chrScroll, "cell 0 1,grow");
      Integer[] chrs = new Integer[27];
      for (int i = 0; i < chrs.length; i++) {
        chrs[i] = i;
      }
      chrList = new JList<>(chrs);
      chrList.setCellRenderer(new DefaultListCellRenderer() {

        @Override
        public Component getListCellRendererComponent(JList<?> list, Object value, int index,
                                                      boolean isSelected, boolean cellHasFocus) {
          // Pad display with spaces:
          Component c = super.getListCellRendererComponent(list, String.format("%1$6d", value),
                                                           index, isSelected, cellHasFocus);
          return c;
        }
      });
      chrList.setFont(COMPONENT_FONT);
      chrList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
      chrList.setVisibleRowCount(3);
      chrScroll.setViewportView(chrList);

      String selChrLbl = "<html>&nbsp;-&nbsp;<u>Select Chromosomes</u></html>";
      JLabel lblSelectChromosomes = new JLabel(selChrLbl);
      lblSelectChromosomes.setFont(HEADER_FONT);
      chrsPanel.topPanel.add(lblSelectChromosomes, "cell 0 0");

      JAccordionPanel otherReqsPanel = new JAccordionPanel();
      panel.add(otherReqsPanel, "cell 0 4 2 1,growx");
      otherReqsPanel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));

      String othReqLbl = "<html>&nbsp;-&nbsp;<u>Other Requirements</u></html>";
      JLabel lblOtherRequirements = new JLabel(othReqLbl);
      lblOtherRequirements.setFont(HEADER_FONT);
      otherReqsPanel.topPanel.add(lblOtherRequirements, "cell 0 0");
      otherReqsPanel.contentPanel.setLayout(new MigLayout("", "[grow]", "[][][][][]"));

      JLabel lblReferenceFile = new JLabel("Reference File:   (required)");
      lblReferenceFile.setFont(COMPONENT_FONT);
      otherReqsPanel.contentPanel.add(lblReferenceFile, "flowx,cell 0 0");

      txtFldRefFile = new JTextField();
      otherReqsPanel.contentPanel.add(txtFldRefFile, "flowx,cell 0 1,growx");
      txtFldRefFile.setColumns(10);

      JButton btnRefFile = new JButton(">");
      btnRefFile.addActionListener(e -> {
        String[] files = new FileChooser(ImputationGUI.this, proj.PROJECT_DIRECTORY.getValue(),
                                         false, false, "Select Reference / Site Map File", log)
                                                                                               .getFiles();
        if (files == null) {
          return;
        } else {
          try {
            ImputationPrep.validateRefFile(files[0], log);
            txtFldRefFile.setText(files[0]);
          } catch (IOException ex) {
            // error reading ref file
            proj.message("Error reading reference file: " + ex.getMessage());
          } catch (IllegalArgumentException ex) {
            // bad ref file
            proj.message("Malformed reference file! - " + ex.getMessage());
          }
        }
      });
      otherReqsPanel.contentPanel.add(btnRefFile, "pad 2 0 -3 0,cell 0 1");

      JLabel lblOutputDirectoryfiles = new JLabel("Output Directory:   (files may be large)");
      lblOutputDirectoryfiles.setFont(COMPONENT_FONT);
      otherReqsPanel.contentPanel.add(lblOutputDirectoryfiles, "cell 0 2");

      txtFldOutDir = new JTextField();
      txtFldOutDir.setColumns(10);
      otherReqsPanel.contentPanel.add(txtFldOutDir, "flowx,cell 0 3,growx");

      JButton btnOutDir = new JButton(">");
      btnOutDir.addActionListener(e -> {
        String dir = new FileChooser(ImputationGUI.this, proj.PROJECT_DIRECTORY.getValue(), false,
                                     true, "Select Output Directory", log).getNavDir();
        if (dir == null || "".equals(dir)) {
          return;
        } else {
          txtFldOutDir.setText(dir);
        }
      });
      otherReqsPanel.contentPanel.add(btnOutDir, "pad 2 0 -3 0,cell 0 3");

      JLabel lblRefHelp = Grafik.getToolTipIconLabel(refFileDesc);
      otherReqsPanel.contentPanel.add(lblRefHelp, "cell 0 0");

      chrsPanel.shrink();
    }
    {
      JPanel buttonPane = new JPanel();
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      buttonPane.setLayout(new MigLayout("ins 5", "[grow][51px][65px]", "[][][23px]"));
      {
        JSeparator separator = new JSeparator();
        buttonPane.add(separator, "cell 0 0 3 1,growx");
      }

      JLabel lblCreate = new JLabel("Create:");
      buttonPane.add(lblCreate, "cell 0 1,alignx right");

      JLabel lblOr = new JLabel("or:");
      buttonPane.add(lblOr, "cell 0 2,alignx right");
      {
        JButton okButton = new JButton("Run");
        okButton.addActionListener(e -> {
          run();
        });
        buttonPane.add(okButton, "cell 1 2,growx,aligny top");
      }
      {
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> {
          close();
        });
        buttonPane.add(cancelButton, "cell 2 2,growx,aligny top");
        getRootPane().setDefaultButton(cancelButton);
      }
      {
        JButton qsubButton = new JButton("QSUB");
        qsubButton.addActionListener(e -> {
          qsub();
        });
        buttonPane.add(qsubButton, "cell 1 1,growx");
      }
      {
        JButton btnScript = new JButton("SCRIPT");
        btnScript.addActionListener(e -> {
          scriptify();
        });
        buttonPane.add(btnScript, "cell 2 1,growx");
      }
    }
    init = false;
  }

  /**
   * Return an enum based on selections
   * 
   * @return One of VCF_ONLY, PLINK_ONLY, PLINK_SHAPEIT, or PLINK_SHAPEIT_MINIMAC. SHAPEIT and
   *         MINIMAC options are not exposed through this gui.
   */
  public IMPUTATION_PIPELINE_PATH getImputationPipeline() {
    boolean selV = rdbtnVcf.isSelected();
    if (selV) {
      return IMPUTATION_PIPELINE_PATH.VCF_ONLY;
    } else {
      boolean sh, mn;
      sh = chckbxShape.isSelected();
      mn = chckbxMini.isSelected();
      if (sh && mn) {
        return IMPUTATION_PIPELINE_PATH.PLINK_SHAPEIT_MINIMAC;
      } else if (sh) {
        return IMPUTATION_PIPELINE_PATH.PLINK_SHAPEIT;
      } else {
        return IMPUTATION_PIPELINE_PATH.PLINK_ONLY;
      }
    }
  }

  public int[] getChromosomes() {
    // indices match chrs, as both start at 0 and incr by 1
    return chrList.getSelectedIndices();
  }

  /**
   * @return A string array of file paths, or blank values, in the following order: MarkerKeeps,
   *         MarkerDrops, SampleKeeps, SampleDrops. Files are not guaranteed to exist or be valid.
   */
  public KeepDrops getKeepDropFiles() {
    String dropSamples = txtFldSampDrop.getText();
    String keepSamples = txtFldSampKeep.getText();
    String dropMarkers = txtFldMkrDrop.getText();
    String keepMarkers = txtFldMkrKeep.getText();

    if (dropSamples.isEmpty()) dropSamples = null;
    if (keepSamples.isEmpty()) keepSamples = null;
    if (dropMarkers.isEmpty()) dropMarkers = null;
    if (keepMarkers.isEmpty()) keepMarkers = null;

    return new KeepDrops(dropSamples, keepSamples, dropMarkers, keepMarkers);
  }

  public String getReferenceFile() {
    return txtFldRefFile.getText();
  }

  public boolean getUseGRC() {
    return chckbxUseGRC.isSelected();
  }

  public String getOutputDirectory() {
    return txtFldOutDir.getText();
  }

  private void run() {
    if (!checkRequirementsOrMessage()) {
      return;
    }
    String propFile = proj.getPropertyFilename();
    int[] chrs = getChromosomes();
    String ref = getReferenceFile();
    String outputDir = getOutputDirectory();

    switch (getImputationPipeline()) {
      case VCF_ONLY:
        ImputationPipeRunner.runVCF(propFile, chrs, ref, getKeepDropFiles(), false,
                                    outputDir + "vcf/" + ext.replaceWithLinuxSafeCharacters(proj.PROJECT_NAME.getValue()),
                                    getUseGRC());
        break;
      case PLINK_ONLY:
        ImputationPipeRunner.runPlink(proj.getPropertyFilename(), getChromosomes(),
                                      getReferenceFile(), getKeepDropFiles(),
                                      outputDir + "/plink/plink");
        break;
      case PLINK_SHAPEIT:
        ImputationPipeRunner.runPlinkAndShapeIt(propFile, chrs, ref, getKeepDropFiles(), outputDir);
        break;
      case PLINK_SHAPEIT_MINIMAC:
        ImputationPipeRunner.runPlinkShapeItAndMinimac(propFile, chrs, ref, getKeepDropFiles(),
                                                       outputDir);
        break;
      default:
        break;
    }
    close();
  }

  private void close() {
    this.setVisible(false);
    this.dispose();
  }

  private String getRunString() {
    String propFile = proj.getPropertyFilename();
    int[] chrs = getChromosomes();
    String ref = getReferenceFile();
    // String plinkSubdir = null; // used for default keep/drop files; ignore here
    String outputDir = getOutputDirectory();
    IMPUTATION_PIPELINE_PATH path = getImputationPipeline();

    StringJoiner imputeStr = new StringJoiner(" ");
    imputeStr.add(Files.getRunString()).add(ImputationPipeline.class.getName());
    imputeStr.add(ImputationPipeline.PROJ_ARG + propFile);
    imputeStr.add(ImputationPipeline.RUN_TYPE_ARG + path.name());
    if (chrs.length > 0) {
      imputeStr.add(ImputationPipeline.CHRS_ARG + ArrayUtils.toStr(chrs, ","));
    }
    imputeStr.add(ImputationPipeline.REF_ARG + ref);
    for (String keepDropArg : generateKeepDropsArgs()) {
      imputeStr.add(keepDropArg);
    }

    switch (path) {
      case VCF_ONLY:
        boolean useGRC = getUseGRC();
        String vcfOutRoot = outputDir + "vcf/"
                            + ext.replaceWithLinuxSafeCharacters(proj.PROJECT_NAME.getValue());
        imputeStr.add(ImputationPipeline.OUT_DIR_AND_ROOT_ARG + vcfOutRoot);
        imputeStr.add(ImputationPipeline.USE_GRC_ARG + useGRC);
        break;
      case PLINK_ONLY:
        imputeStr.add(ImputationPipeline.OUT_DIR_AND_ROOT_ARG + outputDir + "/plink/plink");
        break;
      case PLINK_SHAPEIT:
      case PLINK_SHAPEIT_MINIMAC:
        imputeStr.add(ImputationPipeline.OUT_DIR_ARG + outputDir);
        break;
      default:
        System.err.println("Error - unrecognized imputation path: " + path);
        break;
    }
    return imputeStr.toString();
  }

  private Collection<String> generateKeepDropsArgs() {
    ImmutableSet.Builder<String> keepDropArgsBuilder = ImmutableSet.builder();
    KeepDrops keepDrops = getKeepDropFiles();
    if (keepDrops.getDropSamplesFile() != null) keepDropArgsBuilder.add(ImputationPipeline.DROP_SAMPLES_ARG
                                                                        + keepDrops.getDropSamplesFile());
    if (keepDrops.getKeepSamplesFile() != null) keepDropArgsBuilder.add(ImputationPipeline.KEEP_SAMPLES_ARG
                                                                        + keepDrops.getKeepSamplesFile());
    if (keepDrops.getDropMarkersFile() != null) keepDropArgsBuilder.add(ImputationPipeline.DROP_MARKERS_ARG
                                                                        + keepDrops.getDropMarkersFile());
    if (keepDrops.getKeepMarkersFile() != null) keepDropArgsBuilder.add(ImputationPipeline.KEEP_MARKERS_ARG
                                                                        + keepDrops.getKeepMarkersFile());
    return keepDropArgsBuilder.build();
  }

  private boolean checkRequirementsOrMessage() {
    String refFile = getReferenceFile();
    if (!Files.exists(refFile)) {
      proj.message("Error - a Reference / Site List file is required!");
      return false;
    }
    try {
      ImputationPrep.validateRefFile(refFile, log);
    } catch (IllegalArgumentException e) {
      proj.message("Error - malformed Reference / Site List file: " + e.getMessage());
      return false;
    } catch (IOException e) {
      proj.message("Error - problem reading Reference / Site List file: " + e.getMessage());
      return false;
    }
    String outDir = getOutputDirectory();
    if ("".equals(outDir)) {
      proj.message("Error - output directory not specified!");
      return false;
    }
    if (!Files.exists(outDir)) {
      int opt = JOptionPane.showConfirmDialog(ImputationGUI.this,
                                              "Output directory not found - would you like to create it?",
                                              "Create Output Directory?",
                                              JOptionPane.YES_NO_OPTION);
      if (opt == JOptionPane.YES_OPTION) {
        if (!new File(outDir).mkdirs()) {
          proj.message("Error - couldn't create the desired output directory!");
          return false;
        }
      }
    }
    return true;
  }

  private void qsub() {
    if (checkRequirementsOrMessage()) {
      Qsub.qsubDefaults(proj.PROJECT_DIRECTORY.getValue() + "ImputationPipeline."
                        + ext.replaceWithLinuxSafeCharacters(ext.getDate()) + ".pbs",
                        getRunString());
      close();
    }
  }

  private void scriptify() {
    if (checkRequirementsOrMessage()) {
      Files.write(getRunString(), proj.PROJECT_DIRECTORY.getValue() + "ImputationPipeline."
                                  + ext.replaceWithLinuxSafeCharacters(ext.getDate()) + ".run");
      close();
    }
  }

}
