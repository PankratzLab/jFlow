package org.genvisis.one.ben.imagetag;

import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.border.BevelBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.plaf.basic.BasicTreeUI;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import org.genvisis.one.ben.imagetag.AnnotatedImage.Annotation;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import net.miginfocom.swing.MigLayout;

public class ImageAnnotator {

  private JFrame frmAnnotator;
  private JComboBox<String> sampleCombo;
  private JTree tree;

  IAnnotator annotator = new Annotator();

  BufferedImage selectedImage = AnnotatedImage.createReadyImage();

  private JMenu mnLoadRecent;
  private JTextField fldNewAnnot;
  private JPanel annotPanel;
  private JPanel imagePanel;
  private JPanel controlPanel;
  private JSplitPane splitPane;
  private JMenu mnTravAnn;
  private JMenu mnDelAnn;
  private JRadioButtonMenuItem jrbTravAll;
  private JRadioButtonMenuItem jrbTravAnn;
  private JRadioButtonMenuItem jrbTravNon;
  private ButtonGroup travGrp;
  private HashMap<Annotation, JRadioButtonMenuItem> annTravMap = new HashMap<>();
  private HashMap<Annotation, JMenuItem> annDelMap = new HashMap<>();
  private HashMap<Annotation, JMenuItem> annSaveMap = new HashMap<>();

  private static final String PROP_FILE = ".imageannotator.properties";
  private static final String KEY_LAST_DIR_IMG = "LAST_DIR_IMG";
  private static final String KEY_LAST_DIR_ANN = "LAST_DIR_ANN";
  private static final String KEY_LAST_FILE_ANN = "LAST_FILE_ANN";
  private static final String KEY_LAST_GATE = "LAST_SEL_NODE";
  private static final String KEY_RECENT = "RECENT";

  private String lastOpenedImageDir = null;
  private String lastOpenedAnnFileDir = null;
  private String lastSavedAnnFile = null;
  private String lastSelectedGate = null;
  private String lastSelectedFile = null;
  private HashSet<String> recentAnnotFiles = new HashSet<>();

  private HashMap<Character, Action> mnemonicActions = new HashMap<>();

  private static final int ALL = 0;
  private static final int ANN = 1;
  private static final int NON = 2;
  private static final String BACKUP_DIR = ".backup./";
  private static final String BACKUP_FILE = BACKUP_DIR + "backup.annotations";
  private volatile boolean annotationsChanged = false;
  private volatile boolean checkBxAdvanceAfterAnnotationKeyed = true;

  private void setAnnotationsChanged() {
    annotationsChanged = true;
  }

  private boolean checkAnnotationsChanged() {
    boolean value = annotationsChanged;
    annotationsChanged = false;
    return value;
  }

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    EventQueue.invokeLater(new Runnable() {

      public void run() {
        try {
          ImageAnnotator window = new ImageAnnotator();
          window.frmAnnotator.setVisible(true);
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    });
  }

  /**
   * Create the application.
   */
  public ImageAnnotator() {
    loadProperties();
    initialize();
  }

  private void loadProperties() {
    Properties p = new Properties();
    if (!Files.exists(PROP_FILE)) return;
    try {
      p.load(Files.getAppropriateInputStreamReader(PROP_FILE));
      String lstAnn = p.getProperty(KEY_LAST_DIR_ANN, "");
      String lstImg = p.getProperty(KEY_LAST_DIR_IMG, "");
      String lstFil = p.getProperty(KEY_LAST_FILE_ANN, "");
      String lstGat = p.getProperty(KEY_LAST_GATE, "");
      String rec = p.getProperty(KEY_RECENT, "");
      lastOpenedAnnFileDir = lstAnn;
      lastOpenedImageDir = lstImg;
      lastSavedAnnFile = lstFil;
      lastSelectedGate = lstGat;
      if (!"".equals(rec)) {
        String[] ps = rec.split(";");
        for (String s : ps) {
          if (Files.exists(s)) {
            recentAnnotFiles.add(s);
          }
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void saveProperties() {
    Properties p = new Properties();
    if (lastOpenedImageDir != null && !"".equals(lastOpenedImageDir)
        && Files.exists(lastOpenedImageDir)) {
      p.setProperty(KEY_LAST_DIR_IMG, lastOpenedImageDir);
    }
    if (lastOpenedAnnFileDir != null && !"".equals(lastOpenedAnnFileDir)
        && Files.exists(lastOpenedAnnFileDir)) {
      p.setProperty(KEY_LAST_DIR_ANN, lastOpenedAnnFileDir);
    }
    if (lastSavedAnnFile != null && !"".equals(lastSavedAnnFile)
        && Files.exists(lastSavedAnnFile)) {
      p.setProperty(KEY_LAST_FILE_ANN, lastSavedAnnFile);
    }
    if (lastSelectedGate != null && !"".equals(lastSelectedGate)) {
      p.setProperty(KEY_LAST_GATE, lastSelectedGate);
    }
    if (!recentAnnotFiles.isEmpty()) {
      p.setProperty(KEY_RECENT, ArrayUtils.toStr(recentAnnotFiles, ";"));
    }
    try {
      p.store(Files.getAppropriateWriter(PROP_FILE), "");
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void setLastUsedAnnotationDir(String annDir) {
    this.lastOpenedAnnFileDir = annDir;
  }

  private void setLastSavedAnnotationFile(String annFile) {
    this.lastSavedAnnFile = annFile;
  }

  private String getLastUsedAnnotationDir() {
    if (this.lastOpenedAnnFileDir == null || "".equals(this.lastOpenedAnnFileDir)
        || !Files.exists(this.lastOpenedAnnFileDir)) {
      try {
        this.lastOpenedAnnFileDir = ext.parseDirectoryOfFile(new File(ImageAnnotator.class.getProtectionDomain()
                                                                                          .getCodeSource()
                                                                                          .getLocation()
                                                                                          .toURI()
                                                                                          .getPath()).getAbsolutePath());
      } catch (URISyntaxException e) {
        this.lastOpenedAnnFileDir = "";
      }
    }
    return this.lastOpenedAnnFileDir;
  }

  private void setLastUsedImageDir(String imgDir) {
    this.lastOpenedImageDir = imgDir;
  }

  private String getLastUsedImageDir() {
    if (this.lastOpenedImageDir == null || "".equals(this.lastOpenedImageDir)
        || !Files.exists(this.lastOpenedImageDir)) {
      try {
        this.lastOpenedImageDir = ext.parseDirectoryOfFile(new File(ImageAnnotator.class.getProtectionDomain()
                                                                                        .getCodeSource()
                                                                                        .getLocation()
                                                                                        .toURI()
                                                                                        .getPath()).getAbsolutePath());
      } catch (URISyntaxException e) {
        this.lastOpenedImageDir = "";
      }
    }
    return this.lastOpenedImageDir;
  }

  /**
   * Initialize the contents of the frame.
   */
  private void initialize() {
    frmAnnotator = new JFrame();
    frmAnnotator.setTitle("ImageAnnotator");
    frmAnnotator.setBounds(100, 100, 1200, 800);
    frmAnnotator.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    frmAnnotator.getContentPane().setLayout(new MigLayout("", "[grow]", "[grow]"));
    frmAnnotator.addWindowListener(new WindowAdapter() {

      @Override
      public void windowClosing(WindowEvent e) {
        close();
      }
    });

    splitPane = new JSplitPane();
    splitPane.setResizeWeight(0.94);
    frmAnnotator.getContentPane().add(splitPane, "cell 0 0,grow");

    imagePanel = new JPanel() {

      @Override
      protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Image scaledImage = selectedImage.getScaledInstance(getWidth(), getHeight(),
                                                            Image.SCALE_SMOOTH);
        g.drawImage(scaledImage, 0, 0, getWidth(), getHeight(), null);
      }
    };
    splitPane.setLeftComponent(imagePanel);
    imagePanel.setBorder(new BevelBorder(BevelBorder.RAISED, null, null, null, null));

    JSplitPane splitPane_1 = new JSplitPane();
    splitPane_1.setResizeWeight(0.5);
    splitPane_1.setOrientation(JSplitPane.VERTICAL_SPLIT);
    splitPane.setRightComponent(splitPane_1);

    annotPanel = new JPanel();
    splitPane_1.setRightComponent(annotPanel);
    annotPanel.setLayout(new MigLayout("", "[grow]", "[][][][grow]"));

    JLabel lblAddAnnotation = new JLabel("Add Annotation:");
    annotPanel.add(lblAddAnnotation, "cell 0 0");

    fldNewAnnot = new JTextField();
    ActionMap am = fldNewAnnot.getActionMap();
    InputMap im = fldNewAnnot.getInputMap(JComponent.WHEN_FOCUSED);
    annotPanel.add(fldNewAnnot, "cell 0 1,growx");

    JLabel lblExistingAnnotations = new JLabel("Existing Annotations:");
    annotPanel.add(lblExistingAnnotations, "cell 0 2");

    JScrollPane scrollPane = new JScrollPane();
    annotPanel.add(scrollPane, "cell 0 3,grow");

    annotPanel = new JPanel();
    scrollPane.setViewportView(annotPanel);
    annotPanel.setLayout(new MigLayout("ins 0", "[]", "[]"));

    controlPanel = new JPanel();
    splitPane_1.setLeftComponent(controlPanel);
    controlPanel.setLayout(new MigLayout("", "[110.00,grow]", "[][][grow][][]"));

    sampleCombo = new JComboBox<>();
    sampleCombo.addActionListener(comboListener);

    sampleCombo.setFocusable(false);
    controlPanel.add(sampleCombo, "flowx,cell 0 1,growx");

    JScrollPane scrollPane_1 = new JScrollPane();
    controlPanel.add(scrollPane_1, "cell 0 2,growy, growx, pad 0 0 14 0");

    tree = new JTree();
    scrollPane_1.setViewportView(tree);
    constructTree();
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0), "enter");
    am.put("enter", new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent arg0) {
        createNewAnnotation();
        tree.requestFocusInWindow();
      }
    });

    frmAnnotator.setJMenuBar(createMenuBar());
    checkForBackupFileOrLoadLast();

    startAutoSaveThread();
  }

  public static final String CHR_POS_REGEX = "(.*)_?chr([12]?[0-9[XYM]]+)[:-_]([\\d,]+)[-_]([\\d,]+).*?";

  private int[] parseSampleChrPosIfExists() {
    Matcher m = Pattern.compile(CHR_POS_REGEX).matcher(lastSelectedFile);
    if (m.matches()) {
      byte chr = Positions.chromosomeNumber(m.group(2), false, new Logger());
      int stt = Integer.parseInt(m.group(3));
      int stp = Integer.parseInt(m.group(4));
      return new int[] {chr, stt, stp};
    }
    return null;
  }

  private void constructTree() {
    DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
    DefaultTreeModel dtm = new DefaultTreeModel(rootNode);
    tree.setModel(dtm);
    tree.setShowsRootHandles(false);
    tree.setRootVisible(false);
    tree.addTreeSelectionListener(treeListener);
    tree.setEnabled(false);
    updateTreeKeys(tree);
    expandAllNodes(tree);
  }

  long lastTimeSaved = -1;
  Thread autoSaveThread = null;

  private void checkForBackupFileOrLoadLast() {
    if (Files.exists(BACKUP_FILE)) {
      int opt = JOptionPane.showConfirmDialog(this.frmAnnotator, "Load Auto-saved backup file?",
                                              "Auto-Save File Detected", JOptionPane.YES_NO_OPTION,
                                              JOptionPane.WARNING_MESSAGE);
      if (opt == JOptionPane.YES_OPTION) {
        loadAnnotationFile(BACKUP_FILE);
        new File(BACKUP_FILE).delete();
      } else {
        if (!"".equals(lastSavedAnnFile) && Files.exists(lastSavedAnnFile)) {
          loadAnnotationFile(lastSavedAnnFile);
          selectLastSelectedGate();
        }
      }
    } else {
      if (!"".equals(lastSavedAnnFile) && Files.exists(lastSavedAnnFile)) {
        loadAnnotationFile(lastSavedAnnFile);
        selectLastSelectedGate();
      }
    }
  }

  private void selectLastSelectedGate() {
    if (lastSelectedGate != null && !"".equals(lastSelectedGate)) {
      for (int i = 0; i < tree.getRowCount(); i++) {
        TreePath tp = tree.getPathForRow(i);
        AnnotatedImage ai = ((AnnotatedImage) ((DefaultMutableTreeNode) tp.getLastPathComponent()).getUserObject());
        String gate = ai.getName();
        if (gate.equals(lastSelectedGate)) {
          lastSelectedFile = ai.getImageFile();
          tree.setSelectionPath(tp);
          break;
        }
      }
    }
  }

  private void startAutoSaveThread() {
    autoSaveThread = new Thread(new Runnable() {

      @Override
      public void run() {
        while (true) {
          if (lastTimeSaved == -1 || (System.currentTimeMillis() - lastTimeSaved > 30000)) {
            if (checkAnnotationsChanged()) {
              new File(BACKUP_DIR).mkdir();
              try {
                System.out.println("Auto-saving to backup: "
                                   + new File(BACKUP_FILE).getCanonicalPath());
              } catch (IOException e) {}
              annotator.saveAnnotations(BACKUP_FILE);
            }
            lastTimeSaved = System.currentTimeMillis();
          }
          try {
            Thread.sleep(15000);
          } catch (InterruptedException e) {
            e.printStackTrace();
          }
        }
      }
    }, "AutoSaveThread");
    autoSaveThread.setDaemon(true);
    autoSaveThread.start();
  }

  private JMenuBar createMenuBar() {
    JMenuBar menuBar = new JMenuBar();

    JMenu mnFile = new JMenu("File");
    mnFile.setMnemonic('F');
    menuBar.add(mnFile);

    JMenuItem mntmLoadImgDir = new JMenuItem();
    mntmLoadImgDir.setAction(loadAction);
    mntmLoadImgDir.setText("Load Image Directory");
    mntmLoadImgDir.setMnemonic('L');
    mnFile.add(mntmLoadImgDir);

    JMenuItem mntmLoadAnnotFile = new JMenuItem();
    mntmLoadAnnotFile.setAction(existAction);
    mntmLoadAnnotFile.setText("Load Annotation File");
    mntmLoadAnnotFile.setMnemonic('A');
    mnFile.add(mntmLoadAnnotFile);

    mnLoadRecent = new JMenu("Load Recent Annotation File");
    mnLoadRecent.setMnemonic('R');
    mnFile.add(mnLoadRecent);

    updateRecentFiles();

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmSaveAnnot = new JMenuItem();
    mntmSaveAnnot.setAction(saveAction);
    mntmSaveAnnot.setText("Save Annotations");
    mntmSaveAnnot.setMnemonic('S');
    mnFile.add(mntmSaveAnnot);

    JMenuItem mntmSaveAsAnnot = new JMenuItem();
    mntmSaveAsAnnot.setAction(saveAsAction);
    mntmSaveAsAnnot.setText("Save Annotations to File");
    mntmSaveAsAnnot.setMnemonic('A');
    mnFile.add(mntmSaveAsAnnot);

    mnSaveAnn = new JMenu("Save Annotation Set To File");
    mnFile.add(mnSaveAnn);

    JMenuItem mntmExport = new JMenuItem();
    mntmExport.setAction(exportAction);
    mntmExport.setText("Export Annotated Samples List");
    mntmExport.setMnemonic('E');
    mnFile.add(mntmExport);

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmCopyUCSC = new JMenuItem();
    mntmCopyUCSC.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(Positions.getUCSCformat(parseSampleChrPosIfExists()));
      }
    });
    mntmCopyUCSC.setText("Copy UCSC Location to Clipboard");
    mntmCopyUCSC.setMnemonic('U');
    mnFile.add(mntmCopyUCSC);

    JMenuItem mntmCopyUCSCLink = new JMenuItem();
    mntmCopyUCSCLink.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(Positions.getUCSClink(parseSampleChrPosIfExists()));
      }
    });
    mntmCopyUCSCLink.setText("Copy UCSC Link to Clipboard");
    mntmCopyUCSCLink.setMnemonic('L');
    mnFile.add(mntmCopyUCSCLink);

    JMenuItem mntmCopySample = new JMenuItem();
    mntmCopySample.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(sampleCombo.getSelectedItem().toString());
      }
    });
    mntmCopySample.setText("Copy Sample Name to Clipboard");
    mntmCopySample.setMnemonic('N');
    mnFile.add(mntmCopySample);

    JMenuItem mntmCopySampleFile = new JMenuItem();
    mntmCopySampleFile.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(lastSelectedFile);
      }
    });
    mntmCopySampleFile.setText("Copy File Path to Clipboard");
    mntmCopySampleFile.setMnemonic('N');
    mnFile.add(mntmCopySampleFile);

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    mnDelAnn = new JMenu("Delete Annotation:");
    mnFile.add(mnDelAnn);

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmExit = new JMenuItem();
    mntmExit.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        close();
      }
    });
    mntmExit.setText("Exit");
    mnFile.add(mntmExit);

    JMenu mnNav = new JMenu("Navigation");
    mnNav.setMnemonic('N');
    menuBar.add(mnNav);

    JMenuItem mntmNavTravLbl = new JMenuItem();
    mntmNavTravLbl.setEnabled(false);
    mntmNavTravLbl.setText("Traversal:");
    mnNav.add(mntmNavTravLbl);

    travGrp = new ButtonGroup();

    jrbTravAll = new JRadioButtonMenuItem();
    jrbTravAll.setText("All");
    jrbTravAll.setSelected(true);
    travGrp.add(jrbTravAll);
    mnNav.add(jrbTravAll);
    jrbTravAnn = new JRadioButtonMenuItem();
    jrbTravAnn.setText("Annotated");
    travGrp.add(jrbTravAnn);
    mnNav.add(jrbTravAnn);
    jrbTravNon = new JRadioButtonMenuItem();
    jrbTravNon.setText("Non-Annotated");
    travGrp.add(jrbTravNon);
    mnNav.add(jrbTravNon);
    mnTravAnn = new JMenu("Annotation:");
    mnNav.add(mnTravAnn);

    mnNav.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmNavOptLbl = new JMenuItem();
    mntmNavOptLbl.setEnabled(false);
    mntmNavOptLbl.setText("Other Options:");
    mnNav.add(mntmNavOptLbl);

    JMenuItem mntmNavFind = new JMenuItem();
    mntmNavFind.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        List<String> files = annotator.getRoots();
        if (files.isEmpty()) {
          JOptionPane.showMessageDialog(ImageAnnotator.this.frmAnnotator, "No files to search!");
          return;
        }
        String[] v = files.toArray(new String[files.size()]);
        List<String> sel = FileFinder.showFileFinder(v, false);
        if (!sel.isEmpty()) {
          sampleCombo.setSelectedItem(sel.get(0));
        }
      }
    });
    mntmNavFind.setMnemonic('F');
    mntmNavFind.setText("Find Sample");
    mnNav.add(mntmNavFind);

    JCheckBoxMenuItem mntmNavAdvAfterAnnot = new JCheckBoxMenuItem();
    mntmNavAdvAfterAnnot.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        checkBxAdvanceAfterAnnotationKeyed = mntmNavAdvAfterAnnot.isSelected();
      }
    });
    mntmNavAdvAfterAnnot.setText("Advance after annotating (or SHIFT)");
    mntmNavAdvAfterAnnot.setSelected(true);
    mnNav.add(mntmNavAdvAfterAnnot);

    return menuBar;
  }

  private void addAnnotationToTraversalMenu(Annotation annot) {
    JRadioButtonMenuItem jrb = new JRadioButtonMenuItem();
    jrb.setText(annot.annotation);
    annTravMap.put(annot, jrb);
    travGrp.add(jrb);
    JMenuItem jmn = new JMenuItem();
    jmn.setAction(saveAnnotationAction);
    jmn.setText(annot.annotation);
    annSaveMap.put(annot, jmn);
    JMenuItem jmd = new JMenuItem();
    jmd.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        deleteAnnotation(annot);
      }
    });
    jmd.setText(annot.annotation);
    mnDelAnn.add(jmd);
    annDelMap.put(annot, jmd);
    mnSaveAnn.add(jmn);
  }

  private void deleteAnnotation(Annotation annot) {
    JRadioButtonMenuItem jrb = annTravMap.get(annot);
    travGrp.remove(jrb);
    mnTravAnn.remove(jrb);
    mnSaveAnn.remove(annSaveMap.get(annot));
    mnDelAnn.remove(annDelMap.get(annot));
    annotator.deleteAnnotation(annot);

    setAnnotationsChanged();
    refreshAnnotations();
  }

  private AbstractAction saveAnnotationAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      String annot = e.getActionCommand();
      Annotation ann = null;
      for (Annotation a : annotator.getAnnotations()) {
        if (annot.equals(a.annotation)) {
          ann = a;
          break;
        }
      }
      if (ann == null) {
        return;
      }

      JFileChooser jfc = new JFileChooser(getLastUsedAnnotationDir());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save \"" + ann.annotation + "\" Annotation to File");
      int opt = jfc.showSaveDialog(frmAnnotator);
      if (opt == JFileChooser.APPROVE_OPTION) {
        String annFile = jfc.getSelectedFile().getAbsolutePath();
        setLastUsedAnnotationDir(ext.verifyDirFormat(ext.parseDirectoryOfFile(annFile)));
        annotator.saveAnnotation(ann, annFile);
        saveProperties();
      }
    }
  };

  private void close() {
    if (checkAnnotationsChanged()
        && (lastTimeSaved == -1 || System.currentTimeMillis() - lastTimeSaved > 500)) {
      if (!promptToSaveAnnotations()) {
        setAnnotationsChanged();
        return;
      }
    }
    saveProperties();
    frmAnnotator.setVisible(false);
    frmAnnotator.dispose();
  }

  private void updateRecentFiles() {
    mnLoadRecent.removeAll();
    if (recentAnnotFiles.isEmpty()) {
      JMenuItem mnRecentPlacehldr = new JMenuItem("No Recent Files");
      mnRecentPlacehldr.setEnabled(false);
      mnLoadRecent.add(mnRecentPlacehldr);
    } else {
      for (String s : recentAnnotFiles) {
        JMenuItem mnRecent = new JMenuItem();
        mnRecent.setAction(new AbstractAction() {

          @Override
          public void actionPerformed(ActionEvent e) {
            loadAnnotationFile(s);
          }
        });
        mnRecent.setText(ext.removeDirectoryInfo(s));
        mnLoadRecent.add(mnRecent);
      }
    }
  }

  protected void createNewAnnotation() {
    String newAnnotation = fldNewAnnot.getText().trim();
    fldNewAnnot.setText("");
    if ("".equals(newAnnotation)) {
      return;
    }
    AnnotatedImage.Annotation newAnn = new AnnotatedImage.Annotation(newAnnotation);
    if (!annotator.getAnnotations().contains(newAnn)) {
      addAnnotationToTraversalMenu(newAnn);
      annotator.addNewAnnotation(newAnn);
      setAnnotationsChanged();
      refreshAnnotations();
    }
  }

  private void loadImageFiles() {
    JFileChooser jfc = new JFileChooser(getLastUsedImageDir());
    jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    jfc.setDialogTitle("Open Directory");
    int opt = jfc.showOpenDialog(frmAnnotator);
    if (opt == JFileChooser.APPROVE_OPTION) {
      File f = jfc.getSelectedFile();
      String fS = ext.verifyDirFormat(f.getAbsolutePath());
      setLastUsedImageDir(fS);
      annotator.loadImgDir(fS);
      setAnnotationsChanged();
      reloadControls();
      saveProperties();
    }
  }

  // org.apache.commons.lang3.StringUtils.isEmpty(cs)
  public static boolean isEmpty(final CharSequence cs) {
    return cs == null || cs.length() == 0;
  }

  // adapted from org.apache.commons.lang3.StringUtils.isNumeric(cs)
  public static boolean hasNumeric(final CharSequence cs) {
    if (isEmpty(cs)) {
      return false;
    }
    final int sz = cs.length();
    for (int i = 0; i < sz; i++) {
      if (Character.isDigit(cs.charAt(i))) {
        return true;
      }
    }
    return false;
  }

  private void reloadControls() {
    String[] fcsKeys = annotator.getRoots().toArray(new String[0]);
    // Comparator adapted from https://stackoverflow.com/questions/41085394/java-comparator-alphanumeric-strings
    Arrays.sort(fcsKeys, Comparator.comparingLong(s -> hasNumeric(
                                                                  // if there are any digits in the string, delete all non-digit characters and compare the resulting values
                                                                  s.toString()) ? Long.parseLong(s.toString().replaceAll("\\D", "")) : 0)
                                   // then delete all digit characters and compare the resulting strings
                                   .thenComparing(s -> s.toString().replaceAll("\\d", "")));
    DefaultComboBoxModel<String> dcbm = new DefaultComboBoxModel<>(fcsKeys);
    sampleCombo.setModel(dcbm);
    if (sampleCombo.getModel().getSize() > 0) {
      sampleCombo.setSelectedIndex(0);
    }
    tree.requestFocusInWindow();
  }

  private void refreshAnnotations() {
    mnemonicActions.clear();
    ArrayList<AnnotatedImage.Annotation> allAnnots = annotator.getAnnotations();
    annotPanel.removeAll();
    AnnotatedImage ai = null;
    if (tree.getSelectionModel().getSelectionPath() != null) {
      DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tree.getSelectionModel()
                                                                 .getSelectionPath()
                                                                 .getLastPathComponent();
      ai = (AnnotatedImage) dmtn.getUserObject();
    }
    for (int i = 0; i < allAnnots.size(); i++) {
      JCheckBox annBox = new JCheckBox();
      Annotation ann = allAnnots.get(i);
      List<Annotation> aiAnns = ai == null ? null : ai.getAnnotations();
      if (ai != null && aiAnns.contains(ann)) {
        annBox.setSelected(true);
      } else {
        annBox.setSelected(false);
      }
      final AnnotatedImage ai1 = ai;
      final int ind = i;
      AbstractAction mnemAct = new AbstractAction() {

        @Override
        public void actionPerformed(ActionEvent e) {
          if (ai1 != null) {
            AnnotatedImage.Annotation a = allAnnots.get(ind);
            ArrayList<AnnotatedImage.Annotation> myAnn = ai1.getAnnotations();
            if (myAnn.contains(a)) {
              myAnn.remove(a);
              annBox.setSelected(false);
            } else {
              myAnn.add(a);
              annBox.setSelected(true);
            }
            setAnnotationsChanged();
            annotPanel.repaint();
          }
        }
      };
      annBox.setAction(mnemAct);
      annBox.setText(allAnnots.get(i).annotation);
      annBox.setMnemonic(allAnnots.get(i).mnemonic);
      annBox.setFocusable(false);
      mnemonicActions.put(allAnnots.get(i).mnemonic, mnemAct);
      annotPanel.add(annBox, "cell 0 " + i);
    }
    annotPanel.revalidate();
    annotPanel.repaint();
  }

  private void fireMnem(KeyEvent e) {
    char code = (e.getKeyChar() + "").toUpperCase().charAt(0);
    if (mnemonicActions.containsKey(code)) {
      mnemonicActions.get(code).actionPerformed(null);
      setAnnotationsChanged();
      if (checkBxAdvanceAfterAnnotationKeyed && !e.isShiftDown()) {
        if (tree.getSelectionRows()[0] == tree.getRowCount() - 1) {
          keyRight();
        } else {
          keyDown();
        }
      }
    }
  }

  private int getTraversal() {
    if (jrbTravAll.isSelected()) {
      return ALL;
    } else if (jrbTravAnn.isSelected()) {
      return ANN;
    } else if (jrbTravNon.isSelected()) {
      return NON;
    }
    return -1;
  }

  private int getPrevNodeRow() {
    int trav = getTraversal();
    if (tree.isSelectionEmpty()) {
      return -1;
    }
    int[] rows = tree.getSelectionRows();
    if (rows == null || rows.length == 0) {
      return -1;
    }
    int row = rows[0];
    row--;
    while (row >= 0) {
      TreePath tp = tree.getPathForRow(row);
      DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tp.getLastPathComponent();
      AnnotatedImage ai = (AnnotatedImage) dmtn.getUserObject();
      if (trav == ALL) {
        break;
      } else if (trav == ANN) {
        if (!ai.getAnnotations().isEmpty()) {
          break;
        }
      } else if (trav == NON) {
        if (ai.getAnnotations().isEmpty()) {
          break;
        }
      } else {
        Annotation a = getTravAnnotation();
        if (a == null) return -1;
        if (ai.getAnnotations().contains(a)) {
          break;
        }
      }
      row--;
    }
    return row;
  }

  private Annotation getTravAnnotation() {
    Annotation sel = null;
    for (Entry<Annotation, JRadioButtonMenuItem> entry : annTravMap.entrySet()) {
      if (entry.getValue().isSelected()) {
        sel = entry.getKey();
        break;
      }
    }
    return sel;
  }

  private int getNextNodeRow() {
    int trav = getTraversal();
    if (tree.isSelectionEmpty()) {
      return -1;
    }
    int cnt = getTreeRowCount();
    int[] rows = tree.getSelectionRows();
    if (rows == null || rows.length == 0) {
      return -1;
    }
    int row = rows[0];
    row++;
    while (row < cnt) {
      TreePath tp = tree.getPathForRow(row);
      DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tp.getLastPathComponent();
      AnnotatedImage ai = (AnnotatedImage) dmtn.getUserObject();
      if (trav == ALL) {
        break;
      } else if (trav == ANN) {
        if (!ai.getAnnotations().isEmpty()) {
          break;
        }
      } else if (trav == NON) {
        if (ai.getAnnotations().isEmpty()) {
          break;
        }
      } else {
        Annotation a = getTravAnnotation();
        if (a == null) return -1;
        if (ai.getAnnotations().contains(a)) {
          break;
        }
      }
      row++;
    }
    return row;
  }

  private int[] getPrevFile() {
    int trav = getTraversal();
    int sel = sampleCombo.getSelectedIndex();
    sel--;
    int ind = 0;
    outer: while (sel >= 0) {
      String samp = sampleCombo.getItemAt(sel);
      if (trav == ALL) {
        break;
      } else {
        List<AnnotatedImage> imgs = new ArrayList<>();
        if (annotator.getAnnotationMap().containsKey(samp)) {
          imgs.addAll(annotator.getAnnotationMap().get(samp).values());
          sortAIs(imgs);
        }
        for (int i = imgs.size() - 1; i >= 0; i--) {
          AnnotatedImage ai = imgs.get(i);
          if (trav == ANN) {
            if (ai != null && ai.getAnnotations().size() > 0) {
              ind = i;
              break outer;
            }
          } else if (trav == NON) {
            if (ai == null || ai.getAnnotations().size() == 0) {
              ind = i;
              break outer;
            }
          } else {
            Annotation a = getTravAnnotation();
            if (ai != null && ai.getAnnotations().contains(a)) {
              ind = i;
              break outer;
            }
          }
        }
      }
      sel--;
    }
    return new int[] {sel, ind};
  }

  private int[] getNextFile() {
    int trav = getTraversal();
    int sel = sampleCombo.getSelectedIndex();
    sel++;
    int ind = 0;
    outer: while (sel < sampleCombo.getItemCount()) {
      String samp = sampleCombo.getItemAt(sel);
      if (trav == ALL) {
        break;
      } else {
        List<AnnotatedImage> imgs = new ArrayList<>();
        if (annotator.getAnnotationMap().containsKey(samp)) {
          imgs.addAll(annotator.getAnnotationMap().get(samp).values());
          sortAIs(imgs);
        }

        for (int i = 0; i < imgs.size(); i++) {
          AnnotatedImage ai = imgs.get(i);
          if (trav == ANN) {
            if (ai != null && ai.getAnnotations().size() > 0) {
              ind = i;
              break outer;
            }
          } else if (trav == NON) {
            if (ai == null || ai.getAnnotations().size() == 0) {
              ind = i;
              break outer;
            }
          } else {
            Annotation a = getTravAnnotation();
            if (ai != null && ai.getAnnotations().contains(a)) {
              ind = i;
              break outer;
            }
          }
        }
      }
      sel++;
    }
    return new int[] {sel, ind};
  }

  private void keyUp() { // prev node in tree
    int newRow = getPrevNodeRow();
    if (newRow < 0) return;
    tree.setSelectionRow(newRow);
    tree.scrollRowToVisible(newRow);
  }

  private int getTreeRowCount() {
    int cnt = getChildCount((DefaultMutableTreeNode) ((DefaultTreeModel) tree.getModel()).getRoot());
    return cnt;
  }

  private int getChildCount(DefaultMutableTreeNode parent) {
    int cnt = ((DefaultTreeModel) tree.getModel()).getChildCount(parent);
    Enumeration<?> e = parent.children();
    while (e.hasMoreElements()) {
      Object o = e.nextElement();
      cnt += getChildCount((DefaultMutableTreeNode) o);
    }
    return cnt;
  }

  private void keyDown() { // next node in tree
    int newRow = getNextNodeRow();
    if (newRow >= getTreeRowCount()) return;
    tree.setSelectionRow(newRow);
    tree.scrollRowToVisible(newRow);
  }

  private void keyLeft() { // prev file
    int[] prev = getPrevFile();
    if (prev[0] >= 0) {
      sampleCombo.setSelectedIndex(prev[0]);
    }
    tree.setSelectionRow(prev[1]);
    tree.scrollRowToVisible(prev[1]);
  }

  private void keyRight() { // next file
    int[] next = getNextFile();
    if (next[0] < sampleCombo.getItemCount()) {
      sampleCombo.setSelectedIndex(next[0]);
      tree.setSelectionRow(next[1]);
      tree.scrollRowToVisible(next[1]);
    }
  }

  private void updateTreeKeys(JTree tree) {
    tree.setUI(new BasicTreeUI() {

      protected KeyListener createKeyListener() {
        return new KeyAdapter() {

          @Override
          public void keyPressed(KeyEvent e) {
            e.consume();
            if (e.getKeyCode() == KeyEvent.VK_UP) {
              keyUp();
            } else if (e.getKeyCode() == KeyEvent.VK_DOWN) {
              keyDown();
            } else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
              keyLeft();
            } else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
              keyRight();
            } else {
              if (!(e.isAltDown() && e.getKeyCode() == KeyEvent.VK_F4)) {
                fireMnem(e);
              } else {
                close();
              }
            }
          }
        };
      }
    });
  }

  private void expandAllNodes(JTree tree) {
    int j = tree.getRowCount();
    int i = 0;
    while (i < j) {
      tree.expandRow(i);
      i += 1;
      j = tree.getRowCount();
    }
  }

  private void setSelectedNode(AnnotatedImage node) {
    BufferedImage img = node.getImage();
    selectedImage = img;
    imagePanel.repaint();
    lastSelectedGate = node.getName();
    lastSelectedFile = node.getImageFile();
  }

  private void sortAIs(List<AnnotatedImage> imgs) {
    imgs.sort(new Comparator<AnnotatedImage>() {

      @Override
      public int compare(AnnotatedImage o1, AnnotatedImage o2) {
        String n1 = o1.getName().substring(0, o1.getName().indexOf('~') - 1);
        String n2 = o2.getName().substring(0, o2.getName().indexOf('~') - 1);
        int[] chrPos1 = Positions.parseUCSClocation(n1);
        int[] chrPos2 = Positions.parseUCSClocation(n2);
        int comp;
        comp = Integer.compare(chrPos1[0], chrPos2[0]);
        if (comp != 0) return comp;
        comp = Integer.compare(chrPos1[1], chrPos2[1]);
        if (comp != 0) return comp;
        comp = Integer.compare(chrPos1[2], chrPos2[2]);
        if (comp != 0) return comp;
        comp = Integer.compare(Integer.parseInt(""
                                                + o1.getName().charAt(o1.getName().length() - 1)),
                               Integer.parseInt(""
                                                + o2.getName().charAt(o2.getName().length() - 1)));
        return comp;
      }
    });
  }

  private void updateAvail() {
    HashMap<String, HashMap<String, AnnotatedImage>> map = annotator.getAnnotationMap();
    String sampleName = (String) sampleCombo.getSelectedItem();
    HashMap<String, AnnotatedImage> ann = map.get(sampleName);

    List<AnnotatedImage> imgs = ann == null ? new ArrayList<>() : new ArrayList<>(ann.values());
    sortAIs(imgs);
    DefaultMutableTreeNode root = (DefaultMutableTreeNode) ((DefaultMutableTreeNode) tree.getModel()
                                                                                         .getRoot());
    root.removeAllChildren();
    for (AnnotatedImage img : imgs) {
      DefaultMutableTreeNode dmtn = new DefaultMutableTreeNode();
      dmtn.setUserObject(img);
      root.add(dmtn);
    }

    ((DefaultTreeModel) tree.getModel()).reload();
    expandAllNodes(tree);
    updateTreeKeys(tree);
    tree.repaint();
    tree.setSelectionRow(0);
  }

  private boolean saveAnnotations(boolean prompt) {
    boolean ask = prompt || (null == lastSavedAnnFile || "".equals(lastSavedAnnFile)
                             || !Files.exists(lastSavedAnnFile));
    String file = ask ? null : lastSavedAnnFile;
    if (file == null) {
      JFileChooser jfc = new JFileChooser(getLastUsedAnnotationDir());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save Annotations to File");
      int opt = jfc.showSaveDialog(frmAnnotator);
      if (opt == JFileChooser.APPROVE_OPTION) {
        file = jfc.getSelectedFile().getAbsolutePath();
      }
    }
    if (file == null) {
      return false;
    } else {
      setLastUsedAnnotationDir(ext.verifyDirFormat(ext.parseDirectoryOfFile(file)));
      annotator.saveAnnotations(file);
      setLastSavedAnnotationFile(file);
      lastTimeSaved = System.currentTimeMillis();
      checkAnnotationsChanged();
      if (Files.exists(BACKUP_FILE)) {
        new File(BACKUP_FILE).delete();
      }
      recentAnnotFiles.add(file);
      saveProperties();
      return true;
    }
  }

  private boolean promptToSaveAnnotations() {
    boolean showSave = !(null == lastSavedAnnFile || "".equals(lastSavedAnnFile)
                         || !Files.exists(lastSavedAnnFile));
    final String[] options = showSave ? new String[] {"Save As...", "Save", "Don't Save", "Cancel"}
                                      : new String[] {"Save As...", "Don't Save", "Cancel"};
    int result = JOptionPane.showOptionDialog(this.frmAnnotator, "Save Existing Annotations?",
                                              "Save?", JOptionPane.YES_NO_CANCEL_OPTION,
                                              JOptionPane.QUESTION_MESSAGE, null, options,
                                              "Cancel");
    if (showSave) {
      switch (result) {
        case 0:
          return saveAnnotations(true);
        case 1:
          return saveAnnotations(false);
        case 2:
          return true;
        case 3:
        default:
          return false;
      }
    } else {
      switch (result) {
        case 0:
          return saveAnnotations(true);
        case 1:
          return true;
        case 2:
        default:
          return false;
      }
    }
  }

  private void loadAnnotations() {
    if (this.annotator.getAnnotations().size() > 0) {
      if (!promptToSaveAnnotations()) {
        return;
      }
    }
    JFileChooser jfc = new JFileChooser(getLastUsedAnnotationDir());
    jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    jfc.setDialogTitle("Load Annotations from File");
    int opt = jfc.showOpenDialog(frmAnnotator);
    if (opt == JFileChooser.APPROVE_OPTION) {
      String annFile = jfc.getSelectedFile().getAbsolutePath();
      loadAnnotationFile(annFile);
    }
  }

  private void loadAnnotationFile(String annFile) {
    try {
      annotator = new Annotator();
      annotator.loadAnnotations(annFile);
      setLastSavedAnnotationFile(annFile);
      setLastUsedAnnotationDir(ext.verifyDirFormat(ext.parseDirectoryOfFile(annFile)));
      recentAnnotFiles.add(annFile);
      reloadControls();
      mnTravAnn.removeAll();
      for (Annotation a : annotator.getAnnotations()) {
        addAnnotationToTraversalMenu(a);
      }
      updateAvail();
      tree.setSelectionRow(0);
      tree.scrollRowToVisible(0);
      saveProperties();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  ActionListener comboListener = new ActionListener() {

    @Override
    public void actionPerformed(ActionEvent e) {
      tree.setEnabled(true);
      updateAvail();
    }
  };

  TreeSelectionListener treeListener = new TreeSelectionListener() {

    @Override
    public void valueChanged(TreeSelectionEvent e) {
      if (e.getNewLeadSelectionPath() == null) {
        refreshAnnotations();
        return;
      }
      Object selected = e.getNewLeadSelectionPath().getLastPathComponent();
      if (!(selected instanceof DefaultMutableTreeNode)) {
        System.err.println("Error - selected an invalid item somehow.  Likely programmer error.");
        return;
      }
      DefaultMutableTreeNode selectedNode = (DefaultMutableTreeNode) selected;
      AnnotatedImage ann = ((AnnotatedImage) selectedNode.getUserObject());
      setSelectedNode(ann);
      refreshAnnotations();
    }

  };

  private final Action loadAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      loadImageFiles();
    }
  };

  private final Action existAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      loadAnnotations();
    }
  };

  private final Action exportAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      new AnnotationExportDialog(annotator).setVisible(true);
      // export specific images with specific annotations
      // export sample names with specific annotations
    }
  };

  private final Action saveAsAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      saveAnnotations(true);
    }
  };
  private final Action saveAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      saveAnnotations(false);
    }
  };
  private JMenu mnSaveAnn;

}
