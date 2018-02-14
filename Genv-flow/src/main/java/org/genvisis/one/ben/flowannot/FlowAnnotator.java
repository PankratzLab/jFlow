package org.genvisis.one.ben.flowannot;

import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
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
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Properties;
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
import javax.swing.JRadioButton;
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
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.one.ben.flowannot.AnnotatedImage.Annotation;
import org.genvisis.one.ben.flowannot.IAnnotator.PANEL;
import net.miginfocom.swing.MigLayout;

public class FlowAnnotator {

  private JFrame frmFlowannotator;
  private JComboBox<String> fcsCombo;
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
  private JRadioButtonMenuItem jrbTravAll;
  private JRadioButtonMenuItem jrbTravAnn;
  private JRadioButtonMenuItem jrbTravNon;
  private ButtonGroup travGrp;
  private HashMap<Annotation, JRadioButtonMenuItem> annTravMap = new HashMap<>();

  private static final String PROP_FILE = ".flowannotator.properties";
  private static final String KEY_LAST_DIR_IMG = "LAST_DIR_IMG";
  private static final String KEY_LAST_DIR_ANN = "LAST_DIR_ANN";
  private static final String KEY_LAST_FILE_ANN = "LAST_FILE_ANN";
  private static final String KEY_LAST_GATE = "LAST_SEL_NODE";
  private static final String KEY_RECENT = "RECENT";

  private String lastOpenedImageDir = null;
  private String lastOpenedAnnFileDir = null;
  private String lastSavedAnnFile = null;
  private String lastSelectedGate = null;
  private volatile boolean constructingTree = false;
  private HashSet<String> recentAnnotFiles = new HashSet<>();

  private HashMap<Character, Action> mnemonicActions = new HashMap<>();

  private static final int ALL = 0;
  private static final int ANN = 1;
  private static final int NON = 2;
  private static final String BACKUP_DIR = ".backup./";
  private static final String BACKUP_FILE = BACKUP_DIR + "backup.annotations";
  private volatile boolean annotationsChanged = false;

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
          FlowAnnotator window = new FlowAnnotator();
          window.frmFlowannotator.setVisible(true);
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    });
  }

  /**
   * Create the application.
   */
  public FlowAnnotator() {
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
        this.lastOpenedAnnFileDir = ext.parseDirectoryOfFile(new File(FlowAnnotator.class.getProtectionDomain()
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
        this.lastOpenedImageDir = ext.parseDirectoryOfFile(new File(FlowAnnotator.class.getProtectionDomain()
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
    frmFlowannotator = new JFrame();
    frmFlowannotator.setTitle("FlowAnnotator");
    frmFlowannotator.setBounds(100, 100, 1200, 800);
    frmFlowannotator.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    frmFlowannotator.getContentPane().setLayout(new MigLayout("", "[grow]", "[grow]"));
    frmFlowannotator.addWindowListener(new WindowAdapter() {

      @Override
      public void windowClosing(WindowEvent e) {
        close();
      }
    });

    splitPane = new JSplitPane();
    splitPane.setResizeWeight(0.78);
    frmFlowannotator.getContentPane().add(splitPane, "cell 0 0,grow");

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
    splitPane_1.setResizeWeight(0);
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

    fcsCombo = new JComboBox<String>();
    fcsCombo.addActionListener(comboListener);

    constructingTree = true;
    ButtonGroup panelButtons = new ButtonGroup();
    rdBtnPanel1 = new JRadioButton();
    rdBtnPanel1.addItemListener(panelListener);
    rdBtnPanel1.setText("Panel 1");
    controlPanel.add(rdBtnPanel1, "flowx,cell 0 0,alignx center");
    rdBtnPanel1.setSelected(true);
    panelButtons.add(rdBtnPanel1);
    fcsCombo.setFocusable(false);
    controlPanel.add(fcsCombo, "flowx,cell 0 1,growx");

    JScrollPane scrollPane_1 = new JScrollPane();
    controlPanel.add(scrollPane_1, "cell 0 2,grow");

    tree = new JTree();
    scrollPane_1.setViewportView(tree);
    constructTree(PANEL.PANEL_1);

    rdBtnPanel2 = new JRadioButton();
    rdBtnPanel2.addItemListener(panelListener);
    rdBtnPanel2.setText("Panel 2");
    controlPanel.add(rdBtnPanel2, "cell 0 0,alignx center");
    panelButtons.add(rdBtnPanel2);
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0), "enter");
    am.put("enter", new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent arg0) {
        createNewAnnotation();
        tree.requestFocusInWindow();
      }
    });

    frmFlowannotator.setJMenuBar(createMenuBar());
    checkForBackupFileOrLoadLast();
    constructingTree = false;

    startAutoSaveThread();
  }

  private ItemListener panelListener = new ItemListener() {

    @Override
    public void itemStateChanged(ItemEvent e) {
      if (e.getStateChange() == ItemEvent.SELECTED) {
        PANEL panel = e.getSource() == rdBtnPanel1 ? PANEL.PANEL_1 : PANEL.PANEL_2;
        if (!constructingTree) {
          constructingTree = true;
          constructTree(panel);
          setAnnotationsChanged();
          reloadControls();
          saveProperties();
          constructingTree = false;
        }
      }
    }
  };

  private void constructTree(PANEL panel) {
    DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
    rootNode.add(GateTree.constructTree(panel));
    DefaultTreeModel dtm = new DefaultTreeModel(rootNode);
    tree.setModel(dtm);
    tree.setShowsRootHandles(true);
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
      int opt = JOptionPane.showConfirmDialog(this.frmFlowannotator, "Load Auto-saved backup file?",
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
        String gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tp.getLastPathComponent()).getUserObject()).getGateName();
        if (gate.equals(lastSelectedGate)) {
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
        PANEL p = FlowAnnotator.this.rdBtnPanel1.isSelected() ? PANEL.PANEL_1 : PANEL.PANEL_2;
        List<String> fcs = annotator.getFCSKeys(p);
        if (fcs.isEmpty()) {
          JOptionPane.showMessageDialog(FlowAnnotator.this.frmFlowannotator,
                                        "No files to search in panel " + (p == PANEL.PANEL_1 ? "1"
                                                                                             : "2")
                                                                             + "!");
          return;
        }
        String[] v = fcs.toArray(new String[fcs.size()]);
        List<String> sel = FileFinder.showFileFinder(v, false);
        if (!sel.isEmpty()) {
          fcsCombo.setSelectedItem(sel.get(0));
        }
      }
    });
    mntmNavFind.setMnemonic('F');
    mntmNavFind.setText("Find Sample");
    mnNav.add(mntmNavFind);

    JCheckBoxMenuItem mntmNavKeepGate = new JCheckBoxMenuItem();
    mntmNavKeepGate.setAction(new AbstractAction() {

      @Override
      public void actionPerformed(ActionEvent e) {
        keepGateWhenFileChange = mntmNavKeepGate.isSelected();
      }
    });
    mntmNavKeepGate.setText("Keep Gate When File Changes");
    mntmNavKeepGate.setSelected(true);
    mnNav.add(mntmNavKeepGate);

    return menuBar;
  }

  private void addAnnotationToTraversalMenu(Annotation annot) {
    JRadioButtonMenuItem jrb = new JRadioButtonMenuItem();
    jrb.setText(annot.annotation);
    annTravMap.put(annot, jrb);
    travGrp.add(jrb);
    mnTravAnn.add(jrb);
    JMenuItem jmn = new JMenuItem();
    jmn.setAction(saveAnnotationAction);
    jmn.setText(annot.annotation);
    mnSaveAnn.add(jmn);
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
      int opt = jfc.showSaveDialog(frmFlowannotator);
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
    frmFlowannotator.setVisible(false);
    frmFlowannotator.dispose();
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
    int opt = jfc.showOpenDialog(frmFlowannotator);
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

  private void reloadControls() {
    int panels = 0;
    if (rdBtnPanel1.isSelected()) {
      panels = 1;
    } else if (rdBtnPanel2.isSelected()) {
      panels = 2;
    }
    String[] fcsKeys;
    switch (panels) {
      case 0:
        fcsKeys = new String[0];
        break;
      case 1:
        fcsKeys = annotator.getFCSKeys(PANEL.PANEL_1).toArray(new String[0]);
        break;
      case 2:
        fcsKeys = annotator.getFCSKeys(PANEL.PANEL_2).toArray(new String[0]);
        break;
      case 3:
      default:
        fcsKeys = annotator.getFCSKeys().toArray(new String[0]);
        break;
    }
    DefaultComboBoxModel<String> dcbm = new DefaultComboBoxModel<>(fcsKeys);
    fcsCombo.setModel(dcbm);
    if (fcsCombo.getModel().getSize() > 0) {
      fcsCombo.setSelectedIndex(0);
    }
    tree.requestFocusInWindow();
  }

  private void refreshAnnotations() {
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
      if (ai != null && ai.getAnnotations().contains(allAnnots.get(i))) {
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
  }

  private void fireMnem(char code) {
    if (mnemonicActions.containsKey(code)) {
      mnemonicActions.get(code).actionPerformed(null);
      setAnnotationsChanged();
    }
  }

  boolean keepGateWhenFileChange = true;

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

  private int getPrevFile() {
    int trav = getTraversal();
    int sel = fcsCombo.getSelectedIndex();
    String gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getPathForRow(0)
                                                                  .getLastPathComponent()).getUserObject()).getGateName();
    // default to root, or last leaf?, or any ann/non?
    if (keepGateWhenFileChange && !tree.isSelectionEmpty()) {
      gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getSelectionPath()
                                                             .getLastPathComponent()).getUserObject()).getGateName();
    }
    sel--;
    while (sel >= 0) {
      String fcs = fcsCombo.getItemAt(sel);
      AnnotatedImage ai = annotator.getAnnotationMap().get(fcs).get(gate);
      if (trav == ALL) {
        break;
      } else if (trav == ANN) {
        if (ai != null && ai.getAnnotations().size() > 0) {
          break;
        }
      } else if (trav == NON) {
        if (ai == null || ai.getAnnotations().size() == 0) {
          break;
        }
      } else {
        Annotation a = getTravAnnotation();
        if (a == null) return -1;
        if (ai != null && ai.getAnnotations().contains(a)) {
          break;
        }
      }
      sel--;
    }
    return sel;
  }

  private int getNextFile() {
    int trav = getTraversal();
    int sel = fcsCombo.getSelectedIndex();
    String gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getPathForRow(0)
                                                                  .getLastPathComponent()).getUserObject()).getGateName();
    // default to root, or last leaf?, or any ann/non?
    if (keepGateWhenFileChange && !tree.isSelectionEmpty()) {
      gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getSelectionPath()
                                                             .getLastPathComponent()).getUserObject()).getGateName();
    }
    sel++;
    while (sel < fcsCombo.getItemCount()) {
      String fcs = fcsCombo.getItemAt(sel);
      if (trav == ALL) {
        break;
      } else if (trav == ANN) {
        if (annotator.getAnnotationMap().get(fcs).get(gate).getAnnotations().size() > 0) {
          break;
        }
      } else if (trav == NON) {
        if (annotator.getAnnotationMap().get(fcs).get(gate).getAnnotations().size() == 0) {
          break;
        }
      } else {
        Annotation a = getTravAnnotation();
        if (a == null) return -1;
        if (annotator.getAnnotationMap().get(fcs).get(gate).getAnnotations().contains(a)) {
          break;
        }
      }
      sel++;
    }
    return sel;
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
    int prev = getPrevFile();
    int[] rows = tree.getSelectionRows();
    int row = -1;
    if (rows != null && rows.length > 0) {
      row = rows[0];
    }
    if (prev >= 0) {
      fcsCombo.setSelectedIndex(prev);
    }
    tree.setSelectionRow(keepGateWhenFileChange && row != -1 ? row : 0);
  }

  private void keyRight() { // next file
    int next = getNextFile();
    int[] rows = tree.getSelectionRows();
    int row = -1;
    if (rows != null && rows.length > 0) {
      row = rows[0];
    }
    if (next < fcsCombo.getItemCount()) {
      fcsCombo.setSelectedIndex(next);
      tree.setSelectionRow(keepGateWhenFileChange && row != -1 ? row : 0);
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
                fireMnem((e.getKeyChar() + "").toUpperCase().charAt(0));
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
    refreshAnnotations();
    if (!constructingTree) {
      lastSelectedGate = node.getGateName();
    }
  }

  private void updateAvail() {
    HashMap<String, HashMap<String, AnnotatedImage>> map = annotator.getAnnotationMap();
    String fcsFile = (String) fcsCombo.getSelectedItem();
    HashMap<String, AnnotatedImage> ann = map.get(fcsFile);
    DefaultMutableTreeNode root = (DefaultMutableTreeNode) ((DefaultMutableTreeNode) tree.getModel()
                                                                                         .getRoot()).getChildAt(0);
    updateNode(root, ann);
    ((DefaultTreeModel) tree.getModel()).reload();
    expandAllNodes(tree);
    updateTreeKeys(tree);
    tree.repaint();
  }

  private void updateNode(DefaultMutableTreeNode node, HashMap<String, AnnotatedImage> annMap) {
    AnnotatedImage ai = (AnnotatedImage) node.getUserObject();
    AnnotatedImage newAi = annMap.get(ai.getGateName());
    if (newAi == null) {
      newAi = new AnnotatedImage(ai.getGateName(), ai.isRoot());
      newAi.setMissing(true);
    }
    node.setUserObject(newAi);
    for (int i = 0; i < node.getChildCount(); i++) {
      updateNode((DefaultMutableTreeNode) node.getChildAt(i), annMap);
    }
  }

  private boolean saveAnnotations(boolean prompt) {
    boolean ask = prompt || (null == lastSavedAnnFile || "".equals(lastSavedAnnFile)
                             || !Files.exists(lastSavedAnnFile));
    String file = ask ? null : lastSavedAnnFile;
    if (file == null) {
      JFileChooser jfc = new JFileChooser(getLastUsedAnnotationDir());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save Annotations to File");
      int opt = jfc.showSaveDialog(frmFlowannotator);
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
    int result = JOptionPane.showOptionDialog(this.frmFlowannotator, "Save Existing Annotations?",
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
    int opt = jfc.showOpenDialog(frmFlowannotator);
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
      for (Annotation a : annotator.getAnnotations()) {
        addAnnotationToTraversalMenu(a);
      }
      updateAvail();
      tree.setSelectionRow(0);
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
      if (e.getNewLeadSelectionPath() == null) return;
      Object selected = e.getNewLeadSelectionPath().getLastPathComponent();
      if (!(selected instanceof DefaultMutableTreeNode)) {
        System.err.println("Error - selected an invalid item somehow.  Likely programmer error.");
        return;
      }
      DefaultMutableTreeNode selectedNode = (DefaultMutableTreeNode) selected;
      AnnotatedImage ann = ((AnnotatedImage) selectedNode.getUserObject());
      setSelectedNode(ann);
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
  private JRadioButton rdBtnPanel1;
  private JRadioButton rdBtnPanel2;

}
