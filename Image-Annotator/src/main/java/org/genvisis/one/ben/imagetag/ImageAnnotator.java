package org.genvisis.one.ben.imagetag;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
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
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;
import javax.swing.border.BevelBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.plaf.basic.BasicTreeUI;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;

import org.apache.commons.collections4.MultiSet;
import org.apache.commons.collections4.multiset.HashMultiSet;
import org.genvisis.one.ben.imagetag.AnnotatedImage.Annotation;
import org.genvisis.one.ben.imagetag.IAnnotator.Loaded;
import org.genvisis.one.ben.imagetag.IAnnotator.LoadedSaved;
import org.genvisis.seq.manage.BEDFileReader;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Images;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.Sets;

import net.miginfocom.swing.MigLayout;

public class ImageAnnotator {

  private JFrame frmAnnotator;
  private JComboBox<String> sampleCombo;
  private JTree tree;

  IAnnotator annotator = new Annotator();

  BufferedImage selectedImage = createReadyImage();

  private JMenu mnLoadRecent;
  private JTextField fldNewAnnot;
  private JPanel annotPanel;
  private JPanel imagePanel;
  private JPanel controlPanel;
  private JSplitPane splitPane;
  private JMenu mnTravAnn;
  private JMenu mnDelAnn;
  private JMenu mnSaveAnn;
  private JLabel validationLabel;
  private JLabel validationCountLabel;
  private JMenu mnNavHide;
  private JRadioButtonMenuItem jrbTravAll;
  private JRadioButtonMenuItem jrbTravAnn;
  private JRadioButtonMenuItem jrbTravNon;
  private ButtonGroup travGrp;
  private HashMap<Annotation, JRadioButtonMenuItem> annTravMap = new HashMap<>();
  private HashMap<Annotation, JMenuItem> annDelMap = new HashMap<>();
  private HashMap<Annotation, JMenuItem> annSaveMap = new HashMap<>();
  private HashMap<Annotation, JCheckBoxMenuItem> annHideMap = new HashMap<>();
  private HashSet<Annotation> hiddenAnnotations = new HashSet<>();
  private SORT_TYPE sortOrder = SORT_TYPE.GENOMIC_POSITION;

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

  private static final String BACKUP_DIR = ".backup./";
  private static final String BACKUP_FILE = BACKUP_DIR + "backup.annotations";
  private volatile boolean annotationsChanged = false;
  private volatile boolean checkBxAdvanceAfterAnnotationKeyed = true;
  private volatile boolean checkBxHideMissing = false;

  private void setAnnotationsChanged() {
    annotationsChanged = true;
  }

  private boolean checkAnnotationsChanged() {
    boolean value = annotationsChanged;
    annotationsChanged = false;
    return value;
  }

  /** Launch the application. */
  public static void launch() {

    // set system-wide anti-aliasing
    System.setProperty("awt.useSystemAAFontSettings", "on");
    System.setProperty("swing.aatext", "true");

    ToolTipManager.sharedInstance().setInitialDelay(0);
    ToolTipManager.sharedInstance().setDismissDelay(Integer.MAX_VALUE - 1);
    ToolTipManager.sharedInstance().setReshowDelay(0);
    UIManager.put("ToolTip.background", Color.decode("#F5F5DC"));

    Runnable r = new Runnable() {

      public void run() {
        try {
          ImageAnnotator window = new ImageAnnotator();
          window.frmAnnotator.setVisible(true);
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    };

    EventQueue.invokeLater(r);
  }

  public static void main(String[] args) {
    launch();
  }

  private boolean showDebug = PRELOAD_DEBUG;
  private static final boolean PRELOAD_DEBUG = false;
  private static final int WINDOW_SPEED = 4;
  private static final int WINDOW_BEHAVIOR = 6;
  private static final int WINDOW_DIRECTION = 4;
  private static final int PRELOAD_INIT_WAIT = 3000;
  private static final int PRELOAD_INTERVAL = 1700;
  Timer cacheBalancingTimer = new Timer();

  /** Create the application. */
  public ImageAnnotator() {
    loadProperties();
    initialize();

    cacheBalancingTimer.scheduleAtFixedRate(new TimerTask() {

      @Override
      public void run() {
        SPEED recentSpeed = getRecentUserSpeed(WINDOW_SPEED);
        BEHAVIOR recentBehavior = getRecentUserBehavior(WINDOW_BEHAVIOR);
        DIRECTION recentDirection = getRecentUserDirection(WINDOW_DIRECTION);

        int preload = Math.max(0, ((2 * recentSpeed.modifier) - recentBehavior.adjustment))
                      * recentDirection.sign;

        int curr = tree.getRowForPath(tree.getSelectionPath());
        int treeMax = getTreeRowCount();
        int currSample = sampleCombo.getSelectedIndex();
        int sampMax = sampleCombo.getItemCount();

        int used = 0;
        if (recentDirection.sign == -1) {
          int ind = curr;
          for (int i = curr - 1; i >= Math.max(0, curr + preload); i--) {
            ind = getPrevNodeRow(ind);
            if (ind < 0) break; // TODO preload first images for prev sample
            AnnotatedImage ann = getAnnotatedImage(ind);
            imageCache.getUnchecked(ann);
            used++;
          }
        } else if (recentDirection.sign == 1) {
          int ind = curr;
          for (int i = curr + 1; i <= Math.min(curr + preload, treeMax - 1); i++) {
            ind = getNextNodeRow(ind);
            if (ind >= getTreeRowCount()) break; // TODO preload first images for next sample
            AnnotatedImage ann = getAnnotatedImage(ind);
            imageCache.getUnchecked(ann);
            used++;
          }
        }
        if (showDebug) {
          System.out.println("Preload: " + preload + " | Actual: " + used + " -- basis: "
                             + recentSpeed + ", " + recentBehavior + ", " + recentDirection);
        }
        tree.repaint();
      }
    }, PRELOAD_INIT_WAIT, PRELOAD_INTERVAL);
  }

  CacheLoader<AnnotatedImage, BufferedImage> loader = new CacheLoader<AnnotatedImage, BufferedImage>() {

    @Override
    public BufferedImage load(AnnotatedImage key) throws Exception {
      return createImage(key);
    }
  };
  LoadingCache<AnnotatedImage, BufferedImage> imageCache = CacheBuilder.newBuilder().softValues()
                                                                       .initialCapacity(10)
                                                                       .expireAfterAccess(30,
                                                                                          TimeUnit.SECONDS)
                                                                       .maximumSize(15)
                                                                       .build(loader);

  public BufferedImage getImage(AnnotatedImage ai) {
    return imageCache.getUnchecked(ai);
  }

  private BufferedImage createImage(AnnotatedImage ai) {
    BufferedImage image = null;
    if (ai.getImageFile() != null) {
      if (!ai.getImageFile().contains(";")) {
        try {
          image = ImageIO.read(new File(ai.getImageFile()));
        } catch (IOException e) {
          e.printStackTrace();
          image = createIOExceptionImage(e);
        }
      } else {
        String[] images = ai.getImageFile().split(";");
        image = Images.stitchImages(images, Color.WHITE, false, false);
      }
    } else {
      if (ai.getImageFile() == null) {
        image = createNoFileImage();
      } else if (!Files.exists(ai.getImageFile())) {
        image = createMissingFileImage(ai.getImageFile());
      }
    }
    return image;
  }

  static float fontSize = 16f;

  private static BufferedImage createImage(String msg, String msg2) {
    BufferedImage bi = new BufferedImage(800, 600, BufferedImage.TYPE_INT_RGB);
    Graphics g = bi.getGraphics();
    g.setColor(UIManager.getColor("Panel.background"));
    g.fillRect(0, 0, bi.getWidth(), bi.getHeight());
    g.setColor(Color.BLACK);
    g.setFont(g.getFont().deriveFont(fontSize));
    FontMetrics fm = g.getFontMetrics();
    g.drawString(msg, (bi.getWidth() / 2) - fm.stringWidth(msg) / 2,
                 bi.getHeight() / 2 - fm.getHeight());
    if (msg2 != null) {
      g.drawString(msg2, (bi.getWidth() / 2) - fm.stringWidth(msg2) / 2,
                   bi.getHeight() / 2 + ((int) (fm.getHeight() * 1.5)));
    }
    return bi;
  }

  private BufferedImage createIOExceptionImage(IOException e) {
    return createImage("Exception when loading image:", e.getMessage());
  }

  public static BufferedImage createReadyImage() {
    return createImage("Load Image Directory to Begin.", null);
  }

  private BufferedImage createNoFileImage() {
    return createImage("No image file found!", null);
  }

  private BufferedImage createMissingFileImage(String file) {
    return createImage("File missing:", file);
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

  /** Initialize the contents of the frame. */
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
      private static final long serialVersionUID = 1L;

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
    controlPanel.setLayout(new MigLayout("hidemode 3", "[110.00,grow]", "[][][grow][][]"));

    validationLabel = new JLabel("Validation Results:");
    validationLabel.setFocusable(false);
    validationLabel.setVisible(false);
    controlPanel.add(validationLabel, "cell 0 0, split 2, growx");
    validationCountLabel = new JLabel();
    validationCountLabel.setHorizontalAlignment(SwingConstants.RIGHT);
    validationCountLabel.setFocusable(false);
    validationCountLabel.setVisible(false);
    controlPanel.add(validationCountLabel, "cell 0 0, growx, alignx right");

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
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent arg0) {
        createNewAnnotation();
        tree.requestFocusInWindow();
      }
    });

    frmAnnotator.setJMenuBar(createMenuBar());
    checkForBackupFileOrLoadLast();

    startAutoSaveThread();
    startInteractionProcessingThread();
  }

  public static final String CHR_POS_REGEX = "(.*)_?chr([12]?[0-9[XYM]]+)[:-_]([\\d,]+)[-_]([\\d,]+).*?";

  private int[] parseSampleChrPosIfExists(String nameOrFile) {
    Matcher m = Pattern.compile(CHR_POS_REGEX).matcher(nameOrFile);
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
    tree.setCellRenderer(new DefaultTreeCellRenderer() {
      private static final long serialVersionUID = 1L;

      @Override
      public Component getTreeCellRendererComponent(JTree tree, Object value, boolean sel,
                                                    boolean expanded, boolean leaf, int row,
                                                    boolean hasFocus) {
        Component comp = super.getTreeCellRendererComponent(tree, value, sel, expanded, leaf, row,
                                                            hasFocus);
        Object v = ((DefaultMutableTreeNode) value).getUserObject();
        if (showDebug && v != null && v instanceof AnnotatedImage) {
          comp.setForeground(imageCache.getIfPresent(((AnnotatedImage) ((DefaultMutableTreeNode) value).getUserObject())) == null ? Color.RED
                                                                                                                                  : Color.CYAN);
        } else {
          comp.setForeground(Color.BLACK);
        }
        return comp;
      }
    });
    final JPopupMenu menu = new JPopupMenu();
    tree.addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        super.mousePressed(e);
        if (SwingUtilities.isRightMouseButton(e)) {
          int row = tree.getRowForLocation(e.getX(), e.getY());
          TreePath path = tree.getPathForRow(row);
          DefaultMutableTreeNode last = (DefaultMutableTreeNode) path.getLastPathComponent();
          AnnotatedImage ai = (AnnotatedImage) last.getUserObject();
          String file = ai.getImageFile();

          menu.removeAll();

          addContextItemsToMenu(menu, file);

          menu.add(new JSeparator(SwingConstants.HORIZONTAL));

          JMenuItem closeItem = new JMenuItem();
          closeItem.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent arg0) {
              menu.setVisible(false);
            }
          });
          closeItem.setText("Close");
          menu.add(closeItem);

          menu.show(e.getComponent(), e.getX(), e.getY());
        }
      }

    });
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

  private void startInteractionProcessingThread() {
    interactionProcessingThread.setDaemon(true);
    interactionProcessingThread.start();
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

    JMenuItem mntmValidate = new JMenuItem();
    mntmValidate.setAction(validateAction);
    mntmValidate.setText("Validate Against BED File");
    mntmValidate.setMnemonic('V');
    mnFile.add(mntmValidate);

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    addContextItemsToMenu(mnFile, lastSelectedFile);

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    mnDelAnn = new JMenu("Delete Annotation:");
    mnFile.add(mnDelAnn);

    mnFile.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmExit = new JMenuItem();
    mntmExit.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

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

    JMenuItem mntmSortByLbl = new JMenuItem();
    mntmSortByLbl.setEnabled(false);
    mntmSortByLbl.setText("Sort By:");
    mnNav.add(mntmSortByLbl);

    ButtonGroup sortGroup = new ButtonGroup();
    for (SORT_TYPE sort : SORT_TYPE.values()) {
      JRadioButtonMenuItem mntmSort = new JRadioButtonMenuItem();
      sortGroup.add(mntmSort);
      mntmSort.setAction(new AbstractAction() {
        private static final long serialVersionUID = 1L;

        @Override
        public void actionPerformed(ActionEvent e) {
          ImageAnnotator.this.sortOrder = sort;
          updateAvail();
        }
      });
      mntmSort.setText("Sort by " + sort.dispName);
      if (sort == sortOrder) {
        mntmSort.setSelected(true);
      }
      mnNav.add(mntmSort);
    }

    mnNav.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmNavOptLbl = new JMenuItem();
    mntmNavOptLbl.setEnabled(false);
    mntmNavOptLbl.setText("Other Options:");
    mnNav.add(mntmNavOptLbl);

    JMenuItem mntmNavFind = new JMenuItem();
    mntmNavFind.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

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
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        checkBxAdvanceAfterAnnotationKeyed = mntmNavAdvAfterAnnot.isSelected();
      }
    });
    mntmNavAdvAfterAnnot.setText("Advance after annotating (or SHIFT)");
    mntmNavAdvAfterAnnot.setSelected(true);
    mnNav.add(mntmNavAdvAfterAnnot);

    mnNavHide = new JMenu();
    mnNavHide.setText("Hide:");
    mnNav.add(mnNavHide);

    JCheckBoxMenuItem mntmNavHideMissing = new JCheckBoxMenuItem();
    mntmNavHideMissing.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        checkBxHideMissing = mntmNavHideMissing.isSelected();
        updateAvail();
      }
    });
    mntmNavHideMissing.setText("Hide Missing");
    mnNavHide.add(mntmNavHideMissing);
    mnNavHide.add(new JSeparator(SwingConstants.HORIZONTAL));

    return menuBar;
  }

  private void addContextItemsToMenu(JComponent mnFile, String file) {
    JMenuItem mntmCopyUCSC = new JMenuItem();
    mntmCopyUCSC.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(Positions.getUCSCformat(parseSampleChrPosIfExists(file == null ? lastSelectedFile
                                                                                        : file)));
      }
    });
    mntmCopyUCSC.setText("Copy UCSC Location to Clipboard");
    mntmCopyUCSC.setMnemonic('U');
    mnFile.add(mntmCopyUCSC);

    JMenuItem mntmCopyUCSCLink = new JMenuItem();
    mntmCopyUCSCLink.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(Positions.getUCSClink(parseSampleChrPosIfExists(file == null ? lastSelectedFile
                                                                                      : file)));
      }
    });
    mntmCopyUCSCLink.setText("Copy UCSC Link to Clipboard");
    mntmCopyUCSCLink.setMnemonic('L');
    mnFile.add(mntmCopyUCSCLink);

    JMenuItem mntmCopySample = new JMenuItem();
    mntmCopySample.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

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
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        ext.setClipboard(file == null ? lastSelectedFile : file);
      }
    });
    mntmCopySampleFile.setText("Copy File Path to Clipboard");
    mntmCopySampleFile.setMnemonic('N');
    mnFile.add(mntmCopySampleFile);
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
    annSaveMap.put(annot, jmn);
    JMenuItem jmd = new JMenuItem();
    jmd.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        deleteAnnotation(annot);
      }
    });
    jmd.setText(annot.annotation);
    mnDelAnn.add(jmd);
    annDelMap.put(annot, jmd);
    mnSaveAnn.add(jmn);

    JCheckBoxMenuItem mntmHideAnn = new JCheckBoxMenuItem();
    mntmHideAnn.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        if (mntmHideAnn.isSelected()) {
          hiddenAnnotations.add(annot);
        } else {
          hiddenAnnotations.remove(annot);
        }
        updateAvail();
      }
    });
    mntmHideAnn.setText(annot.annotation);
    annHideMap.put(annot, mntmHideAnn);
    mnNavHide.add(mntmHideAnn);
  }

  private void deleteAnnotation(Annotation annot) {
    JRadioButtonMenuItem jrb = annTravMap.get(annot);
    travGrp.remove(jrb);
    mnTravAnn.remove(jrb);
    mnSaveAnn.remove(annSaveMap.get(annot));
    mnDelAnn.remove(annDelMap.get(annot));
    mnNavHide.remove(annHideMap.get(annot));
    hiddenAnnotations.remove(annot);
    annotator.deleteAnnotation(annot);

    setAnnotationsChanged();
    refreshAnnotations();
    updateAvail();
  }

  private AbstractAction saveAnnotationAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

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
    process = false;
    cacheBalancingTimer.cancel();
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
          private static final long serialVersionUID = 1L;

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
      validationFiles.clear();
      annotator = new Annotator();
      Loaded loaded = annotator.loadImgDir(fS);
      if (loaded.unmatched.size() > 0) {
        JLabel label = new JLabel("<html>" + loaded.getString() + "</html>");
        JScrollPane pane = new JScrollPane(label);
        JPanel panel = new JPanel(new BorderLayout());
        panel.add(pane, BorderLayout.CENTER);
        panel.add(new JLabel("<html>Image names must start with the name of the subdirectory they are in.</html>"),
                  BorderLayout.SOUTH);
        JOptionPane.showMessageDialog(ImageAnnotator.this.frmAnnotator, panel, "Unmatched Files",
                                      JOptionPane.WARNING_MESSAGE);
      }
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
    String[] keys = annotator.getRoots().toArray(new String[0]);
    // Comparator adapted from
    // https://stackoverflow.com/questions/41085394/java-comparator-alphanumeric-strings
    Arrays.sort(keys, Comparator.comparingLong(s -> hasNumeric(
                                                               // if there are any digits in the
                                                               // string, delete all non-digit
                                                               // characters and compare the
                                                               // resulting values
                                                               s.toString()) ? Long.parseLong(s.toString().replaceAll("\\D", "")) : 0)
                                // then delete all digit characters and compare the resulting
                                // strings
                                .thenComparing(s -> s.toString().replaceAll("\\d", "")));
    DefaultComboBoxModel<String> dcbm = new DefaultComboBoxModel<>(keys);
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
        private static final long serialVersionUID = 1L;

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

  private TRAVERSAL getTraversal() {
    if (jrbTravAll.isSelected()) {
      return TRAVERSAL.ALL;
    } else if (jrbTravAnn.isSelected()) {
      return TRAVERSAL.ANN;
    } else if (jrbTravNon.isSelected()) {
      return TRAVERSAL.NON;
    }
    return TRAVERSAL.CUSTOM;
  }

  private int getPrevNodeRow(int start) {
    int row = start;
    TRAVERSAL trav = getTraversal();
    row--;
    while (row >= 0) {
      AnnotatedImage ai = getAnnotatedImage(row);
      if (trav == TRAVERSAL.ALL) {
        break;
      } else if (trav == TRAVERSAL.ANN) {
        if (!ai.getAnnotations().isEmpty()) {
          break;
        }
      } else if (trav == TRAVERSAL.NON) {
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

  private int getPrevNodeRow() {
    if (tree.isSelectionEmpty()) {
      return -1;
    }
    int[] rows = tree.getSelectionRows();
    if (rows == null || rows.length == 0) {
      return -1;
    }
    int row = rows[0];
    return getPrevNodeRow(row);
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

  private int getNextNodeRow(int start) {
    int row = start;
    TRAVERSAL trav = getTraversal();
    int cnt = getTreeRowCount();
    row++;
    while (row < cnt) {
      AnnotatedImage ai = getAnnotatedImage(row);
      if (trav == TRAVERSAL.ALL) {
        break;
      } else if (trav == TRAVERSAL.ANN) {
        if (!ai.getAnnotations().isEmpty()) {
          break;
        }
      } else if (trav == TRAVERSAL.NON) {
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

  private AnnotatedImage getAnnotatedImage(int row) {
    DefaultTreeModel model = (DefaultTreeModel) tree.getModel();
    return (AnnotatedImage) ((DefaultMutableTreeNode) model.getChild(model.getRoot(),
                                                                     row)).getUserObject();
  }

  private int getNextNodeRow() {
    if (tree.isSelectionEmpty()) {
      return -1;
    }
    int[] rows = tree.getSelectionRows();
    if (rows == null || rows.length == 0) {
      return -1;
    }
    int row = rows[0];
    return getNextNodeRow(row);
  }

  private int[] getPrevFile() {
    int sel = sampleCombo.getSelectedIndex();
    return getPrevFile(sel);
  }

  private int[] getPrevFile(int start) {
    int sel = start;
    TRAVERSAL trav = getTraversal();
    sel--;
    int ind = 0;
    outer: while (sel >= 0) {
      String samp = sampleCombo.getItemAt(sel);
      if (trav == TRAVERSAL.ALL) {
        break;
      } else {
        List<AnnotatedImage> imgs = new ArrayList<>();
        if (annotator.getAnnotationMap().containsKey(samp)) {
          imgs.addAll(annotator.getAnnotationMap().get(samp).values());
          sortAIs(imgs);
        }
        for (int i = imgs.size() - 1; i >= 0; i--) {
          AnnotatedImage ai = imgs.get(i);
          if (trav == TRAVERSAL.ANN) {
            if (ai != null && ai.getAnnotations().size() > 0) {
              ind = i;
              break outer;
            }
          } else if (trav == TRAVERSAL.NON) {
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
    int sel = sampleCombo.getSelectedIndex();
    return getNextFile(sel);
  }

  private int[] getNextFile(int start) {
    int sel = start;
    TRAVERSAL trav = getTraversal();
    sel++;
    int ind = 0;
    outer: while (sel < sampleCombo.getItemCount()) {
      String samp = sampleCombo.getItemAt(sel);
      if (trav == TRAVERSAL.ALL) {
        break;
      } else {
        List<AnnotatedImage> imgs = new ArrayList<>();
        if (annotator.getAnnotationMap().containsKey(samp)) {
          imgs.addAll(annotator.getAnnotationMap().get(samp).values());
          sortAIs(imgs);
        }

        for (int i = 0; i < imgs.size(); i++) {
          AnnotatedImage ai = imgs.get(i);
          if (trav == TRAVERSAL.ANN) {
            if (ai != null && ai.getAnnotations().size() > 0) {
              ind = i;
              break outer;
            }
          } else if (trav == TRAVERSAL.NON) {
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

  private void keyUp() { // prev node in tree
    int newRow = getPrevNodeRow();
    if (newRow < 0) return;
    boolean stt = startInteractionEvent();
    if (stt) {
      buildingEvent.direction = DIRECTION.PREV_IMG;
    }
    tree.setSelectionRow(newRow);
    tree.scrollRowToVisible(newRow);
    if (stt) {
      fireInteractionEvent();
    }
  }

  private void keyDown() { // next node in tree
    int newRow = getNextNodeRow();
    if (newRow >= getTreeRowCount()) return;
    boolean stt = startInteractionEvent();
    if (stt) {
      buildingEvent.direction = DIRECTION.NEXT_IMG;
    }
    tree.setSelectionRow(newRow);
    tree.scrollRowToVisible(newRow);
    if (stt) {
      fireInteractionEvent();
    }
  }

  private void keyLeft() { // prev file
    int[] prev = getPrevFile();
    boolean stt = startInteractionEvent();
    if (prev[0] >= 0) {
      if (stt) {
        buildingEvent.direction = DIRECTION.PREV_SMP;
      }
      sampleCombo.setSelectedIndex(prev[0]);
    } else {
      if (stt) {
        buildingEvent.direction = DIRECTION.PREV_IMG;
      }
    }
    tree.setSelectionRow(prev[1]);
    tree.scrollRowToVisible(prev[1]);
    if (stt) {
      fireInteractionEvent();
    }
  }

  private void keyRight() { // next file
    int[] next = getNextFile();
    if (next[0] < sampleCombo.getItemCount()) {
      boolean stt = startInteractionEvent();
      if (stt) {
        buildingEvent.direction = DIRECTION.NEXT_SMP;
      }
      sampleCombo.setSelectedIndex(next[0]);
      tree.setSelectionRow(next[1]);
      tree.scrollRowToVisible(next[1]);
      if (stt) {
        fireInteractionEvent();
      }
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

  private static enum TRAVERSAL {
    ALL, ANN, NON, CUSTOM;
  }

  static enum DIRECTION {
    NEXT_IMG(1), PREV_IMG(-1), NEXT_SMP(1), PREV_SMP(-1), SKIP(0);

    private DIRECTION(int sig) {
      this.sign = sig;
    }

    int sign;
  }

  static enum SPEED {
    REALLY_FAST(8), FAST(5), MODERATE(3), SLOW(1);

    private SPEED(int mod) {
      this.modifier = mod;
    }

    int modifier;
  }

  static enum BEHAVIOR {
    LINEAR(0), // little to no retrieval of past images
    UNSURE(2), // moderate rate of retrieval of past images
    REVIEW(4), // high rate of retrieval of past images
    SEARCH(10); // no pattern of retrieval / not linear

    private BEHAVIOR(int adj) {
      adjustment = adj;
    }

    int adjustment;
  }

  private void setSelectedNode(AnnotatedImage node) {
    BufferedImage img = getImage(node);
    selectedImage = img;
    imagePanel.repaint();
    lastSelectedGate = node.getName();
    lastSelectedFile = node.getImageFile();
  }

  private static enum SORT_TYPE {
    GENOMIC_POSITION("Genomic Position", new Comparator<AnnotatedImage>() {

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
    }), SEGMENT_SIZE("Segment Size", new Comparator<AnnotatedImage>() {

      @Override
      public int compare(AnnotatedImage o1, AnnotatedImage o2) {
        String n1 = o1.getName().substring(0, o1.getName().indexOf('~') - 1);
        String n2 = o2.getName().substring(0, o2.getName().indexOf('~') - 1);
        int[] chrPos1 = Positions.parseUCSClocation(n1);
        int[] chrPos2 = Positions.parseUCSClocation(n2);
        int size1 = chrPos1[2] - chrPos1[1];
        int size2 = chrPos2[2] - chrPos2[1];
        return Integer.compare(size1, size2);
      }
    });

    private SORT_TYPE(String displayName, Comparator<AnnotatedImage> aiComp) {
      this.dispName = displayName;
      this.comparator = aiComp;
    }

    final String dispName;
    final Comparator<AnnotatedImage> comparator;
  }

  private void sortAIs(List<AnnotatedImage> imgs) {
    imgs.sort(sortOrder.comparator);
  }

  class ValidationResult {
    final int countBedSegs;
    final int countFiles;
    final int countOverlap;
    final Set<Segment> extraSegsBed;
    final Set<Segment> extraSegsFile;

    public ValidationResult(int cntBedSegs, int cntFiles, int cntOverlap, Set<Segment> xtraSegsBed,
                            Set<Segment> xtraSegsFile) {
      this.countBedSegs = cntBedSegs;
      this.countFiles = cntFiles;
      this.countOverlap = cntOverlap;
      this.extraSegsBed = xtraSegsBed;
      this.extraSegsFile = xtraSegsFile;
    }

  }

  private void runValidationGUI() {
    JFileChooser jfc = new JFileChooser(getLastUsedAnnotationDir());
    jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    jfc.setDialogTitle("Select BED File");
    int opt = jfc.showOpenDialog(frmAnnotator);
    if (opt == JFileChooser.APPROVE_OPTION) {
      String bedFile = jfc.getSelectedFile().getAbsolutePath();
      validationFiles.put(currSample, bedFile);
      runValidation(bedFile);
    }
  }

  private void runValidation(String bedFile) {
    ValidationResult result = validateAgainstFile(bedFile);
    String text = "";
    int diff1 = result.countBedSegs - result.countOverlap;
    int diff2 = result.countFiles - result.countOverlap;
    boolean valid = diff1 == 0 && diff2 == 0;
    if (diff1 > 0) {
      text += "+" + diff1;
    } else {
      text += diff1;
    }
    text += "/";
    if (diff2 > 0) {
      text += "+" + diff2;
    } else {
      text += diff2;
    }

    String mouseoverText = "<html>";
    mouseoverText += bedFile;
    mouseoverText += "<hr />";
    mouseoverText += "Found [" + result.extraSegsBed.size() + "] extra BED segments.<br />";
    mouseoverText += "Found [" + result.extraSegsFile.size() + "] extra image files.<br />";
    mouseoverText += "Found [" + result.countOverlap + "] overlapping segments/images.<br />";

    mouseoverText += "</html>";

    StringBuilder builder = new StringBuilder();
    for (Segment seg : result.extraSegsBed) {
      builder.append(seg.getChromosomeUCSC() + ":" + (seg.getStart() - 1) + "-" + seg.getStop())
             .append("\n");
    }
    String extraBedSegStr = builder.toString();

    builder = new StringBuilder();
    for (Segment seg : result.extraSegsFile) {
      builder.append(seg.getChromosomeUCSC() + ":" + (seg.getStart() - 1) + "-" + seg.getStop())
             .append("\n");
    }
    String extraFileSegStr = builder.toString();

    new Thread(new Runnable() {
      @Override
      public void run() {
        PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(bedFile, false)
                                                        + "_validation.out");
        writer.println("#Extra BED Segments");
        writer.print(extraBedSegStr);
        writer.println("#Extra IMAGE Segments");
        writer.print(extraFileSegStr);
        writer.close();
      }
    }).start();

    final JPopupMenu menu = new JPopupMenu();
    JMenuItem copyExtraBedSegs = new JMenuItem();
    copyExtraBedSegs.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent arg0) {
        ext.setClipboard(extraBedSegStr);
      }
    });
    copyExtraBedSegs.setText("Copy Extra BED Segments to Clipboard");
    menu.add(copyExtraBedSegs);
    JMenuItem copyExtraImages = new JMenuItem();
    copyExtraImages.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent arg0) {
        ext.setClipboard(extraFileSegStr);
      }
    });
    copyExtraImages.setText("Copy Extra Image Files to Clipboard");
    menu.add(copyExtraImages);

    validationCountLabel.setText(text);
    validationCountLabel.setForeground(valid ? Color.GREEN : Color.RED);
    validationCountLabel.setToolTipText(mouseoverText);
    validationCountLabel.addMouseListener(new MouseAdapter() {
      @Override
      public void mousePressed(MouseEvent e) {
        super.mousePressed(e);
        if (SwingUtilities.isRightMouseButton(e)) {
          menu.show(e.getComponent(), e.getX(), e.getY());
        }
      }

      @Override
      public void mouseReleased(MouseEvent e) {
        super.mouseReleased(e);
        if (SwingUtilities.isRightMouseButton(e) && menu.isVisible()) {
          menu.setVisible(false);
        }
      }
    });
    validationLabel.setVisible(true);
    validationCountLabel.setVisible(true);
    controlPanel.revalidate();
  }

  private ValidationResult validateAgainstFile(String bedFile) {
    Set<Segment> bedSegs = new HashSet<>(BEDFileReader.loadToSegments(bedFile));

    HashMap<String, HashMap<String, AnnotatedImage>> map = annotator.getAnnotationMap();
    String sampleName = (String) sampleCombo.getSelectedItem();
    HashMap<String, AnnotatedImage> ann = map.get(sampleName);
    List<String> files = ann == null ? new ArrayList<>() : new ArrayList<>(ann.keySet());

    Set<Segment> fileSegs = new HashSet<>();
    for (String file : files) {
      int[] v = parseSampleChrPosIfExists(file);
      Segment seg = new Segment((byte) v[0], v[1] + 1, v[2]);
      fileSegs.add(seg);
    }

    int bedSegCnt = bedSegs.size();
    int filesCnt = fileSegs.size();

    int overlapCnt = Sets.intersection(bedSegs, fileSegs).size();
    Set<Segment> extraBedSegs = new HashSet<>();
    Sets.difference(bedSegs, fileSegs).copyInto(extraBedSegs);
    Set<Segment> extraFileSegs = new HashSet<>();
    Sets.difference(fileSegs, bedSegs).copyInto(extraFileSegs);

    return new ValidationResult(bedSegCnt, filesCnt, overlapCnt, extraBedSegs, extraFileSegs);
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
      if (checkBxHideMissing && img.isMissing()) continue;
      boolean hidden = img.getAnnotations().stream().anyMatch(hiddenAnnotations::contains);
      if (hidden) continue;
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
      LoadedSaved saved = annotator.loadAnnotations(annFile);
      setLastSavedAnnotationFile(annFile);
      setLastUsedAnnotationDir(ext.verifyDirFormat(ext.parseDirectoryOfFile(annFile)));
      recentAnnotFiles.add(annFile);
      reloadControls();
      mnTravAnn.removeAll();
      for (JCheckBoxMenuItem mnItm : annHideMap.values()) {
        mnNavHide.remove(mnItm);
      }
      hiddenAnnotations.clear();
      for (Annotation a : annotator.getAnnotations()) {
        addAnnotationToTraversalMenu(a);
      }
      updateAvail();
      tree.setSelectionRow(0);
      tree.scrollRowToVisible(0);
      saveProperties();
      if (saved.missing > 0 || saved.unmatched.size() > 0) {
        JLabel label = new JLabel("<html>" + saved.getString() + "</html>");
        JScrollPane pane = new JScrollPane(label);
        JPanel panel = new JPanel(new BorderLayout());
        panel.add(pane, BorderLayout.CENTER);
        panel.add(new JLabel("<html>Image names must start with the name of the subdirectory they are in."
                             + (saved.missing > 0 ? "<br />Missing files can be shown/hidden in the menu."
                                                  : "")
                             + "</html>"),
                  BorderLayout.SOUTH);
        JOptionPane.showMessageDialog(ImageAnnotator.this.frmAnnotator, panel, "Loading Results",
                                      JOptionPane.INFORMATION_MESSAGE);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private volatile String prevSample;
  private volatile String currSample;

  private Map<String, String> validationFiles = new HashMap<>();

  ActionListener comboListener = new ActionListener() {

    @Override
    public void actionPerformed(ActionEvent e) {
      prevSample = currSample;
      currSample = (String) sampleCombo.getSelectedItem();
      boolean stt = startInteractionEvent();
      if (stt) {
        // either selected by-hand or by search
        if (prevSample == null) {
          buildingEvent.direction = DIRECTION.SKIP;
        } else {
          int prevInd = ((DefaultComboBoxModel<String>) sampleCombo.getModel()).getIndexOf(prevSample);
          int currInd = sampleCombo.getSelectedIndex();
          if (Math.abs(prevInd - currInd) > 1) {
            buildingEvent.direction = DIRECTION.SKIP;
          } else if (currInd > prevInd) {
            buildingEvent.direction = DIRECTION.NEXT_SMP;
          } else if (currInd < prevInd) {
            buildingEvent.direction = DIRECTION.PREV_SMP;
          } else {
            // no change - selected the same sample, which resets the selected image
            buildingEvent.direction = DIRECTION.PREV_IMG;
          }
        }
      }
      tree.setEnabled(true);
      if (prevSample != null && !prevSample.equals(currSample)) {
        validationLabel.setVisible(false);
        validationCountLabel.setVisible(false);
        validationCountLabel.setText("");
        validationCountLabel.setToolTipText("");
        for (MouseListener list : validationCountLabel.getMouseListeners()) {
          validationCountLabel.removeMouseListener(list);
        }
      }
      String file = validationFiles.get(currSample);
      if (file == null
          && Files.exists(lastOpenedImageDir + currSample + "/" + currSample + ".bed")) {
        file = lastOpenedImageDir + currSample + "/" + currSample + ".bed";
        validationFiles.put(currSample, file);
      }
      if (file != null && Files.exists(file)) {
        runValidation(file);
      }
      updateAvail();
      if (stt) {
        fireInteractionEvent();
      }
    }
  };

  TreeSelectionListener treeListener = new TreeSelectionListener() {

    @Override
    public void valueChanged(TreeSelectionEvent e) {
      boolean stt = startInteractionEvent(); // if selecting by hand

      if (e.getNewLeadSelectionPath() == null) {
        // selected new sample, set to first entry
        refreshAnnotations();
      } else {
        Object selected = e.getNewLeadSelectionPath().getLastPathComponent();
        if (!(selected instanceof DefaultMutableTreeNode)) {
          System.err.println("Error - selected an invalid item somehow.  Likely programmer error.");
          cancelInteractionEvent();
          return;
        }

        if (stt) {
          int prev = tree.getRowForPath(e.getOldLeadSelectionPath());
          int curr = tree.getRowForPath(e.getNewLeadSelectionPath());
          if (prev < curr) {
            buildingEvent.direction = DIRECTION.NEXT_IMG;
          } else if (prev > curr) {
            buildingEvent.direction = DIRECTION.PREV_IMG;
          }
        }

        DefaultMutableTreeNode selectedNode = (DefaultMutableTreeNode) selected;
        AnnotatedImage ann = ((AnnotatedImage) selectedNode.getUserObject());
        setSelectedNode(ann);
        refreshAnnotations();
      }
      if (stt) {
        fireInteractionEvent();
      }
    }
  };

  class InteractionEvent {

    long fireTime;
    DIRECTION direction;
    public TRAVERSAL traversal;
    public Annotation traversalAnn;
  }

  private DIRECTION getRecentUserDirection(int win) {
    int sz = interactLog.size();
    int window = (sz < win) ? sz : win;
    List<InteractionEvent> sub = interactLog.subList(sz - window, sz);

    int prev = 0;
    int next = 0;
    int skip = 0;
    for (InteractionEvent e : sub) {
      if (e.direction == DIRECTION.NEXT_IMG || e.direction == DIRECTION.NEXT_SMP) {
        next++;
      } else if (e.direction == DIRECTION.PREV_IMG || e.direction == DIRECTION.PREV_SMP) {
        prev++;
      } else if (e.direction == DIRECTION.SKIP) {
        skip++;
      }
    }
    if (skip > win / 2) {
      return DIRECTION.SKIP;
    } else if (prev > next) {
      return DIRECTION.PREV_IMG;
    } else
      return DIRECTION.NEXT_IMG;
  }

  private BEHAVIOR getRecentUserBehavior(int win) {
    int sz = interactLog.size();
    int window = (sz < win) ? sz : win;
    List<InteractionEvent> sub = interactLog.subList(sz - window, sz);

    MultiSet<BEHAVIOR> behav = new HashMultiSet<>();
    for (int i = sub.size() - 1; i > 0; i--) {
      DIRECTION p1 = sub.get(i).direction;
      DIRECTION p2 = sub.get(i - 1).direction;
      if (p1 == DIRECTION.SKIP || p2 == DIRECTION.SKIP) {
        behav.add(BEHAVIOR.SEARCH);
      } else if ((p1 == DIRECTION.NEXT_IMG || p1 == DIRECTION.NEXT_SMP)
                 && (p2 == DIRECTION.NEXT_IMG || p2 == DIRECTION.NEXT_SMP)) {
        behav.add(BEHAVIOR.LINEAR);
      } else if ((p1 == DIRECTION.PREV_IMG || p1 == DIRECTION.PREV_SMP)
                 && (p2 == DIRECTION.PREV_IMG || p2 == DIRECTION.PREV_SMP)) {
        behav.add(BEHAVIOR.LINEAR);
      } else if ((p1 == DIRECTION.PREV_IMG || p1 == DIRECTION.PREV_SMP)
                 && (p2 == DIRECTION.NEXT_IMG || p2 == DIRECTION.NEXT_SMP)) {
        behav.add(BEHAVIOR.REVIEW);
      } else if ((p1 == DIRECTION.NEXT_IMG || p1 == DIRECTION.NEXT_SMP)
                 && (p2 == DIRECTION.PREV_IMG || p2 == DIRECTION.PREV_SMP)) {
        behav.add(BEHAVIOR.REVIEW);
      }
    }

    int lin = behav.getCount(BEHAVIOR.LINEAR);
    int rev = behav.getCount(BEHAVIOR.REVIEW);
    int sch = behav.getCount(BEHAVIOR.SEARCH);

    if (sch > lin && sch > rev) return BEHAVIOR.SEARCH;
    if (lin > rev) {
      if ((lin - rev) > rev) {
        // twice as many linear movements as reviews
        return BEHAVIOR.LINEAR;
      } else {
        // more linear, but not a lot more
        return BEHAVIOR.UNSURE;
      }
    } else if (rev > lin) {
      if (rev - lin > lin) {
        return BEHAVIOR.REVIEW;
      } else {
        return BEHAVIOR.UNSURE;
      }
    } else {
      return BEHAVIOR.UNSURE;
    }
  }

  private SPEED getRecentUserSpeed(int win) {
    int sz = interactLog.size();
    if (sz == 0) return SPEED.SLOW;
    int window = (sz < win) ? sz : win;
    List<InteractionEvent> sub = interactLog.subList(sz - window, sz);
    double currDelta = System.nanoTime() - sub.get(sub.size() - 1).fireTime;
    if (TimeUnit.NANOSECONDS.toSeconds((long) currDelta) > 15) {
      // if user has stopped interacting with program
      return SPEED.SLOW;
    }
    double[] deltas = new double[sub.size() - 1];
    for (int i = sub.size() - 1; i > 0; i--) {
      deltas[sub.size() - (i + 1)] = sub.get(i).fireTime - sub.get(i - 1).fireTime;
    }

    MultiSet<SPEED> set = new HashMultiSet<>();
    for (int i = 0; i < deltas.length; i++) {
      long secondsDelta = TimeUnit.NANOSECONDS.toSeconds((long) deltas[i]);
      long minutesDelta = TimeUnit.NANOSECONDS.toMinutes((long) deltas[i]);
      long millisDelta = TimeUnit.NANOSECONDS.toMillis((long) deltas[i]);
      if (minutesDelta > 0 || secondsDelta > 2 || millisDelta > 1500) {
        set.add(SPEED.SLOW);
      } else if (millisDelta > 600) {
        set.add(SPEED.MODERATE);
      } else if (millisDelta < 200) {
        set.add(SPEED.REALLY_FAST);
      } else {
        set.add(SPEED.FAST);
      }
    }

    int slow = set.getCount(SPEED.SLOW);
    int mod = set.getCount(SPEED.MODERATE);
    int fast = set.getCount(SPEED.FAST);
    int really = set.getCount(SPEED.FAST);

    int comb = (int) ((slow * -1.5) + fast + really * .5); // weight slow heavier to decay more
    // quickly
    if (mod > Math.abs(comb) || Math.abs(comb) < (win / 2)) {
      return SPEED.MODERATE;
    } else if (comb < 0) {
      return SPEED.SLOW;
    } else {
      if (really > fast) return SPEED.REALLY_FAST;
      return SPEED.FAST;
    }
  }

  volatile boolean process = true;
  Thread interactionProcessingThread = new Thread(new Runnable() {

    @Override
    public void run() {
      while (process) {
        InteractionEvent e;
        try {
          e = interactionQueue.take();
        } catch (InterruptedException e1) {
          System.err.println("Error - interruption in UX processing thread.  Interaction may become stilted.");
          return;
        }
        interactLog.add(e);
      }
    }
  });

  List<InteractionEvent> interactLog = new ArrayList<>();
  LinkedBlockingQueue<InteractionEvent> interactionQueue = new LinkedBlockingQueue<>();

  volatile InteractionEvent buildingEvent = null;

  private boolean startInteractionEvent() {
    if (buildingEvent != null) return false;
    buildingEvent = new InteractionEvent();
    buildingEvent.fireTime = System.nanoTime();
    buildingEvent.traversal = getTraversal();
    if (buildingEvent.traversal == TRAVERSAL.CUSTOM) {
      buildingEvent.traversalAnn = getTravAnnotation();
    }
    return true;
  }

  private void fireInteractionEvent() {
    interactionQueue.offer(buildingEvent);
    buildingEvent = null;
  }

  private void cancelInteractionEvent() {
    buildingEvent = null;
  }

  private final Action loadAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      loadImageFiles();
    }
  };

  private final Action existAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      loadAnnotations();
    }
  };

  private final Action exportAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      new AnnotationExportDialog(annotator).setVisible(true);
      // export specific images with specific annotations
      // export sample names with specific annotations
    }
  };

  private final Action validateAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent arg0) {
      runValidationGUI();
    }
  };

  private final Action saveAsAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      saveAnnotations(true);
    }
  };
  private final Action saveAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      saveAnnotations(false);
    }
  };
}
