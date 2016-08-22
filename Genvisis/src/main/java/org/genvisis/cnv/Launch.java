package org.genvisis.cnv;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import java.util.jar.Attributes;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.DeNovoCNV;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsCrossTabs;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsManhattan;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.ProjectPropertiesEditor;
import org.genvisis.cnv.gui.FileAndOutputSelectorGUI;
import org.genvisis.cnv.gui.ImportProjectGUI;
import org.genvisis.cnv.gui.PlinkExportOptions;
import org.genvisis.cnv.manage.DemoPackage;
import org.genvisis.cnv.manage.ExportCNVsToPedFormat;
import org.genvisis.cnv.manage.GenvisisWorkflow;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.plots.CompPlot;
import org.genvisis.cnv.plots.ForestPlotFrame;
import org.genvisis.cnv.plots.LinePlot;
import org.genvisis.cnv.plots.MosaicPlot;
import org.genvisis.cnv.plots.QQPlotFrame;
import org.genvisis.cnv.plots.ScatterPlot;
import org.genvisis.cnv.plots.SexPlot;
import org.genvisis.cnv.plots.StratPlot;
import org.genvisis.cnv.plots.Trailer;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.common.Aliases;
import org.genvisis.common.CmdLine;
import org.genvisis.common.CurrentManifest;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.HttpUpdate;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.cyto.CytoGUI;

// -XX:+UseConcMarkSweepGC
// -XX:+UseParNewGC

public class Launch extends JFrame implements ActionListener, WindowListener, ItemListener {
  public static final long serialVersionUID = 1L;

  public static final String VERSION = "0.60";

  public static final String EXIT = "Exit";
  public static final String EDIT = "Project Properties Editor";
  public static final String REFRESH = "Refresh";
  public static final String PIPELINE = "Genvisis Project Workflow";
  public static final String NEW_PROJECT = "New Project";
  public static final String IMPORT_PROJECT = "Import Project";

  public static final String MAP_FILES = "Map .csv files to IDs";
  public static final String GENERATE_MARKER_POSITIONS = "Generate marker positions file";
  public static final String PARSE_FILES_CSV = "Parse .csv files";
  public static final String TRANSPOSE_DATA = "Transpose data";
  public static final String MITOPIPELINE = "MitoPipeline";

  public static final String CHECK_SEX = "Check sex";
  public static final String LRR_SD = "LRR Stdevs";
  public static final String CNP_SCAN = "Scan for CNPs";
  public static final String MOSAICISM = "Determine mosaic arms";
  public static final String MARKER_METRICS = "Full QC marker metrics";
  public static final String FILTER_MARKER_METRICS = "Filter marker metrics";
  public static final String TALLY_MARKER_ANNOTATIONS = "Tally marker annotations";
  public static final String TALLY_WITHOUT_DETERMINING_DROPS =
                                                             "Tally without determining dropped markers (much faster)";
  public static final String TALLY_CLUSTER_FILTERS = "Tally all reclustered markers";

  public static final String SCATTER = "Scatter Plot";
  public static final String QQ = "QQ Plot";
  public static final String STRAT = "Stratify Plot";
  public static final String MOSAIC_PLOT = "Mosaic Pot";
  public static final String SEX_PLOT = "Sex Plot";
  public static final String TRAILER = "Trailer Plot";
  public static final String TWOD = "2D Plot";
  public static final String LINE_PLOT = "Line Plot";
  public static final String COMP = "Comp Plot";
  public static final String FOREST_PLOT = "Forest Plot";

  public static final String GENERATE_ABLOOKUP = "Generate AB Lookup";
  public static final String EXPORT_TO_PLINK = "Export to PLINK format";
  public static final String GENERATE_PENNCNV_FILES = "Generate PennCNV files";
  public static final String PARSE_RAW_PENNCNV_RESULTS = "Parse raw PennCNV results files";
  public static final String POPULATIONBAF = "Compute Population BAF file";
  public static final String GCMODEL = "Compute GC model file";
  public static final String CUSTOM_CENTROIDS = "Compute custom centroids file";

  public static final String DENOVO_CNV = "De Novo CNV";
  public static final String EXPORT_CNVS = "Export CNVs to Pedfile format";
  public static final String CYTO_WORKBENCH = "Parse workbench files";
  public static final String PRINCIPAL_COMPONENTS = "Principal Components";
  public static final String GENERATE_DEMO_PACKAGE = "Generate a demo package";
  public static final String ADD_QC_TO_SAMPLE_DATA = "Add sample qc metrics to sample data";
  public static final String CHECK_FOR_UPDATES = "Check for updates";

  public static final String TEST = "Test new program";

  public static final String[][] MENUS =
                                 {{"File", NEW_PROJECT, IMPORT_PROJECT, "Select Project", EDIT,
                                   "Preferences", CHECK_FOR_UPDATES, EXIT},
                                  {"Data", MAP_FILES, GENERATE_MARKER_POSITIONS, PARSE_FILES_CSV,
                                   TRANSPOSE_DATA, PIPELINE}, // ,
                                                              // MITOPIPELINE
                                  {"Quality", CHECK_SEX, LRR_SD, CNP_SCAN, MOSAICISM,
                                   MARKER_METRICS, FILTER_MARKER_METRICS, TALLY_MARKER_ANNOTATIONS,
                                   TALLY_WITHOUT_DETERMINING_DROPS, TALLY_CLUSTER_FILTERS},
                                  {"Plots", SCATTER, QQ, STRAT, MOSAIC_PLOT, SEX_PLOT, TRAILER,
                                   TWOD, LINE_PLOT, COMP, FOREST_PLOT},
                                  {"Tools", GENERATE_ABLOOKUP, EXPORT_TO_PLINK,
                                   GENERATE_PENNCNV_FILES, PARSE_RAW_PENNCNV_RESULTS, POPULATIONBAF,
                                   GCMODEL, CUSTOM_CENTROIDS, DENOVO_CNV, EXPORT_CNVS,
                                   CYTO_WORKBENCH, PRINCIPAL_COMPONENTS, GENERATE_DEMO_PACKAGE,
                                   ADD_QC_TO_SAMPLE_DATA, TEST},
                                  {"Help", "Contents", "Search", "About"}};


  private Project proj;
  private final boolean jar;
  private JComboBox<String> projectsBox;
  private String[] projects;
  private LaunchProperties launchProperties;
  private final String launchPropertiesFile;
  private JTextArea output;
  private JScrollPane scrollPane;
  private final Vector<Thread> threadsRunning;
  private int indexOfCurrentProj;
  private long timestampOfPropertiesFile;
  private long timestampOfSampleDataFile;
  private Logger log;

  private JProgressBar progBar;

  public Launch(String launchPropertiesFile, CurrentManifest currentManifest, boolean jar) {

    super("Genvisis " + currentManifest.getVersion().getVersion());
    this.jar = jar;
    this.launchPropertiesFile = launchPropertiesFile;
    timestampOfPropertiesFile = -1;
    timestampOfSampleDataFile = -1;
    threadsRunning = new Vector<Thread>();
  }

  public void loadProjects() {
    String[] projectNames;

    projects = Files.list(launchProperties.getDirectory(), ".properties", false);
    projectNames = launchProperties.getListOfProjectNames();
    projectsBox.setModel(new DefaultComboBoxModel<String>(projectNames));
  }


  public Project loadProject() {
    proj = new Project(launchProperties.getDirectory() + projects[indexOfCurrentProj], jar);
    proj.setGuiState(true);
    timestampOfPropertiesFile = new Date().getTime();
    timestampOfSampleDataFile = new Date().getTime();
    if (!Files.exists(proj.PROJECT_DIRECTORY.getValue(), proj.JAR_STATUS.getValue())) {
      JOptionPane.showMessageDialog(null,
                                    "Error - the directory ('" + proj.PROJECT_DIRECTORY.getValue()
                                          + "') for project '" + proj.PROJECT_NAME.getValue()
                                          + "' did not exist; creating now. If this was in error, please edit the property file.",
                                    "Error", JOptionPane.ERROR_MESSAGE);
    }

    log = proj.getLog();
    log.linkTextArea(output);

    progBar.setIndeterminate(false);
    progBar.setValue(0);
    progBar.setMaximum(0);
    progBar.setMinimum(0);
    progBar.setString(null);
    progBar.setStringPainted(false);

    proj.initializeProgressMonitor(progBar);
    projectsBox.setSelectedIndex(indexOfCurrentProj);

    return proj;
  }

  public void setIndexOfCurrentProject(String projPropertiesFileName) {
    indexOfCurrentProj = 0;
    for (int i = 0; i < projects.length; i++) {
      if (projects[i].equals(projPropertiesFileName)) {
        indexOfCurrentProj = i;
      }
    }
    if (projects.length > 0) {
      projectsBox.setSelectedIndex(indexOfCurrentProj);
    }
  }

  @SuppressWarnings("rawtypes")
  public static void setUIFont(Font newFont) {
    Enumeration keys = UIManager.getDefaults().keys();
    while (keys.hasMoreElements()) {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if (value != null && value instanceof javax.swing.plaf.FontUIResource) {
        Font oldFont = UIManager.getFont(key);
        UIManager.put(key, newFont.deriveFont(oldFont.getStyle(), oldFont.getSize2D()));
      }
    }
  }

  private static final class ExceptionHandler implements Thread.UncaughtExceptionHandler {
    public Logger log;

    public void setLog(Logger log) {
      this.log = log;
    }

    @Override
    public void uncaughtException(Thread t, Throwable e) {
      if (log != null) {
        log.reportError("Uncaught Exception in Thread {" + t.getName() + "}:");
        log.reportException(e);
      } else {
        System.err.println("Error - Uncaught Exception in Thread {" + t.getName() + "}:");
        e.printStackTrace();
      }
    }
  }

  private static void createAndShowGUI() {
    String launchPropertiesFile;
    final Launch frame;

    try {
      UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
    } catch (Exception e2) {
      System.err.println("Failed loading CrossPlatformLookAndFeel");
      System.err.println(e2);
    }

    ExceptionHandler ueh = new ExceptionHandler();

    Thread.setDefaultUncaughtExceptionHandler(ueh);

    // set system-wide anti-aliasing
    System.setProperty("awt.useSystemAAFontSettings", "on");
    System.setProperty("swing.aatext", "true");

    ToolTipManager.sharedInstance().setInitialDelay(0);
    ToolTipManager.sharedInstance().setDismissDelay(Integer.MAX_VALUE - 1);
    ToolTipManager.sharedInstance().setReshowDelay(0);

    UIManager.put("ToolTip.background", Color.decode("#F5F5DC"));

    // Create and set up the content pane.
    launchPropertiesFile = LaunchProperties.DEFAULT_PROPERTIES_FILE;
    initLaunchProperties(launchPropertiesFile, false, true);
    CurrentManifest manifest = new CurrentManifest(new Attributes());
    try {// try not to break the launch so we will catch anything
      manifest = CurrentManifest.loadGenvisisManifest();
    } catch (Exception e) {
    }
    frame = new Launch(launchPropertiesFile, manifest, false);
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.setJMenuBar(frame.topMenuBar());
    frame.setContentPane(frame.createContentPane());
    frame.createPopupMenu();

    frame.setSize(650, 550);
    frame.setLocation(300, 200);

    frame.addWindowListener(frame);

    javax.swing.SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        frame.setVisible(true);
      }
    });

    // TODO only instantiate when used
    frame.setIndexOfCurrentProject(frame.launchProperties.getProperty(LaunchProperties.LAST_PROJECT_OPENED));
    if (frame.projects.length > 0) {
      frame.loadProject();
    }

    ueh.setLog(frame.log);
    try {
      frame.log.report(HttpUpdate.checkGenvisisVersion(frame.log));
    } catch (Exception e) {

    }
  }

  public static void initLaunchProperties(String launchPropertiesFile, boolean force,
                                          boolean relativePath) {
    String path = LaunchProperties.directoryOfLaunchProperties(launchPropertiesFile);
    String pathToSet;
    if (relativePath) {
      pathToSet = "";
    } else {
      pathToSet = path;
    }
    if (force || !new File(launchPropertiesFile).exists()) {
      System.out.println("creating " + launchPropertiesFile);
      new File(path + "projects/").mkdirs();
      new File(path + "example/").mkdirs();
      Files.writeList(new String[] {"LAST_PROJECT_OPENED=example.properties",
                                    "PROJECTS_DIR=" + pathToSet + "projects/"},
                      launchPropertiesFile);
      if (!new File(path + "projects/example.properties").exists()) {
        Files.writeList(new String[] {"PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/",
                                      "SOURCE_DIRECTORY=sourceFiles/"},
                        path + "projects/example.properties");
      }
    }
    String bat = path + "Launch.bat";
    String sh = path + "Launch.sh";
    String command = path + "Launch.command";

    if (!Files.exists(bat)) {
      Files.write(getLaunchBat(), bat);
    }
    if (!Files.exists(sh)) {
      Files.write(getLaunchSH(), sh);
      Files.chmod(sh);
    }
    if (!Files.exists(command)) {
      Files.write(getLaunchSH(), command);
      Files.chmod(command);
    }

  }

  public Container createContentPane() {
    // Create the content-pane-to-be.
    JPanel contentPane = new JPanel(new BorderLayout());
    contentPane.setOpaque(true);

    contentPane.add(topIconBar(), BorderLayout.NORTH);

    // Create a scrolled text area.
    output = new JTextArea(5, 30);
    output.setEditable(false);
    scrollPane = new JScrollPane(output);

    // Add the text area to the content pane.
    contentPane.add(scrollPane, BorderLayout.CENTER);

    progBar = new JProgressBar();
    contentPane.add(progBar, BorderLayout.SOUTH);
    progBar.setVisible(false);

    return contentPane;
  }

  private JMenuBar topMenuBar() {
    JMenuBar menuBar;
    JMenu menu, submenu;
    JMenuItem menuItem;
    Hashtable<Character, String> hash;

    launchProperties = new LaunchProperties(launchPropertiesFile);
    menuBar = new JMenuBar();
    for (int i = 0; i < MENUS.length; i++) {
      menu = new JMenu(MENUS[i][0]);
      menu.setMnemonic((int) MENUS[i][0].charAt(0));
      menuBar.add(menu);
      hash = new Hashtable<Character, String>();
      for (int j = 1; j < MENUS[i].length; j++) {
        if (MENUS[i][j] == "") {
          break;
        }
        if (MENUS[i][j] == "1") {
          menu.addSeparator();
        } else if (MENUS[i][j].equals(NEW_PROJECT)) {
          menuItem = new JMenuItem(NEW_PROJECT);
          menuItem.addActionListener(this);
          menuItem.setMnemonic(KeyEvent.VK_N);
          menu.add(menuItem);
        } else if (MENUS[i][j].equals(IMPORT_PROJECT)) {
          menuItem = new JMenuItem(IMPORT_PROJECT);
          menuItem.addActionListener(this);
          menuItem.setMnemonic(KeyEvent.VK_I);
          menu.add(menuItem);
        } else if (MENUS[i][j].equals("Select Project")) {
          submenu = new JMenu(MENUS[i][j]);
          projects = Files.list(launchProperties.getDirectory(), ".properties", false);
          for (String project : projects) {
            menuItem = new JMenuItem(ext.rootOf(project, true) + " ");
            menuItem.addActionListener(this);
            submenu.add(menuItem);
          }
          menu.add(submenu);
        } else if (MENUS[i][j].equals(PRINCIPAL_COMPONENTS)) {
          String[] pcSubMenuOptions =
                                    new String[] {PrincipalComponentsManhattan.PRINCIPAL_MANHATTAN_MI,
                                                  PrincipalComponentsCrossTabs.PRINCIPAL_CROSSTABS_MI};
          JMenu pcSubMenu = new JMenu(MENUS[i][j]);
          for (String pcSubMenuOption : pcSubMenuOptions) {
            JMenuItem pcSubItem = new JMenuItem(pcSubMenuOption);
            pcSubItem.addActionListener(this);
            pcSubMenu.add(pcSubItem);
          }
          menu.add(pcSubMenu);

        } else {
          menuItem = new JMenuItem(MENUS[i][j]);
          for (int k = 0; k < MENUS[i][j].length(); k++) {
            if (!hash.containsKey(MENUS[i][j].toLowerCase().charAt(k))) {
              menuItem.setMnemonic((int) MENUS[i][j].toLowerCase().charAt(k));
              hash.put(MENUS[i][j].toLowerCase().charAt(k), "");
              k = MENUS[i][j].length();
            }
          }
          menuItem.addActionListener(this);
          // TODO What is this?
          menuItem.getAccessibleContext().setAccessibleDescription("Under development");
          menu.add(menuItem);
        }
      }
    }
    return menuBar;
  }

  private JPanel topIconBar() {
    JPanel iconBar;
    JButton button;
    String[] icons = null;
    String[] commands = null;

    icons = new String[] {"images/save1.png", "images/edit1.png", "images/refresh.gif",
                          "images/gen_pipe_1.png", "images/scatterPlot2.png",
                          "images/trailerPlot2.png", "images/qqplot.gif", "images/recluster1.png",
                          "images/twoDPlot1.jpg", "images/forestPlot1.png"};
    commands = new String[] {"", EDIT, REFRESH, PIPELINE, SCATTER, TRAILER, QQ, LINE_PLOT, TWOD,
                             FOREST_PLOT};


    iconBar = new JPanel();
    iconBar.setLayout(new FlowLayout(FlowLayout.LEFT));
    for (int i = 0; i < icons.length; i++) {
      button = new JButton(Grafik.getImageIcon(icons[i]));
      button.setActionCommand(commands[i]);
      button.addActionListener(this);
      button.setToolTipText(commands[i]);
      button.setPreferredSize(new Dimension(25, 25));
      button.setBorder(null);
      iconBar.add(button);
    }
    addComponentsToPane(iconBar);

    return iconBar;
  }

  public void addComponentsToPane(final Container pane) {
    launchProperties = new LaunchProperties(launchPropertiesFile);
    projectsBox = new JComboBox();
    loadProjects();
    // In JDK1.4 this prevents action events from being fired when the up/down arrow keys are used
    // on the dropdown menu
    projectsBox.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);

    setIndexOfCurrentProject(launchProperties.getProperty(LaunchProperties.LAST_PROJECT_OPENED));
    if (indexOfCurrentProj > 0 && projectsBox.getItemCount() > 0) {
      projectsBox.setSelectedIndex(indexOfCurrentProj);
    }
    projectsBox.addItemListener(this);
    pane.add(projectsBox);
  }

  public void createPopupMenu() {
    JMenuItem menuItem;

    // Create the popup menu.
    JPopupMenu popup = new JPopupMenu();
    menuItem = new JMenuItem("A popup menu item");
    menuItem.addActionListener(this);
    popup.add(menuItem);
    menuItem = new JMenuItem("Another popup menu item");
    menuItem.addActionListener(this);
    popup.add(menuItem);
  }

  public class IndependentThread implements Runnable {
    private final Project proj;
    private final String command;

    public IndependentThread(Project proj, String command) {
      this.proj = proj;
      this.command = command;
    }

    @Override
    public void run() {
      /*
       * CAUTION/NOTE/TODO: ALL SWING CALLS OR COMPONENT CREATION SHOULD BE WRAPPED IN
       * SwingUtilities.invokeLater();
       */
      if (command.equals(MAP_FILES)) {
        org.genvisis.cnv.manage.SourceFileParser.mapFilenamesToSamples(proj,
                                                                       "filenamesMappedToSamples.txt");
      } else if (command.equals(GENERATE_MARKER_POSITIONS)) {
        org.genvisis.cnv.manage.Markers.generateMarkerPositions(proj,
                                                                proj.getLocationOfSNP_Map(true));
      } else if (command.equals(PARSE_FILES_CSV)) {
        org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, proj.NUM_THREADS.getValue());
      } else if (command.equals(CHECK_SEX)) {
        String blastAnnotationFile = proj.BLAST_ANNOTATION_FILENAME.getValue();
        String nonCrossHybridizingMarkersFile =
                                              MarkerBlastQC.defaultOneHitWondersFilename(blastAnnotationFile);
        if (!Files.exists(nonCrossHybridizingMarkersFile)) {
          if (Files.exists(blastAnnotationFile)) {
            MarkerBlastQC.getOneHitWonders(proj, blastAnnotationFile,
                                           nonCrossHybridizingMarkersFile, 0.8, proj.getLog());
          } else {
            nonCrossHybridizingMarkersFile = null;
          }
        }
        org.genvisis.cnv.qc.SexChecks.sexCheck(proj, true, nonCrossHybridizingMarkersFile);
      } else if (command.equals(TRANSPOSE_DATA)) {
        TransposeData.transposeData(proj, 2000000000, false);
      } else if (command.equals(GENERATE_ABLOOKUP)) {
        ABLookup abLookup;
        String filename;

        filename = proj.PROJECT_DIRECTORY.getValue()
                   + ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
        if (!Files.exists(filename)) {
          abLookup = new ABLookup();
          abLookup.parseFromAnnotationVCF(proj);
          abLookup.writeToFile(filename, proj.getLog());
        }

        ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false);
      } else if (command.equals(EXPORT_TO_PLINK)) {
        PlinkExportOptions peo = new PlinkExportOptions(proj);
        peo.setModal(true);
        peo.setVisible(true);
        if (peo.getCancelled()) {
          return;
        }
        String plinkFileroot = peo.getPlinkRoot();
        if (plinkFileroot == null) {
          return;
        }
        String pedFile = peo.getPedigree();// will change the value of PEDIGREE_FILENAME, so the
                                           // return value here isn't necessary
        if (!new File(pedFile).exists()) {
          log.reportFileNotFound(pedFile);
          return;
        }
        String clusterFiltersFilename = peo.getClusterFilterSelection();
        if (clusterFiltersFilename != null) {
          clusterFiltersFilename = proj.DATA_DIRECTORY.getValue() + clusterFiltersFilename;
          // only care about ab lookup if cluster filters are applied
          String abFile = peo.getABFilename();
          if (abFile == null) {
            ABLookup abLookup;
            String filename;
            filename = proj.PROJECT_DIRECTORY.getValue()
                       + ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
            if (!Files.exists(filename)) {
              abLookup = new ABLookup();
              abLookup.parseFromOriginalGenotypes(proj);
              abLookup.writeToFile(filename, proj.getLog());
            }
            ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false);
            proj.AB_LOOKUP_FILENAME.setValue(filename);
          } else if (!new File(abFile).exists()) {
            log.reportFileNotFound(abFile);
            return;
          } else {
            proj.AB_LOOKUP_FILENAME.setValue(abFile);
          }
        }
        String targetMarkersFilename = peo.getTargetMarkersFile();
        if (peo.getCancelled()) { // getTargetMarkersFile(), if set to CREATE_NEW, can potentially
                                  // be cancelled
          return;
        }

        proj.saveProperties();
        boolean success = false;
        if (peo.exportAsBinary()) {
          success = PlinkData.saveGenvisisToPlinkBedSet(proj, plinkFileroot, clusterFiltersFilename,
                                                        targetMarkersFilename, -1, true);
        } else {
          success = PlinkData.saveGenvisisToPlinkPedSet(proj, plinkFileroot, clusterFiltersFilename,
                                                        targetMarkersFilename);
        }
        if (success) {
          log.report("Success!");
        }
      } else if (command.equals(GENERATE_PENNCNV_FILES)) {
        org.genvisis.cnv.analysis.AnalysisFormats.penncnv(proj, proj.getSampleList().getSamples(),
                                                          null, null, proj.NUM_THREADS.getValue());
      } else if (command.equals(PARSE_RAW_PENNCNV_RESULTS)) {
        // TODO make dialog to ask for filenames with a JCheckBox for denovo parsing
        org.genvisis.cnv.analysis.PennCNV.parseWarnings(proj, "penncnv.log");
        org.genvisis.cnv.analysis.PennCNV.parseResults(proj, "penncnv.rawcnv", false);
      } else if (command.equals(LRR_SD)) {
        org.genvisis.cnv.qc.LrrSd.init(proj, null, null, proj.getProperty(proj.NUM_THREADS));
      } else if (command.equals(CNP_SCAN)) {
        // TODO Genotyping
        // new ScanForCnp(proj, "CNPScanResult.txt");
      } else if (command.equals(DENOVO_CNV)) {
        DeNovoCNV.main("");
      } else if (command.equals(SCATTER)) {
        ScatterPlot.createAndShowGUI(proj, null, null, false);
      } else if (command.equals(QQ)) {
        QQPlotFrame.loadPvals(proj.QQ_FILENAMES.getValue(), "Q-Q Plot",
                              proj.getProperty(proj.DISPLAY_QUANTILES),
                              proj.getProperty(proj.DISPLAY_STANDARD_QQ),
                              proj.getProperty(proj.DISPLAY_ROTATED_QQ), -1, false,
                              proj.QQ_MAX_NEG_LOG10_PVALUE.getValue(), proj.getLog());
      } else if (command.equals(STRAT)) {
        StratPlot.loadStratificationResults(proj);
      } else if (command.equals(MOSAICISM)) {
        Mosaicism.findOutliers(proj);
      } else if (command.equals(MOSAIC_PLOT)) {
        MosaicPlot.loadMosaicismResults(proj);
      } else if (command.equals(SEX_PLOT)) {
        SexPlot.loadSexCheckResults(proj);
      } else if (command.equals(TRAILER)) {
        new Trailer(proj, null, proj.CNV_FILENAMES.getValue(), Trailer.DEFAULT_LOCATION);
      } else if (command.equals(TWOD)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            TwoDPlot.createGUI(proj, true, true);
          }
        });
      } else if (command.equals(LINE_PLOT)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            LinePlot.createAndShowGUI(proj);
          }
        });
      } else if (command.equals(COMP)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            new CompPlot(proj);
          }
        });
      } else if (command.equals(FOREST_PLOT)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            new ForestPlotFrame(proj);
          }
        });
      } else if (command.equals(POPULATIONBAF)) {
        org.genvisis.cnv.analysis.PennCNV.populationBAF(proj);
      } else if (command.equals(CUSTOM_CENTROIDS)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            CentroidCompute.computeAndDumpCentroids(proj);
          }
        });
      } else if (command.equals(EXPORT_CNVS)) {

        String[] inOut =
                       FileAndOutputSelectorGUI.showFileAndOutputSelector(Launch.this, null,
                                                                          JFileChooser.FILES_ONLY,
                                                                          null, null,
                                                                          JFileChooser.FILES_ONLY);
        if (inOut == null) {
          return;
        }

        String cnvFilename = inOut[0];
        String pedFilename = proj.PEDIGREE_FILENAME.getValue();
        String outputRoot = ext.rootOf(inOut[1], false);
        String endOfLine = Files.isWindows() ? "\r\n" : "\n";
        String fileFormat = ExportCNVsToPedFormat.PLINK_BINARY_FORMAT;
        boolean includeDele = true;
        boolean includeDupl = true;
        boolean ordered = false;
        boolean collapsed = false;
        boolean homozygous = false;
        boolean excludeMonomorphicLoci = false;
        int lociPerFile = Integer.MAX_VALUE;
        int window = 0;

        ExportCNVsToPedFormat.export(cnvFilename, pedFilename, outputRoot, endOfLine, fileFormat,
                                     includeDele, includeDupl, ordered, collapsed, homozygous,
                                     excludeMonomorphicLoci, lociPerFile, window, proj.getLog());

      } else if (command.equals(CYTO_WORKBENCH)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            new CytoGUI(proj, proj.PROJECT_DIRECTORY.getValue(), null);
          }
        });
      } else if (command.equals(TEST)) {
        org.genvisis.cnv.qc.SexChecks.sexCheck(proj, true);
        org.genvisis.cnv.qc.LrrSd.init(proj, null, null, proj.getProperty(proj.NUM_THREADS));
        Mosaicism.findOutliers(proj);

        PlinkData.saveGenvisisToPlinkPedSet(proj, "gwas", null,
                                            proj.TARGET_MARKERS_FILENAMES.getValue()[0]);
        CmdLine.run("plink --file gwas --make-bed --out plink", proj.PROJECT_DIRECTORY.getValue());
        new File(proj.PROJECT_DIRECTORY.getValue() + "genome/").mkdirs();
        CmdLine.run("plink --bfile ../plink --freq", proj.PROJECT_DIRECTORY.getValue() + "genome/");
        CmdLine.run("plink --bfile ../plink --missing",
                    proj.PROJECT_DIRECTORY.getValue() + "genome/");


      } else if (command.equals(GCMODEL)) {
        org.genvisis.cnv.analysis.PennCNV.gcModel(proj,
                                                  Files.firstPathToFileThatExists(Aliases.REFERENCE_FOLDERS,
                                                                                  "gc5Base.txt",
                                                                                  true, false, log),
                                                  proj.PROJECT_DIRECTORY.getValue() + "data/custom.gcModel",
                                                  100);
      } else if (command.equals(MARKER_METRICS)) {
        org.genvisis.cnv.qc.MarkerMetrics.fullQC(proj, proj.getSamplesToExclude(), null, true,
                                                 proj.NUM_THREADS.getValue());
      } else if (command.equals(FILTER_MARKER_METRICS)) {
        org.genvisis.cnv.qc.MarkerMetrics.filterMetrics(proj);
      } else if (command.equals(TALLY_MARKER_ANNOTATIONS)) {
        MarkerMetrics.tallyFlaggedReviewedChangedAndDropped(proj, true);
      } else if (command.equals(TALLY_WITHOUT_DETERMINING_DROPS)) {
        MarkerMetrics.tallyFlaggedReviewedChangedAndDropped(proj, false);
      } else if (command.equals(TALLY_CLUSTER_FILTERS)) {
        MarkerMetrics.tallyClusterFilters(proj, proj.getSamplesToInclude(null), null);
      } else if (command.equals(MITOPIPELINE)) {
      } else if (command.equals(PIPELINE)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            GenvisisWorkflow kAndK = new GenvisisWorkflow(proj, Launch.this);
            kAndK.showDialogAndRun();
          }
        });
      } else if (command.equals(PrincipalComponentsManhattan.PRINCIPAL_MANHATTAN_MI)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            PrincipalComponentsManhattan.guiAccess(proj, null);
          }
        });
      } else if (command.equals(PrincipalComponentsCrossTabs.PRINCIPAL_CROSSTABS_MI)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            PrincipalComponentsCrossTabs.guiAccess(proj, null);
          }
        });
      } else if (command.equals(GENERATE_DEMO_PACKAGE)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            DemoPackage.guiAccess(proj);
          }
        });

      } else if (command.equals(ADD_QC_TO_SAMPLE_DATA)) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            SampleQC.parseAndAddToSampleData(proj, 10, -1, false, false, null, false);
          }
        });

      } else {
        log.reportError("Error - unknown command: " + command);
      }
    }
  }

  @Override
  public void actionPerformed(ActionEvent ae) {
    String command = ae.getActionCommand();
    Thread thread;

    if (log == null) {
      log = new Logger();
    }
    log.report("Action performed: " + command + "\n");

    if (proj == null) {
      log.report("Trying again to load project");
      loadProject();
    } else if (timestampOfPropertiesFile < new File(proj.getPropertyFilename()).lastModified()) {
      log.report("Detected a change in the project properties file; reloading from '"
                 + proj.getPropertyFilename() + "'");
      proj = null;
      loadProject();
    } else {
    }

    if (proj != null
        && timestampOfSampleDataFile < new File(proj.SAMPLE_DATA_FILENAME.getValue(false,
                                                                                   false)).lastModified()) {
      log.report("Detected a change in the sampleData file; reloading sample data");
      proj.resetSampleData();
    }

    if (command.equals(EXIT)) {
      System.exit(0);
    } else if (command.equals(EDIT)) {
      log.report("Launching project properties editor...");
      final ProjectPropertiesEditor configurator =
                                                 new ProjectPropertiesEditor(proj,
                                                                             ProjectPropertiesEditor.ALL_PROPERTY_SETS);
      configurator.addWindowListener(new WindowAdapter() {
        @Override
        public void windowClosed(WindowEvent e) {
          Launch.this.requestFocus();
          configurator.dispose();
        };
      });
      configurator.setVisible(true);
    } else if (command.equals(REFRESH)) {
      loadProjects();
      loadProject();
      log.report("Refreshed list of projects");
    } else if (command.equals(PIPELINE)) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          final GenvisisWorkflow kAndK = new GenvisisWorkflow(proj, Launch.this);
          kAndK.showDialogAndRun();
        }
      });
    } else if (command.equals(NEW_PROJECT)) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          final GenvisisWorkflow kAndK = new GenvisisWorkflow(null, Launch.this);
          kAndK.showDialogAndRun();
        }
      });
    } else if (command.equals(IMPORT_PROJECT)) {

      ImportProjectGUI importGUI = new ImportProjectGUI();
      importGUI.setModal(true);
      importGUI.setVisible(true);

      if (!importGUI.getCancelled()) {
        if (importGUI.run()) {
          String newFilename = importGUI.getNewProjectFilename();
          loadProjects();
          setIndexOfCurrentProject(newFilename);
          loadProject();
        }
      }
      importGUI.dispose();

    } else if (command.equals(CHECK_FOR_UPDATES)) {

      HttpUpdate.update("http://genvisis.org/genvisis_dev.jar", "./", log);

    } else if (command.endsWith(" ")) {
      for (int i = 0; i < projects.length; i++) {
        if (command.equals(ext.rootOf(projects[i]) + " ")) {
          projectsBox.setSelectedIndex(i);
          log.report("Selecting: " + projects[i]);
        }
      }
    } else {
      thread = new Thread(new IndependentThread(proj, command));
      thread.start();
      threadsRunning.add(thread);
    }
  }

  @Override
  public void windowOpened(WindowEvent we) {}

  @Override
  public void windowClosing(WindowEvent we) {
  }

  @Override
  public void windowClosed(WindowEvent we) {}

  @Override
  public void windowIconified(WindowEvent we) {}

  @Override
  public void windowDeiconified(WindowEvent we) {}

  @Override
  public void windowActivated(WindowEvent we) {}

  @Override
  public void windowDeactivated(WindowEvent we) {}

  @Override
  public void itemStateChanged(ItemEvent e) {
    if (e.getStateChange() == ItemEvent.SELECTED) {
      indexOfCurrentProj = projectsBox.getSelectedIndex();
      loadProject();
      log.report("\nCurrent project: " + ext.rootOf(projects[indexOfCurrentProj]) + "\n");

      launchProperties.setProperty(LaunchProperties.LAST_PROJECT_OPENED,
                                   projects[projectsBox.getSelectedIndex()]);
      launchProperties.save();
    }
  }

  public static String getDefaultDebugProjectFile(boolean verbose) {
    LaunchProperties launchProperties;
    String dir, filename;

    if (Files.exists(LaunchProperties.DEFAULT_PROPERTIES_FILE)) {
      launchProperties = new LaunchProperties(LaunchProperties.DEFAULT_PROPERTIES_FILE);
      dir = launchProperties.getProperty(LaunchProperties.PROJECTS_DIR);
      filename = launchProperties.getProperty(LaunchProperties.DEBUG_PROJECT_FILENAME);
      if (dir == null || filename == null) {
        if (verbose) {
          System.err.println("Warning - you are trying to access the default debug project properties file, but there is no '"
                             + LaunchProperties.DEBUG_PROJECT_FILENAME + "=' property listed in '"
                             + LaunchProperties.DEFAULT_PROPERTIES_FILE
                             + "'. The default filename is being set to \"default.properties\" in the current directory. However, if that does not exist either, then the program will likely end in an error.");
        }
        dir = "./";
        filename = "default.properties";
      } else if (!Files.exists(dir) || !Files.exists(dir + filename)) {
        if (verbose) {
          System.err.println("Error - default debug project properties file does not exist: " + dir
                             + filename);
        }
      } else {
        if (verbose) {
          System.out.println("The default debug project properties file is currently set to '" + dir
                             + filename + "'");
        }
      }
    } else {
      dir = "./";
      filename = "default.properties";
    }

    return dir + filename;
  }

  public static void main(String[] args) {
    // Check for alternate main requests
    String mainClassName = args.length > 0 ? args[0] : null;
    Class<?> mainClass = null;

    if (mainClassName != null) {
      // Check the given class. If it doesn't exist, prepend org.genvisis package and try again
      try {
        mainClass = Class.forName(mainClassName);
      } catch (ClassNotFoundException exc) {
        mainClassName = "org.genvisis." + mainClassName;
      }

      // Try again with the updated package name
      if (mainClass == null) {
        try {
          mainClass = Class.forName(mainClassName);
        } catch (ClassNotFoundException exc) {
          // Requested class not found
          System.err.println("Requested main class not found: " + mainClassName);
        }
      }

      // If we found a main class, try running it
      if (mainClass != null) {
        Method meth;
        try {
          meth = mainClass.getMethod("main", String[].class);
          String[] params = Arrays.copyOfRange(args, 1, args.length);
          meth.invoke(null, (Object) params);
        } catch (NoSuchMethodException exc) {
          System.err.println("Requested main class does not have main method: " + mainClassName);
        } catch (Exception exc) {
          if (exc instanceof RuntimeException) {
            throw new RuntimeException(exc);
          }
          System.err.println(exc.getMessage());
        }
      }
      // Whether it worked or not, if we were given a main class parameter we skip the UI.
      return;
    }

    try {
      System.out.println(ext.getTime() + "]\tStarting Genvisis...");
      createAndShowGUI();
      System.out.println(ext.getTime() + "]\tGenvisis Loaded.");
    } catch (InternalError e) {
      if (e.getMessage().contains("X11")) {
        System.err.println("Error occurred with X11 forwarding - please install an X11 forwarding server (we recommend Xming - http://sourceforge.net/projects/xming/) or check your X11 forwarding configuration");
      }
    }
  }

  /**
   * @return launch for windows
   */
  public static String getLaunchBat() {
    String bat = "#This script is intended for launch on Windows machines\r\n";
    bat += "#-Xmx2000m indicates 2000 mb of memory, adjust number up or down as needed\r\n";
    bat += "#Script must be in the same directory as vis.jar\r\n";
    bat += "for %%x in (%0) do set BatchPath=%%~dpsx\r\n";
    bat += "for %%x in (%BatchPath%) do set BatchPath=%%~dpsx\r\n";
    bat += "java  -Xmx22000m -cp %BatchPath%/genvisis.jar cnv.Launch  %*\r\n";
    bat += "PAUSE\r\n";
    return bat;
  }

  /**
   * @return launch for linux and mac
   */
  public static String getLaunchSH() {
    String sh = "#!/bin/sh\n";
    sh += "#This script is intended for launch on *nix machines\n";
    sh += "#-Xmx2000m indicates 2000 mb of memory, adjust number up or down as needed\n";
    sh += "#Script must be in the same directory as vis.jar\n";
    sh += "prefix=`dirname $(readlink $0 || echo $0)`\n";
    sh += "exec java -Xmx2000m \\\n";
    sh += "	-cp \"$prefix\"/genvisis.jar cnv.Launch \"$@\"\n";
    return sh;
  }

  public LaunchProperties getLaunchProperties() {
    return launchProperties;
  }

}
