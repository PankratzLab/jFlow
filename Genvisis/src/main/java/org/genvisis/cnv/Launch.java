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
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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

import org.apache.commons.io.FileUtils;
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

public class Launch extends JFrame implements ActionListener, WindowListener {

  public static final long serialVersionUID = 1L;

  public static final String VERSION = "0.60";

  public static final String EXIT = "Exit";
  public static final String EDIT = "Project Properties Editor";
  public static final String REFRESH = "Refresh";
  public static final String PIPELINE = "Genvisis Project Workflow";
  public static final String NEW_PROJECT = "New Project";
  public static final String IMPORT_PROJECT = "Import Project";
  public static final String SELECT_PROJECT = "Select Project";
  public static final String DELETE_PROJECT = "Delete Project";

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

  public static final String[][] MENUS = {
                                          {"File", NEW_PROJECT, IMPORT_PROJECT, SELECT_PROJECT,
                                           DELETE_PROJECT, EDIT, "Preferences", CHECK_FOR_UPDATES,
                                           EXIT},
                                          {"Data", MAP_FILES, GENERATE_MARKER_POSITIONS,
                                           PARSE_FILES_CSV, TRANSPOSE_DATA, PIPELINE}, // ,
                                                                                       // MITOPIPELINE
                                          {"Quality", CHECK_SEX, LRR_SD, CNP_SCAN, MOSAICISM,
                                           MARKER_METRICS, FILTER_MARKER_METRICS,
                                           TALLY_MARKER_ANNOTATIONS,
                                           TALLY_WITHOUT_DETERMINING_DROPS, TALLY_CLUSTER_FILTERS},
                                          {"Plots", SCATTER, QQ, STRAT, MOSAIC_PLOT, SEX_PLOT,
                                           TRAILER, TWOD, LINE_PLOT, COMP, FOREST_PLOT},
                                          {"Tools", GENERATE_ABLOOKUP, EXPORT_TO_PLINK,
                                           GENERATE_PENNCNV_FILES, PARSE_RAW_PENNCNV_RESULTS, POPULATIONBAF, GCMODEL, CUSTOM_CENTROIDS, DENOVO_CNV, EXPORT_CNVS, CYTO_WORKBENCH, PRINCIPAL_COMPONENTS, GENERATE_DEMO_PACKAGE, ADD_QC_TO_SAMPLE_DATA, TEST},
                                          {"Help", "Contents", "Search", "About"}};


  private Project proj;
  private final boolean jar;
  private JComboBox<String> projectsBox;
  private transient List<String> projects;
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
    log = new Logger();
  }

  /**
   * Safely initialize the projects list once, and only once.
   */
  private synchronized void initProjects() {
    if (projects == null) {
      List<String> list = Arrays.asList(launchProperties.getListOfProjectProperties());
      projects = list;
    }
  }

  /**
   * Discover the projects declared in the projects directory and the last opened project.
   */
  public void loadProjects() {
    if (projects == null) {
      initProjects();
    }
    setIndexOfCurrentProject(launchProperties.getProperty(LaunchProperties.LAST_PROJECT_OPENED));
  }

  public Project loadProject() {
    // TODO if indexOfCurrentProj < 0 show dialog asking if user would like to create new project
    proj = new Project(launchProperties.getDirectory() + projects.get(indexOfCurrentProj), jar);
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

  public void setIndexOfCurrentProject(int projectIndex) {
    setIndexOfCurrentProject(projects.get(projectIndex));
  }

  public void setIndexOfCurrentProject(String projPropertiesFileName) {
    indexOfCurrentProj = 0;
    // find the index
    for (int i = 0; i < projects.size(); i++) {
      if (projects.get(i).equals(projPropertiesFileName)) {
        indexOfCurrentProj = i;
      }
    }
    // Update the UI component if available
    if (projectsBox != null && !projects.isEmpty()) {
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

  /**
   * Entry point for graphical use
   */
  private static void createAndShowGUI() {
    String launchPropertiesFile;
    final Launch launchUI;

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
    CurrentManifest manifest = new CurrentManifest(new Attributes());
    try {// try not to break the launch so we will catch anything
      manifest = CurrentManifest.loadGenvisisManifest();
    } catch (Exception e) {
    }
    launchUI = new Launch(launchPropertiesFile, manifest, false);
    launchUI.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    // Have to read which projects are available before creating the menus
    launchUI.initLaunchProperties();

    launchUI.loadProjects();

    launchUI.makeContentPane();
    launchUI.makeTopMenuBar();
    launchUI.makePopupMenu();

    launchUI.setSize(650, 550);
    launchUI.setLocation(300, 200);

    launchUI.addWindowListener(launchUI);

    javax.swing.SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        launchUI.setVisible(true);
      }
    });

    // TODO only instantiate when used
    launchUI.setIndexOfCurrentProject(launchUI.launchProperties.getProperty(LaunchProperties.LAST_PROJECT_OPENED));
    if (launchUI.projects.size() > 0) {
      launchUI.loadProject();
    }

    ueh.setLog(launchUI.log);
    try {
      launchUI.log.report(HttpUpdate.checkGenvisisVersion(launchUI.log));
    } catch (Exception e) {

    }
  }

  private void initLaunchProperties() {
    launchProperties = createLaunchProperties(launchPropertiesFile, false, true);
  }

  private void makeContentPane() {
    // Create the content-pane-to-be.
    JPanel contentPane = new JPanel(new BorderLayout());
    contentPane.setOpaque(true);

    contentPane.add(makeTopIconBar(), BorderLayout.NORTH);

    // Create a scrolled text area.
    output = new JTextArea(5, 30);
    output.setEditable(false);
    log.linkTextArea(output);
    scrollPane = new JScrollPane(output);

    // Add the text area to the content pane.
    contentPane.add(scrollPane, BorderLayout.CENTER);

    progBar = new JProgressBar();
    contentPane.add(progBar, BorderLayout.SOUTH);
    progBar.setVisible(false);

    setContentPane(contentPane);
  }

  private void makeTopMenuBar() {
    JMenuBar menuBar;
    JMenu menu;
    JMenuItem menuItem = null;
    Set<Character> hash;

    menuBar = new JMenuBar();
    // attach mnemonics and actionlisteners to menu elements
    for (String[] element : MENUS) {
      menu = new JMenu(element[0]);
      menu.setMnemonic((int) element[0].charAt(0));
      menuBar.add(menu);
      hash = new HashSet<Character>();
      for (int j = 1; j < element.length && !element[j].isEmpty(); j++) {
        if (element[j] == "1") {
          menu.addSeparator();
          continue;
        } else if (element[j].equals(SELECT_PROJECT)) {
          // Create "select project" submenu
          menuItem = new JMenu(element[j]);
          for (String project : projects) {
            String label = ext.rootOf(project, true) + " ";
            JMenuItem subItem = new JMenuItem(label);
            subItem.addActionListener(this);
            subItem.setMnemonic(getMnemonic(label, hash));
            menuItem.add(subItem);
          }
        } else if (element[j].equals(PRINCIPAL_COMPONENTS)) {
          // Create "principal components" submenu
          menuItem = new JMenu(element[j]);
          for (String pcSubMenuOption : new String[] {PrincipalComponentsManhattan.PRINCIPAL_MANHATTAN_MI,
                                                      PrincipalComponentsCrossTabs.PRINCIPAL_CROSSTABS_MI}) {
            JMenuItem pcSubItem = new JMenuItem(pcSubMenuOption);
            pcSubItem.addActionListener(this);
            menuItem.add(pcSubItem);
          }
        } else {
          // standard menu item
          menuItem = new JMenuItem(element[j]);
          menuItem.addActionListener(this);
        }

        menuItem.setMnemonic(getMnemonic(element[j], hash));
        menu.add(menuItem);
      }
    }
    setJMenuBar(menuBar);
  }

  private int getMnemonic(String string, Set<Character> hash) {
    if (string.equals(NEW_PROJECT)) {
      return KeyEvent.VK_N;
    } else if (string.equals(IMPORT_PROJECT)) {
      return KeyEvent.VK_I;
    } else if (string.equals(DELETE_PROJECT)) {
      return KeyEvent.VK_D;
    } else {
      String s2 = string.toLowerCase();
      for (int k = 0; k < string.length(); k++) {
        char mnemonic = s2.charAt(k);
        if (!hash.contains(mnemonic)) {
          hash.add(mnemonic);
          return mnemonic;
        }
      }
    }
    return 0;
  }

  /**
   * Add icon buttons for various operations
   */
  private JPanel makeTopIconBar() {
    JPanel iconBar = new JPanel();
    iconBar.setLayout(new FlowLayout(FlowLayout.LEFT));

    // Add buttons to left of project selector
    addButtons(iconBar,
               new String[] {"images/save1.png", "images/edit1.png", "images/refresh.gif",
                             "images/gen_pipe_1.png", "images/scatterPlot2.png",
                             "images/trailerPlot2.png", "images/qqplot.gif",
                             "images/recluster1.png", "images/twoDPlot1.jpg",
                             "images/forestPlot1.png"},
               new String[] {"", EDIT, REFRESH, PIPELINE, SCATTER, TRAILER, QQ, LINE_PLOT, TWOD,
                             FOREST_PLOT});

    addProjectSelector(iconBar);

    // Add buttons to right of project selector
    addButtons(iconBar, new String[] {"images/deleteProj.svg.png"}, new String[] {DELETE_PROJECT});

    return iconBar;
  }

  /**
   * Helper method to create a set of standardized buttons to the given pane. {@code icons.length}
   * must equal {@code commands.length}, which is also the number of buttons that will be created.
   */
  private void addButtons(final Container pane, String[] icons, String[] commands) {
    if (icons.length != commands.length) {
      throw new IllegalArgumentException("Error creating topbar buttons. Got " + icons.length
                                         + " icons but " + commands.length + "commands.");
    }

    for (int i = 0; i < icons.length; i++) {
      JButton button = new JButton(Grafik.getImageIcon(icons[i]));
      button.setActionCommand(commands[i]);
      button.addActionListener(this);
      button.setToolTipText(commands[i]);
      button.setPreferredSize(new Dimension(25, 25));
      button.setOpaque(false);
      button.setContentAreaFilled(false);
      button.setBorderPainted(false);
      pane.add(button);
    }
  }

  /**
   * Create project selector combobox and add it to the given pane
   */
  private void addProjectSelector(final Container pane) {

    projectsBox = new JComboBox<String>();
    // In JDK1.4 this prevents action events from being fired when the up/down arrow keys are used
    // on the dropdown menu
    projectsBox.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
    projectsBox.setModel(new DefaultComboBoxModel<String>(launchProperties.getListOfProjectNames()));

    if (indexOfCurrentProj > 0 && projectsBox.getItemCount() > 0) {
      projectsBox.setSelectedIndex(indexOfCurrentProj);
    }

    projectsBox.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent e) {
        if (e.getStateChange() == ItemEvent.SELECTED) {
          indexOfCurrentProj = projectsBox.getSelectedIndex();
          loadProject();
          log.report("\nCurrent project: " + ext.rootOf(projects.get(indexOfCurrentProj)) + "\n");

          launchProperties.setProperty(LaunchProperties.LAST_PROJECT_OPENED,
                                       projects.get(projectsBox.getSelectedIndex()));
          launchProperties.save();
        }
      }
    });
    pane.add(projectsBox);
  }

  private void makePopupMenu() {
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

    log.report("Action performed: " + command + "\n");
    // TODO in Java 7 we can make these a switch statement

    // These options do not require an active project
    if (command.equals(EXIT)) {
      System.exit(0);
    } else if (command.equals(NEW_PROJECT)) {
      createProject();
      return;
    } else if (command.equals(IMPORT_PROJECT)) {
      importProject();
      return;
    } else if (indexOfCurrentProj < 0 || indexOfCurrentProj >= projects.size()) {
      // Command requested requires an active project, so ensure
      // current project is valid
      log.report("No project currently selected. Attempting to create one now.");
      int i = JOptionPane.showOptionDialog(null,
                                           "No projects available. You can create or import one now.",
                                           "Create project?", JOptionPane.OK_CANCEL_OPTION,
                                           JOptionPane.INFORMATION_MESSAGE, null,
                                           new Object[] {"Create", "Import"}, "Create");
      if (i == 0) {
        createProject();
      } else if (i == 1) {
        importProject();
      }
      // TODO it would be nice if we could wait for create/import to finish before deciding if
      // we should return or continue execution
      return;
    }

    refreshTimestamps();

    // These options require an active project
    if (command.equals(EDIT)) {
      log.report("Launching project properties editor...");
      final ProjectPropertiesEditor configurator =
                                                 new ProjectPropertiesEditor(proj,
                                                                             ProjectPropertiesEditor.ALL_PROPERTY_SETS);
      configurator.addWindowListener(new WindowAdapter() {
        @Override
        public void windowClosed(WindowEvent e) {
          Launch.this.requestFocus();
          configurator.dispose();
        }
      });
      configurator.setVisible(true);
    } else if (command.equals(DELETE_PROJECT)) {
      String toDelete = projects.get(indexOfCurrentProj);
      int newIndex = Math.max(0, --indexOfCurrentProj);
      if (new File(launchProperties.getDirectory() + toDelete).delete()) {
        int deleteDir = JOptionPane.showConfirmDialog(null, "Delete source directory as well?",
                                                      "Delete project",
                                                      JOptionPane.YES_NO_CANCEL_OPTION,
                                                      JOptionPane.WARNING_MESSAGE);
        projects = null;

        // Update toDelete to just the project name
        toDelete = ext.rootOf(toDelete, true);
        // Find and remove the project entry in the "select projects" menu
        deleteMenuItem(getJMenuBar(), "File", SELECT_PROJECT, toDelete);

        if (deleteDir == JOptionPane.YES_OPTION) {
          log = new Logger();
          log.linkTextArea(output);
          try {
            FileUtils.deleteDirectory(new File(proj.PROJECT_DIRECTORY.getValue()));
          } catch (IOException e) {
            log.reportTimeWarning("Failed to delete project directory: "
                                + proj.PROJECT_DIRECTORY.getValue());
          }
        }
        loadProjects();
        projectsBox.setModel(new DefaultComboBoxModel<String>(launchProperties.getListOfProjectNames()));
        if (!projects.isEmpty()) {
          loadProject();
          setIndexOfCurrentProject(newIndex);
        }
      } else {
        log.reportTimeWarning("Failed to delete current project: " + toDelete);
      }
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
    } else if (command.equals(CHECK_FOR_UPDATES)) {

      HttpUpdate.update("http://genvisis.org/genvisis_dev.jar", "./", log);

    } else if (command.endsWith(" ")) {
      // FIXME this should be unified with the drop down combobox selector
      for (int i = 0; i < projects.size(); i++) {
        if (command.equals(ext.rootOf(projects.get(i)) + " ")) {
          projectsBox.setSelectedIndex(i);
          log.report("Selecting: " + projects.get(i));
        }
      }
    } else {
      thread = new Thread(new IndependentThread(proj, command));
      thread.start();
      threadsRunning.add(thread);
    }
  }

  /**
   * Recursively search for the menu path specified by the given string array, removing the final
   * entry from the given menu. Each entry except the last is assumed to be a sub menu.
   */
  private void deleteMenuItem(JMenuBar jMenuBar, String... entries) {
    int toDelete = -1;
    for (int i = 0; i < jMenuBar.getMenuCount(); i++) {
      final JMenu menu = jMenuBar.getMenu(i);
      if (menu.getText().equals(entries[0])) {
        if (entries.length > 1) {
          // Enter the JMenu recursion
          deleteMenuItem(menu, Arrays.copyOfRange(entries, 1, entries.length));
          break;
        }

        toDelete = i;
        break;
      }
    }
    if (toDelete >= 0) {
      jMenuBar.remove(toDelete);
    }
    jMenuBar.validate();
  }

  /**
   * @see {@link #deleteMenuItem(JMenuBar, String...)
   */
  private void deleteMenuItem(JMenu menu, String... entries) {
    if (entries.length == 0) {
      return;
    }
    int toDelete = -1;
    for (int i = 0; i < menu.getItemCount(); i++) {
      final JMenuItem item = menu.getItem(i);
      if (item.getText().trim().equals(entries[0])) {
        if (entries.length > 1) {
          // Recursive step
          deleteMenuItem((JMenu) item, Arrays.copyOfRange(entries, 1, entries.length));
        } else {
          log.report("Deleting: " + entries[0] + " from select project menu");
          // End of recursion
          toDelete = i;
        }

        break;
      }
    }
    if (toDelete >= 0) {
      menu.remove(toDelete);
    }
  }

  /**
   * Check the timestamps on the project and sample data. If they are out of date, reload as
   * appropriate.
   */
  private void refreshTimestamps() {
    if (timestampOfPropertiesFile < new File(proj.getPropertyFilename()).lastModified()) {
      log.report("Detected a change in the project properties file; reloading from '"
                 + proj.getPropertyFilename() + "'");
      proj = null;
      loadProject();
    }

    if (proj != null
        && timestampOfSampleDataFile < new File(proj.SAMPLE_DATA_FILENAME.getValue(false,
                                                                                   false)).lastModified()) {
      log.report("Detected a change in the sampleData file; reloading sample data");
      proj.resetSampleData();
    }
  }

  private void importProject() {
    ImportProjectGUI importGUI = new ImportProjectGUI();
    importGUI.setModal(true);
    importGUI.setVisible(true);

    if (!importGUI.getCancelled() && importGUI.run()) {
      String newFilename = importGUI.getNewProjectFilename();
      loadProjects();
      setIndexOfCurrentProject(newFilename);
      loadProject();
    }
    importGUI.dispose();
  }

  private void createProject() {
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        final GenvisisWorkflow kAndK = new GenvisisWorkflow(null, Launch.this);
        kAndK.showDialogAndRun();
      }
    });
  }

  // FIXME factor out to subclass to reduce clutter
  @Override
  public void windowOpened(WindowEvent we) {}

  @Override
  public void windowClosing(WindowEvent we) {}

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


  /**
   * Ensures a minimum launch properties file is available, reading an existing file and creating a
   * new properties file if needed/requested.
   *
   * @param launchPropertiesFile Path to the properties file, e.g. "launch.properties"
   * @param force If true, any existing properties file will be disregarded and potentially
   *        overwritten
   * @param relativePath If true, relative paths will be written instead of fully qualified
   */
  public static LaunchProperties createLaunchProperties(String launchPropertiesFile, boolean force,
                                                        boolean relativePath) {
    Logger log = new Logger();
    LaunchProperties launchProperties;
    String path = relativePath ? ""
                               : LaunchProperties.directoryOfLaunchProperties(launchPropertiesFile);

    String projPath;

    if (force || !new File(launchPropertiesFile).exists()) {
      log.reportTime("No launch properties found. Creating default: " + launchPropertiesFile);
      projPath = path + "projects" + File.separatorChar;
      // Create a default launch.properties
      // Set the "example" project as the default opened project
      Files.writeList(new String[] {"LAST_PROJECT_OPENED=example.properties",
                                    "PROJECTS_DIR=" + projPath},
                      launchPropertiesFile);
    }
    launchProperties = new LaunchProperties(launchPropertiesFile);

    createExampleProject(launchProperties, path);

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

    return launchProperties;
  }

  /**
   * Creates the example project if it doesn't already exist
   *
   * @param properties launch properties file to use to determine project location
   * @param path Base directory for the project data
   */
  private static void createExampleProject(LaunchProperties properties, String path) {
    Logger log = new Logger();
    String examplePath = path + "example" + File.separatorChar;
    String exampleProperties = properties.getDirectory() + "example.properties";

    File f = new File(examplePath);
    if (!f.exists()) {
      f.mkdirs();
      log.reportTime("Creating example project: " + examplePath);
    }

    if (!new File(exampleProperties).exists()) {
      log.reportTime("Creating example project properties: " + exampleProperties);
      Files.writeList(new String[] {"PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/",
                                    "SOURCE_DIRECTORY=sourceFiles/"},
                      exampleProperties);
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
    if (runMainClass(args)) {
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
   * Helper method to execute alternate main classes. Can prepend packages, to handle moved classes.
   *
   * @return true if a main class was executed.
   */
  private static boolean runMainClass(String[] args) {
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
          exc.printStackTrace();
        }
        return true;
      }
    }
    return false;
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
