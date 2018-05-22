package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.stream.Collectors;
import javax.swing.AbstractAction;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.ChromosomeViewer;
import org.genvisis.cnv.gui.CompConfig;
import org.genvisis.cnv.gui.FileNavigator;
import org.genvisis.cnv.gui.ListEditor;
import org.genvisis.cnv.gui.MedianLRRWidget;
import org.genvisis.cnv.gui.RegionNavigator;
import org.genvisis.cnv.gui.RegionNavigator.ChrNavigator;
import org.genvisis.cnv.gui.UITools;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.manage.UCSCtrack;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.var.CNVRectangles;
import org.genvisis.cnv.var.Region;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariantHash;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.Segment;
import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableSet;

/**
 * @author Michael Vieths
 */
public class CompPlot extends JFrame implements ChrNavigator {

  public static final long serialVersionUID = 1L;

  // public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32
  // public static final String DEFAULT_LOCATION = "chr17:15,609,472-40,824,368"; // USP32
  // public static final String DEFAULT_LOCATION = "chr6:161,590,461-163,364,497"; // PARK2
  public static final String DEFAULT_LOCATION = "chr6:161,624,000-163,776,000"; // PARK2 region
  private static final String REGION_LIST_NEW_FILE = "Load Region File";

  public static Color[] colorScheme = {Color.RED, Color.GREEN, Color.BLUE, Color.MAGENTA,
                                       Color.CYAN, Color.ORANGE, Color.YELLOW};

  Project proj;
  private String[] files;
  ArrayList<String> allFiles;
  List<String> filterFiles;
  GeneTrack track;

  // UI Components
  public JPanel compView;
  public RegionNavigator regionNavigator;
  public FileNavigator fileNavigator;
  public CompConfig compConfig;
  public ChromosomeViewer chromosomeViewer;
  public CompPropertyChangeListener cpcl;
  public CompPanel compPanel;

  // Variables configured via subpanels
  // From CompConfig
  int probes;
  int minSize;
  int qualityScore;
  int rectangleHeight;
  String displayMode;
  boolean showExcludes = false;

  CNVRectangles cnvRects;

  // From RegionNavigator
  private Segment location = new Segment((byte) 0, 0, 0);
  private final MarkerDetailSet markerSet;
  String[] allSamples;
  String[] subSamples;
  private final Set<Marker> dropped;

  ArrayList<CNVariantHash> hashes;
  private JMenu delRegionFileMenu;
  private JMenu loadRecentFileMenu;
  private JMenu editRegionFileMenu;
  private ButtonGroup regionButtonGroup;
  private final HashMap<String, JCheckBoxMenuItem> regionFileNameBtn = new HashMap<>();
  private final HashMap<String, String> regionFileNameLoc = new HashMap<>();
  private String[] originalRegionFiles = null;

  public CompPlot(Project proj) {
    super("Genvisis - CompPlot - " + proj.PROJECT_NAME.getValue());
    this.proj = proj;
    markerSet = this.proj.getMarkerSet();
    if (markerSet != null) {
      dropped = proj.getFilteredHash().keySet().stream().map(markerSet.getMarkerNameMap()::get)
                    .filter(Predicates.notNull()).collect(ImmutableSet.toImmutableSet());
    } else {
      dropped = ImmutableSet.of();
    }
    allSamples = proj.getSamples();
    subSamples = ArrayUtils.subArray(proj.getSamples(),
                                     ArrayUtils.booleanNegative(proj.getSamplesToExclude()));
    init();

    addWindowListener(new WindowAdapter() {

      @Override
      public void windowClosing(WindowEvent e) {
        String[] curr = CompPlot.this.proj.REGION_LIST_FILENAMES.getValue();
        HashSet<String> currSet = new HashSet<>();
        for (String s : curr) {
          currSet.add(s);
        }
        for (String s : originalRegionFiles) {
          currSet.remove(s);
        }
        if (currSet.size() > 0) {
          String message = currSet.size() + " files have been added.  ";
          int choice = JOptionPane.showOptionDialog(null,
                                                    message
                                                          + " Would you like to keep this configuration for the next time CompPlot is loaded?",
                                                    "Preserve CompPlot workspace?",
                                                    JOptionPane.YES_NO_CANCEL_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE, null, null, null);
          if (choice == 0) {
            CompPlot.this.proj.saveProperties();
          }
        }
      }
    });
  }

  private void init() {
    // Position

    // Get a list of the .cnv files
    List<String> fileList = Arrays.stream(proj.CNV_FILENAMES.getValue())
                                  .collect(Collectors.toList());
    if (fileList.size() == 0) {
      // CNV_FILENAMES is empty, throw an error and exit
      JOptionPane.showMessageDialog(null,
                                    "Error - CNV_FILENAMES property is empty. No CNV tracks will be displayed.");
    }

    originalRegionFiles = proj.REGION_LIST_FILENAMES.getValue();

    // Get the GeneTrack
    Resource geneTrack = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
                                  .getGTrack();
    if (geneTrack.isAvailable()) {
      track = GeneTrack.load(geneTrack.get());
    } else {
      JOptionPane.showMessageDialog(this,
                                    "Gene track is not installed. Gene boundaries will not be displayed.",
                                    "FYI", JOptionPane.INFORMATION_MESSAGE);
      track = null;
    }

    // Load the variants into memory
    hashes = new ArrayList<>();
    filterFiles = new ArrayList<>();
    Iterator<String> fit = fileList.iterator();
    while (fit.hasNext()) {
      String file = fit.next();
      if (Files.exists(file)) {
        File filename = new File(file);
        filterFiles.add(filename.getName());
        // Load the CNVs out of the files
        CNVariantHash cnvHash = CNVariantHash.load(file, CNVariantHash.CONSTRUCT_ALL,
                                                   proj.getLog());
        hashes.add(cnvHash);
      } else {
        JOptionPane.showMessageDialog(null, "Error - File " + file
                                            + " does not exist, can not display CNV track.");
        fit.remove();
      }
    }

    files = fileList.toArray(new String[fileList.size()]);

    allFiles = new ArrayList<>();
    for (String file2 : files) {
      File file = new File(file2);
      String filename = file.getName();
      allFiles.add(filename);
    }

    setupGUI();
    setJMenuBar(createMenuBar());

    // Initialize the filter attributes
    probes = compConfig.getProbes();
    minSize = compConfig.getMinSize();
    qualityScore = compConfig.getQualityScore();
    rectangleHeight = compConfig.getRectangleHeight();
    setDisplayMode(compConfig.getDisplayMode());

    setPosition(DEFAULT_LOCATION);

    // FIXME update for actual region navigation
    // JCheckBoxMenuItem jcbmi =
    // regionFileNameBtn.get(ext.rootOf(CompPlot.this.regionNavigator.getRegionFile()));
    // if (jcbmi != null) {
    // jcbmi.setSelected(true);
    // }
  }

  public static final int DEFAULT_STARTX = 20;
  public static final int DEFAULT_STARTY = 20;
  public static final int DEFAULT_WIDTH = 1000;
  public static final int DEFAULT_HEIGHT = 720;

  private void setupGUI() {
    // Set the default window size
    setMinimumSize(new Dimension(DEFAULT_STARTX, DEFAULT_STARTY));
    UITools.setSize(this, new Dimension(DEFAULT_WIDTH, DEFAULT_HEIGHT));
    // Close this window but not the entire application on close
    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

    cpcl = new CompPropertyChangeListener(this);

    // Create a new JPanel to contain everything
    compView = new JPanel();
    compView.setLayout(new BorderLayout());

    JPanel topPanel = new JPanel();
    topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));

    regionNavigator = new RegionNavigator(this);
    topPanel.add(regionNavigator);

    fileNavigator = new FileNavigator(files, colorScheme);
    fileNavigator.addPropertyChangeListener(cpcl);
    topPanel.add(fileNavigator);

    compView.add(topPanel, BorderLayout.PAGE_START);

    JPanel viewers = new JPanel();
    viewers.setLayout(new BorderLayout());

    chromosomeViewer = new ChromosomeViewer(location, track);

    viewers.add(chromosomeViewer, BorderLayout.NORTH);
    chromosomeViewer.setPreferredSize(new Dimension(800, 45));

    compPanel = new CompPanel(this);
    compPanel.addPropertyChangeListener(cpcl);
    compPanel.setChromosomeViewer(chromosomeViewer);

    JScrollPane jsp = new JScrollPane(compPanel);
    jsp.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
    viewers.add(jsp, BorderLayout.CENTER);

    compView.add(viewers);

    compConfig = new CompConfig(this);
    compConfig.addPropertyChangeListener(cpcl);

    compView.add(compConfig, BorderLayout.LINE_END);

    add(compView);

    pack();
    // Set the panel visible
    setVisible(true);
  }

  private JMenuBar createMenuBar() {
    JMenuBar menuBar = new JMenuBar();

    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic(KeyEvent.VK_F);

    JMenuItem newRegionFile = new JMenuItem();
    newRegionFile.setAction(new AbstractAction() {

      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        ListEditor le = ListEditor.createRegionListCreator(null,
                                                           proj == null ? null
                                                                        : proj.PROJECT_DIRECTORY.getValue(),
                                                           false);
        le.setModal(true);
        le.setVisible(true);
        if (le.getReturnCode() == JOptionPane.YES_OPTION) {
          final String rgnFile = le.getFileName();
          addFileToList(rgnFile);
          final JMenuItem remove = new JMenuItem();
          remove.setAction(deleteFileAction);
          remove.setText(rgnFile);
          delRegionFileMenu.add(remove);
          // FIXME add actual region navigation
          // regionNavigator.loadRegions();
          // regionNavigator.setRegionFile(rgnFile);
          // regionNavigator.setRegion(0);
          // CompPlot.this.setRegion(regionNavigator.getRegion());
        }
      }
    });
    newRegionFile.setText("New Region List File");
    newRegionFile.setMnemonic(KeyEvent.VK_N);
    fileMenu.add(newRegionFile);

    delRegionFileMenu = new JMenu();
    delRegionFileMenu.setMnemonic(KeyEvent.VK_D);
    delRegionFileMenu.setText("Delete Region List File...");

    for (final String str : proj.REGION_LIST_FILENAMES.getValue()) {
      final JMenuItem remove = new JMenuItem();
      remove.setActionCommand(str);
      remove.setAction(deleteFileAction);
      remove.setText(str);
      delRegionFileMenu.add(remove);
    }

    fileMenu.add(delRegionFileMenu);

    JMenuItem loadRegionFile = new JMenuItem();
    loadRegionFile.setAction(loadNewFileAction);
    loadRegionFile.setText(REGION_LIST_NEW_FILE);
    loadRegionFile.setMnemonic(KeyEvent.VK_L);
    fileMenu.add(loadRegionFile);
    loadRecentFileMenu = new JMenu("Load Recent Region List...");
    loadRecentFileMenu.setMnemonic(KeyEvent.VK_R);
    fileMenu.add(loadRecentFileMenu);

    editRegionFileMenu = new JMenu("Edit Region File");
    editRegionFileMenu.setMnemonic(KeyEvent.VK_E);
    fileMenu.add(editRegionFileMenu);

    menuBar.add(fileMenu);

    regionButtonGroup = new ButtonGroup();
    if (proj != null) {
      String[] files = proj.REGION_LIST_FILENAMES.getValue();
      String name;
      for (String file : files) {
        name = ext.rootOf(file);
        regionFileNameLoc.put(name, file);
        JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
        menuItem.setAction(markerFileSelectAction);
        menuItem.setText(name);
        regionFileNameBtn.put(name, menuItem);
        regionButtonGroup.add(menuItem);
        loadRecentFileMenu.add(menuItem);

        JMenuItem editItem = new JMenuItem();
        editItem.setAction(editFileAction);
        editItem.setText(name);
        editItem.setActionCommand(file);
        editRegionFileMenu.add(editItem);

      }
    }

    JMenu disp = new JMenu("Display");
    disp.setMnemonic(KeyEvent.VK_D);
    menuBar.add(disp);

    JCheckBoxMenuItem displayExcludes = new JCheckBoxMenuItem();
    displayExcludes.setAction(new AbstractAction() {

      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent arg0) {
        setShowExcludes(((JCheckBoxMenuItem) arg0.getSource()).isSelected());
      }
    });
    displayExcludes.setText("Display Excludes?");
    displayExcludes.setSelected(showExcludes);
    displayExcludes.setMnemonic(KeyEvent.VK_E);
    disp.add(displayExcludes);

    JMenu act = new JMenu("Actions");
    act.setMnemonic(KeyEvent.VK_A);
    menuBar.add(act);

    JMenuItem ucsc;
    JMenuItem bedUcsc;
    JMenuItem medianLRR;
    JMenuItem openTrailer;

    ucsc = new JMenuItem();
    ucsc.setAction(ucscAction);
    ucsc.setText("Open Region in UCSC");
    if (Desktop.isDesktopSupported()) {
      ucsc.setToolTipText("View this location on UCSC in a browser");
      ucsc.setEnabled(true);
    } else {
      ucsc.setToolTipText("Browser operations are not supported");
      ucsc.setEnabled(false);
    }
    act.add(ucsc);

    bedUcsc = new JMenuItem();
    bedUcsc.setAction(ucscBedAction);
    bedUcsc.setText("Upload to UCSC");
    if (Desktop.isDesktopSupported()) {
      bedUcsc.setToolTipText("Generate and upload a .BED file to UCSC");
      bedUcsc.setEnabled(true);
    } else {
      bedUcsc.setToolTipText("Browser operations are not supported");
      bedUcsc.setEnabled(false);
    }
    act.add(bedUcsc);

    medianLRR = new JMenuItem();
    medianLRR.setAction(lrrCompAction);
    medianLRR.setText("Median LRR");
    medianLRR.setToolTipText("Compute median Log R Ratios for a region");
    act.add(medianLRR);

    openTrailer = new JMenuItem();
    openTrailer.setAction(openTrailerAction);
    openTrailer.setText("Open all in Trailer");
    openTrailer.setToolTipText("Open the Trailer plot view of all currently visible CNVs");
    act.add(openTrailer);

    return menuBar;
  }

  /**
   * Call {@link #openCNVsInTrailer(List)} with a list of all {@link CNVariant}s in this plot.
   */
  public void openCNVsInTrailer() {
    openCNVsInTrailer(cnvRects.getCNVs());
  }

  /**
   * Open a {@link Trailer} plot with a region list containing an entry for each {@link CNVariant}
   * in the given list.
   */
  public void openCNVsInTrailer(List<CNVariant> selectedCNVs) {
    SampleData sampleData = proj.getSampleData(true);
    int window = proj.getProperty(proj.WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER);

    String[][] sampleRegions = new String[selectedCNVs.size()][];
    for (int i = 0; i < sampleRegions.length; i++) {
      CNVariant cnv = selectedCNVs.get(i);
      String markerPosition = "chr" + location.getChr() + ":" + (cnv.getStart() - window) + "-"
                              + (cnv.getStop() + window);

      // Strip p or q from the end
      if (markerPosition.endsWith("p") || markerPosition.endsWith("q")) {
        markerPosition = markerPosition.substring(0, markerPosition.length() - 1);
      }

      String trailerID = cnv.getFamilyID() + "\t" + cnv.getIndividualID();
      String[] ids = sampleData.lookup(trailerID);
      if (ids != null && ids.length > 0) {
        String id = ids[0];
        sampleRegions[i] = new String[] {id, markerPosition};
      }
    }
    if (sampleRegions.length == 0) {
      proj.message("Error - no valid selected CNVs; cannot launch Trailer without knowing which DNA sample to open");
    } else {
      Trailer t = new Trailer(proj, proj.CNV_FILENAMES.getValue(), sampleRegions);
      t.setVisible(true);
    }
  }

  /**
   * Open a {@link Trailer} plot for each visible CNV.
   */
  private final AbstractAction openTrailerAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      SwingUtilities.invokeLater(new Runnable() {

        @Override
        public void run() {
          openCNVsInTrailer();
        }
      });
    }
  };

  private final AbstractAction lrrCompAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      new Thread(new MedianLRRWidget(proj, regionNavigator.getChrText())).start();
    }
  };

  private final AbstractAction ucscBedAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      // Figure out which files are selected
      // Only allow upload if one file is selected (JDialog warning if multiples)
      List<String> files = getFilterFiles();
      if (files.size() != 1) {
        JOptionPane.showMessageDialog(null,
                                      "One and only one file must be selected before a .BED File can be generated",
                                      "Error", JOptionPane.ERROR_MESSAGE);
      } else {
        // Find the full path to the selected file
        String[] filePaths = proj.CNV_FILENAMES.getValue();
        String compressedFile = "";
        for (String file : filePaths) {
          if (file.endsWith(files.get(0))) {
            System.out.println("File path is " + file);
            compressedFile = file + ".bed..gz"; // TODO this should be .bed.gz?
            // Generate BED file with:
            UCSCtrack.makeTrack(file, file, proj.getLog());
            break;
          }
        }

        // Direct the user to the BED upload page at UCSC Genome Browser
        Desktop desktop = Desktop.getDesktop();
        String URL = Positions.getUCSCUploadLink(Positions.parseUCSClocation(regionNavigator.getChrText()),
                                                 compressedFile);

        // UCSC uses chrX and chrY instead of 23 and 24
        URL = URL.replaceAll("chr23", "chrX");
        URL = URL.replaceAll("chr24", "chrY");
        try {
          URI uri = new URI(URL);
          System.out.println("Browsing to " + URL);
          desktop.browse(uri);
        } catch (URISyntaxException e1) {
          e1.printStackTrace();
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
    }
  };

  private final AbstractAction ucscAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent arg0) {
      Desktop desktop = Desktop.getDesktop();
      String URL = Positions.getUCSClink(Positions.parseUCSClocation(regionNavigator.getChrText()));

      // UCSC uses chrX and chrY instead of 23 and 24
      URL = URL.replaceAll("chr23", "chrX");
      URL = URL.replaceAll("chr24", "chrY");
      try {
        URI uri = new URI(URL);
        System.out.println("Browsing to " + URL);
        desktop.browse(uri);
      } catch (URISyntaxException e) {
        e.printStackTrace();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  };

  private final AbstractAction deleteFileAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      int opt = JOptionPane.showConfirmDialog(CompPlot.this, "Delete file from disk?",
                                              "Delete region file?",
                                              JOptionPane.YES_NO_CANCEL_OPTION,
                                              JOptionPane.QUESTION_MESSAGE, null);
      switch (opt) {
        case JOptionPane.CANCEL_OPTION:
          return;
        case JOptionPane.YES_OPTION:
          boolean deleted = (new File(e.getActionCommand())).delete();
          if (!deleted) {
            JOptionPane.showMessageDialog(CompPlot.this,
                                          "Error - failed to delete file {" + e.getActionCommand()
                                                         + "}",
                                          "Delete File Failed...", JOptionPane.ERROR_MESSAGE);
          }
          break;
        case JOptionPane.NO_OPTION:
          break;
      }
      proj.REGION_LIST_FILENAMES.removeValue(e.getActionCommand());
      proj.REGION_LIST_FILENAMES.getValue();
      // FIXME update with actual region navigation
      // regionNavigator.setRegionFile(val.length > 0 ? val[0] : "");
      // regionNavigator.setRegion(0);
      // CompPlot.this.setRegion(regionNavigator.getRegion());
      delRegionFileMenu.remove((JMenuItem) e.getSource());
      loadRecentFileMenu.remove(regionFileNameBtn.remove(ext.rootOf(e.getActionCommand())));
    }
  };

  private void addFileToList(String rawfile) {
    String file = ext.verifyDirFormat(rawfile);
    file = file.substring(0, file.length() - 1);
    String name = ext.rootOf(file);
    regionFileNameLoc.put(name, file);

    JCheckBoxMenuItem item = new JCheckBoxMenuItem();
    item.setAction(markerFileSelectAction);
    item.setText(name);

    JMenuItem editItem = new JMenuItem();
    editItem.setAction(editFileAction);
    editItem.setText(name);
    editItem.setActionCommand(file);
    editRegionFileMenu.add(editItem);

    regionFileNameBtn.put(name, item);
    regionButtonGroup.add(item);
    loadRecentFileMenu.add(item);

    final JMenuItem remove = new JMenuItem();
    remove.setActionCommand(file);
    remove.setAction(deleteFileAction);
    remove.setText(file);
    delRegionFileMenu.add(remove);

    proj.REGION_LIST_FILENAMES.addValue(file);
    proj.saveProperties(new Property[] {proj.REGION_LIST_FILENAMES});
  }

  private final AbstractAction editFileAction = new AbstractAction() {

    @Override
    public void actionPerformed(ActionEvent e) {
      String file = e.getActionCommand();
      ListEditor editor = ListEditor.createRegionListEditor(null,
                                                            proj == null ? null
                                                                         : proj.PROJECT_DIRECTORY.getValue(),
                                                            false, file);
      editor.setModal(true);
      editor.setVisible(true);

      if (editor.getReturnCode() == JOptionPane.YES_OPTION) {

        String file1 = ext.verifyDirFormat(file);
        file1 = file1.substring(0, file1.length() - 1);
        String name = ext.rootOf(file1);
        regionFileNameBtn.get(name).setSelected(true);
        addFileToList(file1);
      }

    }
  };

  private final AbstractAction markerFileSelectAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      String shortName = ((JCheckBoxMenuItem) e.getSource()).getText();
      // if (!loadingFile) {
      String file = regionFileNameLoc.get(shortName);
      // FIXME update for region files
      if (file == null /* || file.equals(regionNavigator.getRegionFile()) */) {
        return;
      }
      String tempFile = Files.isRelativePath(file) ? proj.PROJECT_DIRECTORY.getValue() + file
                                                   : file;
      if (!Files.exists(tempFile)) {
        proj.message("Error - region file '" + tempFile + "' doesn't exist.");
        regionFileNameBtn.get(shortName).setSelected(true);
      } else {
        // proj.REGION_LIST_FILENAMES.setValue(Array.insertStringAt(file,
        // proj.REGION_LIST_FILENAMES.getValue(), 0));
        // CompPlot.this.regionNavigator.loadRegions();
        // regionNavigator.setRegionFile(file);
        // regionNavigator.setRegion(0);
        // CompPlot.this.setRegion(regionNavigator.getRegion());
        // regionIndex = 0;
        // showRegion();
      }
      /*
       * } else if (loadingFile) { // leave as currently selected marker if
       * (CompPlot.this.regionNavigator.getRegionFile() != "" &&
       * CompPlot.this.regionNavigator.getRegionFile() != null) { String file =
       * ext.rootOf(CompPlot.this.regionNavigator.getRegionFile());
       * regionFileNameBtn.get(file).setSelected(true); } return; }
       */
    }
  };

  AbstractAction loadNewFileAction = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      String newFile = chooseNewFiles();
      if (newFile != null) {
        // FIXME update for region files
        // regionNavigator.loadRegions();
        // regionNavigator.setRegionFile(newFile);
        // regionNavigator.setRegion(0);
        // regionFileNameBtn.get(ext.rootOf(regionNavigator.getRegionFile())).setSelected(true);
        // CompPlot.this.setRegion(regionNavigator.getRegion());
      }
    }
  };

  private String chooseNewFiles() {
    JFileChooser jfc = new JFileChooser((proj != null/* || regionFileName == null */ ? proj.PROJECT_DIRECTORY.getValue()
                                                     : null /*
                                                             * ext.parseDirectoryOfFile(
                                                             * regionFileName)
                                                             */));
    jfc.setMultiSelectionEnabled(true);
    if (jfc.showOpenDialog(CompPlot.this) == JFileChooser.APPROVE_OPTION) {
      File[] files = jfc.getSelectedFiles();
      if (files.length > 0) {
        boolean[] keep = ArrayUtils.booleanArray(files.length, true);
        for (int i = 0; i < files.length; i++) {
          for (String fileName : proj.REGION_LIST_FILENAMES.getValue()) {
            if (ext.rootOf(files[i].toString()).equals(fileName)) {
              keep[i] = false;
            }
          }
        }
        File[] keptFiles = ArrayUtils.subArray(files, keep);
        File[] discards = ArrayUtils.subArray(files, ArrayUtils.booleanNegative(keep));

        if (discards.length > 0) {
          StringBuilder msg = new StringBuilder("The following data file(s) are already present:");
          for (File disc : discards) {
            msg.append("\n").append(disc.getName());
          }
          JOptionPane.showMessageDialog(CompPlot.this, msg.toString());
        }

        for (File kept : keptFiles) {
          if (verifyValidFile(kept.getAbsolutePath())) {
            addFileToList(kept.getAbsolutePath());
          } else {
            proj.getLog().reportError("Error - contents of file {" + kept.getAbsolutePath()
                                      + "} are not valid UCSC regions");
          }
        }
        return keptFiles[0].getAbsolutePath();
      } else {
        File file = jfc.getSelectedFile();
        boolean keep = true;
        for (String fileName : proj.REGION_LIST_FILENAMES.getValue()) {
          if (ext.rootOf(file.toString()).equals(fileName)) {
            keep = false;
          }
        }

        if (!keep) {
          StringBuilder msg = new StringBuilder("The following data file is already present:\n").append(file.getName());
          JOptionPane.showMessageDialog(CompPlot.this, msg.toString());
        } else {
          if (verifyValidFile(file.getAbsolutePath())) {
            addFileToList(file.getAbsolutePath());
            return file.getAbsolutePath();
          } else {
            proj.getLog().reportError("Error - contents of file {" + file.getAbsolutePath()
                                      + "} are not valid UCSC regions");
            return null;
          }
        }

      }

    }
    return null;
  }

  private boolean verifyValidFile(String file) {
    File f = new File(file);
    if (!f.exists()) {
      return false;
    }
    try {
      BufferedReader reader = Files.getAppropriateReader(file);
      String line = null;
      while ((line = reader.readLine()) != null) {
        String[] parts = line.split("\t");
        int[] pos = Positions.parseUCSClocation(parts[0]);
        if (pos == null || pos[0] == -1 || pos[1] == -1 || pos[2] == -1) {
          return false;
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  public void loadCNVs(Segment location) {
    // long startTime = Calendar.getInstance().getTimeInMillis();
    cnvRects = new CNVRectangles(hashes, allFiles, filterFiles, location, probes, minSize,
                                 qualityScore, proj.getSampleData(false),
                                 showExcludes ? allSamples : subSamples);
    cnvRects.setRectangleHeight(rectangleHeight);
    compPanel.setWindow(location.getStart(), location.getStop());
    cnvRects.setScalingFactor(compPanel.getScalingFactor());
    compPanel.setCNVRectangles(cnvRects);

    // long stopTime = Calendar.getInstance().getTimeInMillis();
    //
    // System.out.println("loadCNVs() took " + (stopTime - startTime) + "ms");
  }

  public CNVRectangles getCnvRects() {
    return cnvRects;
  }

  public String[] getFiles() {
    return files;
  }

  public void setFiles(String[] files) {
    this.files = files;
  }

  public void setShowExcludes(boolean show) {
    showExcludes = show;
    loadCNVs(location);
  }

  /*
   * Methods to set values pulled from CompConfig
   */
  public void setProbes(int p) {
    probes = p;
    loadCNVs(location);
  }

  public void setMinSize(int ms) {
    minSize = ms;
    loadCNVs(location);
  }

  public void setQualityScore(int qs) {
    qualityScore = qs;
    loadCNVs(location);
  }

  public void setRectangleHeight(int height) {
    rectangleHeight = height;
    loadCNVs(location);
  }

  public void setDisplayMode(String dm) {
    displayMode = dm;
    compPanel.setDisplayMode(displayMode);
    loadCNVs(location);
  }

  /**
   * Set the selected CNVs from the main UI, e.g. when clicking on a rectangle
   */
  public void setSelectedCNVs(List<CNVariant> cnvs) {
    compConfig.setSelectedCNVs(cnvs);
  }

  /*
   * Method to set values pulled from RegionNavigator
   */
  @Override
  public void setPosition(String region) {
    Region r = new Region(region, track);
    setRegion(r);
  }

  @Override
  public byte getFirstAvailableChr() {
    return markerSet.getChrs()[0];
  }

  @Override
  public byte getLastAvailableChr() {
    return markerSet.getChrs()[markerSet.getChrs().length - 1];
  }

  public void setRegion(Region region) {
    location = new Segment(region.getRegion());

    if (location.getChr() == -1) {
      return;
    }
    if (location.getStart() < 0) {
      location = new Segment(location.getChr(), 1, location.getStop());
    }

    NavigableSet<Marker> chrMarkers = markerSet.getChrMap().get(location.getChr());
    int lastPosition = chrMarkers.isEmpty() ? 0 : chrMarkers.last().getPosition();
    if (location.getStop() == -1 || location.getStop() > lastPosition) {
      location = new Segment(location.getChr(), location.getStart(), lastPosition);
    }
    regionNavigator.setChrFieldText(location);
    chromosomeViewer.updateView(location);
    loadCNVs(location);
    chromosomeViewer.repaint();
  }

  public Region getRegion() {
    return new Region(regionNavigator.getChrText());
  }

  public void setCPLocation(Segment location) {
    this.location = location;
    regionNavigator.setChrFieldText(location);
    chromosomeViewer.updateView(location);
    loadCNVs(location);
    chromosomeViewer.repaint();
  }

  public Segment getCPLocation() {
    return location;
  }

  public Project getProject() {
    return proj;
  }

  public void setFilter(List<String> files) {
    filterFiles = files;
    loadCNVs(location);
  }

  public List<String> getFilterFiles() {
    return filterFiles;
  }

}

class CompPropertyChangeListener implements PropertyChangeListener {

  CompPlot compPlot;

  public CompPropertyChangeListener(CompPlot cp) {
    compPlot = cp;
  }

  @Override
  public void propertyChange(PropertyChangeEvent pve) {
    String propertyName;

    propertyName = pve.getPropertyName();

    if (propertyName.equals("probes")) {
      compPlot.setProbes(Integer.parseInt(pve.getNewValue().toString()));
    } else if (propertyName.equals("minSize")) {
      compPlot.setMinSize(Integer.parseInt(pve.getNewValue().toString()));
    } else if (propertyName.equals("qualityScore")) {
      compPlot.setQualityScore(Integer.parseInt(pve.getNewValue().toString()));
    } else if (propertyName.equals("rectangleHeight")) {
      compPlot.setRectangleHeight(Integer.parseInt(pve.getNewValue().toString()));
    } else if (propertyName.equals("displayMode")) {
      compPlot.setDisplayMode((String) pve.getNewValue());
      // } else if (propertyName.equals("firstRegion")) {
      // } else if (propertyName.equals("previousRegion")) {
      // } else if (propertyName.equals("nextRegion")) {
      // } else if (propertyName.equals("lastRegion")) {
    } else if (propertyName.equals("selectedCNV")) {
      @SuppressWarnings("unchecked")
      ArrayList<CNVariant> cnvs = (ArrayList<CNVariant>) pve.getNewValue();
      compPlot.setSelectedCNVs(cnvs);
    } else if (propertyName.equals("files")) {
      @SuppressWarnings("unchecked")
      ArrayList<String> files = (ArrayList<String>) pve.getNewValue();
      compPlot.setFilter(files);
    } else {
      // System.out.println(pve.getPropertyName() + " changed from " + pve.getOldValue() + " to " +
      // pve.getNewValue());
    }
  }
}
