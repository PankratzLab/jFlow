// replace all exits with a proper disposal of all resources
package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
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
import javax.swing.JRadioButton;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpringLayout;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.annotation.AnalysisParams;
import org.genvisis.cnv.annotation.AnnotationFileLoader.QUERY_ORDER;
import org.genvisis.cnv.annotation.AnnotationParser;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.BlastParams;
import org.genvisis.cnv.annotation.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.MarkerGCAnnotation;
import org.genvisis.cnv.filesys.AnnotationCollection;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.ClusterFilter;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.gui.AnnotationAction;
import org.genvisis.cnv.gui.AutoSaveForScatterPlot;
import org.genvisis.cnv.gui.BlastFrame;
import org.genvisis.cnv.gui.ColorKeyPanel;
import org.genvisis.cnv.gui.CycleRadio;
import org.genvisis.cnv.gui.NewMarkerListDialog;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.PlinkMarkerLoader;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.AlleleFreq;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ProgressMonitor;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.stats.CTable;
import org.genvisis.stats.ContingencyTable;
import org.genvisis.stats.Histogram;
import org.genvisis.stats.ProbDist;

import net.miginfocom.swing.MigLayout;

public class ScatterPlot extends /* JPanel */JFrame implements ActionListener, WindowListener {
  public static final long serialVersionUID = 1L;
  public static final byte DEFAULT_SIZE = 8;
  public static final int DEFAULT_GC_THRESHOLD = 15;
  private static final String ALT_UP = "ALT UP";
  private static final String ALT_DOWN = "ALT DOWN";
  private static final String ALT_LEFT = "ALT LEFT";
  private static final String ALT_RIGHT = "ALT RIGHT";
  private static final String FIRST = "First";
  private static final String PREVIOUS = "Previous";
  private static final String NEXT = "Next";
  private static final String LAST = "Last";
  private static final String SYMMETRY = "Symmetry";
  private static final String EXCLUDE = "Exclude";
  private static final String CORRECTION = "Correction";
  private static final String CLUSTER_FILTER_BACKWARD = "Backward";
  private static final String CLUSTER_FILTER_FORWARD = "Forward";
  private static final String CLUSTER_FILTER_DELETE = "Delete";
  private static final String CAPTURE = "Screen capture";
  private static final String DUMP = "Dump raw data";
  private static final String TRAVERSE_ALL = "Traverse all";
  private static final String TRAVERSE_ANNOTATED_ONLY = "Traverse annotated only";
  private static final String TRAVERSE_UNANNOTATED_ONLY = "Traverse unannotated only";
  private static final String[] RADIOBUTTON_TEXTS = new String[] {TRAVERSE_ALL,
                                                                  TRAVERSE_ANNOTATED_ONLY,
                                                                  TRAVERSE_UNANNOTATED_ONLY};
  public static final String MASK_MISSING = "Mask missing values";
  public static final String UNMASK_MISSING = "Unmask missing values";
  public static final String MENDELIAN_ERROR = "Mendelian Errors";
  public static final Color BACKGROUND_COLOR = Color.WHITE;
  public static final String DEFAULT_MESSAGE = "enter new annotation here";
  public static final String[] GENOTYPE_OPTIONS = new String[] {"-", "A/A", "A/B", "B/B"};
  public static final int NUM_MARKERS_TO_SAVE_IN_HISTORY = 10;
  public static final String[][] TYPES = { /* { "X Raw", "Y Raw" }, */ {"X", "Y"}, {"Theta", "R"},
                                           {"B Allele Freq", "Log R Ratio"}};
  private static final String NEW_LIST_COMMAND = "New List";
  private static final String SAVE_LIST_COMMAND = "Save List";
  private static final String LOAD_LIST_COMMAND = "Load List";
  private static final String BLAST_DETAILS_COMMAND = "BLAST Details";

  private JButton first, previous, next, last;
  private JTextField navigationField;
  // private ScatterPanel selectedScatterPanel;
  private int selectedPanelIndex = 0;
  private ScatterPanel[] scatterPanels;
  private JLabel sizeLabel;
  private JLabel gcLabel;
  private JPanel qcPanel;
  private JPanel blastPanel;
  private JLabel pcLabel;// , nStageStDevLabel, correctionRatioLabel;
  private JScrollPane annotationScrollPane;
  private JPanel annotationPanel;
  private JPanel annotationPanelLowerPart;
  private JCheckBox[] annotationCheckBoxes;
  private JRadioButton[] filterRadioButtons;
  private JTextField addAnnotationField;
  private boolean showAnnotationShortcuts;
  private boolean annotationAutoAdv;
  private boolean showAllMarkersOrNot;
  private boolean showAnnotatedOrUnannotated;
  private boolean isInitilizing;
  private char[] annotationKeys;


  private JComboBox newGenotype;
  private SpringLayout annotationPanelLowerPartLayout;
  private Project proj;
  private String[] masterMarkerList;
  private String[] masterCommentList;
  private String[] markerList;
  private String[] commentList;
  private float[][][][] centroids;
  private String[] centList;
  private String[][] gcCList;
  private GcAdjustorParameters[] gcAdjustorParameters;
  private JCheckBox[] gcAdjustorBoxes;
  private boolean[] displaygcAdjustor;
  private JLabel[] gcAdjustorLabels;
  private JCheckBox[] centBoxes;
  private boolean[] displayCentroids;
  private JLabel[] centroidLabels;
  private int markerIndex;
  private int markerIndexBak;
  private int[] markerIndexHistory;
  private int previousMarkerIndex;
  private JTextField markerName, commentLabel;
  private String[] samples, sampleFIDIIDs;
  private volatile int[] plot_types;
  private volatile int[] classes;
  private byte size;
  private float gcThreshold;
  private long sampleListFingerprint;
  private MarkerLookup markerLookup;
  private boolean jar;
  private SampleData sampleData;
  private JCheckBox symmetryBox;
  // private JCheckBox excludeMarkerBox;
  private JCheckBox excludeSampleBox;
  private JCheckBox correctionBox;
  private JCheckBox maskMissingBox;
  private JCheckBox mendelianErrorBox;
  private volatile boolean[] correction;
  private volatile boolean[] symmetry;
  private volatile boolean[] maskMissing;
  private volatile boolean[] excludeSamples;
  private volatile boolean[] displayMendelianErrors;
  private ClusterFilterCollection clusterFilterCollection;
  private AnnotationCollection annotationCollection;
  private byte currentClusterFilter;
  private JTextField clusterFilterNavigation;
  private boolean isClusterFilterUpdated;
  private boolean isAnnotationUpdated;
  private String sessionID;
  private AutoSaveForScatterPlot autoSave;
  private JRadioButton[] typeRadioButtons;
  private MarkerDataLoader markerDataLoader;
  private ColorKeyPanel colorKeyPanel;
  // private ColorKeyPanel[] colorKeyPanels;
  // private Color[] colorScheme;
  private int indexOfAnnotationUsedAsMarkerList;
  private boolean fail;
  private boolean exitOnClose;
  private boolean showingAll = false;
  private Logger log;
  private double stdevFilter, correctionRatio;
  private int numComponents;
  private PrincipalComponentsResiduals pcResids;
  private Pedigree pedigree;
  private JPanel[] indivPanels;
  private JPanel scatterOverview;
  private JPanel viewPanel;
  private boolean hasAnnotationFile = false;
  private MarkerAnnotationLoader annotationLoader = null;
  private MarkerBlastAnnotation[] blastResults = null;
  private MarkerGCAnnotation[] gcAnnotations = null;
  private ReferenceGenome referenceGenome = null;

  private BlastParams blastParams = null;
  private Hashtable<String, Integer> markerProjectIndices;
  private final HashMap<String, PlinkMarkerLoader> plinkMarkerLoaders =
                                                                      new HashMap<String, PlinkMarkerLoader>();

  private BlastFrame blastFrame;

  private TwoDPlot histFrame;

  private final ItemListener classListener = new ItemListener() {
    @Override
    public void itemStateChanged(ItemEvent ie) {
      JRadioButton jrb = (JRadioButton) ie.getItem();
      if (sampleData != null && jrb.isSelected() && ie.getStateChange() == ItemEvent.SELECTED) {
        for (byte i = 0; i < sampleData.getNumClasses(); i++) {
          if (jrb.getText().equals(sampleData.getClassName(i))) {
            classes[selectedPanelIndex] = i;
            if (colorKeyPanel != null) {
              colorKeyPanel.setCurrentClass(i);
            }
            scatterPanels[selectedPanelIndex].paintAgain();
            scatterOverview.repaint();
            int plinkIndex = i - sampleData.getNumActualClasses() - sampleData.getNumCNVClasses()
                             - sampleData.getBasicClasses().length;
            if (plinkIndex >= 0) {
              loadPlink(plinkIndex);
            }
          }
        }
      }
    }
  };

  // if the user swaps panels too quickly, color codes get mixed up and start changing
  // adding the 'swapping' check and flag eliminates this bug
  private volatile boolean swapping = false;

  public ScatterPlot(Project project, String[] initMarkerList, String[] initCommentList,
                     boolean exitOnClose) {
    super("Genvisis - ScatterPlot - " + project.PROJECT_NAME.getValue());

    String PROG_KEY = "SCATTERPLOT";
    project.getProgressMonitor().beginDeterminateTask(PROG_KEY, "Displaying ScatterPlot...", 10,
                                                      ProgressMonitor.DISPLAY_MODE.GUI_ONLY);

    JPanel scatterPlotPanel = new JPanel();

    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    setContentPane(scatterPlotPanel);
    addWindowListener(this);

    SampleList sampleList;
    // long time;

    // time = new Date().getTime();

    proj = project;
    jar = proj.JAR_STATUS.getValue();
    log = proj.getLog();
    size = DEFAULT_SIZE;
    this.exitOnClose = exitOnClose;
    gcThreshold = DEFAULT_GC_THRESHOLD / 100f;
    stdevFilter = 0;
    correctionRatio = PrincipalComponentsIntensity.DEFAULT_CORRECTION_RATIO;
    markerIndexHistory = Array.intArray(NUM_MARKERS_TO_SAVE_IN_HISTORY, -1);

    if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue(false, true),
                      proj.JAR_STATUS.getValue())) {
      JOptionPane.showMessageDialog(null,
                                    "Directory " + proj.getProperty(proj.MARKER_DATA_DIRECTORY)
                                          + " does not exist; the raw data needs to be parsed and transposed before it can be visualized",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      fail = true;
      proj.getProgressMonitor().endTask(PROG_KEY);
      return;
    }

    proj.getProgressMonitor().updateTask(PROG_KEY);

    if (Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, true), "marker",
                   MarkerData.MARKER_DATA_FILE_EXTENSION, false,
                   proj.JAR_STATUS.getValue()).length == 0) {
      JOptionPane.showMessageDialog(null,
                                    "There is no data in directory "
                                          + proj.getProperty(proj.MARKER_DATA_DIRECTORY)
                                          + "; the raw data needs to be parsed and transposed before it can be visualized",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      fail = true;
      proj.getProgressMonitor().endTask(PROG_KEY);
      return;
    }

    proj.getProgressMonitor().updateTask(PROG_KEY);

    markerLookup = proj.getMarkerLookup();
    if (markerLookup == null) {
      fail = true;
      proj.getProgressMonitor().endTask(PROG_KEY);
      return;
    }

    proj.getProgressMonitor().updateTask(PROG_KEY);

    sampleList = proj.getSampleList();
    samples = sampleList.getSamples();
    sampleListFingerprint = sampleList.getFingerprint();
    sampleData = proj.getSampleData(3, true);
    convertSamples();

    proj.getProgressMonitor().updateTask(PROG_KEY);

    fail = sampleData.failedToLoad();
    if (fail) {
      proj.getLog().reportError("Without a SampleData file, ScatterPlot will not start");
      JOptionPane.showMessageDialog(null, "Without a SampleData file, ScatterPlot will not start",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      proj.getProgressMonitor().endTask(PROG_KEY);
      return;
    }

    proj.getProgressMonitor().updateTask(PROG_KEY);

    String annoFile = proj.BLAST_ANNOTATION_FILENAME.getValue();
    if (Files.exists(annoFile)) {
      hasAnnotationFile = true;
      blastParams = new BlastParams(log);
      annotationLoader = new MarkerAnnotationLoader(proj, new AnalysisParams[] {blastParams},
                                                    annoFile, proj.getMarkerSet(), true);
    }
    String fastaFile = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
    if (Files.exists(fastaFile)) {
      referenceGenome = new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(),
                                            proj.getLog());
    }

    masterMarkerList = initMarkerList;
    masterCommentList = initCommentList;
    if (masterMarkerList == null) {
      String filename = proj.DISPLAY_MARKERS_FILENAMES.getValue()[0];
      loadMarkerListFromFile(filename);
      if (masterMarkerList == null) {
        String[] tmp = Array.subArray(proj.getMarkerNames(), 0,
                                      Math.min(10, proj.getMarkerNames().length));
        Files.writeList(tmp, proj.DISPLAY_MARKERS_FILENAMES.getValue()[0]);
        loadMarkerListFromFile(filename);
        if (masterMarkerList == null) {
          fail = true;
          proj.getProgressMonitor().endTask(PROG_KEY);
          return;
        } else {
          JOptionPane.showMessageDialog(null,
                                        "Generated a temporary display file at '" + filename + "'",
                                        "Created marker list", JOptionPane.INFORMATION_MESSAGE);
        }
      }
      if (masterMarkerList.length == 0) {
        JOptionPane.showMessageDialog(null,
                                      "Error - file '" + filename
                                            + "' was devoid of any valid markers",
                                      "Error", JOptionPane.ERROR_MESSAGE);
        fail = true;
        proj.getProgressMonitor().endTask(PROG_KEY);
        return;
      }
    }
    if (masterCommentList == null) {
      masterCommentList = Array.stringArray(masterMarkerList.length, "");
    }


    proj.getProgressMonitor().updateTask(PROG_KEY);

    resetAfterLoad();

    proj.getProgressMonitor().updateTask(PROG_KEY);

    // Java initializes boolean arrays as false
    maskMissing = new boolean[4];
    symmetry = new boolean[] {true, true, true, true};
    excludeSamples = new boolean[] {true, true, true, true};
    correction = new boolean[4];
    displayMendelianErrors = new boolean[4];
    plot_types = new int[] {0, 0, 0, 0};
    classes = new int[] {2, 2, 2, 2};
    scatterPanels = new ScatterPanel[4];
    for (int i = 0; i < scatterPanels.length; i++) {
      scatterPanels[i] = new ScatterPanel(this, i);
    }
    loadAnnotationCollection();

    setLayout(new BorderLayout());

    setJMenuBar(createJMenuBar());

    proj.getProgressMonitor().updateTask(PROG_KEY);

    scatterOverview = new JPanel(new GridLayout(0, 2));
    indivPanels = new JPanel[4];
    for (int i = 0; i < scatterPanels.length; i++) {
      final int myInd = i;
      indivPanels[i] = new JPanel(new BorderLayout());
      indivPanels[i].setBackground(Color.WHITE);
      MouseAdapter clickAdapter = new MouseAdapter() {
        @Override
        public void mousePressed(MouseEvent e) {
          super.mousePressed(e);
          doPanelSelection(myInd, false);
        }

        @Override
        public void mouseReleased(MouseEvent e) {
          super.mouseReleased(e);
          SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
              updateGUI();
            }
          });
        }
      };
      indivPanels[i].addMouseListener(clickAdapter);
      scatterPanels[i].addMouseListener(clickAdapter);
      indivPanels[i].add(scatterPanels[i], BorderLayout.CENTER);
      // indivPanels[i].setDoubleBuffered(true);
      // scatterPanels[i].setDoubleBuffered(true);
      scatterPanels[i].setSymmetricAxes(true);
      scatterOverview.add(indivPanels[i]);

      proj.getProgressMonitor().updateTask(PROG_KEY);
    }

    colorKeyPanel = new ColorKeyPanel(sampleData, scatterPanels[0], scatterPanels[0].colorScheme,
                                      classListener, 2);

    // scatterOverview.setDoubleBuffered(true);
    viewPanel = new JPanel(new BorderLayout());
    viewPanel.add(indivPanels[selectedPanelIndex], BorderLayout.CENTER);
    viewPanel.add(markerPanel(), BorderLayout.NORTH);
    viewPanel.add(colorKeyPanel, BorderLayout.SOUTH);

    getContentPane().add(viewPanel, BorderLayout.CENTER);
    JScrollPane eastScroll = new JScrollPane();
    JPanel eastPanel = eastPanel();
    eastScroll.setBorder(null);
    eastScroll.setViewportView(eastPanel);
    eastScroll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
    eastScroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
    getContentPane().add(eastScroll, BorderLayout.EAST);

    for (ScatterPanel scatterPanel : scatterPanels) {
      scatterPanel.setPointsGeneratable(true);
      scatterPanel.setUpdateQCPanel(true);
      scatterPanel.setExtraLayersVisible(new byte[] {99});
    }
    displayIndex(navigationField);
    // clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+"
    // of "+clusterFilterCollection.getSize(getMarkerName()));
    // currentClusterFilter=0;
    // setCurrentClusterFilter((byte) (clusterFilterCollection.getSize(getMarkerName())-1));
    // setCurrentClusterFilter();
    if (clusterFilterCollection == null) {
      currentClusterFilter = 0;
    } else {
      currentClusterFilter = (byte) (clusterFilterCollection.getSize(getMarkerName()) - 1);
    }
    // scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
    displayClusterFilterIndex();
    activateAllAnnotationMaps();

    // newGenotype.setSelectedIndex(clusterFilterCollection.getGenotype(getMarkerName(),
    // currentClusterFilter)+1);
    symmetryBox.setSelected(true);
    // excludeBox.setSelected(false);
    // currentStateOfExcludeBox = false;
    if (centList.length > 0) {
      centBoxes[0].setSelected(true);
    }

    inputMapAndActionMap(scatterPlotPanel);
    // inputMapAndActionMap(annotationPanel);
    // inputMapAndActionMap(annotationScrollPane);

    // next.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT);
    // next.setActionMap(actionMap);
    // previous.setActionMap(actionMap);
    scatterPanels[selectedPanelIndex].grabFocus();
    markerProjectIndices = proj.getMarkerIndices();
    updateBLASTPanel();
    updateGUI();

    proj.getProgressMonitor().endTask(PROG_KEY);
  }

  @Override
  public void actionPerformed(ActionEvent ae) {
    String command = ae.getActionCommand();
    if (command.equals(FIRST)) {
      first();
    } else if (command.equals(PREVIOUS)) {
      previous();
    } else if (command.equals(NEXT)) {
      next();
    } else if (command.equals(LAST)) {
      last();
    } else if (command.equals(CLUSTER_FILTER_BACKWARD)) {
      if (clusterFilterCollection.getSize(getMarkerName()) > 0) {
        // seletedScatterPanel.rectangles[currentClusterFilter].setColor((byte)7);
        for (ScatterPanel scatterPanel : scatterPanels) {
          if (scatterPanel.rectangles != null) {
            scatterPanel.rectangles[currentClusterFilter].setColor((byte) 7);
          }
        }
        // scatPanel.setRectangles();
        setCurrentClusterFilter((byte) Math.max(currentClusterFilter - 1, 0));
        // scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
        // clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+"
        // of "+clusterFilterCollection.getSize(getMarkerName()));
        displayClusterFilterIndex();
        // annotationPanelLowerPart();
        // seletedScatterPanel.setPointsGeneratable(true); // why not
        // scatPanel.setPointsGeneratable(false);
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true); // why not scatPanel.setPointsGeneratable(false);
        }
        // scatPanel.setUpdateQcPanel(true);
        updateGUI();
      }
    } else if (command.equals(CLUSTER_FILTER_FORWARD)) {
      if (clusterFilterCollection.getSize(getMarkerName()) > 0) {
        // seletedScatterPanel.rectangles[currentClusterFilter].setColor((byte)7);
        for (ScatterPanel scatterPanel : scatterPanels) {
          if (scatterPanel.rectangles != null) {
            scatterPanel.rectangles[currentClusterFilter].setColor((byte) 7);
          }
        }
        // scatPanel.setRectangles();
        setCurrentClusterFilter((byte) Math.min(currentClusterFilter + 1,
                                                (clusterFilterCollection.getSize(getMarkerName()) == 0 ? 0
                                                                                                       : (clusterFilterCollection.getSize(getMarkerName())
                                                                                                          - 1))));
        // scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
        // clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+"
        // of "+clusterFilterCollection.getSize(getMarkerName()));
        displayClusterFilterIndex();
        // annotationPanelLowerPart();
        // seletedScatterPanel.setPointsGeneratable(true); // why not
        // scatPanel.setPointsGeneratable(false);
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true); // why not scatPanel.setPointsGeneratable(false);
        }
        // scatPanel.setUpdateQcPanel(true);
        updateGUI();
      }
    } else if (command.equals(CLUSTER_FILTER_DELETE)) {
      if (clusterFilterCollection.getSize(getMarkerName()) > 0) {
        // scatPanel.clearRectangles(currentClusterFilter);
        clusterFilterCollection.deleteClusterFilter(getMarkerName(), currentClusterFilter);
        setCurrentClusterFilter((byte) Math.min(currentClusterFilter,
                                                clusterFilterCollection.getSize(getMarkerName())
                                                                      - 1));
        // startAutoSave();
        setClusterFilterUpdated(true);
        // clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+"
        // of "+clusterFilterCollection.getSize(getMarkerName()));
        displayClusterFilterIndex();
        // annotationPanelLowerPart();
        // seletedScatterPanel.setPointsGeneratable(true);
        // seletedScatterPanel.setQcPanelUpdatable(true);
        // seletedScatterPanel.generateRectangles();
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
          scatterPanel.setUpdateQCPanel(true);
          scatterPanel.generateRectangles();
        }
        // if (clusterFilterCollection.getSize(getMarkerName())>0) {
        // scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
        // }
        updateGUI();
      }
    } else if (command.equals(CAPTURE)) {
      doScreenCapture();
    } else if (command.equals(DUMP)) {
      doDataDump();
    } else if (command.equals(MASK_MISSING) || command.equals(UNMASK_MISSING)) {
      maskMissing[selectedPanelIndex] = !maskMissing[selectedPanelIndex];
      // ((AbstractButton)ae.getSource()).setText(maskMissing?UNMASK_MISSING:MASK_MISSING);
      // seletedScatterPanel.setPointsGeneratable(true);
      // seletedScatterPanel.setQcPanelUpdatable(true);
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.setPointsGeneratable(true);
        scatterPanel.setUpdateQCPanel(true);
      }
      updateGUI();
    } else if (command.equals(SYMMETRY)) {
      symmetry[selectedPanelIndex] = !symmetry[selectedPanelIndex];
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.setPointsGeneratable(true);
        scatterPanel.setUpdateQCPanel(true);
      }
      updateGUI();
    } else if (command.equals(CORRECTION)) {
      correction[selectedPanelIndex] = !correction[selectedPanelIndex];
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.setPointsGeneratable(true);
        scatterPanel.setUpdateQCPanel(true);
      }
      updateGUI();
    } else if (command.equals(MENDELIAN_ERROR)) {
      displayMendelianErrors[selectedPanelIndex] = mendelianErrorBox.isSelected();
      updateGUI();
    } else if (command.equals(NEW_LIST_COMMAND)) {
      NewMarkerListDialog newMkrList =
                                     new NewMarkerListDialog(proj.getMarkerNames(),
                                                             proj.getProperty(proj.PROJECT_DIRECTORY));
      newMkrList.setModal(true);
      newMkrList.setVisible(true);
      if (newMkrList.getReturnCode() == JOptionPane.YES_OPTION) {
        String mkrFile = newMkrList.getFileName();
        proj.DISPLAY_MARKERS_FILENAMES.addValue(mkrFile);
        loadMarkerFile(mkrFile);
      }
      // } else if (command.equals(LOAD_LIST_COMMAND + "test")) {
      //
      // final Configurator configurator = new Configurator(proj, new
      // String[][]{Configurator.ALL_PROPERTY_SETS[10]});
      // configurator.setVisible(true);

    } else if (command.equals(LOAD_LIST_COMMAND)) {
      JFileChooser jfc = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
      jfc.setMultiSelectionEnabled(false);
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogType(JFileChooser.OPEN_DIALOG);
      jfc.setDialogTitle("Load marker list file");
      int code = jfc.showDialog(this, "Load File");
      if (code == JFileChooser.APPROVE_OPTION) {
        try {
          String file = jfc.getSelectedFile().getCanonicalPath();
          file = ext.replaceAllWith(file, "\\", "/");
          proj.DISPLAY_MARKERS_FILENAMES.addValue(file);
          loadMarkerFile(file);
        } catch (IOException e) {
          proj.message("Error - invalid file selected");
        }
      }
    } else if (command.equals(SAVE_LIST_COMMAND)) {

      JFileChooser jfc = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
      jfc.setMultiSelectionEnabled(false);
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save Marker List");
      jfc.setDialogType(JFileChooser.SAVE_DIALOG);
      int code = jfc.showSaveDialog(this);
      if (code == JFileChooser.APPROVE_OPTION) {
        Files.writeList(markerList, jfc.getSelectedFile().getAbsolutePath());
      }
    } else if (command.equals(BLAST_DETAILS_COMMAND)) {
      displayBlast();
    } else {
      log.reportError("Error - unknown command '" + command + "'");
    }
  }

  public JPanel annotationPanelLower() {
    JButton removeButton;
    ButtonGroup radioButtonGroup;
    int horizontalMargin;
    int componentHeight;
    int currentHorizontalPos;
    int maxWidth;

    annotationPanelLowerPart.removeAll();
    // annotationPanelLowerPart.repaint();

    horizontalMargin = 2;
    componentHeight = 20;
    currentHorizontalPos = 0;
    radioButtonGroup = new ButtonGroup();
    filterRadioButtons = new JRadioButton[RADIOBUTTON_TEXTS.length];
    for (int i = 0; i < RADIOBUTTON_TEXTS.length; i++) {
      filterRadioButtons[i] = new JRadioButton(RADIOBUTTON_TEXTS[i] + " (n=placeholder)");
      filterRadioButtons[i].setBackground(Color.WHITE);
      filterRadioButtons[i].addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          String commandText;
          commandText = e.getActionCommand();
          if (commandText.startsWith(TRAVERSE_ALL + " (n=")) {
            showAllMarkersOrNot = true;
          } else if (commandText.startsWith(TRAVERSE_ANNOTATED_ONLY + " (n=")) {
            showAllMarkersOrNot = false;
            showAnnotatedOrUnannotated = true;
          } else {
            showAllMarkersOrNot = false;
            showAnnotatedOrUnannotated = false;
          }
        }
      });
      radioButtonGroup.add(filterRadioButtons[i]);
      annotationPanelLowerPart.add(filterRadioButtons[i]);
      annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, filterRadioButtons[i], 5,
                                                   SpringLayout.WEST, annotationPanelLowerPart);
      annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, filterRadioButtons[i],
                                                   currentHorizontalPos, SpringLayout.NORTH,
                                                   annotationPanelLowerPart);
      // horizontalPos += (horizontalMargin + fileterRadioButtons[i].getSize().height);
      currentHorizontalPos += (horizontalMargin + componentHeight);
    }
    filterRadioButtons[0].setSelected(true);

    MouseListener mouseListenerForAnnotationCheckBoxes = new MouseListener() {
      @Override
      public void mouseClicked(MouseEvent e) {
        JCheckBox source;
        JPopupMenu menu;
        String annotation;
        int annotationIndex;

        source = (JCheckBox) e.getSource();
        annotation = source.getText();

        annotationIndex = -1;
        for (int i = 0; i < annotationKeys.length; i++) {
          if (annotationCollection.getDescriptionForComment(annotationKeys[i],
                                                            showAnnotationShortcuts, true)
                                  .equals(annotation)) {
            annotationIndex = i;
            break;
          }
        }

        if (annotation.charAt(0) == '\'' && annotation.charAt(2) == '\'') {
          annotation = annotation.substring(4);
        }
        annotation = annotation.substring(0, annotation.lastIndexOf(" (n="));

        if (e.getButton() == MouseEvent.BUTTON3) {

          menu = new JPopupMenu();
          menu.setName(annotationIndex + "");

          menu.add(new AbstractAction("Display all markers annotated with " + annotation) {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e1) {
              String annotation;
              int annotationIndex;

              annotationIndex = -1;
              annotation = e1.getActionCommand();
              annotation = annotation.substring("Display all markers annotated with ".length());
              for (int i = 0; i < annotationKeys.length; i++) {
                if (annotationCollection.getDescriptionForComment(annotationKeys[i], false, false)
                                        .equals(annotation)) {
                  annotationIndex = i;
                  break;
                }
              }

              changeToNewMarkerList(annotationCollection.getMarkerLists(annotationKeys[annotationIndex]),
                                    markerIndex);
            }

            @Override
            public boolean isEnabled() {
              return indexOfAnnotationUsedAsMarkerList != -2;
            }
          });

          menu.add(new AbstractAction("Revert marker list to original list") {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e2) {
              revertMarkersToOriginalList(markerIndexBak);
            }

            @Override
            public boolean isEnabled() {
              return indexOfAnnotationUsedAsMarkerList != -1;
            }
          });

          menu.add(new AbstractAction((indexOfAnnotationUsedAsMarkerList == -1 ? "Limit"
                                                                               : "Unlimit")
                                      + " Current List to Those with Annotation " + annotation) {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e3) {
              String annotation;
              int annotationIndex;

              annotationIndex = -1;
              annotation = e3.getActionCommand();
              annotation =
                         annotation.substring(((indexOfAnnotationUsedAsMarkerList == -1 ? "Limit"
                                                                                        : "Unlimit")
                                               + " Current List to Those with Annotation ").length());
              for (int i = 0; i < annotationKeys.length; i++) {
                if (annotationCollection.getDescriptionForComment(annotationKeys[i], false, false)
                                        .toLowerCase().equals(annotation)) {
                  annotationIndex = i;
                  break;
                }
              }

              if (indexOfAnnotationUsedAsMarkerList == -1) {
                for (int i = 0; i < annotationIndex; i++) {
                  annotationCheckBoxes[i].setEnabled(false);
                }
                for (int i = annotationIndex + 1; i < annotationCheckBoxes.length; i++) {
                  annotationCheckBoxes[i].setEnabled(false);
                }
                indexOfAnnotationUsedAsMarkerList = annotationIndex;
              } else if (indexOfAnnotationUsedAsMarkerList == annotationIndex) {
                for (JCheckBox annotationCheckBoxe : annotationCheckBoxes) {
                  annotationCheckBoxe.setEnabled(true);
                }
                indexOfAnnotationUsedAsMarkerList = -1;
              }
            }
          });

          menu.show(source, e.getX(), e.getY());
        }
      }

      @Override
      public void mouseEntered(MouseEvent e) {}

      @Override
      public void mouseExited(MouseEvent e) {}

      @Override
      public void mousePressed(MouseEvent e) {}

      // public void mouseClicked(MouseEvent e) {
      // JCheckBox source;
      // String annotation;
      // int annotationIndex;
      //
      // if (e.getButton()==MouseEvent.BUTTON3) {
      // annotationIndex = -1;
      // source = (JCheckBox)e.getSource();
      // annotation = source.getText();
      // for (int i = 0; i < annotationKeys.length; i ++) {
      // if (annotationCollection.getDescriptionForComment(annotationKeys[i],
      // showAnnotationShortcuts, true).equals(annotation)) {
      // annotationIndex = i;
      // break;
      // }
      // }
      //
      // if (indexOfAnnotationUsedAsMarkerList == -1) {
      // for (int i = 0; i < annotationIndex; i ++) {
      // annotationCheckBoxes[i].setEnabled(false);
      // }
      // for (int i = annotationIndex + 1; i < annotationCheckBoxes.length; i ++) {
      // annotationCheckBoxes[i].setEnabled(false);
      // }
      // indexOfAnnotationUsedAsMarkerList = annotationIndex;
      // } else if (indexOfAnnotationUsedAsMarkerList == annotationIndex) {
      // for (int i = 0; i < annotationCheckBoxes.length; i ++) {
      // annotationCheckBoxes[i].setEnabled(true);
      // }
      // indexOfAnnotationUsedAsMarkerList = -1;
      // }
      // }

      @Override
      public void mouseReleased(MouseEvent e) {}
    };

    maxWidth = 200;
    annotationCheckBoxes = new JCheckBox[annotationKeys.length];
    for (int i = 0; i < annotationKeys.length; i++) {
      annotationCheckBoxes[i] =
                              new JCheckBox(annotationCollection.getDescriptionForComment(annotationKeys[i],
                                                                                          showAnnotationShortcuts,
                                                                                          true));
      maxWidth = Math.max(annotationCheckBoxes[i].getPreferredSize().width, maxWidth);
      annotationCheckBoxes[i].setBackground(Color.WHITE);
      if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
        annotationCheckBoxes[i].setSelected(true);
      }
      annotationCheckBoxes[i].addItemListener(new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent itemEvent) {
          if (!isInitilizing) {
            JCheckBox checkBox;
            char currentKey;

            checkBox = (JCheckBox) itemEvent.getSource();
            currentKey = checkBox.getText().charAt(1);

            if (checkBox.getModel().isSelected()) {
              annotationCollection.addAnnotationForMarker(markerList[markerIndex], currentKey);
            } else {
              annotationCollection.removeAnnotationForMarker(markerList[markerIndex], currentKey);
            }
            checkBox.setText(annotationCollection.getDescriptionForComment(currentKey,
                                                                           showAnnotationShortcuts,
                                                                           true));
            updateAnnotationPanelFilterRadioButtons();
            // isAnnotationUpdated = true;
            setAnnotationUpdated(true);
          }
        }
      });
      annotationCheckBoxes[i].addMouseListener(mouseListenerForAnnotationCheckBoxes);

      annotationPanelLowerPart.add(annotationCheckBoxes[i]);
      annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, annotationCheckBoxes[i], 5,
                                                   SpringLayout.WEST, annotationPanelLowerPart);
      annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, annotationCheckBoxes[i],
                                                   currentHorizontalPos, SpringLayout.NORTH,
                                                   annotationPanelLowerPart);
      // horizontalPos += (horizontalMargin + annotationCheckBoxes[i].getSize().height);
      currentHorizontalPos += (horizontalMargin + componentHeight);
      removeButton = new JButton(Grafik.getImageIcon("images/delete2sm.png"));
      removeButton.setActionCommand("removeButton" + i);
      removeButton.setBorder(null);
      removeButton.setSize(6, 6);
      removeButton.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          String commandText;
          int index;

          commandText = e.getActionCommand();
          index = Integer.parseInt(commandText.substring("removeButton".length(),
                                                         commandText.length()));
          annotationCollection.removeAnnotation(proj, annotationKeys[index]);
          removeAnnotationFromMaps(annotationKeys[index]);
          annotationKeys = annotationCollection.getKeys();
          // isAnnotationUpdated = true;
          setAnnotationUpdated(true);
          annotationPanelLower();
        }
      });
      annotationPanelLowerPart.add(removeButton);
      annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, removeButton, 5,
                                                   SpringLayout.EAST, annotationCheckBoxes[i]);
      annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, removeButton, 7,
                                                   SpringLayout.NORTH, annotationCheckBoxes[i]);
    }

    addAnnotationField = new JTextField(DEFAULT_MESSAGE);

    addAnnotationField.addFocusListener(new FocusListener() {
      @Override
      public void focusGained(FocusEvent e) {
        if (addAnnotationField.getText().equals(DEFAULT_MESSAGE)) {
          addAnnotationField.setText("");
          Font font = addAnnotationField.getFont();
          addAnnotationField.setFont(new Font(font.getFontName(), Font.PLAIN, font.getSize()));
        }
        deactivateAllAnnotationMaps();
      }

      @Override
      public void focusLost(FocusEvent e) {
        if (showAnnotationShortcuts) {
          activateAllAnnotationMaps();
        }
        if (addAnnotationField.getText().equals("")) {
          addAnnotationField.setText(DEFAULT_MESSAGE);
          Font font = addAnnotationField.getFont();
          addAnnotationField.setFont(new Font(font.getFontName(), Font.ITALIC, font.getSize()));
        }

      }
    });

    addAnnotationField.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        String str;
        char c;

        str = ((JTextField) e.getSource()).getText();
        // c = str.toLowerCase().charAt(0);
        // str = str.substring(1);
        c = AnnotationCollection.assignKey(str, annotationCollection);
        if (c == 0) {
          log.report("Failed to add the new annotation key '" + str
                     + "', because cannot automatically assign a shortcut for it.");
        } else {
          if (annotationCollection == null) {
            annotationCollection = new AnnotationCollection();
            startAutoSaveToTempFile();
          }
          annotationCollection.addAnnotation(c, str);
          addAnnotationToMaps(c);
          annotationKeys = annotationCollection.getKeys();
          annotationPanelLower();
          // isAnnotationUpdated = true;
          setAnnotationUpdated(true);
        }
      }
    });

    annotationPanelLowerPart.add(addAnnotationField);
    annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, addAnnotationField, 5,
                                                 SpringLayout.WEST, annotationPanelLowerPart);
    annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, addAnnotationField,
                                                 currentHorizontalPos + 3, SpringLayout.NORTH,
                                                 annotationPanelLowerPart);
    // horizontalPos += (horizontalMargin + addAnnotationField.getSize().height);
    // horizontalPos += (horizontalMargin + componentHeight);

    // TODO This could be fixed such that annotationPanel is always a fixed width and we have an
    // inner panel with horizontal scroll bar, but this has proved too complex for the moment
    // annotationPanelLowerPartLayout.getConstraints(annotationPanelLowerPart).setConstraint(SpringLayout.EAST,
    // Spring.constant(250, 250, 250));
    // annotationPanel.setPreferredSize(new Dimension(200, 500));
    annotationPanelLowerPart.setPreferredSize(new Dimension(maxWidth + 42,
                                                            ((annotationKeys == null ? 0
                                                                                     : annotationKeys.length)
                                                             + 4) * 22));

    // annotationPanelLowerPart.invalidate();
    annotationPanel.validate();
    if (annotationScrollPane != null) {
      annotationScrollPane.validate();
      annotationScrollPane.repaint();
    }

    return annotationPanelLowerPart;
  }

  public void changeToNewMarkerList(String[] newMarkerList, int indexOfDefaultFirstMarker) {
    if (!saveClusterFilterAndAnnotationCollection()) {
      return;
    }

    indexOfAnnotationUsedAsMarkerList = -2;
    markerIndexBak = markerIndex;
    // markerList = annotationCollection.getMarkerLists(annotationKeys[annotationIndex]);
    markerList = newMarkerList;
    commentList = new String[markerList.length];
    loadMarkerDataFromList(0);
    displayIndex(navigationField);
    updateGUI();
  }

  public void checkAnnotationBox(char c) {
    int index;

    index = ext.indexOfChar(c, annotationKeys);
    if (!annotationCheckBoxes[index].isSelected()) {
      annotationCheckBoxes[index].setSelected(true);
      if (annotationAutoAdv) {
        actionPerformed(new ActionEvent(next, 0, NEXT));
      }
    }
  }

  public void displayClusterFilterIndex() {
    clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName()) == 0 ? 0
                                                                                           : (currentClusterFilter
                                                                                              + 1))
                                    + " of " + clusterFilterCollection.getSize(getMarkerName()));
  }

  public void displayIndex(JTextField field) {
    field.setText((markerIndex + 1) + " of " + markerList.length);
  }

  public boolean displayMendelianErrors(int index) {
    return displayMendelianErrors[index];
  }

  public boolean failed() {
    return fail;
  }

  public void finishProcessing() {
    updateMarkerIndexHistory();
    displayIndex(navigationField);
    if (showingAll) {
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.setPointsGeneratable(true);
        scatterPanel.setUpdateQCPanel(true);
      }
    } else {
      scatterPanels[selectedPanelIndex].setPointsGeneratable(true);
      scatterPanels[selectedPanelIndex].setUpdateQCPanel(true);
    }
    setCurrentClusterFilter();
    updateBLASTPanel();
    updateGUI();
    if (blastFrame != null && blastFrame.isVisible()) {
      blastFrame.setAnnotations(blastResults[markerIndex], referenceGenome);
      if (histFrame != null/* && histFrame.isVisible() */) {
        double filter = proj.BLAST_PROPORTION_MATCH_FILTER.getValue();
        int probe = proj.ARRAY_TYPE.getValue().getProbeLength();
        setHistogram((int) (filter * probe));
      }
    }
    displayClusterFilterIndex();
  }

  public void first() {
    markerIndex = getAvailableMarker(false, false, showAllMarkersOrNot, showAnnotatedOrUnannotated);
    finishProcessing();
    // SwingUtilities.invokeLater(new Runnable() {
    // @Override
    // public void run() {
    // first.setEnabled(false);
    // previous.setEnabled(false);
    // last.setEnabled(markerIndex < markerList.length - 1);
    // next.setEnabled(markerIndex < markerList.length - 1);
    // }
    // });
  }



  // private JPanel nstageStdevSliderPanel() {
  // JPanel pcSliderPanel = new JPanel();
  //// pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
  // pcSliderPanel.setLayout(new GridLayout(2, 1));
  // pcSliderPanel.setBackground(BACKGROUND_COLOR);
  //
  // JSlider slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
  // // slider.setSize(new Dimension(150, 20));
  // slider.setBackground(BACKGROUND_COLOR);
  // slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 0);
  // slider.setValue(0);
  // slider.setBackground(BACKGROUND_COLOR);
  // nStageStDevLabel = new JLabel("Standard deviation filter "+stdevFilter, JLabel.CENTER);
  // nStageStDevLabel.setFont(new Font("Arial", Font.PLAIN, 16));
  // // tabPanel.add(gcLabel, gbc);
  // pcSliderPanel.add(nStageStDevLabel);
  //
  // slider.addChangeListener(new ChangeListener() {
  // public void stateChanged(ChangeEvent ce) {
  // JSlider slider = (JSlider) ce.getSource();
  // stdevFilter = (double)slider.getValue()/25;
  // nStageStDevLabel.setText("Alt. Genotype StDev Filter: "+stdevFilter);
  //// seletedScatterPanel.setPointsGeneratable(true);
  //// seletedScatterPanel.setQcPanelUpdatable(true);
  //// seletedScatterPanel.paintAgain();
  // for (int i = 0; i < scatterPanels.length; i++) {
  // scatterPanels[i].setPointsGeneratable(true);
  // scatterPanels[i].setQcPanelUpdatable(true);
  // scatterPanels[i].paintAgain();
  // }
  // // qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
  // }
  // });
  // // tabPanel.add(slider, gbc);
  // pcSliderPanel.add(slider);
  //
  // return pcSliderPanel;
  // }

  // private JPanel nstageCorrectionRatio() {
  // JPanel pcSliderPanel = new JPanel();
  //// pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
  // pcSliderPanel.setLayout(new GridLayout(2, 1));
  // pcSliderPanel.setBackground(BACKGROUND_COLOR);
  //
  // JSlider slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
  // // slider.setSize(new Dimension(150, 20));
  // slider.setBackground(BACKGROUND_COLOR);
  // slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 0);
  // slider.setValue(0);
  // slider.setBackground(BACKGROUND_COLOR);
  // correctionRatioLabel = new JLabel("Correction Ratio: " + correctionRatio, JLabel.CENTER);
  // correctionRatioLabel.setFont(new Font("Arial", Font.PLAIN, 16));
  // String usage = "This limits the ratio of PCs to the number of samples in a cluster...";
  //
  // correctionRatioLabel.setToolTipText(usage);
  // // tabPanel.add(gcLabel, gbc);
  // pcSliderPanel.add(correctionRatioLabel);
  //
  // slider.addChangeListener(new ChangeListener() {
  // public void stateChanged(ChangeEvent ce) {
  // JSlider slider = (JSlider) ce.getSource();
  // correctionRatio = (double) slider.getValue() / 100;
  // correctionRatioLabel.setText("Correction Ratio " + correctionRatio);
  //// seletedScatterPanel.setPointsGeneratable(true);
  //// seletedScatterPanel.setQcPanelUpdatable(true);
  //// seletedScatterPanel.paintAgain();
  // for (int i = 0; i < scatterPanels.length; i++) {
  // scatterPanels[i].setPointsGeneratable(true);
  // scatterPanels[i].setQcPanelUpdatable(true);
  // scatterPanels[i].paintAgain();
  // }
  // // qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
  // }
  // });
  // // tabPanel.add(slider, gbc);
  // pcSliderPanel.add(slider);
  //
  // return pcSliderPanel;
  // }

  public float[][][][] getCentroids() {
    return centroids;
  }

  public ClusterFilterCollection getClusterFilterCollection() {
    return clusterFilterCollection;
  }

  public boolean getCorrection(int index) {
    return correction[index];
  }

  public double getCorrectionRatio() {
    return correctionRatio;
  }

  public int getCurrentClass(int index) {
    return classes[index];
  }


  public byte getCurrentClusterFilter() {
    while (currentClusterFilter >= clusterFilterCollection.getSize(getMarkerName())) {
      currentClusterFilter--;
    }
    return currentClusterFilter;
  }

  public MarkerData getCurrentMarkerData() {
    MarkerData markerData;
    int count;

    count = 0;
    while ((markerData = markerDataLoader.getMarkerData(markerIndex)) == null) {
      markerDataLoader.requestIndexBeTheNextFileToLoad(markerIndex);
      try {
        Thread.sleep(250);
        count++;
      } catch (InterruptedException ie) {
      }
      if (count > 8 && count % 8 == 0) {
        log.reportError("Error - ScatterPlot has been waiting on markerDataLoader to load "
                        + markerList[markerIndex] + " for " + (count / 4) + " secounds");
      }
    }

    return markerData;
  }

  public Hashtable<String, String> getDisabledClassValues() {
    return colorKeyPanel.getDisabledClassValues();
  }


  // private int getAvailableMarker(boolean forwardOrBackward, boolean firstOrLastInThatDirection,
  // boolean allMarkersOrOnlyThoseAnnotatedOrUnannotated, boolean annotatedOrUnannotated) {
  // int result;
  //
  // result = markerIndex;
  // if (allMarkersOrOnlyThoseAnnotatedOrUnannotated && indexOfAnnotationControllingMarkerList ==
  // -1) {
  // if (firstOrLastInThatDirection) {
  // if (forwardOrBackward) {
  // result = Math.min(markerIndex + 1, markerList.length - 1);
  // } else {
  // result = Math.max(markerIndex - 1, 0);
  // }
  // } else {
  // if (forwardOrBackward) {
  // result = markerList.length - 1;
  // } else {
  // result = 0;
  // }
  // }
  // } else {
  // if (firstOrLastInThatDirection) {
  // if (forwardOrBackward) {
  // if (indexOfAnnotationControllingMarkerList == -1) {
  // if (annotatedOrUnannotated) {
  // for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
  // if (isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // } else {
  // for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
  // if (! isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // }
  // } else {
  // for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
  // if (annotationCollection.markerHasAnnotation(markerList[i],
  // annotationKeys[indexOfAnnotationControllingMarkerList])) {
  // result = i;
  // break;
  // }
  // }
  // }
  // } else {
  // if (indexOfAnnotationControllingMarkerList == -1) {
  // if (annotatedOrUnannotated) {
  // for (int i = markerIndex - 1; i >= 0; i--) {
  // if (isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // } else {
  // for (int i = markerIndex - 1; i >= 0; i--) {
  // if (! isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // }
  // } else {
  // for (int i = markerIndex - 1; i >= 0; i--) {
  // if (annotationCollection.markerHasAnnotation(markerList[i],
  // annotationKeys[indexOfAnnotationControllingMarkerList])) {
  // result = i;
  // break;
  // }
  // }
  // }
  // }
  // } else {
  // if (forwardOrBackward) {
  // if (indexOfAnnotationControllingMarkerList == -1) {
  // if (annotatedOrUnannotated) {
  // for (int i = isAnnotated.length - 1; i >= 0; i--) {
  // if (isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // } else {
  // for (int i = isAnnotated.length - 1; i >= 0; i--) {
  // if (! isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // }
  // } else {
  // for (int i = isAnnotated.length - 1; i >= 0; i--) {
  // if (annotationCollection.markerHasAnnotation(markerList[i],
  // annotationKeys[indexOfAnnotationControllingMarkerList])) {
  // result = i;
  // break;
  // }
  // }
  // }
  // } else {
  // if (indexOfAnnotationControllingMarkerList == -1) {
  // if (annotatedOrUnannotated) {
  // for (int i = 0; i < isAnnotated.length; i++) {
  // if (isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // } else {
  // for (int i = 0; i < isAnnotated.length; i++) {
  // if (! isAnnotated[i]) {
  // result = i;
  // break;
  // }
  // }
  // }
  // } else {
  // for (int i = 0; i < isAnnotated.length; i++) {
  // if (annotationCollection.markerHasAnnotation(markerList[i],
  // annotationKeys[indexOfAnnotationControllingMarkerList])) {
  // result = i;
  // break;
  // }
  // }
  // }
  // }
  // }
  // }
  //
  // return result;
  // }

  public boolean[] getDisplayCentroids() {
    return displayCentroids;
  }

  public boolean[] getDisplaygcAdjustor() {
    return displaygcAdjustor;
  }

  public GcAdjustorParameters[] getGcAdjustorParameters() {
    return gcAdjustorParameters;
  }

  public String[][] getGcCList() {
    return gcCList;
  }

  public float getGCthreshold() {
    return gcThreshold;
  }

  // public MarkerData[] getMarkerData1() {
  // return markerData;
  // }
  //
  public int getMarkerIndex() {
    return markerIndex;
  }

  public String getMarkerName() {
    return markerList[markerIndex];
  }

  public Hashtable<String, Integer> getMarkerProjectIndices() {
    return markerProjectIndices;
  }

  // public void setCurrentClass(byte newCurrentClass) {
  // currentClass = newCurrentClass;
  // }

  public int getNumComponents() {
    return numComponents;
  }

  public PrincipalComponentsResiduals getPcResids() {
    return pcResids;
  }

  public Pedigree getPedigree() {
    return pedigree;
  }

  public byte getPlinkGenotypeForIndi(String sampleID, int currentClass) {
    int plinkIndex = currentClass - sampleData.getBasicClasses().length
                     - sampleData.getNumActualClasses() - sampleData.getNumCNVClasses();
    String plinkRoot = proj.PLINK_DIR_FILEROOTS.getValue()[plinkIndex];
    String[] lookup = sampleData.lookup(sampleID);
    String fidiid = lookup[1];
    String marker = markerList[markerIndex];
    return plinkMarkerLoaders.get(plinkRoot).getGenotypeForIndi(marker, fidiid);
  }

  public int getPlotType(int index) {
    return plot_types[index];
  }


  public byte getPointSize() {
    return size;
  }

  public Project getProject() {
    return proj;
  }

  public SampleData getSampleData() {
    return sampleData;
  }

  public long getSampleFingerprint() {
    return sampleListFingerprint;
  }

  public String[] getSamples() {
    return samples;
  }

  public double getstdevFilter() {
    return stdevFilter;
  }

  public void goBackToIndexHistory(int indexInHistory) {
    if (markerIndexHistory[indexInHistory] == -1) {
      log.report("Cannot roll back to more than the history has recorded.");
    } else {
      markerIndex = markerIndexHistory[indexInHistory];
      updateMarkerIndexHistory();
    }
  }

  public boolean hideExcludedSamples(int index) {
    return excludeSamples[index];
  }

  public void last() {
    markerIndex = getAvailableMarker(true, false, showAllMarkersOrNot, showAnnotatedOrUnannotated);
    finishProcessing();
    // SwingUtilities.invokeLater(new Runnable() {
    // @Override
    // public void run() {
    // last.setEnabled(false);
    // next.setEnabled(false);
    // first.setEnabled(markerIndex > 0);
    // previous.setEnabled(markerIndex > 0);
    // }
    // });
  }

  public void loadAnnotationCollection() {
    String filename;

    filename = proj.ANNOTATION_FILENAME.getValue(false, false);
    if (new File(filename).exists()) {
      log.report("Loading annotation from: " + filename);
      annotationCollection = (AnnotationCollection) SerializedFiles.readSerial(filename);
    } else {
      log.report("Could not find annotation file: " + filename);
      annotationCollection = new AnnotationCollection();
      annotationCollection.addAnnotation('e', "Extra heterozygote clusters");
      annotationCollection.addAnnotation('m', "Monomorphic");
    }
    startAutoSaveToTempFile();

    annotationKeys = annotationCollection.getKeys();
  }

  public void loadCentroids() {
    String[] markerNames, files;
    Centroids trav;
    float[][][] travCents, targetCents;
    MarkerSet set;
    Vector<float[][][]> v;
    Vector<String> fileList;
    Hashtable<String, String> hash;
    int[] indices;

    files = Files.list(proj.DATA_DIRECTORY.getValue(false, true), ".cent", jar);
    // log.report("Found "+files.length+" .cent files in "+proj.getDir(Project.DATA_DIRECTORY));

    hash = new Hashtable<String, String>();
    for (int i = 0; i < markerList.length; i++) {
      hash.put(markerList[i], i + "");
    }
    set = proj.getMarkerSet();
    markerNames = set.getMarkerNames();
    indices = new int[markerList.length];
    for (int i = 0; i < markerNames.length; i++) {
      if (hash.containsKey(markerNames[i])) {
        indices[Integer.parseInt(hash.get(markerNames[i]))] = i;
      }
    }

    v = new Vector<float[][][]>();
    fileList = new Vector<String>();
    for (String file : files) {
      log.report("<", false, true);
      trav = Centroids.load(proj.DATA_DIRECTORY.getValue(false, true) + file, jar);
      log.report(">", false, true);
      if (trav.getFingerprint() != set.getFingerprint()) {
        log.reportError("Error - Centroids file '" + proj.DATA_DIRECTORY.getValue(false, true)
                        + file
                        + "' does not match up with the fingerprint of the current marker set and therefore will not load; if you don't want to see this error message again, then remove the .cent extension from this file.");
      } else {
        travCents = trav.getCentroids();
        targetCents = new float[markerList.length][][];
        for (int j = 0; j < markerList.length; j++) {
          targetCents[j] = travCents[indices[j]];
        }
        v.add(targetCents);
        fileList.add(ext.rootOf(file));
      }
    }
    log.report("");
    centroids = new float[v.size()][][][];
    for (int i = 0; i < centroids.length; i++) {
      centroids[i] = v.elementAt(i);
    }
    centList = Array.toStringArray(fileList);

  }

  public void loadMarkerListFromFile(String filename) {
    BufferedReader reader;
    Vector<String> markerNames = new Vector<String>();
    Vector<String> markerComments = new Vector<String>();
    String[] line;
    Vector<String> missingMarkers;

    missingMarkers = new Vector<String>();
    try {
      try {
        reader = Files.getReader(filename, jar, true, false);
        if (reader == null) {
          JOptionPane.showMessageDialog(null,
                                        "Failed to load '" + filename
                                              + "'; this is the designated filename in the project properties file. You will need to create a list of the markers that you want to review and place them in this file.",
                                        "Error", JOptionPane.ERROR_MESSAGE);
          return;
        }
        while (reader.ready()) {
          line = reader.readLine().trim().split("\t", -1);
          if (markerLookup.contains(line[0])) {
            markerNames.add(line[0]);
            markerComments.add(line.length > 1 ? line[1] : "");
          } else {
            missingMarkers.add(line[0]);
          }
        }
        reader.close();
        if (missingMarkers.size() > 0) {
          // JOptionPane.showMessageDialog(null, "Error - the following markers were not found in
          // the MarkerSet:\n"+Array.toStr(Array.toStringArray(missingMarkers), "\n"), "Error",
          // JOptionPane.ERROR_MESSAGE);
        }
      } catch (FileNotFoundException fnfe) {
        JOptionPane.showMessageDialog(null, "Error - could not find \"" + filename + "\"", "Error",
                                      JOptionPane.ERROR_MESSAGE);
        // handle error by closing window
      } catch (Exception e) {
        log.reportError("Error reading file \"" + filename + "\"");
        log.reportException(e);
        return;
      }

      masterMarkerList = Array.toStringArray(markerNames);
      masterCommentList = Array.toStringArray(markerComments);

      markerIndex = 0;
      markerIndexHistory = new int[NUM_MARKERS_TO_SAVE_IN_HISTORY];
    } catch (Exception e) {
      log.reportError("Error loading: " + filename);
      log.reportException(e);
    }
  }

  public void loadPlink(int plinkIndex) {
    String plinkFileRoot = proj.PLINK_DIR_FILEROOTS.getValue()[plinkIndex];
    if (plinkMarkerLoaders.containsKey(plinkFileRoot)) {
      return;
    } else {
      PlinkMarkerLoader pml =
                            PlinkMarkerLoader.loadPlinkDataFromListInSeparateThread(proj,
                                                                                    plinkFileRoot,
                                                                                    markerList,
                                                                                    sampleFIDIIDs);
      plinkMarkerLoaders.put(plinkFileRoot, pml);
    }
  }

  public boolean markerDataIsActive() {
    return markerDataLoader != null;
  }

  public boolean maskMissing(int index) {
    return maskMissing[index];
  }

  public void next() {
    markerIndex = getAvailableMarker(true, true, showAllMarkersOrNot, showAnnotatedOrUnannotated);
    finishProcessing();
    // SwingUtilities.invokeLater(new Runnable() {
    // @Override
    // public void run() {
    // last.setEnabled(markerIndex < markerList.length - 1);
    // next.setEnabled(markerIndex < markerList.length - 1);
    // first.setEnabled(markerIndex > 0);
    // previous.setEnabled(markerIndex > 0);
    // }
    // });
  }

  public void prepGCAdjustorLoad() {
    String[] gcParamFiles = proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue();
    ArrayList<String> display = new ArrayList<String>();
    ArrayList<String> actual = new ArrayList<String>();

    for (String gcParamFile : gcParamFiles) {
      if (Files.exists(gcParamFile)) {
        display.add(ext.rootOf(gcParamFile));
        actual.add(gcParamFile);
      }
    }
    gcCList = new String[][] {Array.toStringArray(display), Array.toStringArray(actual)};
    gcAdjustorParameters = new GcAdjustorParameters[gcCList[0].length];
  }

  // private byte qc_prevChr;
  // private int[] qc_prevGeno;
  // private String[] qc_prevSex;
  // private String[] qc_prevOthCls;
  // private int qc_prevIndex;

  public void previous() {
    markerIndex = getAvailableMarker(false, true, showAllMarkersOrNot, showAnnotatedOrUnannotated);
    finishProcessing();
    // SwingUtilities.invokeLater(new Runnable() {
    // @Override
    // public void run() {
    // last.setEnabled(markerIndex < markerList.length - 1);
    // next.setEnabled(markerIndex < markerList.length - 1);
    // first.setEnabled(markerIndex > 0);
    // previous.setEnabled(markerIndex > 0);
    // }
    // });
  }



  // private byte qc_prevChr;
  // private int[] qc_prevGeno;
  // private String[] qc_prevSex;
  // private String[] qc_prevOthCls;
  // private int qc_prevIndex;

  public void revertMarkersToOriginalList(int indexOfDefaultFirstMarker) {
    if (!saveClusterFilterAndAnnotationCollection()) {
      return;
    }

    for (JCheckBox annotationCheckBoxe : annotationCheckBoxes) {
      annotationCheckBoxe.setEnabled(true);
    }
    indexOfAnnotationUsedAsMarkerList = -1;
    markerList = masterMarkerList;
    commentList = masterCommentList;
    loadMarkerDataFromList(indexOfDefaultFirstMarker);
    displayIndex(navigationField);
    updateGUI();
  }

  public boolean saveClusterFilterAndAnnotationCollection() {
    String[] options;
    int choice;
    String annotationFilename;
    String clusterFilterFilename;
    boolean firstTimeClusters;
    String title, message;

    clusterFilterFilename = proj.CLUSTER_FILTER_COLLECTION_FILENAME.getValue(false, false);
    firstTimeClusters = Files.exists(clusterFilterFilename, jar);
    options = new String[] {"Yes, overwrite", "No"};
    if (isClusterFilterUpdated) {
      if (firstTimeClusters) {
        title = "Save cluster filters?";
        message =
                "You first cluster filters for this project have been created. Would you like to save them to a permanent file?";
      } else {
        title = "Update filters?";
        message =
                "New cluster filters have been created. Would you like to append them to the permanent file?";
      }
      choice =
             JOptionPane.showOptionDialog(null, message, title, JOptionPane.YES_NO_CANCEL_OPTION,
                                          JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      if (choice == 0) {
        clusterFilterCollection.serialize(clusterFilterFilename);
        proj.archiveFile(clusterFilterFilename);
        clusterFilterCollection.serialize(clusterFilterFilename);
        log.report("clusterFilters should be backed up now");
      } else if (choice == -1) {
        return false;
      } else {
        // TODO As a double security, move sessionID + ".tempClusterFilters.ser" to
        // BACKUP_DIRECTORY. But then need to delete it from BACKUP_DIRECTORY at some point of time.
        // need a rule for that. also need the code for deletion.
      }
      setClusterFilterUpdated(false);
    }

    annotationFilename = proj.ANNOTATION_FILENAME.getValue(false, false);
    Files.exists(annotationFilename, jar);
    if (isAnnotationUpdated) {
      choice =
             JOptionPane.showOptionDialog(null,
                                          "New Annotations have been generated. Would you like to append them to the permanent annotations file?",
                                          "Update annotations?", JOptionPane.YES_NO_CANCEL_OPTION,
                                          JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      if (choice == 0) {
        annotationCollection.serialize(annotationFilename);
        proj.archiveFile(annotationFilename);
        annotationCollection.serialize(annotationFilename);
      } else if (choice == -1) {
        return false;
      }
      setAnnotationUpdated(false);
    }

    return true;
  }

  public void setAnnotationUpdated(boolean status) {
    isAnnotationUpdated = status;
    autoSave.setAnnotationUpdated(status);
  }

  public void setClusterFilterUpdated(boolean status) {
    isClusterFilterUpdated = status;
    autoSave.setClusterFilterUpdated(status);
  }

  public void setCurrentClusterFilter() {
    setCurrentClusterFilter((byte) (clusterFilterCollection.getSize(getMarkerName()) - 1));
  }

  public void setCurrentClusterFilter(byte currentClusterFilter) {
    ClusterFilter currentFilter;

    this.currentClusterFilter = currentClusterFilter;

    // if (currentClusterFilter >= 0) {
    // log.report("Before: image=="+(scatPanel.image==null?"null":"BImg")+"\t");
    if (clusterFilterCollection.getSize(getMarkerName()) > 0) {
      newGenotype.setSelectedIndex(clusterFilterCollection.getGenotype(markerList[markerIndex],
                                                                       currentClusterFilter)
                                   + 1);
      // if (seletedScatterPanel.getRectangles()!=null) {
      // seletedScatterPanel.generateRectangles();
      // seletedScatterPanel.rectangles[currentClusterFilter].setColor((byte)0);
      // }
      for (ScatterPanel scatterPanel : scatterPanels) {
        if (scatterPanel.getRectangles() != null) {
          scatterPanel.generateRectangles();
          scatterPanel.rectangles[currentClusterFilter].setColor((byte) 0);
        }
      }
      // }
      //
      // if (clusterFilterCollection.getSize(getMarkerName()) > 0) {
      currentFilter = clusterFilterCollection.getClusterFilters(getMarkerName())
                                             .get(currentClusterFilter);
      // seletedScatterPanel.highlightPoints(getCurrentMarkerData().getHighlightStatus(currentFilter));
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.highlightPoints(getCurrentMarkerData().getHighlightStatus(currentFilter));
      }
    }
    // log.report("After: image=="+(scatPanel.image==null?"null":"BImg"));

  }

  public void setNumComponents(int numComponents) {
    this.numComponents = numComponents;
  }

  public void startAutoSaveToTempFile() {
    if (autoSave == null) {
      autoSave = new AutoSaveForScatterPlot(clusterFilterCollection,
                                            proj.DATA_DIRECTORY.getValue(false, true) + sessionID
                                                                     + ".tempClusterFilters.ser",
                                            annotationCollection,
                                            proj.DATA_DIRECTORY.getValue(false, true) + sessionID
                                                                  + ".tempAnnotation.ser",
                                            30);
      new Thread(autoSave).start();
    } else if (clusterFilterCollection != null && autoSave.isClusterFilterNull()) {
      autoSave.addToAutoSave(clusterFilterCollection, proj.DATA_DIRECTORY.getValue(false, true)
                                                      + sessionID + ".tempClusterFilters.ser");
    } else if (annotationCollection != null && autoSave.isAnnotationNull()) {
      autoSave.addToAutoSave(annotationCollection, proj.DATA_DIRECTORY.getValue(false, true)
                                                   + sessionID + ".tempAnnotation.ser");
    }
    autoSave.saveNow();
  }

  public boolean symmetricAxes(int index) {
    return symmetry[index];
  }

  public void uncheckAnnotationBox(char c) {
    int index;

    index = ext.indexOfChar(c, annotationKeys);
    if (annotationCheckBoxes[index].isSelected()) {
      annotationCheckBoxes[index].setSelected(false);
    }
  }

  public void updateAnnotationPanelAnnotationCheckBoxes() {
    isInitilizing = true;
    for (int i = 0; i < annotationCheckBoxes.length; i++) {
      annotationCheckBoxes[i].setText(annotationCollection.getDescriptionForComment(annotationKeys[i],
                                                                                    showAnnotationShortcuts,
                                                                                    true));
      if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
        annotationCheckBoxes[i].setSelected(true);
      } else {
        annotationCheckBoxes[i].setSelected(false);
      }
    }
    isInitilizing = false;
  }

  public void updateAnnotationPanelFilterRadioButtons() {
    int numAnnotated;

    numAnnotated = 0;
    for (String element : markerList) {
      if (annotationCollection.markerHasAnyAnnotation(element)) {
        numAnnotated++;
      }
    }

    for (int i = 0; i < filterRadioButtons.length; i++) {
      filterRadioButtons[i].setText(RADIOBUTTON_TEXTS[i] + " (n="
                                    + (i == 0 ? markerList.length
                                              : (i == 1 ? numAnnotated
                                                        : (markerList.length - numAnnotated)))
                                    + ")");
    }
  }

  public void updateCentLabels() {
    double[] comp;
    String str;
    MarkerData markerData;

    if (markerDataLoader != null) {
      markerData = getCurrentMarkerData();
      if (markerData.getLRRs() != null) {
        for (int i = 0; i < centList.length; i++) {
          comp = markerData.compareLRRs(centroids[i][markerIndex], log);
          str = comp[0] + "";

          for (int j = 2; j < str.length(); j++) {
            if (str.charAt(j) != '9') {
              str = str.substring(0, Math.min(j + 2, str.length()));
            }
          }
          if (centroidLabels != null) {
            centroidLabels[i].setText("LRR corr=" + str + ", err=" + ext.formDeci(comp[1], 3));
          }
        }
      }
    }
  }

  public void updateColorKey(Hashtable<String, String> hash, int index) {
    if (index == selectedPanelIndex) {
      colorKeyPanel.setSisterPanel(scatterPanels[index]);
      colorKeyPanel.updateColorKey(hash);
    }
  }

  public void updateCurrentClusterFilterGenotype(byte newGenotypeSelected,
                                                 boolean updateGenotypeComboBox) {
    clusterFilterCollection.updateGenotype(getMarkerName(), currentClusterFilter,
                                           newGenotypeSelected);
    setClusterFilterUpdated(true);
    // seletedScatterPanel.setPointsGeneratable(true);
    // seletedScatterPanel.setQcPanelUpdatable(true);
    // seletedScatterPanel.paintAgain();
    for (ScatterPanel scatterPanel : scatterPanels) {
      scatterPanel.setPointsGeneratable(true);
      scatterPanel.setUpdateQCPanel(true);
      scatterPanel.paintAgain();
    }
    if (updateGenotypeComboBox) {
      newGenotype.setSelectedIndex(newGenotypeSelected + 1);
    }
  }

  public void updateCurrentClusterFiltersGenotype() {
    byte currentGenotype, newGenotype;

    currentGenotype = clusterFilterCollection.getGenotype(getMarkerName(), currentClusterFilter);
    if (currentGenotype == (ScatterPlot.GENOTYPE_OPTIONS.length - 2)) {
      newGenotype = -1;
    } else {
      newGenotype = (byte) (currentGenotype + 1);
    }
    updateCurrentClusterFilterGenotype(newGenotype, true);
  }


  public synchronized void updateGUI() {
    if (markerDataLoader == null) {
      return;
    }

    // if (markerData[markerIndex]==null) {
    // log.report("Marker data is null at index "+markerIndex);
    //// create visual that we're loading data, (write to scatter panel "Loading data from file")
    // loadMarkerData(markerIndex);
    // }

    if (markerList.length == 0) {
      markerName.setText("Error: marker data was not successfully loaded");
      commentLabel.setText("Check to make sure MarkerLookup is synchronized with the current data");
    } else {
      markerName.setText(markerList[markerIndex]);
      commentLabel.setText(commentList[markerIndex]);
    }
    if (plot_types[selectedPanelIndex] >= 1 && plot_types[selectedPanelIndex] != 4) {
      // if (plot_types[selectedPanelIndex] >= 2 && plot_types[selectedPanelIndex] != 5) {
      symmetryBox.setEnabled(false);
      symmetryBox.setSelected(false);
      // symmetry[selectedPanelIndex] = false;
      scatterPanels[selectedPanelIndex].setSymmetricAxes(false);
    } else {
      symmetryBox.setEnabled(true);
      symmetryBox.setSelected(symmetricAxes(selectedPanelIndex));
      // scatterPanels[selectedPanelIndex].setSymmetricAxes(symmetryBox.isSelected());
      scatterPanels[selectedPanelIndex].setSymmetricAxes(symmetricAxes(selectedPanelIndex));
    }
    // BAF/LRR ???
    if (plot_types[selectedPanelIndex] == 2) {
      // if (plot_types[selectedPanelIndex] == 3) {
      boolean recomputed = false;
      for (int i = 0; i < centBoxes.length; i++) {
        if (centBoxes[i].isSelected()) {
          if (!recomputed) {
            getCurrentMarkerData().recompute(centroids[i][markerIndex]);
            recomputed = true;
          } else {
            centBoxes[i].setSelected(false);
          }
        }

      }
    }
    // if(plot_types[selectedPanelIndex] == 5){
    // symmetryBox.setEnabled(true);
    //// scatterPanels[selectedPanelIndex].setSymmetricAxes(symmetryBox.isSelected());
    // scatterPanels[selectedPanelIndex].setSymmetricAxes(symmetricAxes(selectedPanelIndex));
    // }
    if (filterRadioButtons != null) {
      updateAnnotationPanelFilterRadioButtons();
      updateAnnotationPanelAnnotationCheckBoxes();
    }

    if (markerIndex != previousMarkerIndex) {
      updateCentLabels();
    }
    previousMarkerIndex = markerIndex;

    if (showingAll) {
      for (int i = 0; i < indivPanels.length; i++) {
        String plotTitle = ScatterPlot.TYPES[plot_types[i]][0] + "/"
                           + ScatterPlot.TYPES[plot_types[i]][1] + " - "
                           + sampleData.getClassName(classes[i]);
        if (i == selectedPanelIndex) {
          TitledBorder tb = BorderFactory.createTitledBorder(
                                                             BorderFactory.createEtchedBorder(EtchedBorder.LOWERED,
                                                                                              Color.RED,
                                                                                              Color.GRAY),
                                                             plotTitle, TitledBorder.LEADING,
                                                             TitledBorder.TOP);
          indivPanels[i].setBorder(tb);
        } else {
          TitledBorder tb = BorderFactory.createTitledBorder(
                                                             BorderFactory.createEtchedBorder(EtchedBorder.LOWERED,
                                                                                              Color.GRAY.brighter(),
                                                                                              Color.GRAY),
                                                             plotTitle, TitledBorder.LEADING,
                                                             TitledBorder.TOP);
          indivPanels[i].setBorder(tb);
        }
      }
    } else {
      indivPanels[selectedPanelIndex].setBorder(null);
    }

    // long t2 = System.currentTimeMillis();
    if (showingAll) {
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.createImage();
      }
      for (ScatterPanel scatterPanel : scatterPanels) {
        scatterPanel.repaint();
        // scatterPanels[i].paintAgain();
      }
    } else {
      scatterPanels[selectedPanelIndex].paintAgain();
    }

    // long t3 = System.currentTimeMillis();
    // System.out.println("Updated: ");
    // System.out.println("\tt1: " + t1);
    // System.out.println("\tt2: " + t2);
    // System.out.println("\tt3: " + t3);
  }

  public void updateMarkerIndexHistory() {
    // Oldest history is stored in the end of the array
    for (int i = 0; i < markerIndexHistory.length - 1; i++) {
      markerIndexHistory[i + 1] = markerIndexHistory[i];
    }
    markerIndexHistory[0] = markerIndex;
  }

  public void updateQcPanel(byte chr, int[] genotype, String[] sex, String[] otherClass,
                            int index) {
    // if (qc_prevGeno == null ||
    // qc_prevSex == null ||
    // qc_prevOthCls == null) {
    // qc_prevChr = chr;
    // qc_prevGeno = genotype;
    // qc_prevSex = sex;
    // qc_prevOthCls = otherClass;
    // qc_prevIndex = index;
    // } else {
    // boolean update = false;
    // if (chr != qc_prevChr) {
    // update = true;
    // } else if (index != qc_prevIndex) {
    // update = true;
    // } else if (!Array.equals(genotype, qc_prevGeno)) {
    // update = true;
    // } else if (!Array.equals(sex, qc_prevSex, false)) {
    // update = true;
    // } else if (!Array.equals(otherClass, qc_prevOthCls, false)) {
    // update = true;
    // }
    // if (!update) {
    // System.out.println("Not Updating QC Panel");
    // return;
    // }
    // qc_prevChr = chr;
    // qc_prevGeno = genotype;
    // qc_prevSex = sex;
    // qc_prevOthCls = otherClass;
    // qc_prevIndex = index;
    // }


    int numCalledGenotypes;
    double callrate;
    JLabel qcPanelLabel;
    // JLabel qcCallRateLabel;u
    // JLabel qcHwePvalueLabel;u
    int[] alleleCounts;
    double hweP;
    double minorAlleleFrequency;
    // int[][] sexContingecyTable = new int[2][2];
    CTable classCount;
    String[] called;
    int currentClass;
    int numIncludedSamples;

    numIncludedSamples = 0;
    numCalledGenotypes = 0;
    alleleCounts = new int[3];
    called = new String[genotype.length];
    for (int element : genotype) {
      numIncludedSamples += element == -3 ? 0 : 1;
    }
    if (chr == 23) {
      for (int i = 0; i < genotype.length; i++) {
        if (genotype[i] <= 0) {
          called[i] = "-1";
        } else {
          called[i] = "1";
          if (sex[i].equals("2")) {
            alleleCounts[genotype[i] - 1]++;
          }
          numCalledGenotypes++;
        }
      }
    } else if (chr == 24) {
      for (int i = 0; i < genotype.length; i++) {
        if (genotype[i] <= 0) {
          called[i] = "-1";
        } else {
          called[i] = "1";
          numCalledGenotypes++;
        }
      }
    } else {
      for (int i = 0; i < genotype.length; i++) {
        if (genotype[i] <= 0) {
          called[i] = "-1";
        } else {
          called[i] = "1";
          alleleCounts[genotype[i] - 1]++;
          numCalledGenotypes++;
        }

      }
    }
    hweP = AlleleFreq.HWEsig(alleleCounts);
    callrate = (double) numCalledGenotypes / (double) numIncludedSamples * 100.0;
    if (alleleCounts[0] > alleleCounts[2]) {
      minorAlleleFrequency = (double) (alleleCounts[1] + alleleCounts[2])
                             / (alleleCounts[0] + 2 * alleleCounts[1] + alleleCounts[2]);
    } else {
      minorAlleleFrequency = (double) (alleleCounts[0] + alleleCounts[1])
                             / (alleleCounts[0] + 2 * alleleCounts[1] + alleleCounts[2]);
    }

    qcPanel.removeAll();
    qcPanel.repaint();

    qcPanelLabel = new JLabel("Chromosome:", JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 0");
    qcPanelLabel = new JLabel("" + chr, JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 0");

    qcPanelLabel = new JLabel("Callrate:", JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 1");
    qcPanelLabel = new JLabel(ext.formDeci(callrate, 4) + "%", JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 1");

    qcPanelLabel = new JLabel("Minor Allele Frequency:   ", JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 2");
    qcPanelLabel = new JLabel(ext.prettyP(minorAlleleFrequency + "", 4, 5, 1, false), JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 2");

    qcPanelLabel = new JLabel("HWE p-value:", JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 3");
    qcPanelLabel = new JLabel(ext.prettyP(hweP), JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 3");

    ToolTipManager.sharedInstance().setDismissDelay(100000);

    classCount = new CTable(called, sex);// This is the problem.
    classCount.setCustomNullValues(Array.addStrToArray("-1", CTable.DEFAULT_NULL_VALUES));
    classCount.setCustomLabelsAndOrder(new String[][] {{"-1", "Genotype missing"},
                                                       {"1", "Genotype NOT missing"}},
                                       sampleData.getActualClassColorKey(0));

    // classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype
    // NOT missing"}}, Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null,
    // "missing"}));

    qcPanelLabel = new JLabel("Callrate by sex:", JLabel.LEFT);
    qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 4");
    qcPanelLabel = new JLabel(
                              ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(),
                                                                                      false),
                                                           1)),
                              JLabel.LEFT);
    qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 4");

    classCount = new CTable(CTable.extrapolateCounts(sex, genotype));
    classCount.setCustomNullValues(Array.addStrToArray("-1", CTable.DEFAULT_NULL_VALUES));
    classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(0),
                                                     new String[] {null, "missing"}),
                                       new String[][] {{"A", "Allele A"}, {"B", "Allele B"}});

    // classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(0), new
    // String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"},
    // {".","Missing"}});
    qcPanelLabel = new JLabel("Allele Freq by sex: ", JLabel.LEFT);
    qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 5");
    qcPanelLabel = new JLabel(
                              ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(),
                                                                                      false),
                                                           1)),
                              JLabel.LEFT);
    qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 5");

    currentClass = getCurrentClass(index);
    // currentClass = getCurrentClass();
    if (currentClass > SampleData.BASIC_CLASSES.length
        && currentClass < sampleData.getBasicClasses().length + sampleData.getNumActualClasses()) {
      classCount = new CTable(called, otherClass);// This is the problem.
      classCount.setCustomLabelsAndOrder(new String[][] {{"-1", "Genotype missing"},
                                                         {"1", "Genotype NOT missing"}},
                                         sampleData.getActualClassColorKey(currentClass
                                                                           - SampleData.BASIC_CLASSES.length));
      qcPanelLabel = new JLabel("Callrate by " + sampleData.getClassName(currentClass) + ": ",
                                JLabel.LEFT);
      classCount.setCustomNullValues(Array.addStrToArray("-1",
                                                         Array.addStrToArray("0",
                                                                             CTable.DEFAULT_NULL_VALUES)));
      // classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"},
      // {"1","Genotype NOT missing"}},
      // Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length),
      // new String[] {null, "missing"}));
      qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
      qcPanelLabel.setFont(new Font("Arial", 0, 14));
      qcPanel.add(qcPanelLabel, "cell 0 6");
      // qcPanelLabel = new JLabel("Callrate by "+sampleData.getClassName(currentClass)+":
      // "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(),
      // false), 1) ), JLabel.LEFT);
      // classCount.setCustomNullValues(Array.addStrToArray("-1", Array.addStrToArray("0",
      // CTable.DEFAULT_NULL_VALUES)));
      // //classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"},
      // {"1","Genotype NOT missing"}},
      // Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length),
      // new String[] {null, "missing"}));
      // qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
      // qcPanelLabel.setFont(new Font("Arial", 0, 14));
      // qcPanel.add(qcPanelLabel);
      qcPanelLabel = new JLabel(
                                ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(),
                                                                                        false),
                                                             1)),
                                JLabel.LEFT);
      classCount.setCustomNullValues(Array.addStrToArray("-1",
                                                         Array.addStrToArray("0",
                                                                             CTable.DEFAULT_NULL_VALUES)));
      // classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"},
      // {"1","Genotype NOT missing"}},
      // Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length),
      // new String[] {null, "missing"}));
      qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
      qcPanelLabel.setFont(new Font("Arial", 0, 14));
      qcPanel.add(qcPanelLabel, "cell 1 6");

      classCount = new CTable(CTable.extrapolateCounts(otherClass, genotype));
      classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass
                                                                                         - SampleData.BASIC_CLASSES.length),
                                                       new String[] {null, "missing"}),
                                         new String[][] {{"A", "Allele A"}, {"B", "Allele B"}});
      // classCount.replaceIdWithLabel(SampleData.KEYS_FOR_BASIC_CLASSES[1],sampleData.getActualClassColorKey(0));
      qcPanelLabel = new JLabel("Allele Freq by " + sampleData.getClassName(currentClass) + ": ",
                                JLabel.LEFT);
      classCount.setCustomNullValues(Array.addStrToArray("-1",
                                                         Array.addStrToArray("0",
                                                                             CTable.DEFAULT_NULL_VALUES)));
      // classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length),
      // new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
      qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
      qcPanelLabel.setFont(new Font("Arial", 0, 14));
      qcPanel.add(qcPanelLabel, "cell 0 7");
      // qcPanelLabel = new JLabel("Allele Freq by "+sampleData.getClassName(currentClass)+":
      // "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(),
      // false), 1)), JLabel.LEFT);
      // classCount.setCustomNullValues(Array.addStrToArray("-1", Array.addStrToArray("0",
      // CTable.DEFAULT_NULL_VALUES)));
      // //classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length),
      // new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
      // qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
      // qcPanelLabel.setFont(new Font("Arial", 0, 14));
      // qcPanel.add(qcPanelLabel);
      qcPanelLabel = new JLabel(
                                ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(),
                                                                                        false),
                                                             1)),
                                JLabel.LEFT);
      classCount.setCustomNullValues(Array.addStrToArray("-1",
                                                         Array.addStrToArray("0",
                                                                             CTable.DEFAULT_NULL_VALUES)));
      // classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length),
      // new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
      qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
      qcPanelLabel.setFont(new Font("Arial", 0, 14));
      qcPanel.add(qcPanelLabel, "cell 1 7");
    }
    qcPanelLabel = new JLabel("Marker GC content:", JLabel.LEFT);
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 0 8");
    if (gcAnnotations != null) {
      String gc = gcAnnotations[markerIndex].getAnnotations()[0].getData();
      gc = gc.substring(0, Math.min(gc.length(), 4));
      qcPanelLabel = new JLabel(gc, JLabel.LEFT);
    } else {
      qcPanelLabel = new JLabel("NA", JLabel.LEFT);
    }
    qcPanelLabel.setFont(new Font("Arial", 0, 14));
    qcPanel.add(qcPanelLabel, "cell 1 8");

    // qcPanelLabel = new JLabel("Minor Allele Freq: " + (new
    // DecimalFormat("#.####").format(classCount.getMinorAlleleFrequency())), JLabel.LEFT);
    // qcPanelLabel.setFont(new Font("Arial", 0, 14));
    // qcPanel.add(qcPanelLabel);

    // AlleleFreq.calcFrequency(genotypes);


    /*
     * qcPanelLabel = new JLabel(""+sexContingecyTable[0][1], JLabel.LEFT); qcPanelLabel.setFont(new
     * Font("Arial", 0, 14)); qcPanel.add(qcPanelLabel); qcPanelLabel = new
     * JLabel(""+sexContingecyTable[1][0], JLabel.LEFT); qcPanelLabel.setFont(new Font("Arial", 0,
     * 14)); qcPanel.add(qcPanelLabel); qcPanelLabel = new JLabel(""+sexContingecyTable[1][1],
     * JLabel.LEFT); qcPanelLabel.setFont(new Font("Arial", 0, 14)); qcPanel.add(qcPanelLabel);
     */

    qcPanel.validate();
    qcPanel.repaint();
  }

  @Override
  public void windowActivated(WindowEvent e) {}

  @Override
  public void windowClosed(WindowEvent e) {}

  @Override
  public void windowClosing(WindowEvent e) {
    int count;

    if (!saveClusterFilterAndAnnotationCollection()) {
      return;
    }

    // check blast marker filter change
    if (blastFrame != null) {
      double ratio = ((double) blastFrame.currentAlignFilter)
                     / ((double) proj.ARRAY_TYPE.getValue().getProbeLength());
      if (ratio != proj.BLAST_PROPORTION_MATCH_FILTER.getValue().doubleValue()) {
        proj.BLAST_PROPORTION_MATCH_FILTER.setValue(ratio);
        proj.saveProperties(new Property[] {proj.BLAST_PROPORTION_MATCH_FILTER});
      }
    }

    // TODO notify all threads (e.g., MarkerDataLoader) that they need to close

    if (autoSave != null) {
      autoSave.kill();
    }
    // Killing doesn't actually do anything...
    // for (PlinkMarkerLoader pml : plinkMarkerLoaders.values()) {
    // pml.kill();
    // count = 0;
    // while (!pml.killComplete()) {
    // try {
    // Thread.sleep(250);
    // } catch (InterruptedException ie) {
    // }
    // count++;
    // if (count > 0) {
    // log.reportError("Waiting for PlinkMarkerLoader to wind down...");
    // }
    // }
    // }
    if (markerDataLoader != null) {
      markerDataLoader.kill();
      count = 0;
      while (!markerDataLoader.killComplete()) {
        try {
          Thread.sleep(250);
        } catch (InterruptedException ie) {
        }
        count++;
        if (count > 0) {
          log.reportError("Waiting for markerDataLoader to wind down...");
        }
      }
    }

    if (exitOnClose) {
      System.exit(1);
    } else {
      if (blastFrame != null) {
        blastFrame.setVisible(false);
        blastFrame.dispose();
      }
      ((JFrame) (this)).dispose();
    }
  }

  @Override
  public void windowDeactivated(WindowEvent e) {}

  @Override
  public void windowDeiconified(WindowEvent e) {}

  @Override
  public void windowIconified(WindowEvent e) {}

  @Override
  public void windowOpened(WindowEvent e) {}


  private void activateAllAnnotationMaps() {
    for (char annotationKey : annotationKeys) {
      addAnnotationToMaps(annotationKey);
    }
  }

  private void addAnnotationToMaps(char c) {
    // seletedScatterPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0),
    // 0), "annotation\t"+c);
    // seletedScatterPanel.getActionMap().put("annotation\t"+c, new AnnotationAction(this, c,
    // true));
    //
    // seletedScatterPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0),
    // InputEvent.SHIFT_MASK), "annotation\tShift_"+c);
    // seletedScatterPanel.getActionMap().put("annotation\tShift_"+c, new AnnotationAction(this, c,
    // false));
    for (ScatterPanel scatterPanel : scatterPanels) {
      scatterPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
                  .put(KeyStroke.getKeyStroke((c + "").toUpperCase().charAt(0), 0),
                       "annotation\t" + c);
      scatterPanel.getActionMap().put("annotation\t" + c, new AnnotationAction(this, c, true));

      scatterPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
                  .put(KeyStroke.getKeyStroke((c + "").toUpperCase().charAt(0),
                                              InputEvent.SHIFT_MASK),
                       "annotation\tShift_" + c);
      scatterPanel.getActionMap().put("annotation\tShift_" + c,
                                      new AnnotationAction(this, c, false));
    }
  }

  private void annotationPanel() {
    annotationPanel = new JPanel();
    // annotationScrollPane = new JScrollPane(annotationPanel,
    // JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    // annotationScrollPane.setVisible(true);
    // annotationScrollPane.createVerticalScrollBar();

    annotationPanel.setBackground(Color.WHITE);
    annotationPanel.setLayout(new BoxLayout(annotationPanel, BoxLayout.Y_AXIS));
    // annotationPanelLayout = new SpringLayout();
    // annotationPanel.setLayout(annotationPanelLayout);
    annotationPanel.add(annotationPanelUpper());
    // annotationPanelLayout.putConstraint(SpringLayout.WEST, temp1, 5, SpringLayout.WEST,
    // annotationPanel);
    // annotationPanelLayout.putConstraint(SpringLayout.NORTH, temp1, 5, SpringLayout.NORTH,
    // annotationPanel);
    annotationPanelLowerPart = new JPanel();
    annotationPanelLowerPartLayout = new SpringLayout();
    annotationPanelLowerPart.setLayout(annotationPanelLowerPartLayout);
    annotationPanelLowerPart.setBackground(BACKGROUND_COLOR);
    annotationPanel.add(annotationPanelLower());
    // annotationPanelLayout.putConstraint(SpringLayout.WEST, annotationPanelLowerPart, 5,
    // SpringLayout.WEST, annotationPanel);
    // annotationPanelLayout.putConstraint(SpringLayout.NORTH, annotationPanelLowerPart, 25,
    // SpringLayout.NORTH, annotationPanel);
    annotationScrollPane = new JScrollPane(annotationPanel);
    annotationScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
    annotationScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
  }

  private JPanel annotationPanelUpper() {
    JPanel annotationPanelUpper;
    JCheckBox checkBox;
    JButton button;

    annotationPanelUpper = new JPanel();
    annotationPanelUpper.setLayout(new BoxLayout(annotationPanelUpper, BoxLayout.X_AXIS));
    annotationPanelUpper.setBackground(BACKGROUND_COLOR);

    checkBox = new JCheckBox("Shortcut");
    checkBox.setBackground(Color.WHITE);
    checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
    checkBox.setMnemonic(KeyEvent.VK_S);
    checkBox.setSelected(true);
    checkBox.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent e) {
        JCheckBox checkBox = (JCheckBox) e.getSource();
        if (checkBox.isSelected()) {
          showAnnotationShortcuts = true;
          activateAllAnnotationMaps();
        } else {
          showAnnotationShortcuts = false;
          deactivateAllAnnotationMaps();
        }
      }
    });
    annotationPanelUpper.add(checkBox);

    checkBox = new JCheckBox("Auto Advance");
    checkBox.setBackground(Color.WHITE);
    checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
    checkBox.setMnemonic(KeyEvent.VK_A);
    checkBox.setSelected(true);
    checkBox.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent e) {
        JCheckBox checkBox = (JCheckBox) e.getSource();
        annotationAutoAdv = checkBox.isSelected();
      }
    });
    annotationPanelUpper.add(checkBox);

    // button = new JButton("Import");
    button = new JButton(Grafik.getImageIcon("images/import-icon.png"));
    button.setToolTipText("Import Annotations");
    button.setBackground(Color.WHITE);
    button.setBorder(null);
    button.setSize(10, 11);
    button.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand() != null) {
          String filename, filenameBak;
          JFileChooser fileChooser;
          int returnVal;
          String[] options;
          int choice;

          fileChooser = new JFileChooser(new File(proj.PROJECT_DIRECTORY.getValue()));
          returnVal = fileChooser.showOpenDialog(null);

          if (returnVal == JFileChooser.APPROVE_OPTION) {
            filename = fileChooser.getSelectedFile().getAbsolutePath();
            if (isAnnotationUpdated) {
              options = new String[] {"Yes/Overwrite", "No"};
              choice = JOptionPane.showOptionDialog(null,
                                                    "New Annotations have been generated. Do you want to save them before importing new annotations?",
                                                    "Save new annotations?",
                                                    JOptionPane.YES_NO_CANCEL_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE, null, options,
                                                    options[0]);
              if (choice == 0) {
                // filenameBak = proj.getFilename(proj.getProjectDir(), false, false) +
                // "annotations_bak_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date()))
                // + ".ser";
                filenameBak = proj.PROJECT_DIRECTORY.getValue(false, false) + "annotations_bak_"
                              + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date()))
                              + ".ser";
                log.report("Writing annotations to " + filenameBak);
                annotationCollection.serialize(filenameBak);
              }
            }
            deactivateAllAnnotationMaps();
            if (annotationCollection == null) {
              annotationCollection = AnnotationCollection.loadFromLists(filename, null);
              startAutoSaveToTempFile();
            } else {
              AnnotationCollection.appendFromLists(filename, annotationCollection,
                                                   proj.getMarkerNames(), null);
            }
            setAnnotationUpdated(true);
            annotationKeys = annotationCollection.getKeys();
            activateAllAnnotationMaps();
            annotationPanelLower();
          }
        }
      }
    });
    annotationPanelUpper.add(button);

    // button = new JButton("Export");
    button = new JButton(Grafik.getImageIcon("images/export-icon.png"));
    button.setToolTipText("Export Annotations");
    button.setBackground(Color.WHITE);
    button.setBorder(null);
    button.setSize(10, 11);
    button.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand() != null) {
          annotationCollection.dumpLists(proj);
        }
      }
    });
    annotationPanelUpper.add(button);

    return annotationPanelUpper;
  }

  private JPanel centroidPanel() {
    JPanel centroidPanel;

    centroidPanel = new JPanel();
    // tabPanel.setLayout(new BoxLayout(tabPanel, BoxLayout.Y_AXIS));
    centroidPanel.setLayout(new GridBagLayout());
    centroidPanel.setSize(50, 100);
    centroidPanel.setBackground(BACKGROUND_COLOR);

    GridBagConstraints gbc = new GridBagConstraints();
    gbc.insets = new Insets(1, 3, 0, 30);
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridwidth = GridBagConstraints.REMAINDER;

    ItemListener centListener = new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent ie) {
        int index = ext.indexOfStr(((JCheckBox) ie.getSource()).getText(), centList);
        displayCentroids[index] = ((JCheckBox) ie.getSource()).isSelected();
        centroidLabels[index].setVisible(displayCentroids[index]);
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
        }
        updateGUI();
      }
    };
    JLabel label = new JLabel("  ");
    label.setFont(new Font("Arial", 0, 20));
    centroidPanel.add(label, gbc);
    label = new JLabel("Centroids:");
    label.setFont(new Font("Arial", 0, 20));
    label.setHorizontalAlignment(JLabel.CENTER);
    centroidPanel.add(label, gbc);
    centBoxes = new JCheckBox[centList.length];
    displayCentroids = new boolean[centList.length];
    centroidLabels = new JLabel[centList.length];
    for (int i = 0; i < centList.length; i++) {
      centBoxes[i] = new JCheckBox(centList[i]);
      centBoxes[i].setFont(new Font("Arial", 0, 14));
      centBoxes[i].setSelected(displayCentroids[i]);
      centBoxes[i].addItemListener(centListener);
      centBoxes[i].setBorder(BorderFactory.createLineBorder(ScatterPanel.DEFAULT_COLORS[5 + i], 5));
      centBoxes[i].setBorderPainted(true);
      centBoxes[i].setBackground(BACKGROUND_COLOR);
      centroidPanel.add(centBoxes[i], gbc);

      centroidLabels[i] = new JLabel("LRR correlation not performed");
      centroidLabels[i].setVisible(displayCentroids[i]);
      centroidPanel.add(centroidLabels[i], gbc);
    }

    return centroidPanel;
  }

  private JPanel clusterFilterPanel() {
    JPanel clusterFilterPanel = new JPanel();
    clusterFilterPanel.setBackground(BACKGROUND_COLOR);
    // clusterFilterPanel.setLayout(new BoxLayout(clusterFilterPanel, BoxLayout.PAGE_AXIS));
    clusterFilterPanel.setSize(50, 10);

    // JButton first1 = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
    // first1.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
    // first1.addActionListener(this);
    // first1.setActionCommand(FIRST);
    // first1.setPreferredSize(new Dimension(20, 20));
    JButton backward = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif"));
    backward.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif"));
    backward.addActionListener(this);
    backward.setActionCommand(CLUSTER_FILTER_BACKWARD);
    backward.setPreferredSize(new Dimension(20, 20));
    clusterFilterNavigation = new JTextField("", 5);
    clusterFilterNavigation.setHorizontalAlignment(JTextField.CENTER);
    clusterFilterNavigation.setFont(new Font("Arial", 0, 14));
    // navigationField.setEditable(false);
    clusterFilterNavigation.setBackground(BACKGROUND_COLOR);
    clusterFilterNavigation.addFocusListener(new FocusListener() {
      @Override
      public void focusGained(FocusEvent focusevent) {}

      @Override
      public void focusLost(FocusEvent fe) {
        try {
          int trav = Integer.valueOf(((JTextField) fe.getSource()).getText().split("[\\s]+")[0])
                            .intValue()
                     - 1;
          if (trav >= 0 && trav < clusterFilterCollection.getSize(getMarkerName())) {
            // currentClusterFilter = (byte) trav;
            setCurrentClusterFilter((byte) trav);
          }
        } catch (NumberFormatException nfe) {
        }
        // clusterFilterNavigation.setText((currentClusterFilter+1)+" of
        // "+clusterFilterCollection.getSize(getMarkerName()));
        // displayClusterFilterIndex();
        // seletedScatterPanel.setPointsGeneratable(true);
        // seletedScatterPanel.setQcPanelUpdatable(true);
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
          scatterPanel.setUpdateQCPanel(true);
        }
        // scatPanel.generateRectangles();
        updateGUI();
        displayClusterFilterIndex();
      }
    });

    JButton forward = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif"));
    forward.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif"));
    forward.addActionListener(this);
    forward.setActionCommand(CLUSTER_FILTER_FORWARD);
    forward.setPreferredSize(new Dimension(20, 20));
    // JButton last1 = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
    // last1.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
    // last1.addActionListener(this);
    // last1.setActionCommand(LAST);
    // last1.setPreferredSize(new Dimension(20, 20));
    JButton delete = new JButton(Grafik.getImageIcon("images/delete9.png"));
    delete.setDisabledIcon(Grafik.getImageIcon("images/delete8.png"));
    delete.addActionListener(this);
    delete.setActionCommand(CLUSTER_FILTER_DELETE);
    delete.setPreferredSize(new Dimension(20, 20));
    // clusterFilterPanel.add(first1);
    clusterFilterPanel.add(backward);
    clusterFilterPanel.add(clusterFilterNavigation);
    clusterFilterPanel.add(forward);
    // clusterFilterPanel.add(last1);
    clusterFilterPanel.add(delete);

    newGenotype = new JComboBox(GENOTYPE_OPTIONS);
    ActionListener newGenotypeListener = new ActionListener() {
      @Override
      @SuppressWarnings("unchecked")
      public void actionPerformed(ActionEvent e) {
        byte newGenotypeSelected;
        String newGenoSelected;

        newGenoSelected = (String) ((JComboBox) e.getSource()).getSelectedItem();
        newGenotypeSelected = -2;
        for (int i = 0; i < GENOTYPE_OPTIONS.length; i++) {
          if (newGenoSelected.equals(GENOTYPE_OPTIONS[i])) {
            newGenotypeSelected = (byte) (i - 1);
            break;
          }
        }
        if ((clusterFilterCollection.getGenotype(getMarkerName(),
                                                 currentClusterFilter) != newGenotypeSelected)) {
          updateCurrentClusterFilterGenotype(newGenotypeSelected, false);
        }
      }
    };
    newGenotype.addActionListener(newGenotypeListener);
    // newGenotype.setActionCommand(NEW_GENOTYPE);//???
    clusterFilterPanel.add(newGenotype);
    // typePanel.add(newGenotype, gbc);

    clusterFilterPanel.setBackground(BACKGROUND_COLOR);

    return clusterFilterPanel;
  }

  // private JPanel nstageStdevSliderPanel() {
  // JPanel pcSliderPanel = new JPanel();
  //// pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
  // pcSliderPanel.setLayout(new GridLayout(2, 1));
  // pcSliderPanel.setBackground(BACKGROUND_COLOR);
  //
  // JSlider slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
  // // slider.setSize(new Dimension(150, 20));
  // slider.setBackground(BACKGROUND_COLOR);
  // slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 0);
  // slider.setValue(0);
  // slider.setBackground(BACKGROUND_COLOR);
  // nStageStDevLabel = new JLabel("Standard deviation filter "+stdevFilter, JLabel.CENTER);
  // nStageStDevLabel.setFont(new Font("Arial", Font.PLAIN, 16));
  // // tabPanel.add(gcLabel, gbc);
  // pcSliderPanel.add(nStageStDevLabel);
  //
  // slider.addChangeListener(new ChangeListener() {
  // public void stateChanged(ChangeEvent ce) {
  // JSlider slider = (JSlider) ce.getSource();
  // stdevFilter = (double)slider.getValue()/25;
  // nStageStDevLabel.setText("Alt. Genotype StDev Filter: "+stdevFilter);
  //// seletedScatterPanel.setPointsGeneratable(true);
  //// seletedScatterPanel.setQcPanelUpdatable(true);
  //// seletedScatterPanel.paintAgain();
  // for (int i = 0; i < scatterPanels.length; i++) {
  // scatterPanels[i].setPointsGeneratable(true);
  // scatterPanels[i].setQcPanelUpdatable(true);
  // scatterPanels[i].paintAgain();
  // }
  // // qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
  // }
  // });
  // // tabPanel.add(slider, gbc);
  // pcSliderPanel.add(slider);
  //
  // return pcSliderPanel;
  // }

  // private JPanel nstageCorrectionRatio() {
  // JPanel pcSliderPanel = new JPanel();
  //// pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
  // pcSliderPanel.setLayout(new GridLayout(2, 1));
  // pcSliderPanel.setBackground(BACKGROUND_COLOR);
  //
  // JSlider slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
  // // slider.setSize(new Dimension(150, 20));
  // slider.setBackground(BACKGROUND_COLOR);
  // slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 0);
  // slider.setValue(0);
  // slider.setBackground(BACKGROUND_COLOR);
  // correctionRatioLabel = new JLabel("Correction Ratio: " + correctionRatio, JLabel.CENTER);
  // correctionRatioLabel.setFont(new Font("Arial", Font.PLAIN, 16));
  // String usage = "This limits the ratio of PCs to the number of samples in a cluster...";
  //
  // correctionRatioLabel.setToolTipText(usage);
  // // tabPanel.add(gcLabel, gbc);
  // pcSliderPanel.add(correctionRatioLabel);
  //
  // slider.addChangeListener(new ChangeListener() {
  // public void stateChanged(ChangeEvent ce) {
  // JSlider slider = (JSlider) ce.getSource();
  // correctionRatio = (double) slider.getValue() / 100;
  // correctionRatioLabel.setText("Correction Ratio " + correctionRatio);
  //// seletedScatterPanel.setPointsGeneratable(true);
  //// seletedScatterPanel.setQcPanelUpdatable(true);
  //// seletedScatterPanel.paintAgain();
  // for (int i = 0; i < scatterPanels.length; i++) {
  // scatterPanels[i].setPointsGeneratable(true);
  // scatterPanels[i].setQcPanelUpdatable(true);
  // scatterPanels[i].paintAgain();
  // }
  // // qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
  // }
  // });
  // // tabPanel.add(slider, gbc);
  // pcSliderPanel.add(slider);
  //
  // return pcSliderPanel;
  // }

  private JPanel controlPanel() {
    JPanel controlPanel;

    controlPanel = new JPanel();
    // controlPanel.setLayout(new GridBagLayout());
    // controlPanel.setLayout(new GridLayout(0, 1, 5, 3));
    controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.Y_AXIS));
    controlPanel.setSize(50, 100);
    controlPanel.setBackground(BACKGROUND_COLOR);

    // GridBagConstraints gbc = new GridBagConstraints();
    // gbc.insets = new Insets(1,3,0,30);
    // gbc.weightx = 1.0;
    // gbc.fill = GridBagConstraints.HORIZONTAL;
    // gbc.gridwidth = GridBagConstraints.REMAINDER;
    JPanel plotPanel = new JPanel(new GridLayout(1, 2));
    plotPanel.add(plotTypePanel()); // BorderLayout.EAST);//, gbc);


    JPanel sliderPanel = new JPanel();
    sliderPanel.setLayout(new GridLayout(0, 1));
    sliderPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED),
                                                           "Data Controls", TitledBorder.LEADING,
                                                           TitledBorder.TOP));
    sliderPanel.setBackground(BACKGROUND_COLOR);
    // sliderPanel.add(createSepPanel("Data Controls"));
    sliderPanel.add(sizeSliderPanel());// , gbc);
    sliderPanel.add(gcSliderPanel());// , gbc);
    sliderPanel.add(pcSliderPanel());// , gbc);
    // sliderPanel.add(nstageStdevSliderPanel());//, gbc);
    // sliderPanel.add(nstageCorrectionRatio());//, gbc);

    Font checkBoxFont = new Font("Arial", 0, 14);

    symmetryBox = new JCheckBox("Symmetric axes");
    symmetryBox.setFont(checkBoxFont);
    symmetryBox.setBackground(BACKGROUND_COLOR);
    symmetryBox.setActionCommand(SYMMETRY);
    symmetryBox.addActionListener(this);

    // ItemListener excludeListener = new ItemListener() {
    // public void itemStateChanged(ItemEvent ie) {
    // if (ie.getStateChange() == ItemEvent.SELECTED) {
    // String filename;
    // String[] exclude;
    // Vector<String> newMarkerList;
    // int index;
    //
    // filename = proj.FILTERED_MARKERS_FILENAME.getValue();
    // if (!new File(filename).exists()) {
    // JOptionPane.showOptionDialog(null, "'Exclude Markers' is not activated due to the following
    // file not found:\n " + filename, "File Not Found", JOptionPane.CANCEL_OPTION,
    // JOptionPane.ERROR_MESSAGE, null, new String[] {"Cancel"}, "Cancel");
    // } else {
    // exclude = HashVec.loadFileToStringArray(filename, false, null, false);
    // newMarkerList = new Vector<String>();
    // index = 0;
    // for (int i = 0; i < markerList.length; i++) {
    // for (int j = 0; j < exclude.length; j++) {
    // if (markerList[i].equalsIgnoreCase(exclude[j])) {
    // newMarkerList.add(exclude[j]);
    // if (i == markerIndex) {
    // index = newMarkerList.size() - 1;
    // }
    // break;
    // }
    // }
    // }
    // changeToNewMarkerList(newMarkerList.toArray(new String[0]), index);
    // updateGUI();
    // }
    // } else {
    // int index = markerIndexBak;
    // for (int i = 0; i < masterMarkerList.length; i++) {
    // if(masterMarkerList[i].equalsIgnoreCase(markerList[markerIndex])) {
    // index = i;
    // break;
    // }
    // }
    // revertMarkersToOriginalList(index);
    // updateGUI();
    // }
    // }
    // };
    // excludeMarkerBox = new JCheckBox("<html>Hide Excluded <br />Markers</html>");
    // excludeMarkerBox.setFont(new Font("Arial", 0, 14));
    // excludeMarkerBox.addItemListener(excludeListener);
    // excludeMarkerBox.setBackground(BACKGROUND_COLOR);

    excludeSampleBox = new JCheckBox("<html>Hide Excluded <br />Samples</html>");
    excludeSampleBox.setFont(checkBoxFont);
    excludeSampleBox.setBackground(BACKGROUND_COLOR);
    excludeSampleBox.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent e) {
        excludeSamples[selectedPanelIndex] = e.getStateChange() == ItemEvent.SELECTED;
        updateGUI();
      }
    });
    excludeSampleBox.setSelected(true);

    correctionBox = new JCheckBox("Correct Data");
    correctionBox.setToolTipText("The keyboard shortcut for this function is ctrl+D");
    correctionBox.setFont(new Font("Arial", 0, 14));
    correctionBox.setBackground(BACKGROUND_COLOR);
    correctionBox.setActionCommand(CORRECTION);
    correctionBox.addActionListener(this);


    maskMissingBox = new JCheckBox("Mask Missing");
    // maskMissingBox.setToolTipText("The keyboard shortcut for this function is ctrl+D");
    maskMissingBox.setFont(checkBoxFont);
    // maskMissingBox.addItemListener(correctionListener);
    maskMissingBox.setBackground(BACKGROUND_COLOR);
    maskMissingBox.setActionCommand(MASK_MISSING);
    maskMissingBox.addActionListener(this);

    mendelianErrorBox = new JCheckBox("<html>Display Mendelian <br />Errors</html>");
    mendelianErrorBox.setFont(checkBoxFont);
    mendelianErrorBox.setBackground(BACKGROUND_COLOR);
    mendelianErrorBox.setActionCommand(MENDELIAN_ERROR);
    mendelianErrorBox.addActionListener(this);
    if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
      mendelianErrorBox.setEnabled(false);
    }

    // JPanel boxPanel = new JPanel();
    // boxPanel.setLayout(new GridLayout(3, 2));
    // boxPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED),
    // "View Controls", TitledBorder.LEADING, TitledBorder.TOP));
    // boxPanel.setBackground(BACKGROUND_COLOR);
    // boxPanel.add(symmetryBox);//, gbc);
    // boxPanel.add(correctionBox);//, gbc);
    // boxPanel.add(excludeMarkerBox);//, gbc);
    // boxPanel.add(excludeSampleBox);//, gbc);
    // boxPanel.add(maskMissingBox);//, gbc);
    // controlPanel.add(boxPanel);
    JPanel boxPanel = new JPanel();
    boxPanel.setLayout(new BoxLayout(boxPanel, BoxLayout.Y_AXIS));
    boxPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED),
                                                        "View Controls", TitledBorder.LEADING,
                                                        TitledBorder.TOP));
    boxPanel.setBackground(BACKGROUND_COLOR);
    boxPanel.add(symmetryBox);// , gbc);
    boxPanel.add(correctionBox);// , gbc);
    // boxPanel.add(excludeMarkerBox);//, gbc);
    boxPanel.add(excludeSampleBox);// , gbc);
    boxPanel.add(maskMissingBox);// , gbc);
    boxPanel.add(mendelianErrorBox);
    plotPanel.add(boxPanel, BorderLayout.WEST);
    controlPanel.add(plotPanel);

    controlPanel.add(sliderPanel);

    typeRadioButtons[0].setSelected(true);

    return controlPanel;
  }

  private void convertSamples() {
    sampleFIDIIDs = new String[samples.length];
    for (int i = 0; i < sampleFIDIIDs.length; i++) {
      String[] look = sampleData.lookup(samples[i]);
      if (look == null) {
        log.reportError("Error - no lookup data found for sample {" + samples[i] + "}");
        sampleFIDIIDs[i] = samples[i];
      } else {
        sampleFIDIIDs[i] = sampleData.lookup(samples[i])[1];
      }
    }
  }

  private JMenuBar createJMenuBar() {
    JMenuBar menuBar = new JMenuBar();

    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic(KeyEvent.VK_F);


    JMenuItem newListItem = new JMenuItem();
    newListItem.setMnemonic(KeyEvent.VK_N);
    newListItem.setActionCommand(NEW_LIST_COMMAND);
    newListItem.addActionListener(this);
    newListItem.setText("New Marker List");
    fileMenu.add(newListItem);

    JMenuItem saveItem = new JMenuItem();
    saveItem.setMnemonic(KeyEvent.VK_S);
    saveItem.setActionCommand(SAVE_LIST_COMMAND);
    saveItem.addActionListener(this);
    saveItem.setText("Save Current Marker List");
    fileMenu.add(saveItem);

    JMenuItem loadItem = new JMenuItem();
    loadItem.setMnemonic(KeyEvent.VK_L);
    loadItem.setActionCommand(LOAD_LIST_COMMAND);
    loadItem.addActionListener(this);
    loadItem.setText("Load Marker List");
    fileMenu.add(loadItem);

    // JMenuItem testItem = new JMenuItem();
    // testItem.setActionCommand(LOAD_LIST_COMMAND + "test");
    // testItem.addActionListener(this);
    // testItem.setText("Test Config");
    // fileMenu.add(testItem);

    final JMenu previousListItem = new JMenu();
    previousListItem.setMnemonic(KeyEvent.VK_P);
    previousListItem.setText("Marker Lists");
    fileMenu.add(previousListItem);

    final ActionListener loadFileListener = new ActionListener() {
      @Override
      public void actionPerformed(final ActionEvent e) {
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            String cmd = e.getActionCommand();
            loadMarkerFile(cmd);
          }
        });
      }
    };

    String[] vals = proj.DISPLAY_MARKERS_FILENAMES.getValue();
    HashSet<String> mnemonics = new HashSet<String>();
    for (String val : vals) {
      JMenuItem listEntry = new JMenuItem();
      String[] tmp = val.split("/");
      String mnemChar = tmp[tmp.length - 1].charAt(0) + "";
      int cnt = 1;
      while (mnemonics.contains(mnemChar) && cnt < tmp[tmp.length - 1].length()) {
        mnemChar = tmp[tmp.length - 1].charAt(cnt) + "";
        cnt++;
      }
      if (!mnemonics.contains(mnemChar)) {
        listEntry.setMnemonic(mnemChar.charAt(0));
      }
      listEntry.addActionListener(loadFileListener);
      listEntry.setActionCommand(val);
      listEntry.setText(val);
      previousListItem.add(listEntry);
    }

    JMenu sortItem = new JMenu();
    sortItem.setMnemonic(KeyEvent.VK_S);
    sortItem.setText("Sort Previous Lists");
    fileMenu.add(sortItem);

    ActionListener sortListener = new ActionListener() {
      String prevSort = "list";

      @Override
      public void actionPerformed(ActionEvent ae) {
        if (ae.getActionCommand().equals(prevSort)) {
          return;
        }
        prevSort = ae.getActionCommand();
        previousListItem.removeAll();
        String[] vals = proj.DISPLAY_MARKERS_FILENAMES.getValue();
        HashSet<String> mnemonics = new HashSet<String>();
        if ("alpha".equals(ae.getActionCommand())) {
          int[] indices = Sort.quicksort(vals);
          for (int indice : indices) {
            JMenuItem listEntry = new JMenuItem();
            String[] tmp = vals[indice].split("/");
            String mnemChar = tmp[tmp.length - 1].charAt(0) + "";
            int cnt = 1;
            while (mnemonics.contains(mnemChar) && cnt < tmp[tmp.length - 1].length()) {
              mnemChar = tmp[tmp.length - 1].charAt(cnt) + "";
              cnt++;
            }
            if (!mnemonics.contains(mnemChar)) {
              listEntry.setMnemonic(mnemChar.charAt(0));
            }
            listEntry.addActionListener(loadFileListener);
            listEntry.setActionCommand(vals[indice]);
            listEntry.setText(vals[indice]);
            previousListItem.add(listEntry);
          }
        } else if ("time".equals(ae.getActionCommand())) {
          long[] times = new long[vals.length];
          for (int i = 0; i < vals.length; i++) {
            File f = new File(vals[i]);
            times[i] = f.lastModified();
          }
          int[] indices = Sort.quicksort(times);
          for (int i = indices.length - 1; i >= 0; i--) {
            JMenuItem listEntry = new JMenuItem();
            String[] tmp = vals[indices[i]].split("/");
            String mnemChar = tmp[tmp.length - 1].charAt(0) + "";
            int cnt = 1;
            while (mnemonics.contains(mnemChar) && cnt < tmp[tmp.length - 1].length()) {
              mnemChar = tmp[tmp.length - 1].charAt(cnt) + "";
              cnt++;
            }
            if (!mnemonics.contains(mnemChar)) {
              listEntry.setMnemonic(mnemChar.charAt(0));
            }
            listEntry.addActionListener(loadFileListener);
            listEntry.setActionCommand(vals[indices[i]]);
            listEntry.setText(vals[indices[i]]);
            previousListItem.add(listEntry);
          }
        } else /* if ("list".equals(ae.getActionCommand())) */ {
          for (String val : vals) {
            JMenuItem listEntry = new JMenuItem();
            String[] tmp = val.split("/");
            String mnemChar = tmp[tmp.length - 1].charAt(0) + "";
            int cnt = 1;
            while (mnemonics.contains(mnemChar) && cnt < tmp[tmp.length - 1].length()) {
              mnemChar = tmp[tmp.length - 1].charAt(cnt) + "";
              cnt++;
            }
            if (!mnemonics.contains(mnemChar)) {
              listEntry.setMnemonic(mnemChar.charAt(0));
            }
            listEntry.addActionListener(loadFileListener);
            listEntry.setActionCommand(val);
            listEntry.setText(val);
            previousListItem.add(listEntry);
          }
        }
      }
    };

    JRadioButtonMenuItem sortListOrder = new JRadioButtonMenuItem();
    sortListOrder.setMnemonic(KeyEvent.VK_O);
    sortListOrder.addActionListener(sortListener);
    sortListOrder.setActionCommand("list");
    sortListOrder.setText("Listed Order");
    sortItem.add(sortListOrder);

    JRadioButtonMenuItem sortAlphaNum = new JRadioButtonMenuItem();
    sortAlphaNum.setMnemonic(KeyEvent.VK_A);
    sortAlphaNum.addActionListener(sortListener);
    sortAlphaNum.setActionCommand("alpha");
    sortAlphaNum.setText("Alphanumeric");
    sortItem.add(sortAlphaNum);

    JRadioButtonMenuItem sortTimestamp = new JRadioButtonMenuItem();
    sortTimestamp.setMnemonic(KeyEvent.VK_T);
    sortTimestamp.addActionListener(sortListener);
    sortTimestamp.setActionCommand("time");
    sortTimestamp.setText("Timestamp");
    sortItem.add(sortTimestamp);

    ButtonGroup bg = new ButtonGroup();
    bg.add(sortListOrder);
    bg.add(sortAlphaNum);
    bg.add(sortTimestamp);
    sortListOrder.setSelected(true);
    // sorting will rearrange items on previousListItem menu; likely have to be either final or
    // field

    JMenu delFileMenu = new JMenu();
    delFileMenu.setMnemonic(KeyEvent.VK_D);
    delFileMenu.setText("Delete Previous List");
    fileMenu.add(delFileMenu);

    ActionListener deleteFileListener = new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        String toDelete = e.getActionCommand();
        int opt = JOptionPane.showConfirmDialog(ScatterPlot.this, "Delete file on disk?",
                                                "Delete Marker List File?",
                                                JOptionPane.YES_NO_CANCEL_OPTION);
        boolean deleteFile = false;
        if (opt == JOptionPane.YES_OPTION) {
          deleteFile = true;
        } else if (opt == JOptionPane.NO_OPTION) {
          deleteFile = false;
        } else {
          return; // cancelled
        }
        // String[] vals = proj.DISPLAY_MARKERS_FILENAMES.getValue();
        // String[] newVals = new String[vals.length - 1];
        // int cnt = 0;
        // for (String val : vals) {
        // if (!toDelete.equals(val)) {
        // newVals[cnt] = val;
        // cnt++;
        // }
        // }
        if (deleteFile) {
          (new File(toDelete)).delete();
        }
        proj.DISPLAY_MARKERS_FILENAMES.removeValue(toDelete);
        ScatterPlot.this.setJMenuBar(ScatterPlot.this.createJMenuBar());
        ScatterPlot.this.validate();
        ScatterPlot.this.repaint();
      }
    };

    mnemonics = new HashSet<String>();
    for (String val : vals) {
      JMenuItem listEntry = new JMenuItem();
      String[] tmp = val.split("/");
      String mnemChar = tmp[tmp.length - 1].charAt(0) + "";
      int cnt = 1;
      while (mnemonics.contains(mnemChar) && cnt < tmp[tmp.length - 1].length()) {
        mnemChar = tmp[tmp.length - 1].charAt(cnt) + "";
        cnt++;
      }
      if (!mnemonics.contains(mnemChar)) {
        listEntry.setMnemonic(mnemChar.charAt(0));
      }
      listEntry.addActionListener(deleteFileListener);
      listEntry.setActionCommand(val);
      listEntry.setText(val);
      delFileMenu.add(listEntry);
    }


    JMenu actMenu = new JMenu("Actions");
    actMenu.setMnemonic(KeyEvent.VK_A);

    JMenuItem screenCapItem = new JMenuItem();
    screenCapItem.setActionCommand(CAPTURE);
    screenCapItem.addActionListener(this);
    screenCapItem.setText("Screen Capture");
    screenCapItem.setMnemonic(KeyEvent.VK_S);
    actMenu.add(screenCapItem);

    JMenuItem dumpDataItem = new JMenuItem();
    dumpDataItem.setActionCommand(DUMP);
    dumpDataItem.addActionListener(this);
    dumpDataItem.setText(DUMP);
    dumpDataItem.setMnemonic(KeyEvent.VK_D);
    actMenu.add(dumpDataItem);

    // Mask Missing Values
    // JCheckBoxMenuItem maskMissingItem = new JCheckBoxMenuItem();
    // maskMissingItem.setActionCommand(MASK_MISSING);
    // maskMissingItem.addActionListener(this);
    // maskMissingItem.setText(MASK_MISSING);
    // maskMissingItem.setMnemonic(KeyEvent.VK_M);
    // actMenu.add(maskMissingItem);

    menuBar.add(fileMenu);
    menuBar.add(actMenu);

    return menuBar;
  }

  private void deactivateAllAnnotationMaps() {
    for (char annotationKey : annotationKeys) {
      removeAnnotationFromMaps(annotationKey);
    }
  }

  private void displayBlast() {
    if (blastFrame == null) {
      final ChangeListener alignLengthListener = new ChangeListener() {
        @Override
        public void stateChanged(ChangeEvent e) {
          int value = ((SpinnerNumberModel) ((JSpinner) e.getSource()).getModel()).getNumber()
                                                                                  .intValue();
          blastFrame.updateAnnotations(value - 1);
          blastFrame.updateLabels();
          setHistogram(value);
          updateBLASTPanel();
          SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
              blastFrame.setVisible(true);
              blastFrame.setSecondPanel(histFrame.getPanel());
            }
          });
        }
      };
      blastFrame = new BlastFrame(proj, alignLengthListener);
    }
    blastFrame.setAnnotations(blastResults[markerIndex], referenceGenome);


    if (histFrame == null) {
      histFrame = TwoDPlot.createGUI(getProject(), false, false);
    }

    double filter = proj.BLAST_PROPORTION_MATCH_FILTER.getValue();
    int probe = proj.ARRAY_TYPE.getValue().getProbeLength();
    int alignFilter = blastFrame == null ? (int) (filter * probe) : blastFrame.currentAlignFilter;
    setHistogram(alignFilter);
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        blastFrame.setVisible(true);
        blastFrame.setSecondPanel(histFrame.getPanel());
      }
    });
  }

  private void doDataDump() {
    int count = 1;
    String filename = null;
    do {
      filename = markerList[markerIndex] + "_dump" + (count == 1 ? "" : " v" + count);
      filename = ext.replaceWithLinuxSafeCharacters(filename, true);
      count++;
    } while (new File(proj.PROJECT_DIRECTORY.getValue() + filename + ".xln").exists());
    String outFile = proj.PROJECT_DIRECTORY.getValue() + filename + ".xln";
    getCurrentMarkerData().dump(sampleData, outFile, samples, false, log);
    log.report(ext.getTime() + "]\tWrote data dump to [" + outFile + "]");
  }

  private void doPanelSelection(final int newIndex, final boolean swap) {
    if (swapping) {
      return;
    } else {
      swapping = true;
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          if (selectedPanelIndex == newIndex) {
            if (swap) {
              swapView();
            }
          } else {
            if (!showingAll && swap) {
              swapView();
            }
            selectedPanelIndex = newIndex;
            typeRadioButtons[plot_types[newIndex]].setSelected(true);
            maskMissingBox.setSelected(maskMissing(newIndex));
            excludeSampleBox.setSelected(hideExcludedSamples(newIndex));
            symmetryBox.setSelected(symmetricAxes(newIndex));
            mendelianErrorBox.setSelected(displayMendelianErrors(newIndex));
            correctionBox.setSelected(getCorrection(newIndex));
            colorKeyPanel.setSisterPanel(scatterPanels[newIndex]);
            if (scatterPanels[newIndex].uniqueValueCounts != null) {
              int cls = classes[newIndex];
              colorKeyPanel.setCurrentClass(classes[newIndex]);
              colorKeyPanel.updateColorKey(scatterPanels[newIndex].uniqueValueCounts.convertToHash());
              colorKeyPanel.updateColorVariablePanel(classListener);
              classes[newIndex] = cls;
              colorKeyPanel.setCurrentClass(classes[newIndex]);
            }
          }
          updateGUI();
          swapping = false;
        }
      });
    }
  }

  private void doScreenCapture() {
    int count = 1;
    String filename = null;
    // do {
    // filename = markerList[markerIndex]+"_" + MarkerData.TYPES[plot_type][0] + "-" +
    // MarkerData.TYPES[plot_type][1] + "_" + (count == 1 ? "" : "v" + count);
    // filename = ext.replaceWithLinuxSafeCharacters(filename, true);
    // count++;
    // } while (new File(proj.PROJECT_DIRECTORY.getValue()+filename+".png").exists());
    // for (int i = 0; i < scatterPanels.length; i++) {
    // scatterPanels[i].screenCapture(proj.PROJECT_DIRECTORY.getValue()+filename+"_panel" + i +
    // ".png");
    // log.report(ext.getTime() + "]\tWrote screen capture to
    // ["+proj.PROJECT_DIRECTORY.getValue()+filename+"_panel" + i + ".png]");
    // }
    for (int i = 0; i < scatterPanels.length; i++) {
      do {
        filename = markerList[markerIndex] + "_" + ScatterPlot.TYPES[plot_types[i]][0] + "-"
                   + ScatterPlot.TYPES[plot_types[i]][1] + "_" + (count == 1 ? "" : "v" + count);
        filename = ext.replaceWithLinuxSafeCharacters(filename, true);
        count++;
      } while (new File(proj.PROJECT_DIRECTORY.getValue() + filename + ".png").exists());
      scatterPanels[i].screenCapture(proj.PROJECT_DIRECTORY.getValue() + filename + "_panel" + i
                                     + ".png");
      log.report(ext.getTime() + "]\tWrote screen capture to [" + proj.PROJECT_DIRECTORY.getValue()
                 + filename + "_panel" + i + ".png]");
    }
  }

  private JPanel eastPanel() {
    JPanel eastPanel;
    JTabbedPane tabbedPane;

    tabbedPane = new JTabbedPane();
    tabbedPane.setBackground(BACKGROUND_COLOR);

    tabbedPane.addTab("Control", null, controlPanel(), "Manipulate the chart");
    tabbedPane.setMnemonicAt(0, KeyEvent.VK_1);

    tabbedPane.addTab("Centroid", null, centroidPanel(), "Displays the centroid");
    tabbedPane.setMnemonicAt(1, KeyEvent.VK_2);

    tabbedPane.addTab("GC Adjust", null, gcCorrectionPanel(), "Adjust LRR by gc content");
    // tabbedPane.setMnemonicAt(1, KeyEvent.VK_2);

    annotationPanel();
    tabbedPane.addTab("Annotation", null, annotationScrollPane, "Annotation to this marker");
    tabbedPane.setMnemonicAt(2, KeyEvent.VK_3);

    eastPanel = new JPanel();
    eastPanel.setLayout(new MigLayout("", "[grow]", "[grow,fill][fill][fill][fill]"));
    eastPanel.setBackground(BACKGROUND_COLOR);
    eastPanel.add(tabbedPane, "cell 0 0, grow");

    JTabbedPane qcTabbedPanel = new JTabbedPane();
    qcTabbedPanel.setBackground(BACKGROUND_COLOR);
    qcPanel = new JPanel();
    qcPanel.setLayout(new MigLayout("hidemode 0", "[][grow]", "[][][][][][][]"));
    qcPanel.setBackground(BACKGROUND_COLOR);
    qcTabbedPanel.addTab("QC Metrics", null, qcPanel, "Quality Control Metrics");
    eastPanel.add(qcTabbedPanel, "cell 0 1, grow");

    JTabbedPane blastTabbedPanel = new JTabbedPane();
    blastTabbedPanel.setBackground(BACKGROUND_COLOR);
    blastPanel = new JPanel();
    blastPanel.setLayout(new MigLayout("hidemode 0", "[][grow]", ""));
    blastPanel.setBackground(BACKGROUND_COLOR);
    blastTabbedPanel.addTab("BLAST Metrics", null, blastPanel,
                            "Basic Local Alignment Search Tool Metrics");
    eastPanel.add(blastTabbedPanel, "cell 0 2, grow");

    JTabbedPane cfTabbedPanel = new JTabbedPane();
    cfTabbedPanel.setBackground(BACKGROUND_COLOR);
    JPanel cfPanel = new JPanel();
    cfPanel.setLayout(new GridLayout(0, 1));
    cfPanel.setBackground(BACKGROUND_COLOR);
    cfPanel.add(clusterFilterPanel());
    cfTabbedPanel.addTab("Cluster Filters", null, cfPanel, "Apply Cluster Filters");
    eastPanel.add(cfTabbedPanel, "cell 0 3, grow");

    return eastPanel;
  }

  private JPanel gcCorrectionPanel() {
    JPanel gcAdjustorPanel;

    gcAdjustorPanel = new JPanel();
    // tabPanel.setLayout(new BoxLayout(tabPanel, BoxLayout.Y_AXIS));
    gcAdjustorPanel.setLayout(new GridBagLayout());
    gcAdjustorPanel.setSize(50, 100);
    gcAdjustorPanel.setBackground(BACKGROUND_COLOR);

    GridBagConstraints gbc = new GridBagConstraints();
    gbc.insets = new Insets(1, 3, 0, 30);
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridwidth = GridBagConstraints.REMAINDER;

    ItemListener centListener = new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent ie) {
        int index = ext.indexOfStr(((JCheckBox) ie.getSource()).getText(), gcCList[0]);
        displaygcAdjustor[index] = ((JCheckBox) ie.getSource()).isSelected();

        gcAdjustorLabels[index].setVisible(displaygcAdjustor[index]);
        for (int i = 0; i < gcAdjustorLabels.length; i++) {
          if (i != index) {
            displaygcAdjustor[i] = false;
            gcAdjustorBoxes[i].setSelected(false);
          }
        }
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
        }
        updateGUI();
      }
    };
    JLabel label = new JLabel("  ");
    label.setFont(new Font("Arial", 0, 20));
    gcAdjustorPanel.add(label, gbc);
    label = new JLabel("GC Adjustment (LRR/BAF only) :");
    label.setToolTipText("These adjustments will be applied to LRR values (and BAFs if new centroids were stored");
    label.setFont(new Font("Arial", 0, 20));
    label.setHorizontalAlignment(JLabel.CENTER);
    gcAdjustorPanel.add(label, gbc);
    gcAdjustorBoxes = new JCheckBox[gcCList[0].length];
    displaygcAdjustor = new boolean[gcCList[0].length];
    gcAdjustorLabels = new JLabel[gcCList[0].length];
    for (int i = 0; i < gcCList[0].length; i++) {
      gcAdjustorBoxes[i] = new JCheckBox(gcCList[0][i]);
      gcAdjustorBoxes[i].setFont(new Font("Arial", 0, 14));
      gcAdjustorBoxes[i].setSelected(displaygcAdjustor[i]);
      gcAdjustorBoxes[i].addItemListener(centListener);
      gcAdjustorBoxes[i].setBorder(BorderFactory.createLineBorder(ScatterPanel.DEFAULT_COLORS[5
                                                                                              + i],
                                                                  5));
      gcAdjustorBoxes[i].setBorderPainted(true);
      gcAdjustorBoxes[i].setBackground(BACKGROUND_COLOR);
      gcAdjustorPanel.add(gcAdjustorBoxes[i], gbc);

      gcAdjustorLabels[i] = new JLabel("");
      gcAdjustorLabels[i].setVisible(displaygcAdjustor[i]);
      gcAdjustorPanel.add(gcAdjustorLabels[i], gbc);
    }

    return gcAdjustorPanel;
  }

  private JPanel gcSliderPanel() {
    JPanel gcSliderPanel = new JPanel();
    // gcSliderPanel.setLayout(new BoxLayout(gcSliderPanel, BoxLayout.Y_AXIS));
    gcSliderPanel.setLayout(new GridLayout(2, 1));
    gcSliderPanel.setBackground(BACKGROUND_COLOR);

    JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
    // slider.setSize(new Dimension(150, 20));
    slider.setBackground(BACKGROUND_COLOR);
    slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
    slider.setBackground(BACKGROUND_COLOR);
    gcLabel = new JLabel("GC > 0." + DEFAULT_GC_THRESHOLD, JLabel.CENTER);
    gcLabel.setFont(new Font("Arial", Font.PLAIN, 16));
    // tabPanel.add(gcLabel, gbc);
    gcSliderPanel.add(gcLabel);

    slider.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent ce) {
        JSlider slider = (JSlider) ce.getSource();
        gcThreshold = slider.getValue() / 100f;
        gcLabel.setText("GC > " + ext.formDeci(gcThreshold, 2, true));
        // seletedScatterPanel.setPointsGeneratable(true);
        // seletedScatterPanel.setQcPanelUpdatable(true);
        // seletedScatterPanel.paintAgain();
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
          scatterPanel.setUpdateQCPanel(true);
          scatterPanel.paintAgain();
        }
        // qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
      }
    });
    // tabPanel.add(slider, gbc);
    gcSliderPanel.add(slider);

    return gcSliderPanel;
  }

  private int getAvailableMarker(boolean forwardOrBackward, boolean firstOrLastInThatDirection,
                                 boolean allMarkersOrOnlyThoseAnnotatedOrUnannotated,
                                 boolean annotatedOrUnannotated) {
    int result;
    int begin;
    int end;
    int step;


    result = markerIndex;
    if (forwardOrBackward) {
      if (firstOrLastInThatDirection) {
        begin = markerIndex + 1;
        end = markerList.length;
        step = 1;
      } else {
        begin = markerList.length - 1;
        end = markerIndex;
        step = -1;
      }
    } else {
      if (firstOrLastInThatDirection) {
        begin = markerIndex - 1;
        end = -1;
        step = -1;
      } else {
        begin = 0;
        end = markerIndex;
        step = 1;
      }
    }
    for (int i = begin; i != end; i = i + step) {
      if (indexOfAnnotationUsedAsMarkerList < 0) {
        if (allMarkersOrOnlyThoseAnnotatedOrUnannotated
            || !(annotatedOrUnannotated
                 ^ annotationCollection.markerHasAnyAnnotation(markerList[i]))) {
          result = i;
          break;
        }
      } else if (annotationCollection.markerHasAnnotation(markerList[i],
                                                          annotationKeys[indexOfAnnotationUsedAsMarkerList])) {
        result = i;
        break;
      }
    }

    return result;
  }

  private void inputMapAndActionMap(JComponent comp) {
    InputMap inputMap;
    ActionMap actionMap;

    inputMap = comp.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    actionMap = comp.getActionMap();

    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), ALT_UP);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), ALT_DOWN);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), ALT_LEFT);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), ALT_RIGHT);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_HOME, InputEvent.CTRL_MASK), FIRST);
    // inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, InputEvent.CTRL_MASK), PREVIOUS);
    // inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, InputEvent.CTRL_MASK), NEXT);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, 0), PREVIOUS);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, 0), NEXT);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_END, InputEvent.CTRL_MASK), LAST);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_MASK), SYMMETRY);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_MASK), CORRECTION);

    String SEL_1 = "SEL1";
    String SEL_2 = "SEL2";
    String SEL_3 = "SEL3";
    String SEL_4 = "SEL4";

    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_1, 0/* InputEvent.CTRL_DOWN_MASK */), SEL_1);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_2, 0/* InputEvent.CTRL_DOWN_MASK */), SEL_2);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_3, 0/* InputEvent.CTRL_DOWN_MASK */), SEL_3);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_4, 0/* InputEvent.CTRL_DOWN_MASK */), SEL_4);

    AbstractAction panelSelectionAction = new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        boolean hasFocus = true;
        if (markerName.hasFocus() || commentLabel.hasFocus() || addAnnotationField.hasFocus()
            || navigationField.hasFocus()) {
          hasFocus = false;
        }
        if (!hasFocus) {
          return;
        }
        try {
          int i = Integer.valueOf(e.getActionCommand()) - 1;
          doPanelSelection(i, true);
        } catch (NumberFormatException e2) {
          return;
        }
      }
    };

    actionMap.put(SEL_1, panelSelectionAction);
    actionMap.put(SEL_2, panelSelectionAction);
    actionMap.put(SEL_3, panelSelectionAction);
    actionMap.put(SEL_4, panelSelectionAction);
    //
    // inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_TAB, InputEvent.CTRL_MASK), "SWITCH");
    //
    // actionMap.put("SWITCH", new AbstractAction() {
    // @Override
    // public void actionPerformed(ActionEvent e) {
    // System.out.println("SWAP");
    // swapView();
    // }
    // });

    actionMap.put(ALT_UP, new CycleRadio(typeRadioButtons, -1));
    actionMap.put(ALT_DOWN, new CycleRadio(typeRadioButtons, 1));
    // actionMap.put(ALT_LEFT, new CycleRadio(colorKeyPanel.getClassRadioButtons(), -1));
    // actionMap.put(ALT_RIGHT, new CycleRadio(colorKeyPanel.getClassRadioButtons(), 1));
    actionMap.put(FIRST, new AbstractAction() {
      public static final long serialVersionUID = 4L;

      @Override
      public void actionPerformed(ActionEvent e) {
        first();
      }
    });
    actionMap.put(PREVIOUS, new AbstractAction() {
      public static final long serialVersionUID = 5L;

      @Override
      public void actionPerformed(ActionEvent e) {
        previous();
      }
    });
    actionMap.put(NEXT, new AbstractAction() {
      public static final long serialVersionUID = 6L;

      @Override
      public void actionPerformed(ActionEvent e) {
        next();
      }
    });
    actionMap.put(LAST, new AbstractAction() {
      public static final long serialVersionUID = 7L;

      @Override
      public void actionPerformed(ActionEvent e) {
        last();
      }
    });
    actionMap.put(SYMMETRY, new AbstractAction() {
      public static final long serialVersionUID = 7L;

      @Override
      public void actionPerformed(ActionEvent e) {
        symmetryBox.setSelected(!symmetryBox.isSelected());
      }
    });
    actionMap.put(EXCLUDE, new AbstractAction() {
      public static final long serialVersionUID = 7L;

      @Override
      public void actionPerformed(ActionEvent e) {
        excludeSampleBox.setSelected(!excludeSampleBox.isSelected());
      }
    });
    actionMap.put(CORRECTION, new AbstractAction() {
      public static final long serialVersionUID = 7L;

      @Override
      public void actionPerformed(ActionEvent e) {
        correctionBox.setSelected(!correctionBox.isSelected());
      }
    });
  }

  private boolean loadClusterFilterFiles() {
    String[] otherClusterFilterFiles;
    int choice;
    String[] options = new String[] {"Yes, load and delete old file", "No, just delete old file",
                                     "Cancel and close ScatterPlot"};
    String clusterFilterFilename;

    // clusterFilterFilename = proj.getFilename(proj.CLUSTER_FILTER_COLLECTION_FILENAME);
    clusterFilterFilename = proj.CLUSTER_FILTER_COLLECTION_FILENAME.getValue();
    otherClusterFilterFiles = Files.list(proj.DATA_DIRECTORY.getValue(false, true),
                                         ".tempClusterFilters.ser", jar);
    if (otherClusterFilterFiles.length > 0) {
      choice =
             JOptionPane.showOptionDialog(null,
                                          "Error - either multiple instances of ScatterPlot are running or ScatterPlot failed to close properly\n"
                                                + "last time. The ability to generate new ClusterFilters will be disabled until this file has been\n"
                                                + "removed. Do you want to load the contents of the temporary file into memory before it is deleted?",
                                          "Error", JOptionPane.YES_NO_CANCEL_OPTION,
                                          JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      if (choice == 0) {
        // load the last one in otherClusterFilerFiles[]
        clusterFilterCollection = ClusterFilterCollection.load(
                                                               proj.DATA_DIRECTORY.getValue(false,
                                                                                            true)
                                                               + otherClusterFilterFiles[otherClusterFilterFiles.length
                                                                                         - 1],
                                                               jar);
        for (String otherClusterFilterFile : otherClusterFilterFiles) {
          (new File(proj.DATA_DIRECTORY.getValue(false, true) + otherClusterFilterFile)).delete();
        }
        startAutoSaveToTempFile();
        setClusterFilterUpdated(true);
      } else if (choice == 1) {
        // load permanent
        if (Files.exists(clusterFilterFilename, jar)) {
          clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFilename, jar);
        } else {
          clusterFilterCollection = new ClusterFilterCollection();
        }
        for (String otherClusterFilterFile : otherClusterFilterFiles) {
          (new File(proj.DATA_DIRECTORY.getValue(false, true) + otherClusterFilterFile)).delete();
        }
      } else {
        return false;
      }
    } else if (Files.exists(clusterFilterFilename, jar)) {
      clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFilename, jar);
    } else {
      clusterFilterCollection = new ClusterFilterCollection();
    }

    startAutoSaveToTempFile();

    return true;
  }

  private void loadMarkerDataFromList(int newMarkerIndex) {
    markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
    // markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSameThread(proj, markerList);

    // String[] plinkRoots = proj.PLINK_DIR_FILEROOTS.getValue();
    // String plinkDirFileRoot = proj.DATA_DIRECTORY.getValue() + proj.PLINK_FILEROOT.getValue();
    // plinkGenotypes = (new PlinkMarkerLoader()).run(plinkDirFileRoot, markerList);
    for (String fileRoot : plinkMarkerLoaders.keySet()) {
      plinkMarkerLoaders.put(fileRoot,
                             PlinkMarkerLoader.loadPlinkDataFromListInSeparateThread(proj, fileRoot,
                                                                                     markerList,
                                                                                     sampleFIDIIDs));
    }

    markerIndex = newMarkerIndex;
    updateMarkerIndexHistory();
    previousMarkerIndex = -1;
    // navigationField.getActionListeners()[0].actionPerformed(e);
  }

  private void loadMarkerFile(String file) {
    if (Files.exists(file)) {
      proj.getLog().report("Loading marker list file: " + file);

      loadMarkerListFromFile(file);
      if (masterMarkerList.length == 0) {
        JOptionPane.showMessageDialog(null,
                                      "Error - file '" + file + "' was devoid of any valid markers",
                                      "Error", JOptionPane.ERROR_MESSAGE);
        fail = true;
        return;
      }
      if (masterCommentList == null) {
        masterCommentList = Array.stringArray(masterMarkerList.length, "");
      }
      proj.getProgressMonitor().beginDeterminateTask("SCATTERPLOT",
                                                     "Loading ScatterPlot Marker List...", 10,
                                                     ProgressMonitor.DISPLAY_MODE.GUI_ONLY);
      resetAfterLoad();
      proj.getProgressMonitor().endTask("SCATTERPLOT");
      ScatterPlot.this.setJMenuBar(ScatterPlot.this.createJMenuBar());
      ScatterPlot.this.validate();
      finishProcessing();
    } else {
      proj.getLog().reportError("Error - file " + file + " not found");
    }
  }

  private PrincipalComponentsResiduals loadPcResids() {
    // String pcFile = proj.getFilename(proj.INTENSITY_PC_FILENAME);
    String pcFile = proj.INTENSITY_PC_FILENAME.getValue();
    PrincipalComponentsResiduals pcResids;
    if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + ext.removeDirectoryInfo(pcFile))) {
      proj.getLog().report("Info - loading " + ext.removeDirectoryInfo(pcFile));
      // pcResids = new PrincipalComponentsResiduals(proj, ext.removeDirectoryInfo(pcFile), null,
      // Integer.parseInt(proj.getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), false, 0, false,
      // false, null);
      pcResids = new PrincipalComponentsResiduals(proj, ext.removeDirectoryInfo(pcFile), null,
                                                  proj.getProperty(proj.INTENSITY_PC_NUM_COMPONENTS),
                                                  false, 0, false, false, null);
      setNumComponents(Math.min(numComponents, pcResids.getTotalNumComponents()));
    } else {
      proj.getLog().report("Info - did not find " + proj.PROJECT_DIRECTORY.getValue()
                           + ext.removeDirectoryInfo(pcFile));
      pcResids = null;
    }
    return pcResids;
  }

  private JPanel markerPanel() {
    JPanel descrPanel = new JPanel();
    descrPanel.setLayout(new GridLayout(3, 1));
    // markerName = new JLabel("", JLabel.CENTER);
    markerName = new JTextField("");
    markerName.setBorder(null);
    markerName.setEditable(false);
    markerName.setBackground(BACKGROUND_COLOR);
    markerName.setHorizontalAlignment(JTextField.CENTER);
    markerName.setFont(new Font("Arial", 0, 20));
    descrPanel.add(markerName);

    // commentLabel = new JLabel("", JLabel.CENTER);
    commentLabel = new JTextField("");
    commentLabel.setBorder(null);
    commentLabel.setEditable(false);
    commentLabel.setBackground(BACKGROUND_COLOR);
    commentLabel.setHorizontalAlignment(JTextField.CENTER);
    commentLabel.setFont(new Font("Arial", 0, 14));
    descrPanel.add(commentLabel);

    JPanel navigationPanel = new JPanel();
    first = new JButton(Grafik.getImageIcon("images/firstLast/First.gif"));
    first.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif"));
    first.addActionListener(this);
    first.setActionCommand(FIRST);
    first.setPreferredSize(new Dimension(20, 20));
    // first.setEnabled(false);
    previous = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif"));
    previous.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif"));
    previous.addActionListener(this);
    previous.setActionCommand(PREVIOUS);
    previous.setPreferredSize(new Dimension(20, 20));
    // previous.setEnabled(false);
    navigationField = new JTextField("", 8);
    navigationField.setHorizontalAlignment(JTextField.CENTER);
    navigationField.setFont(new Font("Arial", 0, 14));
    // navigationField.setEditable(false);
    navigationField.setBackground(BACKGROUND_COLOR);
    navigationField.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        try {
          int trav = Integer.valueOf(((JTextField) e.getSource()).getText().split("[\\s]+")[0])
                            .intValue()
                     - 1;
          if (trav >= 0 && trav < markerList.length) {
            markerIndex = trav;
          }
        } catch (NumberFormatException nfe) {
        }
        finishProcessing();
      }
    });

    next = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif"));
    next.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif"));
    next.addActionListener(this);
    next.setActionCommand(NEXT);
    next.setPreferredSize(new Dimension(20, 20));
    last = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif"));
    last.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif"));
    last.addActionListener(this);
    last.setActionCommand(LAST);
    last.setPreferredSize(new Dimension(20, 20));
    navigationPanel.add(first);
    navigationPanel.add(previous);
    navigationPanel.add(navigationField);
    navigationPanel.add(next);
    navigationPanel.add(last);

    navigationPanel.setBackground(BACKGROUND_COLOR);
    descrPanel.add(navigationPanel);
    descrPanel.setBackground(BACKGROUND_COLOR);
    return descrPanel;
  }

  private JPanel pcSliderPanel() {
    JPanel pcSliderPanel = new JPanel();
    // pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
    pcSliderPanel.setLayout(new GridLayout(2, 1));
    pcSliderPanel.setBackground(BACKGROUND_COLOR);

    JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
    // slider.setSize(new Dimension(150, 20));
    slider.setBackground(BACKGROUND_COLOR);
    // slider = new JSlider(JSlider.HORIZONTAL, 0,
    // Integer.parseInt(proj.getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), 0);
    slider = new JSlider(JSlider.HORIZONTAL, 0, proj.getProperty(proj.INTENSITY_PC_NUM_COMPONENTS),
                         0);
    slider.setValue(0);
    slider.setBackground(BACKGROUND_COLOR);
    if (pcResids == null) {
      pcLabel = new JLabel("PC file not detected", JLabel.CENTER);
    } else {
      pcLabel = new JLabel("Total PCs Loaded: " + pcResids.getNumComponents(), JLabel.CENTER);
    }
    pcLabel.setFont(new Font("Arial", Font.PLAIN, 16));
    pcSliderPanel.add(pcLabel);

    slider.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent ce) {
        JSlider slider = (JSlider) ce.getSource();
        numComponents = slider.getValue();
        pcLabel.setText("NumPCs = " + numComponents);
        // seletedScatterPanel.setPointsGeneratable(true);
        // seletedScatterPanel.setQcPanelUpdatable(true);
        // seletedScatterPanel.paintAgain();
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
          scatterPanel.setUpdateQCPanel(true);
          scatterPanel.paintAgain();
        }
      }
    });
    pcSliderPanel.add(slider);
    slider.setEnabled(pcResids != null);

    return pcSliderPanel;
  }

  private JPanel plotTypePanel() {
    JPanel plotTypePanel = new JPanel();
    plotTypePanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED),
                                                             "Plot Type", TitledBorder.LEADING,
                                                             TitledBorder.TOP));

    plotTypePanel.setLayout(new BoxLayout(plotTypePanel, BoxLayout.Y_AXIS));
    // plotTypePanel.setLayout(new GridLayout(3, 1));
    plotTypePanel.setBackground(BACKGROUND_COLOR);

    ItemListener typeListener = new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent ie) {
        JRadioButton jrb = (JRadioButton) ie.getItem();
        if (jrb.isSelected()) {
          // skip X/Y RAW
          for (int i = 0; i < ScatterPlot.TYPES.length; i++) {
            if (jrb.getText().equals(ScatterPlot.TYPES[i][0] + "/" + ScatterPlot.TYPES[i][1])
                || jrb.getText().equals("<html>" + ScatterPlot.TYPES[i][0] + "/<br />"
                                        + ScatterPlot.TYPES[i][1] + "</html>")) {
              plot_types[selectedPanelIndex] = i;
              // seletedScatterPanel.setPointsGeneratable(true);
              // seletedScatterPanel.setQcPanelUpdatable(true);
              for (ScatterPanel scatterPanel : scatterPanels) {
                scatterPanel.setPointsGeneratable(true);
                scatterPanel.setUpdateQCPanel(true);
              }
              updateGUI();
            }
          }
        }
      }
    };
    // --- Beginning of the original block ---
    ButtonGroup typeRadio = new ButtonGroup();
    typeRadioButtons = new JRadioButton[ScatterPlot.TYPES.length];
    // skip X/Y RAW
    for (int i = 0; i < ScatterPlot.TYPES.length; i++) {
      String label = ScatterPlot.TYPES[i][0] + "/" + ScatterPlot.TYPES[i][1];
      if (label.length() > 15) {
        label =
              "<html>" + ScatterPlot.TYPES[i][0] + "/<br />" + ScatterPlot.TYPES[i][1] + "</html>";
      }
      typeRadioButtons[i] = new JRadioButton(label, false);
      typeRadioButtons[i].setFont(new Font("Arial", 0, 14));
      typeRadio.add(typeRadioButtons[i]);
      typeRadioButtons[i].addItemListener(typeListener);
      typeRadioButtons[i].setBackground(BACKGROUND_COLOR);
      // tabPanel.add(typeRadioButtons[i], gbc);
      plotTypePanel.add(typeRadioButtons[i]);
    }
    // typeRadioButtons[0].setSelected(true);

    return plotTypePanel;
  }

  // private byte qc_prevChr;
  // private int[] qc_prevGeno;
  // private String[] qc_prevSex;
  // private String[] qc_prevOthCls;
  // private int qc_prevIndex;

  private void removeAnnotationFromMaps(char c) {
    // seletedScatterPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).remove(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0),
    // 0));
    // seletedScatterPanel.getActionMap().remove("annotation\t"+c);
    for (ScatterPanel scatterPanel : scatterPanels) {
      scatterPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
                  .remove(KeyStroke.getKeyStroke((c + "").toUpperCase().charAt(0), 0));
      scatterPanel.getActionMap().remove("annotation\t" + c);
    }
  }

  private void resetAfterLoad() {
    markerList = masterMarkerList;
    commentList = masterCommentList;
    showAllMarkersOrNot = true;
    showAnnotatedOrUnannotated = true;
    showAnnotationShortcuts = true;
    isInitilizing = true;
    indexOfAnnotationUsedAsMarkerList = -1;

    proj.getProgressMonitor().updateTask("SCATTERPLOT");

    loadMarkerDataFromList(0);

    proj.getProgressMonitor().updateTask("SCATTERPLOT");

    pcResids = loadPcResids();// returns null if not found, marker data should return original x/y
                              // if null
    numComponents = 0;// initialize to 0 PCs
    pedigree = proj.loadPedigree();// returns null if not found

    proj.getProgressMonitor().updateTask("SCATTERPLOT");
    loadCentroids();
    prepGCAdjustorLoad();
    proj.getProgressMonitor().updateTask("SCATTERPLOT");

    sessionID = (new Date().getTime() + "").substring(5);
    isClusterFilterUpdated = false;
    isAnnotationUpdated = false;
    autoSave = null;
    fail = !loadClusterFilterFiles();
    if (fail) {
      proj.getLog()
          .reportError("Chose to ignore prompt for autosaved cluster filters; ScatterPlot will not start");
      return;
    }

    annotationAutoAdv = true;

    if (hasAnnotationFile) {
      try {
        gcAnnotations = MarkerGCAnnotation.initForMarkers(proj, masterMarkerList,
                                                          annotationLoader.getMarkerSet(),
                                                          annotationLoader.getIndices());
        blastResults = MarkerBlastAnnotation.initForMarkers(masterMarkerList);
        ArrayList<AnnotationParser[]> parsers = new ArrayList<AnnotationParser[]>();
        parsers.add(blastResults);
        parsers.add(gcAnnotations);
        annotationLoader.fillAnnotations(masterMarkerList, parsers, QUERY_ORDER.NO_ORDER);
        // for (int i = 0; i < blastResults.length; i++) {
        // // blastResults[i].getMarkerSeqAnnotation().getTopBotProbe();
        // // blastResults[i].getMarkerSeqAnnotation().getTopBotRef();
        // // blastResults[i].getAnnotationLists()[0].get(0).getTag();
        // //
        // System.out.println(masterMarkerList[i]+"\t"+blastResults[i].getMarkerSeqAnnotation().getStrand()+"\t"+blastResults[i].getMarkerSeqAnnotation().getSequence());
        // System.out.println("REF -> " + masterMarkerList[i] + "\t" +
        // blastResults[i].getMarkerSeqAnnotation().getRef().getDisplayString());
        // for (int j = 0; j < blastResults[i].getMarkerSeqAnnotation().getAlts().length; j++) {
        // System.out.println("ALT -> " + masterMarkerList[i] + "\t" +
        // blastResults[i].getMarkerSeqAnnotation().getAlts()[j].getDisplayString());
        // }
        // System.out.println("A -> " + masterMarkerList[i] + "\t" +
        // blastResults[i].getMarkerSeqAnnotation().getA().getDisplayString());
        // System.out.println("B -> " + masterMarkerList[i] + "\t" +
        // blastResults[i].getMarkerSeqAnnotation().getB().getDisplayString());
        // System.out.println("STRAND -> " + masterMarkerList[i] + "\t" +
        // blastResults[i].getMarkerSeqAnnotation().getStrand());
        // }
      } catch (Exception e) {
        blastResults = null;
        proj.getLog().reportException(e);
      }
    }
  }

  private void setFourView() {
    // some roughly correct combination of invalidation calls...
    viewPanel.remove(indivPanels[selectedPanelIndex]);
    scatterOverview.removeAll();
    for (int i = 0; i < indivPanels.length; i++) {
      scatterPanels[i].axisXHeight = 0;
      scatterPanels[i].axisYWidth = 0;
      scatterPanels[i].shrunk = true;
      scatterPanels[i].displayXaxis = false;
      scatterPanels[i].displayYaxis = false;
      scatterOverview.add(indivPanels[i]);
      indivPanels[i].invalidate();
      scatterPanels[i].invalidate();
    }
    scatterOverview.invalidate();
    viewPanel.add(scatterOverview, BorderLayout.CENTER);
    getContentPane().invalidate();
    showingAll = true;
  }

  private void setHistogram(int length) {
    int[] histotemp = blastResults[markerIndex].getAlignmentHistogram(getProject());
    final int[] histogram = new int[histotemp.length];
    for (int i = 0; i < histogram.length; i++) {
      if (i < length) {
        histogram[i] = 0;
      } else {
        histogram[i] = histotemp[i];
      }
    }
    histFrame.removeAllData();
    histFrame.setHistogram(true);
    histFrame.getPanel().overrideAxisLabels("Bins", "");
    histFrame.getPanel().setHistogramOverride(true);
    histFrame.getPanel().setForceXAxisWholeNumbers(true);
    histFrame.getPanel().setForceYAxisWholeNumbers(true);


    final double[] bins = new double[histogram.length];
    int min = 0;
    for (int i = 0; i < histogram.length; i++) {
      bins[i] = i + 1;
      if (histogram[i] == 0) {
        min++;
      }
    }
    final int min2 = min;
    Histogram hist = new Histogram() {
      private static final long serialVersionUID = 1L;
      {
        min = min2;
        max = 50;
        sigfigs = 4;
        counts = histogram;
      }

      @Override
      public double determineStep() {
        return 1d;
      }

      @Override
      public double[] getBins() {
        return bins;
      }
    };
    histFrame.getPanel().setHistogram(hist);
  }

  private void setIndividualView() {
    // some roughly correct combination of invalidation calls...
    viewPanel.remove(scatterOverview);
    viewPanel.add(indivPanels[selectedPanelIndex], BorderLayout.CENTER);
    scatterOverview.removeAll();
    scatterOverview.invalidate();
    for (int i = 0; i < indivPanels.length; i++) {
      scatterPanels[i].axisXHeight = AbstractPanel.HEIGHT_X_AXIS;
      scatterPanels[i].axisYWidth = AbstractPanel.WIDTH_Y_AXIS;
      scatterPanels[i].shrunk = false;
      scatterPanels[i].displayXaxis = true;
      scatterPanels[i].displayYaxis = true;
      indivPanels[i].invalidate();
      scatterPanels[i].invalidate();
    }
    getContentPane().invalidate();
    showingAll = false;
  }

  private JPanel sizeSliderPanel() {
    JPanel sizeSliderPanel = new JPanel();
    // sizeSliderPanel.setLayout(new BoxLayout(sizeSliderPanel, BoxLayout.Y_AXIS));
    sizeSliderPanel.setLayout(new GridLayout(2, 1));
    sizeSliderPanel.setBackground(BACKGROUND_COLOR);

    JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
    // slider.setSize(new Dimension(150, 20));
    slider.setBackground(BACKGROUND_COLOR);
    sizeLabel = new JLabel("Size: " + DEFAULT_SIZE, JLabel.CENTER);
    sizeLabel.setFont(new Font("Arial", Font.PLAIN, 16));
    // tabPanel.add(sizeLabel, gbc);
    sizeSliderPanel.add(sizeLabel);

    slider.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent ce) {
        JSlider slider = (JSlider) ce.getSource();
        sizeLabel.setText("Size = " + slider.getValue());
        size = (byte) slider.getValue();
        // seletedScatterPanel.setPointsGeneratable(true);
        // seletedScatterPanel.setQcPanelUpdatable(false);
        // seletedScatterPanel.paintAgain();
        for (ScatterPanel scatterPanel : scatterPanels) {
          scatterPanel.setPointsGeneratable(true);
          scatterPanel.setUpdateQCPanel(false);
          scatterPanel.paintAgain();
        }
      }
    });
    // tabPanel.add(slider, gbc);
    sizeSliderPanel.add(slider);

    return sizeSliderPanel;
  }

  private void swapView() {
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        if (showingAll) {
          setIndividualView();
        } else {
          setFourView();
        }
        ScatterPlot.this.validate();
        ScatterPlot.this.repaint();
        updateGUI();
      }
    });
  }

  private void updateBLASTPanel() {
    blastPanel.removeAll();
    blastPanel.repaint();

    Font lblFont = new Font("Arial", 0, 14);
    if (!hasAnnotationFile || blastResults == null || blastResults.length == 0
        || markerIndex >= blastResults.length || !blastResults[markerIndex].isFound()) {
      JLabel blastLabel = new JLabel("BLAST Metrics Unavailable", SwingConstants.CENTER);
      blastLabel.setFont(lblFont);
      blastPanel.add(blastLabel, "dock center, grow");
      blastPanel.revalidate();
      blastPanel.repaint();
      return;
    }

    MarkerBlastAnnotation blastResult = blastResults[markerIndex];

    JLabel typeLabel = new JLabel("Has Perfect Match? ", JLabel.LEFT);
    typeLabel.setFont(lblFont);
    blastPanel.add(typeLabel, "cell 0 0");
    boolean has = blastResult.hasPerfectMatch(log);
    String result = has ? "YES" : "NO";
    typeLabel = new JLabel(result, JLabel.LEFT);
    typeLabel.setFont(lblFont);
    typeLabel.setForeground(has ? Color.GREEN : Color.RED);
    blastPanel.add(typeLabel, "cell 1 0, alignx right");

    double filter = proj.BLAST_PROPORTION_MATCH_FILTER.getValue();
    int probe = proj.ARRAY_TYPE.getValue().getProbeLength();
    int alignFilter = blastFrame == null ? (int) (filter * probe) : blastFrame.currentAlignFilter;
    int[] hist = blastResult.getAlignmentHistogram(proj);
    int offTLbls = 0;
    for (int i = hist.length - 1; i >= 0; i--) {
      int len = probe - hist.length + i; // compensate for a histogram that could be of shorter
                                         // length than the probe
      if (len >= alignFilter) {
        if (i == hist.length - 1) {
          if (has && hist[i] > 0) {
            offTLbls--; // remove perfect match
          } else if (has && hist[i] == 0) {
            // error, no record for perfect match in histogram
          }
        }
        offTLbls += hist[i];
      } else {
        break;
      }
    }

    double ratio = blastFrame == null ? filter
                                      : ((double) blastFrame.currentAlignFilter)
                                        / ((double) proj.ARRAY_TYPE.getValue().getProbeLength());
    // int offTLbls = BlastFrame.BlastUtils.filterAnnotations(proj,
    // blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS, log),
    // alignFilter).size();
    typeLabel = new JLabel("# Off-Target Alignments (>" + ratio + "% match): ", JLabel.LEFT);
    typeLabel.setFont(lblFont);
    typeLabel.setForeground(offTLbls > 0 ? Color.RED : Color.BLACK);
    blastPanel.add(typeLabel, "cell 0 1");

    typeLabel = new JLabel("" + offTLbls, JLabel.LEFT);
    typeLabel.setFont(lblFont);
    blastPanel.add(typeLabel, "cell 1 1, alignx right");

    typeLabel = new JLabel("# On-Target (mismatched) Alignments: ", JLabel.LEFT);
    typeLabel.setFont(lblFont);
    blastPanel.add(typeLabel, "cell 0 2");
    ArrayList<BlastAnnotation> onTaligns =
                                         blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT,
                                                                       log);
    String lbl = onTaligns == null ? " 0" : " " + onTaligns.size();
    typeLabel = new JLabel(lbl);
    typeLabel.setFont(lblFont);
    blastPanel.add(typeLabel, "cell 1 2, alignx right");

    JButton blastButton = new JButton();
    blastButton.setActionCommand(BLAST_DETAILS_COMMAND);
    blastButton.setText("Show BLAST Results");
    blastButton.setMnemonic(KeyEvent.VK_B);
    blastButton.addActionListener(this);

    // if (referenceGenome == null) {
    // blastButton.setToolTipText("Error - No reference genome found!");
    // blastButton.setEnabled(false);
    // }

    blastPanel.add(blastButton, "cell 0 3 2 1, center");

    blastPanel.revalidate();
    blastPanel.repaint();
  }

  public static void createAndShowGUI(final Project proj, final String[] markerList,
                                      final String[] commentList, final boolean exitOnClose) {
    new Thread(new Runnable() {
      @Override
      public void run() {
        final ScatterPlot scatterPlot;
        scatterPlot = new ScatterPlot(proj, markerList, commentList, exitOnClose);
        if (!scatterPlot.failed()) {
          SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
              scatterPlot.pack();
              scatterPlot.setSize(1200, 870);
              scatterPlot.setVisible(true);
            }
          });
        } else {
          proj.getLog().reportError("Error - ScatterPlot failed to initialize.");
        }
      }
    }).start();
  }

  public static void main(String[] args) {
    javax.swing.SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        createAndShowGUI(new Project(org.genvisis.cnv.Launch.getDefaultDebugProjectFile(true),
                                     false),
                         null, null, true);
      }
    });
  }
}
