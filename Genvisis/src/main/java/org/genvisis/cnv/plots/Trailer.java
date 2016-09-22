package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.ActionMap;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListCellRenderer;
import javax.swing.InputMap;
import javax.swing.JButton;
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
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.border.EmptyBorder;

import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.analysis.MosaicismDetect;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicBuilder;
import org.genvisis.cnv.analysis.MosaicismQuant;
import org.genvisis.cnv.analysis.MosaicismQuant.MOSAIC_TYPE;
import org.genvisis.cnv.analysis.MosaicismQuant.MosaicQuantResults;
import org.genvisis.cnv.analysis.MosaicismQuant.MosaicQuantWorker;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.gui.ClickListener;
import org.genvisis.cnv.gui.FileActionMenu;
import org.genvisis.cnv.gui.NewRegionListDialog;
import org.genvisis.cnv.gui.RegionNavigator;
import org.genvisis.cnv.gui.SingleClick;
import org.genvisis.cnv.gui.RegionNavigator.ChrNavigator;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.hmm.CNVCaller.CNVCallResult;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.hmm.PFB;
import org.genvisis.cnv.hmm.PennHmm;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.cnv.plots.ColorExt.ColorManager;
import org.genvisis.cnv.prop.StringProperty;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.var.IndiPheno;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.cnv.var.Region;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.TransferableImage;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariant.CNVBuilder;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.mining.Transformations;

import net.miginfocom.swing.MigLayout;

public class Trailer extends JFrame implements ChrNavigator, ActionListener, ClickListener, MouseListener,
                     MouseMotionListener, MouseWheelListener {

  public static final long serialVersionUID = 1L;

  public static final String DEFAULT_LOCATION = "chr1";

  public static final String DEFAULT_SAMPLE = null;
  public static final boolean SHOW_MIDLINE = true;

  public static final double MOUSE_WHEEL_MULTIPLIER = 0.5;
  public static final int WIDTH_BUFFER = 25;
  public static final int HEIGHT_BUFFER = 10;
  // public static final int DYNAMIC_HEIGHT_LIMIT = 2000;
  public static final int DYNAMIC_HEIGHT_LIMIT = 0;
  public static final int DOUBLE_CLICK_INTERVAL = 500;
  public static final int SIZE = 4;
  public static final int MIN_BUFFER = 10000;
  public static final int DEFAULT_STARTX = 20;
  public static final int DEFAULT_STARTY = 20;
  public static final int DEFAULT_WIDTH = 1100;
  public static final int DEFAULT_HEIGHT = 720;
  public static final String REGION_LIST_USE_CNVS = "Use CNVs as Regions...";

  private static final String BLANK_COMMENT = " -- Click here to comment on this region -- ";

  private static final String FIRST_REGION = "First region";
  private static final String PREVIOUS_REGION = "Previous region";
  private static final String NEXT_REGION = "Next region";
  private static final String LAST_REGION = "Last region";
  private static final String TO_SCATTER_PLOT = "To Scatter Plot";
  private static final String TO_COMP_PLOT = "To Comp Plot";
  private static final String REGION_LIST_LOAD_FILE = "Load Region File";
  private static final String REGION_LIST_SAVE_FILE = "Save Regions to File...";
  private static final String REGION_LIST_ADD = "Add current region";
  private static final String REGION_LIST_NEW_FILE = "New Region List File";
  private static final String REGION_LIST_PLACEHOLDER = "Select Region File...";
  private static final String REGION_LIST_PROVIDED = "Region list provided";
  private JComboBox sampleList;
  private RegionNavigator navPanel;
  private String[] samplesPresent;
  private JButton previousRegion, nextRegion, firstRegion, lastRegion;
  private Project proj;
  private String sample;
  private boolean jar;
  private IndiPheno indiPheno;
  private String[] markerNames;
  private PreparedMarkerSet markerSet;
  private long fingerprint;
  private int[] positions;
  private boolean[] dropped;
  private int[][] chrBoundaries;
  private float[] lrrs, lrrValues;
  private float[] bafs, originalBAFs;
  private byte[] genotypes;
  private byte chr;
  private int start, startMarker;
  private int stop, stopMarker;
  // private float lrrMin, lrrMax;
  private boolean inDrag;
  private int startX;
  private String[] cnvLabels;
  private CNVariant[][] cnvs;
  // Variables for whether or not to prompt the user to save regions after given modifications are
  // made.
  private boolean promptCommentSave = true;
  private boolean promptAddRegionSave = true;
  // private String[] regionsList;
  // private int regionsListIndex;

  // FIXME should use org.genvisis.cnv.var.Region class here
  private String[][] regions;
  private int regionIndex;
  private JPanel lrrPanel;
  private JPanel bafPanel;
  private JPanel cnvPanel;
  private SingleClick leftClick;
  private SingleClick rightClick;
  private GeneTrack track;
  private SampleData sampleData;
  private JLabel commentLabel;
  private JTextField commentField;
  private JLabel qcLabel;
  private int transformation_type;
  private boolean transformSeparatelyByChromosome;
  private GcModel gcModel;
  private PennHmm pennHmm;
  private PFB pfb;
  private JCheckBoxMenuItem gcCorrectButton;
  private JCheckBoxMenuItem callCnvsButton;
  private JCheckBoxMenuItem mosaicCallButton;
  private JCheckBoxMenuItem mosaicFButton;
  private JCheckBoxMenuItem beastHButton;

  private String[] qcGenome;
  // private String[] qcChromo;
  private String[] qcRegion;
  private int qcSelection = 0;

  private Hashtable<String, String> namePathMap;
  // private JComboBox<String> centroidsSelection;
  private Logger log;
  private boolean fail;
  boolean isSettingCentroid = false;
  private float[][][] centroids;
  private String currentCentroid;
  private static final String SEX_CENT = "Sex-Specific";
  private JMenu loadRecentFileMenu;
  // private JMenuItem launchScatter;
  // private JMenuItem launchComp;
  private ButtonGroup regionButtonGroup;
  private GCParameterDisplay gcParameterDisplay;
  private Sample samp;
  private ColorManager<String> currentColorManager;
  private ColorManager<String> excludeManager;
  private String[] otherColors;
  private Hashtable<String, ColorManager<String>> previouslyLoadedManagers;
  private Thread updateQCThread = null;
  private StringProperty lastRegionFile;

  private final AbstractAction markerFileSelectAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      String shortName = ((JCheckBoxMenuItem) e.getSource()).getText();
      if (!loadingFile && !REGION_LIST_LOAD_FILE.equals(shortName)
          && !REGION_LIST_PLACEHOLDER.equals(shortName)
          && !REGION_LIST_USE_CNVS.equals(shortName)) {
        String file = regionFileNameLoc.get(shortName);
        if (file == null || file.equals(regionFileName)) {
          return;
        }
        String tempFile = file.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + file : file;
        if (!Files.exists(tempFile)) {
          proj.message("Error - region file '" + shortName + "' doesn't exist.");
          regionFileNameBtn.get(shortName).setSelected(true);
        } else {
          regionFileName = file;
          loadRegions();
          showRegion(0);
        }
      } /*
         * else if (loadingFile && REGION_LIST_PLACEHOLDER.equals(shortName)) { // do nothing }
         */else if (loadingFile || REGION_LIST_PLACEHOLDER.equals(shortName)) {
        // leave as currently selected marker
        if (regionFileName != "" && regionFileName != null) {
          String file = regionFileName;
          if (!REGION_LIST_USE_CNVS.equals(regionFileName)) {
            file = ext.rootOf(regionFileName);
          }
          regionFileNameBtn.get(file).setSelected(true);
        }
        return;
      } else if (REGION_LIST_USE_CNVS.equals(shortName)) {
        if (!shortName.equals(regionFileName)) {
          regionFileName = REGION_LIST_USE_CNVS;
          loadCNVsAsRegions();
          updateCNVs(chr);
          updateGUI();
          showRegion(0);
        }
      }
    }
  };

  /**
   * Menu action for adding current visible area as a region to the current list.
   */
  AbstractAction addRegionAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      if (regions != null) {
        String[][] tmp = new String[regions.length + 1][];
        System.arraycopy(regions, 0, tmp, 0, regions.length);
        tmp[tmp.length - 1] = new String[] {sample, "chr" + chr + ":" + start + "-" + stop};
        regions = tmp;
        promptAddRegionSave = promptAndSaveRegions(promptAddRegionSave);
        showRegion(regions.length - 1);
      }
    }

  };

  /**
   * Menu action for saving the current region list to a file.
   */
  AbstractAction saveRegionFileAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      saveRegionFile();
    }
  };

  AbstractAction loadRegionFileAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      String newFile = chooseNewFiles();
      if (newFile == null) {
        if (regionFileName != null && !"".equals(regionFileName)) {
          if (REGION_LIST_USE_CNVS.equals(regionFileName)) {
            regionFileNameBtn.get(REGION_LIST_USE_CNVS).setSelected(true);
          } else {
            regionFileNameBtn.get(ext.rootOf(regionFileName)).setSelected(true);
          }
        }
      } else {
        String file = ext.verifyDirFormat(newFile);
        file = file.substring(0, file.length() - 1);
        String name = ext.rootOf(file);
        regionFileNameBtn.get(name).setSelected(true);
        regionFileNameBtn.get(name).doClick();
      }
    }
  };

  AbstractAction screencapAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      JFileChooser jfc = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
      jfc.setMultiSelectionEnabled(false);
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save Screen Capture...");
      jfc.setDialogType(JFileChooser.SAVE_DIALOG);
      int code = jfc.showSaveDialog(Trailer.this);
      if (code == JFileChooser.APPROVE_OPTION) {
        String filename = jfc.getSelectedFile().getAbsolutePath();
        doScreenCapture(filename);
      }
    }
  };

  AbstractAction screencapClipboardAction = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      doScreenCapture(null);
    }
  };

  // private Color[] colorScheme = {new Color(33, 31, 53), // dark dark
  // new Color(23, 58, 172), // dark blue
  // new Color(201, 30, 10), // deep red
  // new Color(140, 20, 180), // deep purple
  // new Color(33, 87, 0), // dark green
  // new Color(55, 129, 252), // light blue
  // new Color(217, 109, 194), // pink
  // new Color(94, 88, 214), // light purple
  // new Color(189, 243, 61), // light green
  // };

  static ArrayList<Color[]> getColor() {
    ArrayList<Color[]> colors = new ArrayList<Color[]>();
    colors.add(new Color[] {new Color(23, 58, 172), new Color(55, 129, 252)});
    colors.add(new Color[] {new Color(140, 20, 180), new Color(94, 88, 214)}); // deep/light purple
    colors.add(new Color[] {new Color(33, 87, 0), new Color(189, 243, 61)}); // dark green
    colors.add(new Color[] {new Color(201, 30, 10), new Color(217, 109, 194)}); // deep red/pink
    colors.add(new Color[] {new Color(33, 31, 53), new Color(255, 255, 255)});
    return colors;
  }

  private final ArrayList<Color[]> colorScheme = getColor();
  // private Color[][] colorScheme = {
  // {new Color(23, 58, 172), new Color(55, 129, 252)}, // dark/light blue
  // {new Color(140, 20, 180), new Color(94, 88, 214)}, // deep/light purple
  // {new Color(33, 87, 0), new Color(189, 243, 61)}, // dark green
  // {new Color(201, 30, 10), new Color(217, 109, 194)}, // deep red/pink
  // {new Color(33, 31, 53), new Color(255, 255, 255)} // dark dark/ light light
  // };

  final ArrayList<int[]> activeCNVs = new ArrayList<int[]>();
  volatile int[] selectedCNV = null;
  private final HashMap<String, JCheckBoxMenuItem> regionFileNameBtn =
                                                                     new HashMap<String, JCheckBoxMenuItem>();
  private final HashMap<String, String> regionFileNameLoc = new HashMap<String, String>();
  // private JComboBox<String> regionFileList;
  private String regionFileName;
  private volatile boolean loadingFile = false;
  private JTextField regionField;
  private HashMap<String, JCheckBoxMenuItem> centButtonMap;
  private JCheckBoxMenuItem autoSwitch;

  class CustomCallPopUp extends JPopupMenu {
    /**
     *
     */
    private static final long serialVersionUID = 1L;
    JMenuItem anItem;
    JButton button;

    public CustomCallPopUp() {

      class QuantButton extends JMenuItem {
        /**
         *
         */
        private static final long serialVersionUID = 1L;
        // private Segment toQuant;

        public QuantButton(final Segment toQuant, final CustomCallPopUp callPopUp) {
          super("Quantify Mosaicism for " + toQuant.getUCSClocation());
          // this.toQuant = toQuant;
          addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
              quantHere(toQuant, true);

            }
          });
        }

      }
      class BeastButton extends JMenuItem {
        /**
         *
         */
        private static final long serialVersionUID = 1L;

        // private Segment toQuant;

        public BeastButton(final Segment toQuant, final CustomCallPopUp callPopUp) {
          super("Assign BEAST height for " + toQuant.getUCSClocation());
          // this.toQuant = toQuant;
          addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
              beastHere(toQuant, INTERNAL_CNV_TYPES.BEAST_SCORE);

            }
          });
        }

      }
      class RemoveButton extends JMenuItem {
        /**
         *
         */
        private static final long serialVersionUID = 1L;

        // private Segment toQuant;

        public RemoveButton(final Segment toQuant, final int externalIndex,
                            final INTERNAL_CNV_TYPES type, final CustomCallPopUp callPopUp,
                            final String title, final boolean all) {
          super(title);
          // this.toQuant = toQuant;
          addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
              removeCnvFromPheno(toQuant, externalIndex, all, type);
            }
          });
        }

      }

      if (selectedCNV != null && selectedCNV.length > 1 && cnvs != null && cnvs.length > 0) {
        int externalCnvs = prepInternalClasses();

        Segment toQuant = cnvs[selectedCNV[0]][selectedCNV[1]];
        QuantButton qb = new QuantButton(toQuant, this);
        qb.setFont(new Font("Arial", 0, 12));
        add(qb);
        addSeparator();
        BeastButton bb = new BeastButton(toQuant, this);
        bb.setFont(new Font("Arial", 0, 12));
        add(bb);
        addSeparator();
        INTERNAL_CNV_TYPES type = INTERNAL_CNV_TYPES.getAppropriate(selectedCNV[0], externalCnvs);
        if (type != null) {
          RemoveButton rb = new RemoveButton(toQuant, externalCnvs, type, this,
                                             "Remove " + toQuant.getUCSClocation(), false);
          rb.setFont(new Font("Arial", 0, 12));
          add(rb);
          addSeparator();
          RemoveButton rbAll = new RemoveButton(toQuant, externalCnvs, type, this,
                                                "Remove all from " + type.getName(), true);
          rbAll.setFont(new Font("Arial", 0, 12));
          add(rbAll);
        }
        pack();
      }
    }
  }

  MouseAdapter cnvAdapter = new MouseAdapter() {
    int defaultInitial = ToolTipManager.sharedInstance().getInitialDelay();
    int defaultReshow = ToolTipManager.sharedInstance().getReshowDelay();

    // Object defaultBG = UIManager.get("ToolTip.background");
    @Override
    public void mouseEntered(MouseEvent e) {
      super.mouseEntered(e);
      ToolTipManager.sharedInstance().setReshowDelay(3);
      ToolTipManager.sharedInstance().setInitialDelay(3);
    }

    @Override
    public void mouseExited(MouseEvent e) {
      super.mouseExited(e);
      ToolTipManager.sharedInstance().setReshowDelay(defaultReshow);
      ToolTipManager.sharedInstance().setInitialDelay(defaultInitial);
    }

    @Override
    public void mouseMoved(MouseEvent e) {
      super.mouseMoved(e);
      if (selectedCNV == null) {
        return;
      }
      int x = e.getX();
      int y = e.getY();

      if (cnvs.length <= selectedCNV[0] || cnvs[selectedCNV[0]].length <= selectedCNV[1]) {
        // UIManager.put("ToolTip.background", defaultBG);
        cnvPanel.setToolTipText(null);
        return;
      }

      CNVariant cnv = cnvs[selectedCNV[0]][selectedCNV[1]];
      int cnvX1 = getX(cnv.getStart());
      int cnvX2 = getX(cnv.getStop());
      int cnvY1 = (selectedCNV[0] + 2) * 15;
      int cnvY2 = cnvY1 + 10;
      if (x >= cnvX1 && x <= cnvX2 && y > cnvY1 && y < cnvY2) {
        if (cnvPanel.getToolTipText() != null) {
          // Point locationOnScreen = new Point(e.getXOnScreen(), e.getYOnScreen());
          // Point locationOnComponent = new Point(e.getX(), e.getY());
          // final MouseEvent me = new MouseEvent(cnvPanel, -1, System.currentTimeMillis(), 0,
          // locationOnComponent.x, locationOnComponent.y, locationOnScreen.x, locationOnScreen.y,
          // 0, false, 0);
          // SwingUtilities.invokeLater(new Runnable() {
          // @Override
          // public void run() {
          // ToolTipManager.sharedInstance().mouseMoved(me);
          // }
          // });
        } else {
          StringBuilder txtBld = new StringBuilder();
          txtBld.append("<html>Start: ").append(ext.addCommas(cnv.getStart())).append("<br/>");
          txtBld.append(" Stop: ").append(ext.addCommas(cnv.getStop())).append("<br/>");
          txtBld.append(" Length: ").append(ext.addCommas((cnv.getStop() - cnv.getStart())))
                .append("<br/>");
          txtBld.append("# Mkrs: ").append(cnv.getNumMarkers()).append("<br/>");
          txtBld.append("CN: ").append(cnv.getCN()).append("<br/>");
          txtBld.append("Score: ").append(cnv.getScore()).append("</html>");
          // UIManager.put("ToolTip.background", new
          // javax.swing.plaf.ColorUIResource(colorScheme[selectedCNV[0]][cnv.getCN() < 2 ? 0 :
          // 1]));
          // cnvPanel.createToolTip();
          cnvPanel.setToolTipText(txtBld.toString());
        }
      } else {
        // UIManager.put("ToolTip.background", defaultBG);
        cnvPanel.setToolTipText(null);
      }
    }

    @Override
    public void mousePressed(MouseEvent e) {
      if (e.isPopupTrigger()) {
        if (selectedCNV != null) {
          doPop(e);
        }
      } else {
        startX = e.getPoint().x;
        inDrag = true;
      }
    }

    @Override
    public void mouseReleased(MouseEvent e) {
      if (e.isPopupTrigger()) {
        if (selectedCNV != null) {
          doPop(e);
        }
      }
      inDrag = false;
    }

    @Override
    public void mouseDragged(MouseEvent e) {
      int curX = e.getPoint().x;
      int distance = startX - curX;

      distance *= (stop - start) / (getWidth() - 2 * WIDTH_BUFFER);

      if (distance < 0) {
        distance = Math.max(distance, 1 - start);
      } else {
        distance = Math.min(distance, positions[chrBoundaries[chr][1]] - stop);
      }

      if ((start <= 1 && distance < 0)
          || (stop >= positions[chrBoundaries[chr][1]] && distance > 0)) {

      } else {
        start += distance;
        stop += distance;
      }

      if (inDrag) {
        updateGUI();
        startX = curX;
      }
    }


    private void doPop(MouseEvent e) {
      CustomCallPopUp menu = new CustomCallPopUp();
      menu.show(e.getComponent(), e.getX(), e.getY());

    }

    @Override
    public void mouseClicked(MouseEvent e) {
      super.mouseClicked(e);
      Trailer.this.requestFocusInWindow();
      int x = e.getX();
      int y = e.getY();

      for (int[] cnvInd : activeCNVs) {
        int yMin = (cnvInd[0] + 2) * 15;
        int yMax = yMin + 10;
        int xBegin = getX(cnvs[cnvInd[0]][cnvInd[1]].getStart());
        int xEnd = getX(cnvs[cnvInd[0]][cnvInd[1]].getStop());

        if (x >= xBegin && x <= xEnd && y > yMin && y < yMax) {
          selectedCNV = cnvInd;

          MouseEvent phantom = new MouseEvent(e.getComponent(), MouseEvent.MOUSE_MOVED,
                                              System.currentTimeMillis(), 0, x, y, 0, false);
          ToolTipManager.sharedInstance().mouseMoved(phantom); // order of mouseMoved calls doesn't
                                                               // matter, but both are necessary
          mouseMoved(phantom);
          Trailer.this.repaint();
          return;
        }

      }

      selectedCNV = null;
      Trailer.this.repaint();
    }
  };


  // TODO trailer takes a map<genome : region>
  // builds region array, sets region filename to "From enumerated list"

  // public Trailer(Project proj, String[] cnvFiles, Map<String, List<String>> regionMap) {
  // this(proj, regionMap.keySet().iterator().next(), cnvFiles);
  // cnvLabels = regionMap.keySet().toArray(new String[regionMap.size()]);
  // samplesPresent = cnvLabels;
  // regions = new String[cnvLabels.length][];
  // for (int i=0; i<regions.length; i++) {
  // List<String> list = regionMap.get(cnvLabels[i]);
  // regions[i] = list.toArray(new String[list.size()]);
  // }
  // regionFileName = REGION_LIST_PROVIDED;
  // setBounds(DEFAULT_STARTX, DEFAULT_STARTY, DEFAULT_WIDTH, DEFAULT_HEIGHT);
  // setVisible(true);
  // }

  /**
   * Basic constructor. Use {@link #Trailer(Project, String, String[], String, String[][])} to
   * specify a region list.
   *
   * @see #Trailer(Project, String, String[], String, String[][], int, int, int, int)
   */
  public Trailer(Project proj, String selectedSample, String[] cnvFiles, String startingLocation) {
    this(proj, selectedSample, cnvFiles, startingLocation, null);
  }

  /**
   * Basic constructor allowing specified GUI dimensions. Use
   * {@link #Trailer(Project, String, String[], String, String[][])} to specify a region list.
   *
   * @see #Trailer(Project, String, String[], String, String[][], int, int, int, int)
   */
  public Trailer(Project proj, String selectedSample, String[] cnvFiles,
                 final String startingLocation, final int startX, final int startY, final int width,
                 final int height) {
    this(proj, selectedSample, cnvFiles, startingLocation, null, startX, startY, width, height);
  }

  /**
   * As {@link #Trailer(Project, String, String[], String, String[][])} but uses the first of the
   * specified regions to provide the default selectetd sample/starting location.
   *
   * @see #Trailer(Project, String, String[], String, String[][], int, int, int, int)
   */
  public Trailer(Project proj, String[] cnvFiles, String[][] specifiedRegions) {
    this(proj, specifiedRegions[0][0], cnvFiles, specifiedRegions[0][1], specifiedRegions);
  }

  /**
   * Constructor allowing a custom region list to be defined.
   *
   * @see #Trailer(Project, String, String[], String, String[][], int, int, int, int)
   */
  public Trailer(Project proj, String selectedSample, String[] cnvFiles, String startingLocation,
                 String[][] specifiedRegions) {
    this(proj, selectedSample, cnvFiles, startingLocation, specifiedRegions, DEFAULT_STARTX,
         DEFAULT_STARTY, DEFAULT_WIDTH, DEFAULT_HEIGHT);
  }

  /**
   * Create and show a trailer plot. Trailers visualize the LRR and B allele frequency of a
   * chromsomoal region of one sample genome (DNA).
   *
   * @param proj Project containing sample data of interest
   * @param selectedSample The first sample to be selected
   * @param cnvFiles List of CNV files (e.g. to read regions from)
   * @param startingLocation First chromosomal location to load (of the format:
   *        {@code chr<##>[:startPos-endPos]}
   * @param specifiedRegions (OPTIONAL) a hard-coded enumeration of desired regions. The number of
   *        regions = specifiedRegions.length. Each subarray should have the format:
   *        {@code {sampleID, location, [comment]}}
   * @param startX GUI starting x position
   * @param startY GUI starting y positions
   * @param width GUI width
   * @param height GUI height
   */
  public Trailer(final Project proj, String selectedSample, String[] cnvFiles,
                 final String startingLocation, final String[][] specifiedRegions, final int startX,
                 final int startY, final int width, final int height) {
    // TODO Trailer should have a createAndShowGUI, same as all the other plots, as opposed to being
    // its own frame
    super("Genvisis - Trailer - " + proj.PROJECT_NAME.getValue());
    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent e) {
        if (Trailer.this.proj != null) {
          ArrayList<String> files = new ArrayList<String>(regionFileNameLoc.values());
          String[] currSet = Trailer.this.proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue();

          ArrayList<String> newSet = new ArrayList<String>();
          for (String s : files) {
            boolean found = false;
            for (int i = 0; i < currSet.length; i++) {
              if (currSet[i].equals(s)) {
                found = true;
                break;
              }
            }
            if (!found) {
              newSet.add(s);
            }
          }

          if (newSet.size() > 0) {
            String[] newList = files.toArray(new String[] {});

            String message = newSet.size() + " files have been added.  ";
            int choice =
                       JOptionPane.showOptionDialog(null,
                                                    message
                                                          + " Would you like to keep this configuration for the next time Trailer is loaded?",
                                                    "Preserve Trailer workspace?",
                                                    JOptionPane.YES_NO_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (choice == 0) {
              Trailer.this.proj.INDIVIDUAL_CNV_LIST_FILENAMES.setValue(newList);
            }
          }
          proj.TRAILER_REGION.setValue(regionFileName);
          Trailer.this.proj.saveProperties();
        }
        super.windowClosing(e);
      }
    });

    System.out.println(Array.toStr(cnvFiles, "; "));

    long time;

    this.proj = proj;
    log = proj.getLog();
    jar = proj.JAR_STATUS.getValue();
    fail = false;

    time = new Date().getTime();

    fail = !loadMarkers();
    if (fail) {
      return;
    }
    if (Files.exists(proj.GC_MODEL_FILENAME.getValue(false, false))) {
      gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false),
                                                    true, proj.getLog());
    } else {
      gcModel = null;
    }
    otherColors = Files.list(proj.DATA_DIRECTORY.getValue(), null, ".gcmodel", true, false, true);
    otherColors = Array.concatAll(otherColors, Files.list(proj.PROJECT_DIRECTORY.getValue(), null,
                                                          ".gcmodel", true, false, true));

    generateComponents();
    setJMenuBar(createMenuBar());

    sample = selectedSample == null ? samplesPresent[0] : selectedSample;
    sampleData = proj.getSampleData(2, cnvFiles);
    if (sampleData.failedToLoad()) {
      proj.getLog().reportError("Without a SampleData file, Trailer will not start");
      return;
    }

    time = new Date().getTime();

    Resource gTrack = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), log).getGTrack();
    if (gTrack.isAvailable()) {
      log.report("Loading track from " + gTrack.get());
      track = GeneTrack.load(gTrack.get(), jar);
      log.report("Loaded track in " + ext.getTimeElapsed(time));
    }

    System.out.println("startX: " + startX + "\t startY: " + startY + "\t width: " + width
                       + "\t height: " + height);

    parseLocation(startingLocation);
    addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        Trailer.this.requestFocusInWindow();
      }
    });

    new Thread(new Runnable() {
      @Override
      public void run() {
        final long t = new Date().getTime();


        if (specifiedRegions != null) {
          regions = specifiedRegions;
          regionFileName = REGION_LIST_PROVIDED;
        } else {
          String lastRegionFile = proj.TRAILER_REGION.getValue();
          if (new File(lastRegionFile).exists()) {
            regionFileName = lastRegionFile;
            loadRegionFile(regionFileName);
            loadRegions();
          }
        }

        if (regionFileName == null || REGION_LIST_USE_CNVS.equals(regionFileName)) {
          while (!sampleData.getCNVsLoaded()) {
            try {
              Thread.sleep(30);
            } catch (InterruptedException e) {
            }
          }
          regionFileName = REGION_LIST_USE_CNVS;

          // Select the file from the load recent menu
          regionFileNameBtn.get(REGION_LIST_USE_CNVS).setSelected(true);
        }

        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            cnvLabels = sampleData.getCnvClasses();
            setBounds(startX, startY, width, height);
            setPosition(startingLocation);
            updateSample(sample);
            setVisible(true);
            updateGUI();

            showRegion(0);
            System.out.println("Loaded regions in " + ext.getTimeElapsed(t));
          }
        });
      }
    }).start();
    // cnvLabels = sampleData.getCnvClasses();
    // updateSample(sample);
    // updateLocation(startingLocation);
    // setBounds(startX, startY, width, height);
    // setVisible(true);
    //
    // Trailer.this.regionFileName = REGION_LIST_USE_CNVS;
    // loadCNVsAsRegions();
    // updateCNVs(chr);
    // updateGUI();
    // if (startingLocation.equals(DEFAULT_LOCATION)) {
    // showRegion(0);
    // }

  }

  private void paintLRRPanel(Graphics g) {
    // TODO moving paintComponent code here breaks drawing, and I haven't figured out why.. (cole -
    // 3/6/15)
  }

  private void paintCNVGeneTrackPanel(Graphics g) {

    GeneData[] genes;
    int[][] exons;
    Vector<Segment> v = new Vector<Segment>();
    Segment[] segs;
    int width, begin, end, source;
    Segment currentView;
    String text;

    // g.drawRect(0, 0, this.getWidth()-1, this.getHeight()-1);


    if (track == null) {
      text = "Gene track is not installed";
      width = g.getFontMetrics(g.getFont()).stringWidth(text);
      g.drawString(text, getWidth() / 2 - width / 2, 10);
    } else {
      if (stop - start > 10000000) {
        g.drawString("Zoom in to see genes", 10, 10);
      } else {
        genes = track.getBetween(chr, start, stop, 30);
        g.setColor(Color.BLACK);
        for (GeneData gene : genes) {
          begin = getX(gene.getStart());
          end = getX(gene.getStop());
          g.drawRoundRect(begin, 0 * 15, end - begin, 10, 2, 2);
          v.add(new Segment(begin, end));
          exons = gene.getExonBoundaries();
          for (int j = 0; j < exons.length; j++) {
            begin = getX(exons[j][0]);
            end = getX(exons[j][1]);
            if (j == 0 || j == exons.length - 1) {
              g.fillRoundRect(begin, 0 * 15, end - begin + 1, 10, 2, 2);
            } else {
              g.fillRect(begin, 0 * 15, end - begin + 1, 10);
            }

          }
        }
        Segment.mergeOverlapsAndSort(v);
        segs = Segment.toArray(v);
        g.setFont(new Font("Arial", 0, 14));

        for (GeneData gene : genes) {
          begin = getX(gene.getStart());
          width = g.getFontMetrics(g.getFont()).stringWidth(gene.getGeneName());
          if (!Segment.overlapsAny(new Segment(begin - width - 5, begin - 1), segs)) {
            g.drawString(gene.getGeneName(), begin - width - 3, 0 * 15 + 10);
          }
        }
      }
    }
    currentView = new Segment(chr, start, stop);
    activeCNVs.clear();
    int firstBegin;
    for (int i = 0; cnvs != null && i < cnvs.length; i++) {
      source = i;
      firstBegin = Integer.MAX_VALUE;
      for (int j = 0; j < cnvs[i].length; j++) {
        if (cnvs[i][j].overlaps(currentView)) {
          activeCNVs.add(new int[] {i, j});
          begin = getX(cnvs[i][j].getStart());
          if (begin < firstBegin) {
            firstBegin = begin;
          }
          end = getX(cnvs[i][j].getStop());
          Color[] colors = getAColor(source);
          Color cnvColor =
                         cnvs[i][j].getCN() == PennHmm.LOH_FLAG ? Color.GRAY
                                                                : colors[cnvs[i][j].getCN() < 2 ? 0
                                                                                                : 1];
          g.setColor(cnvColor);
          g.fillRoundRect(begin, (source + 2) * 15, end - begin + 1, 10, 2, 2);
          g.setColor(Color.BLACK);
          if (selectedCNV != null && selectedCNV[0] == i && selectedCNV[1] == j) {
            g.drawRoundRect(begin - 1, (source + 2) * 15 - 1, end - begin + 2, 11, 5, 5);
          }
          // g.drawString(ext.rootOf(cnvLabels[source]), begin+2, (source+1)*15+10);
        }
      }
      if (firstBegin != Integer.MAX_VALUE) {
        g.drawString(ext.rootOf(cnvLabels[source]),
                     firstBegin - Grafik.getTextWidth(cnvLabels[source], g) - 3,
                     (source + 2) * 15 + 10);
      }
    }
  }

  private Color[] getAColor(int index) {
    while (index >= colorScheme.size()) {
      colorScheme.add(ColorExt.generatePastelShades());
    }
    return colorScheme.get(index);
  }

  private int getX(int pos) {
    return (int) ((double) (pos - start) / (double) (stop - start)
                  * (getWidth() - 2 * WIDTH_BUFFER))
           + WIDTH_BUFFER;
  }

  private void saveRegionFile() {
    JFileChooser jfc = new JFileChooser((proj != null ? proj.PROJECT_DIRECTORY.getValue() : ""));
    if (jfc.showSaveDialog(Trailer.this) == JFileChooser.APPROVE_OPTION) {
      String newFile = jfc.getSelectedFile().getAbsolutePath();
      String[] regionLines = new String[regions.length];
      for (int i = 0; i < regionLines.length; i++) {
        regionLines[i] = Array.toStr(regions[i], "\t");
      }
      Files.writeList(regionLines, newFile);
      addFileToList(newFile, true);
      regionFileName = newFile;

      JOptionPane.showMessageDialog(Trailer.this, "Successfully saved " + regions.length
                                                  + " regions to file: " + newFile);
    }
  }

  private String chooseNewFiles() {
    JFileChooser jfc =
                     new JFileChooser((proj != null
                                       || regionFileName == null ? proj.PROJECT_DIRECTORY.getValue()
                                                                 : ext.parseDirectoryOfFile(regionFileName)));
    jfc.setMultiSelectionEnabled(true);
    if (jfc.showOpenDialog(Trailer.this) == JFileChooser.APPROVE_OPTION) {
      File[] files = jfc.getSelectedFiles();
      if (files.length > 0) {
        boolean[] keep = Array.booleanArray(files.length, true);
        for (int i = 0; i < files.length; i++) {
          for (String fileName : regionFileNameLoc.keySet()) {
            if (ext.rootOf(files[i].toString()).equals(fileName)) {
              keep[i] = false;
            }
          }
        }
        File[] keptFiles = Array.subArray(files, keep);
        File[] discards = Array.subArray(files, Array.booleanNegative(keep));

        if (discards.length > 0) {
          StringBuilder msg = new StringBuilder("The following data file(s) are already present:");
          for (File disc : discards) {
            msg.append("\n").append(disc.getName());
          }
          JOptionPane.showMessageDialog(Trailer.this, msg.toString());
        }

        for (File kept : keptFiles) {
          addFileToList(kept.getAbsolutePath());
        }
        return keptFiles[0].getAbsolutePath();
      } else {
        File file = jfc.getSelectedFile();
        boolean keep = true;
        for (String fileName : regionFileNameLoc.keySet()) {
          if (ext.rootOf(file.toString()).equals(fileName)) {
            keep = false;
          }
        }

        if (!keep) {
          StringBuilder msg =
                            new StringBuilder("The following data file is already present:\n").append(file.getName());
          JOptionPane.showMessageDialog(Trailer.this, msg.toString());
          return null;
        } else {
          addFileToList(file.getAbsolutePath());
          return file.getAbsolutePath();
        }

      }

    }
    return null;
  }

  private boolean addFileToList(String rawfile) {
    return addFileToList(rawfile, false);
  }

  /**
   * Adds a file to the "load recent file" menu
   *
   * @param rawfile file to add
   * @param selectFile If true, the file will also be selected in the UI
   * @return true if added successfully
   */
  private boolean addFileToList(String rawfile, boolean selectFile) {
    String file = ext.verifyDirFormat(rawfile);
    file = file.substring(0, file.length() - 1);
    String name = ext.rootOf(file);
    if (regionFileNameLoc.get(name) != null) {
      return false;
    }
    regionFileNameLoc.put(name, file);

    JCheckBoxMenuItem item = new JCheckBoxMenuItem();
    item.setAction(markerFileSelectAction);
    item.setText(name);
    item.setSelected(selectFile);

    regionFileNameBtn.put(name, item);
    regionButtonGroup.add(item);
    loadRecentFileMenu.add(item);
    return true;
  }

  public void generateComponents() {
    JPanel dataPanel = new JPanel();
    dataPanel.setLayout(new GridLayout(3, 1, 5, 5));

    lrrPanel = new JPanel() {
      public static final long serialVersionUID = 2L;

      @Override
      public void paintComponent(Graphics g) {
        paintLRRPanel(g);

        float min, max;

        if (lrrValues != null) {
          min = Array.min(lrrValues);
          max = Array.max(lrrValues);

          // System.out.println("Displaying "+(stopMarker-startMarker)+"
          // markers");
          if (stopMarker - startMarker < DYNAMIC_HEIGHT_LIMIT) {
            min = Float.POSITIVE_INFINITY;
            max = Float.NEGATIVE_INFINITY;
            for (int i = startMarker; i <= stopMarker; i++) {
              if (lrrValues[i] < min) {
                min = lrrValues[i];
              }
              if (lrrValues[i] > max) {
                max = lrrValues[i];
              }
            }
            min = (float) Math.floor(min);
            max = (float) Math.ceil(max);
          }

          // min = -1.5;
          // max = 1.5;

          switch (transformation_type) {
            case -1:
            case 0:
              // Illumina
              min = -3;
              max = 3;
              // Agilent
              // min = -1.5f;
              // max = 1.5f;
              break;
            case 1:
              min = 0;
              max = 1;
              break;
            case 2:
              min = -8;
              max = 8;
              break;
            case 3:
              min = -12;
              max = 12;
              break;

          }

          g.setFont(new Font("Arial", 0, 20));
          g.drawString("Log R Ratio", WIDTH_BUFFER, 20);

          if (SHOW_MIDLINE) {
            g.setColor(Color.LIGHT_GRAY);
            g.drawLine(WIDTH_BUFFER,
                       getHeight() - (int) ((double) (0 - min) / (double) (max - min)
                                            * (getHeight() - 2 * HEIGHT_BUFFER))
                                     - HEIGHT_BUFFER,
                       getWidth() - WIDTH_BUFFER, getHeight()
                                                  - (int) ((double) (0 - min) / (double) (max - min)
                                                           * (getHeight() - 2 * HEIGHT_BUFFER))
                                                  - HEIGHT_BUFFER);
          }

          g.setFont(new Font("Arial", 0, 12));
          for (int i = startMarker; i <= stopMarker; i++) {
            // if (genotypes[i] == 1) {

            if (bafs != null && bafs.length > i && bafs[i] > 0.2 && bafs[i] < 0.8) {
              g.setColor(Color.RED);
              // colorScheme[2]
            } else {
              g.setColor(Color.BLACK);
            }
            if (currentColorManager != null && currentColorManager.hasColorFor(markerNames[i])) {// TODO,
                                                                                                 // all
                                                                                                 // logic
                                                                                                 // with
              Color managed = currentColorManager.getColorItemForVar(markerNames[i]).getColor();
              if (managed != null) {
                g.setColor(managed);
              }
            }
            if (!Float.isNaN(lrrValues[i])) {
              if (dropped[i]
                  || (excludeManager != null && !excludeManager.hasColorFor(markerNames[i]))) {
                // g.drawString("X", getX(positions[i]),
                // getHeight()-(int)((double)(lrrValues[i]-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
              } else if (lrrValues[i] < min) {
                g.drawString("v",
                             Trailer.this.getX(positions[i]), getHeight()
                                                              - (int) ((double) (min - min)
                                                                       / (double) (max - min)
                                                                       * (getHeight()
                                                                          - 2 * HEIGHT_BUFFER))
                                                              - HEIGHT_BUFFER);
              } else if (lrrValues[i] > max) {
                g.drawString("^",
                             Trailer.this.getX(positions[i]), getHeight()
                                                              - (int) ((double) (max - min)
                                                                       / (double) (max - min)
                                                                       * (getHeight()
                                                                          - 2 * HEIGHT_BUFFER))
                                                              - HEIGHT_BUFFER);
              } else {
                g.fillOval(Trailer.this.getX(positions[i]), getHeight()
                                                            - (int) ((double) (lrrValues[i] - min)
                                                                     / (double) (max - min)
                                                                     * (getHeight()
                                                                        - 2 * HEIGHT_BUFFER))
                                                            - HEIGHT_BUFFER,
                           SIZE, SIZE);
              }
            }
          }
        }
      }
    };
    lrrPanel.addMouseListener(this);
    lrrPanel.addMouseMotionListener(this);
    lrrPanel.addMouseWheelListener(this);
    dataPanel.add(lrrPanel);

    cnvPanel = new JPanel() {
      public static final long serialVersionUID = 8L;

      @Override
      public void paintComponent(Graphics g) {
        paintCNVGeneTrackPanel(g);
      }
    };
    cnvPanel.setMaximumSize(new Dimension(getWidth(), 20));
    cnvPanel.addMouseListener(cnvAdapter);
    cnvPanel.addMouseMotionListener(cnvAdapter);
    dataPanel.add(cnvPanel);

    bafPanel = new JPanel() {
      public static final long serialVersionUID = 7L;

      @Override
      public void paintComponent(Graphics g) {
        g.setFont(new Font("Arial", 0, 20));
        g.drawString("B Allele Frequency", WIDTH_BUFFER, 20);

        g.setFont(new Font("Arial", 0, 10));
        if (lrrs != null) {
          for (int i = startMarker; i <= stopMarker; i++) {
            if (!Float.isNaN(lrrs[i])) {
              if (dropped[i]
                  || (excludeManager != null && !excludeManager.hasColorFor(markerNames[i]))) {
                if (excludeManager != null && !excludeManager.hasColorFor(markerNames[i])) {

                } else {
                  g.drawString("X",
                               Trailer.this.getX(positions[i]), getHeight()
                                                                - (int) (bafs[i]
                                                                         * (double) (getHeight() - 2
                                                                                                   * HEIGHT_BUFFER))
                                                                - HEIGHT_BUFFER + 5);
                }
              } else if (genotypes != null && genotypes[i] == -1) {
                g.drawString("+",
                             Trailer.this.getX(positions[i]), getHeight()
                                                              - (int) (bafs[i]
                                                                       * (double) (getHeight() - 4
                                                                                                 * HEIGHT_BUFFER))
                                                              - HEIGHT_BUFFER / 2);
              } else if (bafs != null && bafs.length > i) {
                g.fillOval(Trailer.this.getX(positions[i]), getHeight()
                                                            - (int) (bafs[i]
                                                                     * (double) (getHeight() - 4
                                                                                               * HEIGHT_BUFFER))
                                                            - HEIGHT_BUFFER,
                           SIZE, SIZE);
              } else {
                g.setFont(new Font("Arial", 0, 20));
                g.drawString("Error - no BAF values present for sample.", (getWidth() / 2) - 175,
                             (getHeight() / 2) - 10);
              }
            }
          }
        }
      }
    };
    dataPanel.add(bafPanel);

    getContentPane().add(dataPanel, BorderLayout.CENTER);

    JPanel sampPanel = new JPanel();
    ((FlowLayout) sampPanel.getLayout()).setVgap(0);
    firstRegion = new JButton(Grafik.getImageIcon("images/firstLast/First.gif"));
    firstRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif"));
    firstRegion.addActionListener(this);
    firstRegion.setActionCommand(FIRST_REGION);
    firstRegion.setPreferredSize(new Dimension(25, 25));
    firstRegion.setToolTipText("Go to first region");

    previousRegion = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif"));
    previousRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif"));
    previousRegion.addActionListener(this);
    previousRegion.setActionCommand(PREVIOUS_REGION);
    previousRegion.setPreferredSize(new Dimension(25, 25));
    previousRegion.setToolTipText("Go to previous region");

    sampleList = new JComboBox();
    sampleList.setFont(new Font("Arial", 0, 20));
    createSampleList();

    DefaultListCellRenderer dlcr = new DefaultListCellRenderer();
    dlcr.setHorizontalAlignment(DefaultListCellRenderer.CENTER);
    sampleList.setRenderer(dlcr);
    sampleList.setBorder(BorderFactory.createEtchedBorder());
    sampleList.setEditable(false);
    sampleList.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        @SuppressWarnings("unchecked")
        JComboBox jcb = (JComboBox) e.getSource();
        int index = jcb.getSelectedIndex();
        if (index == samplesPresent.length - 1) {
          createSampleList();
        } else if (sample != samplesPresent[index]) {
          updateSample(samplesPresent[index]);
        }
      }
    });
    sampPanel.add(sampleList);


    JPanel descrPanel = new JPanel();
    descrPanel.setLayout(new MigLayout("gap 0", "[grow, center]", "[]0[]0[]"));

    nextRegion = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif"));
    nextRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif"));
    nextRegion.addActionListener(this);
    nextRegion.setActionCommand(NEXT_REGION);
    nextRegion.setPreferredSize(new Dimension(25, 25));
    nextRegion.setToolTipText("Go to next region");

    lastRegion = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif"));
    lastRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif"));
    lastRegion.addActionListener(this);
    lastRegion.setActionCommand(LAST_REGION);
    lastRegion.setPreferredSize(new Dimension(25, 25));
    lastRegion.setToolTipText("Go to last region");
    sampPanel.setPreferredSize(new Dimension(sampPanel.getPreferredSize().width,
                                             sampleList.getPreferredSize().height + 5));
    descrPanel.add(sampPanel, "cell 0 0");

    JPanel compPanel = new JPanel(new MigLayout("align center, fill, gap 0", "[grow, center]",
                                                "[]5[21:21:21]5[]"));

    JPanel regionPanel = new JPanel();
    ((FlowLayout) regionPanel.getLayout()).setVgap(0);
    regionField = new JTextField("", 8);
    regionField.setHorizontalAlignment(JTextField.CENTER);
    Font font = new Font("Arial", 0, 14);
    regionField.setFont(font);
    regionField.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        try {
          int trav =
                   Integer.valueOf(((JTextField) e.getSource()).getText().trim().split("[\\s]+")[0])
                          .intValue()
                     - 1;
          if (trav >= 0 && trav < regions.length) {
            showRegion(trav);
          }
        } catch (NumberFormatException nfe) {
        }
        updateGUI();
      }
    });
    regionPanel.add(firstRegion);
    regionPanel.add(previousRegion);
    regionPanel.add(regionField);
    regionField.setPreferredSize(new Dimension(regionField.getPreferredSize().width, 26));
    regionPanel.add(nextRegion);
    regionPanel.add(lastRegion);
    compPanel.add(regionPanel, "cell 0 0");

    commentLabel = new JLabel(" ", JLabel.CENTER);
    commentLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
    commentLabel.setFont(font);
    commentLabel.setToolTipText("Click to Edit Comment");
    commentLabel.addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        super.mouseClicked(e);
        commentLabel.setVisible(false);
        commentField.setVisible(true);
        commentField.requestFocusInWindow();
      }
    });
    compPanel.add(commentLabel, "cell 0 1, hidemode 3");

    commentField = new JTextField(20);
    commentField.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        // When enter is pressed, request focus elsewhere to trigger the focus lost event
        // of the comment field.
        Trailer.this.requestFocusInWindow();
      }
    });
    commentField.setAlignmentX(Component.CENTER_ALIGNMENT);
    commentField.setFont(font);
    commentField.addFocusListener(new FocusAdapter() {
      @Override
      public void focusLost(FocusEvent e) {
        super.focusLost(e);
        String newComment = commentField.getText();
        if (regions != null && regions.length > 0) {
          String[] regionDetails = regions[regionIndex];
          if (regionDetails.length >= 3) {
            regionDetails[2] = newComment;
          } else {
            regionDetails = Array.addStrToArray(newComment, regionDetails);
          }
          regions[regionIndex] = regionDetails;
          commentLabel.setText(regions[regionIndex][2].isEmpty() ? BLANK_COMMENT
                                                                 : "region #" + (regionIndex + 1)
                                                                   + ":  "
                                                                   + regions[regionIndex][2]);
          promptCommentSave = promptAndSaveRegions(promptCommentSave);
        }
        commentLabel.setVisible(true);
        commentField.setVisible(false);
      }
    });
    commentField.setVisible(false);
    compPanel.add(commentField, "cell 0 1, hidemode 3");

    qcLabel = new JLabel(" ", JLabel.CENTER);
    qcLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
    qcLabel.setFont(font);
    // qcLabel.setToolTipText("Click for More Details");
    compPanel.add(qcLabel, "cell 0 2");

    descrPanel.add(compPanel, "cell 0 1");
    compPanel.setPreferredSize(new Dimension(compPanel.getPreferredSize().width, 95));

    JPanel navigateChrPanel = new JPanel();

    navPanel = new RegionNavigator(this);
    descrPanel.add(navPanel, "cell 0 2");

    JPanel overPanel = new JPanel();
    overPanel.setLayout(new BoxLayout(overPanel, BoxLayout.LINE_AXIS));

    JSeparator sep = new JSeparator(SwingConstants.VERTICAL);
    sep.setMaximumSize(new Dimension(1, 150));
    overPanel.setBorder(new EmptyBorder(5, 0, 5, 0));

    overPanel.add(Box.createHorizontalGlue());
    overPanel.add(descrPanel);
    overPanel.add(Box.createHorizontalGlue());

    getContentPane().add(overPanel, BorderLayout.NORTH);

    InputMap inputMap = bafPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS_REGION);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT_REGION);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, InputEvent.CTRL_MASK), RegionNavigator.PREVIOUS_CHR);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, InputEvent.CTRL_MASK), RegionNavigator.NEXT_CHR);

    //FIXME this is a duplication of the action listening done by trailer..
    ActionMap actionMap = bafPanel.getActionMap();
    actionMap.put(RegionNavigator.FIRST_CHR, new AbstractAction() {
      public static final long serialVersionUID = 5L;

      @Override
      public void actionPerformed(ActionEvent e) {
        setPosition("chr1");
      }
    });
    actionMap.put(RegionNavigator.PREVIOUS_CHR, new AbstractAction() {
      public static final long serialVersionUID = 6L;

      @Override
      public void actionPerformed(ActionEvent e) {
        setPosition("chr" + Math.max(chr - 1, 1));
      }
    });
    actionMap.put(RegionNavigator.NEXT_CHR, new AbstractAction() {
      public static final long serialVersionUID = 7L;

      @Override
      public void actionPerformed(ActionEvent e) {
        setPosition("chr" + Math.min(chr + 1, 26));
      }
    });
    actionMap.put(RegionNavigator.LAST_CHR, new AbstractAction() {
      public static final long serialVersionUID = 8L;

      @Override
      public void actionPerformed(ActionEvent e) {
        setPosition("chr26");
      }
    });
    actionMap.put(FIRST_REGION, new AbstractAction() {
      public static final long serialVersionUID = 9L;

      @Override
      public void actionPerformed(ActionEvent e) {
        showRegion(0);
      }
    });
    actionMap.put(PREVIOUS_REGION, new AbstractAction() {
      public static final long serialVersionUID = 10L;

      @Override
      public void actionPerformed(ActionEvent e) {
        showRegion(Math.max(regionIndex - 1, 0));
      }
    });
    actionMap.put(NEXT_REGION, new AbstractAction() {
      public static final long serialVersionUID = 11L;

      @Override
      public void actionPerformed(ActionEvent e) {
        showRegion(Math.min(regionIndex + 1, regions.length - 1));
      }
    });
    actionMap.put(LAST_REGION, new AbstractAction() {
      public static final long serialVersionUID = 12L;

      @Override
      public void actionPerformed(ActionEvent e) {
        showRegion(regions.length - 1);
      }
    });
    bafPanel.setActionMap(actionMap);

    // doesn't seem to get captured properly...
    //FIXME move this to RegionNavigator if needed
//    nextChr.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT_REGION);
//    nextChr.setActionMap(actionMap);
//    previousChr.setActionMap(actionMap);
  }

  private void doScreenCapture(String filename) {
    BufferedImage cap = generateScreenshot();

    if (filename == null) {
      TransferableImage ti = new TransferableImage(cap);
      Clipboard c = Toolkit.getDefaultToolkit().getSystemClipboard();
      c.setContents(ti, null);
    } else {
      try {
        ImageIO.write(cap, "png", new File(filename));
      } catch (IOException e) {
        if (proj != null) {
          proj.getLog().reportException(e);
          proj.message("Error occured while writing screen capture to file.  Please check log for more details.");
        }
      }
    }
  }


  private BufferedImage generateScreenshot() {
    int lW = lrrPanel.getWidth();
    int bW = bafPanel.getWidth();
    int cW = cnvPanel.getWidth();
    int lH = lrrPanel.getHeight();
    int bH = bafPanel.getHeight();
    int cH = cnvPanel.getHeight();
    BufferedImage imageLrr = new BufferedImage(lW, lH, BufferedImage.TYPE_INT_RGB);
    BufferedImage imageBaf = new BufferedImage(bW, bH, BufferedImage.TYPE_INT_RGB);
    BufferedImage imageCnv = new BufferedImage(cW, cH, BufferedImage.TYPE_INT_RGB);

    Graphics g = imageLrr.getGraphics();
    g.setColor(lrrPanel.getBackground());
    g.fillRect(0, 0, imageLrr.getWidth(), imageLrr.getHeight());
    lrrPanel.paintAll(g);

    g = imageCnv.getGraphics();
    g.setColor(cnvPanel.getBackground());
    g.fillRect(0, 0, imageCnv.getWidth(), imageCnv.getHeight());
    cnvPanel.paintAll(g);
    if (selectedCNV != null) {
      // CNVariant cnv = cnvs[selectedCNV[0]][selectedCNV[1]];
      // TODO include CNV details in ScreenCapture?
    }

    g = imageBaf.getGraphics();
    g.setColor(bafPanel.getBackground());
    g.fillRect(0, 0, imageBaf.getWidth(), imageBaf.getHeight());
    bafPanel.paintAll(g);

    int w = Math.max(lW, Math.max(bW, cW));
    int h = lH + bH + cH + HEIGHT_BUFFER;
    BufferedImage finalImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);

    g = finalImage.getGraphics();
    g.setColor(lrrPanel.getBackground());
    g.fillRect(0, 0, finalImage.getWidth(), finalImage.getHeight());
    g.drawImage(imageLrr, 0, 0, null);
    g.drawImage(imageCnv, 0, lH + 5, null);
    g.drawImage(imageBaf, 0, lH + cH + 10, null);
    g.dispose();

    return finalImage;
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
        NewRegionListDialog newRgnList = new NewRegionListDialog(
                                                                 proj == null ? null
                                                                              : proj.getSamples(),
                                                                 proj == null ? null
                                                                              : proj.PROJECT_DIRECTORY.getValue(),
                                                                 true);
        newRgnList.setModal(true);
        newRgnList.setVisible(true);
        if (newRgnList.getReturnCode() == JOptionPane.YES_OPTION) {
          String rgnFile = newRgnList.getFileName();
          addFileToList(rgnFile);
          String file = ext.verifyDirFormat(rgnFile);
          file = file.substring(0, file.length() - 1);
          String name = ext.rootOf(file);
          regionFileNameBtn.get(name).setSelected(true);
          regionFileNameBtn.get(name).doClick();
        }
      }
    });
    newRegionFile.setText(REGION_LIST_NEW_FILE);
    newRegionFile.setMnemonic(KeyEvent.VK_N);
    Font font = new Font("Arial", 0, 12);
    newRegionFile.setFont(font);
    fileMenu.add(newRegionFile);

    JMenuItem addRegion = new JMenuItem();
    addRegion.setAction(addRegionAction);
    addRegion.setText(REGION_LIST_ADD);
    addRegion.setFont(font);
    fileMenu.add(addRegion);

    JMenuItem saveRegionFile = new JMenuItem();
    saveRegionFile.setAction(saveRegionFileAction);
    saveRegionFile.setText(REGION_LIST_SAVE_FILE);
    saveRegionFile.setMnemonic(KeyEvent.VK_A);
    saveRegionFile.setFont(font);
    fileMenu.add(saveRegionFile);

    JMenuItem loadRegionFile = new JMenuItem();
    loadRegionFile.setAction(loadRegionFileAction);
    loadRegionFile.setText(REGION_LIST_LOAD_FILE);
    loadRegionFile.setMnemonic(KeyEvent.VK_L);
    loadRegionFile.setFont(font);
    fileMenu.add(loadRegionFile);
    loadRecentFileMenu = new JMenu("Load Recent Region List...");
    loadRecentFileMenu.setMnemonic(KeyEvent.VK_R);
    loadRecentFileMenu.setFont(font);
    fileMenu.add(loadRecentFileMenu);

    JMenuItem screencap1 = new JMenuItem();
    screencap1.setAction(screencapAction);
    screencap1.setMnemonic(KeyEvent.VK_S);
    screencap1.setText("Screen Capture");
    screencap1.setFont(font);
    fileMenu.add(screencap1);

    JMenuItem screencap2 = new JMenuItem();
    screencap2.setAction(screencapClipboardAction);
    screencap2.setMnemonic(KeyEvent.VK_C);
    screencap2.setText("Screen Capture to Clipboard");
    screencap2.setFont(font);
    fileMenu.add(screencap2);

    menuBar.add(fileMenu);

    regionButtonGroup = new ButtonGroup();
    if (proj != null) {
      String[] files = proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue();
      String name;
      for (String file : files) {
        name = ext.rootOf(file);
        regionFileNameLoc.put(name, file);
        JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
        menuItem.setAction(markerFileSelectAction);
        menuItem.setFont(font);
        boolean found = Files.exists(file);
        menuItem.setText(name + (found ? "" : " -- [file not found]"));
        menuItem.setEnabled(found);
        regionFileNameBtn.put(name, menuItem);
        regionButtonGroup.add(menuItem);
        loadRecentFileMenu.add(menuItem);
      }
    }

    JCheckBoxMenuItem item1 = new JCheckBoxMenuItem();
    item1.setAction(markerFileSelectAction);
    item1.setText(REGION_LIST_USE_CNVS);
    item1.setFont(font);
    regionFileNameBtn.put(REGION_LIST_USE_CNVS, item1);
    regionButtonGroup.add(item1);
    loadRecentFileMenu.add(item1);

    final JRadioButtonMenuItem[] transformBtns =
                                               new JRadioButtonMenuItem[Transforms.TRANFORMATIONS.length];
    gcCorrectButton = new JCheckBoxMenuItem(GcAdjustor.GC_ADJUSTOR_TITLE[0], false);// stays hidden
                                                                                    // if gcModel is
                                                                                    // not detected

    {
      JMenu transMenu = new JMenu("Transformation");

      JMenuItem lbl1 = transMenu.add("LRR Transforms:");
      lbl1.setEnabled(false);
      lbl1.setFont(font);

      ButtonGroup lrrBtnGrp = new ButtonGroup();
      ButtonGroup transBtnGrp = new ButtonGroup();
      ItemListener typeListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JRadioButtonMenuItem jrb = (JRadioButtonMenuItem) ie.getItem();
          if (jrb.isSelected()) {
            centroids = null;
            bafs = originalBAFs;
            for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
              if (jrb.getText().equals(Transforms.TRANFORMATIONS[i])) {
                transformation_type = i;
                lrrValues = getNewLRRs(proj, lrrs, transformation_type,
                                       transformSeparatelyByChromosome, markerSet, gcModel,
                                       gcCorrectButton.isSelected(), true, log);
                updateGUI();
              }
            }
          }
        }
      };
      HashSet<String> mnems = new HashSet<String>();
      for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
        transformBtns[i] = new JRadioButtonMenuItem(Transforms.TRANFORMATIONS[i]);
        transformBtns[i].addItemListener(typeListener);
        transformBtns[i].setFont(font);
        String[] wds = Transforms.TRANFORMATIONS[i].split("[\\s]+");
        int ind = 0;
        String mnem = wds[ind].substring(0, 1);
        while (mnems.contains(mnem)) {
          mnem = wds[++ind].substring(0, 1);
        }
        mnems.add(mnem);
        transformBtns[i].setMnemonic(mnem.charAt(0));
        lrrBtnGrp.add(transformBtns[i]);
        transMenu.add(transformBtns[i]);
      }
      transformBtns[0].setSelected(true);

      transMenu.addSeparator();
      JMenuItem lbl2 = transMenu.add("Transform By:");
      lbl2.setEnabled(false);
      lbl2.setFont(font);

      JRadioButtonMenuItem[] scopeBtns = new JRadioButtonMenuItem[Transforms.SCOPES.length];
      ItemListener scopeListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JRadioButtonMenuItem jrb = (JRadioButtonMenuItem) ie.getItem();
          if (jrb.isSelected()) {
            transformSeparatelyByChromosome = jrb.getText().equals(Transforms.SCOPES[1]);
            lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome,
                                   markerSet, gcModel, gcCorrectButton.isSelected(), true, log);
            updateGUI();
          }
        }
      };
      mnems = new HashSet<String>();
      for (int i = 0; i < Transforms.SCOPES.length; i++) {
        scopeBtns[i] = new JRadioButtonMenuItem(Transforms.SCOPES[i]);
        scopeBtns[i].addItemListener(scopeListener);
        scopeBtns[i].setFont(font);

        String[] wds = Transforms.SCOPES[i].split("[\\s]+");
        int ind = 0;
        String mnem = wds[ind].substring(0, 1);
        while (mnems.contains(mnem)) {
          mnem = wds[++ind].substring(0, 1);
        }
        mnems.add(mnem);
        scopeBtns[i].setMnemonic(mnem.charAt(0));
        transBtnGrp.add(scopeBtns[i]);
        transMenu.add(scopeBtns[i]);
      }
      scopeBtns[0].setSelected(true);

      menuBar.add(transMenu);
    }
    {
      gcParameterDisplay = new GCParameterDisplay(proj, this, proj.getLog());

      JMenu adjMenu = new JMenu("Adjustments");
      adjMenu.setMnemonic(KeyEvent.VK_J);
      ItemListener gcListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          if (jrb.isSelected()) {
            gcParameterDisplay.getButtonMap().get("None").setSelected(true);
          }
          lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome,
                                 markerSet, gcModel, jrb.isSelected(), jrb.isSelected(), log);
          updateGUI();
        }
      };

      gcCorrectButton.setToolTipText("GC correction will be applied prior to any transformation");
      gcCorrectButton.addItemListener(gcListener);
      gcCorrectButton.setFont(font);
      // act.addSeparator();
      adjMenu.add(gcCorrectButton).setEnabled(gcModel != null);
      // adjMenu.addSeparator();
      adjMenu.add(gcParameterDisplay.getActionMenu());

      menuBar.add(adjMenu);
    }

    {
      JMenu qcMenu = new JMenu("Show QC");
      String[] opts = new String[] {"Hide QC", "Genome", /* "Chromosome", */ "Region"};
      ItemListener qcListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent arg0) {
          if (arg0.getStateChange() == ItemEvent.SELECTED) {
            String cmd = ((AbstractButton) arg0.getSource()).getActionCommand();
            if ("Hide QC".equals(cmd)) {
              qcSelection = 0;
            } else if ("Genome".equals(cmd)) {
              qcSelection = 1;
              // } else if ("Chromosome".equals(cmd)) {
              // qcSelection = 2;
            } else if ("Region".equals(cmd)) {
              qcSelection = 2;
            }
            updateQCDisplay();
          }
        }
      };
      ButtonGroup qcBtnGrp = new ButtonGroup();
      for (int i = 0; i < opts.length; i++) {
        JRadioButtonMenuItem qcButton = new JRadioButtonMenuItem(opts[i]);
        qcButton.setActionCommand(opts[i]);
        qcButton.addItemListener(qcListener);
        qcButton.setFont(font);
        qcBtnGrp.add(qcButton);
        qcMenu.add(qcButton);
        if (i == 0) {
          qcButton.setSelected(true);
        }
      }
      menuBar.add(qcMenu);

    }
    // TODO


    {
      JMenu centMenu = new JMenu("Centroids");
      autoSwitch = new JCheckBoxMenuItem();
      autoSwitch.setText("Auto-Select Sex Centroid by Sample Sex");
      autoSwitch.setMnemonic(KeyEvent.VK_A);
      autoSwitch.setFont(font);
      centMenu.add(autoSwitch);

      ItemListener centListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          if (jrb.isSelected() && !"None".equals(jrb.getText())) {
            if (namePathMap == null || namePathMap.isEmpty()) {
              jrb.setSelected(false);
              if (transformBtns != null && transformBtns[0] != null) {
                transformBtns[0].setSelected(true);
              }
              return;
            }
            if (!isSettingCentroid) {
              isSettingCentroid = true;
              transformation_type = -1;
              setCentroid(jrb.getText());
              loadValues();
              updateGUI();
              repaint();
              isSettingCentroid = false;
            } else {
              setCentroid(jrb.getText());
            }
          } else {
            currentCentroid = null;
            centroids = null;
            bafs = originalBAFs;
            loadValues();
            updateGUI();
            repaint();
          }
        }
      };
      namePathMap = new Hashtable<String, String>();
      Vector<String> centFiles = new Vector<String>();
      centFiles.add(proj.ORIGINAL_CENTROIDS_FILENAME.getValue());
      centFiles.add(proj.GENOTYPE_CENTROIDS_FILENAME.getValue());
      centFiles.add(proj.CUSTOM_CENTROIDS_FILENAME.getValue());
      centFiles.add(proj.CHIMERA_CENTROIDS_FILENAME.getValue());

      String tempFileMale = proj.SEX_CENTROIDS_MALE_FILENAME.getValue();
      String tempFileFemale = proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue();
      if (tempFileMale != null && tempFileFemale != null && !"".equals(tempFileMale)
          && !"".equals(tempFileFemale)) {
        centFiles.add(SEX_CENT);
      }

      for (String file : centFiles) {
        if (SEX_CENT.equals(file)) {
          namePathMap.put(SEX_CENT + " - Male", tempFileMale);
          namePathMap.put(SEX_CENT + " - Female", tempFileFemale);
        } else if (Files.exists(file)) {
          String name = file.substring(file.lastIndexOf("/") + 1);
          name = name.substring(0, name.lastIndexOf("."));
          namePathMap.put(name, file);
        }
      }

      centButtonMap = new HashMap<String, JCheckBoxMenuItem>();
      ButtonGroup centButtons = new ButtonGroup();
      JCheckBoxMenuItem centBox = new JCheckBoxMenuItem("None");
      centBox.addItemListener(centListener);
      centBox.setFont(font);
      centBox.setSelected(true);
      centButtons.add(centBox);

      JMenu lbl3 = new JMenu("Derive from Centroids");
      lbl3.setMnemonic(KeyEvent.VK_D);
      lbl3.setFont(font);
      centMenu.add(lbl3);
      lbl3.add(centBox);
      String[] centKeys = namePathMap.keySet().toArray(new String[] {});
      for (String key : centKeys) {
        centBox = new JCheckBoxMenuItem(key);
        centBox.addItemListener(centListener);
        centBox.setFont(font);
        centButtons.add(centBox);
        lbl3.add(centBox);
        centButtonMap.put(key, centBox);
      }

      menuBar.add(centMenu);
    }

    {
      // act.addSeparator();
      JMenu cnvMenu = new JMenu("CNV Calls");

      ItemListener cnvListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          if (pennHmm == null) {
            if (Files.exists(proj.HMM_FILENAME.getValue())) {
              pennHmm = PennHmm.loadPennHmm(proj.HMM_FILENAME.getValue(), proj.getLog());

            } else {
              pennHmm = null;
            }
          }
          if (pfb == null) {
            if (Files.exists(proj.CUSTOM_PFB_FILENAME.getValue())) {
              pfb = PFB.loadPFB(proj);

            } else {
              pfb = null;
            }
          }
          if (pfb == null || pennHmm == null) {
            if (pennHmm == null) {
              proj.getLog().reportTimeError("Could not load " + proj.HMM_FILENAME.getName()
                                            + " defined by " + proj.HMM_FILENAME.getValue());
            } else if (pfb == null) {
              proj.getLog().reportTimeError("Could not load " + proj.CUSTOM_PFB_FILENAME.getName()
                                            + " defined by " + proj.CUSTOM_PFB_FILENAME.getValue());
            }
          } else {
            CNVCallResult callResult = CNVCaller.callCNVsFor(proj, pennHmm, sample,
                                                             Array.toDoubleArray(lrrs),
                                                             Array.toDoubleArray(bafs), gcModel,
                                                             pfb, markerSet, new int[] {chr}, null,
                                                             false, CNVCaller.DEFUALT_MIN_SITES,
                                                             CNVCaller.DEFUALT_MIN_CONF,
                                                             PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT,
                                                             proj.NUM_THREADS.getValue(), true);
            int externalCNVs = prepInternalClasses();
            addCnvsToPheno(callResult.getChrCNVs().getLoci(), externalCNVs,
                           INTERNAL_CNV_TYPES.CNV_CALLER);
            // addCnvsToPheno(callResult.getChrCNVsReverse().getLoci(),
            // externalCNVs,INTERNAL_CNV_TYPES.REV_CNV_CALLER);
            // addCnvsToPheno(callResult.getChrCNVsReverseConsensus().getLoci(), externalCNVs,
            // INTERNAL_CNV_TYPES.CONSENSUS);
            sampleData.getSampleHash().put(sample.toLowerCase(), indiPheno);
            updateCNVs(chr);
          }
          jrb.setSelected(false);
          updateGUI();
        }
      };
      callCnvsButton = new JCheckBoxMenuItem("Call CNVs", false);// stays hidden if gcModel is not
                                                                 // detected
      callCnvsButton.addItemListener(cnvListener);
      callCnvsButton.setFont(font);
      cnvMenu.add(callCnvsButton);
      // cnvMenu.addSeparator();

      ItemListener mosaicListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          MosaicBuilder builder = new MosaicBuilder();
          builder.verbose(true);
          MosaicismDetect md = builder.build(proj, sample, markerSet, Array.toDoubleArray(bafs));
          Segment seg = new Segment(chr, 0, Integer.MAX_VALUE);
          LocusSet<MosaicRegion> mosSet = md.callMosaic(seg, false);
          int externalCNVs = prepInternalClasses();
          addCnvsToPheno(mosSet.getLoci(), externalCNVs, INTERNAL_CNV_TYPES.MOSAIC_CALLER);
          sampleData.getSampleHash().put(sample.toLowerCase(), indiPheno);
          updateCNVs(chr);

          jrb.setSelected(false);
          updateGUI();
        }
      };
      mosaicCallButton = new JCheckBoxMenuItem("Call Mosaicism (extra beta)", false);
      mosaicCallButton.addItemListener(mosaicListener);
      mosaicCallButton.setFont(font);
      cnvMenu.add(mosaicCallButton);

      ItemListener mosaicFListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          Segment quantSeg = new Segment(chr, positions[startMarker], positions[stopMarker]);
          quantHere(quantSeg, false);
          jrb.setSelected(false);
        }



      };
      mosaicFButton = new JCheckBoxMenuItem("Quantify Mosaicism (extra extra beta)", false);
      mosaicFButton.addItemListener(mosaicFListener);
      mosaicFButton.setFont(font);
      cnvMenu.add(mosaicFButton);

      ItemListener beastHListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          Segment quantSeg = new Segment(chr, positions[startMarker], positions[stopMarker]);
          beastHere(quantSeg, INTERNAL_CNV_TYPES.BEAST_SCORE_CUSTOM);
          jrb.setSelected(false);
        }

      };

      beastHButton = new JCheckBoxMenuItem("Compute BEAST Height", false);
      beastHButton.addItemListener(beastHListener);
      beastHButton.setFont(font);
      cnvMenu.add(beastHButton);

      // cnvMenu.addSeparator();
      menuBar.add(cnvMenu);
    }

    JMenu act = new JMenu("Actions");
    act.setMnemonic(KeyEvent.VK_A);
    menuBar.add(act);

    JMenuItem launchScatter = new JMenuItem();
    launchScatter.setText(TO_SCATTER_PLOT);
    launchScatter.setMnemonic(KeyEvent.VK_S);
    launchScatter.setFont(font);
    launchScatter.addActionListener(this);
    act.add(launchScatter);
    JMenuItem launchComp = new JMenuItem();
    launchComp.setText(TO_COMP_PLOT);
    launchComp.setMnemonic(KeyEvent.VK_C);
    launchComp.setFont(font);
    launchComp.addActionListener(this);
    act.add(launchComp);
    // act.addSeparator();
    {

      previouslyLoadedManagers = new Hashtable<String, ColorExt.ColorManager<String>>();
      JMenu colorMenu = new JMenu("Colors");
      ArrayList<String[]> optsTmp = new ArrayList<String[]>();
      optsTmp.add(new String[] {"Default", "Default", "Set color scheme to the Genvisis default"});
      if (gcModel != null) {
        optsTmp.add(new String[] {"GC content", "GC content", "Color by GC content"});
      }
      if (proj.MARKER_COLOR_KEY_FILENAMES.getValue() != null) {
        for (int i = 0; i < proj.MARKER_COLOR_KEY_FILENAMES.getValue().length; i++) {
          optsTmp.add(new String[] {ext.rootOf(proj.MARKER_COLOR_KEY_FILENAMES.getValue()[i]),
                                    proj.MARKER_COLOR_KEY_FILENAMES.getValue()[i],
                                    "Use this custom color file"});
        }
      }

      for (String otherColor : otherColors) {
        if (Files.exists(otherColor)) {
          optsTmp.add(new String[] {ext.rootOf(otherColor), otherColor,
                                    "Use this custom color file"});
        }
      }

      ItemListener colorListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent arg0) {
          if (arg0.getStateChange() == ItemEvent.SELECTED) {
            String cmd = ((AbstractButton) arg0.getSource()).getActionCommand();
            if ("Default".equals(cmd)) {
              currentColorManager = null;
            } else if ("GC content".equals(cmd)) {
              if (gcModel == null) {
                log.reportTimeError("Internal error, null gc model");
              } else {
                currentColorManager = gcModel.getColorManager();// stored within, doesent regenerate
              }
            } else if (proj.MARKER_COLOR_KEY_FILENAMES.getValue() != null
                       && ext.indexOfStr(cmd, proj.MARKER_COLOR_KEY_FILENAMES.getValue()) >= 0) {
              int index = ext.indexOfStr(cmd, proj.MARKER_COLOR_KEY_FILENAMES.getValue());
              if (previouslyLoadedManagers.containsKey(cmd)) {
                currentColorManager = previouslyLoadedManagers.get(cmd);
              } else {
                ColorManager<String> tmp =
                                         ColorExt.getColorManager(proj,
                                                                  proj.MARKER_COLOR_KEY_FILENAMES.getValue()[index]);
                previouslyLoadedManagers.put(cmd, tmp);
                currentColorManager = previouslyLoadedManagers.get(cmd);
              }
            } else if (ext.indexOfStr(cmd, otherColors) >= 0) {
              int index = ext.indexOfStr(cmd, otherColors);
              if (previouslyLoadedManagers.containsKey(cmd)) {
                currentColorManager = previouslyLoadedManagers.get(cmd);
              } else {
                try {
                  ColorManager<String> tmp = GcModel.populateFromFile(otherColors[index], true, log)
                                                    .getColorManager();
                  previouslyLoadedManagers.put(cmd, tmp);
                  currentColorManager = previouslyLoadedManagers.get(cmd);
                } catch (Exception e) {
                  log.reportTimeWarning("Could not load additional gc-model file "
                                        + otherColors[index]);
                }
              }
            }


            else {
              log.reportTimeError("Internal error, Invalid color command");
            }
            updateQCDisplay();
          }
        }
      };
      ButtonGroup colorBtnGroup = new ButtonGroup();
      for (int i = 0; i < optsTmp.size(); i++) {

        JRadioButtonMenuItem colorButton = new JRadioButtonMenuItem(optsTmp.get(i)[0]);
        colorButton.setActionCommand(optsTmp.get(i)[1]);
        colorButton.setToolTipText(optsTmp.get(i)[2]);
        colorButton.addItemListener(colorListener);
        colorButton.setFont(font);
        colorBtnGroup.add(colorButton);
        colorMenu.add(colorButton);
        if (i == 0) {
          colorButton.setSelected(true);
        }
      }
      menuBar.add(colorMenu);

    }
    {

      JMenu excudeMenu = new JMenu("ExcludeBy");
      ArrayList<String[]> optsTmp = new ArrayList<String[]>();
      optsTmp.add(new String[] {"None", "None", "No Exclusions"});
      if (proj.MARKER_COLOR_KEY_FILENAMES.getValue() != null) {
        for (int i = 0; i < proj.MARKER_COLOR_KEY_FILENAMES.getValue().length; i++) {
          optsTmp.add(new String[] {ext.rootOf(proj.MARKER_COLOR_KEY_FILENAMES.getValue()[i]),
                                    proj.MARKER_COLOR_KEY_FILENAMES.getValue()[i],
                                    "Use this custom color file"});
        }
      }
      ItemListener excludeListener = new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent arg0) {
          if (arg0.getStateChange() == ItemEvent.SELECTED) {
            String cmd = ((AbstractButton) arg0.getSource()).getActionCommand();
            if ("None".equals(cmd)) {
              excludeManager = null;
            } else if (proj.MARKER_COLOR_KEY_FILENAMES.getValue() != null
                       && ext.indexOfStr(cmd, proj.MARKER_COLOR_KEY_FILENAMES.getValue()) >= 0) {
              int index = ext.indexOfStr(cmd, proj.MARKER_COLOR_KEY_FILENAMES.getValue());
              if (previouslyLoadedManagers.containsKey(cmd)) {
                excludeManager = previouslyLoadedManagers.get(cmd);
              } else {
                ColorManager<String> tmp =
                                         ColorExt.getColorManager(proj,
                                                                  proj.MARKER_COLOR_KEY_FILENAMES.getValue()[index]);
                previouslyLoadedManagers.put(cmd, tmp);
                excludeManager = previouslyLoadedManagers.get(cmd);
              }
            }

            else {
              log.reportTimeError("Internal error, Invalid color command");
            }
            updateQCDisplay();
          }
        }
      };
      ButtonGroup excludeBtnGroup = new ButtonGroup();
      for (int i = 0; i < optsTmp.size(); i++) {

        JRadioButtonMenuItem excludeButton = new JRadioButtonMenuItem(optsTmp.get(i)[0]);
        excludeButton.setActionCommand(optsTmp.get(i)[1]);
        excludeButton.setToolTipText(optsTmp.get(i)[2]);
        excludeButton.addItemListener(excludeListener);
        excludeButton.setFont(font);
        excludeBtnGroup.add(excludeButton);
        excudeMenu.add(excludeButton);
        if (i == 0) {
          excludeButton.setSelected(true);
        }
      }
      menuBar.add(excudeMenu);

    }
    return menuBar;
  }

  @Override
  public void actionPerformed(ActionEvent ae) {
    String command = ae.getActionCommand();
    // String[] filenames;

    if (command.equals(FIRST_REGION)) {
      if (regions == null || regions.length == 0) {
        JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error",
                                      JOptionPane.ERROR_MESSAGE);
        return;
      }
      showRegion(0);
    } else if (command.equals(PREVIOUS_REGION)) {
      if (regions == null || regions.length == 0) {
        JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error",
                                      JOptionPane.ERROR_MESSAGE);
        return;
      }
      showRegion(Math.max(regionIndex - 1, 0));
    } else if (command.equals(NEXT_REGION)) {
      if (regions == null || regions.length == 0) {
        JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error",
                                      JOptionPane.ERROR_MESSAGE);
        return;
      }
      showRegion(Math.min(regionIndex + 1, regions.length - 1));
    } else if (command.equals(LAST_REGION)) {
      if (regions == null || regions.length == 0) {
        JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error",
                                      JOptionPane.ERROR_MESSAGE);
        return;
      }
      showRegion(regions.length - 1);
    } else if (command.equals(TO_SCATTER_PLOT)) {
      if (proj == null) {
        JOptionPane.showConfirmDialog(this, "Error - a Project is required to open ScatterPlot",
                                      "Error - no Project", JOptionPane.CANCEL_OPTION,
                                      JOptionPane.ERROR_MESSAGE);
        return;
      }
      String[] listOfMarkers;

      listOfMarkers = new String[stopMarker - startMarker];
      if (listOfMarkers.length == 0) {
        JOptionPane.showMessageDialog(null,
                                      "There are no markers within this region; ScatterPlot will not bother launching",
                                      "Error", JOptionPane.ERROR_MESSAGE);
      } else {
        for (int i = startMarker; i < stopMarker; i++) {
          listOfMarkers[i - startMarker] = markerNames[i];
        }
        ScatterPlot.createAndShowGUI(proj, listOfMarkers, null, false);
      }
    } else if (command.equals(TO_COMP_PLOT)) {
      if (proj == null) {
        JOptionPane.showConfirmDialog(this, "Error - a Project is required to open CompPlot",
                                      "Error - no Project", JOptionPane.CANCEL_OPTION,
                                      JOptionPane.ERROR_MESSAGE);
        return;
      }
      final Region toRegion = new Region(new int[] {chr, start, stop});
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          CompPlot cp = new CompPlot(proj);
          cp.setRegion(toRegion);
        }
      });
    } else {
      System.err.println("Error - unknown command '" + command + "'");
    }
  }

  @Override
  public void mouseClicked(MouseEvent e) {
    Trailer.this.requestFocusInWindow();
    if (e.getButton() == MouseEvent.BUTTON1) {
      if (leftClick != null && leftClick.isAlive()) {
        leftClick.cancel();
        leftClick = null;
        zoomProportionally(false, e.getPoint(), true);
      } else {
        new Thread(leftClick = new SingleClick(this, MouseEvent.BUTTON1,
                                               DOUBLE_CLICK_INTERVAL)).start();
      }

    } else if (e.getButton() == MouseEvent.BUTTON3) {
      if (rightClick != null && rightClick.isAlive()) {
        rightClick.cancel();
        rightClick = null;
        zoom(1, 1);
      } else {
        new Thread(rightClick = new SingleClick(this, MouseEvent.BUTTON3,
                                                DOUBLE_CLICK_INTERVAL)).start();
      }
    }
  }

  @Override
  public void singleLeftClick() {
    // System.out.println("Single left click");
  }

  @Override
  public void singleRightClick() {
    // System.out.println("Single right click");
  }

  @Override
  public void mouseEntered(MouseEvent e) {}

  @Override
  public void mouseExited(MouseEvent e) {}

  @Override
  public void mousePressed(MouseEvent e) {
    startX = e.getPoint().x;
    inDrag = true;
  }

  @Override
  public void mouseReleased(MouseEvent e) {
    inDrag = false;
  }

  @Override
  public void mouseDragged(MouseEvent e) {
    int curX = e.getPoint().x;
    int distance = startX - curX;

    distance *= (stop - start) / (getWidth() - 2 * WIDTH_BUFFER);

    if (distance < 0) {
      distance = Math.max(distance, 1 - start);
    } else {
      distance = Math.min(distance, positions[chrBoundaries[chr][1]] - stop);
    }

    if ((start <= 1 && distance < 0)
        || (stop >= positions[chrBoundaries[chr][1]] && distance > 0)) {

    } else {
      start += distance;
      stop += distance;
    }

    if (inDrag) {
      updateGUI();
      startX = curX;
    }
  }

  @Override
  public void mouseMoved(MouseEvent e) {}

  @Override
  public void mouseWheelMoved(MouseWheelEvent e) {
    zoomProportionally(e.getWheelRotation() > 0, e.getPoint(), false);
  }

  public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
    int width = lrrPanel.getWidth() - 2 * WIDTH_BUFFER;
    double x = p.getX() - WIDTH_BUFFER;
    double xHat = width - x;
    double left = x / width;
    double right = xHat / width;
    double multiplier = MOUSE_WHEEL_MULTIPLIER / (outNotIn ? 1 : -2);

    if (!outNotIn && center) {
      right = 0.25 - right;
      left = 0.25 - left;
      multiplier = MOUSE_WHEEL_MULTIPLIER;
    }

    zoom(left * multiplier, right * multiplier);

  }

  public void zoom(double leftProportion, double rightProportion) {
    int dist = stop - start;
    start = start - (int) (leftProportion * dist);
    stop = stop + (int) (rightProportion * dist);
    updateGUI();
  }

  class MessageOfEncouragment implements Runnable {
    private final String message;
    private final Project proj;
    private boolean noLongerNecessary;

    public MessageOfEncouragment(String message, Project proj) {
      this.message = message;
      this.proj = proj;
      noLongerNecessary = false;
    }

    @Override
    public void run() {
      int count;

      count = 0;
      while (!noLongerNecessary && count < 6) {
        try {
          Thread.sleep(500);
        } catch (InterruptedException ie) {
        }
        count++;
      }

      if (!noLongerNecessary) {
        proj.message(message, "Patience...", JOptionPane.INFORMATION_MESSAGE);
      }
    }

    public void disregard() {
      noLongerNecessary = true;
    }

  }

  public void createSampleList() {
    long time;
    String[] filesPresent;
    FontMetrics fontMetrics;
    String refresh;
    int maxWidth;
    MessageOfEncouragment mess;

    time = new Date().getTime();
    log.report("  Getting a list of all files with extension " + Sample.SAMPLE_FILE_EXTENSION
               + " (if the process hangs here the first time after reverse transposing, please be patient, the operating system is busy indexing the new files) ...");
    mess =
         new MessageOfEncouragment("Getting a list of all sample files is taking longer than usual and probably means that your recently created files are still being indexed on the hard drive. Please be patient...",
                                   proj);
    new Thread(mess).start();
    filesPresent = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true),
                              Sample.SAMPLE_FILE_EXTENSION, jar);
    log.report("Getting list of files took " + ext.getTimeElapsed(time));
    time = new Date().getTime();
    fontMetrics = sampleList.getFontMetrics(sampleList.getFont());
    refresh = "refresh list";
    maxWidth = fontMetrics.stringWidth(refresh);
    System.out.println("  computing font metrics took " + ext.getTimeElapsed(time));
    System.out.println("Determined sample list in " + ext.getTimeElapsed(time));
    mess.disregard();

    if (filesPresent == null || filesPresent.length == 0) {
      // samplesPresent = new String[] {proj.get(Project.SAMPLE_DIRECTORY)+" directory is empty",
      // refresh};
      samplesPresent = new String[] {proj.SAMPLE_DIRECTORY.getValue(false, true)
                                     + " directory is empty", refresh};
      maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(samplesPresent[0]));
    } else {
      samplesPresent = new String[filesPresent.length + 1];
      for (int i = 0; i < filesPresent.length; i++) {
        samplesPresent[i] = filesPresent[i].substring(0, filesPresent[i].lastIndexOf("."));
        maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(samplesPresent[i]));
      }
      samplesPresent[filesPresent.length] = refresh;
    }

    sampleList.setModel(new DefaultComboBoxModel(samplesPresent));
    sampleList.setPreferredSize(new Dimension(maxWidth + 50, 30));

  }

  public boolean loadMarkers() {
    Hashtable<String, String> hash;
    byte[] chrs;
    long time;
    int chr;

    time = new Date().getTime();

    hash = proj.getFilteredHash();
    markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
    if (markerSet == null) {
      JOptionPane.showMessageDialog(null,
                                    "Error - Failed to load the MarkerSet file; make sure the raw data is parsed",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      log.reportError("Error - failed to load MarkerSet for project " + proj.PROJECT_NAME.getValue()
                      + "; make sure the raw data is parsed");
      return false;
    }
    fingerprint = markerSet.getFingerprint();
    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();

    dropped = new boolean[markerNames.length];
    chrBoundaries = new int[27][2];
    for (int i = 0; i < chrBoundaries.length; i++) {
      // chrBoundaries[i][0] = chrBoundaries[i][1] = -1;
      chrBoundaries[i][0] = chrBoundaries[i][1] = 0;
    }
    chr = 0;
    for (int i = 0; i < markerNames.length; i++) {
      dropped[i] = hash.containsKey(markerNames[i]);
      if (chrs[i] > chr) {
        if (chr != 0) {
          chrBoundaries[chr][1] = i - 1;
        }
        chr = chrs[i];
        chrBoundaries[chr][0] = i;
      }
    }
    chrBoundaries[chr][1] = markerNames.length - 1;
    chrBoundaries[0][0] = 0;
    chrBoundaries[0][1] = markerNames.length - 1;

    System.out.println("Read in data for " + markerNames.length + " markers in "
                       + ext.getTimeElapsed(time));
    return true;
  }

  public void loadValues() {
    // Need to load every time?
    if (samp == null || !samp.getSampleName().equals(sample)) {
      samp = proj.getPartialSampleFromRandomAccessFile(sample, false, true, true, true, false);
    }
    if (samp == null) {
      System.err.println("Error - unspecified sample from "
                         + proj.SAMPLE_DIRECTORY.getValue(false, true));
    } else if (samp.getFingerprint() != fingerprint) {
      System.err.println("Error - Sample " + proj.SAMPLE_DIRECTORY.getValue(false, true) + sample
                         + Sample.SAMPLE_FILE_EXTENSION + " has a different fingerprint ("
                         + samp.getFingerprint() + ") than the MarkerSet (" + fingerprint + ")");
    } else {

      if (currentCentroid != null && currentCentroid.startsWith(SEX_CENT)
          && autoSwitch.isSelected()) {
        SampleData sampleData = proj.getSampleData(0, false);
        int sex = sampleData.getSexForIndividual(samp.getSampleName());
        if (sex == 1) {
          if (currentCentroid.endsWith("Female") && !isSettingCentroid) {
            centButtonMap.get(SEX_CENT + " - Male").setSelected(true);
            log.report("Switching to specified male centroid file");
          }
        } else if (sex == 2) {
          if (currentCentroid.endsWith("Male") && !isSettingCentroid) {
            centButtonMap.get(SEX_CENT + " - Female").setSelected(true);
            log.report("Switching to specified female centroid file");
          }
        } else {
          log.report("Warning - no sex specified for sample " + samp.getSampleName()
                     + "; using currently selected centroid file");
        }
      }

      lrrs = samp.getLRRs();
      bafs = samp.getBAFs();
      genotypes = samp.getAB_Genotypes();

      if (transformation_type > 0) {
        lrrValues = Transforms.transform(lrrs, transformation_type, transformSeparatelyByChromosome,
                                         markerSet);
      } else if (transformation_type < 0 && centroids != null) {
        lrrValues = samp.getLRRs(centroids);
        originalBAFs = bafs;
        bafs = samp.getBAFs(centroids);
      } else if (gcParameterDisplay.getCurrentParamIndex() >= 0) {
        int sampleIndex = gcParameterDisplay.getSampleIndex().get(sample);
        GcAdjustorParameters tmp =
                                 gcParameterDisplay.getParams()
                                                   .get(gcParameterDisplay.getCurrentParamIndex());
        lrrValues = samp.getGCCorrectedLRR(tmp, sampleIndex, log);
        bafs = samp.getBAF(tmp, sampleIndex, log);
      } else {
        lrrValues = lrrs;
        originalBAFs = bafs;
      }
    }
    // lrrMin = Math.floor(lrrMin);
    // lrrMax = Math.ceil(lrrMax);
  }


  /**
   * Updates the {@link #chr}, {@link #start} and {@link #stop} fields, based on the provided
   * location.
   *
   * @param location Either a valid gene name recognized by the installed gene track, or a UCSC
   *        location of the format "chr#:start-stop"
   */
  private void parseLocation(String location) {
    byte oldChr;
    int[] loc;

    oldChr = chr;

    if (location == null) {
      System.err.println("Error - null location");
      return;
    }

    if (!location.startsWith("chr")) {
      if (track == null) {
        JOptionPane.showMessageDialog(this,
                                      "Cannot parse '" + location
                                            + "' since the gene track has either not been installed or did not load properly.",
                                      "Error", JOptionPane.ERROR_MESSAGE);
        return;
      } else {
        loc = track.lookupPosition(location);
        if (loc[0] == -1) {
          JOptionPane.showMessageDialog(this,
                                        "'" + location
                                              + "' is not a valid gene name and is not a valid UCSC location.",
                                        "Error", JOptionPane.ERROR_MESSAGE);
          return;
        }
      }
    } else {
      loc = Positions.parseUCSClocation(location);
    }

    if (loc == null) {
      JOptionPane.showMessageDialog(this, "'" + location + "' is not a valid UCSC location.",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    chr = (byte) loc[0];
    if (chr == -1) {
      chr = oldChr;
      return;
    }
    if (chr != oldChr) {
      start = stop = -1;
      updateCNVs(chr);
    }
    start = loc[1];
    stop = loc[2];
    if (start == -1 || start < 0) {
      start = 1;
    }
    if (stop == -1 || stop > positions[chrBoundaries[chr][1]]) {
      stop = positions[chrBoundaries[chr][1]];
    }
    if (start >= stop) {
      start = stop - 1;
    }
  }

  /**
   * Does {@link #parseLocation(String)} and also refreshes the UI.
   *
   * @param location Either a valid gene name recognized by the installed gene track, or a UCSC
   *        location of the format "chr#:start-stop"
   */
  @Override
  public void setPosition(String location) {
    parseLocation(location);
    updateGUI();
  }

  public void updateSample(String newSample) {
    boolean found;
    long time;

    found = sampleList.getSelectedItem().equals(newSample);
    if (found) {
      time = new Date().getTime();
      // System.out.println("Loading CNVs...");
      // System.out.println("took "+ext.getTimeElapsed(time));
      sample = newSample;
      time = new Date().getTime();
      System.out.print("Found " + sample + "...");
      indiPheno = sampleData.getIndiPheno(sample.toLowerCase());

      if (indiPheno == null) {
        // if (!sample.equals(proj.get(Project.SAMPLE_DIRECTORY)+" directory is empty")) {
        if (!sample.equals(proj.SAMPLE_DIRECTORY.getValue(false, true) + " directory is empty")) {
          JOptionPane.showMessageDialog(this,
                                        "Sample '" + sample
                                              + "' was not present in the SampleData file",
                                        "Error", JOptionPane.ERROR_MESSAGE);
        }
        return;
      } else {
        if (proj.CNV_FILENAMES.getValue() == null || proj.CNV_FILENAMES.getValue().length == 0
            || indiPheno.getCnvClasses().size() <= proj.CNV_FILENAMES.getValue().length) {
          for (int i = 0; i < INTERNAL_CNV_TYPES.values().length; i++) {
            indiPheno.getCnvClasses().add(new Hashtable<String, CNVariant[]>());
          }

        }
      }
      loadValues();

      if (REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
        loadCNVsAsRegions();
      }

      updateCNVs(chr);
      updateGUI();

      if (REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
        showRegion(0);
      }
      System.out.println("updated in " + ext.getTimeElapsed(time));
    } else {
      for (int i = 0; i < samplesPresent.length && !found; i++) {
        if (samplesPresent[i].equals(newSample)) {
          sampleList.setSelectedIndex(i);
          found = true;
        }
      }
      if (!found) {
        if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, true) + newSample
                         + Sample.SAMPLE_FILE_EXTENSION, jar)) {
          createSampleList();
          updateSample(newSample);
        } else {
          JOptionPane.showMessageDialog(this,
                                        "Data was not found for the next sample in the list ("
                                              + newSample + ").",
                                        "Error", JOptionPane.ERROR_MESSAGE);
        }
      }
    }

  }

  public void updateGUI() {
    if (start <= 1) {
      startMarker = chrBoundaries[chr][0];
      start = 1;
    } else {
      startMarker = Array.binarySearch(positions, start, chrBoundaries[chr][0],
                                       chrBoundaries[chr][1], false);
    }

    if (stop >= positions[chrBoundaries[chr][1]]) {
      stop = positions[chrBoundaries[chr][1]];
      stopMarker = chrBoundaries[chr][1];
    } else {
      stopMarker = Array.binarySearch(positions, stop, chrBoundaries[chr][0], chrBoundaries[chr][1],
                                      false);
    }

    if (startMarker == -1) {
      System.err.println("Error - failed to find startMarker");
      startMarker = chrBoundaries[chr][0];
    }

    if (stopMarker == -1) {
      System.err.println("Error - failed to find stopMarker");
      stopMarker = chrBoundaries[chr][1];
    }

    displayIndex();

    updateQC(true, true, true);

    repaint();
  }

  private void updateQC(final boolean updateGenome, final boolean updateChr,
                        final boolean updateRegion) {
    if (updateQCThread != null) {
      return;
    }
    updateQCThread = new Thread(new Runnable() {
      @Override
      public void run() {
        if (samp != null) {
          while (inDrag) {
            Thread.yield();
          }
          boolean[] markersForCallrate = null;
          GcModel gcModelToUse = gcModel != null && gcCorrectButton.isSelected() ? gcModel : null;
          // boolean fastQC = false;//gcFastButton.isSelected();
          GC_CORRECTION_METHOD correctionMethod = GC_CORRECTION_METHOD.GENVISIS_GC;// !gcFastButton.isSelected();
          if (updateGenome) {
            boolean[] markersForEverythingElseGenome = Array.booleanNegative(dropped);
            qcGenome = LrrSd.LrrSdPerSample(proj, markerSet, sample, samp, centroids,
                                            markersForCallrate, markersForEverythingElseGenome,
                                            gcModelToUse, correctionMethod, log);
          }
          // if (updateChr) {
          // boolean[] markersForEverythingElseChromosome = Array.booleanNegative(dropped);
          // byte[] chrs = markerSet.getChrs();
          // // TODO check array lengths are the same
          // for (int i = 0; i < chrs.length; i++) {
          // if (chrs[i] != chr) {
          // markersForEverythingElseChromosome[i] = false;
          // }
          // }
          // qcChromo = LrrSd.LrrSdPerSample(proj, sample, samp, centroids, markersForCallrate,
          // markersForEverythingElseChromosome, gcModelToUse, fastQC, log);
          // }
          if (updateRegion) {
            boolean[] markersForEverythingElseRegion = Array.booleanNegative(dropped);
            for (int i = 0; i < markersForEverythingElseRegion.length; i++) {
              if (i < startMarker || i > stopMarker) {
                markersForEverythingElseRegion[i] = false;
              }
            }
            // GC correction not valid for a region
            qcRegion = LrrSd.LrrSdPerSample(proj, markerSet, sample, samp, centroids,
                                            markersForCallrate, markersForEverythingElseRegion,
                                            null, correctionMethod, log);
          }
          updateQCDisplay();
        }
        updateQCThread = null;
      }
    });
    updateQCThread.start();
  }

  private void updateQCDisplay() {
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        String[] qcDetails = null;
        switch (qcSelection) {
          case 0:
            qcDetails = null;
            break;
          case 1:
            qcDetails = qcGenome;
            break;
          // case 2:
          // qcDetails = qcChromo;
          // break;
          case 2:
            qcDetails = qcRegion;
            break;
        }
        // 0 sampleID
        // 1 Array.mean(lrrs, true)
        // 2 Array.stdev(lrrs, true)
        // 3 lrrsdBound
        // 4 Array.stdev(bafs, true)
        // 5 (abCallRate > 0 ? abCallRate : forwardCallRate)
        // 6 (abCallRate > 0 ? abHetRate : forwardHetRate)
        // 7 wfPrior
        // 8 gcwfPrior
        // 9 wfPost
        // 10 gcwfPost
        // 11 lrrsdPost
        // 12 lrrsdPostBound
        // 13 multimodal
        // 14 Array.toStr(bafBinCounts)
        if (qcDetails == null) {
          qcLabel.setText("");
        } else {
          String lrrSd = "LRR SD: " + qcDetails[2];
          String wf = "WF: " + qcDetails[9]; // + qcDetails[7];
          String gcWF = "GC WF: " + qcDetails[10];// + qcDetails[8];
          qcLabel.setText(lrrSd + " | " + wf + " | " + gcWF);
        }
        repaint();
      }
    });
  }

  /**
   * Puts the specified file in the list of region files. Note that {@link #loadRegions()} must be
   * called independently.
   */
  public void loadRegionFile(String filename) {
    addFileToList(filename);
    String file = ext.verifyDirFormat(filename);
    file = file.substring(0, file.length() - 1);
    String name = ext.rootOf(file);
    regionFileNameBtn.get(name).setSelected(true);
    regionFileNameBtn.get(name).doClick();
  }

  public void loadRegions() {
    BufferedReader reader;
    Vector<String[]> v;
    String line;
    String[] parts;
    int ignoredLines, countMissingRegions, invalidSamples;

    try {
      String file = regionFileName.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue()
                                                      + regionFileName
                                                    : regionFileName;
      reader = Files.getAppropriateReader(file);// Files.getReader(file, jar, false, false);
      System.out.print("Loading regions from " + regionFileName + "...");
      v = new Vector<String[]>();
      ignoredLines = countMissingRegions = invalidSamples = 0;
      line = null;
      // FIXME here is the logic for turning strings to regions
      while ((line = reader.readLine()) != null) {
        parts = line.trim().split("\t");
        if (sampleData.lookup(parts[0]) == null) {
          log.reportError("Error - '" + parts[0] + "' is not a valid sample id");
          invalidSamples++;
        } else if (parts.length == 1) {
          v.add(new String[] {parts[0], "chr1"});
          countMissingRegions++;
        } else if (parts.length > 1 && parts[1].startsWith("chr")) {
          v.add(parts);
        } else {
          ignoredLines++;
        }
      }
      System.out.println(" loaded " + v.size() + " regions");
      regions = Matrix.toStringArrays(v);
      if (invalidSamples > 0) {
        JOptionPane.showMessageDialog(null,
                                      "Error - there were " + invalidSamples
                                            + " invalid samples in '" + regionFileName
                                            + "' that were ignored because they could not be found",
                                      "Error", JOptionPane.ERROR_MESSAGE);
      }
      if (countMissingRegions > 0) {
        JOptionPane.showMessageDialog(null,
                                      "Warning - there were " + countMissingRegions + " lines in '"
                                            + regionFileName
                                            + "' without a chromosomal region listed; using \"chr1\" for all missing values",
                                      "Warning", JOptionPane.ERROR_MESSAGE);
      }
      if (ignoredLines > 1) {
        JOptionPane.showMessageDialog(null,
                                      "Error - there were " + ignoredLines + " regions in '"
                                            + regionFileName
                                            + "' that were ignored due to improper formatting",
                                      "Error", JOptionPane.ERROR_MESSAGE);
      }
      reader.close();


    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + regionFileName + "\" not found in data directory");
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + regionFileName + "\"");
    }
  }

  public void loadCNVsAsRegions() {
    int buffer;
    Segment[][] segs;
    int count;

    count = 0;
    segs = new Segment[26][];
    for (int i = 0; i < segs.length; i++) {
      updateCNVs((byte) i);
      segs[i] = findUniqueRegions(cnvs);
      count += segs[i].length;
    }

    regions = new String[count][2];
    count = 0;
    for (Segment[] seg : segs) {
      for (Segment element : seg) {
        regions[count][0] = sample;
        buffer = Math.max(element.getSize() / 2, MIN_BUFFER);
        regions[count][1] = Positions.getUCSCformat(new int[] {element.getChr(),
                                                               element.getStart() - buffer,
                                                               element.getStop() + buffer});
        count++;
      }
    }
  }

  public static Segment[] findUniqueRegions(CNVariant[][] cnvs) { // haven't actually coded the
                                                                  // collapsing of the segments yet
    Segment[] segs;
    int count;

    count = 0;
    for (CNVariant[] cnv : cnvs) {
      count += cnv.length;
    }

    segs = new Segment[count];
    count = 0;
    for (CNVariant[] cnv : cnvs) {
      for (CNVariant element : cnv) {
        segs[count] = element;
        count++;
      }
    }

    return segs;
  }

  /**
   * Updates the active region and sample (if necessary).
   */
  public void showRegion(int regionIndex) {
    Trailer.this.regionIndex = regionIndex;
    if (regions == null || regions.length == 0) {
      regionField.setText("");
      commentLabel.setText(" ");
      return;
    }
    setPosition(regions[regionIndex][1]);
    if (regions[regionIndex].length > 2) {
      commentField.setText(regions[regionIndex][2]);
      commentLabel.setText("region #" + (regionIndex + 1) + ":  " + regions[regionIndex][2]);
    } else {
      commentField.setText("");
      commentLabel.setText(BLANK_COMMENT);
    }

    if (!regions[regionIndex][0].equals(sample)) {
      updateSample(regions[regionIndex][0]);
    }

    regionField.setText((regionIndex + 1) + " of " + regions.length);
  }

  // public static Vector<String[]> loadFileToVec(String filename, boolean demo, boolean
  // ignoreFirstLine, int[] cols, boolean onlyIfAbsent) {
  // BufferedReader reader = null;
  // Vector<String[]> v = new Vector<String[]>();
  // String[] trav;
  // String[] line;
  //
  // try {
  // if (demo) {
  // reader = new BufferedReader(new
  // InputStreamReader(ClassLoader.getSystemResourceAsStream(filename)));
  // } else {
  // reader = new BufferedReader(new FileReader(filename));
  // }
  // if (ignoreFirstLine) {
  // reader.readLine();
  // }
  // while (reader.ready()) {
  // trav = line = reader.readLine().trim().split("\t", -1);
  // if (cols!=null) {
  // trav = new String[cols.length];
  // for (int i = 0; i<cols.length; i++) {
  // trav[i] = line[cols[i]];
  // }
  // }
  // if (onlyIfAbsent) {
  // if (!v.contains(trav)) {
  // v.add(trav);
  // }
  // } else {
  // v.add(trav);
  // }
  // }
  // reader.close();
  // } catch (FileNotFoundException fnfe) {
  // System.err.println("Error: file \""+filename+"\" not found in current directory");
  // System.exit(1);
  // } catch (Exception e) {
  // System.err.println("Error reading file \""+filename+"\"");
  // e.printStackTrace();
  // System.exit(2);
  // }
  //
  // return v;
  // }

  public void displayIndex() {
    if (start < 0) {
      start = 1;
    }
    if (stop > positions[chrBoundaries[chr][1]]) {
      stop = positions[chrBoundaries[chr][1]];
    }

    navPanel.setChrFieldText(chr, start, stop);
  }

  /**
   * Updates the {@link #cnvs} array for the specified chromosome. This first dimension of the
   * resulting array = the number of CNV labels, while the second dimension is the CNV array for the
   * current {@link IndiPheno} at that label index and chromosome.
   */
  private void updateCNVs(byte chr) {
    if (cnvLabels == null) {
      cnvs = new CNVariant[0][];
      return;
    }
    cnvs = new CNVariant[cnvLabels.length][];
    if (indiPheno != null) {
      for (int i = 0; i < cnvLabels.length; i++) {
        cnvs[i] = indiPheno.getCNVs(i, chr);
        if (cnvs[i] == null) {
          cnvs[i] = new CNVariant[0];
        }
      }
    }
  }

  /**
   * @param proj
   * @param lrrsToTransform the log R ratio input array
   * @param transformation_type transformation type for
   *        {@link Transforms#transform(float[], int, boolean, MarkerSet)}
   * @param transformSeparatelyByChromosome transform log R ratios separately
   * @param markerSet
   * @param gcModel a gc model, if the gcmodel is null and correctGC is true, we report an error
   * @param correctGC whether to perform gc correction, must have a valid {@link GcModel} to correct
   * @param correctGCFirst perform the gc correction first, and then any transformations
   * @param log
   * @return
   */
  private static float[] getNewLRRs(Project proj, float[] lrrsToTransform, int transformation_type,
                                    boolean transformSeparatelyByChromosome,
                                    PreparedMarkerSet markerSet, GcModel gcModel, boolean correctGC,
                                    boolean correctGCFirst, Logger log) {
    float[] tmpLrrs = lrrsToTransform; // make sure not to modify
    if (gcModel == null && correctGC) {
      log.reportError("Error - gc Correction was flagged and the model was null, this should not happen...skipping gc correction");
      correctGC = false;
    }
    if (correctGC && correctGCFirst) {
      tmpLrrs = Array.toFloatArray(GcAdjustor
                                             .getComputedAdjustor(proj, markerSet, tmpLrrs, gcModel,
                                                                  GC_CORRECTION_METHOD.GENVISIS_GC,
                                                                  false, false, true)
                                             .getCorrectedIntensities());
    }

    if (transformation_type > 0) {
      tmpLrrs = Transforms.transform(tmpLrrs, transformation_type, transformSeparatelyByChromosome,
                                     markerSet);
    }
    if (correctGC && !correctGCFirst) {
      tmpLrrs = Array.toFloatArray(GcAdjustor
                                             .getComputedAdjustor(proj, markerSet, tmpLrrs, gcModel,
                                                                  GC_CORRECTION_METHOD.GENVISIS_GC,
                                                                  false, false, true)
                                             .getCorrectedIntensities());
    }
    if (Transforms.TRANSFORMATION_TYPES[transformation_type] == Transformations.MAD_SCALED) {
      for (int i = 0; i < tmpLrrs.length; i++) {
        tmpLrrs[i] *= 2;// for now
      }
    }
    return tmpLrrs;
  }

  private void setCentroid(String name) {
    // String name = centroidsSelection.getSelectedItem().toString();
    if (!name.equals(currentCentroid)) {
      String path = namePathMap.get(name);
      Centroids cent = Centroids.load(path, false);
      centroids = cent.getCentroids();
      currentCentroid = name;
    }
  }

  private int prepInternalClasses() {
    int externalCNVs = 0;
    if (proj.CNV_FILENAMES.getValue() != null) {
      for (int i = 0; i < proj.CNV_FILENAMES.getValue().length; i++) {
        if (Files.exists(proj.CNV_FILENAMES.getValue()[i])) {
          externalCNVs++;
        }
      }
    }
    String[] internals = new String[INTERNAL_CNV_TYPES.values().length];
    for (int i = 0; i < INTERNAL_CNV_TYPES.values().length; i++) {
      internals[INTERNAL_CNV_TYPES.values()[i].getIndex()] =
                                                           INTERNAL_CNV_TYPES.values()[i].getName();
    }
    if (sampleData.getCnvClasses().length <= externalCNVs) {
      sampleData.setCnvClasses(Array.concatAll(sampleData.getCnvClasses(), internals));
      cnvLabels = Array.concatAll(sampleData.getCnvClasses());
    }
    if (indiPheno.getCnvClasses().size() <= externalCNVs) {
      for (int i = 0; i < INTERNAL_CNV_TYPES.values().length; i++) {
        indiPheno.getCnvClasses().add(new Hashtable<String, CNVariant[]>());
      }
    }
    return externalCNVs;
  }

  private void beastHere(Segment beastSeg, INTERNAL_CNV_TYPES type) {
    if (type != INTERNAL_CNV_TYPES.BEAST_SCORE && type != INTERNAL_CNV_TYPES.BEAST_SCORE_CUSTOM) {
      throw new IllegalArgumentException("Method only for internal beast types");
    }
    CNVBuilder builder = new CNVBuilder();
    builder.chr(beastSeg.getChr());
    builder.start(beastSeg.getStart());
    builder.stop(beastSeg.getStop());
    builder.cn(2);
    builder.familyID(sample);
    builder.individualID(sample);
    int[] indices = markerSet.getIndicesOfMarkersIn(beastSeg, null, log);
    builder.numMarkers(indices.length);
    BeastScore beastScore =
                          new BeastScore(lrrValues, new int[][] {markerSet.getIndicesByChr()[chr]},
                                         new int[][] {indices}, log);
    beastScore.computeBeastScores();
    builder.score(beastScore.getBeastHeights()[0]);
    int externalCNVs = prepInternalClasses();

    addCnvsToPheno(new CNVariant[] {builder.build()}, externalCNVs, type);
    sampleData.getSampleHash().put(sample.toLowerCase(), indiPheno);
    updateCNVs(chr);
    updateGUI();
  }

  private void quantHere(Segment quantSeg, boolean checkAlreadyCalled) {
    MosaicQuantWorker worker = new MosaicQuantWorker(new Segment[] {quantSeg}, proj, sample,
                                                     MOSAIC_TYPE.values(), 5);
    CNVBuilder builder = new CNVBuilder();
    builder.chr(quantSeg.getChr());
    builder.start(quantSeg.getStart());
    builder.stop(quantSeg.getStop());
    builder.cn(2);
    builder.familyID(sample);
    builder.individualID(sample);
    builder.numMarkers(stopMarker - startMarker);
    builder.score(Double.NaN);
    CNVariant[] tmp = new CNVariant[MOSAIC_TYPE.values().length];
    WorkerHive<MosaicQuantResults[]> hive =
                                          new WorkerHive<MosaicismQuant.MosaicQuantResults[]>(1, 10,
                                                                                              proj.getLog());
    hive.addCallable(worker);
    hive.execute(true);
    for (int i = 0; i < tmp.length; i++) {
      tmp[i] = builder.build();
    }
    ArrayList<MosaicQuantResults[]> mqrs = hive.getResults();
    MosaicQuantResults[] mqr = mqrs.get(0);
    for (int i = 0; i < mqr.length; i++) {
      builder.score(mqr[i].getFs()[0]);
      builder.numMarkers(mqr[i].getNumMarkers()[0]);
      tmp[i] = builder.build();
    }
    int externalCNVs = prepInternalClasses();
    // private static final String[] INTERNAL_CNV_CLASSES = new String[] { "CNVCaller",
    // "RevCNVCaller", "Consensus", "MosaicCaller", "MONOSOMY_DISOMYF", "CUSTOMF", "BEAST_SCORE",
    // "BEAST_SCORE_CUSTOM" };
    // private static final int[] INTERNAL_CNV_CLASSES_INDICES = new int[] { 0, 1, 2, 3, 4, 5, 6, 7
    // };
    addCnvsToPheno(tmp, externalCNVs, INTERNAL_CNV_TYPES.MONSOMOMY_DYSOMYF);
    if (!checkAlreadyCalled || selectedCNV == null
        || selectedCNV[0] != externalCNVs + INTERNAL_CNV_TYPES.MOSAIC_CALLER.getIndex()) {
      MosaicBuilder builderMosaic = new MosaicBuilder();
      builderMosaic.verbose(true);
      MosaicismDetect md = builderMosaic.build(proj, sample, markerSet, Array.toDoubleArray(bafs));
      LocusSet<MosaicRegion> mosSet = md.callMosaic(quantSeg, true);

      if (mosSet.getLoci().length != 1) {
        proj.getLog().reportTimeError("Mosaic caller not in force call mode");
      } else {
        addCnvsToPheno(new CNVariant[] {mosSet.getLoci()[0]}, externalCNVs,
                       INTERNAL_CNV_TYPES.CUSTOMF);
      }
    }
    sampleData.getSampleHash().put(sample.toLowerCase(), indiPheno);
    updateCNVs(chr);
    updateGUI();
  }



  /**
   * Remove a cnv, currently operates on segments
   */
  private void removeCnvFromPheno(Segment remove, int externalCNVs, boolean all,
                                  INTERNAL_CNV_TYPES type) {
    int key = externalCNVs + type.getIndex();
    CNVariant[] tmpCurrent = indiPheno.getCnvClasses().get(key).get(chr + "");
    if (!all) {
      if (tmpCurrent != null && tmpCurrent.length > 0) {
        ArrayList<CNVariant> retain = new ArrayList<CNVariant>();
        for (int i = 0; i < tmpCurrent.length; i++) {
          if (!remove.equals(tmpCurrent[i])) {
            retain.add(tmpCurrent[i]);
          }
        }
        tmpCurrent = retain.toArray(new CNVariant[retain.size()]);
      }
    } else {
      tmpCurrent = new CNVariant[] {};
    }
    indiPheno.getCnvClasses().get(key).put(chr + "", tmpCurrent);
    updateCNVs(chr);
    updateGUI();

  }

  /**
   * @param tmp
   * @param externalCNVs the number of external (file based) cnv classes
   * @param internalIndex the index of the internal index
   */
  private void addCnvsToPheno(CNVariant[] tmp, int externalCNVs, INTERNAL_CNV_TYPES type) {
    int key = externalCNVs + type.getIndex();
    CNVariant[] tmpCurrent = indiPheno.getCnvClasses().get(key).get(chr + "");

    if (tmpCurrent != null && tmpCurrent.length > 0) {
      boolean[] use = Array.booleanArray(tmp.length, true);
      for (CNVariant element : tmpCurrent) {
        for (int j = 0; j < tmp.length; j++) {
          if (use[j]) {
            if (element.equals(tmp[j])) {
              use[j] = false;
            }
          }
        }
      }
      ArrayList<CNVariant> uniqAdd = new ArrayList<CNVariant>();
      for (int i = 0; i < use.length; i++) {
        if (use[i]) {
          uniqAdd.add(tmp[i]);
        }
      }
      CNVariant[] uniq = uniqAdd.toArray(new CNVariant[uniqAdd.size()]);
      tmpCurrent = Array.concatAll(uniq, tmpCurrent);
    } else {
      tmpCurrent = tmp;
    }

    indiPheno.getCnvClasses().get(key).put(chr + "", tmpCurrent);
  }

  /**
   * Updates the region file with the current region state. If no valid region file is currently
   * being tracked, prompts the user if {@code doPrompt} is set. The return value can be assigned to
   * the input variable to avoid future prompts.
   *
   * @param doPrompt Whether or not to prompt the user to create a region file if it doesn't already
   *        exist.
   * @return Updated value of the {@code doPrompt} param
   */
  private boolean promptAndSaveRegions(boolean doPrompt) {
    if (doPrompt && (regionFileName == null || !new File(regionFileName).exists())) {
      int createRegion = JOptionPane.showConfirmDialog(null,
                                                       "No region file found - changes will not persist.\nWould you like to create a region file now?",
                                                       "Warning - temporary changes",
                                                       JOptionPane.YES_NO_OPTION);
      if (JOptionPane.YES_OPTION == createRegion) {
        saveRegionFile();
      } else {
        doPrompt = false;
      }
    }

    if (regionFileName != null && new File(regionFileName).exists()) {
      Files.writeMatrix(regions, regionFileName, "\t");
    }

    return doPrompt;
  }

  // public static CNVariant[] loadCNVfiles(Project proj, String[] filenames) {
  // BufferedReader reader;
  // Vector<CNVariant> v = null;
  // String[] line;
  // String id = null;
  // boolean jar;
  //
  // jar = proj.getJarStatus();
  //
  // v = new Vector<CNVariant>();
  // for (int i = 0; i<filenames.length; i++) {
  // try {
  // reader = Files.getReader(filenames[i], jar, true, false);
  // if (reader!=null) {
  // reader.mark(1000);
  // line = reader.readLine().trim().split("[\\s]+");
  // if (!line[2].toLowerCase().equals("chr")&&ext.chromosomeNumber(line[2])!=-1) {
  // reader.reset();
  // }
  // while (reader.ready()) {
  // line = reader.readLine().trim().split("[\\s]+");
  // if (line[1].equals(id)) {
  // v.add(new CNVariant(line, i));
  // }
  // }
  // reader.close();
  // }
  // } catch (IOException ioe) {
  // System.err.println("Error reading file \""+filenames[i]+"\"");
  // ioe.printStackTrace();
  // }
  //
  // }
  //
  // return CNVariant.sortCNVs(CNVariant.toArray(v));
  // }
  /**
   * Handle the actions for gc parameter files ~ fast version of gc-adjustment
   *
   */
  private static class GCParameterDisplay extends FileActionMenu {
    private final ArrayList<GcAdjustorParameters> params;
    private final Trailer trailer;
    private ItemListener itemListener;
    private final Hashtable<String, Integer> sampleIndex;
    private int currentParamIndex;
    private final Logger log;

    private GCParameterDisplay(Project proj, Trailer trailer, Logger log) {
      super("Adjust with gc parameter file", proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue());

      String[] samps = proj.getSamples();
      sampleIndex = new Hashtable<String, Integer>();
      for (int i = 0; i < samps.length; i++) {
        sampleIndex.put(samps[i], i);
      }
      this.trailer = trailer;
      currentParamIndex = -1;
      params = new ArrayList<GcAdjustorParameters>();
      this.log = log;
      developItemListener();
      developMenu();
      getActionMenu().setToolTipText("This uses a precomputed GC correction to adjust the sample by GC content");
      getActionMenu().setMnemonic(KeyEvent.VK_A);
    }

    public Hashtable<String, Integer> getSampleIndex() {
      return sampleIndex;
    }

    public int getCurrentParamIndex() {
      return currentParamIndex;
    }

    public ArrayList<GcAdjustorParameters> getParams() {
      return params;
    }

    @Override
    public ItemListener getListener() {
      return itemListener;
    }

    private void developItemListener() {
      int num = getNamePathMap().keySet().size();
      for (int i = 0; i < num; i++) {
        params.add(null);
      }

      itemListener = new ItemListener() {

        @Override
        public void itemStateChanged(ItemEvent ie) {
          JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
          if (jrb.isSelected() && !"None".equals(jrb.getText())) {
            if (getNamePathMap() == null || getNamePathMap().isEmpty()) {
              jrb.setSelected(false);
              currentParamIndex = -1;
              return;
            }
            if (getNamePathMap().containsKey(jrb.getText())) {
              currentParamIndex = getNamePathMap().get(jrb.getText());
              if (params.get(currentParamIndex) == null) {// lazy load the params
                String fileToLoad = getExistingFiles().get(currentParamIndex);
                log.reportTimeInfo("Loading GC parameter file " + fileToLoad);
                params.set(currentParamIndex, GcAdjustorParameters.readSerial(fileToLoad, log));
              }
              trailer.gcCorrectButton.setSelected(false);
            } else {
              currentParamIndex = -1;
            }
            trailer.loadValues();
            trailer.updateGUI();
            trailer.repaint();
          } else {
            currentParamIndex = -1;
            trailer.loadValues();
            trailer.updateGUI();
            trailer.repaint();
          }
        }

      };
    }
  }

  /**
   * Organizes dynamic calling of cnvs and similar metrics
   *
   */
  private enum INTERNAL_CNV_TYPES {
                                   CNV_CALLER("CNVCaller", 0),
                                   // REV_CNV_CALLER("RevCNVCaller", 1),
                                   // CONSENSUS("Consensus", 2),
                                   MOSAIC_CALLER("MosaicCaller",
                                                 1), MONSOMOMY_DYSOMYF("MONOSOMY_DISOMYF",
                                                                       2), CUSTOMF("CUSTOMF",
                                                                                   3), BEAST_SCORE("BEAST_Score",
                                                                                                   4), BEAST_SCORE_CUSTOM("BEAST_Score_Custom",
                                                                                                                          5);

    private String name;
    private int index;

    private INTERNAL_CNV_TYPES(String name, int index) {
      this.name = name;
      this.index = index;
    }

    public String getName() {
      return name;
    }

    public int getIndex() {
      return index;
    }

    private static INTERNAL_CNV_TYPES getAppropriate(int selectedCnvs, int externalIndex) {
      if (selectedCnvs < externalIndex) {
        return null;
      } else {
        int key = selectedCnvs - externalIndex;
        for (int i = 0; i < INTERNAL_CNV_TYPES.values().length; i++) {
          if (key == INTERNAL_CNV_TYPES.values()[i].getIndex()) {
            return INTERNAL_CNV_TYPES.values()[i];
          }
        }
        return null;
      }
    }
  }

  public static void main(String[] args) {
    Project proj = new Project(args[0], false);

    // Project proj = new Project(cnv.Launch.getDefaultDebugProjectFile(true), false);
    // Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
    new Trailer(proj, DEFAULT_SAMPLE, proj.CNV_FILENAMES.getValue(), DEFAULT_LOCATION);
  }
}


