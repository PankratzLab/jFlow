package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Set;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListModel;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
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
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.border.SoftBevelBorder;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.gui.JCheckBoxListCellRenderer;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.Transforms.TRANSFORMATION;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.miginfocom.swing.MigLayout;

public class RegionReviewer extends JFrame {

  private static final long serialVersionUID = 1L;

  enum SORT {
    SAMPLE("Sample Order"), LRRMEDIAN("LRR Median");

    private SORT(String d) {
      this.desc = d;
    }

    private final String desc;

    public String getDescription() {
      return desc;
    }
  }

  private Project proj;
  private MarkerDetailSet mds;
  private NavigableMap<Byte, NavigableSet<Marker>> markerChrMap;

  private Segment seg;
  private TRANSFORMATION transform;
  private SORT sorting = SORT.SAMPLE;

  private Dimension panelSize = new Dimension(100, 50);
  private int size = 4;

  private Map<String, JPanel> indivPanels;
  private Map<String, JSeparator> indivSeparators;
  private Map<String, LRRPanel> lrrPanels;
  private Map<String, BAFPanel> bafPanels;
  private Map<String, Double> lrrMedians;

  private List<String> samples;
  private List<String> hidden;
  private String selectedSample;
  private Set<Annotation> annotations;
  private Multimap<String, Annotation> annotationMap;
  private HashMap<Character, Action> mnemonicActions;
  private JMenu mnTravAnn;
  private JMenu mnSaveAnn;
  private JMenu mnDelAnn;
  private JRadioButtonMenuItem jrbTravAll;
  private JRadioButtonMenuItem jrbTravAnn;
  private JRadioButtonMenuItem jrbTravNon;
  private ButtonGroup travGrp;
  private HashMap<Annotation, JRadioButtonMenuItem> annTravMap;
  private HashMap<Annotation, JMenuItem> annDelMap;
  private HashMap<Annotation, JMenuItem> annSaveMap;

  private volatile boolean annotationsChanged = false;
  private volatile boolean checkBxAdvanceAfterAnnotationKeyed = true;

  public static final int WIDTH_BUFFER = 0;
  public static final int HEIGHT_BUFFER = 2;

  private JPanel samplePanel;
  private JScrollPane scrollPane;
  private JLabel loadingLabel;
  private JTextField fldAddAnnot;

  private JPanel annotPanel;
  private JTextField fldRegion;
  private JComboBox<String> cmboSample;
  private AbstractAction selectPrev = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      if (selectedSample == null) return;
      int ind = samples.indexOf(selectedSample);
      if (jrbTravAll.isSelected()) {
        ind--;
        if (ind == -1) return;
      } else if (jrbTravAnn.isSelected()) {
        do {
          ind--;
          if (ind == -1) return;
        } while (annotationMap.get(samples.get(ind)).size() > 0);
      } else if (jrbTravNon.isSelected()) {
        do {
          ind--;
          if (ind == -1) return;
        } while (annotationMap.get(samples.get(ind)).size() == 0);
      } else {
        Annotation sel = null;
        for (Entry<Annotation, JRadioButtonMenuItem> entry : annTravMap.entrySet()) {
          if (entry.getValue().isSelected()) {
            sel = entry.getKey();
            break;
          }
        }
        do {
          ind--;
          if (ind == -1) return;
        } while (!annotationMap.get(samples.get(ind)).contains(sel));
      }
      selectSample(samples.get(ind));
    }
  };

  private AbstractAction selectNext = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      if (selectedSample == null) return;
      int ind = samples.indexOf(selectedSample);
      if (jrbTravAll.isSelected()) {
        ind++;
        if (ind == samples.size()) return; // no
                                           // sample
                                           // found
      } else if (jrbTravAnn.isSelected()) {
        do {
          ind++;
          if (ind == samples.size()) return; // no sample found
        } while (annotationMap.get(samples.get(ind)).size() > 0);
      } else if (jrbTravNon.isSelected()) {
        do {
          ind++;
          if (ind == samples.size()) return; // no
                                             // sample
                                             // found
        } while (annotationMap.get(samples.get(ind)).size() == 0);
      } else {
        Annotation sel = null;
        for (Entry<Annotation, JRadioButtonMenuItem> entry : annTravMap.entrySet()) {
          if (entry.getValue().isSelected()) {
            sel = entry.getKey();
            break;
          }
        }
        do {
          ind++;
          if (ind == samples.size()) return; // no
                                             // sample
                                             // found
        } while (!annotationMap.get(samples.get(ind)).contains(sel));
      }
      selectSample(samples.get(ind));
    }
  };

  private Comparator<String> sampleComparator = new Comparator<String>() {
    @Override
    public int compare(String o1, String o2) {
      return Integer.compare(proj.getSampleIndices().get(o1), proj.getSampleIndices().get(o2));
    }
  };

  private Comparator<String> lrrComparator = new Comparator<String>() {
    @Override
    public int compare(String o1, String o2) {
      return Double.compare(lrrMedians.get(o2) == null ? 0 : lrrMedians.get(o2),
                            lrrMedians.get(o1) == null ? 0 : lrrMedians.get(o1));
    }
  };

  private KeyAdapter mnemonicKeyAdapter = new KeyAdapter() {

    @Override
    public void keyPressed(KeyEvent e) {
      if (e.getKeyCode() != KeyEvent.VK_UP && e.getKeyCode() != KeyEvent.VK_DOWN
          && e.getKeyCode() != KeyEvent.VK_SPACE) {
        e.consume();
        if (!(e.isAltDown() && e.getKeyCode() == KeyEvent.VK_F4)) {
          fireMnem(e);
        } else {
          close();
        }
      }
    }
  };

  private MouseAdapter createSelectMouseAdapter(String sample, JPanel indivPanel) {
    return new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e)) {
          selectSample(sample);
        } else if (SwingUtilities.isRightMouseButton(e)) {
          JPopupMenu menu = new JPopupMenu();
          JMenuItem copySample = new JMenuItem();
          copySample.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e) {
              ext.setClipboard(sample);
            }
          });
          copySample.setText("Copy Sample Name to Clipboard");
          menu.add(copySample);

          JMenuItem hide = new JMenuItem();
          hide.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e) {
              hideSample(sample);
            }
          });
          hide.setText("Hide Sample");
          menu.add(hide);

          JMenuItem trailer = new JMenuItem();
          trailer.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e) {
              new Thread(new Runnable() {
                @Override
                public void run() {
                  Trailer t = new Trailer(proj, sample, proj.CNV_FILENAMES.getValue(),
                                          seg.getUCSClocation());
                  t.setVisible(true, false);
                }
              }).start();
            }
          });
          trailer.setText("Open Sample/Region in Trailer");
          menu.add(trailer);

          menu.show(indivPanel, e.getX(), e.getY());
        }
        indivPanel.requestFocus();
      }
    };
  }

  private AbstractAction actionSetRegion = new AbstractAction() {

    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      int[] pos = Positions.parseUCSClocation(fldRegion.getText());
      if (pos == null) {
        fldRegion.setText(Positions.getUCSCformat(new int[] {seg.getChr(), seg.getStart(),
                                                             seg.getStart()}));
        samplePanel.requestFocus();
        return;
      }
      if (pos[1] == -1) {
        pos[1] = 0;
      }
      if (pos[2] == -1) {
        switch (proj.GENOME_BUILD_VERSION.getValue()) {
          case HG18:
            pos[2] = Positions.CHROMOSOME_LENGTHS_B36_HG18[pos[0]];
            break;
          case HG19:
            pos[2] = Positions.CHROMOSOME_LENGTHS_B37_HG19[pos[0]];
            break;
          case HG38:
            pos[2] = Positions.CHROMOSOME_LENGTHS_B38_HG38[pos[0]];
            break;
        }
      }
      new Thread(new Runnable() {
        @Override
        public void run() {
          setRegion(new Segment((byte) pos[0], pos[1], pos[2]));
          samplePanel.requestFocus();
        }
      }).start();
    }
  };

  private AbstractAction actionJumpSample = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      if (cmboSample.getSelectedIndex() == 0) return;
      selectSample((String) cmboSample.getSelectedItem());
      cmboSample.setSelectedIndex(0);
      samplePanel.requestFocus();
    }
  };

  private AbstractAction actionCreateAnnotation = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent arg0) {
      createNewAnnotation();
      samplePanel.requestFocus();
    }
  };

  private AbstractAction actionShowHideSamples = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      showHideSamples();
    }
  };

  private AbstractAction actionPanelHeight = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      int hgt = -1;
      do {
        Object val = JOptionPane.showInputDialog(RegionReviewer.this, "Input Panel Height",
                                                 "Panel Height", JOptionPane.PLAIN_MESSAGE, null,
                                                 null, (int) panelSize.getHeight());
        if (val == null) break;
        try {
          hgt = Integer.parseInt(val.toString());
        } catch (NumberFormatException e1) {
          hgt = -1;
        }
      } while (hgt == -1);
      if (hgt == -1) return;
      panelSize = new Dimension(100, hgt);
      for (String s : samples) {
        lrrPanels.get(s).setMinimumSize(panelSize);
        bafPanels.get(s).setMinimumSize(panelSize);
        lrrPanels.get(s).revalidate();
        bafPanels.get(s).revalidate();
      }
      scrollPane.getViewport().validate();
      scrollPane.getViewport().repaint();
      samplePanel.requestFocus();
    }
  };

  private AbstractAction actionLoad = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      JFileChooser jfc = new JFileChooser(proj.PROJECT_DIRECTORY.getValue());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Load Annotation File");
      int opt = jfc.showOpenDialog(RegionReviewer.this);
      if (opt == JFileChooser.APPROVE_OPTION) {
        String annFile = jfc.getSelectedFile().getAbsolutePath();
        try {
          load(annFile);
        } catch (IOException e1) {
          proj.getLog().reportException(e1);
        }
      }
    }
  };

  private AbstractAction actionSaveAll = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      JFileChooser jfc = new JFileChooser(proj.PROJECT_DIRECTORY.getValue());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save All Samples / Annotations to File");
      int opt = jfc.showSaveDialog(RegionReviewer.this);
      if (opt == JFileChooser.APPROVE_OPTION) {
        String annFile = jfc.getSelectedFile().getAbsolutePath();
        saveAll(annFile);
      }
    }
  };

  private AbstractAction actionSaveAllAnnotated = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      JFileChooser jfc = new JFileChooser(proj.PROJECT_DIRECTORY.getValue());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save Annotated Samples to File");
      int opt = jfc.showSaveDialog(RegionReviewer.this);
      if (opt == JFileChooser.APPROVE_OPTION) {
        String annFile = jfc.getSelectedFile().getAbsolutePath();
        saveAllAnnotated(annFile);
      }
    }
  };

  private AbstractAction actionSaveAnnotated = new AbstractAction() {
    private static final long serialVersionUID = 1L;

    @Override
    public void actionPerformed(ActionEvent e) {
      String annot = e.getActionCommand();
      Annotation ann = null;
      for (Annotation a : annotations) {
        if (annot.equals(a.annotation)) {
          ann = a;
          break;
        }
      }
      if (ann == null) {
        return;
      }

      JFileChooser jfc = new JFileChooser(proj.PROJECT_DIRECTORY.getValue());
      jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
      jfc.setDialogTitle("Save \"" + ann.annotation + "\" Annotation to File");
      int opt = jfc.showSaveDialog(RegionReviewer.this);
      if (opt == JFileChooser.APPROVE_OPTION) {
        String annFile = jfc.getSelectedFile().getAbsolutePath();
        saveAnnotated(annFile, ann);
      }
    }
  };

  public RegionReviewer(Project proj) {
    setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
    setTitle("Region Review");

    createJMenuBar();

    this.proj = proj;
    this.mds = proj.getMarkerSet();
    this.markerChrMap = mds.getChrMap();
    this.samples = new ArrayList<String>(Arrays.asList(proj.getSamples())) {
      private static final long serialVersionUID = 1L;

      @Override
      public boolean add(String e) {
        boolean b = super.add(e);
        switch (sorting) {
          case SAMPLE:
            Collections.sort(samples, sampleComparator);
            break;
          case LRRMEDIAN:
            Collections.sort(samples, lrrComparator);
        }
        return b;
      }
    };
    this.hidden = new ArrayList<>();
    this.annotations = new LinkedHashSet<>();
    this.annotationMap = HashMultimap.create();
    this.mnemonicActions = new HashMap<>();
    this.annTravMap = new HashMap<>();
    this.annDelMap = new HashMap<>();
    this.annSaveMap = new HashMap<>();

    this.transform = TRANSFORMATION.RAW;

    indivPanels = new HashMap<>();
    indivSeparators = new HashMap<>();
    lrrPanels = new HashMap<>();
    lrrMedians = new HashMap<>();
    bafPanels = new HashMap<>();

    JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
    splitPane.setOneTouchExpandable(false);
    splitPane.setResizeWeight(0.85);
    this.setContentPane(splitPane);

    samplePanel = new JPanel(new MigLayout("flowy, fill, gap 1, hidemode 3", "[]", ""));
    samplePanel.addKeyListener(mnemonicKeyAdapter);
    scrollPane = new JScrollPane(samplePanel);
    scrollPane.getVerticalScrollBar().setUnitIncrement(16);
    splitPane.setLeftComponent(scrollPane);
    String[] samples = proj.getSamples();
    for (int s = 0; s < samples.length; s++) {
      String sample = samples[s];
      lrrPanels.put(sample, new LRRPanel());
      bafPanels.put(sample, new BAFPanel());
      lrrPanels.get(sample).setMinimumSize(panelSize);
      bafPanels.get(sample).setMinimumSize(panelSize);
      JPanel indivPanel = new JPanel(new MigLayout("flowy, fill, gap 1", "[]", ""));
      indivPanel.setBorder(null);
      indivPanel.addMouseListener(createSelectMouseAdapter(sample, indivPanel));
      indivPanel.addKeyListener(mnemonicKeyAdapter);
      JLabel sampleLabel = new JLabel(sample);
      sampleLabel.setFont(new Font("Arial", Font.BOLD, 9));
      indivPanel.add(sampleLabel);
      indivPanel.add(lrrPanels.get(sample), "grow");
      indivPanel.add(bafPanels.get(sample), "grow");
      setMappings(indivPanel);
      setMappings(lrrPanels.get(sample));
      setMappings(bafPanels.get(sample));
      indivPanels.put(sample, indivPanel);
      samplePanel.add(indivPanel, "grow");
      indivSeparators.put(sample, new JSeparator(SwingConstants.HORIZONTAL));
      samplePanel.add(indivSeparators.get(sample), "growx");
    }

    JPanel controlPanel = new JPanel(new MigLayout("hidemode 3", "[grow, fill][]",
                                                   "[][][][][][][][grow]"));
    splitPane.setRightComponent(controlPanel);

    Font labelFont = new Font("Arial", 0, 13);

    int rowIndex = 0;
    JLabel lblRegion = new JLabel("Region:");
    lblRegion.setFont(labelFont);
    controlPanel.add(lblRegion, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;
    fldRegion = new JTextField(25);
    fldRegion.addActionListener(actionSetRegion);
    controlPanel.add(fldRegion, "cell 0 " + rowIndex + " 2 1");
    loadingLabel = new JLabel(Grafik.getImageIcon("images/load.gif"));
    loadingLabel.setVisible(false);
    controlPanel.add(loadingLabel, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;

    JLabel lblSkip = new JLabel("Skip to Sample:");
    lblSkip.setFont(labelFont);
    controlPanel.add(lblSkip, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;
    cmboSample = new JComboBox<>(ArrayUtils.addStrToArray("", proj.getSamples(), 0));
    cmboSample.addActionListener(actionJumpSample);
    controlPanel.add(cmboSample, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;

    JLabel lblAddAnnot = new JLabel("Add Annotation:");
    lblAddAnnot.setFont(labelFont);
    controlPanel.add(lblAddAnnot, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;
    fldAddAnnot = new JTextField(25);
    InputMap im = fldAddAnnot.getInputMap();
    ActionMap am = fldAddAnnot.getActionMap();
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0), "enter");
    am.put("enter", actionCreateAnnotation);
    controlPanel.add(fldAddAnnot, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;

    JLabel lblAnnot = new JLabel("Annotations:");
    lblAnnot.setFont(labelFont);
    controlPanel.add(lblAnnot, "cell 0 " + rowIndex + " 2 1");
    rowIndex++;

    annotPanel = new JPanel();
    JScrollPane scrollPane2 = new JScrollPane(annotPanel);
    scrollPane2.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
    annotPanel.setLayout(new MigLayout("ins 0", "[grow]", "[grow]"));
    controlPanel.add(scrollPane2, "cell 0 " + rowIndex + " 2 1, grow");
    rowIndex++;

    setMappings((JComponent) this.getContentPane());
    setMappings(fldRegion);
    setMappings(fldAddAnnot);
    setMappings(splitPane);
    setMappings(scrollPane);
    setMappings(scrollPane.getViewport());
    setMappings(samplePanel);
    setMappings(controlPanel);

    this.setSize(600, 800);
    setVisible(true);
    repaint();
  }

  private void createJMenuBar() {
    JMenuBar menubar = new JMenuBar();

    JMenu menuOptions = new JMenu();
    menuOptions.setText("Options");
    menuOptions.setMnemonic('O');
    menubar.add(menuOptions);

    JMenuItem itemLoad = new JMenuItem();
    itemLoad.setAction(actionLoad);
    itemLoad.setText("Load File");
    menuOptions.add(itemLoad);

    JMenuItem itemSaveAll = new JMenuItem();
    itemSaveAll.setAction(actionSaveAll);
    itemSaveAll.setText("Save All");
    menuOptions.add(itemSaveAll);

    JMenuItem itemSaveAnn = new JMenuItem();
    itemSaveAnn.setAction(actionSaveAllAnnotated);
    itemSaveAnn.setText("Save All Annotated");
    menuOptions.add(itemSaveAnn);

    mnSaveAnn = new JMenu("Save Annotated:");
    menuOptions.add(mnSaveAnn);

    menuOptions.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem itemShowHide = new JMenuItem();
    itemShowHide.setAction(actionShowHideSamples);
    itemShowHide.setText("Show / Hide Samples");
    itemShowHide.setMnemonic('H');
    menuOptions.add(itemShowHide);

    JMenuItem itemPanelHgt = new JMenuItem();
    itemPanelHgt.setAction(actionPanelHeight);
    itemPanelHgt.setText("Change Panel Height");
    itemPanelHgt.setMnemonic('C');
    menuOptions.add(itemPanelHgt);

    menuOptions.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmNavTravLbl = new JMenuItem();
    mntmNavTravLbl.setEnabled(false);
    mntmNavTravLbl.setText("Traversal:");
    menuOptions.add(mntmNavTravLbl);

    travGrp = new ButtonGroup();

    jrbTravAll = new JRadioButtonMenuItem();
    jrbTravAll.setText("All");
    jrbTravAll.setSelected(true);
    travGrp.add(jrbTravAll);
    menuOptions.add(jrbTravAll);
    jrbTravAnn = new JRadioButtonMenuItem();
    jrbTravAnn.setText("Annotated");
    travGrp.add(jrbTravAnn);
    menuOptions.add(jrbTravAnn);
    jrbTravNon = new JRadioButtonMenuItem();
    jrbTravNon.setText("Non-Annotated");
    travGrp.add(jrbTravNon);
    menuOptions.add(jrbTravNon);
    mnTravAnn = new JMenu("Annotation:");
    menuOptions.add(mnTravAnn);

    JCheckBoxMenuItem mntmNavAdvAfterAnnot = new JCheckBoxMenuItem();
    mntmNavAdvAfterAnnot.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        checkBxAdvanceAfterAnnotationKeyed = mntmNavAdvAfterAnnot.isSelected();
      }
    });
    mntmNavAdvAfterAnnot.setText("Advance after Annotating");
    mntmNavAdvAfterAnnot.setSelected(checkBxAdvanceAfterAnnotationKeyed);
    menuOptions.add(mntmNavAdvAfterAnnot);

    menuOptions.add(new JSeparator(SwingConstants.HORIZONTAL));

    JMenuItem mntmSortLbl = new JMenuItem();
    mntmSortLbl.setEnabled(false);
    mntmSortLbl.setText("Sort By:");
    menuOptions.add(mntmSortLbl);

    ButtonGroup sortGrp = new ButtonGroup();
    for (SORT sort : SORT.values()) {
      JRadioButtonMenuItem jmrb = new JRadioButtonMenuItem();
      jmrb.setAction(new AbstractAction() {
        private static final long serialVersionUID = 1L;

        @Override
        public void actionPerformed(ActionEvent e) {
          RegionReviewer.this.sorting = sort;
          reorder();
        }
      });
      jmrb.setText(sort.getDescription());
      if (this.sorting == sort) {
        jmrb.setSelected(true);
      }
      sortGrp.add(jmrb);
      menuOptions.add(jmrb);
    }

    menuOptions.add(new JSeparator(SwingConstants.HORIZONTAL));

    mnDelAnn = new JMenu("Delete Annotation:");
    menuOptions.add(mnDelAnn);

    this.setJMenuBar(menubar);
  }

  private void close() {
    if (checkAnnotationsChanged()) {
      // TODO
    }
    setVisible(false);
    dispose();
  }

  private void setMappings(JComponent comp) {
    InputMap im = comp.getInputMap(JComponent.WHEN_FOCUSED);
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, 0, false), "UP");
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, 0, false), "DOWN");
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, KeyEvent.SHIFT_DOWN_MASK, false), "UP");
    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0, false), "DOWN");
    ActionMap am = comp.getActionMap();
    am.put("UP", selectPrev);
    am.put("DOWN", selectNext);
  }

  private void showHideSamples() {
    DefaultListModel<JCheckBox> listModel = new DefaultListModel<>();
    for (int i = 0; i < proj.getNumberOfParsedSamples(); i++) {
      listModel.addElement(new JCheckBox(proj.getSamples()[i],
                                         hidden.contains(proj.getSamples()[i])));
    }
    JList<JCheckBox> samplesList = new JList<JCheckBox>();
    JScrollPane scroll = new JScrollPane(samplesList);
    samplesList.setCellRenderer(new JCheckBoxListCellRenderer(getFont()));
    samplesList.addMouseListener(new MouseAdapter() {

      int prevIndex = -1;

      public void mousePressed(MouseEvent e) {
        int index = samplesList.locationToIndex(e.getPoint());
        if (index != -1) {
          JCheckBox checkbox = (JCheckBox) samplesList.getModel().getElementAt(index);
          if (checkbox.isEnabled()) {
            boolean cl = e.getClickCount() >= 2;
            boolean li = prevIndex != -1 && index == prevIndex;
            if (cl || li) {
              checkbox.setSelected(!checkbox.isSelected());
            }
            prevIndex = index;
          }
          e.getComponent().repaint();
          checkbox.repaint();
          scroll.getViewport().repaint();
        }
      }
    });
    samplesList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    samplesList.setModel(listModel);

    JDialog dialog = new JDialog(this, true);
    JPanel pnlBtns = new JPanel();
    JButton btnOk = new JButton();
    btnOk.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        Map<String, Boolean> showHide = new HashMap<>();
        for (int i = 0; i < listModel.getSize(); i++) {
          showHide.put(proj.getSamples()[i], listModel.get(i).isSelected());
        }
        showHide.entrySet().stream().forEach(ent -> {
          if (ent.getValue().booleanValue()) {
            hideSample(ent.getKey());
          } else {
            showSample(ent.getKey());
          }
        });
        dialog.setVisible(false);
        dialog.dispose();
      }
    });
    btnOk.setText("OK");
    JButton btnCancel = new JButton();
    btnCancel.setAction(new AbstractAction() {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        dialog.setVisible(false);
        dialog.dispose();
      }
    });
    btnCancel.setText("Cancel");
    pnlBtns.add(btnOk);
    pnlBtns.add(btnCancel);

    dialog.getContentPane().setLayout(new MigLayout());
    dialog.getContentPane().add(new JLabel("Select Samples to Hide:"), "dock north");
    dialog.getContentPane().add(scroll, "dock center");
    dialog.getContentPane().add(pnlBtns, "dock south");
    dialog.setLocation(fldAddAnnot.getLocationOnScreen());
    dialog.setSize(400, 600);
    dialog.setVisible(true);
  }

  private void hideSample(String sample) {
    hidden.add(sample);
    samples.remove(sample);
    indivPanels.get(sample).setVisible(false);
    indivSeparators.get(sample).setVisible(false);
    if (selectedSample != null && selectedSample.equals(sample)) {
      indivPanels.get(sample).setBorder(null);
      this.selectedSample = null;
    }
  }

  private void showSample(String sample) {
    hidden.remove(sample);
    samples.add(sample);
    indivPanels.get(sample).setVisible(true);
    indivSeparators.get(sample).setVisible(true);
  }

  private void selectSample(String sample) {
    if (selectedSample != null && selectedSample.equals(sample)) {
      return;
    }
    if (selectedSample != null) {
      indivPanels.get(selectedSample).setBorder(null);
    }
    this.selectedSample = sample;
    indivPanels.get(this.selectedSample).setBorder(new SoftBevelBorder(BevelBorder.LOWERED));
    int h = indivPanels.get(this.selectedSample).getBounds().height;
    Rectangle rect = new Rectangle(0, 0, indivPanels.get(RegionReviewer.this.selectedSample)
                                                    .getBounds().width,
                                   h + 10);
    indivPanels.get(RegionReviewer.this.selectedSample).scrollRectToVisible(rect);
    refreshAnnotations();
  }

  public void displayRegion(Segment seg) {
    fldRegion.setText(seg.getUCSClocation());
    new Thread(new Runnable() {
      @Override
      public void run() {
        setRegion(seg);
        samplePanel.requestFocus();
      }
    }).start();
    repaint();
  }

  private void setRegion(Segment seg) {
    this.seg = seg;

    try {
      SwingUtilities.invokeAndWait(new Runnable() {
        @Override
        public void run() {
          loadingLabel.setVisible(true);
        }
      });
    } catch (InvocationTargetException | InterruptedException e) {
      loadingLabel.setVisible(true);
    }

    redraw();
  }

  private void redraw() {
    if (this.seg == null) return;
    List<Marker> markersInSeg = this.mds.getMarkersInSeg(seg).asList();
    MDL mdl = new MDL(proj, markersInSeg);
    MarkerData[] data = new MarkerData[markersInSeg.size()];
    for (int i = 0; i < data.length; i++) {
      if (!mdl.hasNext()) {
        // TODO error!
      }
      data[i] = mdl.next();
    }
    Map<String, Integer> inds = proj.getSampleIndices();
    for (int s = 0; s < samples.size(); s++) {
      String sample = samples.get(s);
      Map<Marker, Float> lrrs = new HashMap<>();
      Map<Marker, Float> bafs = new HashMap<>();
      for (int m = 0; m < data.length; m++) {
        lrrs.put(markersInSeg.get(m), data[m].getLRRs()[inds.get(sample)]);
        bafs.put(markersInSeg.get(m), data[m].getBAFs()[inds.get(sample)]);
      }
      lrrMedians.put(sample, ArrayUtils.median(lrrs.values()));
      lrrPanels.get(sample).setData(lrrs, transform);
      bafPanels.get(sample).setData(bafs, transform);
    }
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        scrollPane.getViewport().repaint();
        loadingLabel.setVisible(false);
        samplePanel.requestFocus();
      }
    });
  }

  private void reorder() {
    samplePanel.removeAll();
    switch (sorting) {
      case SAMPLE:
        Collections.sort(samples, sampleComparator);
        break;
      case LRRMEDIAN:
        Collections.sort(samples, lrrComparator);
    }
    for (String sample : samples) {
      samplePanel.add(indivPanels.get(sample), "grow");
      samplePanel.add(indivSeparators.get(sample), "growx");
    }
    samplePanel.revalidate();
    samplePanel.repaint();
  }

  public void setReviewPanelSize(int width, int height) {
    this.panelSize = new Dimension(width, height);
  }

  private void createNewAnnotation() {
    String newAnnotation = fldAddAnnot.getText().trim();
    fldAddAnnot.setText("");
    if ("".equals(newAnnotation)) {
      return;
    }
    addAnnotation(newAnnotation);
  }

  private Annotation addAnnotation(String newAnnotation) {
    Annotation ann = new Annotation(newAnnotation);
    if (!annotations.contains(ann)) {
      annotations.add(ann);
      addAnnotationToTraversalMenu(ann);
      refreshAnnotations();
    }
    return ann;
  }

  private void refreshAnnotations() {
    mnemonicActions.clear();
    ArrayList<Annotation> allAnnots = new ArrayList<>(annotations);
    annotPanel.removeAll();
    for (int i = 0; i < allAnnots.size(); i++) {
      JCheckBox annBox = new JCheckBox();
      Annotation ann = allAnnots.get(i);
      Collection<Annotation> aiAnns = selectedSample != null
                                      && annotationMap.containsKey(selectedSample) ? annotationMap.get(selectedSample)
                                                                                   : null;
      if (aiAnns != null && aiAnns.contains(ann)) {
        annBox.setSelected(true);
      } else {
        annBox.setSelected(false);
      }
      final int ind = i;
      AbstractAction mnemAct = new AbstractAction() {
        private static final long serialVersionUID = 1L;

        @Override
        public void actionPerformed(ActionEvent e) {
          if (selectedSample == null) return;
          Annotation a = allAnnots.get(ind);
          if (annotationMap.containsEntry(selectedSample, a)) {
            annotationMap.remove(selectedSample, a);
            annBox.setSelected(false);
          } else {
            annotationMap.put(selectedSample, a);
            annBox.setSelected(true);
          }
          setAnnotationsChanged();
          annotPanel.repaint();
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

  private void addAnnotationToTraversalMenu(Annotation annot) {
    JRadioButtonMenuItem jrb = new JRadioButtonMenuItem();
    jrb.setText(annot.annotation);
    annTravMap.put(annot, jrb);
    travGrp.add(jrb);
    mnTravAnn.add(jrb);
    JMenuItem jmn = new JMenuItem();
    jmn.setAction(actionSaveAnnotated);
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
  }

  private void deleteAnnotation(Annotation annot) {
    JRadioButtonMenuItem jrb = annTravMap.get(annot);
    travGrp.remove(jrb);
    mnTravAnn.remove(jrb);
    mnSaveAnn.remove(annSaveMap.get(annot));
    mnDelAnn.remove(annDelMap.get(annot));
    annotationMap.asMap().entrySet().stream().forEach(e -> e.getValue().remove(annot));
    setAnnotationsChanged();
    refreshAnnotations();
  }

  private void setAnnotationsChanged() {
    annotationsChanged = true;
  }

  private boolean checkAnnotationsChanged() {
    boolean value = annotationsChanged;
    annotationsChanged = false;
    return value;
  }

  private void saveAll(String file) {
    PrintWriter writer = Files.getAppropriateWriter(file);
    for (String s : proj.getSamples()) {
      writer.print(s + "\t" + seg.getUCSClocation() + "\t");
      Collection<Annotation> annots = annotationMap.get(s);
      int index = 0;
      for (Annotation annot : annots) {
        writer.print(annot.annotation);
        if (index < annots.size() - 1) {
          writer.print("|");
        }
        index++;
      }
      writer.println();
    }
    writer.close();
  }

  private void saveAllAnnotated(String file) {
    PrintWriter writer = Files.getAppropriateWriter(file);
    for (String s : proj.getSamples()) {
      Collection<Annotation> annots = annotationMap.get(s);
      if (annots.size() == 0) continue;
      writer.print(s + "\t" + seg.getUCSClocation() + "\t");
      int index = 0;
      for (Annotation annot : annots) {
        writer.print(annot.annotation);
        if (index < annots.size() - 1) {
          writer.print("|");
        }
        index++;
      }
      writer.println();
    }
    writer.close();
  }

  private void saveAnnotated(String file, Annotation annot) {
    PrintWriter writer = Files.getAppropriateWriter(file);
    for (String s : proj.getSamples()) {
      Collection<Annotation> annots = annotationMap.get(s);
      if (annots.size() == 0 || !annots.contains(annot)) continue;
      writer.print(s + "\t" + seg.getUCSClocation() + "\t");
      int index = 0;
      for (Annotation annot1 : annots) {
        writer.print(annot1.annotation);
        if (index < annots.size() - 1) {
          writer.print("|");
        }
        index++;
      }
      writer.println();
    }
    writer.close();
  }

  private void load(String file) throws IOException {
    BufferedReader reader = Files.getAppropriateReader(file);
    String line = null;
    String[] parts;
    Segment seg = null;
    while ((line = reader.readLine()) != null) {
      parts = line.split("\t", -1);
      String sample = parts[0];
      if (!samples.contains(sample) && !hidden.contains(sample)) {
        proj.getLog().reportTimeWarning("Unknown sample found, not present in project: " + sample);
        continue;
      }
      String region = parts[1];
      String annots = parts.length > 2 ? parts[2] : "";
      if (seg == null) {
        seg = new Segment(region);
      } else if (!seg.matches(new Segment(region))) {
        proj.getLog().reportError("Mismatch in regions; expected " + seg.getUCSClocation()
                                  + "; found " + region);
        continue;
      }
      if (!annots.isEmpty()) {
        int ind = annots.lastIndexOf('^');
        if (ind >= 0) {
          annots = annots.substring(ind);
        }
        String[] anns = annots.split("\\|");
        for (String ann : anns) {
          annotationMap.put(sample, addAnnotation(ann));
        }
      }
    }
    reader.close();
    displayRegion(seg);
  }

  private void fireMnem(KeyEvent e) {
    char code = (e.getKeyChar() + "").toUpperCase().charAt(0);
    if (mnemonicActions.containsKey(code)) {
      mnemonicActions.get(code).actionPerformed(null);
      setAnnotationsChanged();
      if (checkBxAdvanceAfterAnnotationKeyed && !e.isShiftDown()) {
        selectNext.actionPerformed(null);
      }
    }
  }

  class LRRPanel extends JPanel {
    private static final long serialVersionUID = 1L;
    private Map<Marker, Float> data;
    private int dispMin, dispMax;

    public void setData(Map<Marker, Float> data, TRANSFORMATION transform) {
      this.data = data;
      this.dispMin = transform.getDisplayMin();
      this.dispMax = transform.getDisplayMax();
    }

    @Override
    protected void paintComponent(Graphics g) {
      g.setColor(Color.LIGHT_GRAY);
      g.drawLine(WIDTH_BUFFER,
                 getHeight() - (((int) (.5 * (double) getHeight())) - (2 * HEIGHT_BUFFER))
                               - HEIGHT_BUFFER,
                 getWidth() - WIDTH_BUFFER, getHeight() - (((int) (.5 * (double) getHeight()))
                                                           - (2 * HEIGHT_BUFFER))
                                            - HEIGHT_BUFFER);
      if (data != null) {
        g.setColor(Color.black);
        for (Marker m : data.keySet()) {
          int x = getX(m.getPosition());
          if (data.get(m) < dispMin) {
            g.drawString("v", x,
                         getHeight()
                                 - (int) ((double) (dispMin - dispMin)
                                          / (double) (dispMax - dispMin)
                                          * (getHeight() - 2 * HEIGHT_BUFFER - 1))
                                 - HEIGHT_BUFFER - 1);
          } else if (data.get(m) > dispMax) {
            g.drawString("^", x,
                         getHeight()
                                 - (int) ((double) (dispMax - dispMin)
                                          / (double) (dispMax - dispMin)
                                          * (getHeight() - 2 * HEIGHT_BUFFER - 1))
                                 - HEIGHT_BUFFER - 1);
          } else {
            int y = getHeight() - (int) ((double) (data.get(m) - dispMin) / (dispMax - dispMin)
                                         * (getHeight() - 2 * HEIGHT_BUFFER - 1))
                    - HEIGHT_BUFFER - 1;
            g.fillOval(x, y, size, size);
          }
        }
      }
    }

    private int getX(int pos) {
      byte chr = seg.getChr();
      int start = seg.getStart();
      int stop = seg.getStop();
      if (start == -1 || start < 0) {
        start = 1;
      }
      if (stop == -1 || stop > markerChrMap.get(chr).last().getPosition()) {
        stop = markerChrMap.get(chr).last().getPosition();
      }
      if (start >= stop) {
        start = stop - 1;
      }
      return (int) ((double) (pos - start) / (double) (stop - start)
                    * (getWidth() - 2 * WIDTH_BUFFER))
             + WIDTH_BUFFER;
    }
  }

  class BAFPanel extends JPanel {
    private static final long serialVersionUID = 1L;
    private Map<Marker, Float> data;

    public void setData(Map<Marker, Float> data, TRANSFORMATION transform) {
      this.data = data;
    }

    @Override
    protected void paintComponent(Graphics g) {
      if (data != null) {
        g.setColor(Color.black);
        for (Marker marker : data.keySet()) {
          int x = getX(marker.getPosition());
          float yv = data.get(marker);
          int y = getHeight() - (int) (yv * (double) (getHeight() - 4 * HEIGHT_BUFFER - 1))
                  - 2 * HEIGHT_BUFFER;
          g.fillOval(x, y, size, size);
        }
      }
    }

    private int getX(int pos) {
      byte chr = seg.getChr();
      int start = seg.getStart();
      int stop = seg.getStop();
      if (start == -1 || start < 0) {
        start = 1;
      }
      if (stop == -1 || stop > markerChrMap.get(chr).last().getPosition()) {
        stop = markerChrMap.get(chr).last().getPosition();
      }
      if (start >= stop) {
        start = stop - 1;
      }
      return (int) ((double) (pos - start) / (double) (stop - start)
                    * (getWidth() - 2 * WIDTH_BUFFER))
             + WIDTH_BUFFER;
    }
  }

  public static class Annotation {

    final String annotation;
    final char mnemonic;

    public Annotation(String ann) {
      this.annotation = ann;
      this.mnemonic = determineMnemonic(ann);
    }

    private static final HashMap<String, Character> mnemonicMap = new HashMap<>();

    public static Character determineMnemonic(String ann) {
      Character mnemonic;
      if (mnemonicMap.containsKey(ann)) {
        mnemonic = mnemonicMap.get(ann);
      } else {
        int ind = 0;
        Character c = ann.toUpperCase().charAt(ind);
        while (mnemonicMap.containsValue(c)) {
          ind++;
          c = ann.toUpperCase().charAt(ind);
        }
        if (ind == ann.length()) {
          System.err.println("Error - all possible mnemonic characters already used for annotation {"
                             + ann + "}.  Using alphanumerics instead.");
          String alphanum = "abcdefghijklmnopqrstuvwxyz0123456789";
          ind = 0;
          while (mnemonicMap.containsValue(alphanum.charAt(ind))) {
            ind++;
            if (ind == alphanum.length()) {
              System.err.println("Error - ran out of alphanumeric mnemonics too!");
              break;
            }
          }
          if (ind < alphanum.length()) {
            mnemonic = alphanum.charAt(ind);
            mnemonicMap.put(ann, mnemonic);
          } else {
            mnemonic = '-';
          }
        } else {
          mnemonicMap.put(ann, c);
          mnemonic = c;
        }
      }
      return mnemonic;
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((annotation == null) ? 0 : annotation.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      Annotation other = (Annotation) obj;
      if (annotation == null) {
        if (other.annotation != null) return false;
      } else if (!annotation.equals(other.annotation)) return false;
      return true;
    }
  }

  public static void main(String[] args) {
    Project proj = new Project("G:\\WorkFiles\\projects\\Ovation-P2-MSI.properties");
    RegionReviewer rr = new RegionReviewer(proj);
    rr.displayRegion(new Segment("chr10:135096576-135406504"));
    rr.repaint();
  }

}
