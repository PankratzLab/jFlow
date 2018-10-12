package org.genvisis.link;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.BorderFactory;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import org.genvisis.cnv.gui.JPanelFlowLayoutComponentListener;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ProgressBarDialog;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.PSF.Colors.BLUES;
import org.pankratzlab.common.PSF.Colors.GREENS;
import org.pankratzlab.common.PSF.Colors.REDS;
import org.pankratzlab.common.PSF.Colors.VIOLETS;
import org.pankratzlab.phenoprep.LinkageMap;
import org.pankratzlab.shared.gui.UITools;

public class PlotResults extends JFrame implements ActionListener {

  public static final long serialVersionUID = 1L;
  public static final boolean PLOT_INFO = true;

  public static final String[] RESULT_TYPES = {"Allegro", "Merlin", "Merlin-VC", "Merlin-Regress"};
  public static final String[][] RESULT_TYPE_NAMES = {{"chrom\\d\\d.lin.out",
                                                       "chrom\\d\\d.exp.out"},
                                                      {"merlin-chr\\d\\d-nonparametric.tbl",
                                                       "merlin-chr\\d\\d-info.tbl"},
                                                      {"vc-chr\\d\\d-vc-chr\\d\\d.tbl",
                                                       "vc-chr\\d\\d-info.tbl"},
                                                      {"regress-chr\\d\\d-regress-chr\\d\\d.tbl",
                                                       "regress-chr\\d\\d-regress-chr\\d\\d.tbl"}};
  public static final String[][][] RESULT_TYPE_HEADERS = {{{"location", "LOD", "dhat", "NPL", "Zlr",
                                                            "marker"},
                                                           {"location", "LOD", "dhat", "NPL", "Zlr",
                                                            "info", "marker"}},
                                                          {{"CHR", "POS", "LABEL", "ANALYSIS",
                                                            "ZSCORE", "DELTA", "LOD", "PVALUE"},
                                                           {"CHR", "POS", "LABEL", "INFO"}},
                                                          {{"CHR", "POS", "LABEL", "TRAIT", "H2",
                                                            "LOD", "PVALUE"},
                                                           {"CHR", "POS", "LABEL", "INFO"}},
                                                          {{"CHR", "POS", "PHENOTYPE", "H2", "SD",
                                                            "INFO", "LOD", "PVALUE"},
                                                           {"CHR", "POS", "PHENOTYPE", "H2", "SD",
                                                            "INFO", "LOD", "PVALUE"}}};
  public static final int[][] RESULT_TYPE_LOD_INDICES = {{0, 1, 5}, {1, 6, 2}, {1, 5, 2},
                                                         {1, 6, 1}};
  public static final int[][] RESULT_TYPE_INFO_INDICES = {{0, 5}, {1, 3}, {1, 3}, {1, 5}};
  public static final int[] RESULT_TYPE_SKIP_EXTRA_LINES = {0, 2, 0, 0};
  public static final Color BACKGROUND_COLOR = Color.WHITE;

  public static final int HEAD_BUFFER = 25;
  public static final int HEIGHT_X_AXIS = 55;
  public static final int WIDTH_Y_AXIS = 50;
  public static final int WIDTH_INFO_AXIS = 70;
  public static final int SIZE = 8;
  public static final int AXIS_THICKNESS = 4;
  public static final int TICK_THICKNESS = 3;
  public static final int TICK_LENGTH = 15;
  public static final int MARKER_SIZE = 6;
  public static final int NUM_CHR = 23;
  public static final String[] MARKER_POS_HEADER = {"Chr", "Marker", "Sex_Ave_cM", "Used"};
  private static final String FIRST = "First";
  private static final String PREVIOUS = "Previous";
  private static final String NEXT = "Next";
  private static final String LAST = "Last";
  private static final String CHECK_ALL = "Check all";
  private static final String UNCHECK_ALL = "Uncheck all";

  private static final Color[] DEFAULT_COLOR_SCHEME = {BLUES.MIDNIGHT_EXPRESS, // dark dark
                                                       BLUES.PERSIAN_BLUE, // dark blue
                                                       REDS.VENETIAN_RED, // deep red
                                                       VIOLETS.BLUE_VIOLET, // deep purple
                                                       GREENS.GREEN, // dark green
                                                       BLUES.DODGER_BLUE, // light blue
                                                       VIOLETS.ORCHID, // pink
                                                       BLUES.SLATE_BLUE, // light purple
                                                       GREENS.GREEN_YELLOW, // light green

                                                       BLUES.MIDNIGHT_EXPRESS, // dark dark
                                                       BLUES.PERSIAN_BLUE, // dark blue
                                                       REDS.VENETIAN_RED, // deep red
                                                       VIOLETS.BLUE_VIOLET, // deep purple
                                                       GREENS.GREEN, // dark green
                                                       BLUES.DODGER_BLUE, // light blue
                                                       VIOLETS.ORCHID, // pink
                                                       BLUES.SLATE_BLUE, // light purple
                                                       GREENS.GREEN_YELLOW, // light green

                                                       BLUES.MIDNIGHT_EXPRESS, // dark dark
                                                       BLUES.PERSIAN_BLUE, // dark blue
                                                       REDS.VENETIAN_RED, // deep red
                                                       VIOLETS.BLUE_VIOLET, // deep purple
                                                       GREENS.GREEN, // dark green
                                                       BLUES.DODGER_BLUE, // light blue
                                                       VIOLETS.ORCHID, // pink
                                                       BLUES.SLATE_BLUE, // light purple
                                                       GREENS.GREEN_YELLOW, // light green
  };

  private JButton first, previous, next, last;
  private double[][][][] data;
  private double[][][][] info;
  private String[][][] markerNames;
  private double[][][] markerPositions;
  private int[][][] markerUsage;
  private double[][] minMaxes;
  private int chr;
  private boolean[] inUse;
  private String[] dirNames;
  private JCheckBox[] dirBoxes;
  private JCheckBox alignBox;
  private JCheckBox infoBox;
  private JLabel chrNum;
  private JLabel commentLabel;
  private boolean alignMaps;
  private boolean plotInformativeness;

  private final Color[] colorScheme = DEFAULT_COLOR_SCHEME;

  public PlotResults(String title, String[] dirs, double[][][][] newData, double[][][][] newInfo,
                     String[][][][] newMarkerInfo) {
    super(title == null ? "Linkage Plots" : title);
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    JPanel panel;
    int[] keys;
    double[] positions;

    if (dirs.length == 0) {
      System.err.println("Error - nothing to plotl directory is empty");
      return;
    }

    data = newData;
    info = newInfo;
    dirNames = dirs;
    inUse = ArrayUtils.booleanArray(dirs.length, true);

    markerNames = new String[dirs.length][NUM_CHR][];
    markerPositions = new double[dirs.length][NUM_CHR][];
    markerUsage = new int[dirs.length][NUM_CHR][];
    for (int i = 0; i < dirs.length; i++) {
      for (int j = 0; j < NUM_CHR; j++) {
        if (newMarkerInfo[i][j] != null) {
          positions = new double[newMarkerInfo[i][j].length];
          for (int k = 0; k < newMarkerInfo[i][j].length; k++) {
            positions[k] = Double.parseDouble(newMarkerInfo[i][j][k][1]);
          }
          keys = Sort.getSortedIndices(positions);
          markerNames[i][j] = new String[positions.length];
          markerPositions[i][j] = new double[positions.length];
          markerUsage[i][j] = new int[positions.length];
          for (int k = 0; k < positions.length; k++) {
            markerNames[i][j][k] = newMarkerInfo[i][j][keys[k]][0];
            markerPositions[i][j][k] = positions[keys[k]] - positions[keys[0]];
            markerUsage[i][j][k] = Integer.parseInt(newMarkerInfo[i][j][keys[k]][2]);
          }
        }
      }
    }

    minMaxes = new double[NUM_CHR][4];
    for (int i = 0; i < NUM_CHR; i++) {
      for (int j = 0; j < data.length; j++) {
        if (markerPositions[j][i] != null) {
          minMaxes[i][1] = Math.max(markerPositions[j][i][markerPositions[j][i].length - 1],
                                    minMaxes[i][1]);
        }
        for (int k = 0; data[j][i] != null && k < data[j][i].length; k++) {
          minMaxes[i][0] = Math.min(data[j][i][k][0], minMaxes[i][0]);
          minMaxes[i][1] = Math.max(data[j][i][k][0], minMaxes[i][1]);
          minMaxes[i][2] = Math.min(data[j][i][k][1], minMaxes[i][2]);
          minMaxes[i][3] = Math.max(data[j][i][k][1], minMaxes[i][3]);
        }
      }
      minMaxes[i][3] = Math.max(Math.ceil(minMaxes[i][3]), 3);
    }

    JPanel graphPanel = new JPanel();
    graphPanel.setLayout(new BorderLayout());
    graphPanel.setBackground(BACKGROUND_COLOR);

    panel = new JPanel() {

      public static final long serialVersionUID = 3L;

      private int xMin = 0;
      private int xMax = getWidth();
      private int yMin = 0;
      private int yMax = getWidth();
      private double plotXmax;

      @Override
      public void paintComponent(Graphics g) {
        String str;
        int step = 10;
        int indexSet;

        while (minMaxes[chr - 1][1] / step > 10) {
          step += 10;
        }
        plotXmax = Math.ceil(minMaxes[chr - 1][1] / step) * step;

        // plot
        xMin = WIDTH_Y_AXIS;
        xMax = getWidth() - WIDTH_INFO_AXIS;
        yMin = HEIGHT_X_AXIS;
        yMax = getHeight() - HEAD_BUFFER;
        for (int i = 0; i < data.length; i++) {
          if (inUse[i] && data[i][chr - 1] != null) {
            g.setColor(colorScheme[i]);
            for (int j = 0; j < data[i][chr - 1].length - 1; j++) {
              Grafik.drawThickLine(g,
                                   getX(data[alignMaps && data[0][chr - 1] != null ? 0
                                                                                   : i][chr
                                                                                        - 1][j][0]),
                                   getY(data[i][chr - 1][j][1]),
                                   getX(data[alignMaps
                                             && data[0][chr - 1] != null ? 0 : i][chr - 1][j
                                                                                           + 1][0]),
                                   getY(data[i][chr - 1][j + 1][1]), SIZE / 2, colorScheme[i]);
            }
          }
          if (plotInformativeness && inUse[i] && info[i][chr - 1] != null) {
            g.setColor(colorScheme[i]);
            for (int j = 0; j < info[i][chr - 1].length - 1; j++) {
              Grafik.drawThickLine(g,
                                   getX(info[alignMaps && info[0][chr - 1] != null ? 0
                                                                                   : i][chr
                                                                                        - 1][j][0]),
                                   getNormY(info[i][chr - 1][j][1]),
                                   getX(info[alignMaps
                                             && info[0][chr - 1] != null ? 0 : i][chr - 1][j
                                                                                           + 1][0]),
                                   getNormY(info[i][chr - 1][j + 1][1]), SIZE / 4, colorScheme[i]);
            }
          }
        }

        // x-Axis
        xMin = WIDTH_Y_AXIS;
        xMax = getWidth() - WIDTH_INFO_AXIS;
        yMin = 0;
        yMax = HEIGHT_X_AXIS;
        g.setFont(new Font("Arial", 0, 28));
        for (int i = 0; i <= plotXmax / step; i++) {
          Grafik.drawThickLine(g, getX(i * step), getHeight() - yMax, getX(i * step),
                               getHeight() - (yMax - TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
          g.drawString(i * step + "", getX(i * step) - (i * step + "").length() * 8,
                       getHeight() - (yMax - TICK_LENGTH - 30));
        }
        Grafik.drawThickLine(g, xMin - (int) Math.ceil(AXIS_THICKNESS / 2.0), getHeight() - yMax,
                             xMax + (int) Math.ceil(AXIS_THICKNESS / 2.0), getHeight() - yMax,
                             AXIS_THICKNESS, Color.BLACK);

        // y-axis
        xMin = 0;
        xMax = WIDTH_Y_AXIS;
        yMin = HEIGHT_X_AXIS;
        yMax = getHeight() - HEAD_BUFFER;
        g.setFont(new Font("Arial", 0, 28));
        for (int i = 0; i <= minMaxes[chr - 1][3]; i++) {
          Grafik.drawThickLine(g, xMax - TICK_LENGTH, getY(i), xMax, getY(i), TICK_THICKNESS,
                               Color.BLACK);
          g.drawString(i + "", xMax - TICK_LENGTH - (i + "").length() * 15 - 5, getY(i) + 9);
        }
        Grafik.drawThickLine(g, xMax, getY(minMaxes[chr - 1][2]), xMax,
                             getY(minMaxes[chr - 1][3]) - (int) Math.ceil(TICK_THICKNESS / 2.0),
                             AXIS_THICKNESS, Color.BLACK);

        // info axis
        xMin = getWidth() - WIDTH_INFO_AXIS;
        xMax = getWidth();
        yMin = HEIGHT_X_AXIS;
        yMax = getHeight() - HEAD_BUFFER;
        g.setFont(new Font("Arial", 0, 28));
        for (double i = 0; i <= 1; i += 0.2) {
          Grafik.drawThickLine(g, xMin, getNormY(i), xMin + TICK_LENGTH, getNormY(i),
                               TICK_THICKNESS, Color.BLACK);
          str = ext.formDeci(i, 1, true);
          g.drawString(str, xMin + TICK_LENGTH + 5, getNormY(i) + 9);
        }
        Grafik.drawThickLine(g, xMin, getNormY(0), xMin,
                             getNormY(1) - (int) Math.ceil(TICK_THICKNESS / 2.0), AXIS_THICKNESS,
                             Color.BLACK);

        // markers
        xMin = WIDTH_Y_AXIS;
        xMax = getWidth() - WIDTH_INFO_AXIS;
        indexSet = -1;
        for (int i = 0; i < dirNames.length; i++) {
          if (inUse[i] && indexSet == -1 && markerPositions[i][chr - 1] != null) {
            indexSet = i;
          }
        }
        if (indexSet != -1) {
          g.setColor(colorScheme[indexSet]);
          for (int i = 0; i < markerUsage[indexSet][chr - 1].length; i++) {
            if (markerUsage[indexSet][chr - 1][i] == 1) {
              g.fillOval(getX(markerPositions[indexSet][chr - 1][i]), HEAD_BUFFER, MARKER_SIZE,
                         MARKER_SIZE);
            } else if (markerUsage[indexSet][chr - 1][i] == -1) {
              g.drawOval(getX(markerPositions[indexSet][chr - 1][i]), HEAD_BUFFER, MARKER_SIZE,
                         MARKER_SIZE);
            } else {
              g.setFont(new Font("Arial", Font.BOLD, 15));
              g.drawString("X", getX(markerPositions[indexSet][chr - 1][i]) - 3, HEAD_BUFFER + 9);
            }
          }

        }
      }

      public int getX(double x) {
        return (int) ((x - minMaxes[chr - 1][0]) / (plotXmax - minMaxes[chr - 1][0])
                      * (xMax - xMin))
               + xMin;
      }

      public int getY(double y) {
        return getHeight()
               - (int) ((y - minMaxes[chr - 1][2]) / (minMaxes[chr - 1][3] - minMaxes[chr - 1][2])
                        * (yMax - yMin) + yMin);
      }

      public int getNormY(double y) {
        return getHeight() - (int) (y * (yMax - yMin) + yMin);
      }
    };
    graphPanel.add(panel, BorderLayout.CENTER);

    graphPanel.add(new JLabel(Grafik.getImageIcon("images/distance28.gif")), BorderLayout.SOUTH);
    graphPanel.add(new JLabel(Grafik.getImageIcon("images/LODscore28+73.gif")), BorderLayout.WEST);
    graphPanel.add(new JLabel(Grafik.getImageIcon("images/information20.gif")), BorderLayout.EAST);
    getContentPane().add(graphPanel, BorderLayout.CENTER);

    JPanel descrPanel = new JPanel();
    descrPanel.setLayout(new GridLayout(2, 1));

    JPanel navigationPanel = new JPanel();
    // navigationPanel.setLayout(new GridLayout(1,5));

    first = new JButton(Grafik.getImageIcon("images/firstLast/First.gif"));
    first.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif"));
    first.addActionListener(this);
    first.setActionCommand(FIRST);
    first.setPreferredSize(new Dimension(20, 20));

    previous = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif"));
    previous.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif"));
    previous.addActionListener(this);
    previous.setActionCommand(PREVIOUS);
    previous.setPreferredSize(new Dimension(20, 20));

    chrNum = new JLabel("Chromosome #", JLabel.CENTER);
    chrNum.setFont(new Font("Arial", 0, 20));
    chrNum.setPreferredSize(new Dimension(160, 20));

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
    navigationPanel.add(chrNum);
    navigationPanel.add(next);
    navigationPanel.add(last);
    navigationPanel.setBackground(BACKGROUND_COLOR);
    descrPanel.add(navigationPanel);

    commentLabel = new JLabel("", JLabel.CENTER);
    commentLabel.setFont(new Font("Arial", 0, 14));
    descrPanel.add(commentLabel);
    descrPanel.setBackground(BACKGROUND_COLOR);

    getContentPane().add(descrPanel, BorderLayout.NORTH);

    JPanel classPanel = new JPanel();
    classPanel.setLayout(new FlowLayout());
    classPanel.addComponentListener(new JPanelFlowLayoutComponentListener());

    // filterPanel.setLayout(new BoxLayout(filterPanel,
    // BoxLayout.PAGE_AXIS));
    JLabel label = new JLabel("Directories with data:");
    label.setFont(new Font("Arial", 0, 20));
    classPanel.add(label);

    ItemListener classListener = new ItemListener() {

      @Override
      public void itemStateChanged(ItemEvent ie) {
        inUse[ext.indexOfStr(((JCheckBox) ie.getSource()).getText(),
                             dirNames)] = ((JCheckBox) ie.getSource()).isSelected();
        updateGUI();
      }
    };
    dirBoxes = new JCheckBox[dirNames.length];
    boolean alignable = true;
    int[] chrLengths = ArrayUtils.intArray(NUM_CHR, -1);
    for (int i = 0; i < dirNames.length; i++) {
      dirBoxes[i] = new JCheckBox(dirNames[i]);
      dirBoxes[i].setFont(new Font("Arial", 0, 14));
      dirBoxes[i].setSelected(inUse[i]);
      dirBoxes[i].addItemListener(classListener);
      dirBoxes[i].setBorder(BorderFactory.createLineBorder(colorScheme[i], 5));
      dirBoxes[i].setBorderPainted(true);
      dirBoxes[i].setBackground(BACKGROUND_COLOR);
      classPanel.add(dirBoxes[i]);
      for (int j = 0; j < NUM_CHR; j++) {
        if (data[i][j] != null) {
          if (chrLengths[j] == -1) {
            chrLengths[j] = data[i][j].length;
          } else if (data[i][j].length != chrLengths[j]) {
            alignable = false;
          }
        }
      }
    }
    if (data.length > 1 && alignable) {
      label = new JLabel("Check to align maps:");
      label.setFont(new Font("Arial", 0, 20));
      classPanel.add(label);
      alignBox = new JCheckBox("");
      alignBox.setSelected(true);
      alignBox.setBackground(BACKGROUND_COLOR);
      alignMaps = true;
      alignBox.addItemListener(new ItemListener() {

        @Override
        public void itemStateChanged(ItemEvent ie) {
          alignMaps = ((JCheckBox) ie.getSource()).isSelected();
          updateGUI();
        }
      });
      classPanel.add(alignBox);
    }
    classPanel.setBackground(BACKGROUND_COLOR);

    label = new JLabel("Plot informativeness: ");
    label.setFont(new Font("Arial", 0, 20));
    classPanel.add(label);
    infoBox = new JCheckBox("");
    // infoBox.setSelected(true);
    // plotInformativeness = true;
    infoBox.setSelected(false);
    plotInformativeness = false;
    infoBox.setBackground(BACKGROUND_COLOR);
    infoBox.addItemListener(new ItemListener() {

      @Override
      public void itemStateChanged(ItemEvent ie) {
        plotInformativeness = ((JCheckBox) ie.getSource()).isSelected();
        updateGUI();
      }
    });
    classPanel.add(infoBox);

    JButton button = new JButton(CHECK_ALL);
    button.addActionListener(this);
    button.setActionCommand(CHECK_ALL);
    button.setPreferredSize(new Dimension(120, 25));
    classPanel.add(button);

    button = new JButton(UNCHECK_ALL);
    button.addActionListener(this);
    button.setActionCommand(UNCHECK_ALL);
    button.setPreferredSize(new Dimension(120, 25));
    classPanel.add(button);

    getContentPane().add(classPanel, BorderLayout.SOUTH);

    InputMap inputMap = panel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_HOME, InputEvent.CTRL_MASK), FIRST);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, InputEvent.CTRL_MASK), PREVIOUS);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, InputEvent.CTRL_MASK), NEXT);
    inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_END, InputEvent.CTRL_MASK), LAST);
    ActionMap actionMap = panel.getActionMap();
    actionMap.put(FIRST, new AbstractAction() {

      public static final long serialVersionUID = 4L;

      @Override
      public void actionPerformed(ActionEvent e) {
        chr = 1;
        checkButtons();
        updateGUI();
      }
    });
    actionMap.put(PREVIOUS, new AbstractAction() {

      public static final long serialVersionUID = 5L;

      @Override
      public void actionPerformed(ActionEvent e) {
        chr = Math.max(chr - 1, 1);
        checkButtons();
        updateGUI();
      }
    });
    actionMap.put(NEXT, new AbstractAction() {

      public static final long serialVersionUID = 6L;

      @Override
      public void actionPerformed(ActionEvent e) {
        chr = Math.min(chr + 1, NUM_CHR);
        checkButtons();
        updateGUI();
      }
    });
    actionMap.put(LAST, new AbstractAction() {

      public static final long serialVersionUID = 7L;

      @Override
      public void actionPerformed(ActionEvent e) {
        chr = NUM_CHR;
        checkButtons();
        updateGUI();
      }
    });
    panel.setActionMap(actionMap);

    next.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT);
    next.setActionMap(actionMap);
    previous.setActionMap(actionMap);

    chr = 16;
    updateGUI();

    setMinimumSize(new Dimension(20, 20));

    UITools.setSize(this, new Dimension(1000, 720));
    pack();
    setVisible(true);
  }

  public void updateGUI() {
    chrNum.setText("Chromosome " + chr);
    for (int i = 0; i < data.length; i++) {
      dirBoxes[i].setEnabled(data[i][chr - 1] != null);
    }
    repaint();
  }

  @Override
  public void actionPerformed(ActionEvent ae) {
    String command = ae.getActionCommand();

    if (command.equals(FIRST)) {
      chr = 1;
      updateGUI();
    } else if (command.equals(PREVIOUS)) {
      chr = Math.max(chr - 1, 1);
      updateGUI();
    } else if (command.equals(NEXT)) {
      chr = Math.min(chr + 1, NUM_CHR);
      updateGUI();
    } else if (command.equals(LAST)) {
      chr = NUM_CHR;
      updateGUI();
    } else if (command.equals(CHECK_ALL)) {
      for (int i = 0; i < dirNames.length; i++) {
        inUse[i] = true;
        dirBoxes[i].setSelected(inUse[i]);
      }
      updateGUI();
    } else if (command.equals(UNCHECK_ALL)) {
      for (int i = 0; i < dirNames.length; i++) {
        inUse[i] = false;
        dirBoxes[i].setSelected(inUse[i]);
      }
      updateGUI();
    } else {
      System.err.println("Error - unknown command '" + command + "'");
    }

    checkButtons();
  }

  public void checkButtons() {
    first.setEnabled(chr != 1);
    previous.setEnabled(chr != 1);
    next.setEnabled(chr != NUM_CHR);
    last.setEnabled(chr != NUM_CHR);
  }

  public static void loadResults(String dir) {
    BufferedReader reader;
    String[] line, dirs, files;
    double[][][][] data, info;
    String[][][][] markerData;
    IntVector dataIV, infoIV, mapIV;
    Vector<double[]> result;
    Hashtable<String, Vector<String>> hash;
    LinkageMap map;
    String[] markerNames;
    double[] dists, infoPair;
    Vector<String> v;
    int resultType;
    String filename;
    ProgressBarDialog prog;
    Dimension screenSize;
    String title;

    dirs = Files.listDirectories(dir);
    screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    prog = new ProgressBarDialog("Loading results...", 0, dirs.length * 23, screenSize.width,
                                 screenSize.height, 0);

    if (dirs.length > DEFAULT_COLOR_SCHEME.length) {
      System.err.println("Error - there are currently more data sources than there are colors defined!");
      return;
    }

    data = new double[dirs.length][NUM_CHR][][];
    info = new double[dirs.length][NUM_CHR][][];
    markerData = new String[dirs.length][NUM_CHR][][];
    for (int i = 0; i < dirs.length; i++) {
      files = Files.list(dir + dirs[i], "");
      resultType = -1;
      for (int j = 0; resultType < 0 && j < files.length; j++) {
        for (int k = 0; resultType < 0 && k < RESULT_TYPES.length; k++) {
          Pattern pattern = Pattern.compile(RESULT_TYPE_NAMES[k][0]);
          Matcher matcher = pattern.matcher(files[j]);
          if (matcher.find()) {
            resultType = k;
          }
        }
      }
      if (resultType == -1) {
        System.err.println("Error - could not determine which program generated the results in "
                           + dirs[i] + "/");
      } else {
        System.out.println("Parsing " + RESULT_TYPES[resultType] + " data from " + dirs[i] + "/");
        dataIV = new IntVector();
        infoIV = new IntVector();
        mapIV = new IntVector();
        hash = new Hashtable<>();
        for (int chr = 1; chr <= NUM_CHR; chr++) {
          prog.setProgress(i * 23 + chr);

          filename = ext.replaceAllWith(RESULT_TYPE_NAMES[resultType][0], "\\d\\d",
                                        ext.chrome(chr));
          try {
            reader = Files.getReader(dir + dirs[i] + "/" + filename, false, false);
            if (reader == null) {
              throw new FileNotFoundException(filename);
            }
            result = new Vector<>();
            ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                            RESULT_TYPE_HEADERS[resultType][0], true);
            for (int j = 0; j < RESULT_TYPE_SKIP_EXTRA_LINES[resultType]; j++) {
              reader.readLine();
            }
            while (reader.ready()) {
              line = reader.readLine().trim()
                           .split(filename.endsWith(".tbl") ? "\t" : PSF.Regex.GREEDY_WHITESPACE);
              result.add(new double[] {Double.parseDouble(line[RESULT_TYPE_LOD_INDICES[resultType][0]]),
                                       Double.parseDouble(line[RESULT_TYPE_LOD_INDICES[resultType][1]])});
              if (!line[RESULT_TYPE_LOD_INDICES[resultType][2]].equals("-")
                  && !line[RESULT_TYPE_LOD_INDICES[resultType][2]].equals(line[RESULT_TYPE_LOD_INDICES[resultType][0]])) {
                HashVec.addToHashVec(hash, chr + "",
                                     line[RESULT_TYPE_LOD_INDICES[resultType][2]] + "\t"
                                                     + line[RESULT_TYPE_LOD_INDICES[resultType][0]]
                                                     + "\t" + "1",
                                     true);
              }
            }
            reader.close();
            data[i][chr - 1] = Matrix.toDoubleArrays(result);
          } catch (FileNotFoundException fnfe) {
            dataIV.add(chr);
          } catch (NullPointerException fnfe) {
            System.err.println("Error reading file \"" + dirs[i] + "/" + filename + "\"");
            System.exit(2);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + dirs[i] + "/" + filename + "\"");
            System.exit(2);
          }
          filename = ext.replaceAllWith(RESULT_TYPE_NAMES[resultType][1], "\\d\\d",
                                        ext.chrome(chr));
          try {
            reader = Files.getReader(dir + dirs[i] + "/" + filename, false, false);
            result = new Vector<>();
            ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                            RESULT_TYPE_HEADERS[resultType][1], true);
            while (reader.ready()) {
              line = reader.readLine().trim()
                           .split(filename.endsWith(".tbl") ? "\t" : PSF.Regex.GREEDY_WHITESPACE);
              infoPair = new double[] {Double.parseDouble(line[RESULT_TYPE_INFO_INDICES[resultType][0]]),
                                       Double.parseDouble(line[RESULT_TYPE_INFO_INDICES[resultType][1]])};
              if (infoPair[1] > 1) {
                infoPair[1] /= 100;
              }
              result.add(infoPair);
            }
            reader.close();
            info[i][chr - 1] = Matrix.toDoubleArrays(result);
          } catch (FileNotFoundException fnfe) {
            infoIV.add(chr);
          } catch (NullPointerException fnfe) {
            infoIV.add(chr);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + dirs[i] + "/" + filename + "\"");
            System.exit(2);
          }
        }
        if (dataIV.size() > 0) {
          System.out.println("  Missing results for chromosomes: "
                             + ArrayUtils.toStr(dataIV.toArray(), ", "));
        }
        if (infoIV.size() > 0) {
          System.out.println("  Missing info for chromosomes: "
                             + ArrayUtils.toStr(infoIV.toArray(), ", "));
        }

        if (Files.exists(dir + dirs[i] + "/markerPositions.dat")) {
          hash = new Hashtable<>();
          try {
            reader = Files.getReader(dir + dirs[i] + "/markerPositions.dat", true, false);
            ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                            MARKER_POS_HEADER, true);
            while (reader.ready()) {
              line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
              HashVec.addToHashVec(hash, line[0], line[1] + "\t" + line[2] + "\t" + line[3], true);
            }
            reader.close();
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + dirs[i] + "/markerPositions.dat"
                               + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + dirs[i] + "/markerPositions.dat" + "\"");
            System.exit(2);
          }
        } else {
          for (int chr = 1; chr <= NUM_CHR; chr++) {
            if (new File(dir + dirs[i] + "/map" + ext.chrome(chr) + ".dat").exists()) {
              map = new LinkageMap(dir + dirs[i] + "/", chr);
              markerNames = map.getMarkerNames();
              dists = map.getCumulativePositions(false);
              for (int j = 0; j < markerNames.length; j++) {
                HashVec.addToHashVec(hash, chr + "", markerNames[j] + "\t" + dists[j] + "\t" + "1",
                                     true);
              }
            } else {
              mapIV.add(chr);
            }
          }
        }
        for (int chr = 1; chr <= NUM_CHR; chr++) {
          v = hash.get(chr + "");
          if (v != null) {
            markerData[i][chr - 1] = new String[v.size()][];
            for (int j = 0; j < v.size(); j++) {
              markerData[i][chr - 1][j] = v.elementAt(j).split(PSF.Regex.GREEDY_WHITESPACE);
            }
          }
        }
      }
    }
    prog.close();

    if (Files.exists(dir + "title.txt")) {
      title = HashVec.loadFileToStringArray(dir + "title.txt", false, new int[] {0}, true, false,
                                            "\t")[0];
    } else {
      title = null;
    }
    new PlotResults(title, dirs, data, info, markerData);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\hearing\\04_hearing_dx\\results\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\results\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\6K
    // Screens\\08_AltPheno\\_PDish_LargerFams\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\6K
    // Screens\\08_AltPheno\\_PDish_SmallFriesi\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\6K
    // Screens\\08_AltPheno\\_VPD\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\6K
    // Screens\\08_AltPheno\\allFams\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\UMN\\Pankow\\Linkage\\results\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\BOSS\\Linkage\\results\\";
    String dir = "D:/BOSS/Linkage/results/";
    // String dir = "D:/BOSS/Linkage/PCA_all_files/results/";

    String usage = "\\n" + "link.PlotResults requires 0-1 arguments\n"
                   + "   (0) directory (i.e. dir=" + dir + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-notJar")) {
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      loadResults(dir);
      // loadAllegroResults(dir, jar);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
