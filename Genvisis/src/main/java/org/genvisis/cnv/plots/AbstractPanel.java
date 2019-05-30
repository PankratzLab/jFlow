package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.Timer;

import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ProgressBarDialog;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.mining.Distance;
import org.pankratzlab.common.stats.Maths;

import com.google.common.primitives.Bytes;

public abstract class AbstractPanel extends JPanel implements MouseListener, MouseMotionListener,
                                    MouseWheelListener, ComponentListener, ActionListener {

  public static final long serialVersionUID = 1L;

  public static final boolean DEBUGGING = false;

  public static final int HEAD_BUFFER = 25;
  public static final int HEIGHT_X_AXIS = 90;
  public static final int WIDTH_Y_AXIS = 120;
  public static final int WIDTH_BUFFER = 50;
  public static final int AXIS_THICKNESS = 4;
  public static final int TICK_THICKNESS = 3;
  public static final int TICK_LENGTH = 15;
  public static final int DEFAULT_LOOKUP_RESOLUTION = 20;
  public static final int AXIS_FONT_SIZE = 28;
  public static final int DEFAULT_X_AXIS_LABEL_PAD = 11;
  public static final double DOUBLE_INACCURACY_HEDGE = 0.00001;
  public static final double MINIMUM_ZOOM_PROPORTION_WINDOW = 0.0001;
  public static final float DEFAULT_MOUSE_WHEEL_MULTIPLIER = 0.5f;
  public static final int DEFAULT_PLOTPOINTSET_SIZE = 10;
  public static final int SIZE = 12;
  public static final double HIGHLIGHT_DISTANCE = 20;// = Math.sqrt(SIZE*SIZE/2);
  public static final int DELAY = 0; // A control variable to reduce the repaint() operations during
                                     // component resizing;
  public static final int SCATTER_PLOT_TYPE = 1;
  public static final int HEAT_MAP_TYPE = 2;

  public static final int DEFAULT_TYPE = SCATTER_PLOT_TYPE;

  public static final int IMAGE_NULL = 0;
  public static final int IMAGE_STARTED = 1;
  public static final int IMAGE_COMPLETE = 2;

  protected Color[] colorScheme;
  protected int canvasSectionMinimumX;
  protected int canvasSectionMaximumX;
  protected int canvasSectionMinimumY;
  protected int canvasSectionMaximumY;
  private int axisYWidth = WIDTH_Y_AXIS;
  private int axisXHeight = HEIGHT_X_AXIS;
  protected int titleHeight = 0;
  protected double plotXmax, plotYmax;
  protected double plotXmin, plotYmin;
  protected volatile BufferedImage image;
  protected String prevPos = "";
  protected Hashtable<String, IntVector> locLookup;
  protected PlotPoint[] points; // make private when worked out
  protected GenericLine[] lines;
  protected GenericRectangle[] rectangles;
  protected GenericRectangle highlightRectangle;
  protected String xAxisLabel;
  protected String yAxisLabel;
  protected int xAxisLabelPad = DEFAULT_X_AXIS_LABEL_PAD;
  protected String title;
  protected boolean displayXAxis;
  protected boolean displayXLabel;
  protected boolean displayXAxisScale;
  protected boolean displayYAxis;
  protected boolean displayYLabel;
  protected boolean displayYAxisScale;
  protected boolean displayGrid;
  protected boolean displayTitle;
  protected boolean xAxisWholeNumbers;
  protected boolean yAxisWholeNumbers;
  protected int missingWidth;
  protected int nanWidth;
  private int axisFontSize = AXIS_FONT_SIZE;
  protected float forcePlotXmax, forcePlotYmax;
  protected float forcePlotXmin, forcePlotYmin;
  protected boolean createLookup;
  protected volatile boolean invertX;
  protected volatile boolean invertY;
  protected boolean makeSymmetric;
  protected String errorMessage;
  protected float mouseWheelMultiplier;
  protected boolean zoomable;
  protected boolean swapable; // 4/27/2012
  protected boolean invertable; // 4/27/2012
  protected boolean truncate;
  protected float[][] zoomSubsets;
  protected IntVector prox;
  protected int chartType;

  private Logger log = new Logger();

  private boolean inDrag;
  private volatile int startX, startY;
  // private int titleX, titleY;
  private int titleLocation;
  private final int plotPointSetSize;
  private final int totalNumPlotPointSets;
  private final int currentPlotPointSet;
  private final int lastIndexInPlotPointSet;
  private int currentIndexInPlotPointSet;
  private int lookupResolution;
  private boolean flow; // A control variable. If resizing is not yet done, don't start
                        // generatePoints() or drawAll();
  private volatile int imageStatus; // A control variable. If drawAll() is not yet done, don't start
                                    // paintComponent();
  private byte[] layersInBase;
  private byte[] extraLayersVisible;
  // private boolean pointsGeneratable;
  protected Timer waitingTimer; // A control variable to reduce the repaint() operations during
                                // component resizing;
  private String nullMessage;
  private final boolean randomTest;
  private int numberOfNaNSamples;
  private final boolean antiAlias = true;
  private volatile boolean beEfficient;
  private HashSet<Integer> pointsPlotted;

  public AbstractPanel() {
    canvasSectionMinimumX = 0;
    canvasSectionMaximumX = getWidth();
    canvasSectionMinimumY = 0;
    canvasSectionMaximumY = getWidth();
    displayXAxis = true;
    displayXLabel = true;
    displayXAxisScale = true;
    displayYAxis = true;
    displayYLabel = true;
    displayYAxisScale = true;
    displayGrid = false;
    createLookup = true;
    missingWidth = -1;
    nanWidth = -1;
    forcePlotXmax = forcePlotYmax = forcePlotXmin = forcePlotYmin = Float.NaN;
    mouseWheelMultiplier = DEFAULT_MOUSE_WHEEL_MULTIPLIER;
    zoomable = false;
    swapable = true;
    resetZoomProportions();
    plotPointSetSize = DEFAULT_PLOTPOINTSET_SIZE;
    points = new PlotPoint[plotPointSetSize];
    totalNumPlotPointSets = 1;
    currentPlotPointSet = 0;
    lastIndexInPlotPointSet = -1;
    currentIndexInPlotPointSet = -1;
    randomTest = false;
    setChartType(DEFAULT_TYPE);

    layersInBase = null;
    extraLayersVisible = null;

    image = null;
    locLookup = new Hashtable<>();
    imageStatus = IMAGE_NULL;
    flow = true;
    // pointsGeneratable = true;
    beEfficient = true;

    colorScheme = new Color[] {Color.BLACK, Color.GRAY};
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
    addComponentListener(this);
  }

  public void setLog(Logger log) {
    this.log = log;
  }

  public void createLookup(boolean value) {
    createLookup = value;
  }

  public void setImageStatus(int status) {
    if (DEBUGGING) {
      log.report("Set image status to " + status);
    }
    imageStatus = status;
  }

  protected boolean imageIsFinal() {
    return imageStatus == IMAGE_COMPLETE;
  }

  public void resetCurrentIndexInPlotPointSet() {
    currentIndexInPlotPointSet = -1;
  }

  public boolean morePlotPointsExist() {
    return currentPlotPointSet < totalNumPlotPointSets
           && currentIndexInPlotPointSet <= lastIndexInPlotPointSet;
  }

  public void setNullMessage(String str) {
    nullMessage = str;
  }

  public void setColorScheme(Color[] scheme) {
    colorScheme = scheme;
  }

  public void setZoomable(boolean zoomable, boolean truncate) {
    this.zoomable = zoomable;
    this.truncate = truncate;
  }

  public void beEfficient(boolean value) {
    beEfficient = value;
  }

  public void paintAgain() {
    image = null;
    setImageStatus(IMAGE_NULL);
    repaint();
  }

  public boolean isRandomTest() {
    return randomTest;
  }

  public String getNullMessage() {
    return nullMessage;
  }

  public boolean isFlow() {
    return flow;
  }

  public int getNumberOfNaNSamples() {
    return numberOfNaNSamples;
  }

  public void setNumberOfNaNSamples(int numberOfNaNSamples) {
    this.numberOfNaNSamples = numberOfNaNSamples;
  }

  public int getLookupResolution() {
    return lookupResolution;
  }

  public byte[] getLayersInBase() {
    return layersInBase;
  }

  public byte[] getExtraLayersVisible() {
    return extraLayersVisible;
  }

  public void setLayersInBase(byte[] layers) {
    layersInBase = layers;
  }

  public void setExtraLayersVisible(byte[] layers) {
    extraLayersVisible = layers;
  }

  @Override
  public void paintComponent(final Graphics g) {
    // this either doesn't affect anything or gets caught in an infinite loop while the lookup is
    // being created
    // while (imageStatus == IMAGE_STARTED) {
    // if (DEBUGGING) {
    // log.report("Additional call to paint before the first was completed; sleeping
    // 100ms");
    // }
    // try {
    // Thread.sleep(100);
    // } catch (InterruptedException ie) {
    // }
    // }

    if (image == null) {
      if (DEBUGGING) {
        log.report("createImage() being called from paintComponent()");
      }
      createImage(); // if you remove this, you get a blank screen and at least QQPlot ends up with
                     // a double title panel
    } else if (DEBUGGING) {
      log.report("Skipping image creation");
    }

    g.drawImage(image, 0, 0, AbstractPanel.this);
    if (extraLayersVisible != null && extraLayersVisible.length > 0) {
      drawAll(g, false);
    }
  }

  public void screenCapture(String filename) {
    boolean headless = GraphicsEnvironment.isHeadless();
    String imgDir = ext.parseDirectoryOfFile(filename);
    boolean mkdirs = new File(imgDir).mkdirs();
    File imgFile = new File(filename);
    if (mkdirs || Files.exists(ext.parseDirectoryOfFile(filename))) {
      if (image == null) {
        do {
          createImage();
          Thread.yield();
        } while (image == null || imageStatus != IMAGE_COMPLETE);
      }

      Runnable screenFunc = () -> {
        try {
          ImageIO.write(image, "png", imgFile);
        } catch (IOException ie) {
          if (headless) {
            log.reportError("Error while trying to save the plot");
          } else {
            JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
          }
        }
      };
      if (!SwingUtilities.isEventDispatchThread()) {
        try {
          SwingUtilities.invokeAndWait(screenFunc);
        } catch (InvocationTargetException | InterruptedException e) {
          log.reportException(e);
        }
      } else {
        screenFunc.run();
      }
    } else {
      if (headless) {
        log.reportError("Error creating directory in which to save the plot");
      } else {
        JOptionPane.showMessageDialog(null, "Error creating directory in which to save the plot");
      }
    }
  }

  public void setForcePlotXmin(float forcePlotXmin) {
    this.forcePlotXmin = forcePlotXmin;
  }

  public void setForcePlotXmax(float forcePlotXmax) {
    this.forcePlotXmax = forcePlotXmax;
  }

  public void setForcePlotYmin(float forcePlotYmin) {
    this.forcePlotYmin = forcePlotYmin;
  }

  public void setForcePlotYmax(float forcePlotYmax) {
    this.forcePlotYmax = forcePlotYmax;
  }

  public void setForceXAxisWholeNumbers(boolean whole) {
    xAxisWholeNumbers = whole;
  }

  public void setForceYAxisWholeNumbers(boolean whole) {
    yAxisWholeNumbers = whole;
  }

  public abstract void generatePoints();

  public abstract void highlightPoints();

  public abstract void assignAxisLabels();

  public int getAxisYWidth() {
    return axisYWidth;
  }

  public void setAxisYWidth(int axisYWidth) {
    this.axisYWidth = axisYWidth;
  }

  public int getAxisXHeight() {
    return axisXHeight;
  }

  public void setAxisXHeight(int axisXHeight) {
    this.axisXHeight = axisXHeight;
  }

  public void setAxisFontSize(int sz) {
    if (sz > 0) {
      axisFontSize = sz;
    }
  }

  public int getAxisFontSize() {
    return axisFontSize;
  }

  public void setChartType(int chartType) {
    this.chartType = chartType;
  }

  public void setXinversion(boolean b) {
    invertX = b;
  }

  public void setYinversion(boolean b) {
    invertY = b;
  }

  public void toggleXinversion() {
    invertX = !invertX;
  }

  public void toggleYinversion() {
    invertY = !invertY;
  }

  public void setSymmetricAxes(boolean b) {
    makeSymmetric = b;
  }

  public void drawAll(Graphics g, boolean base) {
    float minimumObservedRawX, maximumObservedRawX, minimumObservedRawY, maximumObservedRawY;
    double[] plotMinMaxStep; // needs to be double, else x <= plotXmax can be inexact and
                             // leave off
                             // the last tick mark
    String pos;
    int xLook, yLook;
    FontMetrics fontMetrics = null;
    Hashtable<String, Vector<PlotPoint>> layers;
    Vector<PlotPoint> layer;
    String trav;
    String[] keys;
    int step;
    long time;
    ProgressBarDialog prog;
    int rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel;

    setImageStatus(IMAGE_STARTED);

    long fullTime = System.currentTimeMillis();

    pointsPlotted = new HashSet<>();

    if (g instanceof Graphics2D) {
      ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                        antiAlias ? RenderingHints.VALUE_ANTIALIAS_ON
                                                  : RenderingHints.VALUE_ANTIALIAS_OFF);
      ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                                        antiAlias ? RenderingHints.VALUE_TEXT_ANTIALIAS_ON
                                                  : RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
    }

    // Set control variables; Generate data for the plot; set Lookup Resolution; Prepare AxisLabels.
    if (randomTest) {
      points = new PlotPoint[1000000];
      for (int i = 0; i < points.length; i++) {
        points[i] = new PlotPoint("", PointType.FILLED_CIRCLE, (float) Math.random(),
                                  (float) Math.random(), (byte) 5, (byte) 0, (byte) 0);
      }
    } else {
      generatePoints();
    }
    highlightPoints();

    if ((points == null || points.length == 0) && (rectangles == null || rectangles.length == 0)) {
      locLookup.clear();
      g.setColor(Color.WHITE);
      g.fillRect(0, 0, getWidth(), getHeight());
      if (nullMessage != null) {
        g.setColor(Color.BLACK);
        g.drawString(nullMessage,
                     getWidth() / 2 - g.getFontMetrics(g.getFont()).stringWidth(nullMessage) / 2,
                     getHeight() / 2);
      }
      setImageStatus(IMAGE_COMPLETE);
      return;
    }

    setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);
    assignAxisLabels();

    // Scan for rawX, rawY range of the data points
    minimumObservedRawX = Float.MAX_VALUE;
    maximumObservedRawX = Float.MIN_VALUE;
    minimumObservedRawY = Float.MAX_VALUE;
    maximumObservedRawY = Float.MIN_VALUE;
    for (int i = 0; i < points.length && flow; i++) {
      if (points[i] != null && Float.isFinite(points[i].getRawX())
          && Float.isFinite(points[i].getRawY())) {
        minimumObservedRawX = Maths.min(minimumObservedRawX, points[i].getRawX());
        maximumObservedRawX = Maths.max(maximumObservedRawX, points[i].getRawX());
        minimumObservedRawY = Maths.min(minimumObservedRawY, points[i].getRawY());
        maximumObservedRawY = Maths.max(maximumObservedRawY, points[i].getRawY());
      }
    }

    for (int i = 0; lines != null && i < lines.length && flow; i++) {
      if (lines[i] != null && lines[i].getLayer() == 0 && Float.isFinite(lines[i].getStartX())
          && Float.isFinite(lines[i].getStopX()) && Float.isFinite(lines[i].getStartY())
          && Float.isFinite(lines[i].getStopY())) {
        minimumObservedRawX = Maths.min(minimumObservedRawX, lines[i].getStartX());
        maximumObservedRawX = Maths.max(maximumObservedRawX, lines[i].getStartX());
        minimumObservedRawY = Maths.min(minimumObservedRawY, lines[i].getStartY());
        maximumObservedRawY = Maths.max(maximumObservedRawY, lines[i].getStartY());

        minimumObservedRawX = Maths.min(minimumObservedRawX, lines[i].getStopX());
        maximumObservedRawX = Maths.max(maximumObservedRawX, lines[i].getStopX());
        minimumObservedRawY = Maths.min(minimumObservedRawY, lines[i].getStopY());
        maximumObservedRawY = Maths.max(maximumObservedRawY, lines[i].getStopY());
      }
    }
    for (int i = 0; rectangles != null && i < rectangles.length && flow; i++) {
      if (rectangles[i] != null && rectangles[i].getLayer() == 0
          && Float.isFinite(rectangles[i].getStartXValue())
          && Float.isFinite(rectangles[i].getStopXValue())
          && Float.isFinite(rectangles[i].getStartYValue())
          && Float.isFinite(rectangles[i].getStopYValue())) {
        minimumObservedRawX = Maths.min(minimumObservedRawX, rectangles[i].getStartXValue());
        maximumObservedRawX = Maths.max(maximumObservedRawX, rectangles[i].getStartXValue());
        minimumObservedRawY = Maths.min(minimumObservedRawY, rectangles[i].getStartYValue());
        maximumObservedRawY = Maths.max(maximumObservedRawY, rectangles[i].getStartYValue());

        minimumObservedRawX = Maths.min(minimumObservedRawX, rectangles[i].getStopXValue());
        maximumObservedRawX = Maths.max(maximumObservedRawX, rectangles[i].getStopXValue());
        minimumObservedRawY = Maths.min(minimumObservedRawY, rectangles[i].getStopYValue());
        maximumObservedRawY = Maths.max(maximumObservedRawY, rectangles[i].getStopYValue());
      }
    }

    minimumObservedRawX = minimumObservedRawX == Float.MAX_VALUE ? 0 : minimumObservedRawX;
    maximumObservedRawX = maximumObservedRawX == Float.MIN_VALUE ? 1 : maximumObservedRawX;
    minimumObservedRawY = minimumObservedRawY == Float.MAX_VALUE ? 0 : minimumObservedRawY;
    maximumObservedRawY = maximumObservedRawY == Float.MIN_VALUE ? 1 : maximumObservedRawY;

    // otherwise step is off
    minimumObservedRawX = minimumObservedRawX > 0 ? 0 : minimumObservedRawX;
    minimumObservedRawY = minimumObservedRawY > 0 ? 0 : minimumObservedRawY;

    if (!Float.isNaN(forcePlotXmin)) {
      if (forcePlotXmin > minimumObservedRawX) {
        if (DEBUGGING) {
          log.reportError("WARNING - specified [minimum X boundary : " + forcePlotXmin
                          + "] is higher than the data point with the [lowest X value : "
                          + minimumObservedRawX + "]");
        }
      }
      minimumObservedRawX = forcePlotXmin;
    }
    if (Float.isNaN(forcePlotXmax)) {
      maximumObservedRawX = (maximumObservedRawX
                             + (maximumObservedRawX - minimumObservedRawX) * (float) 0.01);
    } else {
      if (forcePlotXmax < maximumObservedRawX) {
        if (DEBUGGING) {
          log.reportError("WARNING - specified [maximum X boundary : " + forcePlotXmax
                          + "] is lower than the data point with the [highest X value : "
                          + maximumObservedRawX + "]");
        }
      }
      maximumObservedRawX = forcePlotXmax;
    }
    if (!Float.isNaN(forcePlotYmin)) {
      if (forcePlotYmin > minimumObservedRawY) {
        if (DEBUGGING) {
          log.reportError("WARNING - specified [minimum Y boundary : " + forcePlotYmin
                          + "] is higher than the data point with the [lowest Y value : "
                          + minimumObservedRawY + "]");
        }
      }
      minimumObservedRawY = forcePlotYmin;
    }
    if (Float.isNaN(forcePlotYmax)) {
      maximumObservedRawY = (maximumObservedRawY
                             + (maximumObservedRawY - minimumObservedRawY) * (float) 0.01);
    } else {
      if (forcePlotYmax < maximumObservedRawY) {
        if (DEBUGGING) {
          log.reportError("WARNING - specified [maximum Y boundary : " + forcePlotYmax
                          + "] is lower than the data point with the [highest Y value : "
                          + maximumObservedRawY + "]");
        }
      }
      maximumObservedRawY = forcePlotYmax;
    }

    if (makeSymmetric) {
      maximumObservedRawX = Math.max(maximumObservedRawX, maximumObservedRawY);
      maximumObservedRawY = maximumObservedRawX;
      minimumObservedRawX = Math.min(minimumObservedRawX, minimumObservedRawY);
      minimumObservedRawY = minimumObservedRawX;
    }

    numberOfNaNSamples = 0;

    if (base) {
      if (DEBUGGING) {
        log.report("Drawing base image.");
      }
      g.fillRect(0, 0, getWidth(), getHeight());
      g.setFont(new Font("Arial", 0, getAxisFontSize()));

      titleHeight = calcTitleHeight(g, base, fontMetrics);

      drawTitle(g, base, fontMetrics);

      fontMetrics = g.getFontMetrics(g.getFont());
      missingWidth = fontMetrics.stringWidth("X");
      missingWidth = fontMetrics.stringWidth("X");

      // Calculate the plot area's range (X-axis, Y-axis)
      temporarilySetCanvasForXAxis();
      plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
      if (xAxisWholeNumbers) {
        if (plotMinMaxStep[2] < 1) {
          plotMinMaxStep[2] = 1;
        } else {
          plotMinMaxStep[2] = Math.round(plotMinMaxStep[2]);
        }
        if (plotMinMaxStep[3] >= (minimumObservedRawX - plotMinMaxStep[2])) {
          plotMinMaxStep[3] = plotMinMaxStep[3] - plotMinMaxStep[2];
        }
      }
      plotXmin = Float.isNaN(forcePlotXmin) ? plotMinMaxStep[0] : forcePlotXmin;
      plotXmax = Float.isNaN(forcePlotXmax) ? plotMinMaxStep[1] : forcePlotXmax;

      temporarilySetCanvasForYAxis();
      plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
      if (yAxisWholeNumbers) {
        if (plotMinMaxStep[2] < 1) {
          plotMinMaxStep[2] = 1;
        } else {
          plotMinMaxStep[2] = Math.round(plotMinMaxStep[2]);
        }
        if (plotMinMaxStep[3] >= (minimumObservedRawX - plotMinMaxStep[2])) {
          plotMinMaxStep[3] = plotMinMaxStep[3] - plotMinMaxStep[2];
        }
      }
      plotYmin = plotMinMaxStep[0];
      plotYmax = plotMinMaxStep[1];

    }

    // TODO outercoordinates
    resetCanvasToPlotArea();

    g.setClip(canvasSectionMinimumX, HEAD_BUFFER, canvasSectionMaximumX - canvasSectionMinimumX + 1,
              getHeight() - axisXHeight - 24);

    // Draw the lines
    for (int i = 0; lines != null && i < lines.length && flow; i++) {
      if ((base && (layersInBase == null || Bytes.indexOf(layersInBase, lines[i].getLayer()) >= 0))
          || (!base && Bytes.indexOf(extraLayersVisible, lines[i].getLayer()) >= 0)) {
        double x1raw = lines[i].getStartX();
        double y1raw = lines[i].getStartY();
        double x2raw = lines[i].getStopX();
        double y2raw = lines[i].getStopY();

        double x1final = x1raw, y1final = y1raw, x2final = x2raw, y2final = y2raw;

        int x1 = getXPixel(x1final);
        int y1 = getYPixel(y1final);
        int x2 = getXPixel(x2final);
        int y2 = getYPixel(y2final);

        Grafik.drawThickLine(g, x1, y1, x2, y2, lines[i].getThickness(),
                             colorScheme[lines[i].getColor()], lines[i].getDirection());
      }
    }

    // Draw the rectangles for clusterFilters
    for (int i = 0; rectangles != null && i < rectangles.length && flow; i++) {
      if ((base
           && (layersInBase == null || Bytes.indexOf(layersInBase, rectangles[i].getLayer()) >= 0))
          || (!base && Bytes.indexOf(extraLayersVisible, rectangles[i].getLayer()) >= 0)) {
        rectangleXPixel = Math.min(getXPixel(rectangles[i].getStartXValue()),
                                   getXPixel(rectangles[i].getStopXValue()));
        rectangleYPixel = Math.min(getYPixel(rectangles[i].getStartYValue()),
                                   getYPixel(rectangles[i].getStopYValue()));
        rectangleWidthPixel = Math.abs(getXPixel(rectangles[i].getStartXValue())
                                       - getXPixel(rectangles[i].getStopXValue()));
        rectangleHeightPixel = (Math.abs(getYPixel(rectangles[i].getStartYValue())
                                         - getYPixel(rectangles[i].getStopYValue())));
        g.setColor(colorScheme[rectangles[i].getColor()]);

        if (rectangles[i].getFill()) {
          if (rectangles[i].getColor() != rectangles[i].getFillColor()) {
            if (rectangles[i].getRoundedCorners()) {
              g.drawRoundRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel,
                              rectangleHeightPixel, 2, 2);
              g.setColor(colorScheme[rectangles[i].getFillColor()]);
              g.fillRoundRect(rectangleXPixel + rectangles[i].getThickness(),
                              rectangleYPixel + rectangles[i].getThickness(),
                              rectangleWidthPixel - (rectangles[i].getThickness() * 2) + 1,
                              rectangleHeightPixel - (rectangles[i].getThickness() * 2) + 1, 2, 2);
              g.setColor(colorScheme[rectangles[i].getColor()]);
            } else {
              drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel,
                            rectangleHeightPixel, rectangles[i].getThickness());
              g.setColor(colorScheme[rectangles[i].getFillColor()]);
              g.fillRect(rectangleXPixel + rectangles[i].getThickness(),
                         rectangleYPixel + rectangles[i].getThickness(),
                         rectangleWidthPixel - (rectangles[i].getThickness() * 2) + 1,
                         rectangleHeightPixel - (rectangles[i].getThickness() * 2) + 1);
              g.setColor(colorScheme[rectangles[i].getColor()]);
            }
          } else {
            if (rectangles[i].getRoundedCorners()) {
              g.fillRoundRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel,
                              rectangleHeightPixel, 2, 2);
            } else {
              g.fillRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel,
                         rectangleHeightPixel);
            }
          }
        } else {
          if (rectangles[i].getRoundedCorners()) {
            g.drawRoundRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel,
                            rectangleHeightPixel, 2, 2);
          } else {
            drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel,
                          rectangleHeightPixel, rectangles[i].getThickness());
          }
        }
      }
    }

    g.setClip(null);
    // Draw the rectangle outlined by dragging the mouse
    if (highlightRectangle != null) {
      rectangleXPixel = Math.min(getXPixel(highlightRectangle.getStartXValue()),
                                 getXPixel(highlightRectangle.getStopXValue()));
      rectangleYPixel = Math.min(getYPixel(highlightRectangle.getStartYValue()),
                                 getYPixel(highlightRectangle.getStopYValue()));
      rectangleWidthPixel = Math.abs(getXPixel(highlightRectangle.getStartXValue())
                                     - getXPixel(highlightRectangle.getStopXValue()));
      rectangleHeightPixel = (Math.abs(getYPixel(highlightRectangle.getStartYValue())
                                       - getYPixel(highlightRectangle.getStopYValue())));
      g.setColor(colorScheme[0]);
      drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel,
                    (byte) 1);
    }

    // Draw data points, also build the lookup matrix for nearby points.
    locLookup.clear();
    prog = null;
    time = new Date().getTime();
    step = Math.max((points.length) / 100, 1);
    layers = new Hashtable<>();

    if (chartType == HEAT_MAP_TYPE) {
      drawHeatMap(g);
    } else if (chartType == SCATTER_PLOT_TYPE) {
      boolean headless = GraphicsEnvironment.isHeadless();
      for (int i = 0; i < points.length && flow; i++) {
        if (!headless && base && i % step == 0) {
          if (new Date().getTime() - time > 1000) {
            if (prog == null) {
              prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(),
                                           getHeight(), 5000);
            }
            prog.setProgress(i);
          }
        }
        if (points[i] == null || points[i].getColor() == -1 || !points[i].isVisible()) {

        } else if (truncate && (points[i].getRawX() < plotXmin
                                || points[i].getRawX() - plotXmax > plotXmax / 1000.0
                                || points[i].getRawY() < plotYmin
                                || points[i].getRawY() > plotYmax)) {} else {
          trav = points[i].getLayer() + "";
          if (points[i].isHighlighted()
              || (base && (layersInBase == null
                           || Bytes.indexOf(layersInBase, points[i].getLayer()) >= 0))
              || (!base && Bytes.indexOf(extraLayersVisible, points[i].getLayer()) >= 0)) {
            if (trav.equals("0")) {
              if (points[i].getType() != PointType.NOT_A_NUMBER) {
                drawPoint(g, points[i]);
              } else if (base) {
                numberOfNaNSamples++;
              }
            } else {
              if (layers.containsKey(trav)) {
                layer = layers.get(trav);
              } else {
                layers.put(trav, layer = new Vector<>(points.length));
              }
              layer.add(points[i]);
            }
          }
          if (createLookup && points[i] != null) {
            xLook = (int) Math.floor(getXPixel(points[i].getRawX()) / lookupResolution);
            yLook = (int) Math.floor(getYPixel(points[i].getRawY()) / lookupResolution);
            for (int j = xLook - 1; j <= xLook + 1; j++) {
              for (int k = yLook - 1; k <= yLook + 1; k++) {
                pos = j + "x" + k;
                if (locLookup.containsKey(pos)) {
                  locLookup.get(pos).add(i);
                } else {
                  locLookup.put(pos, new IntVector(new int[] {i}));
                }
              }
            }
          }
        }
      }

      // Draw those points with layer > 0.
      keys = HashVec.getKeys(layers);
      for (int i = 0; i < keys.length && flow; i++) {
        layer = layers.get(keys[i]);
        for (int j = 0; j < layer.size(); j++) {
          if (layer.elementAt(j).getType() != PointType.NOT_A_NUMBER) {
            drawPoint(g, layer.elementAt(j));
          } else {
            numberOfNaNSamples++;
          }
        }
      }
    } else {
      log.reportError("Error - invalid chart type: " + chartType);
    }

    if (numberOfNaNSamples > 0) {
      g.drawString(PlotPoint.NAN_STR + " (n=" + numberOfNaNSamples + ")",
                   getXPixel(0) - nanWidth / 2, getYPixel(0) + 60 + points[0].getSize() / 2);
    }

    if (base && displayGrid) {
      for (double d = 0; d < 1.0; d += 0.1) {
        g.drawLine(getXPixel(d), getYPixel(0), getXPixel(d), getYPixel(canvasSectionMaximumY));
      }
      for (double d = -0.5; d < 0.5; d += 0.1) {
        g.drawLine(getXPixel(0), getYPixel(d), getXPixel(canvasSectionMaximumX), getYPixel(d));
      }
    }

    if (base) {
      temporarilySetCanvasForXAxis();
      plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
      drawXAxis(g, plotMinMaxStep, fontMetrics);

      temporarilySetCanvasForYAxis();
      plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
      drawYAxis(g, plotMinMaxStep);
      if (errorMessage != null) {
        g.drawString(errorMessage,
                     (getWidth() - getAxisYWidth()) / 2 - fontMetrics.stringWidth(errorMessage) / 2
                                   + getAxisYWidth(),
                     (getHeight() - HEAD_BUFFER - axisXHeight) / 2 - 20 + HEAD_BUFFER);
      }
      resetCanvasToPlotArea();
    }

    setImageStatus(IMAGE_COMPLETE);

    if (base && prog != null) {
      prog.close();
    }

    refreshOtherComponents();

    if (DEBUGGING) {
      log.report("Took " + ext.getTimeElapsed(fullTime) + " to draw "
                 + (createLookup ? "(and create lookup for) " : "") + points.length + " points");
    }
  }

  private void resetCanvasToPlotArea() {
    canvasSectionMinimumX = getAxisYWidth();// WIDTH_Y_AXIS;
    canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
    canvasSectionMinimumY = axisXHeight;// HEIGHT_X_AXIS;
    canvasSectionMaximumY = getHeight() - (HEAD_BUFFER + titleHeight);
  }

  /**
   * {@link #resetCanvasToPlotArea()} MUST be called after this method to ensure canvas is set
   * correctly for future responses to plot interaction
   */
  private void temporarilySetCanvasForXAxis() {
    canvasSectionMinimumX = getAxisYWidth();
    canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
    canvasSectionMinimumY = titleHeight;
    canvasSectionMaximumY = axisXHeight;
  }

  /**
   * {@link #resetCanvasToPlotArea()} MUST be called after this method to ensure canvas is set
   * correctly for future responses to plot interaction
   */
  private void temporarilySetCanvasForYAxis() {
    canvasSectionMinimumX = 0;
    canvasSectionMaximumX = getAxisYWidth();
    canvasSectionMinimumY = axisXHeight;
    canvasSectionMaximumY = getHeight() - (HEAD_BUFFER + titleHeight);
  }

  private void drawYAxis(Graphics g, double[] plotMinMaxStep) {
    int sigFigs;
    String str;
    BufferedImage yLabel;
    Graphics gfx;
    FontMetrics fontMetrics;
    sigFigs = ext.getNumSigFig(plotMinMaxStep[2]);
    float minSize = g.getFont().getSize2D();
    Font prevFont = g.getFont();
    for (double y = plotMinMaxStep[3]; y <= plotYmax; y += plotMinMaxStep[2]) {
      if (y >= plotYmin || !truncate) {
        str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
        if (str.length() == 5) {
          minSize = Math.min(minSize, prevFont.getSize2D() - 5);
        } else if (str.length() >= 6) {
          minSize = Math.min(minSize, prevFont.getSize2D() - 10);
        }
      }
    }
    Font minFont = prevFont.deriveFont(minSize);
    g.setFont(minFont);
    fontMetrics = g.getFontMetrics();
    if (displayYAxisScale) {
      for (double y = plotMinMaxStep[3]; y <= plotYmax; y += plotMinMaxStep[2]) {
        if (y >= plotYmin || !truncate) {
          Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(y),
                               canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
          str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
          g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 5 - fontMetrics.stringWidth(str),
                       getYPixel(y) + fontMetrics.getHeight() / 2);
        }
      }
    }
    g.setFont(prevFont);
    fontMetrics = g.getFontMetrics();
    if (displayYAxis) {
      Grafik.drawThickLine(g, canvasSectionMaximumX, getYPixel(plotYmin), canvasSectionMaximumX,
                           getYPixel(plotYmax) - (int) Math.ceil(TICK_THICKNESS / 2.0),
                           AXIS_THICKNESS, Color.BLACK);
    }
    if (displayYLabel) {
      int strWidth = fontMetrics.stringWidth(yAxisLabel);
      if (strWidth > 0) {
        yLabel = new BufferedImage(strWidth, 36, BufferedImage.TYPE_INT_RGB);
        gfx = yLabel.createGraphics();
        gfx.setFont(new Font("Arial", 0, getAxisFontSize()));
        gfx.setColor(Color.WHITE);
        gfx.fillRect(0, 0, getWidth(), getHeight());
        gfx.setColor(Color.BLACK);
        gfx.drawString(yAxisLabel, 0, yLabel.getHeight() - 6);

        int leftPad = 5; // TODO scale with window size if very small
        g.drawImage(Grafik.rotateImage(yLabel,
                                       true),
                    leftPad, (getHeight() + titleHeight - axisXHeight) / 2
                             - fontMetrics.stringWidth(yAxisLabel) / 2,
                    this);
      }
    }
  }

  protected void drawXAxis(Graphics g, double[] plotMinMaxStep, FontMetrics fontMetrics) {
    int sigFigs;
    String str;
    sigFigs = ext.getNumSigFig(plotMinMaxStep[2]);
    if (displayXAxisScale) {
      for (double x = plotMinMaxStep[3]; x <= plotXmax; x += plotMinMaxStep[2]) {
        int x1 = getXPixel(x);
        if ((x >= plotXmin && x1 >= canvasSectionMinimumX) || !truncate) {
          Grafik.drawThickLine(g, x1, getHeight() - canvasSectionMaximumY, x1,
                               getHeight() - (canvasSectionMaximumY - TICK_LENGTH), TICK_THICKNESS,
                               Color.BLACK);
          str = ext.formDeci(Math.abs(x) < DOUBLE_INACCURACY_HEDGE ? 0 : x, sigFigs, true);
          g.drawString(str, x1 - str.length() * 8,
                       getHeight() - (canvasSectionMaximumY - TICK_LENGTH - (2 * xAxisLabelPad)
                                      - 5));
        }
      }
    }
    if (displayXAxis) {
      Grafik.drawThickLine(g, canvasSectionMinimumX - (int) Math.ceil(AXIS_THICKNESS / 2.0),
                           getHeight() - canvasSectionMaximumY,
                           canvasSectionMaximumX + (int) Math.ceil(AXIS_THICKNESS / 2.0),
                           getHeight() - canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
    }
    if (xAxisLabel != null && !"".equals(xAxisLabel) && displayXLabel) {
      int yPad = getHeight() - xAxisLabelPad; // TODO scale with window size if very small
      g.drawString(xAxisLabel, (getWidth() - getAxisYWidth()) / 2
                               - fontMetrics.stringWidth(xAxisLabel) / 2 + getAxisYWidth(),
                   yPad);
    }
  }

  public void drawTitle(Graphics g, boolean base, FontMetrics fontMetrics) {
    if (base && displayTitle && title != null && !title.equals("")) {
      int titleX, titleY;
      int PAD = 5;
      int fontHeight = (fontMetrics == null ? 25 : fontMetrics.getHeight());
      int titleWidth = (fontMetrics == null ? title.length() * PAD
                                            : fontMetrics.stringWidth(title));

      switch (titleLocation) {
        default: // DEFAULT TO TOP
        case SwingConstants.NORTH:
          titleX = (getWidth() - getAxisYWidth()) / 2 + getAxisYWidth() - titleWidth / 2;
          titleY = PAD + fontHeight;
          break;
        case SwingConstants.NORTH_EAST:
          titleX = getWidth() - 2 * PAD - titleWidth;
          titleY = PAD + fontHeight;
          break;
        case SwingConstants.NORTH_WEST:
          titleX = getAxisYWidth() + 2 * PAD;
          titleY = PAD + fontHeight;
          break;
      }
      Color currColor = g.getColor();
      g.setColor(Color.black);
      g.drawString(title, titleX, titleY);
      g.setColor(currColor);
    }
  }

  public int calcTitleHeight(Graphics g, boolean base, FontMetrics fontMetrics) {
    return 0;
  }

  /**
   * Currently works with:<br />
   * <ul>
   * <li>SwingConstants.NORTH</li>
   * <li>SwingConstants.NORTH_EAST</li>
   * <li>SwingConstants.NORTH_WEST</li>
   * </ul>
   *
   * @param loc
   */
  public void setTitleLocation(int loc) {
    titleLocation = loc;
  }

  public void setTitle(String title) {
    this.title = title;
  }

  public void setDisplayTitle(boolean show) {
    displayTitle = show;
  }

  public boolean getDisplayTitle() {
    return displayTitle;
  }

  public String getTitle() {
    return title;
  }

  public void refreshOtherComponents() {}

  public void drawHeatMap(Graphics g) {
    int nRows, nColumns;
    int[][] gridIntensities;
    int[][][] gridColors;

    int width = 1, height = 1, radius = 2;

    nRows = getHeight();
    nColumns = getWidth();
    gridIntensities = getGridIntensityForHeapMapGrid(nRows, nColumns, width, height, radius);
    gridColors = getColorFromIntensityForHeapMapGrid(gridIntensities);
    for (int i = 0; i < nColumns; i++) {
      for (int j = 0; j < nRows; j++) {
        if (gridIntensities[i][j] != 0) {
          g.setColor(new Color(gridColors[i][j][0], gridColors[i][j][1], gridColors[i][j][2]));
          g.drawRect(i + canvasSectionMinimumX, j - canvasSectionMinimumY, 1, 1);
          // g.fillRect(i + canvasSectionMinimumX, j - canvasSectionMinimumY, 2, 2);
        }
      }
    }
  }

  public int[][] getGridIntensityForHeapMapGrid(int nRows, int nColumns, int cellWidth,
                                                int cellHeight, int neighbor) {
    int xPixel, yPixel;
    int[][] intensities;
    boolean zoomedIn;
    int[] origin;

    origin = new int[] {0, 0};

    zoomedIn = (Math.abs(getXValueFromXPixel(5) - getXValueFromXPixel(0)) < 0.002)
               || (Math.abs(getYValueFromYPixel(5) - getYValueFromYPixel(0)) < 0.002);
    intensities = new int[nColumns][nRows];
    for (PlotPoint point : points) {
      if (point != null) {
        xPixel = getXPixel(point.getRawX()) - canvasSectionMinimumX;
        yPixel = getYPixel(point.getRawY()) + canvasSectionMinimumY;

        for (int j = -neighbor; j <= neighbor; j++) {
          for (int k = -neighbor; k <= neighbor; k++) {
            if ((xPixel + j) >= 0 && (xPixel + j) < nColumns && (yPixel + k) >= 0
                && (yPixel + k) < nRows) { // && Distance.euclidean(new int[] {Math.abs(j),
                                           // Math.abs(k)}, origin) < neighbor*(neighbor-1)
              if (zoomedIn) {
                intensities[xPixel + j][yPixel + k]++;
              } else {
                intensities[xPixel
                            + j][yPixel
                                 + k] += neighbor * neighbor
                                         - Distance.euclidean(new int[] {Math.abs(j), Math.abs(k)},
                                                              origin);
              }
            }
          }
        }
      }
    }
    return intensities;
  }

  public int[][][] getColorFromIntensityForHeapMapGrid(int[][] intensities) {
    int[][][] color;
    int max;

    max = 0;
    for (int[] intensitie : intensities) {
      for (int element : intensitie) {
        if (max < element) {
          max = element;
        }
      }
    }

    color = new int[intensities.length][intensities[0].length][3];
    for (int i = 0; i < intensities.length; i++) {
      for (int j = 0; j < intensities[i].length; j++) {
        if (intensities[i][j] != 0) {
          color[i][j] = Grafik.getHeatmapColor((double) intensities[i][j] / (double) max);
        }
      }
    }

    return color;
  }

  @Override
  public void mouseClicked(MouseEvent e) {}

  @Override
  public void mouseEntered(MouseEvent e) {}

  @Override
  public void mouseExited(MouseEvent e) {}

  // public void mouseMoved(MouseEvent e) {}

  @Override
  public void mouseMoved(MouseEvent event) {
    Graphics g = getGraphics();
    IntVector indicesOfNearbyPoints;
    String pos;
    int x, y, dataPointIndex;
    byte size;

    if (imageIsFinal() && chartType != HEAT_MAP_TYPE) {
      x = event.getX();
      y = event.getY();
      resetCanvasToPlotArea();
      pos = (int) Math.floor(x / DEFAULT_LOOKUP_RESOLUTION) + "x"
            + (int) Math.floor(y / DEFAULT_LOOKUP_RESOLUTION);
      if (!pos.equals(prevPos)) {
        repaint();
      }

      indicesOfNearbyPoints = lookupNearbyPoints(x, y, pos);
      prox = new IntVector();

      size = SIZE * 2;
      g.setColor(Color.GRAY);
      for (int i = 0; indicesOfNearbyPoints != null && i < indicesOfNearbyPoints.size(); i++) {
        dataPointIndex = indicesOfNearbyPoints.elementAt(i);
        if (Distance.euclidean(new int[] {x, y},
                               new int[] {getXPixel(points[dataPointIndex].getRawX()),
                                          getYPixel(points[dataPointIndex].getRawY())}) < HIGHLIGHT_DISTANCE) {
          prox.add(dataPointIndex);
          // g.setColor(Color.YELLOW);
          // g.fillOval(getX(points[dataPointIndex].getRawX()) - size/2,
          // getY(points[dataPointIndex].getRawY()) - size/2, size, size);
          g.drawOval(getXPixel(points[dataPointIndex].getRawX()) - size / 2,
                     getYPixel(points[dataPointIndex].getRawY()) - size / 2, size, size);

        }
      }

      prevPos = pos;

      // TODO indicesOfNearbySamples = lookupNearbyPoints(x, y, pos);
    }
  }

  @Override
  public void mousePressed(MouseEvent e) {
    startX = e.getPoint().x;
    startY = e.getPoint().y;
    inDrag = true;
  }

  @Override
  public void mouseReleased(MouseEvent e) {
    inDrag = false;
  }

  @Override
  public void mouseDragged(MouseEvent e) {
    int curX, curY;
    double distance = -1;

    curX = e.getPoint().x;
    curY = e.getPoint().y;
    // if (invertX) {
    // curX = (2 * startX) - curX;
    // }
    // if (invertY) {
    // curY = (2 * startY) - curY;
    // }

    for (int i = 0; i < 2; i++) {
      if (i == 0) {
        distance = (startX - curX) * (zoomSubsets[0][1] - zoomSubsets[0][0])
                   / (getWidth() - WIDTH_BUFFER - getAxisYWidth()/* WIDTH_Y_AXIS */);
      } else {
        distance = (curY - startY) * (zoomSubsets[1][1] - zoomSubsets[1][0])
                   / (getHeight() - HEAD_BUFFER - axisXHeight/* HEIGHT_X_AXIS */);
      }

      if (distance < 0) {
        distance = Math.max(distance, -1 * zoomSubsets[i][0]);
      } else {
        distance = Math.min(distance, 1 - zoomSubsets[i][1]);
      }

      if ((zoomSubsets[i][0] <= 0 && distance < 0) || (zoomSubsets[i][1] >= 1 && distance > 0)) {

      } else {
        // if ((invertX && i == 0) || (invertY && i == 1)) {
        // zoomSubsets[i][0] -= distance;
        // zoomSubsets[i][1] -= distance;
        // } else {
        // zoomSubsets[i][0] += distance;
        // zoomSubsets[i][1] += distance;
        // }
        zoomSubsets[i][0] += distance;
        zoomSubsets[i][1] += distance;
      }
    }
    if (inDrag) {
      paintAgain();
      startX = curX;
      startY = curY;
    }
  }

  @Override
  public void mouseWheelMoved(MouseWheelEvent e) {
    if (zoomable) {
      if (e.getWheelRotation() < 0
          && zoomSubsets[0][1] - zoomSubsets[0][0] < MINIMUM_ZOOM_PROPORTION_WINDOW) {
        return;
      }
      zoomProportionally(e.getWheelRotation() > 0, e.getPoint(), false);
    }
  }

  @Override
  public void componentHidden(ComponentEvent e) {}

  @Override
  public void componentMoved(ComponentEvent e) {}

  @Override
  public void componentResized(ComponentEvent e) {
    int w = (int) e.getComponent().getSize().getWidth();
    int h = (int) e.getComponent().getSize().getHeight();
    if (trackedW == -1 || trackedH == -1) {
      setTrackedSize(w, h);
      return;
    }
    if (w != trackedW || h != trackedH) {
      setTrackedSize(w, h);
    }

    if (waitingTimer == null) {
      /* Start waiting for DELAY to elapse. */
      waitingTimer = new Timer(DELAY, this);
      waitingTimer.start();
    } else {
      /* Event came too soon, swallow it by resetting the timer.. */
      waitingTimer.restart();
    }
  }

  private volatile int trackedW = -1;
  private volatile int trackedH = -1;

  public void setTrackedSize(int w, int h) {
    trackedW = w;
    trackedH = h;
  }

  @Override
  public void componentShown(ComponentEvent e) {}

  @Override
  public void actionPerformed(ActionEvent ae) {
    /* Timer finished? */
    if (ae.getSource() == waitingTimer) {
      /* Stop timer */
      waitingTimer.stop();
      waitingTimer = null;
      /* Resize */
      setFlow(true);
      if (DEBUGGING) {
        log.reportError("Action performed in AbstractPanel");
      }
      createImage();
      repaint();
    }

  }

  public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
    float x, y, multiplier, dist;
    float[][] proportions;
    int width, height;
    boolean changed;

    proportions = new float[2][2];

    width = getWidth() - /* WIDTH_Y_AXIS */getAxisYWidth() - WIDTH_BUFFER;
    x = (float) p.getX() - getAxisYWidth();// WIDTH_Y_AXIS;
    proportions[0][0] = x / width;
    proportions[0][1] = (width - x) / width;

    height = getHeight() - /* HEIGHT_X_AXIS */axisXHeight - HEAD_BUFFER;
    y = (float) p.getY() - HEAD_BUFFER; // could be HEAD_BUFFER
    proportions[1][0] = (height - y) / height; // reversed because of top down
    proportions[1][1] = y / height;

    multiplier = mouseWheelMultiplier / (outNotIn ? 1 : -2);

    if (!outNotIn && center) {
      for (int i = 0; i < proportions.length; i++) {
        for (int j = 0; j < proportions[i].length; j++) {
          proportions[i][j] = 0.25f - proportions[i][j];
        }
      }
      multiplier = mouseWheelMultiplier;
    }

    for (int i = 0; i < proportions.length; i++) {
      dist = zoomSubsets[i][1] - zoomSubsets[i][0];
      for (int j = 0; j < proportions[i].length; j++) {
        zoomSubsets[i][j] = zoomSubsets[i][j]
                            + (j == 0 ? -1 : 1) * (proportions[i][j] * multiplier * dist);
      }
    }

    if (zoomSubsets[0][1] - zoomSubsets[0][0] > 1) {
      resetZoomProportions();
    }

    // necessary to prevent distortion in X versus Y when zooming in/out near an edge
    do {
      changed = false;
      for (int i = 0; i < 2; i++) {
        if (zoomSubsets[i][0] < -0.0001) {
          zoomSubsets[i][1] -= zoomSubsets[i][0];
          zoomSubsets[i][0] = 0;
          changed = true;
        }
        if (zoomSubsets[i][1] > 1.0001) {
          zoomSubsets[i][0] -= zoomSubsets[i][1] - 1;
          zoomSubsets[i][1] = 1;
          changed = true;
        }
      }
    } while (changed);

    paintAgain();
  }

  public void resetZoomProportions() {
    zoomSubsets = new float[2][2];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        zoomSubsets[j][i] = i;
      }
    }
  }

  public int getNumPointsPlottedEfficiently() {
    if (beEfficient) {
      return pointsPlotted.size();// +"";
    } else {
      return -1;// "[not set to be efficient]";
    }
  }

  private int getEfficientPointCode(int x, int y, int sz, int clr) {
    // Standard hashCode implementation fails here (fun behavior though!)
    // final int prime = 31;
    // int result = 17;
    // result = prime * result + x;
    // result = prime * result + y;
    // result = prime * result + sz;
    // result = prime * result + clr;
    // return result;

    // VERY slow, possibly not performing properly:
    // return (new int[]{x, y, sz, clr}).hashCode();

    // Szudzik: (unfortunately only applicable for two integers - should we multiplex (i.e. x/y &
    // sz/clr?)
    // https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
    int coord = szudzikCode(x, y);
    int meta = szudzikCode(sz, clr);
    int ret = szudzikCode(coord, meta);
    return ret;
  }

  private int szudzikCode(int x, int y) {
    int A = x >= 0 ? 2 * x : -2 * x - 1;
    int B = y >= 0 ? 2 * y : -2 * y - 1;
    int C = A >= B ? A * A + A + B : A + B * B;
    return x < 0 && y < 0 || x >= 0 && y >= 0 ? C : -C - 1;
  }

  public void drawPoint(Graphics g, PlotPoint point) {
    int x, y, size, color;

    x = getXPixel(point.getRawX());
    y = getYPixel(point.getRawY());
    size = point.getSize();
    color = point.getColor();

    g.setColor(colorScheme[color]);

    if (size == 0) {
      return;
    }
    if (beEfficient) {
      int code = getEfficientPointCode(x, y, size, color);
      if (pointsPlotted.contains(code)) {
        return;
      } else {
        pointsPlotted.add(code);
      }
      // Files.appendStringToFile("listOfPoints.out", x+":"+y+":"+size+":"+color);
    }

    switch (point.getType()) {
      case FILLED_CIRCLE:
        g.fillOval(x - size / 2, y - size / 2, size, size);
        break;
      case OPEN_CIRCLE:
        g.drawOval(x - size / 2, y - size / 2, size, size);
        break;
      case FILLED_SQUARE:
        g.fillPolygon(new int[] {x - size / 2, x - size / 2, x + size / 2, x + size / 2},
                      new int[] {y - size / 2, y + size / 2, y + size / 2, y - size / 2}, 4);
        break;
      case OPEN_SQUARE:
        g.drawPolygon(new int[] {x - size / 2, x - size / 2, x + size / 2, x + size / 2},
                      new int[] {y - size / 2, y + size / 2, y + size / 2, y - size / 2}, 4);
        break;
      case FILLED_TRIANGLE:
        Grafik.drawTriangle(g, x, y, size, true);
        // g.drawPolygon(new int[] {x-size/2, x, +size/2},
        // new int[] {y+size/2, y-size/2, y+size/2},
        // 3);
        break;
      case MISSING:
        // g.drawString(PlotPoint.MISSING_STR, getX(point.getRawX())-missingWidth/2,
        // getY(point.getRawY())+size/2);
        if (PlotPoint.MISSING_STR.equals("X") || PlotPoint.MISSING_STR.equals("x")) {
          g.drawLine(x - size / 4, y - size / 4, x + size / 4, y + size / 4);
          g.drawLine(x - size / 4, y + size / 4, x + size / 4, y - size / 4);
        } else {
          g.setFont(new Font("Arial", 0, size));
          g.drawString(PlotPoint.MISSING_STR, -size / 2, y + size / 2);
          g.setFont(new Font("Arial", 0, getAxisFontSize()));
        }
        break;
      case NOT_A_NUMBER:
        // g.drawString(PlotPoint.NAN_STR, getX(point.getRawX())-nanWidth/2,
        // getY(point.getRawY())-30+size/2);
        break;
      default:
        log.reportError("Error - invalid PlotPoint type");
    }
  }

  public double calcStepStep(double range) {
    String[] line;

    try {
      line = new DecimalFormat("0.0E0").format(range).split("E");
      return Math.pow(10, Integer.parseInt(line[1]) - 1)
             * (Double.parseDouble(line[0]) > 2.0 ? 5 : 1);
    } catch (Exception e) {
      log.reportError("Error - could not parse stepStep from range '" + range + "'");
      return Double.NaN;
    }
  }

  public double[] getPlotMinMaxStep(double min, double max, Graphics g, boolean xAxis) {
    double range, plotStep, stepStep, plotMin, plotMax;
    double zoomMin, zoomMax, dist, tempD;
    int numHashes, wid, hgt, canvasRange;
    FontMetrics fontMetrics;
    int sf, temp;

    range = max - min;

    plotStep = stepStep = calcStepStep(range);
    sf = ext.getNumSigFig(stepStep);

    fontMetrics = g.getFontMetrics(g.getFont());
    if (xAxis) {
      wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)),
                     fontMetrics.stringWidth(ext.formDeci(max, sf)));
      numHashes = 12;
      canvasRange = canvasSectionMaximumX - canvasSectionMinimumX;
      temp = (wid + 30) * numHashes;
      while (temp > canvasRange) {
        numHashes -= 1;
        temp = (wid + 30) * numHashes;
      }
    } else {
      hgt = fontMetrics.getHeight();
      numHashes = 10;
      canvasRange = canvasSectionMaximumY - canvasSectionMinimumY;
      temp = (hgt + 10) * numHashes;
      while (temp > canvasRange) {
        numHashes -= 1;
        temp = (hgt + 10) * numHashes;
      }
    }

    numHashes = Math.max(numHashes, 1);
    while (range / plotStep > numHashes) {
      plotStep += stepStep;
    }

    plotMin = plotMax = 0;
    while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
      plotMax += plotStep;
    }
    if (min > plotMin) {
      tempD = min - (plotMin + 2 * plotStep);
      while (tempD > DOUBLE_INACCURACY_HEDGE) {
        plotMin += plotStep;
        tempD = min - (plotMin + 2 * plotStep);
      }
    } else {
      tempD = min - plotMin;
      while (tempD < -DOUBLE_INACCURACY_HEDGE) {
        plotMin -= plotStep;
        tempD = min - plotMin;
      }
    }

    // plotMin = plotMax = 0;
    // while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
    // plotMax += plotStep;
    // }
    // while (min - plotMin < -1*DOUBLE_INACCURACY_HEDGE) { // double check this, untested
    // plotMin -= plotStep;
    // }
    // // if (min >= 0 && plotMin < 0) {
    // // plotMin = 0;
    // // }
    //
    // // log.report(Float.parseFloat(ext.formDeci(plotMin,
    // sf))+"\t"+Float.parseFloat(ext.formDeci(plotMax,
    // sf))+"\t"+Float.parseFloat(ext.formDeci(plotStep, sf)));
    //
    if (zoomable) {
      dist = plotMax - plotMin;
      zoomMin = plotMin + zoomSubsets[xAxis ? 0 : 1][0] * dist;
      zoomMax = plotMax - (1 - zoomSubsets[xAxis ? 0 : 1][1]) * dist;

      range = zoomMax - zoomMin;
      plotStep = stepStep = calcStepStep(range);
      sf = ext.getNumSigFig(stepStep);

      if (xAxis) {
        fontMetrics = g.getFontMetrics(g.getFont());
        wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)),
                       fontMetrics.stringWidth(ext.formDeci(max, sf)));
        numHashes = 12;
        while ((wid + 30) * numHashes > canvasSectionMaximumX - canvasSectionMinimumX) {
          numHashes -= 2;
        }
      } else {
        numHashes = 10;
      }

      numHashes = Math.max(numHashes, 1);
      while (range / plotStep > numHashes) {
        plotStep += stepStep;
      }
      return new double[] {zoomMin, zoomMax, Double.parseDouble(ext.formDeci(plotStep, sf)),
                           Double.parseDouble(ext.formDeci(plotMin, sf))};
    } else {
      return new double[] {Double.parseDouble(ext.formDeci(plotMin, sf)),
                           Double.parseDouble(ext.formDeci(plotMax, sf)),
                           Double.parseDouble(ext.formDeci(plotStep, sf)),
                           Double.parseDouble(ext.formDeci(plotMin, sf))};
    }
  }

  public int getXPixel(double x) {
    if (invertX) {
      return (int) ((plotXmax - x) / (plotXmax - plotXmin)
                    * (canvasSectionMaximumX - canvasSectionMinimumX))
             + canvasSectionMinimumX;
    } else {
      return (int) ((x - plotXmin) / (plotXmax - plotXmin)
                    * (canvasSectionMaximumX - canvasSectionMinimumX))
             + canvasSectionMinimumX;
    }
  }

  public int getYPixel(double y) {
    if (invertY) {
      return getHeight()
             - (int) ((plotYmax - y) / (plotYmax - plotYmin)
                      * (canvasSectionMaximumY - canvasSectionMinimumY) + canvasSectionMinimumY);
    } else {
      return getHeight()
             - (int) ((y - plotYmin) / (plotYmax - plotYmin)
                      * (canvasSectionMaximumY - canvasSectionMinimumY) + canvasSectionMinimumY);
    }
  }

  /**
   * Converts mouse location int X,Y into data points' value double rawX,rawY. The control variable
   * invertX is assigned elsewhere in the class.
   *
   * @param mouseX the mouse location X
   * @return the rawX value of the corresponding data point.
   */
  public double getXValueFromXPixel(int mouseX) {
    if (invertX) {
      return plotXmax
             - ((double) (mouseX - canvasSectionMinimumX)
                / (double) (canvasSectionMaximumX - canvasSectionMinimumX) * (plotXmax - plotXmin));
    } else {
      return plotXmin
             + ((double) (mouseX - canvasSectionMinimumX)
                / (double) (canvasSectionMaximumX - canvasSectionMinimumX) * (plotXmax - plotXmin));
    }
  }

  /**
   * Converts mouse location int X,Y into data points' value double rawX,rawY The control variable
   * invertY is assigned elsewhere in the class.
   *
   * @param mouseY the mouse location Y
   * @return the rawY value of the corresponding data point.
   */
  public double getYValueFromYPixel(int mouseY) {
    if (invertY) {
      return plotYmax
             + ((double) (mouseY + canvasSectionMinimumY - getHeight())
                / (double) (canvasSectionMaximumY - canvasSectionMinimumY) * (plotYmax - plotYmin));
    } else {
      return plotYmin
             - ((double) (mouseY + canvasSectionMinimumY - getHeight())
                / (double) (canvasSectionMaximumY - canvasSectionMinimumY) * (plotYmax - plotYmin));
    }
  }

  public double getXIntercept(double x1, double y1, double x2, double y2, double yNew) {
    // return (x1 * (y1 - 2*y2 + yNew) + x2 * (y2 - yNew)) / (y1 - y2);

    return ((x2 - x1) * (yNew - y2 + x2 * ((y2 - y1) / (x2 - x1)))) / (y2 - y1);

  }

  public double getYIntercept(double x1, double y1, double x2, double y2, double xNew) {
    double s = ((y2 - y1) / (x2 - x1));
    double b = y2 - s * x2;
    return s * xNew + b;
  }

  // /**
  // * Screens the data points to find out those that fall into the range of (double rawXmin, double
  // rawXmax, double rawYmin, double rawYmax)
  // * @param rawXmin
  // * @param rawXmax
  // * @param rawYmin
  // * @param rawYmax
  // */
  // public void highlightPoints(double rawXmin, double rawYmin, double rawXmax, double rawYmax) {
  // if (points.length==100){
  // for (int i=0; i<points.length; i++) {
  // if (points[i].getRawX()>=rawXmin
  // && points[i].getRawY()>=rawYmin
  // && points[i].getRawX()<=rawXmax
  // && points[i].getRawY()<=rawYmax) {
  // points[i].setHighlighted(true);
  // log.report("Highlighting: "+points[i].getRawX()+","+points[i].getRawY());
  // } else {
  // points[i].setHighlighted(false);
  // }
  //
  // }
  // }
  // }
  //
  /**
   * Highlights those points that need to be highlighted
   *
   * @param array
   */
  public void highlightPoints(boolean[] array) {
    if (points.length != array.length) {
      if (DEBUGGING) {
        log.reportError("Error - mismatched array size when highlighting");
      }
    } else {
      for (int i = 0; i < points.length; i++) {
        if (points[i] != null && array[i]) {
          points[i].setHighlighted(true);
        } else if (points[i] != null) {
          points[i].setHighlighted(false);
        }

      }
    }
  }

  public void createImage() {
    setImageStatus(IMAGE_STARTED);
    image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
    flow = true;
    if (DEBUGGING) {
      log.report("Drawing base image");
    }
    drawAll(image.createGraphics(), true);
  }

  public void setLookupResolution(int lookupResolution) {
    this.lookupResolution = lookupResolution;
  }

  public IntVector lookupNearbyPoints(int x, int y, String pos) {
    IntVector iv = locLookup.get(pos);
    IntVector indicesOfDataPoints = new IntVector();

    for (int i = 0; iv != null && i < iv.size(); i++) {
      if (Distance.euclidean(new int[] {x, y},
                             new int[] {getXPixel(points[iv.elementAt(i)].getRawX()),
                                        getYPixel(points[iv.elementAt(i)].getRawY())}) < HIGHLIGHT_DISTANCE) {
        indicesOfDataPoints.add(iv.elementAt(i));
      }
    }

    return indicesOfDataPoints;
  }

  public byte lookupNearbyRectangles(int xPixel, int yPixel) {
    int rectangleStartXPixel, rectangleStartYPixel, rectangleStopXPixel, rectangleStopYPixel;
    byte indicesRectangle;
    int distance;
    int minDistance;
    int[] rectangleEdges;

    indicesRectangle = -1;
    minDistance = Integer.MAX_VALUE;
    for (byte i = 0; rectangles != null && i < rectangles.length; i++) {
      rectangleStartXPixel = getXPixel(rectangles[i].getStartXValue());
      rectangleStartYPixel = getYPixel(rectangles[i].getStartYValue());
      rectangleStopXPixel = getXPixel(rectangles[i].getStopXValue());
      rectangleStopYPixel = getYPixel(rectangles[i].getStopYValue());

      rectangleEdges = new int[] {Math.min(rectangleStartXPixel, rectangleStopXPixel),
                                  Math.max(rectangleStartXPixel, rectangleStopXPixel),
                                  Math.min(rectangleStartYPixel, rectangleStopYPixel),
                                  Math.max(rectangleStartYPixel, rectangleStopYPixel)};
      for (int j = 0; j < 2; j++) {
        distance = Math.abs(xPixel - rectangleEdges[j]);
        if (distance < HIGHLIGHT_DISTANCE && distance < minDistance
            && yPixel > rectangleEdges[2] - HIGHLIGHT_DISTANCE
            && yPixel < rectangleEdges[3] + HIGHLIGHT_DISTANCE) {
          minDistance = distance;
          indicesRectangle = i;
        }
      }
      for (int j = 2; j < 4; j++) {
        distance = Math.abs(yPixel - rectangleEdges[j]);
        if (distance < HIGHLIGHT_DISTANCE && distance < minDistance
            && xPixel > rectangleEdges[0] - HIGHLIGHT_DISTANCE
            && xPixel < rectangleEdges[1] + HIGHLIGHT_DISTANCE) {
          minDistance = distance;
          indicesRectangle = i;
        }
      }
    }

    return indicesRectangle;
  }

  public byte lookupResidingRectangles(int xPixel, int yPixel) {
    int rectangleStartXPixel, rectangleStartYPixel, rectangleStopXPixel, rectangleStopYPixel;
    byte indicesRectangle;
    int distance;
    int minDistance;
    int[] rectangleEdges;
    int[] rectangleSearchRange;

    indicesRectangle = -1;
    minDistance = Integer.MAX_VALUE;
    for (byte i = 0; rectangles != null && i < rectangles.length; i++) {
      rectangleStartXPixel = getXPixel(rectangles[i].getStartXValue());
      rectangleStartYPixel = getYPixel(rectangles[i].getStartYValue());
      rectangleStopXPixel = getXPixel(rectangles[i].getStopXValue());
      rectangleStopYPixel = getYPixel(rectangles[i].getStopYValue());

      rectangleEdges = new int[] {Math.min(rectangleStartXPixel, rectangleStopXPixel),
                                  Math.max(rectangleStartXPixel, rectangleStopXPixel),
                                  Math.min(rectangleStartYPixel, rectangleStopYPixel),
                                  Math.max(rectangleStartYPixel, rectangleStopYPixel)};
      rectangleSearchRange = new int[] {(int) (rectangleEdges[0] - HIGHLIGHT_DISTANCE),
                                        (int) (rectangleEdges[1] + HIGHLIGHT_DISTANCE),
                                        (int) (rectangleEdges[2] - HIGHLIGHT_DISTANCE),
                                        (int) (rectangleEdges[3] + HIGHLIGHT_DISTANCE)};
      if (xPixel >= rectangleSearchRange[0] && xPixel <= rectangleSearchRange[1]
          && yPixel >= rectangleSearchRange[2] && yPixel <= rectangleSearchRange[3]) {
        for (int j = 0; j < 2; j++) {
          distance = Math.abs(xPixel - rectangleEdges[j]);
          if (minDistance > distance) {
            minDistance = distance;
            indicesRectangle = i;
          }
        }
        for (int j = 2; j < 4; j++) {
          distance = Math.abs(yPixel - rectangleEdges[j]);
          if (minDistance > distance) {
            minDistance = distance;
            indicesRectangle = i;
          }
        }
      }
    }

    return indicesRectangle;
  }

  public void setFlow(boolean flow) {
    this.flow = flow;
  }

  public boolean getFlow() {
    return flow;
  }

  public void setPointsGeneratable(boolean pointsGeneratable) {
    // this.pointsGeneratable = pointsGeneratable;
  }

  public boolean isPointsGeneratable() {
    return true /* pointsGeneratable */;
  }

  public void setSwapable(boolean swapable) {
    this.swapable = swapable;
  }

  public boolean isSwapable() {
    return swapable;
  }

  public void drawRectThick(Graphics g, int x, int y, int width, int height, byte thickness) {
    for (byte i = 0; i < thickness; i++) {
      g.drawRect(x + i, y + i, width - 2 * i, height - 2 * i);
    }
  }
}
