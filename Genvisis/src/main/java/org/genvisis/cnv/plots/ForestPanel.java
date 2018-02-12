package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.RenderingHints;
import java.text.DecimalFormat;
import java.util.ArrayList;
import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Grafik;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF.Colors.BLUES;
import org.genvisis.common.PSF.Colors.GREENS;
import org.genvisis.common.PSF.Colors.REDS;
import org.genvisis.common.PSF.Colors.VIOLETS;
import org.genvisis.common.ext;
import org.genvisis.stats.Maths;
import com.google.common.primitives.Bytes;

/**
 * Forest Panel
 */
public class ForestPanel extends AbstractPanel {

  // FORESTPLOT
  static final long serialVersionUID = 4023579107L;

  public static final Color[] DEFAULT_COLORS = {BLUES.MIDNIGHT_EXPRESS, // dark dark
                                                BLUES.PERSIAN_BLUE, // dark blue
                                                REDS.VENETIAN_RED, // deep red
                                                VIOLETS.BLUE_VIOLET, // deep purple
                                                GREENS.GREEN, // dark green
                                                BLUES.DODGER_BLUE, // light blue
                                                BLUES.SLATE_BLUE, // light purple
                                                GREENS.GREEN_YELLOW, // light green
                                                VIOLETS.ORCHID, // pink
                                                BLUES.NAVY, BLUES.CORNFLOWER_BLUE,
                                                BLUES.DARK_SLATE_BLUE, BLUES.SLATE_BLUE,
                                                BLUES.MEDIUM_SLATE_BLUE, BLUES.LIGHT_SLATE_BLUE,
                                                BLUES.MEDIUM_BLUE, BLUES.ROYAL_BLUE, Color.BLUE,
                                                BLUES.DODGER_BLUE, BLUES.DEEP_SKY_BLUE,
                                                BLUES.LIGHT_SKY_BLUE, BLUES.LIGHT_SKY_BLUE,
                                                BLUES.STEEL_BLUE, BLUES.LIGHT_STEEL_BLUE,
                                                BLUES.LIGHT_BLUE, BLUES.POWDER_BLUE,
                                                BLUES.PALE_TURQUOISE, BLUES.DARK_TURQUOISE,
                                                BLUES.MEDIUM_TURQUOISE, BLUES.TURQUOISE, BLUES.AQUA,
                                                BLUES.LIGHT_CYAN,

  };

  public static final String META_LABEL = "Overall";
  public static final Color META_COLOR = Color.BLACK;

  private final ForestPlot forestPlot;
  private final Logger log;
  private boolean rectangleGeneratable;
  private final DecimalFormat precision2Decimal;
  private final boolean antiAlias = true;
  boolean oddsDisplay = false;

  public ForestPanel(ForestPlot forestPlot, Logger log) {
    super();
    rectangleGeneratable = true;

    this.forestPlot = forestPlot;
    this.log = log;
    setColorScheme(DEFAULT_COLORS);
    precision2Decimal = new DecimalFormat("#0.00");

    setZoomable(false, true);
    // setZoomable(true, true);
    setColorScheme(DEFAULT_COLORS);
  }

  public boolean isRectangleGeneratable() {
    return rectangleGeneratable;
  }

  public void setRectangleGeneratable(boolean rectangleGeneratable) {
    this.rectangleGeneratable = rectangleGeneratable;
  }

  @Override
  public void generatePoints() {
    if (forestPlot.getCurrentMetaStudy() == null) {
      ForestInput input;
      String errorMessage;

      if (forestPlot.getDataIndices().size() == 0) {
        errorMessage = "No data file selected, or no data found in input file.  Please select a data file.";
      } else {
        input = forestPlot.getDataIndices().get(forestPlot.getCurrentDataIndex());
        errorMessage = "Cannot generate points for marker " + input.marker
                       + " because the data did not load; check to see if file \"" + input.file
                       + "\" actually exists and if the beta/se columns are named as expected (e.g., expecting beta.Study1 and not Study1.beta; overall results need to be exactly beta/se or effect/stderr)";
      }
      setNullMessage(errorMessage);
      if (log != null) {
        log.reportError("Error - " + errorMessage);
      }

      points = new PlotPoint[0];
      setPointsGeneratable(false);
      setRectangleGeneratable(false);

      return;
    }
    @SuppressWarnings("unchecked")
    ArrayList<StudyData> currentData = (ArrayList<StudyData>) forestPlot.getCurrentMetaStudy()
                                                                        .getStudies().clone();
    // PlotPoint[] tempPoints = new PlotPoint[currentData.size()];
    ArrayList<PlotPoint> tempPoints = new ArrayList<PlotPoint>();
    ArrayList<GenericLine> linesData = new ArrayList<GenericLine>();

    float xAxisValue, yAxisValue;
    lines = new GenericLine[0];

    for (int i = 0; i < currentData.size(); i++) {
      StudyData currStudy = currentData.get(i);
      // if(currentData.get(i).getBeta() != 0.0 && currentData.get(i).getStderr() != 0.0){
      if (currStudy instanceof StudyBreak) {
        yAxisValue = (float) i + 1;
        PlotPoint leftEnd = new PlotPoint("", PointType.FILLED_CIRCLE, oddsDisplay ? 1 : 0,
                                          yAxisValue, (byte) 5, (byte) 0, (byte) 0);
        yAxisValue = (float) i + 1;
        PlotPoint rightEnd = new PlotPoint("", PointType.FILLED_CIRCLE, oddsDisplay ? 1 : 0,
                                           yAxisValue, (byte) 5, (byte) 0, (byte) 0);

        linesData.add(new GenericLine(leftEnd, rightEnd, (byte) 1, (byte) 0, (byte) 0, false));

        yAxisValue = (float) i + 1;

        tempPoints.add(new PlotPoint(" | ", PointType.FILLED_CIRCLE, oddsDisplay ? 1 : 0,
                                     yAxisValue, (byte) 3, (byte) 0, (byte) 0));
        tempPoints.get(tempPoints.size() - 1).setVisible(false);
        // tempPoints[i] = new PlotPoint(" | ", (byte) 0, oddsDisplay ? 1 : 0, yAxisValue, (byte) 3,
        // (byte) 0, (byte) 0);
        // tempPoints[i].setVisible(false);
      } else {
        xAxisValue = currStudy.getConfInterval(oddsDisplay)[0];
        yAxisValue = (float) i + 1;
        PlotPoint leftEnd = new PlotPoint(currStudy.getDisplayLabel(), currStudy.getShape(),
                                          xAxisValue, yAxisValue, (byte) 5, (byte) 0, (byte) 0);

        xAxisValue = currStudy.getConfInterval(oddsDisplay)[1];
        yAxisValue = (float) i + 1;
        PlotPoint rightEnd = new PlotPoint(currStudy.getDisplayLabel(), currStudy.getShape(),
                                           xAxisValue, yAxisValue, (byte) 5, (byte) 0, (byte) 0);

        linesData.add(new GenericLine(leftEnd, rightEnd, (byte) 1, (byte) 0, (byte) 0, false));

        xAxisValue = currStudy.getBeta(oddsDisplay);
        yAxisValue = (float) i + 1;
        tempPoints.add(new PlotPoint(currStudy.getDisplayLabel() + "|"
                                     + prepareRightMarkers(currStudy), currStudy.getShape(),
                                     xAxisValue, yAxisValue, (byte) 3, (byte) 0, (byte) 0));
        tempPoints.get(tempPoints.size() - 1).setVisible(false);
      }
    }
    // points = tempPoints;
    points = tempPoints.toArray(new PlotPoint[tempPoints.size()]);
    lines = ArrayUtils.concatAll(lines, linesData.toArray(new GenericLine[linesData.size()]));
  }

  private String prepareRightMarkers(float beta, float conf0, float conf1) {
    if ((beta == 0.0f && conf0 == 0.0f && conf1 == 0.0f)
        || (oddsDisplay && beta == 1.0f && conf0 == 1.0f && conf1 == 1.0f)) {
      return " monomorphic ";
    }
    return String.format("%1$4s (%2$4s, %3$4s)", precision2Decimal.format(beta),
                         precision2Decimal.format(conf0), precision2Decimal.format(conf1));
  }

  private String prepareRightMarkers(StudyData forestTree) {
    return prepareRightMarkers(forestTree.getBeta(oddsDisplay),
                               forestTree.getConfInterval(oddsDisplay)[0],
                               forestTree.getConfInterval(oddsDisplay)[1]);
  }

  private void generateRectangles(Graphics g) {
    ArrayList<StudyData> currentData = forestPlot.getCurrentMetaStudy().getStudies();
    int breakCnt = 0;
    for (StudyData sd : currentData) {
      if (sd instanceof StudyBreak) {
        breakCnt++;
      }
    }
    ArrayList<GenericRectangle> rectData = new ArrayList<GenericRectangle>();
    rectangles = new GenericRectangle[0];

    float xAxisValue, yAxisValue;

    float xAxisStep = (float) calcStepStep(plotXmax - plotXmin);
    xAxisStep = Math.max(xAxisStep, 0.1f);
    // float yAxisStep = (float) calcStepStep(plotYmax - plotYmin);
    float yAxisStep = (float) calcStepStep(currentData.size() - breakCnt);
    yAxisStep = Math.max(yAxisStep, 0.1f);

    for (int i = 0; i < currentData.size(); i++) {
      if (currentData.get(i) instanceof StudyBreak) {
        continue;
      }
      if (currentData.get(i).getBeta(oddsDisplay) != 0
          && currentData.get(i).getStderr(oddsDisplay) != 0) {
        xAxisValue = currentData.get(i).getBeta(oddsDisplay);
        yAxisValue = (float) i + 1;
        // max scale = .25
        // float scale = currentData.get(i).getZScore() / forestPlot.getMaxZScore() / 4;

        float scale = (float) (Math.sqrt(currentData.get(i).getZScore())
                               / Math.sqrt(forestPlot.getMaxZScore()) / 4);
        float xDelta = xAxisStep * scale;
        float yDelta = yAxisStep * scale;
        xDelta = Math.max(xDelta, 0.05f);
        yDelta = Math.max(yDelta, 0.05f);
        rectData.add(new GenericRectangle(xAxisValue - xDelta, yAxisValue - yDelta,
                                          xAxisValue + xDelta, yAxisValue + yDelta, (byte) 5, true,
                                          false, (byte) 0, (byte) 0, false));
      }
    }

    rectangles = ArrayUtils.concatAll(rectangles,
                                      rectData.toArray(new GenericRectangle[rectData.size()]));
  }

  @Override
  public void highlightPoints() {
    byte defaultSize;
    defaultSize = 8;
    for (PlotPoint point : points) {
      if (point.isHighlighted()) {
        point.setSize((byte) (defaultSize * 1.5));
      } else {
        point.setSize((defaultSize));
      }

    }
  }

  @Override
  public void assignAxisLabels() {
    displayXAxis = displayYAxis = true;
    xAxisLabel = (oddsDisplay ? "Odds Ratio" : "Relative Risk") + " (" + forestPlot.getPlotLabel()
                 + ")";
    yAxisLabel = " ";
  }

  @Override
  public double[] getPlotMinMaxStep(double min, double max, Graphics g, boolean xAxis) {
    double range, plotStep, stepStep, plotMin, plotMax;
    double zoomMin, zoomMax, dist, tempD;
    int numHashes, textWidth, canvasRange;
    FontMetrics fontMetrics = g.getFontMetrics(g.getFont());
    int sf, temp;

    range = max - min;

    plotStep = stepStep = calcStepStep(range);
    sf = ext.getNumSigFig(stepStep);

    if (xAxis) {
      textWidth = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf, true)),
                           fontMetrics.stringWidth(ext.formDeci(max, sf, true)));
      numHashes = 12;
      canvasRange = canvasSectionMaximumX - canvasSectionMinimumX;
      temp = (textWidth + 30) * numHashes;
      while (temp > canvasRange) {
        numHashes -= 1;
        temp = (textWidth + 20) * numHashes;
      }
    } else {
      numHashes = 10;
    }
    numHashes = Math.max(numHashes, 1);

    while (range / plotStep >= numHashes) {
      plotStep += stepStep;
    }

    plotMin = plotMax = 0;
    while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
      plotMax += xAxis ? plotStep : 1;
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

    dist = plotMax - plotMin;
    zoomMin = plotMin + zoomSubsets[xAxis ? 0 : 1][0] * dist;
    zoomMax = plotMax - (1 - zoomSubsets[xAxis ? 0 : 1][1]) * dist;

    range = zoomMax - zoomMin;
    plotStep = stepStep = calcStepStep(range);
    sf = ext.getNumSigFig(stepStep);

    if (xAxis) {
      fontMetrics = g.getFontMetrics(g.getFont());
      textWidth = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf, true)),
                           fontMetrics.stringWidth(ext.formDeci(max, sf, true)));
      numHashes = 12;
      canvasRange = canvasSectionMaximumX - canvasSectionMinimumX;
      temp = (textWidth + 30) * numHashes;
      while (temp > canvasRange) {
        numHashes -= 1;
        temp = (textWidth + 20) * numHashes;
      }
    } else {
      numHashes = (canvasSectionMaximumY - canvasSectionMinimumY) / (fontMetrics.getHeight());
    }
    numHashes = Math.max(numHashes, 1);

    while (range / plotStep >= numHashes) {
      plotStep += stepStep;
    }

    double minBnd = Double.parseDouble(ext.formDeci(plotMin, sf, false));
    if (minBnd > zoomMin) {
      minBnd = zoomMin - (plotStep / 2d);
    }

    double[] retArr = new double[] {zoomMin, zoomMax,
                                    Double.parseDouble(ext.formDeci(plotStep, sf, true)), minBnd,
                                    sf};
    return retArr;
  }

  public double calcStepStep(double range) {
    String[] line;
    int temp;
    // double val;
    double val2;
    try {
      line = new DecimalFormat("0.0E0").format(range).split("E");
      temp = Integer.parseInt(line[1]) - 1;
      // val = Math.pow(10, temp); // ticks at 0.01
      val2 = Math.pow(10, temp) * (Double.parseDouble(line[0]) > 2.0 ? 5 : 1); // ticks at 0.05
      return val2;
    } catch (Exception e) {
      System.err.println("Error - could not parse stepStep from range '" + range + "'");
      return Double.NaN;
    }
  }

  @Override
  public void drawAll(Graphics g, boolean base) {
    float minimumObservedRawX, maximumObservedRawX, minimumObservedRawY, maximumObservedRawY;
    double[] plotMinMaxStep; // needs to be double, else x <= plotXmax can be inexact and leave off
                            // the last tick mark
    int sigFigs;
    String str;
    FontMetrics fontMetrics = g.getFontMetrics();
    // int[] order = null;
    int rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel;

    if (g instanceof Graphics2D) {
      ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                        antiAlias ? RenderingHints.VALUE_ANTIALIAS_ON
                                                  : RenderingHints.VALUE_ANTIALIAS_OFF);
      ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                                        antiAlias ? RenderingHints.VALUE_TEXT_ANTIALIAS_ON
                                                  : RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
    }

    setImageStatus(IMAGE_STARTED);
    generatePoints();

    if (/* !isPointsGeneratable() || */points.length == 0) {
      g.setColor(Color.WHITE);
      g.fillRect(0, 0, getWidth(), getHeight());
      if (getNullMessage() != null) {
        g.setColor(Color.BLACK);
        String msg = getNullMessage();
        boolean skip = false;
        int width = g.getFontMetrics().stringWidth(msg);
        int index = msg.length() / 2;
        if (width >= getWidth() - 20) {
          // search for a space
          while (msg.charAt(index) != ' ' && index < msg.length()) {
            index++;
          }
          if (index == msg.length()) {
            index = msg.length() / 2;
            // search for a space other direction
            while (msg.charAt(index) != ' ' && index >= 0) {
              index--;
            }
          }
          if (index != 0) {
            skip = true;
          }
        }

        int x = getWidth() / 2 - g.getFontMetrics(g.getFont()).stringWidth(msg) / 2;
        int y = getHeight() / 2;

        if (!skip) {
          g.drawString(msg, x, y);
        } else {
          String s1, s2;
          int x1, x2, y1, y2;
          s1 = msg.substring(0, index);
          s2 = msg.substring(index);
          x1 = getWidth() / 2 - (g.getFontMetrics().stringWidth(s1) / 2);
          x2 = getWidth() / 2 - (g.getFontMetrics().stringWidth(s2) / 2);
          y1 = getHeight() / 2;
          y2 = getHeight() / 2 + g.getFontMetrics().getHeight() + 3;
          g.drawString(s1, x1, y1);
          g.drawString(s2, x2, y2);
        }
      }

      setImageStatus(IMAGE_COMPLETE);
      return;
    }

    setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);
    assignAxisLabels();

    // TODO simplify and condense:
    // TODO THESE NEVER CHANGE - DON'T RE-CALC EVERY TIME!

    // Scan for rawX, rawY range of the data points
    minimumObservedRawX = Float.MAX_VALUE;
    maximumObservedRawX = Float.MIN_VALUE;
    minimumObservedRawY = Float.MAX_VALUE;
    maximumObservedRawY = Float.MIN_VALUE;
    for (int i = 0; i < points.length && isFlow(); i++) {
      if (points[i] != null) {
        minimumObservedRawX = Maths.min(minimumObservedRawX, points[i].getRawX());
        maximumObservedRawX = Maths.max(maximumObservedRawX, points[i].getRawX());
        minimumObservedRawY = Maths.min(minimumObservedRawY, points[i].getRawY());
        maximumObservedRawY = Maths.max(maximumObservedRawY, points[i].getRawY());
      }
    }

    for (int i = 0; lines != null && i < lines.length && isFlow(); i++) {
      if (lines[i] != null) {
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

    minimumObservedRawX = minimumObservedRawX == Float.MAX_VALUE ? 0 : minimumObservedRawX;
    maximumObservedRawX = maximumObservedRawX == Float.MIN_VALUE ? 1 : maximumObservedRawX;
    minimumObservedRawY = minimumObservedRawY == Float.MAX_VALUE ? 0 : minimumObservedRawY;
    maximumObservedRawY = maximumObservedRawY == Float.MIN_VALUE ? 1 : maximumObservedRawY;

    // otherwise step is off
    minimumObservedRawX = !oddsDisplay ? (minimumObservedRawX > 0 ? 0 : minimumObservedRawX)
                                       : minimumObservedRawX;
    minimumObservedRawY = minimumObservedRawY > 0 ? 0 : minimumObservedRawY;

    minimumObservedRawX = Float.isNaN(forcePlotXmin) ? minimumObservedRawX : forcePlotXmin;
    maximumObservedRawX = Float.isNaN(forcePlotXmax) ? (maximumObservedRawX
                                                        + (maximumObservedRawX
                                                           - minimumObservedRawX)
                                                          * (float) 0.01)
                                                     : forcePlotXmax;
    minimumObservedRawY = Float.isNaN(forcePlotYmin) ? minimumObservedRawY : forcePlotYmin;
    maximumObservedRawY = Float.isNaN(forcePlotYmax) ? (maximumObservedRawY
                                                        + (maximumObservedRawY
                                                           - minimumObservedRawY)
                                                          * (float) 0.01)
                                                     : forcePlotYmax;

    setNumberOfNaNSamples(0);
    int leftsize = 0;
    int rightsize = 0;
    if (base) {
      setAxisFontSize(22);
      g.fillRect(0, 0, getWidth(), getHeight());
      g.setFont(new Font("Arial", 0, axisFontSize));

      fontMetrics = g.getFontMetrics(g.getFont());
      missingWidth = fontMetrics.stringWidth("X");
      missingWidth = fontMetrics.stringWidth("X");

      // Calculate the plot area's range (X-axis, Y-axis)
      plotMinMaxStep = null;

      // float xRange = getXPixel(maximumObservedRawX) - getXPixel(minimumObservedRawX);
      // float yRange = getYPixel(minimumObservedRawY) - getYPixel(maximumObservedRawY);

      leftsize = determineLongestLeft(g, getMarkerFontSize(g));
      rightsize = determineRightBorder(g, getMarkerFontSize(g));

      if (displayXAxis) {
        canvasSectionMinimumX = WIDTH_BUFFER + leftsize;
        canvasSectionMaximumX = getWidth() - rightsize;
        canvasSectionMinimumY = 0;
        canvasSectionMaximumY = axisXHeight;// HEIGHT_X_AXIS;
        plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
        plotXmin = plotMinMaxStep[0];
        plotXmax = plotMinMaxStep[1];

        sigFigs = (int) plotMinMaxStep[4];
        // for (double x = plotMinMaxStep[3]; x <= plotXmax; x += plotMinMaxStep[2]) {
        double loopMax = plotXmin;
        while (loopMax < plotXmax) {
          loopMax += plotMinMaxStep[2];
        }
        for (double x = plotXmin; x <= plotXmax; x += plotMinMaxStep[2]) {
          if (x >= plotXmin || !truncate) {
            Grafik.drawThickLine(g, getXPixel(x), getHeight() - canvasSectionMaximumY, getXPixel(x),
                                 getHeight() - (canvasSectionMaximumY - TICK_LENGTH),
                                 TICK_THICKNESS, Color.BLACK);
            str = ext.formDeci(Math.abs(x) < DOUBLE_INACCURACY_HEDGE ? 0 : x, sigFigs, true);
            g.drawString(str, getXPixel(x) - str.length() * 8,
                         getHeight() - (canvasSectionMaximumY - TICK_LENGTH - 30));
          }
        }
        if (Math.abs(loopMax - plotXmax) > DOUBLE_INACCURACY_HEDGE) {
          Grafik.drawThickLine(g, getXPixel(plotXmax), getHeight() - canvasSectionMaximumY,
                               getXPixel(plotXmax),
                               getHeight() - (canvasSectionMaximumY - TICK_LENGTH), TICK_THICKNESS,
                               Color.BLACK);
          str = ext.formDeci(Math.abs(plotXmax) < DOUBLE_INACCURACY_HEDGE ? 0 : plotXmax, sigFigs,
                             true);
          g.drawString(str, getXPixel(plotXmax) - str.length() * 8,
                       getHeight() - (canvasSectionMaximumY - TICK_LENGTH - 30));
        }
        Grafik.drawThickLine(g, canvasSectionMinimumX - (int) Math.ceil(AXIS_THICKNESS / 2.0),
                             getHeight() - canvasSectionMaximumY,
                             canvasSectionMaximumX + (int) Math.ceil(AXIS_THICKNESS / 2.0),
                             getHeight() - canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
        g.drawString(xAxisLabel,
                     // graph-centered plot name:
                     canvasSectionMinimumX + ((canvasSectionMaximumX - canvasSectionMinimumX) / 2
                                              - (fontMetrics.stringWidth(xAxisLabel) / 2)),
                     // window-centered plot name:
                     // (getWidth() - WIDTH_Y_AXIS) / 2 - fontMetrics.stringWidth(xAxisLabel) / 2 +
                     // WIDTH_Y_AXIS,
                     getHeight() - 20);
      }

      if (displayYAxis) {
        g.setFont(new Font("Arial", 0, (int) getMarkerFontSize(g)));
        fontMetrics = g.getFontMetrics();
        canvasSectionMinimumX = 0;
        canvasSectionMaximumX = WIDTH_BUFFER + rightsize;
        canvasSectionMinimumY = /* HEIGHT_X_AXIS */axisXHeight + (fontMetrics.getHeight() * 2);
        canvasSectionMaximumY = getHeight() - (4 * HEAD_BUFFER);

        plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
        plotYmin = plotMinMaxStep[0];
        plotYmax = plotMinMaxStep[1];

        g.drawString((oddsDisplay ? "OR" : "RR") + " (95% CI)", getWidth() - rightsize + 15,
                     3 * HEAD_BUFFER + 14);

        sigFigs = (int) plotMinMaxStep[4];
        double step = 1;// Math.max(1, Math.round(plotMinMaxStep[2] * 2) / 2.0f);
        for (double y = plotMinMaxStep[3]; y <= plotYmax; y += step) {
          if ((y >= plotYmin && y == (int) y && y <= points.length && y > 0) || !truncate) {
            // Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(y),
            // canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
            str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
            str = str.split("\\.")[0];
            int index = Integer.parseInt(str) - 1;
            String left = points[index].getId().split("\\|")[0];
            g.drawString(left, WIDTH_BUFFER + leftsize - fontMetrics.stringWidth(left) - 15,
                         getYPixel(y) + 7);
            String right = points[index].getId().split("\\|")[1];
            g.drawString(right, getWidth() - rightsize + 15, getYPixel(y) + 7);
          }
        }
        g.drawString(META_LABEL, WIDTH_BUFFER + leftsize - fontMetrics.stringWidth(META_LABEL) - 15,
                     getHeight() - axisXHeight/* HEIGHT_X_AXIS */ - fontMetrics.getHeight() - 10);

        g.drawString(prepareRightMarkers(forestPlot.getCurrentMetaStudy().getMetaBeta(oddsDisplay),
                                         forestPlot.getCurrentMetaStudy()
                                                   .getMetaConf(oddsDisplay)[0],
                                         forestPlot.getCurrentMetaStudy()
                                                   .getMetaConf(oddsDisplay)[1]),
                     getWidth() - rightsize + 15,
                     getHeight() - /* HEIGHT_X_AXIS */axisXHeight - fontMetrics.getHeight() - 10);
        // Grafik.drawThickLine(g, canvasSectionMaximumX, getYPixel(plotYmin),
        // canvasSectionMaximumX, getYPixel(plotYmax) - (int) Math.ceil((double) TICK_THICKNESS /
        // 2.0), AXIS_THICKNESS, Color.BLACK);
        // g.setFont(new Font("Arial", 0, AXIS_FONT_SIZE));
        // yLabel = new BufferedImage(fontMetrics.stringWidth(yAxisLabel), 36,
        // BufferedImage.TYPE_INT_RGB);
        // gfx = yLabel.createGraphics();
        // gfx.setFont(new Font("Arial", 0, 28));
        // gfx.setColor(Color.WHITE);
        // gfx.fillRect(0, 0, getWidth(), getHeight());
        // gfx.setColor(Color.BLACK);
        // gfx.drawString(yAxisLabel, 0, yLabel.getHeight() - 6);
        //
        // g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight() - HEIGHT_X_AXIS) / 2 -
        // fontMetrics.stringWidth(yAxisLabel) / 2, this);
      }
      //
      // if (errorMessage != null) {
      // g.drawString(errorMessage, (getWidth() - WIDTH_Y_AXIS) / 2 -
      // fontMetrics.stringWidth(errorMessage) / 2 + WIDTH_Y_AXIS, (getHeight() - HEAD_BUFFER -
      // HEIGHT_X_AXIS) / 2 - 20 + HEAD_BUFFER);
      // }

    }

    // TODO outercoordinates
    canvasSectionMinimumX = WIDTH_BUFFER + leftsize;
    canvasSectionMaximumX = getWidth() - rightsize;
    canvasSectionMinimumY = /* HEIGHT_X_AXIS */axisXHeight + (fontMetrics.getHeight() * 2);
    canvasSectionMaximumY = getHeight() - (4 * HEAD_BUFFER);

    // Draw the lines
    for (int i = 0; lines != null && i < lines.length && isFlow(); i++) {
      if ((base && (getLayersInBase() == null
                    || Bytes.indexOf(getLayersInBase(), lines[i].getLayer()) >= 0))
          || (!base && Bytes.indexOf(getExtraLayersVisible(), lines[i].getLayer()) >= 0)) {
        Grafik.drawThickLine(g, getXPixel(lines[i].getStartX()), getYPixel(lines[i].getStartY()),
                             getXPixel(lines[i].getStopX()), getYPixel(lines[i].getStopY()),
                             lines[i].getThickness(), colorScheme[lines[i].getColor()]);
      }
    }

    if (isRectangleGeneratable()) {
      generateRectangles(g);
    }

    // Draw the rectangles for clusterFilters
    int actWidth, actHeight, wDiff, hDiff;
    for (int i = 0; rectangles != null && i < rectangles.length && isFlow(); i++) {
      if ((base && (getLayersInBase() == null
                    || Bytes.indexOf(getLayersInBase(), rectangles[i].getLayer()) >= 0))
          || (!base && Bytes.indexOf(getExtraLayersVisible(), rectangles[i].getLayer()) >= 0)) {
        rectangleXPixel = Math.min(getXPixel(rectangles[i].getStartXValue()),
                                   getXPixel(rectangles[i].getStopXValue()));
        rectangleYPixel = Math.min(getYPixel(rectangles[i].getStartYValue()),
                                   getYPixel(rectangles[i].getStopYValue()));
        rectangleWidthPixel = Math.abs(getXPixel(rectangles[i].getStartXValue())
                                       - getXPixel(rectangles[i].getStopXValue()));
        rectangleHeightPixel = Math.abs(getYPixel(rectangles[i].getStartYValue())
                                        - getYPixel(rectangles[i].getStopYValue()));

        actWidth = Math.min(rectangleWidthPixel, rectangleHeightPixel);
        actHeight = Math.min(rectangleWidthPixel, rectangleHeightPixel);

        wDiff = rectangleWidthPixel - actWidth;
        hDiff = rectangleHeightPixel - actHeight;

        rectangleXPixel += wDiff / 2;
        rectangleYPixel += hDiff / 2;

        g.setColor(colorScheme[rectangles[i].getColor()]);

        if (rectangles[i].getFill()) {
          if (rectangles[i].getRoundedCorners()) {
            g.fillRoundRect(rectangleXPixel, rectangleYPixel, actWidth, actHeight, 2, 2);
          } else {
            g.fillRect(rectangleXPixel, rectangleYPixel, actWidth, actHeight);
          }
        } else {
          if (rectangles[i].getRoundedCorners()) {
            g.drawRoundRect(rectangleXPixel, rectangleYPixel, actWidth, actHeight, 2, 2);
          } else {
            drawRectThick(g, rectangleXPixel, rectangleYPixel, actWidth, actHeight,
                          rectangles[i].getThickness());
          }
        }
      }
    }

    if (base) {
      g.setColor(Color.BLACK);
      double val = forestPlot.getCurrentMetaStudy().getMetaBeta(false)
                   - 1.96 * forestPlot.getCurrentMetaStudy().getMetaStderr(false);
      if (oddsDisplay) {
        val = Math.exp(val);
      }
      int xL = getXPixel(val);
      val = forestPlot.getCurrentMetaStudy().getMetaBeta(false);
      if (oddsDisplay) {
        val = Math.exp(val);
      }
      int xM = getXPixel(val);
      val = forestPlot.getCurrentMetaStudy().getMetaBeta(false)
            + 1.96 * forestPlot.getCurrentMetaStudy().getMetaStderr(false);
      if (oddsDisplay) {
        val = Math.exp(val);
      }
      int xR = getXPixel(val);
      // int xL = getXPixel(forestPlot.getCurrentMetaStudy().getMetaBeta(oddsDisplay) - 1.96 *
      // forestPlot.getCurrentMetaStudy().getMetaStderr(oddsDisplay));
      // int xM = getXPixel(forestPlot.getCurrentMetaStudy().getMetaBeta(oddsDisplay));
      // int xR = getXPixel(forestPlot.getCurrentMetaStudy().getMetaBeta(oddsDisplay) + 1.96 *
      // forestPlot.getCurrentMetaStudy().getMetaStderr(oddsDisplay));

      int yM = getHeight() - /* HEIGHT_X_AXIS */axisXHeight - fontMetrics.getHeight() - 15;
      int yU = yM - (fontMetrics.getHeight() / 2) - 1;
      int yD = yM + (fontMetrics.getHeight() / 2) + 1;

      Grafik.drawThickLine(g, xL, yM, xM, yU, 2, META_COLOR);
      Grafik.drawThickLine(g, xM, yU, xR, yM, 2, META_COLOR);
      Grafik.drawThickLine(g, xL, yM, xM, yD, 2, META_COLOR);
      Grafik.drawThickLine(g, xM, yD, xR, yM, 2, META_COLOR);

      int yMin = (4 * HEAD_BUFFER) - 5;
      int yMax = getHeight() - /* HEIGHT_X_AXIS */axisXHeight;

      Grafik.drawThickLine(g, getXPixel(oddsDisplay ? 1.0 : 0.0), yMin,
                           getXPixel(oddsDisplay ? 1.0 : 0.0), yMax, 3, Color.BLACK);

      int dashSize = 10;
      int dashSpacing = 5;
      int yStart;
      for (int i = 0; i < ((yMax - yMin) / (dashSize + dashSpacing)) + 1; i++) {
        yStart = yMin + (dashSize * i) + (dashSpacing * i);
        Grafik.drawThickLine(g, xM + 1, yStart, xM + 1, yStart + dashSize, 2, Color.GRAY);
      }

      g.setColor(Color.BLACK);
      g.setFont(new Font("Arial", Font.PLAIN, 22));
      if (forestPlot.getDataIndices().size() > 0) {
        String comm = forestPlot.getDataIndices().get(forestPlot.getCurrentDataIndex()).comment;
        if (!"".equals(comm)) {
          int w = g.getFontMetrics().stringWidth(comm) / 2;
          // TODO handle this better than checking GraphicsEnvironment.isHeadless() to determine
          // where to put title
          g.drawString(comm, getWidth() / 2 - w,
                       (GraphicsEnvironment.isHeadless() ? 1 : 3) * HEAD_BUFFER + 10);
        }
      }
    }

    setImageStatus(IMAGE_COMPLETE);
    refreshOtherComponents();
  }

  @Override
  public int getYPixel(double y) {
    double yMax = plotYmax; // # of studies
    double yMin = plotYmin; // always 0
    return getHeight()
           - (int) ((y - yMin) / (yMax - yMin) * (canvasSectionMaximumY - canvasSectionMinimumY)
                    + canvasSectionMinimumY);
  }

  private int determineRightBorder(Graphics g, double markerFontSize) {
    Font defaultFont = new Font("Arial", 0, (int) markerFontSize);
    FontMetrics fmt = g.getFontMetrics(defaultFont);
    int width = fmt.stringWidth("2222   2222  2222 2222");
    return (int) Math.floor(width);
  }

  private int determineLongestLeft(Graphics g, double markerFontSize) {
    Font defaultFont = new Font("Arial", 0, (int) markerFontSize);
    FontMetrics fmt = g.getFontMetrics(defaultFont);
    int width = fmt.stringWidth(forestPlot.getLongestStudyName());
    return (int) Math.floor(width);
  }

  private double getMarkerFontSize(Graphics g) {
    // double scale = 0.60;
    // double mkSz = ((getWidth() - WIDTH_BUFFER) - WIDTH_Y_AXIS) / (points.length * 2) * scale;
    // double mkSz = 20 - Math.atan(points.length - 15);
    // return mkSz;
    String largestLbl = forestPlot.getLongestStudyName();
    double yRange = getHeight() - (2 * HEAD_BUFFER) - (2 * axisXHeight/* HEIGHT_X_AXIS */);
    double studies = points.length;
    double xMax = getWidth() / 4;
    double yMax = yRange / studies;

    int fontSz = 20;
    FontMetrics fm = g.getFontMetrics(new Font("Arial", 0, fontSz));
    while (fm.getHeight() > yMax || fm.stringWidth(largestLbl) > xMax) {
      fontSz--;
      fm = g.getFontMetrics(new Font("Arial", 0, fontSz));
    }

    return fontSz;
  }

}
