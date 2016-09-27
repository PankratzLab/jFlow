package org.genvisis.one.ben.fcs.sub;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import org.genvisis.cnv.plots.GenericLine;
import org.genvisis.cnv.plots.GenericPath;
import org.genvisis.cnv.plots.GenericRectangle;
import org.genvisis.cnv.plots.PlotPoint;
import org.genvisis.common.Array;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2;
import org.genvisis.stats.LeastSquares;

public class OneDPanel extends AbstractPanel2 {
  public static final long serialVersionUID = 3L;
  public static final int LOOKUP_RESOLUTION = 20;
  public static final Color[] DEFAULT_COLORS = {
      new Color(33, 31, 53), // dark dark
      new Color(201, 30, 10), // deep red
      new Color(182, 182, 182), // light grey
      new Color(94, 88, 214), // light purple
      new Color(182, 182, 182, 180),
      new Color(189, 243, 61), // light green
      new Color(217, 109, 194), // pink
      new Color(33, 87, 0), // dark green
      new Color(23, 58, 172), // dark blue
      new Color(140, 20, 180), // deep purple
      new Color(220, 220, 220), // very light grey
      new Color(0, 0, 128), // ALL KINDS OF BLUES
      new Color(55, 129, 252), // light blue
      new Color(100, 149, 237), new Color(72, 61, 139), new Color(106, 90, 205),
      new Color(123, 104, 238), new Color(132, 112, 255), new Color(0, 0, 205),
      new Color(65, 105, 225), new Color(0, 0, 255), new Color(30, 144, 255),
      new Color(0, 191, 255), new Color(135, 206, 250), new Color(135, 206, 250),
      new Color(70, 130, 180), new Color(176, 196, 222), new Color(173, 216, 230),
      new Color(176, 224, 230), new Color(175, 238, 238), new Color(0, 206, 209),
      new Color(72, 209, 204), new Color(64, 224, 208), new Color(0, 255, 255),
      new Color(224, 255, 255),

  };
  private static final byte POINT_SIZE = 5;
  private static final int MISSING_SIZE_MULT = 3;

  public static enum PLOT_TYPE {
    BOX_PLOT, DOT_LINE_PLOT,
  }

  private PLOT_TYPE currentPlot;

  public void setPlotType(PLOT_TYPE type) {
    this.currentPlot = type;
    switch (type) {
      case BOX_PLOT:
        this.setAxisFontSize(12);
        rectangles = new GenericRectangle[0];

        setForcePlotXMax(20);
        setDisplayXAxis(false);
        break;
      case DOT_LINE_PLOT:

        setAxisFontSize(20);

        setForcePlotXMax(Float.NaN);
        setDisplayXAxis(true);

        break;
    }

    // for all:
    setDoubleBuffered(false);
    setSymmetricAxes(false);
    setYAxis(AXIS_SCALE.LIN);
    setXAxis(AXIS_SCALE.LIN);
    setZoomable(false, true);
    setColorScheme(DEFAULT_COLORS);
    createLookup(true);
    paintAgain();
  }

  public OneDPanel() {
    super();
    setPlotType(PLOT_TYPE.BOX_PLOT);
  }

  double[][] data;// = {11.8, 0.93, 1.76, 14, 16.5, 17.1, 32.5, 33.4, 16.8, 21.5, 13.1, 22.2, 22.2,
                  // 16, 16.2};
  String[][] dataLabels;
  String plotLabel;
  
  HashMap<String, int[][]> regressionLimits = new HashMap<String, int[][]>();
  HashMap<String, ArrayList<String>> locallyDroppedPoints = new HashMap<String, ArrayList<String>>();
  HashSet<String> globallyDroppedPoints = new HashSet<String>();
  
  public void setShowMean15Line(boolean showMean15Line) {
    this.showMean15Line = showMean15Line;
  }

  public void setShowRegressionLine(boolean showRegressionLine) {
    this.showRegressionLine = showRegressionLine;
  }

  public void setShow1SDLines(boolean show1sdLines) {
    show1SDLines = show1sdLines;
  }

  public void setShow2SDLines(boolean show2sdLines) {
    show2SDLines = show2sdLines;
  }

  private volatile boolean showMean15Line = true;
  private volatile boolean showRegressionLine = true;
  private volatile boolean show1SDLines = true;
  private volatile boolean show2SDLines = true;

  public void setData(String dataName, String[] files, double[] data) {
    this.plotLabel = dataName;
    this.dataLabels = new String[][] {files};
    this.data = new double[][] {data};
    if (!regressionLimits.containsKey(dataName)) {
      int[][] regLimits = new int[1][4];
      regLimits[0][0] = 0;
      regLimits[0][1] = data.length - 1;
      regLimits[0][2] = 0;
      regLimits[0][3] = data.length - 1;
      regressionLimits.put(dataName, regLimits);
    }
    if (!locallyDroppedPoints.containsKey(dataName)) {
      locallyDroppedPoints.put(dataName, new ArrayList<String>());
      for (String s : globallyDroppedPoints) {
        locallyDroppedPoints.get(dataName).add(s);
      }
    }
  }

  public void setData(String dataName, String[][] files, double[][] data) {
    this.plotLabel = dataName;
    this.dataLabels = files;
    this.data = data;
    if (!regressionLimits.containsKey(dataName)) {
      int[][] regLimits = new int[data.length][4];
      for (int d = 0; d < data.length; d++) {
        regLimits[d][0] = 0;
        regLimits[d][1] = data[d].length - 1;
        regLimits[d][2] = 0;
        regLimits[d][3] = data[d].length - 1;
        regressionLimits.put(dataName, regLimits);
      }
    }
    if (!locallyDroppedPoints.containsKey(dataName)) {
      locallyDroppedPoints.put(dataName, new ArrayList<String>());
      for (String s : globallyDroppedPoints) {
        locallyDroppedPoints.get(dataName).add(s);
      }
    }
  }

  public void generatePointsRectanglesAndLines() {
    if (data == null) {
      setNullMessage("Error, no data set..");
      points = new PlotPoint[0];
      rectangles = new GenericRectangle[0];
      lines = new GenericLine[0];
      polygons = new GenericPath[0];
      return;
    }

    switch (currentPlot) {
      case BOX_PLOT:
        generateBoxPlot();
        break;
      case DOT_LINE_PLOT:
        generateDotLinePlot();
        break;
    }

  }

  private void generateDotLinePlot() {
    byte type;
    float xAxisValue, yAxisValue;
    byte size = POINT_SIZE;

    int numPoints = 0;
    for (double[] dataArr : data) {
      numPoints += dataArr.length;
    }
    points = new PlotPoint[numPoints];
    ArrayList<GenericLine> lineList = new ArrayList<GenericLine>();
    ArrayList<GenericRectangle> rects = new ArrayList<GenericRectangle>();

    int ind = 0;
    // int lInd = 0;
    byte color;

    for (int d = 0; d < data.length; d++) {
      for (int i = 0; i < data[d].length; i++) {
        xAxisValue = (float) ind;
        yAxisValue = (float) data[d][i];
        if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
          type = PlotPoint.NOT_A_NUMBER;
        } else {
          if (showRegressionLine) {
            if (locallyDroppedPoints.get(plotLabel).contains(dataLabels[d][i])) {
              type = PlotPoint.MISSING;
              size = (byte) ((int) POINT_SIZE * MISSING_SIZE_MULT);
            } else {
              type = PlotPoint.FILLED_CIRCLE;
              size = POINT_SIZE;
            }
          } else {
            type = PlotPoint.FILLED_CIRCLE;
            size = POINT_SIZE;
          }
        }
        
        color = (byte) d; // TODO apply gating for colors
        points[ind] = new PlotPoint(dataLabels[d][i], type, xAxisValue, yAxisValue, size, color, (byte) 0);
        if (i < data[d].length - 1) {
          lineList.add(new GenericLine(xAxisValue, yAxisValue, (float) ind + 1, (float) data[d][i + 1], (byte) 1, (byte) d, (byte) 0));
        }
        ind++;
      }
    }

    float mean = (float) Array.mean(data[0]);
    float mean15 = mean * .15f;
    float sd = (float) Array.stdev(data[0], true);
    lineList.add(new GenericLine(0, mean, numPoints + 1, mean, (byte) 1, (byte) 0, (byte) 99));
    if (showMean15Line) {
      lineList.add(new GenericLine(0, mean - mean15, numPoints + 1, mean - mean15, (byte) 1, (byte) getMeanColor(), (byte) 0));
      lineList.add(new GenericLine(0, mean + mean15, numPoints + 1, mean + mean15, (byte) 1, (byte) getMeanColor(), (byte) 0));
    }
    if (show1SDLines) {
      lineList.add(new GenericLine(0, mean - sd, numPoints + 1, mean - sd, (byte) 1, (byte) get1SDColor(), (byte) 0));
      lineList.add(new GenericLine(0, mean + sd, numPoints + 1, mean + sd, (byte) 1, (byte) get1SDColor(), (byte) 0));
    }
    if (show2SDLines) {
      lineList.add(new GenericLine(0, mean - 2 * sd, numPoints + 1, mean - 2 * sd, (byte) 1, (byte) get2SDColor(), (byte) 0));
      lineList.add(new GenericLine(0, mean + 2 * sd, numPoints + 1, mean + 2 * sd, (byte) 1, (byte) get2SDColor(), (byte) 0));
    }

    double dataMin = Math.min(mean - mean15, Array.min(data[0])), dataMax = Math.max(mean + mean15, Array.max(data[0]));
    for (int i = 1; i < data.length; i++) {
      if (data[i].length == 0)
        continue;
      dataMin = Math.min(dataMin, Array.min(data[i]));
      dataMax = Math.max(dataMax, Array.max(data[i]));
    }

    double rngY = dataMax - dataMin;
    double pad = 1;
    while ((rngY = rngY / 10) > 1) {
      pad += 1;
    }
    
    if (showRegressionLine) {
      int prevSum = 0;
      for (int d = 0; d < data.length; d++) {
        int[] limits = regressionLimits.get(plotLabel)[d];
        ArrayList<String> dropped = locallyDroppedPoints.get(plotLabel);
        int count = limits[1] - limits[0];
        for (int i = limits[0]; i <= limits[1]; i++) {
          if (dropped.contains(dataLabels[d][i])) {
            count--;
          }
        }
        
        if (count >= 3) {
          double[] deps = new double[count + 1];
          double[] indeps = new double[count + 1];
          int dropCnt = 0;
          for (int i = limits[0]; i <= limits[1]; i++) {
            if (dropped.contains(dataLabels[d][i])) {
              dropCnt++;
            } else {
              indeps[i - dropCnt - limits[0]] = i - dropCnt - limits[0];
              deps[i - dropCnt - limits[0]] = data[d][i];
            }
          }
          try {
            LeastSquares ls = new LeastSquares(deps, indeps);
            if (!ls.analysisFailed()) {
              double b = ls.getBetas()[0];
              double m = ls.getBetas()[1];
              double yStart = b;
              double yEnd = b + /*(limits[1] - limits[0] + 1) +*/ m;
              lineList.add(new GenericLine(prevSum, (float) yStart, prevSum + data[d].length - 1, (float) yEnd, (byte) 2, (byte) d, (byte) 0));
            } else {
              System.err.println("Error - regression failed");
            }
          } catch (Exception e) {
            System.err.println("Error - regression failed with exception: " + e.getMessage());
          }

          rects.add(new GenericRectangle((float) (prevSum + limits[2]), (float) (dataMin - pad), (float) (prevSum + limits[0]), (float) (dataMax + pad), (byte) 1, true, false, (byte) 0, (byte) 4, (byte) 99, false));
          rects.add(new GenericRectangle((float) (prevSum + limits[1]), (float) (dataMin - pad), (float) (prevSum + limits[3]), (float) (dataMax + pad), (byte) 1, true, false, (byte) 0, (byte) 4, (byte) 99, false));
        }
        prevSum += data[d].length;
      }
    }

    setForcePlotYMin((float) (dataMin - pad));
    setForcePlotYMax((float) (dataMax + pad));
    setForcePlotXMin(0);

    setYAxis(AXIS_SCALE.LIN);
    setXAxis(AXIS_SCALE.LIN);
    
    lines = lineList.toArray(new GenericLine[lineList.size()]);
    rectangles = new GenericRectangle[rects.size()];
    rectangles = rects.toArray(rectangles);
  }

  private void generateBoxPlot() {
    // points for any data above/below wiskLow/wiskHigh
    ArrayList<GenericLine> lns = new ArrayList<GenericLine>();
    ArrayList<PlotPoint> pts = new ArrayList<PlotPoint>();

    double xMin = Double.MAX_VALUE, xMax = Double.MIN_VALUE;
    double min = Double.MAX_VALUE, max = Double.MIN_VALUE;
    for (int i = 0; i < data.length; i++) {
      if (data[i].length == 0)
        continue;
      if (data[i].length <= 2) {
        float xLow = 20 * (i - 1) + 2;
        float xHigh = xLow + 18;

        byte col = (byte) i;

        for (int d = 0; d < data[i].length; d++) {
          min = Math.min(min, data[i][d]);
          max = Math.max(max, data[i][d]);
          lns.add(new GenericLine(xLow, (float) data[i][d], xHigh, (float) data[i][d], (byte) 4,
              col, (byte) 0));
        }
      } else if (data[i].length > 2) {
        float xLow = 20 * i + 2;
        xMin = Math.min(xLow - 5, xMin);
        float xHigh = xLow + 18;
        xMax = Math.max(xHigh + 5, xMax);
        float xMed = xLow + (xHigh - xLow) / 2;

        byte col = (byte) i;

        double med = Array.median(data[i]);
        double qr25 = Array.quantExclusive(data[i], 0.25);
        double qr75 = Array.quantExclusive(data[i], 0.75);
        double iqr = Array.iqrExclusive(data[i]);
        double wiskLow = qr25 - 1.5 * iqr;
        double wiskHigh = qr75 + 1.5 * iqr;
        min = Math.min(min, Math.min(wiskLow, Array.min(data[i])));
        max = Math.max(max, Math.max(wiskHigh, Array.max(data[i])));

        // line @ med
        lns.add(new GenericLine(xLow, (float) med, xHigh, (float) med, (byte) 4, col, (byte) 0));
        // line @ qr25
        lns.add(new GenericLine(xLow, (float) qr25, xHigh, (float) qr25, (byte) 2, col, (byte) 0));
        // line @ qr75
        lns.add(new GenericLine(xLow, (float) qr75, xHigh, (float) qr75, (byte) 2, col, (byte) 0));
        // small line at wiskLow
        lns.add(new GenericLine(xMed - (xMed - xLow) / 2, (float) wiskLow,
            xMed + (xMed - xLow) / 2, (float) wiskLow, (byte) 2, col, (byte) 0));
        // small line at wiskHigh
        lns.add(new GenericLine(xMed - (xMed - xLow) / 2, (float) wiskHigh, xMed + (xMed - xLow)
            / 2, (float) wiskHigh, (byte) 2, col, (byte) 0));
        // line from qr25 -> wiskLow
        lns.add(new GenericLine(xMed, (float) qr25, xMed, (float) wiskLow, (byte) 1, col, (byte) 0));
        // line from qr75 -> wiskHigh
        lns.add(new GenericLine(xMed, (float) qr75, xMed, (float) wiskHigh, (byte) 1, col, (byte) 0));
        // two lines vert, from qr25 to qr75
        lns.add(new GenericLine(xLow, (float) qr25, xLow, (float) qr75, (byte) 1, col, (byte) 0));
        lns.add(new GenericLine(xHigh, (float) qr25, xHigh, (float) qr75, (byte) 1, col, (byte) 0));

        for (int j = 0; j < data[i].length; j++) {
          if (data[i][j] < wiskLow || data[i][j] > wiskHigh) {
            pts.add(new PlotPoint(dataLabels[i][j], PlotPoint.FILLED_CIRCLE, xMed,
                (float) data[i][j], (byte) POINT_SIZE, col, (byte) 0));
          }
        }
      }
    }

    double rngY = max - min;
    double pad = 1;
    while ((rngY = rngY / 10) > 1) {
      pad += 1;
    }

    setForcePlotYMin((float) (min - pad));
    setForcePlotYMax((float) (max + pad));
    setPlotYMin((float) (min - pad));
    setPlotYMax((float) (max + pad));

    setForcePlotXMin((float) (xMin));
    setForcePlotXMax((float) (xMax));
    setPlotXMin((float) (xMin));
    setPlotXMax((float) (xMax));


    lines = lns.toArray(new GenericLine[lns.size()]);
    points = pts.toArray(new PlotPoint[pts.size()]);

  }
  
  volatile int dragInd = -1;
  volatile int drag = 0;
  
  @Override
  public void mouseReleased(MouseEvent e) {
    if (currentPlot == PLOT_TYPE.BOX_PLOT) return;
    dragInd = -1;
    drag = 0;
  }
  
  @Override
  public void mousePressed(MouseEvent e) {
    if (currentPlot == PLOT_TYPE.BOX_PLOT) return;
    
    dragInd = -1;
    drag = 0;
    int prevSum = 0;
    if (!showRegressionLine) return;
    for (int d = 0; d < data.length; d++) {
      double diffStart = Math.abs(getXValueFromXPixel(e.getX()) - regressionLimits.get(plotLabel)[d][0] - prevSum);
      double diffEnd = Math.abs(getXValueFromXPixel(e.getX()) - regressionLimits.get(plotLabel)[d][1] - prevSum);
    
      if (diffStart < 1) {
        dragInd = d;
        drag = 1;
      } else if (diffEnd < 1) {
        dragInd = d;
        drag = 2;
      }
      prevSum += data[d].length;
    }
    
    paintAgain();
  }
  
  @Override
  public void mouseDragged(MouseEvent e) {
    if (currentPlot == PLOT_TYPE.BOX_PLOT) return;
    if (dragInd == -1) return;
    int sum = 0;
    if (dragInd > 0) {
      for (int i = 0; i < dragInd; i++) {
        sum += data[i].length;
      }
    }
    switch(drag) {
      case 0:
        return;
      case 1:
        regressionLimits.get(plotLabel)[dragInd][0] = Math.round((float) getXValueFromXPixel(e.getX()) - sum);
        regressionLimits.get(plotLabel)[dragInd][0] = Math.max(regressionLimits.get(plotLabel)[dragInd][0], regressionLimits.get(plotLabel)[dragInd][2]);
        break;
      case 2:
        regressionLimits.get(plotLabel)[dragInd][1] = Math.round((float) getXValueFromXPixel(e.getX()) - sum);
        regressionLimits.get(plotLabel)[dragInd][1] = Math.min(regressionLimits.get(plotLabel)[dragInd][1], regressionLimits.get(plotLabel)[dragInd][3]);
        break;
    }
    paintAgain();
  }

  public void mouseClicked(MouseEvent e) {
    JPopupMenu menu;
    byte maxNumPoints;

    if (prox != null && prox.size() > 0) {
      menu = new JPopupMenu();
      maxNumPoints = (byte) Math.min(20, prox.size());
      for (int i = 0; i < maxNumPoints; i++) {
        final int ind = prox.elementAt(i);
        JMenuItem jmi = new JMenuItem(points[ind].getId() + " -- " + points[ind].getRawY());
        jmi.addActionListener(new ActionListener() {
          @Override
          public void actionPerformed(ActionEvent e) {
            ext.setClipboard(points[ind].getId());
          }
        });
        menu.add(jmi);
        if (showRegressionLine) {
          
          boolean localDrop, globalDrop;
          localDrop = (locallyDroppedPoints.get(plotLabel).contains(points[ind].getId()));
          globalDrop = (globallyDroppedPoints.contains(points[ind].getId()));
          
          String promptLocal = null, promptGlobal = null;
          int caseVal = -1;
          if (globalDrop && localDrop) {
            promptLocal = "Re-add data-point to this regression only";
            promptGlobal = "Re-add data-point to all regressions";
            caseVal = 0;
          } else if (globalDrop && !localDrop) {
            promptLocal = "Drop data-point from this regression only";
            promptGlobal = "Re-add data-point to all regressions";
            caseVal = 1;
          } else if (!globalDrop && localDrop) {
            promptLocal = "Re-add data-point to this regression only";
            promptGlobal = "Drop data-point from all regressions";
            caseVal = 2;
          } else if (!globalDrop && !localDrop) {
            promptLocal = "Drop data-point from this regression only";
            promptGlobal = "Drop data-point from all regressions";
            caseVal = 3;
          }
          
          if (caseVal >= 0) { 
            final int caseValF = caseVal;
          
            JMenuItem jmi2 = new JMenuItem(promptLocal);
            jmi2.addActionListener(new ActionListener() {
              @Override
              public void actionPerformed(ActionEvent e) {
                switch (caseValF) {
                  case 0:
                    locallyDroppedPoints.get(plotLabel).remove(points[ind].getId());
                    break;
                  case 1:
                    locallyDroppedPoints.get(plotLabel).add(points[ind].getId());
                    break;
                  case 2:
                    locallyDroppedPoints.get(plotLabel).remove(points[ind].getId());
                    break;
                  case 3:
                    locallyDroppedPoints.get(plotLabel).add(points[ind].getId());
                    break;
                }
                firePropertyChange("REGRESSION", null, null);
                paintAgain();
              }
            });
            menu.add(jmi2);
            
            JMenuItem jmi3 = new JMenuItem(promptGlobal);
            jmi3.addActionListener(new ActionListener() {
              @Override
              public void actionPerformed(ActionEvent e) {
                switch (caseValF) {
                  case 0:
                    globallyDroppedPoints.remove(points[ind].getId());
                    for (ArrayList<String> s : locallyDroppedPoints.values()) {
                      s.remove(points[ind].getId());
                    }
                    break;
                  case 1:
                    globallyDroppedPoints.remove(points[ind].getId());
                    for (ArrayList<String> s : locallyDroppedPoints.values()) {
                      s.remove(points[ind].getId());
                    }
                    break;
                  case 2:
                    globallyDroppedPoints.add(points[ind].getId());
                    for (ArrayList<String> s : locallyDroppedPoints.values()) {
                      s.add(points[ind].getId());
                    }
                    break;
                  case 3:
                    globallyDroppedPoints.add(points[ind].getId());
                    for (ArrayList<String> s : locallyDroppedPoints.values()) {
                      s.add(points[ind].getId());
                    }
                    break;
                }
                firePropertyChange("REGRESSION", null, null);
                paintAgain();
              }
            });
            menu.add(jmi3);
          }
        }
      }
      menu.show(this, e.getX(), e.getY());
    }
  }

  @Override
  public void highlightPoints() {
    byte defaultSize;

    defaultSize = POINT_SIZE;
    for (int i = 0; i < points.length; i++) {
      points[i].setSize((byte) ((points[i].getType() == PlotPoint.MISSING ? defaultSize * MISSING_SIZE_MULT : defaultSize) * (points[i].isHighlighted() ? 1.5 : 1)));
    }
  }

  @Override
  public void assignAxisLabels() {
    // String[] pts = plotLabel.split("\\|");
    // setXAxisLabel("");//pts[0].trim().replaceAll("/", " /\n");
    // setYAxisLabel(pts[1].trim());

    // setXAxisLabel("File by Date");
    // setYAxisLabel("Mean - " + plotLabel);
  }

  public boolean isShowMean15Line() {
    return showMean15Line;
  }

  public boolean isShowRegressionLine() {
    return showRegressionLine;
  }

  public boolean isShow1SDLines() {
    return show1SDLines;
  }

  public boolean isShow2SDLines() {
    return show2SDLines;
  }

  public int get1SDColor() {
    return 6;
  }

  public int get2SDColor() {
    return 3;
  }

  public int getMeanColor() {
    return 2;
  }

}
