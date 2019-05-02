package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.AbstractAction;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.LaunchAction;
import org.genvisis.cnv.plots.ManhattanPlot.ManhattanDataPoint;
import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.pankratzlab.common.Grafik;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;

public class ManhattanPanel extends AbstractPanel {

  ManhattanPlot mp;
  byte pointSize = 2;
  byte layer = 1;
  int numPointColors = 3;

  double[] lineValuesToDraw = new double[] {.00001, 0.00000005};
  int[] sizeMult = new int[] {2, 4};
  int[] lineColors = new int[] {3, 3};
  int[] aboveLineColors = new int[] {4, 5};

  HashMap<Integer, int[]> linearizedChrBnds = new HashMap<>();

  public ManhattanPanel(ManhattanPlot parent) {
    super();
    this.mp = parent;
    setZoomable(true, true);
    setSymmetricAxes(false);

    setColorScheme(new Color[] {Color.GRAY, Color.BLACK, Color.LIGHT_GRAY, Color.RED, Color.MAGENTA,
                                Color.GREEN});

  }

  @Override
  protected void drawXAxis(Graphics g, double[] plotMinMaxStep, FontMetrics fontMetrics) {
    String str;
    if (displayXAxisScale) {
      int prevEnd = -1;
      for (int i : linearizedChrBnds.keySet()) {
        int[] bnds = linearizedChrBnds.get(i);
        int x = (int) (bnds[0] + .5 * (bnds[1] - bnds[0]));
        if (x >= plotXmin || !truncate) {
          str = Positions.CHR_CODES[i];
          int xLoc = getXPixel(x) - fontMetrics.stringWidth(str) / 2;
          int len = TICK_LENGTH;
          if (xLoc <= prevEnd) {
            len -= TICK_LENGTH / 3;
          }
          Grafik.drawThickLine(g, getXPixel(x), getHeight() - canvasSectionMaximumY, getXPixel(x),
                               getHeight() - (canvasSectionMaximumY - len), TICK_THICKNESS,
                               Color.BLACK);
          if (xLoc > prevEnd) {
            g.drawString(str, xLoc, getHeight() - (canvasSectionMaximumY - TICK_LENGTH - 30));
            prevEnd = xLoc + fontMetrics.stringWidth(str) + 1;
          }
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
      g.drawString(xAxisLabel,
                   (getWidth() - getAxisYWidth()/* WIDTH_Y_AXIS */)
                               / 2 - fontMetrics.stringWidth(xAxisLabel) / 2
                               + getAxisYWidth()/* WIDTH_Y_AXIS */,
                   getHeight() - 20);
    }
  }

  @Override
  public void generatePoints() {
    ArrayList<ManhattanDataPoint> dataPoints = mp.getData();
    if (dataPoints == null) {
      points = new PlotPoint[0];
      lines = new GenericLine[0];
      setNullMessage("Please Load a Data File.");
      return;
    } else if (dataPoints.isEmpty()) {
      points = new PlotPoint[0];
      lines = new GenericLine[0];
      if (mp.isDataLoaded()) {
        setNullMessage("No Data Passed Filters.");
      } else {
        setNullMessage("Data Loading, Please Wait...");
      }
      return;
    } else {
      setNullMessage(null);
    }
    setForcePlotXmin(dataPoints.get(0).linearLoc - ManhattanPlot.LIN_CHR_BUFFER);
    setForcePlotXmax(dataPoints.get(dataPoints.size() - 1).linearLoc
                     + ManhattanPlot.LIN_CHR_BUFFER);
    points = new PlotPoint[dataPoints.size()];
    for (int i = 0, count = dataPoints.size(); i < count; i++) {
      ManhattanDataPoint mdp = dataPoints.get(i);
      int[] bnds = linearizedChrBnds.get(mdp.chr);
      if (bnds == null) {
        bnds = new int[] {Integer.MAX_VALUE, Integer.MIN_VALUE};
        linearizedChrBnds.put(mdp.chr, bnds);
      }
      if (mdp.linearLoc < bnds[0]) {
        bnds[0] = mdp.linearLoc;
      }
      if (mdp.linearLoc > bnds[1]) {
        bnds[1] = mdp.linearLoc;
      }
      points[i] = new PlotPoint(mdp.mkr == null ? mdp.chr + ":" + mdp.pos : mdp.mkr,
                                PointType.FILLED_CIRCLE, (float) mdp.linearLoc,
                                (float) mdp.transformedPVal, getSize(mdp.transformedPVal),
                                getPointColor(mdp.chr, mdp.transformedPVal), layer);
    }

    lines = new GenericLine[lineValuesToDraw.length];
    for (int i = 0; i < lineValuesToDraw.length; i++) {
      float v = (float) -Math.log10(lineValuesToDraw[i]);
      lines[i] = new GenericLine(Integer.MIN_VALUE, v, (float) Integer.MAX_VALUE, v, (byte) 1,
                                 (byte) lineColors[i], (byte) 99);
    }
  }

  private byte getPointColor(int chr, double transP) {
    int c = (chr % numPointColors);
    for (int i = 0; i < lineValuesToDraw.length; i++) {
      float v = (float) -Math.log10(lineValuesToDraw[i]);
      if (transP < v) {
        break;
      }
      c = aboveLineColors[i];
    }
    return (byte) c;
  }

  private byte getSize(double pVal) {
    int sz = pointSize;
    for (int i = 0; i < lineValuesToDraw.length; i++) {
      float v = (float) -Math.log10(lineValuesToDraw[i]);
      if (pVal < v) {
        break;
      }
      sz = pointSize * sizeMult[i];
    }
    return (byte) sz;
  }

  @Override
  public void highlightPoints() {
    // TODO Auto-generated method stub

  }

  @Override
  public void mouseClicked(MouseEvent e) {
    JPopupMenu menu;

    if (SwingUtilities.isRightMouseButton(e)) {
      if (prox != null && prox.size() > 0) {
        if (prox.size() <= 10) {
          menu = new JPopupMenu();
          Project proj = mp.getProject();
          for (int i = 0; i < prox.size(); i++) {
            ManhattanDataPoint mdp = mp.getData().get(prox.get(i));

            JMenu subMen = new JMenu(mdp.mkr);
            JMenuItem jmi;

            jmi = new JMenuItem();
            jmi.setEnabled(false);
            jmi.setText("P-val: " + mdp.originalPVal + " | " + mdp.transformedPVal);
            subMen.add(jmi);

            jmi = new JMenuItem(new AbstractAction() {

              @Override
              public void actionPerformed(ActionEvent arg0) {
                ext.setClipboard(mdp.mkr);
              }
            });
            jmi.setText("Copy Name to Clipboard");
            subMen.add(jmi);

            jmi = new JMenuItem(new AbstractAction() {

              @Override
              public void actionPerformed(ActionEvent e) {
                ext.setClipboard(mdp.chr + ":" + mdp.pos);
              }
            });
            jmi.setText("Copy Position to Clipboard");
            subMen.add(jmi);

            jmi = new JMenuItem(new AbstractAction() {

              @Override
              public void actionPerformed(ActionEvent e) {
                ext.setClipboard(mdp.mkr + "\t" + mdp.chr + ":" + mdp.pos);
              }
            });
            jmi.setText("Copy Name & Position to Clipboard");
            subMen.add(jmi);

            if (proj != null) {
              LaunchAction act = new LaunchAction(proj, mdp.mkr,
                                                  colorScheme[getPointColor(mdp.chr,
                                                                            mdp.transformedPVal)]);
              subMen.add(act);
            }

            for (Entry<String, String> other : mdp.otherData.entrySet()) {
              jmi = new JMenuItem();
              jmi.setEnabled(false);
              jmi.setText(other.getKey() + ": " + other.getValue());
              subMen.add(jmi);
            }

            menu.add(subMen);
          }
          menu.show(this, e.getX(), e.getY());
        }
      }
    }
  }

  @Override
  public void mouseMoved(MouseEvent event) {
    super.mouseMoved(event);

    ToolTipManager.sharedInstance().setReshowDelay(3);
    ToolTipManager.sharedInstance().setInitialDelay(3);

    StringBuilder sb = new StringBuilder("<html>");
    if (prox != null && prox.size() > 0) {
      if (prox.size() == 1) {
        sb.append(points[prox.get(0)].getId());
      } else if (prox.size() <= 10) {
        for (int i = 0; i < prox.size(); i++) {
          ManhattanDataPoint mdp = mp.getData().get(prox.get(i));
          sb.append(mdp.mkr);
          if (i < prox.size() - 1) {
            sb.append("<br/>");
          }
        }
      } else {
        sb = null;
      }
      if (sb != null) {
        sb.append("</html>");
        setToolTipText(sb.toString());
      } else {
        setToolTipText(null);
      }
    } else {
      setToolTipText(null);
    }

  };

  @Override
  public void assignAxisLabels() {
    xAxisLabel = "Chromosome";
    yAxisLabel = "-Log10(pVal)";
  }

}
