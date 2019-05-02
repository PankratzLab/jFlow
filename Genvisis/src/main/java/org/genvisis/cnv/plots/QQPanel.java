package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;

import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.ext;

public class QQPanel extends AbstractPanel {

  public static final long serialVersionUID = 1L;

  protected static final Font DESCR_FONT = new Font("Arial", 0, 20);
  protected static final int PAD = 5;

  private final double[][] pvals;
  private final boolean log10;
  private final boolean rotated;
  private final float maxValue;
  private final String[] descriptions;

  public QQPanel(String[] labels, double[][] pvals, boolean[][] useForLambda, boolean log10,
                 boolean rotated, float maxValue, Color[] colorScheme, Logger log) {
    super();

    this.setLog(log);
    this.pvals = pvals;
    this.log10 = log10;
    this.rotated = rotated;
    this.maxValue = maxValue;
    this.colorScheme = colorScheme;

    descriptions = new String[pvals.length];

    log.report("File\tTrait\tLambda");
    for (int i = 0; i < pvals.length; i++) {
      String lamb = ext.formDeci(ArrayUtils.lambda(ArrayUtils.subArray(pvals[i], useForLambda[i])),
                                 4);
      descriptions[i] = "lambda = " + lamb + " (" + labels[i] + ")";
      log.report(ArrayUtils.toStr(ext.replaceAllWith(labels[i], "'", "")
                                     .split(PSF.Regex.GREEDY_WHITESPACE))
                 + "\t" + lamb);
    }

    createLookup(false);

    setZoomable(true, true);
  }

  protected double[][] getPvals() {
    return pvals;
  }

  @Override
  public void drawTitle(Graphics g, boolean base, FontMetrics fontMetrics) {
    if (base) {
      Color currColor = g.getColor();
      fontMetrics = g.getFontMetrics(DESCR_FONT);
      int fontHeight = (fontMetrics == null ? 25 : fontMetrics.getHeight());
      int descrY = 0;

      for (int i = 0; i < descriptions.length; i++) {
        int descrX = PAD + getAxisYWidth();
        descrY += calcSingleDescrHeight(fontHeight);

        g.setColor(descriptions.length == 1 ? colorScheme[0] : colorScheme[i + 2]);
        g.drawString(descriptions[i], descrX, descrY);
      }
      g.setColor(currColor);
    }
  }

  @Override
  public int calcTitleHeight(Graphics g, boolean base, FontMetrics fontMetrics) {
    if (base) {
      fontMetrics = g.getFontMetrics(DESCR_FONT);
      int fontHeight = (fontMetrics == null ? 25 : fontMetrics.getHeight());
      return calcSingleDescrHeight(fontHeight) * descriptions.length;
    }
    return 0;
  }

  private static int calcSingleDescrHeight(int fontHeight) {
    return fontHeight + PAD * 2;
  }

  @Override
  public void assignAxisLabels() {
    if (log10) {
      if (rotated) {
        xAxisLabel = "-log10(rank/n)";
        yAxisLabel = "-log10(p-value) - -log10(rank/n)";
      } else {
        xAxisLabel = "-log10(rank/n)";
        yAxisLabel = "-log10(p-value)";
      }
      plotXmin = 0;
      plotYmin = 0;
    } else {
      xAxisLabel = "Expected quantiles";
      yAxisLabel = "Observed quantiles";
    }
  }

  public boolean invertX() {
    return false;
  }

  public boolean invertY() {
    return false;
  }

  @Override
  public void highlightPoints() {}

  @Override
  public void generatePoints() {
    int[] keys;
    int count;
    int max;

    count = 0;
    for (double[] pval : pvals) {
      count += pval.length;
    }

    // if (lines != null && points != null && points.length == count && count > 0) {
    // return;
    // }
    lines = new GenericLine[1];
    max = 0;
    for (double[] pval : pvals) {
      max = Math.max(max, pval.length);
    }
    max = (int) Math.ceil((-1 * Math.log10(1.0 / max)));
    if (rotated) {
      lines[0] = new GenericLine(0, 0, max, 0, (byte) 2, (byte) 1, (byte) 0);
    } else if (log10) {
      // lines[0] = new PlotLine(0, 0, (float)(-1*Math.log10((1.0/pvals.length))),
      // (float)(-1*Math.log10((1.0/pvals.length))), (byte)2, (byte)1);
      lines[0] = new GenericLine(0, 0, max, max, (byte) 2, (byte) 1, (byte) 0);
    } else {
      lines[0] = new GenericLine(0, 0, 1, 1, (byte) 2, (byte) 1, (byte) 0);
    }

    points = new PlotPoint[count];

    count = 0;
    for (int i = 0; i < pvals.length; i++) {
      keys = Sort.getSortedIndices(pvals[i]);

      for (int j = 0; j < pvals[i].length; j++) {
        if (rotated) {
          points[count] = new PlotPoint(keys[j] + "", PointType.FILLED_CIRCLE,
                                        (float) (-1 * Math.log10(((double) keys[j] + 1)
                                                                 / pvals[i].length)),
                                        Math.min(maxValue,
                                                 (float) (-1 * Math.log10(pvals[i][j]))
                                                           - (float) (-1
                                                                      * Math.log10(((double) keys[j]
                                                                                    + 1)
                                                                                   / pvals[i].length))),
                                        (byte) 6, (byte) (pvals.length == 1 ? 0 : i + 2), (byte) 0);
        } else if (log10) {
          points[count] = new PlotPoint(keys[j] + "", PointType.FILLED_CIRCLE,
                                        (float) (-1 * Math.log10(((double) keys[j] + 1)
                                                                 / pvals[i].length)),
                                        Math.min(maxValue, (float) (-1 * Math.log10(pvals[i][j]))),
                                        (byte) 6, (byte) (pvals.length == 1 ? 0 : i + 2), (byte) 0);
        } else {
          points[count] = new PlotPoint(keys[j] + "", PointType.FILLED_CIRCLE,
                                        (float) (((double) keys[j] + 1) / pvals[i].length),
                                        (float) pvals[i][j], (byte) 6,
                                        (byte) (pvals.length == 1 ? 0 : i + 2), (byte) 0);
        }
        count++;
      }
    }

  }

  public static void main(String[] args) {
    QQPlot.main(new String[] {});
  }
}
