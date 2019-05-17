package org.genvisis.one.ben;

import java.awt.BorderLayout;
import java.awt.Color;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.WindowConstants;
import org.genvisis.cnv.plots.AbstractPanel;
import org.genvisis.cnv.plots.PlotPoint;
import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.pankratzlab.common.PSF.Colors;

public class AlleleCountGenotyper {

  JFrame frame = new JFrame();

  public static void main(String[] args) {
    AlleleCountGenotyper acgt = new AlleleCountGenotyper();
    acgt.frame.setSize(600, 600);
    acgt.frame.setVisible(true);
  }

  JSlider sliderAMax;
  JSlider sliderBMax;

  public AlleleCountGenotyper() {
    frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
    frame.setLayout(new BorderLayout());
    frame.add(new AbstractPanel() {

      {
        setColorScheme(new Color[] {Color.BLACK, Color.RED, Colors.VIOLETS.ORCHID, Color.BLUE});
      }

      @Override
      public void highlightPoints() {

      }

      @Override
      public void generatePoints() {
        int aMax = sliderAMax.getValue();
        int bMax = sliderBMax.getValue();
        points = new PlotPoint[aMax * bMax];
        for (int a = 0; a < aMax; a++) {
          for (int b = 0; b < bMax; b++) {
            int ind = a * bMax + b;
            byte color = compute(a, b);
            points[ind] = new PlotPoint("", PointType.FILLED_SQUARE, a, b, (byte) 10, color,
                                        (byte) 1);
          }
        }
      }

      @Override
      public void assignAxisLabels() {
        xAxisLabel = "A";
        yAxisLabel = "B";
      }
    }, BorderLayout.CENTER);

    JPanel comp = new JPanel();
    frame.add(comp, BorderLayout.SOUTH);
    sliderAMax = new JSlider(JSlider.HORIZONTAL, 0, 100, 10);
    sliderAMax.setPaintLabels(true);
    comp.add(sliderAMax);
    sliderBMax = new JSlider(JSlider.HORIZONTAL, 0, 100, 10);
    sliderBMax.setPaintLabels(true);
    comp.add(sliderBMax);
  }

  byte compute(int a, int b) {
    if (b == 0) return 1; // pure homoA
    if (a == b) return 2; // pure hetero
    if (a == 0) return 3; // pure homoB

    double a1 = a / (double) (a + b);
    double b1 = b / (double) (a + b);

    if (a1 > .85) return 1;
    if (b1 > .85) return 3;

    return 2;
  }

}
