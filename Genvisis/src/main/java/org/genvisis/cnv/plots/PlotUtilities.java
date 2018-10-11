package org.genvisis.cnv.plots;

import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.utils.gwas.windows.HitWindows;

public final class PlotUtilities {

  private PlotUtilities() {}

  public static void runHitWindows(String inputFile) {
    HitWindows.main(new String[] {"file=" + inputFile,
                                  "out=" + ext.rootOf(inputFile, false) + ".hits"});
  }

  public static void createManPlotScreenshot(String inputFile) {
    ManhattanPlot mp = new ManhattanPlot(null);
    mp.loadFileAuto(inputFile);
    mp.waitForData();
    mp.getManPan().setSize(800, 600);
    mp.screenshot(ext.rootOf(inputFile, false) + "_manPlot.png");
  }

  public static void createQQPlotScreenshot(String inputFile) {
    QQPlot qqPlot = QQPlot.loadPvals(new String[] {inputFile}, "Q-Q Plot", false, true, false, -1,
                                     -1, false, Float.MAX_VALUE, new Logger());
    qqPlot.screenCap(ext.rootOf(inputFile, false) + "_qqPlot.png");
  }

  public static void createAFPlotScreenshot(String inputFile) {
    AFPlot afPlot = new AFPlot(null);
    afPlot.loadFromFile(inputFile, null);
    afPlot.waitForData();
    afPlot.screenshot(ext.rootOf(inputFile, false) + "_afPlot.png");
  }

}
