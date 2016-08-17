package org.genvisis.cnv.plots;

import java.io.Serializable;

import org.genvisis.common.SerializedFiles;

public class PlotPointSet implements Serializable {
  public static final long serialVersionUID = 1L;

  public static PlotPointSet load(String filename, boolean jar) {
    return (PlotPointSet) SerializedFiles.readSerial(filename, jar, true);
  }

  private final PlotPoint[] points;

  public PlotPointSet(PlotPoint[] points) {
    this.points = points;
  }

  public PlotPoint[] getPlotPoints() {
    return points;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }
}
