package org.genvisis.cnv.var;

public class IndiPoints {
  private final String id;
  private final double[][] datapoints;

  public IndiPoints(String id, double[][] datapoints) {
    this.id = id;
    this.datapoints = datapoints;
  }

  public double[][] getDatapoints() {
    return datapoints;
  }

  public String getId() {
    return id;
  }
}
