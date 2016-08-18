package org.genvisis.mining;

public class Point {
  public double value;

  public Point(double value) {
    this.value = value;
  }

  @Override
  public Point clone() {
    return new Point(value);
  }
}
