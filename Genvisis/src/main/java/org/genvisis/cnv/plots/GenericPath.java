package org.genvisis.cnv.plots;

import java.awt.geom.Path2D;
import java.util.Arrays;

public class GenericPath {

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + color;
    result = prime * result + (editable ? 1231 : 1237);
    result = prime * result + (fill ? 1231 : 1237);
    result = prime * result + fillColor;
    result = prime * result + Arrays.hashCode(foci);
    result = prime * result + layer;
    result = prime * result + ((myPath == null) ? 0 : myPath.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    GenericPath other = (GenericPath) obj;
    if (color != other.color) {
      return false;
    }
    if (editable != other.editable) {
      return false;
    }
    if (fill != other.fill) {
      return false;
    }
    if (fillColor != other.fillColor) {
      return false;
    }
    if (!Arrays.deepEquals(foci, other.foci)) {
      return false;
    }
    if (layer != other.layer) {
      return false;
    }
    if (myPath == null) {
      if (other.myPath != null) {
        return false;
      }
    } else if (!myPath.equals(other.myPath)) {
      return false;
    }
    return true;
  }

  private final String label;
  private final Path2D myPath;
  private final byte color;
  private final byte fillColor;
  private final byte layer;
  private final boolean fill;
  private boolean editable;
  private double[][] foci;

  public GenericPath(String label, Path2D path, byte color, byte fillColor, byte layer,
                     boolean fill, boolean editable) {
    this.label = label;
    myPath = path;
    this.color = color;
    this.fillColor = fillColor;
    this.layer = layer;
    this.fill = fill;
    this.editable = editable;
  }

  // public Ellipse2D getPath() {
  // return myShape;
  // }

  public byte getColor() {
    return color;
  }

  public byte getFillColor() {
    return fillColor;
  }

  public byte getLayer() {
    return layer;
  }

  public boolean getFill() {
    return fill;
  }

  public boolean getEditable() {
    return editable;
  }

  public void setEditable(boolean b) {
    editable = b;
  }

  public String getLabel() {
    return label;
  }

  public Path2D getPath() {
    return myPath;
  }



}
