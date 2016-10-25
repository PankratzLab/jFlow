package org.genvisis.one.ben.fcs.gating;

import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;

public class GateDimension {
  String paramName;
  Gate owner;
  AXIS_SCALE scale;

  public GateDimension(Gate gate, String param, AXIS_SCALE scale) {
    this.paramName = param;
    this.owner = gate;
    this.scale = scale;
  }

  public static class RectangleGateDimension extends GateDimension {
    public RectangleGateDimension(RectangleGate gate, String param, AXIS_SCALE scale) {
      super(gate, param, scale);
    }

    public RectangleGateDimension(RectangleGate gate, String param, AXIS_SCALE scale, float min,
        float max) {
      super(gate, param, scale);
      this.min = Math.min(min, max);
      this.max = Math.max(min, max);
    }

    private float min, max;

    public float getMin() {
      return min;
    }

    public float getMax() {
      return max;
    }

    public void setMin(float min2) {
      min = min2;
      owner.parentGating = null;
    }

    public void setMax(float max2) {
      max = max2;
      owner.parentGating = null;
    }
  }

  public String getParam() {
    return paramName;
  }


}
