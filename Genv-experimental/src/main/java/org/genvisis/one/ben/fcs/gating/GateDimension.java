package org.genvisis.one.ben.fcs.gating;

import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;

public class GateDimension {
    String paramName;
    Gate owner;
    
    public GateDimension(Gate gate, String param) {
        this.paramName = param;
        this.owner = gate;
    }
    
    public static class RectangleGateDimension extends GateDimension {
        public RectangleGateDimension(RectangleGate gate, String param) {
            super(gate, param);
        }
        
        public RectangleGateDimension(RectangleGate gate, String param, float min, float max) {
            super(gate, param);
            this.min = Math.min(min, max);
            this.max = Math.max(min, max);
        }
        
        float min, max;

        public float getMin() {
            return min;
        }
        
        public float getMax() {
            return max;
        }

        public void setMin(float min2) {
            this.min = min2;
            this.owner.parentGating = null;
        }

        public void setMax(float max2) {
            this.max = max2;
            this.owner.parentGating = null;
        }
    }

    public String getParam() {
        return paramName;
    }
    
    
}

