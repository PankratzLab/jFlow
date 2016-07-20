package org.genvisis.one.ben.fcs.gating;

public class GateDimension {
    String paramName;
    
    public GateDimension(String param) {
        this.paramName = param;
    }
    
    public static class RectangleGateDimension extends GateDimension {
        public RectangleGateDimension(String param) {
            super(param);
        }
        
        public RectangleGateDimension(String param, float min, float max) {
            super(param);
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
        }

        public void setMax(float max2) {
            this.max = max2;
        }
    }

    public String getParam() {
        return paramName;
    }
    
    
}

