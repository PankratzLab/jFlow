package one.ben.fcs.gating;

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
            this.min = min;
            this.max = max;
        }
        
        float min, max;

        public float getMin() {
            return min;
        }
        
        public float getMax() {
            return max;
        }
    }
    
    
}

