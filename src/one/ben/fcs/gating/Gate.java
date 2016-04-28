package one.ben.fcs.gating;

import java.awt.geom.Path2D;
import java.util.ArrayList;

import one.ben.fcs.FCSDataLoader;
import one.ben.fcs.gating.GateDimension.RectangleGateDimension;

public abstract class Gate {
    String id;
    String parentID;
    Gate parentGate;
    protected ArrayList<Gate> children = new ArrayList<Gate>();
    protected ArrayList<GateDimension> dimensions = new ArrayList<GateDimension>();

    public void addDimension(GateDimension gd) {
        this.dimensions.add(gd);
    }
    
    public abstract boolean[] gate(FCSDataLoader dataLoader);
    
    public static class RectangleGate extends Gate {

        @Override
        public void addDimension(GateDimension gd) {
            if (!(gd instanceof RectangleGateDimension)) {
                System.err.println("Error - only RectangleGateDimensions can be added to a RectangleGate");
                return;
            }
            this.dimensions.add(gd);
        }
        
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] : this.parentGate.gate(dataLoader);
            boolean[][] paramIncludes = new boolean[dimensions.size()][dataLoader.getCount()];
            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                RectangleGateDimension rgd = (RectangleGateDimension) dimensions.get(p);
                float[] paramData = dataLoader.getData(rgd.paramName, true);
                for (int i = 0; i < dataLoader.getCount(); i++) {
                    paramIncludes[p][i] = rgd.min <= paramData[i] && rgd.max > paramData[i]; 
                }
            }
            for (int i = 0; i < dataLoader.getCount(); i++) {
                boolean include = true;
                if (this.parentGate != null && !includes[i]) {
                    continue;
                }
                for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                    if (!paramIncludes[p][i]) {
                        include = false;
                        break;
                    }
                }
                includes[i] = include;
            }
            return includes;
        }
        
    }
    
    public static class PolygonGate extends Gate {
        ArrayList<Float> verticesX = new ArrayList<Float>(); 
        ArrayList<Float> verticesY = new ArrayList<Float>(); 
        Path2D myPath;

        @Override
        public void addDimension(GateDimension gd) {
            if (this.dimensions.size() == 2) {
                System.err.println("Error - cannot add more than two dimensions to a PolygonGate");
                return;
            }
            this.dimensions.add(gd);
        }
        
        private Path2D constructPath() {
            Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
            path.moveTo(verticesX.get(0), verticesY.get(0));
            for(int i = 1; i < verticesX.size(); ++i) {
               path.lineTo(verticesX.get(i), verticesY.get(i));
            }
            path.closePath();
            return path;
        }
        
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] : this.parentGate.gate(dataLoader);
            if (myPath == null) {
                myPath = constructPath();
            }
            float[][] paramData = new float[dimensions.size()][];
            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                GateDimension gd = dimensions.get(p);
                paramData[p] = dataLoader.getData(gd.paramName, true);
            }
            for (int i = 0; i < dataLoader.getCount(); i++) {
                if (this.parentGate != null && !includes[i]) {
                    continue;
                }
                includes[i] = myPath.contains(paramData[0][i], paramData[1][i]);
            }
            
            return includes;
        }
        
    }
    
    public static class EllipsoidGate extends Gate {
        public double[][] foci;
        public double[][] edges;
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            throw new UnsupportedOperationException();
        }
    }
    
    public static class QuadrantGate extends Gate {
        // UNSUPPORTED
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            throw new UnsupportedOperationException();
        }
    }
    
    public static class BooleanGate extends Gate {
        // UNSUPPORTED
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            throw new UnsupportedOperationException();
        }
    }
    
}
