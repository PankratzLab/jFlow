package one.ben.fcs.gating;

import java.awt.geom.Path2D;
import java.util.ArrayList;
import java.util.Random;

import common.CmdLineProcess.OUTPUT_Mode;
import one.ben.fcs.FCSDataLoader;
import one.ben.fcs.FCSDataLoader.DATA_SET;
import one.ben.fcs.gating.GateDimension.RectangleGateDimension;

public abstract class Gate {
    protected String popName;
    protected String id;
    protected String parentID;
    protected Gate parentGate;
    protected ArrayList<Gate> children = new ArrayList<Gate>();
    protected ArrayList<GateDimension> dimensions = new ArrayList<GateDimension>();
//    protected HashMap<String, boolean[]> gatingCache = new HashMap<String, boolean[]>();
    
    static final Random rand = new Random();
    
    private static String generateID() {
        return "ID" + (rand.nextInt(35000000) + 35000000);
    }
    
    public Gate(Gate parentGate2) {
        this(parentGate2, "PopulationName", generateID());
    }

    public Gate(Gate parentGate2, String popName, String id) {
        if (parentGate2 != null) {
            this.parentGate = parentGate2;
            this.parentID = parentGate2.id;
        }
        this.id = id;
        this.popName = popName;
    }

    public Gate getParentGate() {
        return parentGate;
    }
    
    protected boolean[] parentGating = null;
    
    public boolean[] getParentGating(FCSDataLoader dataLoader) {
//        if (parentGating == null) {
//            if (this.parentGate == null) {
//                return null;
//            } else {
//                return parentGating = this.parentGate.gate(dataLoader);
//            }
//        }
//        return parentGating;
        return this.parentGate == null ? null : this.parentGate.gate(dataLoader);
    }
    
//    public void clearCache() {
//        gatingCache.clear();
//    }
    
    public String getID() {
        return id;
    }
    
    public String getFullNameAndGatingPath() {
        StringBuilder full = new StringBuilder();
        
        if (!"".equals(getName())) {
            full.append(getName());
            full.append(" (");
            ArrayList<GateDimension> dims = getDimensions();
            for (int i = 0; i < dims.size(); i++) {
                full.append(dims.get(i).paramName);
                if (i < dims.size() - 1) {
                    full.append(" v ");
                }
            }
            full.append(")");
        }
        
        if (this.parentGate != null) {
            String p = this.parentGate.getFullNameAndGatingPath();
            full.insert(0, " / ");
            full.insert(0, p);
        }
        
        return full.toString();
    }
    
    public String getName() {
        return popName;
    }

    public ArrayList<Gate> getChildGates() {
        return children;
    }
    
    public ArrayList<GateDimension> getDimensions() {
        return dimensions;
    }
    
    public void addDimension(GateDimension gd) {
        this.dimensions.add(gd);
    }
    
    public abstract boolean[] gate(FCSDataLoader dataLoader);
    
    public static class RectangleGate extends Gate {
        
        public RectangleGate(Gate parentGate) {
            super(parentGate);
        }
        
        public RectangleGate(Gate parentGate, String popName, String id) {
            super(parentGate, popName, id);
        }
        
        public RectangleGateDimension getDimension(String param) {
            for (GateDimension gd : this.dimensions) {
                if (gd.paramName.equals(param)) {
                    return (RectangleGateDimension) gd;
                }
            }
            return null;
        }
        
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
//            if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
//                return gatingCache.get(dataLoader.getLoadedFile());
//            }
            boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] : (this.parentGating = this.parentGate.gate(dataLoader));
            boolean[][] paramIncludes = new boolean[dimensions.size()][dataLoader.getCount()];
            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                RectangleGateDimension rgd = (RectangleGateDimension) dimensions.get(p);
                if (!dataLoader.getAllDisplayableNames(DATA_SET.ALL).contains(rgd.paramName)) {
                    return null;
                }
                float[] paramData = dataLoader.getData(rgd.paramName, true);
                for (int i = 0; i < dataLoader.getCount(); i++) {
                    paramIncludes[p][i] = rgd.min <= paramData[i] && rgd.max > paramData[i]; 
                }
            }
            for (int i = 0; i < dataLoader.getCount(); i++) {
                boolean include = true;
                if (includes == null || (this.parentGate != null && !includes[i])) {
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
//            gatingCache.put(dataLoader.getLoadedFile(), includes);
            return includes;
        }
        
    }
    
    public static class PolygonGate extends Gate {
        ArrayList<Float> verticesX = new ArrayList<Float>(); 
        ArrayList<Float> verticesY = new ArrayList<Float>(); 
        Path2D myPath;
        
        public PolygonGate(Gate parentGate, String name, String id) {
            super(parentGate, name, id);
        }
        
        public PolygonGate(Gate parentGate) {
            super(parentGate);
        }
        
        @Override
        public void addDimension(GateDimension gd) {
            if (this.dimensions.size() == 2) {
                System.err.println("Error - cannot add more than two dimensions to a PolygonGate");
                return;
            }
            this.dimensions.add(gd);
        }
        
        public void setPath(Path2D pth) {
            this.myPath = pth;
//            clearCache();
        }
        
        public Path2D getPath() {
            if (myPath == null) {
                myPath = constructPath();
            }
            return myPath;
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
//            if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
//                return gatingCache.get(dataLoader.getLoadedFile());
//            }
            boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] : (this.parentGating = this.parentGate.gate(dataLoader));
            if (myPath == null) {
                myPath = constructPath();
            }
            float[][] paramData = new float[dimensions.size()][];
            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                GateDimension gd = dimensions.get(p);
                if (!dataLoader.getAllDisplayableNames(DATA_SET.ALL).contains(gd.paramName)) {
                    return null;
                }
                paramData[p] = dataLoader.getData(gd.paramName, true);
            }
            for (int i = 0; i < dataLoader.getCount(); i++) {
                if (includes == null || (this.parentGate != null && !includes[i])) {
                    continue;
                }
                includes[i] = myPath.contains(paramData[0][i], paramData[1][i]);
            }
//            gatingCache.put(dataLoader.getLoadedFile(), includes);
            return includes;
        }
        
    }
    
    public static class EllipsoidGate extends Gate {
        public EllipsoidGate() {
            super(null);
            throw new UnsupportedOperationException();
        }
        
//        public double[][] foci;
//        public double[][] edges;
//        
//        public Shape createEllipse(double[][] points) {
//            double minX = Math.min(Math.min(points[0][0], points[1][0]), Math.min(points[2][0], points[3][0]));
//            double minY = Math.min(Math.min(points[0][1], points[1][1]), Math.min(points[2][1], points[3][1]));
//            double maxY = Math.max(Math.max(points[0][1], points[1][1]), Math.max(points[2][1], points[3][1]));
//
//            double width = Math.sqrt((points[0][0] - points[1][0]) * (points[0][0] - points[1][0]) + (points[0][1] - points[1][1]) * (points[0][1] - points[1][1]));
//            double height = Math.sqrt((points[3][0] - points[2][0]) * (points[3][0] - points[2][0]) + (points[3][1] - points[2][1]) * (points[3][1] - points[2][1]));
//
//            double yD, xD;
//            yD = points[1][1] - points[0][1];
//            xD = points[1][0] - points[0][0];
//            double radAng = Math.atan2(yD, xD);
//            AffineTransform trans = AffineTransform.getRotateInstance(radAng, minX + width / 2, maxY + height / 2);
//            Ellipse2D ell = new Ellipse2D.Double(minX, maxY - height, width, height);
//            Shape rotEll = trans.createTransformedShape(ell);
//            return rotEll;
//        }
//        
//        Shape myShape; 
        
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            throw new UnsupportedOperationException();
//            boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] : this.parentGate.gate(dataLoader);
//            if (myShape == null) {
//                myShape = createEllipse(this.edges);
//            }
//            float[][] paramData = new float[dimensions.size()][];
//            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
//                GateDimension gd = dimensions.get(p);
//                paramData[p] = dataLoader.getData(gd.paramName, true);
//            }
//            for (int i = 0; i < dataLoader.getCount(); i++) {
//                if (this.parentGate != null && !includes[i]) {
//                    continue;
//                }
//                includes[i] = myShape.contains(paramData[0][i], paramData[1][i]);
//            }
//            
//            return includes;
        }
    }
    
//    public static class EllipseCalc {
//        
//        public static double[][] calcAll(double[] xy1, double[] xy2, double[] xy3, double[] xy4) {
//            double m1 = calcSlope(xy1[0], xy1[1], xy2[0], xy2[1]);
//            double m2 = calcSlope(xy3[0], xy3[1], xy4[0], xy4[1]);
//            
//            double[] newPt1, newPt2, newPt3, newPt4;
//            
//            newPt1 = calc(m1, m2, xy1[0], xy1[1], xy3[0], xy3[1]);
//            newPt2 = calc(m1, m2, xy1[0], xy1[1], xy4[0], xy4[1]);
//            newPt3 = calc(m1, m2, xy2[0], xy2[1], xy3[0], xy3[1]);
//            newPt4 = calc(m1, m2, xy2[0], xy2[1], xy4[0], xy4[1]);
//            
//            return new double[][]{newPt1, newPt2, newPt3, newPt4};
//        }
//        
//        private static double calcSlope(double x1, double y1, double x2, double y2) {
//            return (y2 - y1) / (x2 - x1);
//        }
//        
//        private static double[] calc(double m1, double m2, double x1, double y1, double x3, double y3) {
//            double xNew = ((y3 - m1 * x3) - (y1 - m2 * x1)) / (m2 - m1);
//            double yNew = (m2 * (y3 - m1 * x3) - m1 * (y1 - m2 * x1)) / (m2 - m1);
//            return new double[]{xNew, yNew};
//        }
//        
//    }
    
    public static class QuadrantGate extends Gate {
        public QuadrantGate() {
            super(null);
            throw new UnsupportedOperationException();
        }
        // UNSUPPORTED
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            throw new UnsupportedOperationException();
        }
    }
    
    public static class BooleanGate extends Gate {
        public BooleanGate() {
            super(null);
            throw new UnsupportedOperationException();
        }
        // UNSUPPORTED
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            throw new UnsupportedOperationException();
        }
    }
    
}
