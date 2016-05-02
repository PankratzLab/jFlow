package one.ben.fcs.gating;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import cnv.plots.GenericPath;
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
        
        public void setPath(Path2D pth) {
            this.myPath = pth;
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
