package org.genvisis.one.ben.fcs.gating;

import java.awt.Rectangle;
import java.awt.geom.Area;
import java.awt.geom.IllegalPathStateException;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import org.apache.poi.sl.usermodel.Shape;
import org.genvisis.common.Array;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;


public abstract class Gate {
    
    static class GateUtils {
        
        /**
        * Tests if the specified coordinates are inside the closed
        * boundary of the specified {@link PathIterator}.
        * <p>
        * This method provides a basic facility for implementors of
        * the {@link Shape} interface to implement support for the
        * {@link Shape#contains(double, double)} method.
        *
        * @param pi the specified {@code PathIterator}
        * @param x the specified X coordinate
        * @param y the specified Y coordinate
        * @return {@code true} if the specified coordinates are inside the
        *         specified {@code PathIterator}; {@code false} otherwise
        * @since 1.6
        */
       public static boolean contains(PathIterator pi, double x, double y) {
           if (x * 0.0 + y * 0.0 == 0.0) {
               /* N * 0.0 is 0.0 only if N is finite.
                * Here we know that both x and y are finite.
                */
               int mask = (pi.getWindingRule() == Path2D.WIND_NON_ZERO ? -1 : 1);
               int cross = pointCrossingsForPath(pi, x, y);
               return ((cross & mask) != 0);
           } else {
               /* Either x or y was infinite or NaN.
                * A NaN always produces a negative response to any test
                * and Infinity values cannot be "inside" any path so
                * they should return false as well.
                */
               return false;
           }
       }

        /**
         * Calculates the number of times the given path
         * crosses the ray extending to the right from (px,py).
         * If the point lies on a part of the path,
         * then no crossings are counted for that intersection.
         * +1 is added for each crossing where the Y coordinate is increasing
         * -1 is added for each crossing where the Y coordinate is decreasing
         * The return value is the sum of all crossings for every segment in
         * the path.
         * The path must start with a SEG_MOVETO, otherwise an exception is
         * thrown.
         * The caller must check p[xy] for NaN values.
         * The caller may also reject infinite p[xy] values as well.
         */
        public static int pointCrossingsForPath(PathIterator pi,
                                                double px, double py)
        {
            if (pi.isDone()) {
                return 0;
            }
            double coords[] = new double[6];
            if (pi.currentSegment(coords) != PathIterator.SEG_MOVETO) {
                throw new IllegalPathStateException("missing initial moveto "+
                                                    "in path definition");
            }
            pi.next();
            double movx = coords[0];
            double movy = coords[1];
            double curx = movx;
            double cury = movy;
            double endx, endy;
            int crossings = 0;
            while (!pi.isDone()) {
                switch (pi.currentSegment(coords)) {
                case PathIterator.SEG_MOVETO:
                    if (cury != movy) {
                        crossings += pointCrossingsForLine(px, py,
                                                           curx, cury,
                                                           movx, movy);
                    }
                    movx = curx = coords[0];
                    movy = cury = coords[1];
                    break;
                case PathIterator.SEG_LINETO:
                    endx = coords[0];
                    endy = coords[1];
                    crossings += pointCrossingsForLine(px, py,
                                                       curx, cury,
                                                       endx, endy);
                    curx = endx;
                    cury = endy;
                    break;
                case PathIterator.SEG_QUADTO:
                    endx = coords[2];
                    endy = coords[3];
                    crossings += pointCrossingsForQuad(px, py,
                                                       curx, cury,
                                                       coords[0], coords[1],
                                                       endx, endy, 0);
                    curx = endx;
                    cury = endy;
                    break;
                case PathIterator.SEG_CUBICTO:
                    endx = coords[4];
                    endy = coords[5];
                    crossings += pointCrossingsForCubic(px, py,
                                                        curx, cury,
                                                        coords[0], coords[1],
                                                        coords[2], coords[3],
                                                        endx, endy, 0);
                    curx = endx;
                    cury = endy;
                    break;
                case PathIterator.SEG_CLOSE:
                    if (cury != movy) {
                        crossings += pointCrossingsForLine(px, py,
                                                           curx, cury,
                                                           movx, movy);
                    }
                    curx = movx;
                    cury = movy;
                    break;
                }
                pi.next();
            }
            if (cury != movy) {
                crossings += pointCrossingsForLine(px, py,
                                                   curx, cury,
                                                   movx, movy);
            }
            return crossings;
        }

        /**
         * Calculates the number of times the line from (x0,y0) to (x1,y1)
         * crosses the ray extending to the right from (px,py).
         * If the point lies on the line, then no crossings are recorded.
         * +1 is returned for a crossing where the Y coordinate is increasing
         * -1 is returned for a crossing where the Y coordinate is decreasing
         */
        public static int pointCrossingsForLine(double px, double py,
                                                double x0, double y0,
                                                double x1, double y1)
        {
            if (py <  y0 && py <  y1) return 0;
            if (py >= y0 && py >= y1) return 0;
            // assert(y0 != y1);
            if (px >= x0 && px >= x1) return 0;
            if (px <  x0 && px <  x1) return (y0 < y1) ? 1 : -1;
            double xintercept = x0 + (py - y0) * (x1 - x0) / (y1 - y0);
            if (px >= xintercept) return 0;
            return (y0 < y1) ? 1 : -1;
        }
        

        /**
         * Calculates the number of times the quad from (x0,y0) to (x1,y1)
         * crosses the ray extending to the right from (px,py).
         * If the point lies on a part of the curve,
         * then no crossings are counted for that intersection.
         * the level parameter should be 0 at the top-level call and will count
         * up for each recursion level to prevent infinite recursion
         * +1 is added for each crossing where the Y coordinate is increasing
         * -1 is added for each crossing where the Y coordinate is decreasing
         */
        public static int pointCrossingsForQuad(double px, double py,
                                                double x0, double y0,
                                                double xc, double yc,
                                                double x1, double y1, int level)
        {
            if (py <  y0 && py <  yc && py <  y1) return 0;
            if (py >= y0 && py >= yc && py >= y1) return 0;
            // Note y0 could equal y1...
            if (px >= x0 && px >= xc && px >= x1) return 0;
            if (px <  x0 && px <  xc && px <  x1) {
                if (py >= y0) {
                    if (py < y1) return 1;
                } else {
                    // py < y0
                    if (py >= y1) return -1;
                }
                // py outside of y01 range, and/or y0==y1
                return 0;
            }
            // double precision only has 52 bits of mantissa
            if (level > 52) return pointCrossingsForLine(px, py, x0, y0, x1, y1);
            double x0c = (x0 + xc) / 2;
            double y0c = (y0 + yc) / 2;
            double xc1 = (xc + x1) / 2;
            double yc1 = (yc + y1) / 2;
            xc = (x0c + xc1) / 2;
            yc = (y0c + yc1) / 2;
            if (Double.isNaN(xc) || Double.isNaN(yc)) {
                // [xy]c are NaN if any of [xy]0c or [xy]c1 are NaN
                // [xy]0c or [xy]c1 are NaN if any of [xy][0c1] are NaN
                // These values are also NaN if opposing infinities are added
                return 0;
            }
            return (pointCrossingsForQuad(px, py,
                                          x0, y0, x0c, y0c, xc, yc,
                                          level+1) +
                    pointCrossingsForQuad(px, py,
                                          xc, yc, xc1, yc1, x1, y1,
                                          level+1));
        }

        /**
         * Calculates the number of times the cubic from (x0,y0) to (x1,y1)
         * crosses the ray extending to the right from (px,py).
         * If the point lies on a part of the curve,
         * then no crossings are counted for that intersection.
         * the level parameter should be 0 at the top-level call and will count
         * up for each recursion level to prevent infinite recursion
         * +1 is added for each crossing where the Y coordinate is increasing
         * -1 is added for each crossing where the Y coordinate is decreasing
         */
        public static int pointCrossingsForCubic(double px, double py,
                                                 double x0, double y0,
                                                 double xc0, double yc0,
                                                 double xc1, double yc1,
                                                 double x1, double y1, int level)
        {
            if (py <  y0 && py <  yc0 && py <  yc1 && py <  y1) return 0;
            if (py >= y0 && py >= yc0 && py >= yc1 && py >= y1) return 0;
            // Note y0 could equal yc0...
            if (px >= x0 && px >= xc0 && px >= xc1 && px >= x1) return 0;
            if (px <  x0 && px <  xc0 && px <  xc1 && px <  x1) {
                if (py >= y0) {
                    if (py < y1) return 1;
                } else {
                    // py < y0
                    if (py >= y1) return -1;
                }
                // py outside of y01 range, and/or y0==yc0
                return 0;
            }
            // double precision only has 52 bits of mantissa
            if (level > 52) return pointCrossingsForLine(px, py, x0, y0, x1, y1);
            double xmid = (xc0 + xc1) / 2;
            double ymid = (yc0 + yc1) / 2;
            xc0 = (x0 + xc0) / 2;
            yc0 = (y0 + yc0) / 2;
            xc1 = (xc1 + x1) / 2;
            yc1 = (yc1 + y1) / 2;
            double xc0m = (xc0 + xmid) / 2;
            double yc0m = (yc0 + ymid) / 2;
            double xmc1 = (xmid + xc1) / 2;
            double ymc1 = (ymid + yc1) / 2;
            xmid = (xc0m + xmc1) / 2;
            ymid = (yc0m + ymc1) / 2;
            if (Double.isNaN(xmid) || Double.isNaN(ymid)) {
                // [xy]mid are NaN if any of [xy]c0m or [xy]mc1 are NaN
                // [xy]c0m or [xy]mc1 are NaN if any of [xy][c][01] are NaN
                // These values are also NaN if opposing infinities are added
                return 0;
            }
            return (pointCrossingsForCubic(px, py,
                                           x0, y0, xc0, yc0,
                                           xc0m, yc0m, xmid, ymid, level+1) +
                    pointCrossingsForCubic(px, py,
                                           xmid, ymid, xmc1, ymc1,
                                           xc1, yc1, x1, y1, level+1));
        }
        
    }
    
    
    protected String popName;
    protected String id;
    protected String parentID;
    protected Gate parentGate;
    protected ArrayList<Gate> children = new ArrayList<Gate>();
    protected ArrayList<GateDimension> dimensions = new ArrayList<GateDimension>();
//    protected HashMap<String, boolean[]> gatingCache = new HashMap<String, boolean[]>();
    protected int displayLevel = 0;
    
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
        if (parentGating == null) {
            if (this.parentGate == null) {
                return parentGating = Array.booleanArray(dataLoader.getCount(), true);
            } else {
                return parentGating = this.parentGate.gate(dataLoader);
            }
        }
        return parentGating;
//        return this.parentGate == null ? null : this.parentGate.gate(dataLoader);
    }
    
//    public void clearCache() {
//        gatingCache.clear();
//    }
    
    public String getID() {
        return id;
    }
    
    public void setLevel(int lvl) {
        this.displayLevel = lvl;
    }

    public int getLevel() {
        return this.displayLevel;
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
            this.parentGating = null;
        }
        
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
//            if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
//                return gatingCache.get(dataLoader.getLoadedFile());
//            }
            boolean[] includes = this.parentGate == null ? Array.booleanArray(dataLoader.getCount(), true) : this.parentGate.gate(dataLoader);
            this.parentGating = includes == null ? null : Arrays.copyOf(includes, includes.length);
            boolean[][] paramIncludes = new boolean[dimensions.size()][dataLoader.getCount()];
            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                RectangleGateDimension rgd = (RectangleGateDimension) dimensions.get(p);
                if (!dataLoader.getAllDisplayableNames(DATA_SET.ALL).contains(rgd.paramName)) {
                    return null;
                }
                double[] paramData = dataLoader.getData(rgd.paramName, true);
                for (int i = 0; i < dataLoader.getCount(); i++) {
                    paramIncludes[p][i] = rgd.min <= paramData[i] && rgd.max >= paramData[i]; 
                }
            }
            if (includes == null) {
                return includes;
            }
            for (int i = 0; i < dataLoader.getCount(); i++) {
                boolean include = true;
                if (this.parentGate != null && !includes[i]) {
//                    includes[i] = false;
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
        private ArrayList<Double> verticesX = new ArrayList<Double>(); 
        private ArrayList<Double> verticesY = new ArrayList<Double>(); 
        private Path2D myPath;
        private int gateResolution;
        private boolean mimicFlowJo = false;
        
        public PolygonGate(Gate parentGate, String name, String id) {
            super(parentGate, name, id);
        }
        
        public PolygonGate(Gate parentGate) {
            super(parentGate);
        }
        
        public boolean getMimicsFlowJo() {
            return mimicFlowJo;
        }
        
        public void setShouldMimicFlowJoGating(boolean mimic) {
            this.mimicFlowJo = mimic;
            prepGating();
            this.parentGating = null;
        }
        
        @Override
        public void addDimension(GateDimension gd) {
            if (this.dimensions.size() == 2) {
                System.err.println("Error - cannot add more than two dimensions to a PolygonGate");
                return;
            }
            this.dimensions.add(gd);
            this.parentGating = null;
        }
        
        public void setPath(Path2D pth) {
            this.myPath = pth;
            this.parentGating = null;
            resetVertices();
            prepGating();
//            clearCache();
        }
        
        private void resetVertices() {
            double[] coords = new double[6];
            PathIterator pi = myPath.getPathIterator(null);
            while (!pi.isDone()) {
                pi.currentSegment(coords);
                addVertex(coords[0], coords[1]);
                pi.next(); // TODO may require a 'next()' to start?
            }
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
            for (int i = 1; i < verticesX.size(); ++i) {
               path.lineTo(verticesX.get(i), verticesY.get(i));
            }
            path.lineTo(verticesX.get(0), verticesY.get(0));
            path.closePath();
            return path;
        }
        
        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
//            if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
//                return gatingCache.get(dataLoader.getLoadedFile());
//            }
            boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] : this.parentGate.gate(dataLoader);
            this.parentGating = includes == null ? null : includes.clone();
            if (myPath == null) {
                myPath = constructPath();
            }
            double[][] paramData = new double[dimensions.size()][];
            for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
                GateDimension gd = dimensions.get(p);
                if (!dataLoader.getAllDisplayableNames(DATA_SET.ALL).contains(gd.paramName)) {
                    return null;
                }
                paramData[p] = dataLoader.getData(gd.paramName, true);
            }
            if (includes == null) {
                return includes;
            }
            for (int i = 0; i < dataLoader.getCount(); i++) {
                if (this.parentGate != null && !includes[i]) {
//                    includes[i] = false;
                    continue;
                }
                if (mimicFlowJo) {
                    for (Rectangle rect : myRects) {
                        if (rect.contains(paramData[0][i], paramData[1][i])) {
                            includes[i] = true;
                            break;
                        }
                    }
                } else {
                    if (myPath.contains(paramData[0][i], paramData[1][i])) {
                        includes[i] = true;
                    }
                }
                
            }            

//            gatingCache.put(dataLoader.getLoadedFile(), includes);
            return includes;
        }
        
        int range = 262144;
        int numBins = 256;
        int binStep = range / numBins; // 1024
        ArrayList<Rectangle> myRects = new ArrayList<Rectangle>();
        
        void prepGating() {
            ArrayList<Rectangle> vertexRects = new ArrayList<Rectangle>();
            ArrayList<Rectangle> rects = new ArrayList<Rectangle>();
            for (int i = 0; i < numBins; i++) {
                for (int j = 0; j < numBins; j++) {
                    rects.add(new Rectangle(i * binStep + binStep / 2, j * binStep + binStep / 2, binStep, binStep));
                }
            }
            for (Rectangle r : rects) {
                for (int i = 0; i < verticesX.size(); i++) {
                    double x, y;
                    x = verticesX.get(i);
                    y = verticesY.get(i);
                    if (r.contains(x, y)) {
                        vertexRects.add(r);
                    }
                }
            }
            
            Area a;
            Path2D path = new Path2D.Double();
            path.moveTo(vertexRects.get(0).getCenterX(), vertexRects.get(0).getCenterY());
            for (int i = 1; i < vertexRects.size(); i++) {
                path.lineTo(vertexRects.get(i).getCenterX(), vertexRects.get(i).getCenterY());
            }
            path.lineTo(vertexRects.get(0).getCenterX(), vertexRects.get(0).getCenterY());
            path.closePath();
            a = new Area(path);
            
            int index = 1;
            while (!a.isSingular()) {
                Collections.swap(vertexRects, 0, index);
                path = new Path2D.Double();
                path.moveTo(vertexRects.get(0).getCenterX(), vertexRects.get(0).getCenterY());
                for (int i = 1; i < vertexRects.size(); i++) {
                    path.lineTo(vertexRects.get(i).getCenterX(), vertexRects.get(i).getCenterY());
                }
                path.lineTo(vertexRects.get(0).getCenterX(), vertexRects.get(0).getCenterY());
                path.closePath();
                a = new Area(path);
                index = (index + 1) % vertexRects.size();
            }
            
            
            for (Rectangle rect : rects) {
                if (vertexRects.contains(rect) || path.contains(rect) || (path.intersects(rect) && path.contains(rect.getCenterX(), rect.getCenterY()))) {
                    myRects.add(rect);
                }
            }
            
            rects.clear();
            rects = null;
        }

        public void setGateResolution(int res) {
            this.gateResolution = res;
        }

        public void addVertex(Double fX, Double fY) {
            this.verticesX.add(fX);
            this.verticesY.add(fY);
            this.parentGating = null;
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