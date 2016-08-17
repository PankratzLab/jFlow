package org.genvisis.one.ben.fcs.gating;

import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import org.genvisis.common.Array;
import org.genvisis.common.Numbers;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;

public abstract class Gate {

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

    @Override
    public String getXMLTag() {
      return "gating:BooleanGate";
    }
  }

  public static class EllipsoidGate extends Gate {
    public EllipsoidGate() {
      super(null);
      throw new UnsupportedOperationException();
    }

    // public double[][] foci;
    // public double[][] edges;
    //
    // public Shape createEllipse(double[][] points) {
    // double minX = Math.min(Math.min(points[0][0], points[1][0]), Math.min(points[2][0],
    // points[3][0]));
    // double minY = Math.min(Math.min(points[0][1], points[1][1]), Math.min(points[2][1],
    // points[3][1]));
    // double maxY = Math.max(Math.max(points[0][1], points[1][1]), Math.max(points[2][1],
    // points[3][1]));
    //
    // double width = Math.sqrt((points[0][0] - points[1][0]) * (points[0][0] - points[1][0]) +
    // (points[0][1] - points[1][1]) * (points[0][1] - points[1][1]));
    // double height = Math.sqrt((points[3][0] - points[2][0]) * (points[3][0] - points[2][0]) +
    // (points[3][1] - points[2][1]) * (points[3][1] - points[2][1]));
    //
    // double yD, xD;
    // yD = points[1][1] - points[0][1];
    // xD = points[1][0] - points[0][0];
    // double radAng = Math.atan2(yD, xD);
    // AffineTransform trans = AffineTransform.getRotateInstance(radAng, minX + width / 2, maxY +
    // height / 2);
    // Ellipse2D ell = new Ellipse2D.Double(minX, maxY - height, width, height);
    // Shape rotEll = trans.createTransformedShape(ell);
    // return rotEll;
    // }
    //
    // Shape myShape;

    @Override
    public boolean[] gate(FCSDataLoader dataLoader) {
      throw new UnsupportedOperationException();
      // boolean[] includes = this.parentGate == null ? new boolean[dataLoader.getCount()] :
      // this.parentGate.gate(dataLoader);
      // if (myShape == null) {
      // myShape = createEllipse(this.edges);
      // }
      // float[][] paramData = new float[dimensions.size()][];
      // for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
      // GateDimension gd = dimensions.get(p);
      // paramData[p] = dataLoader.getData(gd.paramName, true);
      // }
      // for (int i = 0; i < dataLoader.getCount(); i++) {
      // if (this.parentGate != null && !includes[i]) {
      // continue;
      // }
      // includes[i] = myShape.contains(paramData[0][i], paramData[1][i]);
      // }
      //
      // return includes;
    }

    @Override
    public String getXMLTag() {
      return "gating:EllipsoidGate";
    }
  }

  public static class PolygonGate extends Gate {
    private final ArrayList<Double> verticesX = new ArrayList<Double>();
    private final ArrayList<Double> verticesY = new ArrayList<Double>();
    private Path2D myPath;
    private int gateResolution = 256; // default
    private boolean mimicFlowJo = false;

    int range = 262144; // TODO make this dynamic based on data loader param range (set when dims
                        // are set) - may have different ranges for each dim

    ArrayList<Rectangle> myRects = new ArrayList<Rectangle>();

    public PolygonGate(Gate parentGate) {
      this(parentGate, null, generateID(), false);
    }

    public PolygonGate(Gate parentGate, String name) {
      this(parentGate, name, generateID(), false);
    }

    public PolygonGate(Gate parentGate, String name, boolean flowjo) {
      this(parentGate, name, generateID(), flowjo);
    }

    public PolygonGate(Gate parentGate, String name, String id) {
      this(parentGate, name, id, false);
    }

    public PolygonGate(Gate parentGate, String name, String id, boolean flowjo) {
      super(parentGate, name, id);
      setShouldMimicFlowJoGating(flowjo);
    }

    @Override
    public void addDimension(GateDimension gd) {
      if (dimensions.size() == 2) {
        System.err.println("Error - cannot add more than two dimensions to a PolygonGate");
        return;
      }
      dimensions.add(gd);
      parentGating = null;
    }

    public void addVertex(Double fX, Double fY) {
      verticesX.add(fX);
      verticesY.add(fY);
      parentGating = null;
    }

    private Path2D constructPath() {
      Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
      path.moveTo(verticesX.get(0), verticesY.get(0));
      for (int i = 1; i < verticesX.size(); ++i) {
        path.lineTo(verticesX.get(i), verticesY.get(i));
      }
      path.closePath();
      return path;
    }

    @Override
    public boolean[] gate(FCSDataLoader dataLoader) {
      // if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
      // return gatingCache.get(dataLoader.getLoadedFile());
      // }
      boolean[] includes =
          parentGate == null ? new boolean[dataLoader.getCount()] : parentGate.gate(dataLoader);
      parentGating = parentGate == null ? Array.booleanArray(dataLoader.getCount(), true)
          : Arrays.copyOf(includes, includes.length);
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
        if (parentGate != null && !includes[i]) {
          continue;
        }
        boolean include = false;
        if (mimicFlowJo) {
          for (Rectangle rect : myRects) {
            if (rect.contains(paramData[0][i], paramData[1][i])) {
              include = true;
              break;
            }
          }
        } else {
          if (myPath.contains(paramData[0][i], paramData[1][i])) {
            include = true;
          }
        }
        includes[i] = include;
      }

      // gatingCache.put(dataLoader.getLoadedFile(), includes);
      return includes;
    }

    public boolean getMimicsFlowJo() {
      return mimicFlowJo;
    }

    public Path2D getPath() {
      if (myPath == null) {
        myPath = constructPath();
      }
      return myPath;
    }

    @Override
    public String getXMLTag() {
      return "gating:PolygonGate";
    }

    void prepGating() {
      ArrayList<Rectangle> vertexRects = new ArrayList<Rectangle>();
      ArrayList<Rectangle> rects = new ArrayList<Rectangle>();
      int binStep = range / gateResolution;
      for (int i = 0; i < gateResolution; i++) {
        for (int j = 0; j < gateResolution; j++) {
          rects.add(new Rectangle(i * binStep + binStep / 2, j * binStep + binStep / 2, binStep,
              binStep));
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


      double xSum = 0, ySum = 0;
      for (Rectangle r : vertexRects) {
        xSum += r.getCenterX();
        ySum += r.getCenterY();
      }
      final double x = xSum / vertexRects.size();
      final double y = ySum / vertexRects.size();

      // taken from:
      // https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
      Collections.sort(vertexRects, new Comparator<Rectangle>() {
        @Override
        public int compare(Rectangle o1, Rectangle o2) {
          if (o1.getCenterX() - x >= 0 && o2.getCenterX() - x < 0) {
            return -1;
          }
          if (o1.getCenterX() - x < 0 && o2.getCenterX() - x >= 0) {
            return 1;
          }
          if (o1.getCenterX() - x == 0 && o2.getCenterX() - x == 0) {
            if (o1.getCenterY() - y >= 0 || o2.getCenterY() - y >= 0) {
              return o1.getCenterY() > o2.getCenterY() ? -1 : 1;
            }
            return o2.getCenterY() > o1.getCenterY() ? -1 : 1;
          }

          // compute the cross product of vectors (center -> a) x (center -> b)
          double det = (o1.getCenterX() - x) * (o2.getCenterY() - y)
              - (o2.getCenterX() - x) * (o1.getCenterY() - y);
          if (det < 0) {
            return -1;
          }
          if (det > 0) {
            return 1;
          }

          // points a and b are on the same line from the center
          // check which point is closer to the center
          double d1 = (o1.getCenterX() - x) * (o1.getCenterX() - x)
              + (o1.getCenterY() - y) * (o1.getCenterY() - y);
          double d2 = (o2.getCenterX() - x) * (o2.getCenterX() - x)
              + (o2.getCenterY() - y) * (o2.getCenterY() - y);
          return d1 > d2 ? -1 : d1 < d2 ? 1 : 0;
        }
      });

      Path2D path = new Path2D.Double();
      path.moveTo(vertexRects.get(0).getCenterX(), vertexRects.get(0).getCenterY());
      for (int i = 1; i < vertexRects.size(); i++) {
        path.lineTo(vertexRects.get(i).getCenterX(), vertexRects.get(i).getCenterY());
      }
      path.closePath();

      myRects.clear();
      for (Rectangle rect : rects) {
        if (vertexRects.contains(rect) || path.contains(rect)
            || (path.intersects(rect) && path.contains(rect.getCenterX(), rect.getCenterY()))) {
          myRects.add(rect);
        }
      }

      rects.clear();
      rects = null;
    }

    private void resetVertices() {
      verticesX.clear();
      verticesY.clear();
      double[] coords = new double[6];
      PathIterator pi = myPath.getPathIterator(null);
      while (!pi.isDone()) {
        pi.currentSegment(coords);
        addVertex(coords[0], coords[1]);
        pi.next(); // TODO may require a 'next()' to start?
      }
    }

    public void setGateResolution(int res) {
      gateResolution = res;
    }

    public void setPath(Path2D pth) {
      myPath = pth;
      parentGating = null;
      resetVertices();
      if (mimicFlowJo) {
        prepGating();
      }
      // clearCache();
    }

    public void setShouldMimicFlowJoGating(boolean mimic) {
      mimicFlowJo = mimic;
      if (mimicFlowJo) {
        prepGating();
      }
      parentGating = null;
    }

    public void transform(AffineTransform at) {
      myPath.transform(at);
      resetVertices();
      if (mimicFlowJo) {
        prepGating();
      }
      parentGating = null;
    }


  }
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

    @Override
    public String getXMLTag() {
      return "gating:RectangleGate"; // TODO stored as four RectangleGate tags
    }
  }
  public static class RectangleGate extends Gate {

    public RectangleGate(Gate parentGate) {
      super(parentGate);
    }

    public RectangleGate(Gate parentGate, String popName) {
      super(parentGate, popName);
    }

    public RectangleGate(Gate parentGate, String popName, String id) {
      super(parentGate, popName, id);
    }

    @Override
    public void addDimension(GateDimension gd) {
      if (!(gd instanceof RectangleGateDimension)) {
        System.err.println("Error - only RectangleGateDimensions can be added to a RectangleGate");
        return;
      }
      dimensions.add(gd);
      parentGating = null;
    }

    @Override
    public boolean[] gate(FCSDataLoader dataLoader) {
      // if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
      // return gatingCache.get(dataLoader.getLoadedFile());
      // }
      boolean[] includes =
          parentGate == null ? new boolean[dataLoader.getCount()] : parentGate.gate(dataLoader);
      parentGating = parentGate == null ? Array.booleanArray(dataLoader.getCount(), true)
          : Arrays.copyOf(includes, includes.length);
      boolean[][] paramIncludes = new boolean[dimensions.size()][dataLoader.getCount()];
      for (int p = 0, pCount = dimensions.size(); p < pCount; p++) {
        RectangleGateDimension rgd = (RectangleGateDimension) dimensions.get(p);
        if (!dataLoader.getAllDisplayableNames(DATA_SET.ALL).contains(rgd.paramName)) {
          return null;
        }
        double[] paramData = dataLoader.getData(rgd.paramName, true);
        for (int i = 0; i < dataLoader.getCount(); i++) {
          // inclusive min, exclusive max - see gating-ml spec
          paramIncludes[p][i] = (!Numbers.isFinite(rgd.getMin()) || rgd.getMin() <= paramData[i])
              && (!Numbers.isFinite(rgd.getMax()) || rgd.getMax() > paramData[i]);
        }
      }
      if (includes == null) {
        return includes;
      }
      for (int i = 0; i < dataLoader.getCount(); i++) {
        boolean include = true;
        if (parentGate != null && !includes[i]) {
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
      // gatingCache.put(dataLoader.getLoadedFile(), includes);
      return includes;
    }

    public RectangleGateDimension getDimension(String param) {
      for (GateDimension gd : dimensions) {
        if (gd.paramName.equals(param)) {
          return (RectangleGateDimension) gd;
        }
      }
      return null;
    }

    @Override
    public String getXMLTag() {
      return "gating:RectangleGate";
    }

  }

  static final Random rand = new Random();

  private static String generateID() {
    return "ID" + (rand.nextInt(35000000) + 35000000);
  }

  protected int insideColorCode;
  protected String xmlTag;

  protected String popName;

  protected String id;

  protected String parentID;

  protected Gate parentGate;

  protected ArrayList<Gate> children = new ArrayList<Gate>();

  protected ArrayList<GateDimension> dimensions = new ArrayList<GateDimension>();

  // protected HashMap<String, boolean[]> gatingCache = new HashMap<String, boolean[]>();
  protected int displayLevel = 0;

  protected boolean[] parentGating = null;

  // public void clearCache() {
  // gatingCache.clear();
  // }

  public Gate(Gate parentGate2) {
    this(parentGate2, null, generateID());
  }

  public Gate(Gate parentGate2, String popName) {
    this(parentGate2, popName, generateID());
  }

  public Gate(Gate parentGate2, String popName, String id) {
    if (parentGate2 != null) {
      parentGate = parentGate2;
      parentID = parentGate2.id;
    }
    this.id = id;
    this.popName = popName;
  }

  public void addChildGate(Gate rg) {
    children.add(rg);
  }

  public void addDimension(GateDimension gd) {
    dimensions.add(gd);
  }

  public int countAllChildren() {
    int cnt = getChildGates().size();
    for (Gate g2 : getChildGates()) {
      cnt += g2.countAllChildren();
    }
    return cnt;
  }

  public abstract boolean[] gate(FCSDataLoader dataLoader);

  public ArrayList<Gate> getChildGates() {
    return children;
  }

  public int getColorCode() {
    return insideColorCode;
  }

  public ArrayList<GateDimension> getDimensions() {
    return dimensions;
  }

  public String getFullNameAndGatingPath() {
    StringBuilder full = new StringBuilder();

    if (!"".equals(getName())) {
      full.append(getName());
    }

    full.append(" (");
    ArrayList<GateDimension> dims = getDimensions();
    for (int i = 0; i < dims.size(); i++) {
      full.append(dims.get(i).paramName);
      if (i < dims.size() - 1) {
        full.append(" v ");
      }
    }
    full.append(")");

    if (parentGate != null) {
      String p = parentGate.getFullNameAndGatingPath();
      full.insert(0, " / ");
      full.insert(0, p);
    }

    return full.toString();
  }

  public String getID() {
    return id;
  }

  public int getLevel() {
    return displayLevel;
  }

  public String getName() {
    return popName;
  }

  public Gate getParentGate() {
    return parentGate;
  }

  public boolean[] getParentGating(FCSDataLoader dataLoader) {
    if (parentGating == null) {
      if (parentGate == null) {
        return parentGating = Array.booleanArray(dataLoader.getCount(), true);
      } else {
        return parentGating = parentGate.gate(dataLoader);
      }
    }
    return parentGating;
    // return this.parentGate == null ? null : this.parentGate.gate(dataLoader);
  }

  public abstract String getXMLTag();

  public void setColor(int colorCode) {
    insideColorCode = colorCode;
  }

  // public static class EllipseCalc {
  //
  // public static double[][] calcAll(double[] xy1, double[] xy2, double[] xy3, double[] xy4) {
  // double m1 = calcSlope(xy1[0], xy1[1], xy2[0], xy2[1]);
  // double m2 = calcSlope(xy3[0], xy3[1], xy4[0], xy4[1]);
  //
  // double[] newPt1, newPt2, newPt3, newPt4;
  //
  // newPt1 = calc(m1, m2, xy1[0], xy1[1], xy3[0], xy3[1]);
  // newPt2 = calc(m1, m2, xy1[0], xy1[1], xy4[0], xy4[1]);
  // newPt3 = calc(m1, m2, xy2[0], xy2[1], xy3[0], xy3[1]);
  // newPt4 = calc(m1, m2, xy2[0], xy2[1], xy4[0], xy4[1]);
  //
  // return new double[][]{newPt1, newPt2, newPt3, newPt4};
  // }
  //
  // private static double calcSlope(double x1, double y1, double x2, double y2) {
  // return (y2 - y1) / (x2 - x1);
  // }
  //
  // private static double[] calc(double m1, double m2, double x1, double y1, double x3, double y3)
  // {
  // double xNew = ((y3 - m1 * x3) - (y1 - m2 * x1)) / (m2 - m1);
  // double yNew = (m2 * (y3 - m1 * x3) - m1 * (y1 - m2 * x1)) / (m2 - m1);
  // return new double[]{xNew, yNew};
  // }
  //
  // }

  public void setLevel(int lvl) {
    displayLevel = lvl;
  }

  public void setName(String newName) {
    popName = newName;
  }

}
