package org.genvisis.one.ben.fcs.gating;

import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Random;

import org.genvisis.common.Array;
import org.genvisis.common.Numbers;
import org.genvisis.one.ben.fcs.AbstractPanel2;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;
//import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.genvisis.one.ben.fcs.sub.EMModel;

import edu.stanford.facs.logicle.Logicle;

public abstract class Gate {

  protected int insideColorCode;

  protected String xmlTag;

  protected String popName;
  protected String id;
  protected String parentID;
  protected Gate parentGate;
  protected ArrayList<Gate> children = new ArrayList<Gate>();
  protected GateDimension xDim, yDim;
//  protected ArrayList<GateDimension> dimensions = new ArrayList<GateDimension>();
  // protected HashMap<String, boolean[]> gatingCache = new HashMap<String, boolean[]>();
  protected int displayLevel = 0;
  
  static final Random rand = new Random();

  private static String generateID() {
    String id;
    id = "ID" + (rand.nextInt(35000000) + 35000000);
    return id;
  }

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

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((id == null) ? 0 : id.hashCode());
    result = prime * result + ((popName == null) ? 0 : popName.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    Gate other = (Gate) obj;
    if (id == null) {
      if (other.id != null)
        return false;
    } else if (!id.equals(other.id))
      return false;
    if (popName == null) {
      if (other.popName != null)
        return false;
    } else if (!popName.equals(other.popName))
      return false;
    return true;
  }

  public Gate getParentGate() {
    return parentGate;
  }

  protected boolean[] parentGating = null;

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

  // public void clearCache() {
  // gatingCache.clear();
  // }

  public String getID() {
    return id;
  }

  public void setLevel(int lvl) {
    displayLevel = lvl;
  }

  public int getLevel() {
    return displayLevel;
  }

  public int getColorCode() {
    return insideColorCode;
  }

  public void setColor(int colorCode) {
    insideColorCode = colorCode;
  }

  public String getFullNameAndGatingPath() {
    StringBuilder full = new StringBuilder();

    if (!"".equals(getName())) {
      full.append(getName());
    }

    full.append(" (");
    full.append(getXDimension().paramName);
    if (getYDimension() != null) {
      full.append(" v ").append(getYDimension().paramName);
    }
    full.append(")");

    if (parentGate != null) {
      String p = parentGate.getFullNameAndGatingPath();
      full.insert(0, " / ");
      full.insert(0, p);
    }

    return full.toString();
  }

  public String getName() {
    return popName;
  }

  public void setName(String newName) {
    popName = newName;
  }

  public ArrayList<Gate> getChildGates() {
    return children;
  }

  public int countAllChildren() {
    int cnt = getChildGates().size();
    for (Gate g2 : getChildGates()) {
      cnt += g2.countAllChildren();
    }
    return cnt;
  }

  public void addChildGate(Gate rg) {
    children.add(rg);
  }

  public GateDimension getXDimension() {
    return xDim;
  }

  public void setXDimension(GateDimension gd) {
    this.xDim = gd;
  }
  
  public GateDimension getYDimension() {
    return yDim;
  }
  
  public void setYDimension(GateDimension gd) {
    this.yDim = gd;
  }
  
//  public ArrayList<GateDimension> getDimensions() {
//    return dimensions;
//  }

//  public void addDimension(GateDimension gd) {
//    dimensions.add(gd);
//  }

  public abstract boolean[] gate(FCSDataLoader dataLoader);
  public abstract Gate copy(Gate parentGate);

  public abstract String getXMLTag();

  public static class RectangleGate extends Gate {

    public RectangleGate(Gate parentGate) {
      super(parentGate);
    }

    public RectangleGate(Gate parentGate, String popName, String id) {
      super(parentGate, popName, id);
    }

    public RectangleGate(Gate parentGate, String popName) {
      super(parentGate, popName);
    }
    
    @Override
    public Gate copy(Gate parentGate) {
      RectangleGate rg = new RectangleGate(parentGate, this.popName);
      rg.setXDimension(new RectangleGateDimension(rg, this.getXDimension().paramName, ((RectangleGateDimension) this.getXDimension()).getMin(), ((RectangleGateDimension) this.getXDimension()).getMax()));
      if (this.getYDimension() != null) {
        rg.setYDimension(new RectangleGateDimension(rg, this.getYDimension().paramName, ((RectangleGateDimension) this.getYDimension()).getMin(), ((RectangleGateDimension) this.getYDimension()).getMax()));
      }
      rg.setLevel(this.getLevel());
      for (Gate c : children) {
        rg.children.add(c.copy(rg));
      }
      return rg;
    }

    @Override
    public String getXMLTag() {
      return "gating:RectangleGate";
    }

    @Override
    public void setXDimension(GateDimension gd) {
      if (!(gd instanceof RectangleGateDimension)) {
        System.err.println("Error - only RectangleGateDimensions can be added to a RectangleGate");
        return;
      }
      super.setXDimension(gd);
      parentGating = null;
    }

    @Override
    public void setYDimension(GateDimension gd) {
      if (!(gd instanceof RectangleGateDimension)) {
        System.err.println("Error - only RectangleGateDimensions can be added to a RectangleGate");
        return;
      }
      super.setYDimension(gd);
      parentGating = null;
    }

    @Override
    public boolean[] gate(FCSDataLoader dataLoader) {
      // if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
      // return gatingCache.get(dataLoader.getLoadedFile());
      // }
      boolean[] includes = parentGate == null ? new boolean[dataLoader.getCount()] : parentGate.gate(dataLoader);
      if (includes == null) {
        return includes;
      }
      parentGating = parentGate == null ? Array.booleanArray(dataLoader.getCount(), true) : Arrays.copyOf(includes, includes.length);
      boolean[][] paramIncludes = new boolean[getYDimension() == null ? 1 : 2][dataLoader.getCount()];
      
      RectangleGateDimension rgd = (RectangleGateDimension) getXDimension();
      if (!dataLoader.containsParam(rgd.paramName)) {
        return null;
      }
      double[] paramData = dataLoader.getData(rgd.paramName, true);
      for (int i = 0; i < dataLoader.getCount(); i++) {
        // inclusive min, exclusive max - see gating-ml spec
        paramIncludes[0][i] = (!Numbers.isFinite(rgd.getMin()) || rgd.getMin() <= paramData[i]) && (!Numbers.isFinite(rgd.getMax()) || rgd.getMax() > paramData[i]);
      }
       
      if ((rgd = (RectangleGateDimension) getYDimension()) != null) {
        if (!dataLoader.containsParam(rgd.paramName)) {
          return null;
        }
        paramData = dataLoader.getData(rgd.paramName, true);
        for (int i = 0; i < dataLoader.getCount(); i++) {
          // inclusive min, exclusive max - see gating-ml spec
          paramIncludes[1][i] = (!Numbers.isFinite(rgd.getMin()) || rgd.getMin() <= paramData[i]) && (!Numbers.isFinite(rgd.getMax()) || rgd.getMax() > paramData[i]);
        }
        
      }
      
      for (int i = 0; i < dataLoader.getCount(); i++) {
        boolean include = true;
        if (parentGate != null && !includes[i]) {
          continue;
        }
        for (int p = 0, pCount = paramIncludes.length; p < pCount; p++) {
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

  }

  public static class PolygonGate extends Gate {
    private static int DEFAULT_GATE_RESOLUTION = 256; // default
    private static int gateResolution = DEFAULT_GATE_RESOLUTION; // default
    private static int range = 262144; // TODO make this dynamic based on data loader param range
                                       // (set when dims are set) - may have different ranges for
                                       // each dim
//    private static ArrayList<Rectangle> DEFAULT_RECTS = new ArrayList<Rectangle>();
//    private static ArrayList<Rectangle> rects = new ArrayList<Rectangle>();
    private static int binStep = range / gateResolution;
    private static Rectangle[][] rectArray;
    static {
//      rects = prepRects();
//      DEFAULT_RECTS = rects;
      setupRects();
    }
    static void setupRects() {
      rectArray = new Rectangle[gateResolution * 2][gateResolution * 2];
      for (int i = -gateResolution; i < gateResolution; i++) {
        for (int j = -gateResolution; j < gateResolution; j++) {
          rectArray[i + gateResolution][j + gateResolution] = new Rectangle(i * binStep + binStep / 2, j * binStep + binStep / 2, binStep, binStep);
        }
      }
    }
//    static ArrayList<Rectangle> prepRects() {
//      ArrayList<Rectangle> rects = new ArrayList<Rectangle>();
//      for (int i = -gateResolution; i < gateResolution; i++) {
//        for (int j = -gateResolution; j < gateResolution; j++) {
//          rects.add(new Rectangle(i * binStep + binStep / 2, j * binStep + binStep / 2, binStep, binStep));
//        }
//      }
//      return rects;
//    }

    private final ArrayList<Double> verticesX = new ArrayList<Double>();
    private final ArrayList<Double> verticesY = new ArrayList<Double>();
    private Path2D myPath;
    private Path2D transformedPath;
    ArrayList<Rectangle> myRects = new ArrayList<Rectangle>();
    private boolean mimicFlowJo = false;

    public PolygonGate(Gate parentGate, String name, String id, boolean flowjo) {
      super(parentGate, name, id);
      setShouldMimicFlowJoGating(flowjo);
    }

    public PolygonGate(Gate parentGate, String name, String id) {
      this(parentGate, name, id, false);
    }

    public PolygonGate(Gate parentGate, String name, boolean flowjo) {
      this(parentGate, name, generateID(), flowjo);
    }

    public PolygonGate(Gate parentGate, String name) {
      this(parentGate, name, generateID(), false);
    }

    public PolygonGate(Gate parentGate) {
      this(parentGate, null, generateID(), false);
    }

    @Override
    public String getXMLTag() {
      return "gating:PolygonGate";
    }

    public boolean getMimicsFlowJo() {
      return mimicFlowJo;
    }

    public void setShouldMimicFlowJoGating(boolean mimic) {
      mimicFlowJo = mimic;
//      if (mimicFlowJo) {
//        prepGating();
//      }
      parentGating = null;
      transformedPath = null;
    }

    @Override
    public void setXDimension(GateDimension gd) {
      super.setXDimension(gd);
      parentGating = null;
      transformedPath = null;
    }
    
    @Override
    public void setYDimension(GateDimension gd) {
      super.setYDimension(gd);
      parentGating = null;
      transformedPath = null;
    }

    public void setPath(Path2D pth) {
      myPath = pth;
      parentGating = null;
      resetVertices();
//      if (mimicFlowJo) {
//        prepGating();
//      }
      transformedPath = null;
      // clearCache();
    }

    private void resetVertices() {
      verticesX.clear();
      verticesY.clear();
      double[] coords = new double[6];
      PathIterator pi = myPath.getPathIterator(null);
      while (!pi.isDone()) {
        int type = pi.currentSegment(coords);
        if (type != PathIterator.SEG_CLOSE) {
          addVertex(coords[0], coords[1]);
        }
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
      for (int i = 1; i < verticesX.size(); i++) {
        path.lineTo(verticesX.get(i), verticesY.get(i));
      }
      path.closePath();
      return path;
    }
    
    private Path2D transformPath(FCSDataLoader dataLoader) {
      Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
      AxisTransform xTr, yTr;
      xTr = dataLoader.getParamTransform(getXDimension().paramName);
      yTr = dataLoader.getParamTransform(getYDimension().paramName);
      if (verticesX.isEmpty() || verticesY.isEmpty()) {
        resetVertices();
      }
//      path.moveTo(x ? verticesX.get(0) : (xTr.scaleX(verticesX.get(0)) * gateResolution), y ? verticesY.get(0) : (yTr.scaleY(verticesY.get(0)) * gateResolution));
      path.moveTo(xTr.scaleX(verticesX.get(0)),yTr.scaleY(verticesY.get(0)));
      for (int i = 1; i < verticesX.size(); i++) {
        path.lineTo(xTr.scaleX(verticesX.get(i)),  yTr.scaleY(verticesY.get(i)));
      }
      path.closePath();
      return path;
    }
    

    @Override
    public Gate copy(Gate parentGate) {
      PolygonGate pg = new PolygonGate(parentGate, this.popName, this.mimicFlowJo);
      pg.setXDimension(new GateDimension(pg, getXDimension().paramName));
      pg.setYDimension(new GateDimension(pg, getYDimension().paramName));
      pg.setLevel(this.getLevel());
      pg.setPath((Path2D) this.myPath.clone());
      for (Gate c : children) {
        pg.children.add(c.copy(pg));
      }
      return pg;
    }
    
    @Override
    public boolean[] gate(FCSDataLoader dataLoader) {
      // if (gatingCache.containsKey(dataLoader.getLoadedFile())) {
      // return gatingCache.get(dataLoader.getLoadedFile());
      // }
      boolean[] includes = parentGate == null ? new boolean[dataLoader.getCount()] : parentGate.gate(dataLoader);
      if (includes == null) {
        return includes;
      }
      parentGating = parentGate == null ? Array.booleanArray(dataLoader.getCount(), true) : Arrays.copyOf(includes, includes.length);
      if (myPath == null) {
        myPath = constructPath();
      }
      if (transformedPath == null) {
        transformedPath = transformPath(dataLoader);
        prepGating(dataLoader);
      }
      if (!dataLoader.containsParam(getXDimension().paramName) || !dataLoader.containsParam(getYDimension().paramName)) {
        return null;
      }
      double[][] paramData = {
          dataLoader.getData(getXDimension().paramName, true),
          dataLoader.getData(getYDimension().paramName, true)
      };

      boolean xT, yT;
      AxisTransform xTr, yTr;
      xTr = dataLoader.getParamTransform(getXDimension().paramName);
      yTr = dataLoader.getParamTransform(getYDimension().paramName);
      for (int i = 0; i < dataLoader.getCount(); i++) {
        if (parentGate != null && !includes[i]) {
          continue;
        }
        boolean include = false;
        if (mimicFlowJo) {
//          int xInd, yInd;
//          xInd = (int) (dataLoader.getScaleForParam(getDimensions().get(0).paramName) == AXIS_SCALE.LIN ? ((paramData[0][i] / binStep) + gateResolution) : xTr.scaleX(paramData[0][i]));
//          yInd = (int) (dataLoader.getScaleForParam(getDimensions().get(1).paramName) == AXIS_SCALE.LIN ? ((paramData[1][i] / binStep) + gateResolution) : yTr.scaleY(paramData[1][i]));
//          if (myRects.contains(rectArray[xInd][yInd])) {
//            include = true;
//          }
          
          double x, y;
          x = xTr.scaleX(paramData[0][i]);
          y = yTr.scaleY(paramData[1][i]);
          for (Rectangle rect : myRects) {
//            if (rect.contains(xT ? paramData[0][i] : xTr.scaleX(paramData[0][i]) * gateResolution * binStep, yT ? paramData[1][i] : yTr.scaleY(paramData[1][i]) * gateResolution * binStep)) {
            if (rect.contains(paramData[0][i], paramData[1][i])) {
//            if (rect.contains(x, y)) {
              include = true;
              break;
            }
          }
        } else {
//          double x = paramData[0][i];
//          double y = paramData[1][i];
          double x = xTr.scaleX(paramData[0][i]);
          double y = yTr.scaleY(paramData[1][i]);
//          if (myPath.contains(x, y)) {
//            include = true;
//          }
          if (transformedPath.contains(x, y)) {
            include = true;
          }
        }
        includes[i] = include;
      }

      // gatingCache.put(dataLoader.getLoadedFile(), includes);
      return includes;
    }

    void prepGating(FCSDataLoader dataLoader) {

      AxisTransform xTr, yTr;
//      boolean xT, yT;
      xTr = dataLoader.getParamTransform(getXDimension().paramName);
      yTr = dataLoader.getParamTransform(getYDimension().paramName);
      
      ArrayList<Rectangle> vertexRects = new ArrayList<Rectangle>();
      for (int i = 0; i < verticesX.size(); i++) {
//        double x, y;
        int xInd, yInd;
//        xInd = (int) (dataLoader.getScaleForParam(getDimensions().get(0).paramName) == AXIS_SCALE.LIN ? ((verticesX.get(i) / binStep) + gateResolution) : xTr.scaleX(verticesX.get(i)));
//        yInd = (int) (dataLoader.getScaleForParam(getDimensions().get(1).paramName) == AXIS_SCALE.LIN ? ((verticesY.get(i) / binStep) + gateResolution) : yTr.scaleY(verticesY.get(i)));
//        xInd = (int) ((xTr.scaleX(verticesX.get(i)) / binStep) + gateResolution);
//        yInd = (int) ((yTr.scaleY(verticesY.get(i)) / binStep) + gateResolution);
        xInd = (int) ((verticesX.get(i) / binStep) + gateResolution);
        yInd = (int) ((verticesY.get(i) / binStep) + gateResolution);
//        x = verticesX.get(i);
//        y = verticesY.get(i);
//        for (int j = 0; j < rectArray.length; j++) {
//          for (int k = 0; k < rectArray[j].length; k++) {
//            if (rectArray[j][k].contains(x, y)) {
//              vertexRects.add(rectArray[j][k]);
//            }
//          }
//        }
        vertexRects.add(rectArray[xInd][yInd]);
      }
//      for (Rectangle r : rects) {
//        for (int i = 0; i < verticesX.size(); i++) {
//          double x, y;
//          x = xT ? verticesX.get(i) : xTr.scaleX(verticesX.get(i)) * gateResolution * binStep;
//          y = yT ? verticesY.get(i) : yTr.scaleY(verticesY.get(i)) * gateResolution * binStep;
////          x = verticesX.get(i);
////          y = verticesY.get(i);
//          if (r.contains(x, y)) {
//            vertexRects.add(r);
//          }
//        }
//      }

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
          double det =
              (o1.getCenterX() - x) * (o2.getCenterY() - y) - (o2.getCenterX() - x)
                  * (o1.getCenterY() - y);
          if (det < 0) {
            return -1;
          }
          if (det > 0) {
            return 1;
          }

          // points a and b are on the same line from the center
          // check which point is closer to the center
          double d1 =
              (o1.getCenterX() - x) * (o1.getCenterX() - x) + (o1.getCenterY() - y)
                  * (o1.getCenterY() - y);
          double d2 =
              (o2.getCenterX() - x) * (o2.getCenterX() - x) + (o2.getCenterY() - y)
                  * (o2.getCenterY() - y);
          return d1 > d2 ? -1 : d1 < d2 ? 1 : 0;
        }
      });

      double x1 = vertexRects.get(0).getCenterX();
      double y1 = vertexRects.get(0).getCenterY();
      Path2D path = new Path2D.Double();
//      path.moveTo(xT ? x1 : xTr.scaleX(x1) * gateResolution * binStep, yT ? y1 : yTr.scaleY(y1) * gateResolution * binStep);
      path.moveTo(x1, y1);
      for (int i = 1; i < vertexRects.size(); i++) {
        x1 = vertexRects.get(i).getCenterX();
        y1 = vertexRects.get(i).getCenterY();
//        path.lineTo(xT ? x1 : xTr.scaleX(x1) * gateResolution * binStep, yT ? y1 : yTr.scaleY(y1) * gateResolution * binStep);
        path.lineTo(x1, y1);
      }
      path.closePath();

      myRects.clear();
      for (int i = 0; i < rectArray.length; i++) {
        for (int j = 0; j < rectArray[i].length; j++) {
          Rectangle rect = rectArray[i][j];
          if (vertexRects.contains(rect) || path.contains(rect)
              || (path.intersects(rect) && path.contains(rect.getCenterX(), rect.getCenterY()))) {
            myRects.add(rect);
          }
        }
      }
//      for (Rectangle rect : rects) {
//        if (vertexRects.contains(rect) || path.contains(rect)
//            || (path.intersects(rect) && path.contains(rect.getCenterX(), rect.getCenterY()))) {
//          myRects.add(rect);
//        }
//      }
    }

    public void setGateResolution(int res) {
      if (res == gateResolution) {
        return;
      }
      gateResolution = res;
      binStep = range / gateResolution;
      if (res == DEFAULT_GATE_RESOLUTION) {
//        rects = DEFAULT_RECTS;
//      } else {
//        rects = prepRects();
        setupRects();
      }
    }

    public void addVertex(Double fX, Double fY) {
      verticesX.add(fX);
      verticesY.add(fY);
      parentGating = null;
      transformedPath = null;
    }

    public void transform(AffineTransform at) {
      myPath.transform(at);
      resetVertices();
//      if (mimicFlowJo) {
//        prepGating();
//      }
      parentGating = null;
      transformedPath = null;
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
    public String getXMLTag() {
      return "gating:EllipsoidGate";
    }

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
    public Gate copy(Gate parentGate) {
      throw new UnsupportedOperationException();
    }
    
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

  public static class QuadrantGate extends Gate {
    public QuadrantGate() {
      super(null);
      throw new UnsupportedOperationException();
    }

    @Override
    public String getXMLTag() {
      return "gating:RectangleGate"; // TODO stored as four RectangleGate tags
    }

    // UNSUPPORTED
    @Override
    public boolean[] gate(FCSDataLoader dataLoader) {
      throw new UnsupportedOperationException();
    }

    @Override
    public Gate copy(Gate parentGate) {
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

    @Override
    public String getXMLTag() {
      return "gating:BooleanGate";
    }

    @Override
    public Gate copy(Gate parentGate) {
      throw new UnsupportedOperationException();
    }
  }

}
