package org.genvisis.one.ben.fcs;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.flowcyt.cfcs.CFCSAbstractData;
import org.flowcyt.cfcs.CFCSData;
import org.flowcyt.cfcs.CFCSDataSet;
import org.flowcyt.cfcs.CFCSError;
import org.flowcyt.cfcs.CFCSKeywords;
import org.flowcyt.cfcs.CFCSListModeData;
import org.flowcyt.cfcs.CFCSParameter;
import org.flowcyt.cfcs.CFCSParameters;
import org.flowcyt.cfcs.CFCSSpillover;
import org.flowcyt.cfcs.CFCSSystem;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;

import scala.annotation.target.getter;
import edu.stanford.facs.logicle.Logicle;

public class FCSDataLoader {

  private static final String COMPENSATED_PREPEND = "Comp-";
  private static final int COMP_LEN = COMPENSATED_PREPEND.length();

  //
  // public static void main(String[] args) {
  // // String fcsFilename = "F:\\Flow\\P1-B&C-CD3-APC-Cy7 or CD4-APC-Cy7_ULTRA BRIGHT RAINBOW
  // BEADS_URB_001.fcs";
  // // String fcsFilename = "F:\\Flow\\P1- PBMC-A&C rest_panel one_PBMC-C P1 1HR rest_003.fcs";
  // String fcsFilename = "F:\\Flow\\P1- PBMC-A&C rest_panel one_PBMC-A P1 1HR rest_002.fcs";
  // try {
  // (new FCSDataLoader()).loadData(fcsFilename);
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  // }

  public static enum LOAD_STATE {
    LOADED, LOADING, PARTIALLY_LOADED, LOADING_REMAINDER, UNLOADED;
  }

  private volatile LOAD_STATE state = LOAD_STATE.UNLOADED;

  ArrayList<String> paramNamesInOrder;
  ArrayList<String> paramShortNamesInOrder;
  ArrayList<String> paramLongNamesInOrder;
  private HashMap<String, AXIS_SCALE> paramScales;
  public HashMap<String, AxisTransform> paramTransforms;
  private final ArrayList<Integer> ranges;
  LinkedHashSet<String> compensatedNames;
  HashMap<String, Integer> compensatedIndices;
  int eventCount = -1;
  int loadedCount = 0;
  double[][] allData;
  double[][] compensatedData;
  private String loadedFile = null;
  CFCSData dataObj = null;
  CFCSSpillover spillObj = null;
  int paramsCount = -1;
  boolean isTransposed = false;
  ArrayList<Integer> indicesToLoad = new ArrayList<Integer>();
  Date lastModified;
  Date runDate;
  Date beginTime;
  Date endTime;
  double[] paramScaling;

  public enum DATA_SET {
    ALL, COMPENSATED, UNCOMPENSATED;
  }

  public FCSDataLoader() {
    paramNamesInOrder = new ArrayList<String>();
    paramShortNamesInOrder = new ArrayList<String>();
    paramLongNamesInOrder = new ArrayList<String>();
    paramScales = new HashMap<>();
    paramTransforms = new HashMap<String, AxisTransform>();
    ranges = new ArrayList<Integer>();
    compensatedNames = new LinkedHashSet<String>();
    compensatedIndices = new HashMap<String, Integer>();
  }

  public void emptyAndReset() {
    loadInBGThread.interrupt();
    try {
      loadInBGThread.join();
    } catch (InterruptedException e) {
      /* wouldn't this be what we want? */}
    setState(LOAD_STATE.UNLOADED);
    paramNamesInOrder = new ArrayList<String>();
    paramShortNamesInOrder = new ArrayList<String>();
    paramLongNamesInOrder = new ArrayList<String>();
    paramScales = new HashMap<>();
    paramTransforms = new HashMap<String, AxisTransform>();
    compensatedNames = new LinkedHashSet<String>();
    compensatedIndices = new HashMap<String, Integer>();
    eventCount = -1;
    loadedCount = 0;
    allData = null;
    compensatedData = null;
    loadedFile = null;
    dataObj = null;
    spillObj = null;
    paramsCount = -1;
    isTransposed = false;
    indicesToLoad = new ArrayList<Integer>();
    paramScaling = null;
    runDate = null;
    beginTime = null;
    endTime = null;
    System.gc();
  }

  public String getLoadedFile() {
    return loadedFile;
  }

  private/* synchronized */void setState(LOAD_STATE state) {
    this.state = state;
  }

  public/* synchronized */LOAD_STATE getLoadState() {
    return state;
  }

  public int[] getLoadedStatus() {
    if (getLoadState() == LOAD_STATE.LOADED) {
      return new int[] {1, 1}; // complete
    } else if (eventCount >= 0) {
      return new int[] {loadedCount, eventCount}; // in progress
    } else {
      return new int[] {-1, -1}; // indeterminate
    }
  }

  public int getCount() {
    LOAD_STATE currState = getLoadState();
    if (currState == LOAD_STATE.LOADED) {
      return eventCount;
    } else if (currState == LOAD_STATE.UNLOADED || currState == LOAD_STATE.LOADING) {
      return 0;
    } else {
      return loadedCount;
    }
  }

  public Date getLastModified() {
    return lastModified;
  }

  public Date getRunDate() {
    return runDate;
  }

  public Date getBeginTime() {
    return beginTime;
  }

  public Date getEndTime() {
    return endTime;
  }

  private String[] loadTimeData(CFCSKeywords keys) {
    String dKey = "$DATE";
    String bKey = "$BTIM";
    String eKey = "$ETIM";

    String dVal = null, eVal = null, bVal = null;
    try {
      dVal = keys.getKeyword(dKey).getKeywordValue();
    } catch (Exception e) {
    }
    try {
      bVal = keys.getKeyword(bKey).getKeywordValue();
    } catch (Exception e) {
    }
    try {
      eVal = keys.getKeyword(eKey).getKeywordValue();
    } catch (Exception e) {
    }

    return new String[] {dVal, bVal, eVal};
  }

  private double[] scaling(CFCSKeywords keys, int cnt) {
    double[] scal = new double[cnt];
    for (int i = 0; i < cnt; i++) {
      String k = "$P" + (i + 1) + "G";
      scal[i] = 1;
      try {
        scal[i] = keys.getKeyword(k).getKeywordDoubleValue();
      } catch (Exception e) {
      }
    }
    return scal;
  }

  public void loadData(String fcsFilename) throws IOException {
    // synchronized(this) {
    if (getLoadState() != LOAD_STATE.UNLOADED) {
      return;
    }
    setState(LOAD_STATE.LOADING);
    // }
    loadedFile = fcsFilename;

    CFCSSystem syst = new CFCSSystem();
    File sysFile = new File(fcsFilename);
    URL fileURL = (sysFile).toURI().toURL();
    syst.open(fileURL);

    CFCSDataSet dset = syst.getDataSet(0);
    dataObj = dset.getData();
    eventCount = ((CFCSAbstractData) dataObj).getCount();
    CFCSParameters params = dset.getParameters();
    paramsCount = params.getCount();
    CFCSKeywords keys = dset.getKeywords();

    String[] dateTimes = loadTimeData(keys);
    SimpleDateFormat sdfDate = new SimpleDateFormat("dd-MMM-yyyy");
    SimpleDateFormat sdfTime = new SimpleDateFormat("HH:mm:ss");
    if (dateTimes[0] != null) {
      try {
        runDate = sdfDate.parse(dateTimes[0]);
      } catch (ParseException e1) {
        System.out.println("Unable to parse date: " + dateTimes[0]);
      }
    }
    if (dateTimes[1] != null) {
      try {
        beginTime = sdfTime.parse(dateTimes[1]);
      } catch (ParseException e1) {
        System.out.println("Unable to parse time: " + dateTimes[1]);
      }
    }
    if (dateTimes[2] != null) {
      try {
        endTime = sdfTime.parse(dateTimes[2]);
      } catch (ParseException e1) {
        System.out.println("Unable to parse time: " + dateTimes[2]);
      }
    }
    paramScaling = scaling(keys, paramsCount);

    spillObj = keys.getSpillover();
    lastModified = keys.getLastModified();
    if (lastModified == null) {
      System.err.println("Warning - FCS file " + fcsFilename
          + " does NOT contain a last modified date - using the last-modified system file date.");
      lastModified = new Date(sysFile.lastModified());
    }

    HashMap<String, String> names = new HashMap<String, String>();
    for (int i = 0; i < paramsCount; i++) {
      CFCSParameter param = params.getParameter(i);
      String name = null;
      String shortName = null;
      try {
        name = param.getFullName();
      } catch (Exception e) {
      }
      try {
        shortName = param.getShortName();
      } catch (Exception e) {
      }
      // if (name == null) {
      // // TODO Error, no name set (included in spillover matrix?) -- will this scenario ever
      // happen?
      // name = "P" + i;
      // }
      String actName = shortName;
      if (name != null) {
        actName = shortName + " (" + name + ")";
      }
      paramNamesInOrder.add(actName);
      paramShortNamesInOrder.add(shortName);
      if (name != null) {
        paramLongNamesInOrder.add(name);
      }
      names.put(shortName, actName);

      String axisKeywork = "P" + (i + 1) + "DISPLAY";
      AXIS_SCALE scale = AXIS_SCALE.LIN;
      try {
        scale = AXIS_SCALE.valueOf(keys.getKeyword(axisKeywork).getKeywordValue());
      } catch (Exception e) {
        System.err.println("Warning - no axis scale set for parameter " + paramNamesInOrder.get(i)
            + "; assuming a linear scale.");
      };
      paramScales.put(actName, scale);
      paramTransforms.put(actName, getDefaultTransform(scale));

      String rangeKeyword = "$P" + (i + 1) + "R";
      int rng = -1;
      try {
        rng = keys.getKeyword(rangeKeyword).getKeywordIntegerValue();
      } catch (Exception e) {
        System.err.println("Warning - no parameter range value for parameter "
            + paramNamesInOrder.get(i) + "; assuming standard of 262144");
        rng = 262144;
      }
      ranges.add(rng);

    }

    String[] arr = spillObj.getParameterNames();
    for (int i = 0, count = arr.length; i < count; i++) {
      compensatedNames.add(names.get(arr[i]));
      compensatedIndices.put(names.get(arr[i]), i);
      paramScales.put(COMPENSATED_PREPEND + names.get(arr[i]), getScaleForParam(names.get(arr[i])));
      paramTransforms.put(COMPENSATED_PREPEND + names.get(arr[i]), getDefaultTransform(getScaleForParam(names.get(arr[i]))));
    }

    if (dataObj.getType() == CFCSData.LISTMODE) {
      CFCSListModeData listData = ((CFCSListModeData) dataObj);
      allData = new double[listData.getCount()][];
      setState(LOAD_STATE.PARTIALLY_LOADED);
      for (int i = 0; i < listData.getCount(); i++) {
        double[] newData = Array.doubleArray(paramsCount, Double.NaN);
        try {
          listData.getEvent/* AsInTheFile */(i, newData); // should be getEventAsInTheFile???
          if (Double.isNaN(newData[0])) {
            indicesToLoad.add(i);
            loadedCount--;
          }
        } catch (CFCSError e) {
          indicesToLoad.add(i);
          loadedCount--;
        }
        for (int p = 0; p < paramsCount; p++) {
          newData[p] *= paramScaling[p];
        }
        allData[i] = newData;
        loadedCount++;
      }
      setState(LOAD_STATE.LOADING_REMAINDER);
      loadInBGThread.start();
    } else {
      System.err.println("Error - UNSUPPORTED DATA TYPE.");
    }
  }

  private AxisTransform getDefaultTransform(AXIS_SCALE scale) {
    switch (scale) {
      case LOG:
        System.err.println("Error - NO DEFINED AXIS TRANSFORM FOR LOG AXES!");
      case BIEX:
        return createBiexAxisTransform();
      case LIN:
      default:
        return createPassthroughTransform();
    }
  }

  private static AxisTransform createPassthroughTransform() {
    return new AxisTransform(null) {
      @Override
      public double scaleY(double val) { return val; }
      @Override
      public double scaleX(double val) { return val; }
      @Override
      public double inverseY(double val) { return val; }
      @Override
      public double inverseX(double val) { return val; }
    };
  }
  
  static int biexResMinX = 0;
  static int biexResMaxX = 256;
  static int biexResMinY = 0;
  static int biexResMaxY = 256;
  public void setBiexRangeX(int min, int max) {
    biexResMinX = min;
    biexResMaxX = max;
  }
  public void setBiexRangeY(int min, int max) {
    biexResMinY = min;
    biexResMaxY = max;
  }
  
  private static AxisTransform createBiexAxisTransform() {
    double DEFAULT_T = 262144;
    double DEFAULT_W = Math.log10(Math.abs(-100));
    double DEFAULT_M = Math.log10(DEFAULT_T);
    double DEFAULT_A = Math.min(0, 1);
    Logicle lgl = new Logicle(DEFAULT_T, DEFAULT_W, DEFAULT_M, DEFAULT_A);
    
    int DEFAULT_RESOLUTION = 256;
    int resolution = DEFAULT_RESOLUTION;
    
    return new AxisTransform(null) {
      
      private double scaleLin(double val, int biexMin, int biexMax) {
        int[] screenMinMax = {biexMin, biexMax};
        double[] plotMinMax = {lgl.inverse(0), lgl.inverse(1)};
        return (int) ((val - plotMinMax[0]) / (plotMinMax[1] - plotMinMax[0]) * (screenMinMax[1] - screenMinMax[0])) + screenMinMax[0];
      }
      
      @Override
      public double scaleY(double val) {
        return scaleLin(lgl.scale(val) * lgl.inverse(1), biexResMinY, biexResMaxY);
//        return lgl.scale(val);
      }
      
      @Override
      public double scaleX(double val) {
        return scaleLin(lgl.scale(val) * lgl.inverse(1), biexResMinX, biexResMaxX);
//        return lgl.scale(val);
      }
      
      @Override
      public double inverseY(double val) {
        return lgl.inverse(val / lgl.inverse(1));
//        return lgl.inverse(val);
      }
      
      @Override
      public double inverseX(double val) {
        return lgl.inverse(val / lgl.inverse(1));
//        return lgl.inverse(val);
      }
    };
  }
  
  public AxisTransform getParamTransform(String param) {
    return paramTransforms.get(getInternalParamName(param));
  }
  
  Thread loadInBGThread = new Thread(new Runnable() {
    @Override
    public void run() {
      if (indicesToLoad.size() > 0) {
        CFCSListModeData listData = ((CFCSListModeData) dataObj);
        while (!listData.isLoaded() && !Thread.currentThread().isInterrupted()) {
          Thread.yield();
        }
        while (indicesToLoad.size() > 0 && !Thread.currentThread().isInterrupted()) {
          int indToLoad = indicesToLoad.get(indicesToLoad.size() - 1);
          double[] newData = Array.doubleArray(paramsCount, Double.NaN);
          try {
            listData.getEvent/* AsInTheFile */(indToLoad, newData); // should be
                                                                    // getEventAsInTheFile???
            indicesToLoad.remove(indicesToLoad.size() - 1);
          } catch (CFCSError e) {
            continue;
          }
          for (int p = 0; p < paramsCount; p++) {
            newData[p] *= paramScaling[p];
          }
          allData[indToLoad] = newData;
          loadedCount++;
        }
      }
      if (Thread.currentThread().isInterrupted()) {
        return;
      }
      compensatedData =
          compensateSmall(paramNamesInOrder, allData,
              compensatedNames.toArray(new String[compensatedNames.size()]),
              getInvertedSpilloverMatrix(spillObj.getSpilloverCoefficients()));
      if (Thread.currentThread().isInterrupted()) {
        return;
      }
      allData = Matrix.transpose(allData);
      if (Thread.currentThread().isInterrupted()) {
        return;
      }
      compensatedData = Matrix.transpose(compensatedData);
      if (Thread.currentThread().isInterrupted()) {
        return;
      }
      isTransposed = true;
      dataObj = null;
      setState(LOAD_STATE.LOADED);
      System.gc();
    }
  });

  public void setSpilloverMatrix(double[][] newCoeffs) {
    String[] names = compensatedNames.toArray(new String[compensatedNames.size()]);
    compensatedData =
        compensateSmall(paramNamesInOrder, allData, names, getInvertedSpilloverMatrix(newCoeffs));
    spillObj = CFCSSpillover.createSpilloverMatrix(names, newCoeffs);
    // TODO set as internal obj to CFCS objects???

  }
  
  public double[] getData(String colName, boolean waitIfNecessary) {
    String columnName = colName.startsWith(COMPENSATED_PREPEND) ? COMPENSATED_PREPEND + getInternalParamName(colName.substring(COMP_LEN)) : getInternalParamName(colName);
    LOAD_STATE currState = getLoadState();
    if (columnName.startsWith(COMPENSATED_PREPEND)) {
      if (currState == LOAD_STATE.LOADED) {
        if (isTransposed) {
          return compensatedData[compensatedIndices.get(columnName.substring(COMP_LEN))];
        } else {
          return Matrix.extractColumn(compensatedData, compensatedIndices.get(columnName.substring(COMP_LEN)));
        }
      } else {
        if (currState != LOAD_STATE.UNLOADED && waitIfNecessary) {
          while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
            Thread.yield();
          }
          if (isTransposed) {
            return compensatedData[compensatedIndices.get(columnName.substring(COMP_LEN))];
          } else {
            return Matrix.extractColumn(compensatedData,
                compensatedIndices.get(columnName.substring(COMP_LEN)));
          }
        } else {
          int len = eventCount == -1 ? 0 : eventCount;
          return Array.doubleArray(len, Double.NaN);
        }
      }
    } else {
      if (currState == LOAD_STATE.UNLOADED) {
        return new double[0];
      } else {
        if (currState == LOAD_STATE.LOADING && !waitIfNecessary) {
          int len = eventCount == -1 ? 0 : eventCount;
          return Array.doubleArray(len, Double.NaN);
        } else {
          if (currState != LOAD_STATE.LOADED && waitIfNecessary) {
            while ((currState = getLoadState()) != LOAD_STATE.LOADED) { // TODO wait for complete
                                                                        // data, or at least some?
              Thread.yield();
            }
          }
          // TODO WARNING, RETURNED DATA /MAY/ BE INCOMPLETE
          if (isTransposed) {
            return allData[paramNamesInOrder.indexOf(columnName)];
          } else {
            return Matrix.extractColumn(allData, paramNamesInOrder.indexOf(columnName));
          }
        }
      }
    }
  }

  private String getInternalParamName(String name) {
    String nm = name;
    boolean prepend = nm.startsWith(COMPENSATED_PREPEND); 
    if (prepend) {
      nm = nm.substring(COMP_LEN);
    }
    int ind = paramNamesInOrder.indexOf(nm);
    if (ind != -1) {
      return prepend ? COMPENSATED_PREPEND + nm : nm;
    } else if ((ind = paramShortNamesInOrder.indexOf(nm)) != -1) {
      return prepend ? COMPENSATED_PREPEND + paramNamesInOrder.get(ind) : paramNamesInOrder.get(ind);
    } else if ((ind = paramLongNamesInOrder.indexOf(nm)) != -1) {
      return prepend ? COMPENSATED_PREPEND + paramNamesInOrder.get(ind) : paramNamesInOrder.get(ind);
    }
    return prepend ? COMPENSATED_PREPEND + nm : nm;
  }
  
  public double[] getDataLine(ArrayList<String> params, int ind) {
    double[] line = new double[params.size()];
    
    LOAD_STATE currState = getLoadState();
    for (int i = 0; i < params.size(); i++) {
      String columnName = getInternalParamName(params.get(i));
      if (columnName.startsWith(COMPENSATED_PREPEND)) {
        if (currState == LOAD_STATE.LOADED) {
          if (isTransposed) {
            line[i] = compensatedData[compensatedIndices.get(columnName.substring(COMP_LEN))][ind];
          } else {
            line[i] = compensatedData[ind][compensatedIndices.get(columnName.substring(COMP_LEN))];
          }
        } else {
          if (currState != LOAD_STATE.UNLOADED) {
            while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
              Thread.yield();
            }
            if (isTransposed) {
              line[i] = compensatedData[compensatedIndices.get(columnName.substring(COMP_LEN))][ind];
            } else {
              line[i] = compensatedData[ind][compensatedIndices.get(columnName.substring(COMP_LEN))];
            }
          } else {
            line[i] = Double.NaN;
          }
        }
      } else {
        if (currState == LOAD_STATE.UNLOADED) {
          line[i] = Double.NaN;
        } else {
          if (currState == LOAD_STATE.LOADING) {
            line[i] = Double.NaN;
          } else {
            if (currState != LOAD_STATE.LOADED) {
              while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
                Thread.yield();
              }
            }
            // TODO WARNING, RETURNED DATA /MAY/ BE INCOMPLETE
            if (isTransposed) {
              line[i] = allData[paramNamesInOrder.indexOf(columnName)][ind];
            } else {
              line[i] = allData[ind][paramNamesInOrder.indexOf(columnName)];
            }
          }
        }
      }
    }
    return line;
  }

  public boolean containsParam(String paramName) {
    String nm = paramName;
    if (paramName.startsWith(COMPENSATED_PREPEND)) {
      nm = paramName.substring(COMP_LEN);
    }
    return (compensatedNames.contains(nm) || paramNamesInOrder.contains(nm)) || paramShortNamesInOrder.contains(nm) || paramLongNamesInOrder.contains(nm);
  }
  
  public ArrayList<String> getAllDisplayableNames(DATA_SET set) {
    ArrayList<String> names = new ArrayList<String>();
    switch (set) {
      case ALL:
        names = new ArrayList<String>(paramNamesInOrder);
        for (int i = 0; i < names.size(); i++) {
          if (compensatedNames.contains(names.get(i))) {
            names.add(COMPENSATED_PREPEND + names.get(i));
          }
        }
        break;
      case COMPENSATED:
        names = new ArrayList<String>(paramNamesInOrder);
        for (int i = 0; i < names.size(); i++) {
          if (compensatedNames.contains(names.get(i))) {
            names.set(i, COMPENSATED_PREPEND + names.get(i));
          }
        }
        break;
      case UNCOMPENSATED:
        names = new ArrayList<String>(paramNamesInOrder);
        break;
      default:
        break;
    }
    return names;
  }

  public void exportData(String outputFilename, DATA_SET set) {
    PrintWriter writer = Files.getAppropriateWriter(outputFilename);
    ArrayList<String> colNames = getAllDisplayableNames(set);
    for (int i = 0, count = colNames.size(); i < count; i++) {
      writer.print(colNames.get(i));
      if (i < count - 1) {
        writer.print(",");
      }
    }
    writer.println();

    switch (set) {
      case ALL:
        for (int i = 0; i < eventCount; i++) {
          for (double[] element : allData) {
            writer.print(ext.formDeci(element[i], 10));
            writer.print(",");
          }
          for (int p = allData.length, count = colNames.size(); p < count; p++) {
            writer.print(ext.formDeci(compensatedData[p - allData.length][i], 10));
            if (p < count - 1) {
              writer.print(",");
            }
          }
          writer.println();
        }

        break;
      case COMPENSATED:
        for (int i = 0; i < eventCount; i++) {
          for (int p = 0, c = 0; p < allData.length; p++) {
            if (compensatedNames.contains(paramNamesInOrder.get(p))) {
              writer.print(ext.formDeci(compensatedData[c++][i], 10));
            } else {
              writer.print(ext.formDeci(allData[p][i], 10));
            }
            if (p < allData.length - 1) {
              writer.print(",");
            }
          }
          writer.println();
        }

        break;
      case UNCOMPENSATED:
        for (int i = 0; i < eventCount; i++) {
          for (int p = 0; p < allData.length; p++) {
            writer.print(ext.formDeci(allData[p][i], 10));
            if (p < allData.length - 1) {
              writer.print(",");
            }
          }
          writer.println();
        }

        break;
      default:
        break;

    }
    writer.flush();
    writer.close();
  }

  private static DenseMatrix64F getInvertedSpilloverMatrix(double[][] coeffs) {
    DenseMatrix64F spillMatrix = new DenseMatrix64F(coeffs);
    CommonOps.invert(spillMatrix);
    return spillMatrix;
  }

  private static double[][] compensateSmall(ArrayList<String> dataColNames, double[][] data,
      String[] spillColNames, DenseMatrix64F spillMatrix) {
    int[] spillLookup = new int[dataColNames.size()];
    for (int i = 0; i < spillLookup.length; i++) {
      spillLookup[i] = ext.indexOfStr(dataColNames.get(i), spillColNames);
    }

    double[][] compensated = new double[data.length][];
    for (int r = 0; r < data.length; r++) {
      compensated[r] = new double[spillColNames.length];
      for (int c = 0, cS = 0; c < data[r].length; c++) {
        if (spillLookup[c] == -1) {
          continue;
        }
        float sum = 0f;
        for (int c1 = 0; c1 < data[r].length; c1++) {
          if (spillLookup[c1] == -1) {
            continue;
          }
          sum += data[r][c1] * spillMatrix.get(spillLookup[c1], spillLookup[c]);
        }
        compensated[r][cS++] = sum;
      }
    }

    return compensated;
  }

  //
  // private static float[][] compensate(ArrayList<String> dataColNames, float[][] data, String[]
  // spillColNames, DenseMatrix64F spillMatrix) {
  // int[] spillLookup = new int[dataColNames.size()];
  // for (int i = 0; i < spillLookup.length; i++) {
  // spillLookup[i] = ext.indexOfStr(dataColNames.get(i), spillColNames);
  // }
  //
  // float[][] compensated = new float[data.length][];
  // for (int r = 0; r < data.length; r++) {
  // compensated[r] = new float[data[r].length];
  // for (int c = 0; c < data[r].length; c++) {
  // if (spillLookup[c] == -1) {
  // compensated[r][c] = data[r][c];
  // continue;
  // }
  // float sum = 0f;
  // for (int c1 = 0; c1 < data[r].length; c1++) {
  // if (spillLookup[c1] == -1) continue;
  // sum += data[r][c1] * spillMatrix.get(spillLookup[c1], spillLookup[c]);
  // }
  // compensated[r][c] = sum;
  // }
  // }
  //
  // return compensated;
  // }

  public AXIS_SCALE getScaleForParam(String string) {
    AXIS_SCALE scale = paramScales.get(getInternalParamName(string));
    return scale == null ? AXIS_SCALE.LIN : scale;
  }

  public void setScaleForParam(String dataName, AXIS_SCALE scale) {
    paramScales.put(getInternalParamName(dataName), scale);
    paramTransforms.put(getInternalParamName(dataName), getDefaultTransform(scale));
  }

}
