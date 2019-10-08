package org.genvisis.fcs;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import org.genvisis.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.fcs.AbstractPanel2.AxisTransform;
import org.genvisis.jfcs.FCSKeywords;
import org.genvisis.jfcs.FCSReader;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;

public class FCSDataLoader {

  private static final String COMPENSATED_PREPEND = "Comp-";
  private static final int COMP_LEN = COMPENSATED_PREPEND.length();
  static final String GATING_KEY = "GATING";

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
  String[] presetGating;
  private long startLoadTime;
  private String loadedFile = null;
  int paramsCount = -1;
  Date lastModified;
  Date runDate;
  Date beginTime;
  Date endTime;
  double[] paramScaling;
  private FCSReader reader;

  public enum DATA_SET {
    ALL, COMPENSATED, UNCOMPENSATED;
  }

  public FCSDataLoader() {
    paramNamesInOrder = new ArrayList<>();
    paramShortNamesInOrder = new ArrayList<>();
    paramLongNamesInOrder = new ArrayList<>();
    paramScales = new HashMap<>();
    paramTransforms = new HashMap<>();
    ranges = new ArrayList<>();
    compensatedNames = new LinkedHashSet<>();
    compensatedIndices = new HashMap<>();
  }

  public void emptyAndReset() {
    setState(LOAD_STATE.UNLOADED);
    paramNamesInOrder = new ArrayList<>();
    paramShortNamesInOrder = new ArrayList<>();
    paramLongNamesInOrder = new ArrayList<>();
    paramScales = new HashMap<>();
    paramTransforms = new HashMap<>();
    compensatedNames = new LinkedHashSet<>();
    compensatedIndices = new HashMap<>();
    eventCount = -1;
    loadedCount = 0;
    presetGating = null;
    loadedFile = null;
    paramsCount = -1;
    paramScaling = null;
    runDate = null;
    beginTime = null;
    endTime = null;
    if (reader == null) {
      System.out.println("Reader was already null here!");
    }
    reader.dispose();
    reader = null;
    System.gc();
  }

  public String getLoadedFile() {
    return loadedFile;
  }

  private synchronized void setState(LOAD_STATE state) {
    this.state = state;
  }

  public synchronized LOAD_STATE getLoadState() {
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

  private void loadTimeData(FCSKeywords keys) {
    String dKey = "$DATE";
    String bKey = "$BTIM";
    String eKey = "$ETIM";

    String dVal = null, eVal = null, bVal = null;
    try {
      dVal = keys.getKeyword(dKey);
    } catch (Exception e) {}
    try {
      bVal = keys.getKeyword(bKey);
    } catch (Exception e) {}
    try {
      eVal = keys.getKeyword(eKey);
    } catch (Exception e) {}

    SimpleDateFormat sdfDate = new SimpleDateFormat("dd-MMM-yyyy");
    SimpleDateFormat sdfTime = new SimpleDateFormat("HH:mm:ss");
    if (dVal != null) {
      try {
        runDate = sdfDate.parse(dVal);
      } catch (ParseException e1) {
        System.out.println("Unable to parse date: " + dVal);
      }
    }
    if (bVal != null) {
      try {
        beginTime = sdfTime.parse(bVal);
      } catch (ParseException e1) {
        System.out.println("Unable to parse time: " + bVal);
      }
    }
    if (eVal != null) {
      try {
        endTime = sdfTime.parse(eVal);
      } catch (ParseException e1) {
        System.out.println("Unable to parse time: " + eVal);
      }
    }
  }

  public void loadData(String fcsFilename) throws IOException {
    // synchronized(this) {
    if (getLoadState() != LOAD_STATE.UNLOADED) {
      return;
    }
    setState(LOAD_STATE.LOADING);
    // }
    loadedFile = fcsFilename;

    startLoadTime = System.nanoTime();
    FCSReader reader = FCSReader.open(fcsFilename);

    FCSKeywords keys = reader.getKeywords();
    eventCount = keys.getEventCount();
    paramsCount = keys.getParameterCount();

    loadTimeData(keys);

    lastModified = null;// keys.getLastModified();
    if (lastModified == null) {
      lastModified = new Date(new File(fcsFilename).lastModified());
    }

    String gating = null;
    try {
      gating = keys.getKeyword(FCSDataLoader.GATING_KEY);
    } catch (Exception e) {
      // System.err.println("Info - no precreated gating info available.");
    }
    if (gating != null) {
      String[] sp = gating.split(",");
      presetGating = sp;
    } else {
      presetGating = null;
    }

    HashMap<String, String> names = new HashMap<>();
    for (int i = 0; i < paramsCount; i++) {
      String name = null;
      String shortName = null;
      try {
        name = keys.getParameterLongName(i);
      } catch (Exception e) {}
      try {
        shortName = keys.getParameterShortName(i);
      } catch (Exception e) {}
      // if (name == null) {
      // // TODO Error, no name set (included in spillover matrix?) -- will this scenario ever
      // happen?
      // name = "P" + i;
      // }
      String actName = shortName;
      if (name != null && !name.equals(shortName)) {
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
        scale = AXIS_SCALE.valueOf(keys.getKeyword(axisKeywork));
      } catch (Exception e) {
        // System.err.println("Warning - no axis scale set for parameter " +
        // paramNamesInOrder.get(i)
        // + "; assuming a linear scale.");
      }
      ;
      if (scale == AXIS_SCALE.LOG) {
        scale = AXIS_SCALE.BIEX;
      }
      paramScales.put(actName, scale);
      paramTransforms.put(actName, getDefaultTransform(scale));

      String rangeKeyword = "$P" + (i + 1) + "R";
      int rng = -1;
      try {
        rng = keys.getKeywordInt(rangeKeyword);
      } catch (Exception e) {
        System.err.println("Warning - no parameter range value for parameter "
                           + paramNamesInOrder.get(i) + "; assuming standard of 262144");
        rng = 262144;
      }
      ranges.add(rng);
    }

    String[] arr = reader.getSpillover().getParameterNames();
    for (int i = 0, count = arr.length; i < count; i++) {
      compensatedNames.add(names.get(arr[i]));
      compensatedIndices.put(names.get(arr[i]), i);
      paramScales.put(COMPENSATED_PREPEND + names.get(arr[i]), getScaleForParam(names.get(arr[i])));
      paramTransforms.put(COMPENSATED_PREPEND + names.get(arr[i]),
                          getDefaultTransform(getScaleForParam(names.get(arr[i]))));
    }

    this.reader = reader;
    setState(LOAD_STATE.LOADED);
  }

  private AxisTransform getDefaultTransform(AXIS_SCALE scale) {
    switch (scale) {
      case LOG:
        System.err.println("Error - NO DEFINED AXIS TRANSFORM FOR LOG AXES!");
      case BIEX:
        return AxisTransform.createBiexTransform();
      case LIN:
      default:
        return AxisTransform.createLinearTransform(0, 262144, 1);
    }
  }

  public AxisTransform getParamTransform(String param) {
    return paramTransforms.get(getInternalParamName(param));
  }

  public void setTransformMap(HashMap<String, AxisTransform> map) {
    HashMap<String, AxisTransform> newMap = new HashMap<>();
    for (String key : map.keySet()) {
      String nm = getInternalParamName(key);
      newMap.put(nm, map.get(key));
    }
    for (String key : this.paramTransforms.keySet()) {
      if (newMap.containsKey(key)) {
        continue;
      }
      newMap.put(key, this.paramTransforms.get(key));
    }
    this.paramTransforms = newMap;
  }

  public double[] getData(String colName, boolean waitIfNecessary) {
    boolean compensated = false;
    String columnName = (compensated = colName.startsWith(COMPENSATED_PREPEND)) ? COMPENSATED_PREPEND
                                                                                  + getInternalParamName(colName.substring(COMP_LEN))
                                                                                : getInternalParamName(colName);
    LOAD_STATE currState = getLoadState();
    double[] data;
    if (currState == LOAD_STATE.LOADED) {
      int index = paramNamesInOrder.indexOf(compensated ? columnName.substring(COMP_LEN)
                                                        : columnName);
      if (index < 0) {
        System.err.println("Error - gate index was " + index + " for " + colName + " | "
                           + columnName + " -- gates: " + getAllDisplayableNames(DATA_SET.ALL));
      }
      data = reader.getParamAsDoubles(paramShortNamesInOrder.get(index), compensated);
    } else {
      if (currState != LOAD_STATE.UNLOADED && waitIfNecessary) {
        while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
          Thread.yield();
        }
        data = reader.getParamAsDoubles(columnName, compensated);
      } else {
        int len = eventCount == -1 ? 0 : eventCount;
        return ArrayUtils.doubleArray(len, Double.NaN);
      }
    }
    return data;
  }

  public String getInternalParamName(String name) {
    String nm = name;
    boolean prepend = nm.startsWith(COMPENSATED_PREPEND);
    if (prepend) {
      nm = nm.substring(COMP_LEN);
    }
    int ind = paramNamesInOrder.indexOf(nm);
    if (ind != -1) {
      return prepend ? COMPENSATED_PREPEND + nm : nm;
    } else if ((ind = paramShortNamesInOrder.indexOf(nm)) != -1) {
      return prepend ? COMPENSATED_PREPEND + paramNamesInOrder.get(ind)
                     : paramNamesInOrder.get(ind);
    } else if ((ind = paramLongNamesInOrder.indexOf(nm)) != -1) {
      return prepend ? COMPENSATED_PREPEND + paramNamesInOrder.get(ind)
                     : paramNamesInOrder.get(ind);
    }
    return prepend ? COMPENSATED_PREPEND + nm : nm;
  }

  public String getPresetGateAssignment(int line) {
    if (presetGating == null || presetGating.length <= line || line < 0) {
      return "";
    }
    return presetGating[line];
  }

  public double[] getDataLine(List<String> params, int ind) {
    double[] line = new double[params.size()];
    double[] data = reader.getEventAsDoubles(ind, false);
    double[] compData = reader.getEventAsDoubles(ind, true);

    LOAD_STATE currState = getLoadState();
    if (currState != LOAD_STATE.UNLOADED) {
      while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
        Thread.yield();
      }
    }
    for (int i = 0; i < params.size(); i++) {
      String columnName = getInternalParamName(params.get(i));
      if (columnName.startsWith(COMPENSATED_PREPEND)) {
        line[i] = compData[compensatedIndices.get(columnName.substring(COMP_LEN))];
      } else {
        line[i] = data[paramNamesInOrder.indexOf(columnName)];
      }
    }
    return line;
  }

  public boolean containsParam(String paramName) {
    String nm = paramName;
    if (paramName.startsWith(COMPENSATED_PREPEND)) {
      nm = paramName.substring(COMP_LEN);
    }
    return (compensatedNames.contains(nm) || paramNamesInOrder.contains(nm))
           || paramShortNamesInOrder.contains(nm) || paramLongNamesInOrder.contains(nm);
  }

  public ArrayList<String> getAllDisplayableNames(DATA_SET set) {
    ArrayList<String> names = new ArrayList<>();
    switch (set) {
      case ALL:
        names = new ArrayList<>(paramNamesInOrder);
        for (int i = 0; i < names.size(); i++) {
          if (compensatedNames.contains(names.get(i))) {
            names.add(COMPENSATED_PREPEND + names.get(i));
          }
        }
        break;
      case COMPENSATED:
        names = new ArrayList<>(paramNamesInOrder);
        for (int i = 0; i < names.size(); i++) {
          if (compensatedNames.contains(names.get(i))) {
            names.set(i, COMPENSATED_PREPEND + names.get(i));
          }
        }
        break;
      case UNCOMPENSATED:
        names = new ArrayList<>(paramNamesInOrder);
        break;
      default:
        break;
    }
    return names;
  }

  public void exportData(String outputFilename, DATA_SET set) {}

  public AXIS_SCALE getScaleForParam(String string) {
    AXIS_SCALE scale = paramScales.get(getInternalParamName(string));
    return scale == null ? AXIS_SCALE.LIN : scale;
  }

  public void setScaleForParam(String dataName, AXIS_SCALE scale) {
    paramScales.put(getInternalParamName(dataName), scale);
    paramTransforms.put(getInternalParamName(dataName), getDefaultTransform(scale));
  }

  // tacked on functionality:

  private Map<String, boolean[]> gateOverride = new HashMap<>();
  private Map<String, List<String>> gateOverrideMatch = new HashMap<>();

  public void loadGateOverrides(String file, String match) {
    String[][] strData = HashVec.loadFileToStringMatrix(file, false, null);
    String[] hdr = strData[0];

    gateOverride = new HashMap<>();
    int count = strData.length - 1;
    for (int i = 0; i < hdr.length; i++) {
      gateOverride.put(hdr[i], new boolean[count]);
    }
    for (int i = 0; i < count; i++) {
      for (int h = 0; h < hdr.length; h++) {
        gateOverride.get(hdr[h])[i] = Boolean.parseBoolean(strData[i + 1][h]);
      }
    }

    gateOverrideMatch = new HashMap<>();
    try {
      BufferedReader reader = Files.getAppropriateReader(match);
      String l = null;
      while ((l = reader.readLine()) != null) {
        if (l.trim().equals("")) continue;
        String[] pts = l.trim().split("\t");
        String key = pts[0].trim();
        if (key.contains("(")) {
          key = key.substring(0, key.indexOf('(')).trim();
        }
        gateOverrideMatch.put(key, new ArrayList<>());
        for (int i = 1; i < pts.length; i++) {
          gateOverrideMatch.get(key).add(pts[i]);
        }
      }
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public boolean[] getOverrideGating(String gateName) {
    String key = gateName.contains("(") ? gateName.substring(0, gateName.indexOf('(')).trim()
                                        : gateName;
    if (gateOverrideMatch.containsKey(key)) {
      List<String> ovvr = gateOverrideMatch.get(key);
      boolean[] start = Arrays.copyOf(gateOverride.get(ovvr.get(0)), getCount());
      if (start == null) return null;
      for (int i = 1; i < ovvr.size(); i++) {
        if (gateOverride.get(ovvr.get(i)) == null) return null;
        start = ArrayUtils.booleanArrayAnd(start, gateOverride.get(ovvr.get(i)));
      }
      return start;
    } else if (gateOverride.containsKey(key)) {
      return gateOverride.get(key);
    } else
      return null;
  }

}
