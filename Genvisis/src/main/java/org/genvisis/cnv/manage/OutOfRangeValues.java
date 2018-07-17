package org.genvisis.cnv.manage;

import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.Map.Entry;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;

public class OutOfRangeValues {

  public enum TYPE {
    X, Y, LRR;
  }

  /**
   * @param markerName
   * @return Map<String, Map<OutOfRangeValues.TYPE, Float>> map of sample names to out-of-range
   *         values
   */
  public Map<String, Map<OutOfRangeValues.TYPE, Float>> getMarkerOutliers(String markerName) {
    return mkrMap.containsKey(markerName) ? mkrMap.get(markerName) : new HashMap<>();
  }

  /**
   * @param sampleName
   * @return Map<String, Map<OutOfRangeValues.TYPE, Float>> map of marker names to out-of-range
   *         values
   */
  public Map<String, Map<OutOfRangeValues.TYPE, Float>> getSampleOutliers(String sampleName) {
    return smpMap.containsKey(sampleName) ? smpMap.get(sampleName) : new HashMap<>();
  }

  public boolean hasSample(String sampleName) {
    return smpMap.containsKey(sampleName);
  }

  public Hashtable<String, Float> getSampleOutliersForFile(Project proj, String sampleName) {
    Hashtable<String, Float> table = new Hashtable<>();
    Map<String, Marker> markerNameMap = proj.getMarkerSet().getMarkerNameMap();
    Map<Marker, Integer> markerIndexMap = proj.getMarkerSet().getMarkerIndexMap();
    Map<String, Map<OutOfRangeValues.TYPE, Float>> outs = getSampleOutliers(sampleName);
    for (Entry<String, Map<OutOfRangeValues.TYPE, Float>> e : outs.entrySet()) {
      for (Entry<OutOfRangeValues.TYPE, Float> e1 : e.getValue().entrySet()) {
        table.put(markerIndexMap.get(markerNameMap.get(e.getKey())).intValue() + "\t"
                  + e1.getKey().name().toLowerCase(), e1.getValue());
      }
    }
    return table;
  }

  public static OutOfRangeValues construct(Project proj) {
    Hashtable<String, Float> set = MarkerDataLoader.loadOutliers(proj);
    OutOfRangeValues outs = new OutOfRangeValues();
    String[] allMarkers = proj.getMarkerNames();
    for (Entry<String, Float> out : set.entrySet()) {
      String[] k = out.getKey().split("\t");
      int mInd = Integer.parseInt(k[0]);
      String sName = k[1];
      OutOfRangeValues.TYPE typ = TYPE.valueOf(k[2].toUpperCase());
      outs.add(allMarkers[mInd], sName, typ, out.getValue());
    }
    return outs;
  }

  private OutOfRangeValues() {
    mkrMap = new HashMap<>();
    smpMap = new HashMap<>();
  }

  Map<String, Map<String, Map<OutOfRangeValues.TYPE, Float>>> mkrMap;
  Map<String, Map<String, Map<OutOfRangeValues.TYPE, Float>>> smpMap;

  private void add(String mkr, String samp, OutOfRangeValues.TYPE typ, Float value) {
    if (!mkrMap.containsKey(mkr)) {
      mkrMap.put(mkr, new HashMap<>());
    }
    Map<String, Map<OutOfRangeValues.TYPE, Float>> outMap = mkrMap.get(mkr);
    if (!outMap.containsKey(samp)) {
      outMap.put(samp, new HashMap<>());
    }
    outMap.get(samp).put(typ, value);

    if (!smpMap.containsKey(samp)) {
      smpMap.put(samp, new HashMap<>());
    }
    outMap = smpMap.get(samp);
    if (!outMap.containsKey(mkr)) {
      outMap.put(mkr, new HashMap<>());
    }
    outMap.get(mkr).put(typ, value);
  }

}
