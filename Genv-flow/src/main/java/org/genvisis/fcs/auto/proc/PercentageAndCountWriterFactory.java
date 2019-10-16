package org.genvisis.fcs.auto.proc;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import org.pankratzlab.common.Files;

public class PercentageAndCountWriterFactory implements ProcessorFactory<PercentageAndCountWriter> {

  ConcurrentHashMap<Object, Map<String, Map<String, Double>>> resultMapPcts = new ConcurrentHashMap<>();
  ConcurrentHashMap<Object, Map<String, Map<String, Integer>>> resultMapCnts = new ConcurrentHashMap<>();
  String outputDir;
  int panel;
  final String ovvrDir;
  final String ovvrSfx;
  final String ovvrMatch;

  public PercentageAndCountWriterFactory(String outDir, int panel, String ovvrDir, String ovvrSfx,
                                         String ovvrMatch) {
    super();
    this.outputDir = outDir;
    this.panel = panel;
    this.ovvrDir = ovvrDir;
    this.ovvrSfx = ovvrSfx;
    this.ovvrMatch = ovvrMatch;
  }

  @Override
  public PercentageAndCountWriter createProcessor(Object owner, int index) {
    Map<String, Map<String, Double>> map;
    Map<String, Map<String, Integer>> map2;
    synchronized (resultMapPcts) {
      map = resultMapPcts.get(owner);
      if (map == null) {
        map = new ConcurrentHashMap<>();
        resultMapPcts.put(owner, map);
      }
      map2 = resultMapCnts.get(owner);
      if (map2 == null) {
        map2 = new ConcurrentHashMap<>();
        resultMapCnts.put(owner, map2);
      }
    }
    return new PercentageAndCountWriter(map, map2, ovvrDir, ovvrSfx, ovvrMatch);
  }

  private <T> void writeMap(Map<String, Map<String, T>> resultMap, String fileName) {
    PrintWriter writer = Files.getAppropriateWriter(outputDir + fileName);
    boolean header = false;
    HashSet<String> headerSet = new HashSet<>();
    for (Entry<String, Map<String, T>> entry : resultMap.entrySet()) {
      headerSet.addAll(entry.getValue().keySet());
    }
    ArrayList<String> headers = new ArrayList<>(headerSet);
    Collections.sort(headers, new Comparator<String>() {

      @Override
      public int compare(String o1, String o2) {
        return new Integer(o1.length()).compareTo(o2.length());
      }
    });
    for (Entry<String, Map<String, T>> entry : resultMap.entrySet()) {
      if (!header) {
        StringBuilder sb = new StringBuilder();
        for (String s : headers) {
          sb.append("\t").append(s);
        }
        writer.println(sb.toString());
        header = true;
      }
      StringBuilder sb = new StringBuilder(entry.getKey());
      for (String s : headers) {
        sb.append("\t").append(entry.getValue().get(s));
      }
      writer.println(sb.toString());
    }
    writer.flush();
    writer.close();
  }

  @Override
  public void cleanup(Object owner) {
    Map<String, Map<String, Double>> pctMap = resultMapPcts.get(owner);
    Map<String, Map<String, Integer>> cntMap = resultMapCnts.get(owner);
    if (pctMap != null) {
      writeMap(pctMap, "p" + panel + ".pcts.xln");
    }
    if (cntMap != null) {
      writeMap(cntMap, "p" + panel + ".cnts.xln");
    }
  }
}
