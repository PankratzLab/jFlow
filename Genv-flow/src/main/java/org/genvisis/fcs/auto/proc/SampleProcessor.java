package org.genvisis.fcs.auto.proc;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.genvisis.fcs.sub.EMInitializer;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.PlainTextExport;

@FunctionalInterface
public interface SampleProcessor {

  public void processSample(SampleNode sn, Logger log) throws IOException;

}

abstract class AbstractSampleProcessor implements SampleProcessor {

  private final Map<String, String> dimSwitch = new HashMap<>();

  protected AbstractSampleProcessor() {}

  public void addDimensionNameOverride(String from, String to) {
    dimSwitch.put(from, to);
  }

  public String replaceName(String name) {
    return dimSwitch.containsKey(name) ? dimSwitch.get(name) : name;
  }

  public String getReplacedName(String name) {
    String newName = name;
    for (Entry<String, String> dim : dimSwitch.entrySet()) {
      newName = newName.replaceAll(dim.getKey(), dim.getValue());
    }
    return newName;
  }

}

class PercentageAndCountWriter extends AbstractLoadingSampleProcessor {

  final Map<String, Map<String, Double>> pctMap;
  final Map<String, Map<String, Integer>> cntMap;
  final String ovvrDir;
  final String ovvrSfx;
  final String ovvrMatch;

  public PercentageAndCountWriter(Map<String, Map<String, Double>> pctMap,
                                  Map<String, Map<String, Integer>> cntMap, String ovvrDir,
                                  String ovvrSfx, String ovvrMatch) {
    this.pctMap = pctMap;
    this.cntMap = cntMap;
    this.ovvrDir = ovvrDir;
    this.ovvrSfx = ovvrSfx;
    this.ovvrMatch = ovvrMatch;
  }

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    if (!Files.exists(sn.fcsFile)) {
      return;
    }
    loadPopsAndGates(sn);
    loadData(sn);

    if (ovvrDir != null && Files.exists(ovvrDir + ext.removeDirectoryInfo(sn.fcsFile) + ovvrSfx)) {
      d.loadGateOverrides(ovvrDir + ext.removeDirectoryInfo(sn.fcsFile) + ovvrSfx, ovvrMatch);
    }

    Map<String, Double> pcts = pctMap.get(sn.fcsFile);
    if (pcts == null) {
      pcts = new ConcurrentHashMap<>();
      pctMap.put(sn.fcsFile, pcts);
    }
    Map<String, Integer> cnts = cntMap.get(sn.fcsFile);
    if (cnts == null) {
      cnts = new ConcurrentHashMap<>();
      cntMap.put(sn.fcsFile, cnts);
    }
    HashSet<Gate> allGates = new HashSet<>(sn.gating.gateMap.values());
    for (Gate g : allGates) {
      boolean[] gating = g.gate(d);
      boolean[] parent = g.getParentGating(d);
      int g1 = ArrayUtils.booleanArraySum(gating);
      int g2 = ArrayUtils.booleanArraySum(parent);

      String name = g.getFullNameAndGatingPath();
      name = getReplacedName(name);
      pcts.put(name, 100 * ((double) g1) / (double) g2);
      cnts.put(name, g1);
    }

    d.emptyAndReset();
    d = null;
  }
}

class LeafDataSamplerFactory implements ProcessorFactory<LeafDataSampler> {

  final List<String> params;
  static final int SAMPLING = 500;
  final String outDir;
  ConcurrentHashMap<Object, Map<String, PrintWriter>> writers;
  ConcurrentHashMap<Object, Map<String, AtomicInteger>> counts;
  ConcurrentHashMap<Object, String> outRoots;
  Logger log;

  public LeafDataSamplerFactory(String outputDir, Logger log) {
    this.params = new ArrayList<>();
    for (String s : EMInitializer.DATA_COLUMNS) {
      params.add(s);
    }
    writers = new ConcurrentHashMap<>();
    counts = new ConcurrentHashMap<>();
    outRoots = new ConcurrentHashMap<>();
    outDir = outputDir;
    this.log = log;
  }

  @Override
  public LeafDataSampler createProcessor(Object owner, int index) {
    Map<String, PrintWriter> gateWriters;
    Map<String, AtomicInteger> writeCounts;
    String outputFileRoot;

    synchronized (writers) {
      gateWriters = writers.get(owner);
      if (gateWriters == null) {
        gateWriters = new ConcurrentHashMap<>();
        writers.put(owner, gateWriters);
      }
    }
    synchronized (counts) {
      writeCounts = counts.get(owner);
      if (writeCounts == null) {
        writeCounts = new ConcurrentHashMap<>();
        counts.put(owner, writeCounts);
      }
    }
    synchronized (outRoots) {
      outputFileRoot = outRoots.get(owner);
      if (outputFileRoot == null) {
        outputFileRoot = outDir + "p" + index + ".";
        outRoots.put(owner, outputFileRoot);
      }
    }

    return new LeafDataSampler(SAMPLING, outputFileRoot, gateWriters, writeCounts, params);
  }

  @Override
  public void cleanup(Object owner) {
    for (Entry<String, PrintWriter> entry : writers.get(owner).entrySet()) {
      entry.getValue().flush();
      entry.getValue().close();
      log.reportTime("Wrote " + counts.get(owner).get(entry.getKey()).get() + " for gate "
                     + entry.getKey());
    }
  }
}

class LeafDataSampler extends AbstractLoadingSampleProcessor {

  private static final String FILE_EXT = ".data";
  private int sampleSize = 1000;
  private Map<String, PrintWriter> writers;
  private Map<String, AtomicInteger> writeCounts;
  private List<String> params;
  private String fileRoot;

  public LeafDataSampler(int sampleSize, String outputFileRoot,
                         Map<String, PrintWriter> gateWriters,
                         Map<String, AtomicInteger> writeCounts, List<String> paramsInOrder) {
    this.sampleSize = sampleSize;
    this.writers = gateWriters;
    this.writeCounts = writeCounts;
    this.params = paramsInOrder;
    this.fileRoot = outputFileRoot;
  }

  private void setupNewGateWriter(PrintWriter writer) {
    for (int i = 0; i < params.size(); i++) {
      writer.print(params.get(i));
      if (i < params.size() - 1) {
        writer.print("\t");
      }
    }
    writer.println();
    writer.flush();
  }

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    if (!Files.exists(sn.fcsFile)) {
      return;
    }
    loadPopsAndGates(sn);
    loadData(sn);

    HashSet<Gate> leafGates = sn.gating.getAllLeafGates();
    HashSet<Gate> map = new HashSet<>();
    for (Gate g : leafGates) {
      map.add(g);
    }
    for (Gate g : map) {
      PrintWriter writer;
      AtomicInteger counter;
      synchronized (writers) {
        writer = writers.get(g.getName());
        if (writer == null) {
          String filename = fileRoot
                            + ext.replaceWithLinuxSafeCharacters(g.getName() + "_"
                                                                 + g.getXDimension().getParam()
                                                                 + "_"
                                                                 + g.getYDimension().getParam()
                                                                 + FILE_EXT, false);
          writer = Files.getAppropriateWriter(filename);
          setupNewGateWriter(writer);
          writers.put(g.getName(), writer);
          writeCounts.put(g.getName(), counter = new AtomicInteger());
        } else {
          counter = writeCounts.get(g.getName());
        }
      }
      boolean[] incl = g.gate(d);
      int[] indices = ArrayUtils.booleanArrayToIndices(incl);
      Random rand = new Random();
      for (int i = 0; i < sampleSize; i++) {
        int ind = indices[rand.nextInt(indices.length)];
        double[] line = d.getDataLine(params, ind);
        synchronized (writer) {
          for (int l = 0; l < line.length; l++) {
            writer.print(line[l]);
            writer.print(l < line.length - 1 ? "\t" : "");
          }
          writer.println();
        }
        counter.incrementAndGet();
      }
      synchronized (writer) {
        writer.flush();
      }
    }

    d.emptyAndReset();
    d = null;
  }

  static class GateAssignmentFactory implements ProcessorFactory<GateAssignmentProcessor> {

    ConcurrentHashMap<Object, ConcurrentHashMap<String, Integer>> ownerMaps = new ConcurrentHashMap<>();
    String outputDir = "/scratch.global/cole0482/FCS/";

    @Override
    public GateAssignmentProcessor createProcessor(Object owner, int index) {
      ConcurrentHashMap<String, Integer> map;
      synchronized (ownerMaps) {
        map = ownerMaps.get(owner);
        if (map == null) {
          map = new ConcurrentHashMap<>();
          ownerMaps.put(owner, map);
        }
      }
      return new GateAssignmentProcessor(map, outputDir);
    }

    @Override
    public void cleanup(Object owner) {
      // unused
    }

  }

  static class GateAssignmentProcessor extends AbstractLoadingSampleProcessor {

    Map<String, Integer> gateCoding;
    String outputDir;

    public GateAssignmentProcessor(Map<String, Integer> gateCoding, String outputDir) {
      this.gateCoding = gateCoding;
      this.outputDir = outputDir;
    }

    private void setup(SampleNode sn) {
      for (String s : sn.gating.getAllGateNames()) {
        Gate g = sn.gating.gateMap.get(s);
        Integer code = gateCoding.get(g.getName());
        if (code == null) {
          code = 1;
          Gate g1 = g;
          while (g1 != null) {
            g1 = g1.getParentGate();
            if (g1 != null) {
              code += g1.getChildGates().size();
            }
          }
          while (gateCoding.containsValue(code)) {
            code++;
          }
          gateCoding.put(g.getName(), code);
        }
      }
    }

    @Override
    public void processSample(SampleNode sn, Logger log) throws IOException {
      if (!Files.exists(sn.fcsFile)) {
        return;
      }
      loadPopsAndGates(sn);
      loadData(sn);

      int[] coding = ArrayUtils.intArray(d.getCount(), 0);
      synchronized (gateCoding) {
        setup(sn);
      }

      HashSet<Gate> leaves = sn.gating.getAllLeafGates();
      for (Gate g : leaves) {
        assign(sn, g, coding);
      }

      CodingData cd = new CodingData();
      cd.codingLookup = gateCoding;
      cd.coding = coding;
      cd.file = sn.fcsFile;
      SerializedFiles.writeSerial(cd, outputDir + ext.rootOf(sn.fcsFile, true) + "_coding.ser");
      cd.exportToText(outputDir + ext.rootOf(sn.fcsFile, true) + "_coding.xln", log);
    }

    private void assign(SampleNode sn, Gate g, int[] coding) {
      boolean[] gating = g.gate(d);
      Integer code = gateCoding.get(g.getName());
      if (code == null) {
        code = gateCoding.size() + 1;
        gateCoding.put(g.getName(), code);
      }
      for (int i = 0; i < coding.length; i++) {
        if (gating[i]) {
          if (coding[i] == 0) {
            coding[i] = code.intValue();
          } else {
            int c = coding[i];
            Gate g2 = null;
            for (Entry<String, Integer> ent : gateCoding.entrySet()) {
              if (ent.getValue() == c) {
                g2 = sn.gating.gateMap.get(ent.getKey());
              }
            }
            if (g2 == null || g.getGateTreeLevel() > g2.getGateTreeLevel()) {
              coding[i] = code.intValue();
            }
          }
        }
      }
      if (g.getParentGate() != null) {
        assign(sn, g.getParentGate(), coding);
      }
    }

  }

  static class CodingData implements Serializable, PlainTextExport {

    private Map<String, Integer> codingLookup;
    int[] coding;
    String file;

    @Override
    public void exportToText(String outputFile, Logger log) {
      Files.writeArray(ArrayUtils.toStringArray(coding), outputFile);
      PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(outputFile, false) + "_map.txt");
      for (Entry<String, Integer> look : codingLookup.entrySet()) {
        writer.println(look.getKey() + "=" + look.getValue());
      }
      writer.flush();
      writer.close();
    }
  }

}
