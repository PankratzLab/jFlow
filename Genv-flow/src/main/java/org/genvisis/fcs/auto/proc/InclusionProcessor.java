package org.genvisis.fcs.auto.proc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class InclusionProcessor extends AbstractLoadingSampleProcessor {

  String outDir;
  final String ovvrDir;
  final String ovvrSfx;
  final String ovvrMatch;

  public InclusionProcessor(String o, String ovvrDir, String ovvrSuff, String ovvrMatch) {
    super();
    outDir = o;
    this.ovvrDir = ovvrDir;
    this.ovvrSfx = ovvrSuff;
    this.ovvrMatch = ovvrMatch;
  }

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    if (!Files.exists(sn.fcsFile)) {
      return;
    }
    loadPopsAndGates(sn);
    loadData(sn);

    String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
    String outFile = outDir + cleanedName + "/" + cleanedName + ".incl.xln.gz";

    if (Files.exists(outFile) && Files.countLines(outFile, 1) == d.getCount()) {
      log.reportTime("Found full output for sample " + sn.fcsFile);
      return;
    }

    if (ovvrDir != null && Files.exists(ovvrDir + ext.removeDirectoryInfo(sn.fcsFile) + ovvrSfx)) {
      d.loadGateOverrides(ovvrDir + ext.removeDirectoryInfo(sn.fcsFile) + ovvrSfx, ovvrMatch);
    }

    long t1 = System.nanoTime();
    HashSet<Gate> allGates1 = new HashSet<>(sn.gating.gateMap.values());
    ArrayList<Gate> gateList = new ArrayList<>();
    for (Gate g : allGates1) {
      gateList.add(g);
    }
    HashMap<String, boolean[]> incl = new HashMap<>();

    for (int i = 0; i < gateList.size(); i++) {
      Gate g = gateList.get(i);
      boolean[] gating = g.gate(d);
      String name = g.getFullNameAndGatingPath();
      name = getReplacedName(name);
      incl.put(name, gating);
    }

    new File(outDir + cleanedName + "/").mkdirs();

    PrintWriter writer = Files.getAppropriateWriter(outFile);
    writer.print("Index");
    for (int i = 0; i < gateList.size(); i++) {
      Gate g = gateList.get(i);
      String name = g.getFullNameAndGatingPath();
      name = getReplacedName(name);
      writer.print("\t" + name);
    }
    writer.println();

    for (int i = 0; i < d.getCount(); i++) {
      writer.print(Integer.toString(i));
      for (int gi = 0; gi < gateList.size(); gi++) {
        Gate g = gateList.get(gi);
        String name = g.getFullNameAndGatingPath();
        name = getReplacedName(name);
        writer.print("\t" + Boolean.toString(incl.get(name)[i]));
      }
      writer.println();
    }
    writer.close();
    log.reportTime("Processed " + sn.fcsFile + " in " + ext.getTimeElapsedNanos(t1));
  }
}
