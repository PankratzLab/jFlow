package org.genvisis.one.ben.fcs.auto.proc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class InclusionProcessor extends AbstractSampleProcessor {

  String outDir;
  final String ovvrDir;
  final String ovvrSfx;
  final String ovvrMatch;
  static final Map<String, String> dimSwitch = new HashMap<>();

  {
    dimSwitch.put("Comp-BV 605-A (CD95)", "Comp-BV605-A (CD95)");
    dimSwitch.put("Comp-BV 510-A (CD28)", "Comp-BV510-A (CD28)");
    dimSwitch.put("Comp-BB 515-A (CD27)", "Comp-BB515-A (CD27)");
    dimSwitch.put("Comp-BB515-A (CD27)", "Comp-FITC-A (CD27)");
    dimSwitch.put("Comp-BV 421-A (CCR7)", "Comp-BV421-A (CCR7)");
    dimSwitch.put("Comp-BV 711-A (CD45RA)", "Comp-BV711-A (CD45RA)");
  }

  public InclusionProcessor(String o, String ovvrDir, String ovvrSuff, String ovvrMatch) {
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

    if (ovvrDir != null) {
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
      for (Entry<String, String> dim : dimSwitch.entrySet()) {
        name.replaceAll(dim.getKey(), dim.getValue());
      }
      incl.put(name, gating);
    }

    String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
    String outFile = outDir + cleanedName + "/" + cleanedName + ".incl.xln.gz";
    new File(outDir + cleanedName + "/").mkdirs();

    PrintWriter writer = Files.getAppropriateWriter(outFile);
    writer.print("Index");
    for (int i = 0; i < gateList.size(); i++) {
      Gate g = gateList.get(i);
      String name = g.getFullNameAndGatingPath();
      for (Entry<String, String> dim : dimSwitch.entrySet()) {
        name.replaceAll(dim.getKey(), dim.getValue());
      }
      writer.print("\t" + name);
    }
    writer.println();

    for (int i = 0; i < d.getCount(); i++) {
      writer.print(Integer.toString(i));
      for (int gi = 0; gi < gateList.size(); gi++) {
        Gate g = gateList.get(gi);
        String name = g.getFullNameAndGatingPath();
        for (Entry<String, String> dim : dimSwitch.entrySet()) {
          name.replaceAll(dim.getKey(), dim.getValue());
        }
        writer.print("\t" + Boolean.toString(incl.get(name)[i]));
      }
      writer.println();
    }
    writer.close();
    log.reportTime("Processed " + sn.fcsFile + " in " + ext.getTimeElapsedNanos(t1));
  }
}
