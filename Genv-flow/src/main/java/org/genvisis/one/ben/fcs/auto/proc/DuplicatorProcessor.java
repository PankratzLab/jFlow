package org.genvisis.one.ben.fcs.auto.proc;

import java.io.IOException;
import java.util.HashSet;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.FCSFileDuplicator;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

class DuplicatorProcessor extends AbstractSampleProcessor {

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    loadPopsAndGates(sn);
    loadData(sn);

    String[] gatings = ArrayUtils.stringArray(d.getCount(), "");

    HashSet<Gate> leaves = sn.gating.getAllLeafGates();
    for (Gate g : leaves) {
      boolean[] b = g.gate(d);
      for (int i = 0; i < b.length; i++) {
        if (b[i]) {
          if ("".equals(gatings[i])) {
            gatings[i] = g.getName();
          } else {
            int g1 = countDepth(g);
            int g2 = countDepth(sn.gating.gateMap.get(gatings[i]));
            if (g1 > g2) {
              gatings[i] = g.getName();
            }
          }
        }
      }
    }

    String file = sn.fcsFile;
    String newFile = ext.rootOf(file, false) + "_gated.fcs";
    FCSFileDuplicator.createFrom(file, newFile, log, gatings);

  }

  private int countDepth(Gate g) {
    int cnt = 0;
    Gate g1 = g;
    while (g1.getParentGate() != null) {
      cnt++;
      g1 = g1.getParentGate();
    }
    return cnt;
  }
}
