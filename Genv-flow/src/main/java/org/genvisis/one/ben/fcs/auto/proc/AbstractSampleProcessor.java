package org.genvisis.one.ben.fcs.auto.proc;

import java.io.IOException;
import java.util.HashMap;
import org.genvisis.common.Logger;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.gating.GateFileUtils;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public abstract class AbstractSampleProcessor implements SampleProcessor {

  NodeList popList;
  HashMap<String, Element> popMap = new HashMap<>();
  HashMap<String, Element> gateMap = new HashMap<>();
  FCSDataLoader d;

  protected AbstractSampleProcessor() {}

  void loadPopsAndGates(SampleNode sn) {
    popList = sn.sampleNode.getElementsByTagName("Population");
    // annoOffsetX/annoOffsetY ???
    // panelState ???
    for (int i = 0; i < popList.getLength(); i++) {
      Element e = (Element) popList.item(i);
      Element g = (Element) e.getElementsByTagName("Gate").item(0);
      popMap.put(g.getAttribute("gating:id"), e);
      gateMap.put(g.getAttribute("gating:id"), g);
    }
  }

  void loadData(SampleNode sn) throws IOException {
    long t1 = System.currentTimeMillis();
    d = new FCSDataLoader();
    d.loadData(sn.fcsFile);
    (new Logger()).reportTimeElapsed("Loaded FCS ... ", t1);
    d.setTransformMap(sn.savedTransforms);
    GateFileUtils.updateGateParams(d, sn.gating.gateRoots);
    sn.gating.paramGateMap = GateFileUtils.parameterizeGates(sn.gating.gateMap);
  }

}
