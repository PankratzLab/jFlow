package org.genvisis.fcs.auto.proc;

import java.io.IOException;
import java.util.HashMap;

import org.genvisis.fcs.FCSDataLoader;
import org.genvisis.fcs.gating.GateFileUtils;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.Logger;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public abstract class AbstractLoadingSampleProcessor extends AbstractSampleProcessor {

  NodeList popList;
  HashMap<String, Element> popMap = new HashMap<>();
  HashMap<String, Element> gateMap = new HashMap<>();
  FCSDataLoader d;

  protected AbstractLoadingSampleProcessor(String dimensionOverrideFile) {
    super();
    loadDimOverrides(dimensionOverrideFile);
  }

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
