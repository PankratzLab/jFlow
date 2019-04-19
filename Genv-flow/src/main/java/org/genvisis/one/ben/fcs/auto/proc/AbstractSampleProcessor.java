package org.genvisis.one.ben.fcs.auto.proc;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.gating.GateFileUtils;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.Logger;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public abstract class AbstractSampleProcessor implements SampleProcessor {

  NodeList popList;
  HashMap<String, Element> popMap = new HashMap<>();
  HashMap<String, Element> gateMap = new HashMap<>();
  FCSDataLoader d;

  static final Map<String, String> dimSwitch = new HashMap<>();

  {
    dimSwitch.put("Comp-BV 605-A \\(CD95\\)", "Comp-BV605-A \\(CD95\\)");
    dimSwitch.put("Comp-BV 510-A \\(CD28\\)", "Comp-BV510-A \\(CD28\\)");
    dimSwitch.put("Comp-BB 515-A \\(CD27\\)", "Comp-BB515-A \\(CD27\\)");
    dimSwitch.put("Comp-BB515-A \\(CD27\\)", "Comp-FITC-A \\(CD27\\)");
    dimSwitch.put("Comp-BV 421-A \\(CCR7\\)", "Comp-BV421-A \\(CCR7\\)");
    dimSwitch.put("Comp-BV 711-A \\(CD45RA\\)", "Comp-BV711-A \\(CD45RA\\)");
    dimSwitch.put("Comp-BUV 395-A \\(CD8\\)", "Comp-BUV396-A \\(CD8\\)");
    dimSwitch.put("LIVE/DEAD", "L/D");
    dimSwitch.put("Comp-PE-Cy7 \\(blue\\)-A \\(CD19\\)", "Comp-PE-Cy7-A \\(CD19\\)");
    dimSwitch.put("Comp-BUV 737-A \\(IgD\\)", "Comp-BUV737-A \\(IgD\\)");
  }

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
