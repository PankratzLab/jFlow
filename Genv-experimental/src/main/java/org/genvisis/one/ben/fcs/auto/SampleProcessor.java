package org.genvisis.one.ben.fcs.auto;

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

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.PlainTextExport;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.GateFileUtils;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.genvisis.one.ben.fcs.sub.EMInitializer;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

@FunctionalInterface
interface SampleProcessor {
  public void processSample(SampleNode sn, Logger log) throws IOException;
}

interface ProcessorFactory<T extends SampleProcessor> {
	T createProcessor(Object owner, int index);
	void cleanup(Object owner);
}

abstract class AbstractSampleProcessor implements SampleProcessor {
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
    while(d.getLoadState() != LOAD_STATE.LOADED) {
      try {
        Thread.sleep(100);
      } catch (InterruptedException e) {
      }
    }
    (new Logger()).reportTimeElapsed("Loaded FCS ... ", t1);
    d.setTransformMap(sn.savedTransforms);
    GateFileUtils.updateGateParams(d, sn.gating.gateRoots);
    sn.gating.paramGateMap = GateFileUtils.parameterizeGates(sn.gating.gateMap);
  }
  
}

class PercentageWriterFactory implements ProcessorFactory<PercentageWriter> {
	ConcurrentHashMap<Object, Map<String, Map<String, Double>>> resultMap = new ConcurrentHashMap<>();
	
	@Override
	public PercentageWriter createProcessor(Object owner, int index) {
		Map<String, Map<String, Double>> map;
		synchronized(resultMap) {
			map = resultMap.get(owner);
			if (map == null) {
				map = new ConcurrentHashMap<>();
				resultMap.put(owner, map);
			}
		}
		return new PercentageWriter(map);
	}
	
	@Override
	public void cleanup(Object owner) {
		Map<String, Map<String, Double>> resMap = resultMap.get(owner);
		if (resMap == null) {
			return;
		}
		PrintWriter writer = Files.getAppropriateWriter("F:/Flow/gatingTest/p1.pcts.xln");
		boolean header = false;
		for (Entry<String, Map<String, Double>> entry : resMap.entrySet()) {
			if (!header) {
	  		StringBuilder sb = new StringBuilder();
	  		for (String s : order) {
	  			sb.append("\t").append(s);
	  		}
	  		writer.println(sb.toString());
	  		header = true;
			}
			StringBuilder sb = new StringBuilder(entry.getKey());
			for (String s : order) {
				sb.append("\t").append(entry.getValue().get(s));
			}
			writer.println(sb.toString());
		}
		writer.flush();
		writer.close();
	}

	private static final String[] order =
																				{
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A)",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H)",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A)",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / IgD+ memory Bcells (CD27+) (Comp-BUV 737-A (IgD) v Comp-BB515-A (CD27))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / IgD- memory Bcells (CD27+) (Comp-BUV 737-A (IgD) v Comp-BB515-A (CD27))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / naive Bcells (CD27- IgD+) (Comp-BUV 737-A (IgD) v Comp-BB515-A (CD27))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / activated cytotoxic Tcells (CD8+ HLA-DR+) (Comp-PE-CF594-A (HLA-DR) v Comp-BUV 395-A (CD8))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / cytotoxic Tcells CD27- , CD28+ (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / pE cytotoxic Tcells (CD27-  CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / pE1 cytotoxic Tcells (CD27+  CD28+) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / pE2 cytotoxic Tcells (CD27+ , CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM1 cytotoxic Tcells (CD27+  CD28+) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM2 cytotoxic Tcells (CD27+  CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM3 cytotoxic Tcells (CD27-  CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM4 cytotoxic Tcells (CD27-  CD28+) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / activated helper Tcells (CD4+ HLA-DR+) (Comp-PE-CF594-A (HLA-DR) v Comp-APC-Cy7-A (CD4))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
																					"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",};
}

class PercentageWriter extends AbstractSampleProcessor {
	final Map<String, Map<String, Double>> pctMap;
	
	public PercentageWriter(Map<String, Map<String, Double>> pctMap) {
		this.pctMap = pctMap;
	}
	
	@Override
	public void processSample(SampleNode sn, Logger log) throws IOException {
    if (!Files.exists(sn.fcsFile)) { 
    	return;
    }
    loadPopsAndGates(sn);
    loadData(sn);
		
    Map<String, Double> pcts = pctMap.get(sn.fcsFile);
    if (pcts == null) {
    	pcts = new HashMap<>();
    	pctMap.put(sn.fcsFile, pcts);
    }
    HashSet<Gate> allGates = new HashSet<>(sn.gating.gateMap.values());
    for (Gate g : allGates) {
    	boolean[] gating = g.gate(d);
    	boolean[] parent = g.getParentGating(d);
    	int g1 = Array.booleanArraySum(gating);
    	int g2 = Array.booleanArraySum(parent);
    	pcts.put(g.getFullNameAndGatingPath(), ((double) g1) / (double) g2);
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
		
		synchronized(writers) {
			gateWriters = writers.get(owner);
			if (gateWriters == null) {
				gateWriters = new ConcurrentHashMap<>();
				writers.put(owner, gateWriters);
			}
		}
		synchronized(counts) {
			writeCounts = counts.get(owner);
			if (writeCounts == null) {
				writeCounts = new ConcurrentHashMap<>();
				counts.put(owner, writeCounts);
			}
		}
		synchronized(outRoots) {
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
			log.reportTime("Wrote " + counts.get(owner).get(entry.getKey()).get() + " for gate " + entry.getKey());
		}
	}
}

class LeafDataSampler extends AbstractSampleProcessor {
	private static final String FILE_EXT = ".data";
	private int sampleSize = 1000;
	private Map<String, PrintWriter> writers;
	private Map<String, AtomicInteger> writeCounts;
	private List<String> params;
	private String fileRoot;
	
	public LeafDataSampler(int sampleSize, String outputFileRoot, Map<String, PrintWriter> gateWriters, Map<String, AtomicInteger> writeCounts, List<String> paramsInOrder) { 
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
    	synchronized(writers) {
	    	writer = writers.get(g.getName());
	    	if (writer == null) {
	    		String filename = fileRoot + ext.replaceWithLinuxSafeCharacters(g.getName() + "_" + g.getXDimension().getParam() + "_" + g.getYDimension().getParam() + FILE_EXT, false);
	    		writer = Files.getAppropriateWriter(filename);
	    		setupNewGateWriter(writer);
	    		writers.put(g.getName(), writer);
	    		writeCounts.put(g.getName(), counter = new AtomicInteger());
	    	} else {
	    		counter = writeCounts.get(g.getName());
	    	}
    	}
      boolean[] incl = g.gate(d);
      int[] indices = Array.booleanArrayToIndices(incl);
      Random rand = new Random();
      for (int i = 0; i < sampleSize; i++) {
        int ind = indices[rand.nextInt(indices.length)];
        double[] line = d.getDataLine(params, ind);
        synchronized(writer) {
	        for (int l = 0; l < line.length; l++) {
	          writer.print(line[l]);
	          writer.print(l < line.length - 1 ? "\t" : "");
	        }
	        writer.println();
        }
        counter.incrementAndGet();
      }
      synchronized(writer) {
      	writer.flush();
      }
    }
    
    d.emptyAndReset();
    d = null;
  }


	static class GateAssignmentFactory implements ProcessorFactory<GateAssignmentProcessor> {
		
		ConcurrentHashMap<Object, ConcurrentHashMap<String, Integer>> ownerMaps = new ConcurrentHashMap<>();
		String outputDir = "";
		@Override
		public GateAssignmentProcessor createProcessor(Object owner, int index) {
			ConcurrentHashMap<String, Integer> map;
			synchronized(ownerMaps) {
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
	
	static class GateAssignmentProcessor extends AbstractSampleProcessor {
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
			
	    int[] coding = Array.intArray(d.getCount(), 0);
	    synchronized(gateCoding) {
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
			Files.writeArray(Array.toStringArray(coding), outputFile);
			PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(outputFile, false) + "_map.txt");
			for (Entry<String, Integer> look : codingLookup.entrySet()) {
				writer.println(look.getKey() + "=" + look.getValue());
			}
			writer.flush();
			writer.close();
		}
	}
	
  
}

