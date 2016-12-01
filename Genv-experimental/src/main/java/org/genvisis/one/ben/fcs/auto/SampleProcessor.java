package org.genvisis.one.ben.fcs.auto;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
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
  
  default void updateGateParams(FCSDataLoader dataLoader, ArrayList<Gate> gates) {
  	for (Gate g : gates) {
  		String p = g.getXDimension().getParam();
  		String d = dataLoader.getInternalParamName(p);
  		if (!p.equals(d)) {
  			g.getXDimension().setParam(d);
  		}
  		p = g.getYDimension().getParam();
  		d = dataLoader.getInternalParamName(p);
  		if (!p.equals(d)) {
  			g.getYDimension().setParam(d);
  		}
  		updateGateParams(dataLoader, g.getChildGates());
  	}
  }
  
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
    updateGateParams(d, sn.gating.gateRoots);
    sn.gating.paramGateMap = GateFileUtils.parameterizeGates(sn.gating.gateMap);
  }
  
}


class LeafDataSampler extends AbstractSampleProcessor {
	private static final String FILE_EXT = ".data";
	private int sampleSize = 1000;
	private Map<String, PrintWriter> writers;
	private List<String> params;
	private String fileRoot;
	
	public LeafDataSampler(int sampleSize, String outputFileRoot, Map<String, PrintWriter> gateWriters, List<String> paramsInOrder) { 
		this.sampleSize = sampleSize;
		this.writers = gateWriters;
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
    	synchronized(writers) {
	    	writer = writers.get(g.getName());
	    	if (writer == null) {
	    		log.reportError("No output writer found for gate " + g.getName());
	    		String filename = fileRoot + ext.replaceWithLinuxSafeCharacters(g.getID() + "_" + g.getXDimension().getParam() + "_" + g.getYDimension().getParam() + FILE_EXT, false);
	    		writer = Files.getAppropriateWriter(filename);
	    		setupNewGateWriter(writer);
	    		writers.put(g.getName(), writer);
	    	}
    	}
      boolean[] incl = g.gate(d);
      int[] indices = Array.booleanArrayToIndices(incl);
      log.report("Extracting " + sampleSize + " data samples from " + Array.booleanArraySum(incl) + " possible data points in file " + sn.fcsFile + " for gate " + g.getName());
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
      }
      synchronized(writer) {
      	writer.flush();
      }
    }
    
    d.emptyAndReset();
    d = null;
  }
  
}

