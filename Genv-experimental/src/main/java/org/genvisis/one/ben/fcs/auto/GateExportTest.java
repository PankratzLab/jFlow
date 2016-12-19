package org.genvisis.one.ben.fcs.auto;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.PlainTextExport;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class GateExportTest {

//	private static final String WSP = "F:/Flow/gatingTest/";
//	private static final String FCS = "F:/Flow/gatingTest/2016-07-13_PANEL 1_ZF_Group one_F1631196_002.fcs";
	
	private static final String FCS = "/home/pankrat2/shared/flow/firstSetOf102/";
	private static final String WSP = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
	
	
	private void run() throws IOException {

//		Class<? extends SampleProcessor> processorClass = LeafDataSampler.class;
//		Class<? extends SampleProcessor> processorClass = PercentageWriter.class;

		SamplingPipeline sp = new SamplingPipeline(null, WSP, FCS, null, "/scratch.global/cole0482/FCS/", ProcClass.class);
		
		sp.run();
		
		while (!sp.checkFinished()) {
			Thread.yield();
		}
		
	}
	
	static class ProcClass extends AbstractSampleProcessor {
		Map<String, Integer> gateCoding;
		String outputDir;
		
		public ProcClass(Map<String, Integer> gateCoding, String outputDir) {
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
		Map<String, Integer> codingLookup;
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
	
	public static void main(String[] args) throws IOException {
		new GateExportTest().run();
	}
	
}
