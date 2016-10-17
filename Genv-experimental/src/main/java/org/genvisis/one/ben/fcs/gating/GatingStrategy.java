package org.genvisis.one.ben.fcs.gating;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class GatingStrategy {
	HashSet<String> allNames = new HashSet<String>();
	HashMap<String, Gate> gateMap = new HashMap<String, Gate>(); // ID -> Gate
	HashMap<String, ArrayList<Gate>> paramGateMap = new HashMap<String, ArrayList<Gate>>(); // paramName -> list of applicable Gates
	ArrayList<Gate> gateRoots = new ArrayList<Gate>();

	private String fileName;

    private void addGateToLookups(Gate g2) {
      allNames.add(g2.getName() == null || "".equals(g2.getName()) ? g2.getID() : g2.getName());
      gateMap.put(g2.getID(), g2);
      if (g2.getName() != null && !g2.getName().equals("")) {
        gateMap.put(g2.getName(), g2);
      }
      for (GateDimension gd : g2.dimensions) {
        ArrayList<Gate> gates = paramGateMap.get(gd.paramName);
        if (gates == null) {
          gates = new ArrayList<Gate>();
          paramGateMap.put(gd.paramName, gates);
        }
        gates.add(g2);
      }
      
      for (Gate c : g2.children) {
        addGateToLookups(c);
      }
    }
	
	public GatingStrategy copy(String newFileName) { // create copy where gates are same params but different IDs
	  GatingStrategy gs = new GatingStrategy();
	  gs.fileName = newFileName;
	  
	  for (Gate g : gateRoots) {
	    Gate g2 = g.copy(null);
	    gs.addGateToLookups(g2);
	    gs.gateRoots.add(g);
	  }
	  
	  return gs;
	}
	
	public ArrayList<Gate> getGatesForParam(String paramName) {
		ArrayList<Gate> gates = paramGateMap.get(paramName);
		if (gates == null) {
			gates = new ArrayList<Gate>();
		}
		return gates;
	}

	public ArrayList<Gate> getGatesForParamOnly(String paramName) {
		ArrayList<Gate> gates = getGatesForParam(paramName);
		if (gates == null) {
			return new ArrayList<Gate>();
		}
		ArrayList<Gate> ret = new ArrayList<Gate>();
		for (Gate g : gates) {
			if (g.dimensions.size() == 1) {
				ret.add(g);
			}
		}
		return ret;
	}

	public boolean gateNameExists(String name) {
		return allNames.contains(name);
	}

	public void addRootGate(Gate g) {
		gateRoots.add(g);
		for (GateDimension gd : g.dimensions) {
			ArrayList<Gate> gates = paramGateMap.get(gd.paramName);
			if (gates == null) {
				gates = new ArrayList<Gate>();
				paramGateMap.put(gd.paramName, gates);
			}
			gates.add(g);
		}
		gateMap.put(g.getID(), g);
		if (g.getName() != null && !g.getName().equals("")) {
		  gateMap.put(g.getName(), g);
		}
		allNames.add(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
	}

	public ArrayList<Gate> getRootGates() {
		return gateRoots;
	}


	public void deleteGate(Gate g) {
		deleteGate(g, gateRoots);
	}

	private boolean deleteGate(Gate g, ArrayList<Gate> gates) {
		if (gates.contains(g)) {
			gates.remove(g);
			allNames.remove(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
			gateMap.remove(g.getID());
	        if (g.getName() != null && !g.getName().equals("")) {
	          gateMap.remove(g.getName());
	        }
	        for (Gate c : g.children) {
	          deleteGate(c);
	        }
			// TODO remove from paramGateMap?
			return true;
		} else {
			for (Gate g2 : gates) {
				return deleteGate(g, g2.getChildGates());
			}
		}
		return false;
	}

	public String getFile() {
		return fileName;
	}

	public void setFile(String filename2) {
		fileName = filename2;
	}

}
