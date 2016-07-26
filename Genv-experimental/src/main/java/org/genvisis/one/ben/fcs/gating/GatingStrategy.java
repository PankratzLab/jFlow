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
    
    public ArrayList<Gate> getGatesForParam(String paramName) {
        ArrayList<Gate> gates = paramGateMap.get(paramName);
        if (gates == null) gates = new ArrayList<Gate>();
        return gates;
    }
    
    public ArrayList<Gate> getGatesForParamOnly(String paramName) {
        ArrayList<Gate> gates = getGatesForParam(paramName);
        if (gates == null) return new ArrayList<Gate>();
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
        allNames.add(g.getName());
    }

    @SuppressWarnings("unchecked")
    public ArrayList<Gate> getRootGates() {
        return (ArrayList<Gate>) gateRoots.clone();
    }

    public String getFile() {
        return fileName;
    }

    public void setFile(String filename2) {
        this.fileName = filename2;
    }
    
}