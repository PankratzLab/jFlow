package one.ben.fcs.gating;

import java.util.ArrayList;
import java.util.HashMap;

public class GatingStrategy {
    HashMap<String, Gate> gateMap; // ID -> Gate
    HashMap<String, ArrayList<Gate>> paramGateMap; // paramName -> list of applicable Gates 
    ArrayList<Gate> gateRoots;
    
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

    public ArrayList<Gate> getGatesForParams(String xCol, String yCol) {
        ArrayList<Gate> gates = paramGateMap.get(xCol);
        if (gates == null) return new ArrayList<Gate>();
        ArrayList<Gate> ret = new ArrayList<Gate>();
        gt : for (Gate g : gates) {
            if (g.dimensions.size() == 2) {
                for (GateDimension gd : g.dimensions) {
                    if (!gd.paramName.equals(yCol) && !gd.paramName.equals(xCol)) {
                        continue gt;
                    }
                }
                ret.add(g);
            }
        }
        return ret;
    }
    
}
