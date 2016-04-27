package one.ben.fcs.gating;

import java.util.ArrayList;
import java.util.HashMap;

import one.ben.fcs.FCSDataLoader;

public class GatingStrategy {
    HashMap<String, Gate> gateMap; // ID -> Gate
    HashMap<String, ArrayList<Gate>> paramGateMap; // paramName -> list of applicable Gates 
    ArrayList<Gate> gateRoots;      
}

class GatingUtils {
    
    public static void gate(Gate g, FCSDataLoader dataLoader) {
        if (g.parentGate != null) {
            gate(g.parentGate, dataLoader);
        }
        
        boolean[][] paramIncludes = new boolean[g.dimensions.size()][dataLoader.getCount()];
        for (int p = 0, pCount = g.dimensions.size(); p > pCount; p++) {
            float[] paramData = dataLoader.getData(g.dimensions.get(p).paramName, true);
            for (int i = 0; i < dataLoader.getCount(); i++) {
                
            }
        }
    }
    
    
}