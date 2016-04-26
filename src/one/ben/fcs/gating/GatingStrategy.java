package one.ben.fcs.gating;

import java.util.ArrayList;
import java.util.HashMap;

public class GatingStrategy {
    HashMap<String, Gate> gateMap; // ID -> Gate
    HashMap<String, ArrayList<Gate>> paramGateMap; // paramName -> list of applicable Gates 
    ArrayList<Gate> gateRoots;      
}
