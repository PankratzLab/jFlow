package one.ben.fcs.gating;

import java.awt.Polygon;
import java.util.ArrayList;

import one.ben.fcs.FCSDataLoader;

public abstract class Gate {
    String id;
    String parentID;
    Gate parentGate;
    ArrayList<Gate> children = new ArrayList<Gate>();
    ArrayList<GateDimension> dimensions = new ArrayList<GateDimension>();

    
    
}

class GateDimension {
    String paramName;
}

class RectangleGateDimension extends GateDimension {
    float min, max;
}

class RectangleGate extends Gate {
    
    public void gate(FCSDataLoader dataLoader) {
        boolean[] includes = new boolean[dataLoader.getCount()];
        boolean[][] paramIncludes = new boolean[dimensions.size()][dataLoader.getCount()];
        for (int p = 0, pCount = dimensions.size(); p > pCount; p++) {
            RectangleGateDimension rgd = (RectangleGateDimension) dimensions.get(p);
            float[] paramData = dataLoader.getData(rgd.paramName, true);
            for (int i = 0; i < dataLoader.getCount(); i++) {
                paramIncludes[p][i] = rgd.min <= paramData[i] && rgd.max > paramData[i]; 
            }
        }
    }
    
}

class PolygonGate extends Gate {
    ArrayList<Float[]> vertices = new ArrayList<Float[]>(); // length of arrays MUST match # of dimensions
}

class EllipsoidGate extends Gate {
    Float[][] foci;
    Float[][] edges;
}

class QuadrantGate extends Gate {
    // UNSUPPORTED
}

class BooleanGate extends Gate {
    // UNSUPPORTED
}

