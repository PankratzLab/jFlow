package one.ben.fcs.gating;

import java.util.ArrayList;

public class Gate {
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

