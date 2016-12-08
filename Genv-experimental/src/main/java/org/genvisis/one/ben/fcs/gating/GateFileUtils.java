package org.genvisis.one.ben.fcs.gating;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.gating.Gate.BooleanGate;
import org.genvisis.one.ben.fcs.gating.Gate.EllipsoidGate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.Gate.QuadrantGate;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class GateFileUtils {

	public static void updateGateParams(FCSDataLoader dataLoader, List<Gate> gates) {
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
	
	public static HashMap<String, ArrayList<Gate>> parameterizeGates(HashMap<String, Gate> gateMap) {
    HashMap<String, ArrayList<Gate>> paramGates = new HashMap<String, ArrayList<Gate>>();
    for (Gate g : gateMap.values()) {
      ArrayList<Gate> gates = paramGates.get(g.getXDimension().paramName);
      if (gates == null) {
        gates = new ArrayList<Gate>();
        paramGates.put(g.getXDimension().paramName, gates);
      }
      gates.add(g);
      if (g.getYDimension() != null) {
        gates = paramGates.get(g.getYDimension().paramName);
        if (gates == null) {
          gates = new ArrayList<Gate>();
          paramGates.put(g.getYDimension().paramName, gates);
        }
        gates.add(g);
      }
    }
    return paramGates;
  }

  public static HashMap<String, Gate> buildPopGraph(NodeList allGates, boolean flowJo) {
    HashMap<String, Gate> gateMap = new HashMap<String, Gate>();
    for (int i = 0, count = allGates.getLength(); i < count; i++) {
      Element popNode = (Element) allGates.item(i);
      String popName = popNode.getAttribute("name");

      Element gateNode = null;
      for (int n = 0, cnt = popNode.getChildNodes().getLength(); n < cnt; n++) {
        if (popNode.getChildNodes().item(n).getNodeName().equals("Gate")) {
          gateNode = (Element) popNode.getChildNodes().item(n);
          break;
        }
      }

      String id = gateNode.getAttribute("gating:id");
      String parentID = gateNode.getAttribute("gating:parent_id");
      Node actualGateNode = null;
      for (int n = 0, cnt = gateNode.getChildNodes().getLength(); n < cnt; n++) {
        if (gateNode.getChildNodes().item(n).getNodeName().startsWith("gating:")) {
          actualGateNode = gateNode.getChildNodes().item(n);
          break;
        }
      }
      if (actualGateNode != null) {
        Gate newGate = buildGate(popName, id, parentID, actualGateNode, flowJo);
        gateMap.put(id, newGate);
        if (popName != null && !popName.equals("")) {
          gateMap.put(popName, newGate);
        }
      }

    }
    return gateMap;
  }

  public static ArrayList<Gate> connectGates(HashMap<String, Gate> gateMap) {
    HashSet<Gate> rootGates = new HashSet<Gate>();
    for (Gate g : gateMap.values()) {
      if (null != g.parentID && !"".equals(g.parentID)) {
        g.parentGate = gateMap.get(g.parentID);
        if (!g.parentGate.children.contains(g)) {
        	g.parentGate.children.add(g);
        }
      } else {
        rootGates.add(g);
      }
    }
    return new ArrayList<Gate>(rootGates);
  }
  
  private static HashSet<String> LIN_PARAMS = new HashSet<>();
  {
    LIN_PARAMS.add("time");
    LIN_PARAMS.add("fsc-a");
    LIN_PARAMS.add("fsc-w");
    LIN_PARAMS.add("fsc-h");
    LIN_PARAMS.add("ssc-a");
    LIN_PARAMS.add("ssc-w");
    LIN_PARAMS.add("ssc-h");
  }
  private static AXIS_SCALE getAxisScaleHardcoded(String param) {
    if (LIN_PARAMS.contains(param.toLowerCase())) {
      return AXIS_SCALE.LIN;
    } else {
      return AXIS_SCALE.BIEX;
    }
  }
  
  static Gate buildGate(String popName, String id, String parentID, Node gateNode, boolean flowJo) {
    String gateType = gateNode.getNodeName().substring(7); // 7 = "gating:".length()
    Gate gate = null;
    if ("RectangleGate".equals(gateType)) {
      gate = new RectangleGate(null, popName, id);
      ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
      for (int i = 0; i < dimNodes.size(); i++) {
        Node dimNode = dimNodes.get(i);
        String param = ((Element) getFirstChild(dimNode, "data-type:fcs-dimension")).getAttribute("data-type:name");
        RectangleGateDimension gd = new RectangleGateDimension((RectangleGate) gate, param);
        String min = ((Element) dimNode).getAttribute("gating:min");
        String max = ((Element) dimNode).getAttribute("gating:max");
        // ((Element) dimNode).getAttribute("yRatio"); // TODO dunno what yRatio is used for yet
        gd.paramName = param;
        gd.setMin("".equals(min) ? Float.NEGATIVE_INFINITY : Float.parseFloat(min));
        gd.setMax("".equals(max) ? Float.POSITIVE_INFINITY : Float.parseFloat(max));
        if (i == 0) {
          gate.setXDimension(gd);
        } else if (i == 1) {
          gate.setYDimension(gd);
        }
      }
    } else if ("EllipsoidGate".equals(gateType)) {
      gate = new EllipsoidGate();
      // // String gatingDistance = ((Element) gateNode).getAttribute("gating:distance"); // TODO
      // dunno what gating:distance is for yet
      // ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
      // for (int i = 0; i < dimNodes.size(); i++) {
      // Node dimNode = dimNodes.get(i);
      // String param = ((Element) getFirstChild(dimNode,
      // "data-type:fcs-dimension")).getAttribute("data-type:name");
      // GateDimension gd = new GateDimension(param);
      // gd.paramName = param;
      // gate.dimensions.add(gd);
      // }
      // ((EllipsoidGate) gate).foci = new double[2][2];
      // ArrayList<Node> fociNodes = getChildNodes(getFirstChild(gateNode, "gating:foci"),
      // "gating:vertex");
      // // TODO check that fociNodes.size() == 2??
      // for (int i = 0; i < ((EllipsoidGate) gate).foci.length; i++) {
      // ArrayList<Node> coordNodes = getChildNodes(fociNodes.get(i), "gating:coordinate");
      // ((EllipsoidGate) gate).foci[i][0] = Double.parseDouble(((Element)
      // coordNodes.get(0)).getAttribute("data-type:value"));
      // ((EllipsoidGate) gate).foci[i][1] = Double.parseDouble(((Element)
      // coordNodes.get(1)).getAttribute("data-type:value"));
      // }
      // ArrayList<Node> edgeNodes = getChildNodes(getFirstChild(gateNode, "gating:edge"),
      // "gating:vertex");
      // ((EllipsoidGate) gate).edges = new double[4][2];
      // for (int i = 0; i < ((EllipsoidGate) gate).edges.length; i++) {
      // ArrayList<Node> coordNodes = getChildNodes(edgeNodes.get(i), "gating:coordinate");
      // ((EllipsoidGate) gate).edges[i][0] = Double.parseDouble(((Element)
      // coordNodes.get(0)).getAttribute("data-type:value"));
      // ((EllipsoidGate) gate).edges[i][1] = Double.parseDouble(((Element)
      // coordNodes.get(1)).getAttribute("data-type:value"));
      // }
    } else if ("PolygonGate".equals(gateType)) {
      gate = new PolygonGate(null, popName, id);
      String resStr = gateNode.getAttributes().getNamedItem("gateResolution").getNodeValue();
      int res = -1;
      try {
        res = Integer.parseInt(resStr);
      } catch (NumberFormatException e) {
      }
      ((PolygonGate) gate).setGateResolution(res);
      ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
      for (int i = 0; i < dimNodes.size(); i++) {
        Node dimNode = dimNodes.get(i);
        String param =
            ((Element) getFirstChild(dimNode, "data-type:fcs-dimension"))
                .getAttribute("data-type:name");
        GateDimension gd = new GateDimension(gate, param);
        gd.paramName = param;
        if (i == 0) {
          gate.setXDimension(gd);
        } else if (i == 1) {
          gate.setYDimension(gd);
        }
      }
      ArrayList<Node> vertexNodes = getChildNodes(gateNode, "gating:vertex");
      for (Node n : vertexNodes) {
        ArrayList<Node> coordNodes = getChildNodes(n, "gating:coordinate");
        Double fX =
            Double.parseDouble(((Element) coordNodes.get(0)).getAttribute("data-type:value"));
        Double fY =
            Double.parseDouble(((Element) coordNodes.get(1)).getAttribute("data-type:value"));
        ((PolygonGate) gate).addVertex(fX, fY);
      }
//      ((PolygonGate) gate).prepGating();
      ((PolygonGate) gate).setShouldMimicFlowJoGating(flowJo);
    } else if ("QuadrantGate".equals(gateType)) {
      gate = new QuadrantGate();

    } else if ("BooleanGate".equals(gateType)) {
      gate = new BooleanGate();

    }
    gate.id = id;
    gate.parentID = parentID;
    return gate;
  }

  private static Node getFirstChild(Node nd, String name) {
    NodeList children = nd.getChildNodes();
    for (int i = 0; i < children.getLength(); i++) {
      if (children.item(i).getNodeName().equals(name)) {
        return children.item(i);
      }
    }
    return null;
  }

  private static ArrayList<Node> getChildNodes(Node nd, String name) {
    ArrayList<Node> retNodes = new ArrayList<Node>();
    NodeList children = nd.getChildNodes();
    for (int i = 0; i < children.getLength(); i++) {
      if (children.item(i).getNodeName().equals(name)) {
        retNodes.add(children.item(i));
      }
    }
    return retNodes;
  }
  
}