package org.genvisis.one.ben.fcs.gating;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.one.ben.fcs.gating.Gate.BooleanGate;
import org.genvisis.one.ben.fcs.gating.Gate.EllipsoidGate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.Gate.QuadrantGate;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class GateFileReader {

    public static GatingStrategy readGateFile(String filename) throws ParserConfigurationException, SAXException, IOException {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse(new File(filename));
        doc.getDocumentElement().normalize();
        
        // TODO read Population nodes instead (one gate per population, pop. defines name, gives count, etc
        NodeList allPops = doc.getElementsByTagName("Population");
        
        GatingStrategy gs = new GatingStrategy();
        gs.setFile(filename);
        gs.gateMap = buildPopGraph(allPops);
        gs.gateRoots = connectGates(gs.gateMap);
        gs.paramGateMap = parameterizeGates(gs.gateMap);
        
        return gs;
    }
    
    private static HashMap<String, ArrayList<Gate>> parameterizeGates(HashMap<String, Gate> gateMap) {
        HashMap<String, ArrayList<Gate>> paramGates = new HashMap<String, ArrayList<Gate>>();
        for (Gate g : gateMap.values()) {
            for (GateDimension gd : g.dimensions) {
                ArrayList<Gate> gates = paramGates.get(gd.paramName);
                if (gates == null) {
                    gates = new ArrayList<Gate>();
                    paramGates.put(gd.paramName, gates);
                }
                gates.add(g);
            }
        }
        return paramGates;
    }
    
    private static HashMap<String, Gate> buildPopGraph(NodeList allGates) {
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
                Gate newGate = buildGate(popName, id, parentID, actualGateNode); 
                gateMap.put(id, newGate);
            }
            
        }
        return gateMap;
    }
    
    private static ArrayList<Gate> connectGates(HashMap<String, Gate> gateMap) {
        ArrayList<Gate> rootGates = new ArrayList<Gate>();
        for (Gate g : gateMap.values()) {
            if (null != g.parentID && !"".equals(g.parentID)) {
                g.parentGate = gateMap.get(g.parentID);
                g.parentGate.children.add(g);
            } else {
                rootGates.add(g);
            }
        }
        return rootGates;
    }
    
    private static Gate buildGate(String popName, String id, String parentID, Node gateNode) {
        String gateType = gateNode.getNodeName().substring(7); // 7 = "gating:".length()
        Gate gate = null;
        if ("RectangleGate".equals(gateType)) {
            gate = new RectangleGate(null, popName, id);
            ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
            for (int i = 0; i < dimNodes.size(); i++) {
                Node dimNode = dimNodes.get(i);
                String param = ((Element) getFirstChild(dimNode, "data-type:fcs-dimension")).getAttribute("data-type:name");
                RectangleGateDimension gd = new RectangleGateDimension(param);
                String min = ((Element) dimNode).getAttribute("gating:min");
                String max = ((Element) dimNode).getAttribute("gating:max");
//                ((Element) dimNode).getAttribute("yRatio"); // TODO dunno what yRatio is used for yet
                gd.paramName = param;
                gd.min = "".equals(min) ? Float.NEGATIVE_INFINITY : Float.parseFloat(min);
                gd.max = "".equals(max) ? Float.POSITIVE_INFINITY : Float.parseFloat(max);
                gate.dimensions.add(gd);
            }
        } else if ("EllipsoidGate".equals(gateType)) {
            gate = new EllipsoidGate();
////            String gatingDistance = ((Element) gateNode).getAttribute("gating:distance"); // TODO dunno what gating:distance is for yet
//            ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
//            for (int i = 0; i < dimNodes.size(); i++) {
//                Node dimNode = dimNodes.get(i);
//                String param = ((Element) getFirstChild(dimNode, "data-type:fcs-dimension")).getAttribute("data-type:name");
//                GateDimension gd = new GateDimension(param);
//                gd.paramName = param;
//                gate.dimensions.add(gd);
//            }
//            ((EllipsoidGate) gate).foci = new double[2][2];  
//            ArrayList<Node> fociNodes = getChildNodes(getFirstChild(gateNode, "gating:foci"), "gating:vertex");
//            // TODO check that fociNodes.size() == 2??
//            for (int i = 0; i < ((EllipsoidGate) gate).foci.length; i++) {
//                ArrayList<Node> coordNodes = getChildNodes(fociNodes.get(i), "gating:coordinate");
//                ((EllipsoidGate) gate).foci[i][0] = Double.parseDouble(((Element) coordNodes.get(0)).getAttribute("data-type:value"));
//                ((EllipsoidGate) gate).foci[i][1] = Double.parseDouble(((Element) coordNodes.get(1)).getAttribute("data-type:value"));
//            }
//            ArrayList<Node> edgeNodes = getChildNodes(getFirstChild(gateNode, "gating:edge"), "gating:vertex");
//            ((EllipsoidGate) gate).edges = new double[4][2];
//            for (int i = 0; i < ((EllipsoidGate) gate).edges.length; i++) {
//                ArrayList<Node> coordNodes = getChildNodes(edgeNodes.get(i), "gating:coordinate");
//                ((EllipsoidGate) gate).edges[i][0] = Double.parseDouble(((Element) coordNodes.get(0)).getAttribute("data-type:value"));
//                ((EllipsoidGate) gate).edges[i][1] = Double.parseDouble(((Element) coordNodes.get(1)).getAttribute("data-type:value"));
//            }
        } else if ("PolygonGate".equals(gateType)) {
            gate = new PolygonGate(null, popName, id);
            String resStr = gateNode.getAttributes().getNamedItem("gateResolution").getNodeValue();
            int res = -1;
            try { res = Integer.parseInt(resStr); } catch (NumberFormatException e) {}
            ((PolygonGate)gate).setGateResolution(res);
            ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
            for (int i = 0; i < dimNodes.size(); i++) {
                Node dimNode = dimNodes.get(i);
                String param = ((Element) getFirstChild(dimNode, "data-type:fcs-dimension")).getAttribute("data-type:name");
                GateDimension gd = new GateDimension(param);
                gd.paramName = param;
                gate.dimensions.add(gd);
            }
            ArrayList<Node> vertexNodes = getChildNodes(gateNode, "gating:vertex");
            for (Node n : vertexNodes) {
                ArrayList<Node> coordNodes = getChildNodes(n, "gating:coordinate");
                Double fX = Double.parseDouble(((Element) coordNodes.get(0)).getAttribute("data-type:value"));
                Double fY = Double.parseDouble(((Element) coordNodes.get(1)).getAttribute("data-type:value"));
                ((PolygonGate) gate).addVertex(fX, fY);
            }
            ((PolygonGate)gate).prepGating();
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