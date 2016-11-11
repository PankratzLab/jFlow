package org.genvisis.one.ben.fcs.gating;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.gating.Gate.BooleanGate;
import org.genvisis.one.ben.fcs.gating.Gate.EllipsoidGate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.Gate.QuadrantGate;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class GateFileReader {

  // TODO remove all dependence on this variable, which may not be possibly due to info limitations
  // in files
  private static final AXIS_SCALE DEFAULT_SCALE = AXIS_SCALE.LIN;

  public static Gating readGateFile(String gateFile) throws ParserConfigurationException,
      SAXException, IOException {
    return gateFile.toLowerCase().endsWith(".wsp") || gateFile.toLowerCase().endsWith(".wspt") ? readFlowJoGatingFile(gateFile)
        : readGatingMLFile(gateFile);
  }

  public static Workbench loadWorkspace(String file) throws ParserConfigurationException, SAXException, IOException {
    Workbench wb = new Workbench();
    
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    Document doc = builder.parse(new File(file));
    doc.getDocumentElement().normalize();

    // TODO load templateGating
    NodeList groups = doc.getElementsByTagName("GroupNode");
    for (int i = 0, count = groups.getLength(); i < count; i++) {
      Element groupNode = (Element) groups.item(i);
      
    }
    
    
    
    
    NodeList samples = doc.getElementsByTagName("Sample");

    for (int i = 0, count = samples.getLength(); i < count; i++) {
      Element e = (Element) samples.item(i);
      String fcsFile = ((Element) e.getElementsByTagName("DataSet").item(0)).getAttribute("uri");
      
      Element transforms = ((Element) e.getElementsByTagName("Transformations").item(0));
      parseTransforms(transforms);
      
      Element sampleNode = (Element) e.getElementsByTagName("SampleNode").item(0);
      String id = sampleNode.getAttribute("sampleID");
      SampleNode sn = new SampleNode();
      if (fcsFile.startsWith("file:\\")) {
        fcsFile = fcsFile.substring("file:\\".length());
      } else if (fcsFile.startsWith("file:/")) {
        fcsFile = fcsFile.substring(6);
      } 
      sn.id = sampleNode.getAttribute("name");
      sn.fcsFile = fcsFile;
      sn.sampleNode = sampleNode;
      sn.doc = doc;
      Gating gs = new Gating();
      gs.setFile(file);
      NodeList nodes = sampleNode.getElementsByTagName("Population");
      gs.gateMap = GateFileUtils.buildPopGraph(nodes, true);
      gs.gateRoots = GateFileUtils.connectGates(gs.gateMap);
      gs.paramGateMap = GateFileUtils.parameterizeGates(gs.gateMap);
      for (Gate g : gs.gateMap.values()) {
          gs.allNames.add(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
      }
      sn.gating = gs;
      
      wb.samples.put(sn.id, sn);
    }
    
    return wb;
  }

  private static void parseTransforms(Element transforms) {
    NodeList nList = transforms.getChildNodes();
    for (int i = 0, count = nList.getLength(); i < count; i++) {
      Node n = nList.item(i);
      if (!n.getNodeName().startsWith("transforms:")) continue; 
      
      String param = ((Element) ((Element) n).getElementsByTagName("data-type:parameter").item(0)).getAttribute("data-type:name");
      if (n.getNodeName().endsWith("linear")) {
        String min = ((Element) n).getAttribute("transforms:minRange");
        String max = ((Element) n).getAttribute("transforms:maxRange");
        String gain = ((Element) n).getAttribute("gain");
//        System.out.println();
      } else if (n.getNodeName().endsWith("biex")) {
        String len = ((Element) n).getAttribute("transforms:length");
        String rng = ((Element) n).getAttribute("transforms:maxRange");
        String neg = ((Element) n).getAttribute("transforms:neg");
        String wid = ((Element) n).getAttribute("transforms:width");
        String pos = ((Element) n).getAttribute("transforms:pos");
//        System.out.println();
      } else if (n.getNodeName().endsWith("log")) {
        ((Element) n).getAttribute("transforms:offset");
        ((Element) n).getAttribute("transforms:decades");
//        System.out.println();
      }
    }
  }
  
  public static Gating readFlowJoGatingFile(String filename)
      throws ParserConfigurationException, SAXException, IOException {
    return readFile(filename, true);
  }

  public static Gating readGatingMLFile(String filename)
      throws ParserConfigurationException, SAXException, IOException {
    return readFile(filename, false);
  }

  private static Gating readFile(String filename, boolean flowjo)
      throws ParserConfigurationException, SAXException, IOException {
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    Document doc = builder.parse(new File(filename));
    doc.getDocumentElement().normalize();
    NodeList nodes =
        flowjo ? doc.getElementsByTagName("Population") : doc
            .getElementsByTagName("gating:Gating-ML");

    // TODO load <Transforms> and subelements from FlowJo files
    // TODO write transforms to gating-ml file
    // TODO why aren't transforms exported to Gating-ML by FlowJo?

    Gating gs = new Gating();
    gs.setFile(filename);
    gs.gateMap = flowjo ? buildPopGraph(nodes) : buildGateGraph(nodes);
    gs.gateRoots = connectGates(gs.gateMap);
    gs.paramGateMap = parameterizeGates(gs.gateMap);
    for (Gate g : gs.gateMap.values()) {
      gs.allNames.add(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
    }


    return gs;
  }

  private static HashMap<String, ArrayList<Gate>> parameterizeGates(HashMap<String, Gate> gateMap) {
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

  private static HashMap<String, Gate> buildGateGraph(NodeList topGate) {
    Element topLevel = (Element) topGate.item(0);
    NodeList allGates = topLevel.getChildNodes();
    HashMap<String, Gate> gateMap = new HashMap<String, Gate>();
    for (int i = 0, count = allGates.getLength(); i < count; i++) {
      Node gateNode = allGates.item(i);
      if (!gateNode.getNodeName().startsWith("gating:")) {
        continue;
      }
      NamedNodeMap attrs = gateNode.getAttributes();
      String id = attrs.getNamedItem("gating:id").getNodeValue();
      Node parentAttr = attrs.getNamedItem("gating:parent_id");
      String parentID = parentAttr == null ? null : parentAttr.getNodeValue();
      Gate newGate = buildGate(null, id, parentID, gateNode);
      gateMap.put(id, newGate);
    }
    return gateMap;
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
        if (popName != null && !popName.equals("")) {
          gateMap.put(popName, newGate);
        }
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
        RectangleGateDimension gd =
            new RectangleGateDimension((RectangleGate) gate, param);
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
      // // dunno what gating:distance is for yet
      // ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
      // for (int i = 0; i < dimNodes.size(); i++) {
      // Node dimNode = dimNodes.get(i);
      // String param = ((Element) getFirstChild(dimNode,
      // // "data-type:fcs-dimension")).getAttribute("data-type:name");
      // GateDimension gd = new GateDimension(param);
      // gd.paramName = param;
      // gate.dimensions.add(gd);
      // }
      // ((EllipsoidGate) gate).foci = new double[2][2];
      // ArrayList<Node> fociNodes = getChildNodes(getFirstChild(gateNode, "gating:foci"),
      // // "gating:vertex");
      // // TODO check that fociNodes.size() == 2??
      // for (int i = 0; i < ((EllipsoidGate) gate).foci.length; i++) {
      // ArrayList<Node> coordNodes = getChildNodes(fociNodes.get(i), "gating:coordinate");
      // ((EllipsoidGate) gate).foci[i][0] = Double.parseDouble(((Element)
      // // coordNodes.get(0)).getAttribute("data-type:value"));
      // ((EllipsoidGate) gate).foci[i][1] = Double.parseDouble(((Element)
      // // coordNodes.get(1)).getAttribute("data-type:value"));
      // }
      // ArrayList<Node> edgeNodes = getChildNodes(getFirstChild(gateNode, "gating:edge"),
      // // "gating:vertex");
      // ((EllipsoidGate) gate).edges = new double[4][2];
      // for (int i = 0; i < ((EllipsoidGate) gate).edges.length; i++) {
      // ArrayList<Node> coordNodes = getChildNodes(edgeNodes.get(i), "gating:coordinate");
      // ((EllipsoidGate) gate).edges[i][0] = Double.parseDouble(((Element)
      // // coordNodes.get(0)).getAttribute("data-type:value"));
      // ((EllipsoidGate) gate).edges[i][1] = Double.parseDouble(((Element)
      // // coordNodes.get(1)).getAttribute("data-type:value"));
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
