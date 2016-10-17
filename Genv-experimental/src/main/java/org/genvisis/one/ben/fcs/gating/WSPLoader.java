package org.genvisis.one.ben.fcs.gating;

import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Result;
import javax.xml.transform.Source;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Numbers;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.gating.Gate.BooleanGate;
import org.genvisis.one.ben.fcs.gating.Gate.EllipsoidGate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.Gate.QuadrantGate;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class WSPLoader {

  private static final String[][] PANELS = {
    {"panel 1", "p1"},
    {"panel 2", "p2"},
  };
  
  private ArrayList<String> autoGates = new ArrayList<String>();
  
  public void load(String file) throws ParserConfigurationException, SAXException, IOException {
    long t1 = System.currentTimeMillis();
    long t2 = System.currentTimeMillis();
    
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    Document doc = builder.parse(new File(file));
    doc.getDocumentElement().normalize();
    
    NodeList nList = doc.getElementsByTagName("GroupNode");
    Element panel1 = null;
    Element panel2 = null;
    nodeLoop : for (int i = 0; i < nList.getLength(); i++) {
      Node n = nList.item(i);
      if (n.getNodeType() == Node.ELEMENT_NODE) {
        Element e = ((Element) n);
        String name = e.getAttributes().getNamedItem("name").getNodeValue();
        for (String chk : PANELS[0]) {
          if (name.toLowerCase().equals(chk)) {
            panel1 = e;
            continue nodeLoop;
          }
        }
        for (String chk : PANELS[1]) {
          if (name.toLowerCase().equals(chk)) {
            panel2 = e;
            continue nodeLoop;
          }
        }
      }
    }
    
    // TODO parse template gating
    
    ArrayList<String> panel1IDs;
    ArrayList<String> panel2IDs;
    final HashMap<String, SampleNode> panel1Nodes = new HashMap<String, SampleNode>();
    HashMap<String, SampleNode> panel2Nodes = new HashMap<String, SampleNode>();
    
    panel1IDs = panel1 != null ? loadSampleList(panel1) : new ArrayList<String>();
    panel2IDs = panel2 != null ? loadSampleList(panel2) : new ArrayList<String>();

    NodeList samples = doc.getElementsByTagName("Sample");

    String s = ext.getTimeElapsed(t1);
    System.out.println("Read: " + s);
    t1 = System.currentTimeMillis();
    
    for (int i = 0, count = samples.getLength(); i < count; i++) {
      Element e = (Element) samples.item(i);
      String fcsFile = ((Element) e.getElementsByTagName("DataSet").item(0)).getAttribute("uri");
      
      Element sampleNode = (Element) e.getElementsByTagName("SampleNode").item(0);
      String id = sampleNode.getAttribute("sampleID");
      SampleNode sn = new SampleNode();
      sn.id = id;
      if (fcsFile.startsWith("file:\\")) {
        fcsFile = fcsFile.substring("file:\\".length());
      } else if (fcsFile.startsWith("file:/")) {
        fcsFile = fcsFile.substring(6);
      } 
      sn.fcsFile = fcsFile;
      sn.sampleNode = sampleNode;
      sn.doc = doc;
      GatingStrategy gs = new GatingStrategy();
      gs.setFile(file);
      NodeList nodes = sampleNode.getElementsByTagName("Population");
      gs.gateMap = GateFileUtils.buildPopGraph(nodes);
      gs.gateRoots = GateFileUtils.connectGates(gs.gateMap);
      gs.paramGateMap = GateFileUtils.parameterizeGates(gs.gateMap);
      for (Gate g : gs.gateMap.values()) {
          gs.allNames.add(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
      }
      sn.gating = gs;
      if (panel1IDs.contains(id)) {
        panel1Nodes.put(id, sn);
      }
      if (panel2IDs.contains(id)) {
        panel2Nodes.put(id, sn);
      }
    }
    
    s = ext.getTimeElapsed(t1);
    System.out.println("Parse " + samples.getLength() + ": " + s);
    t1 = System.currentTimeMillis();
    
    int proc = Runtime.getRuntime().availableProcessors();
    final ThreadPoolExecutor threadPool = new ThreadPoolExecutor(proc, proc, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
    final ConcurrentLinkedQueue<SampleNode> queue = new ConcurrentLinkedQueue<SampleNode>();
    
    for (final String s1 : panel1IDs) {
      final SampleNode sn = panel1Nodes.get(s1);
      if (Files.exists(sn.fcsFile)) {
        queue.add(sn);
      } else {
        System.err.println("Error - file not found: " + sn.fcsFile);
      }
    }
    new Thread(new Runnable() {
      @Override
      public void run() {
        while (!threadPool.isShutdown()) {
          final SampleNode sn = queue.poll();
          if (sn == null) {
            if (threadPool.getQueue().size() == 0 && threadPool.getActiveCount() == 0) {
              threadPool.shutdown();
            }
            try {
              Thread.sleep(100);
            } catch (InterruptedException e) {}
            continue;
          }
          Runnable run = new Runnable() {
            @Override
            public void run() {
              try {
                
                updateGating(sn);
                
              } catch (IOException e) {
                System.err.println("Error - " + e.getMessage());
                // do not re-add to queue
              } catch (OutOfMemoryError e) {
                System.gc();
                queue.add(sn);
                // cleanup and re-add to queue
              }
            }
          };
          threadPool.execute(run);
        }
      }
    }).run();
    while (!threadPool.isShutdown()) {
      try {
        Thread.sleep(1000);
      } catch (InterruptedException e) {}
    }
    try {
      threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    } catch (InterruptedException e1) {}
    
    s = ext.getTimeElapsed(t1);
    System.out.println("Update: " + s);
    t1 = System.currentTimeMillis();
    
    doc.normalize();
    try {
      Transformer transformer = TransformerFactory.newInstance().newTransformer();
      transformer.setOutputProperty(OutputKeys.INDENT, "yes");
      transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
      Result output = new StreamResult(new File(ext.verifyDirFormat(ext.parseDirectoryOfFile(file)) + "outputChanged.xml"));
      Source input = new DOMSource(doc);
      transformer.transform(input, output);
    } catch (TransformerException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    s = ext.getTimeElapsed(t1);
    System.out.println("Write: " + s);
    System.out.println("Total: " + ext.getTimeElapsed(t2));
  }
  
  private void updateGating(SampleNode sn) throws IOException {
    NodeList nList = sn.sampleNode.getElementsByTagName("Population"); 
    // annoOffsetX/annoOffsetY ???
    // panelState ???
    HashMap<String, Element> map = new HashMap<String, Element>();
    HashMap<String, Element> mapG = new HashMap<String, Element>();
    for (int i = 0; i < nList.getLength(); i++) {
      Element e = (Element) nList.item(i);
      Element g = (Element) e.getElementsByTagName("Gate").item(0);
      map.put(g.getAttribute("gating:id"), e);
      mapG.put(g.getAttribute("gating:id"), g);
    }
    
    // TODO load data file
    FCSDataLoader d = new FCSDataLoader();
    d.loadData(sn.fcsFile);
    
//    HashMap<String, boolean[]> gatingMap = new HashMap<String, boolean[]>();
    
    for (String gateName : sn.gating.allNames) {
      Gate gate = sn.gating.gateMap.get(gateName);
      
//      if (gate.children.isEmpty()) {
//        continue;
//      }
//      if (autoGates.contains(gateName)) {
//        // TODO adjust gate
//      } else {
//        continue;
//      }
      
      String id = gate.getID();
      Element p = map.get(id);
      
      if (d != null) {
        boolean[] gating = gate.gate(d);
//        gatingMap.put(gate.getName() == null || gate.getName().equals("") ? gate.getID() : gate.getName(), gating);
        p.setAttribute("count", "" + Array.booleanArraySum(gating));
      }
      
      Element g = mapG.get(gate.getID());

      if (gate instanceof PolygonGate) {
        Element actGate = (Element) g.getElementsByTagName("gating:PolygonGate").item(0);
        NodeList nList2 = g.getElementsByTagName("data-type:fcs-dimension");
        int first = ((Element) nList.item(0)).getAttribute("data-type:name").equals(gate.getDimensions().get(0).paramName) ? 0 : 1;
        int second = first == 0 ? 1 : 0;
        nList2 = g.getElementsByTagName("gating:vertex");
        for (int i = nList2.getLength() - 1; i >= 0; i--) {
          actGate.removeChild(nList2.item(i));
        }
        Path2D path = ((PolygonGate) gate).getPath();
        PathIterator pi = path.getPathIterator(null);
        while (!pi.isDone()) {
          double[] coords = new double[6];
          int type = pi.currentSegment(coords);
          if (type == PathIterator.SEG_CLOSE) {
            pi.next();
            continue;
          }
          Element vert = sn.doc.createElement("gating:vertex");
          Element vert1 = sn.doc.createElement("gating:coordinate");
          Element vert2 = sn.doc.createElement("gating:coordinate");

          vert1.setAttribute("data-type:value", coords[first] + "");
          vert2.setAttribute("data-type:value", coords[second] + "");

          vert.appendChild(vert2);
          vert.appendChild(vert1);
          actGate.appendChild(vert);
          pi.next();
        }

      } else if (gate instanceof RectangleGate) {
        Element actGate = (Element) g.getElementsByTagName("gating:RectangleGate").item(0);
        NodeList nList2 = actGate.getElementsByTagName("gating:dimension");
        for (int i = nList2.getLength() - 1; i >= 0; i--) {
          actGate.removeChild(nList2.item(i));
        }
        
        for (GateDimension gd : gate.getDimensions()) {
          Element dim1 = sn.doc.createElement("gating:dimension");
          
          if (gd instanceof RectangleGateDimension) {
              RectangleGateDimension rgd = (RectangleGateDimension) gd;
              if (Numbers.isFinite(rgd.getMax())) {
                  dim1.setAttribute("gating:max", "" + Math.max(((RectangleGateDimension) gd).getMin(), ((RectangleGateDimension) gd).getMax()));
              }
              if (Numbers.isFinite(rgd.getMin())) {
                  dim1.setAttribute("gating:min", "" + Math.min(((RectangleGateDimension) gd).getMin(), ((RectangleGateDimension) gd).getMax()));
              }
          }
  
          Element dim2 = sn.doc.createElement("data-type:fcs-dimension"); 
          dim2.setAttribute("data-type:name", gd.getParam());
          dim1.appendChild(dim2);
          actGate.appendChild(dim1);
        }
        
      }
      
    }
  }
  
  private ArrayList<String> loadSampleList(Element panelNode) {
    NodeList nList = ((Element) panelNode.getElementsByTagName("SampleRefs").item(0)).getElementsByTagName("SampleRef");
    ArrayList<String> sampleRef = new ArrayList<String>();
    for (int i = 0; i < nList.getLength(); i++) {
      sampleRef.add(((Element) nList.item(i)).getAttribute("sampleID"));
    }
    return sampleRef;
  }
  
  
  public static void main(String[] args) {
    try {
      new WSPLoader().load("F:/Flow/713-715.wsp");
    } catch (ParserConfigurationException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (SAXException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
}



class GateFileUtils {

  static HashMap<String, ArrayList<Gate>> parameterizeGates(HashMap<String, Gate> gateMap) {
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

  static HashMap<String, Gate> buildPopGraph(NodeList allGates) {
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

  static ArrayList<Gate> connectGates(HashMap<String, Gate> gateMap) {
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

  static Gate buildGate(String popName, String id, String parentID, Node gateNode) {
    String gateType = gateNode.getNodeName().substring(7); // 7 = "gating:".length()
    Gate gate = null;
    if ("RectangleGate".equals(gateType)) {
      gate = new RectangleGate(null, popName, id);
      ArrayList<Node> dimNodes = getChildNodes(gateNode, "gating:dimension");
      for (int i = 0; i < dimNodes.size(); i++) {
        Node dimNode = dimNodes.get(i);
        String param =
            ((Element) getFirstChild(dimNode, "data-type:fcs-dimension"))
                .getAttribute("data-type:name");
        RectangleGateDimension gd =
            new RectangleGateDimension((RectangleGate) gate, param, AXIS_SCALE.LIN);
        String min = ((Element) dimNode).getAttribute("gating:min");
        String max = ((Element) dimNode).getAttribute("gating:max");
        // ((Element) dimNode).getAttribute("yRatio"); // TODO dunno what yRatio is used for yet
        gd.paramName = param;
        gd.setMin("".equals(min) ? Float.NEGATIVE_INFINITY : Float.parseFloat(min));
        gd.setMax("".equals(max) ? Float.POSITIVE_INFINITY : Float.parseFloat(max));
        gate.dimensions.add(gd);
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
        GateDimension gd = new GateDimension(gate, param, AXIS_SCALE.LIN);
        gd.paramName = param;
        gate.dimensions.add(gd);
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
      ((PolygonGate) gate).prepGating();
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
