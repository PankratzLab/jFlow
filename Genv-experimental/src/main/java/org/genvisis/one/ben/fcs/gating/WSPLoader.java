package org.genvisis.one.ben.fcs.gating;

import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ConcurrentLinkedQueue;
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
import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.gating.Gate.BooleanGate;
import org.genvisis.one.ben.fcs.gating.Gate.EllipsoidGate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.Gate.QuadrantGate;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.genvisis.one.ben.fcs.sub.EMInitializer;
import org.genvisis.one.ben.fcs.sub.EMModel;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import edu.stanford.facs.logicle.Logicle;

public class WSPLoader {

  private static final String[][] PANELS = {
    {"panel 1", "p1"},
    {"panel 2", "p2"},
  };
  
  private void updateGateParams(FCSDataLoader dataLoader, ArrayList<Gate> gates) {
  	for (Gate g : gates) {
  		String p = g.getXDimension().paramName;
  		String d = dataLoader.getInternalParamName(p);
  		if (!p.equals(d)) {
  			g.getXDimension().paramName = d;
  		}
  		p = g.getYDimension().paramName;
  		d = dataLoader.getInternalParamName(p);
  		if (!p.equals(d)) {
  			g.getYDimension().paramName = d;
  		}
  		updateGateParams(dataLoader, g.children);
  	}
  }
  
  public <T extends SampleProcessor> void load(String file, Class<T> processorClass) throws ParserConfigurationException, SAXException, IOException {
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

    ArrayList<SampleNode> allSamples = new ArrayList<>();
    
    NodeList samples = doc.getElementsByTagName("Sample");

    String s = ext.getTimeElapsed(t1);
    System.out.println("Read: " + s);
    t1 = System.currentTimeMillis();
    
    for (int i = 0, count = samples.getLength(); i < count; i++) {
      Element e = (Element) samples.item(i);
      String fcsFile = ((Element) e.getElementsByTagName("DataSet").item(0)).getAttribute("uri");
      
      Element transforms = ((Element) e.getElementsByTagName("Transformations").item(0));
      HashMap<String, AxisTransform> transformMap = parseTransforms(transforms);
      
      Element sampleNode = (Element) e.getElementsByTagName("SampleNode").item(0);
      String id = sampleNode.getAttribute("sampleID");
      SampleNode sn = new SampleNode();
      sn.id = id;
      if (fcsFile.startsWith("file:\\")) {
        fcsFile = fcsFile.substring("file:\\".length());
      } else if (fcsFile.startsWith("file:/")) {
        fcsFile = fcsFile.substring(6);
      } 
      try {
        sn.fcsFile = URLDecoder.decode(fcsFile, "utf-8");
      } catch (UnsupportedEncodingException e2) {
        System.err.println("Error - " + e2.getMessage());
        sn.fcsFile = fcsFile;
      }
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
      sn.savedTransforms = transformMap;
      
      if (panel1IDs.contains(id)) {
        panel1Nodes.put(id, sn);
      }
      if (panel2IDs.contains(id)) {
        panel2Nodes.put(id, sn);
      }
      allSamples.add(sn);
    }
    
    s = ext.getTimeElapsed(t1);
    System.out.println("Parse " + samples.getLength() + ": " + s);
    t1 = System.currentTimeMillis();
    
    int proc = Runtime.getRuntime().availableProcessors();
    final ThreadPoolExecutor threadPool = new ThreadPoolExecutor(proc, proc, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
    final ConcurrentLinkedQueue<SampleNode> queue = new ConcurrentLinkedQueue<SampleNode>();
    
    for (SampleNode sn : allSamples) {
    	if (Files.exists(sn.fcsFile)) {
    		queue.add(sn);
    	}
    }
//    for (final String s1 : panel1IDs) {
//      final SampleNode sn = panel1Nodes.get(s1);
//      if (Files.exists(sn.fcsFile)) {
//        queue.add(sn);
//      } else {
//        System.err.println("Error - file not found: " + sn.fcsFile);
//      }
//    }
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
                  processorClass.getConstructor(WSPLoader.class).newInstance(WSPLoader.this).processSample(sn);
                  
              } catch (IOException e) {
                System.err.println("Error - " + e.getMessage());
                // do not re-add to queue
              } catch (OutOfMemoryError e) {
                System.gc();
                queue.add(sn);
                // cleanup and re-add to queue
              } catch (InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException
                  | NoSuchMethodException | SecurityException  e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
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
  
  private HashMap<String, AxisTransform> parseTransforms(Element transforms) {
  	HashMap<String, AxisTransform> map = new HashMap<>();
    NodeList nList = transforms.getChildNodes();
    for (int i = 0, count = nList.getLength(); i < count; i++) {
      Node n = nList.item(i);
      if (!n.getNodeName().startsWith("transforms:")) continue; 
      
      String param = ((Element) ((Element) n).getElementsByTagName("data-type:parameter").item(0)).getAttribute("data-type:name");
      AxisTransform at = null;
      if (n.getNodeName().endsWith("linear")) {
        String min = ((Element) n).getAttribute("transforms:minRange");
        String max = ((Element) n).getAttribute("transforms:maxRange");
        String gain = ((Element) n).getAttribute("gain");
        double m1, m2, g;
        m1 = Double.parseDouble(min);
        m2 = Double.parseDouble(max);
        g = Double.parseDouble(gain);
        at = AxisTransform.createLinearTransform(m1, m2, g);
      } else if (n.getNodeName().endsWith("biex")) {
        String len = ((Element) n).getAttribute("transforms:length");
        String rng = ((Element) n).getAttribute("transforms:maxRange");
        String neg = ((Element) n).getAttribute("transforms:neg");
        String wid = ((Element) n).getAttribute("transforms:width");
        String pos = ((Element) n).getAttribute("transforms:pos");
        double T = Double.parseDouble(rng);
        double W = Math.log10(Math.abs(Double.parseDouble(wid)));
        double M = Double.parseDouble(pos);
        double A = Double.parseDouble(neg);
        at = AxisTransform.createBiexTransform(T, W, M, A);
      } else if (n.getNodeName().endsWith("log")) {
        ((Element) n).getAttribute("transforms:offset");
        ((Element) n).getAttribute("transforms:decades");
        // TODO log scale
        at = AxisTransform.createBiexTransform();
      }
      map.put(param, at);
    }
    return map;
  }

  abstract class SampleProcessor {
    public abstract void processSample(SampleNode sn) throws IOException;
  }
  
  abstract class AbstractSampleProcessor extends SampleProcessor {
    NodeList popList;
    HashMap<String, Element> popMap = new HashMap<String, Element>();
    HashMap<String, Element> gateMap = new HashMap<String, Element>();
    FCSDataLoader d;
    
    private AbstractSampleProcessor() {}
    
    void loadPopsAndGates(SampleNode sn) {
      popList = sn.sampleNode.getElementsByTagName("Population"); 
      // annoOffsetX/annoOffsetY ???
      // panelState ???
      for (int i = 0; i < popList.getLength(); i++) {
        Element e = (Element) popList.item(i);
        Element g = (Element) e.getElementsByTagName("Gate").item(0);
        popMap.put(g.getAttribute("gating:id"), e);
        gateMap.put(g.getAttribute("gating:id"), g);
      }
    }
    
    void loadData(SampleNode sn) throws IOException {
      long t1 = System.currentTimeMillis();
      d = new FCSDataLoader();
      d.loadData(sn.fcsFile);
      while(d.getLoadState() != LOAD_STATE.LOADED) {
        try {
          Thread.sleep(100);
        } catch (InterruptedException e) {
        }
      }
      (new Logger()).reportTimeElapsed("Loaded FCS ... ", t1);
      d.setTransformMap(sn.savedTransforms);
      updateGateParams(d, sn.gating.gateRoots);
      sn.gating.paramGateMap = GateFileUtils.parameterizeGates(sn.gating.gateMap);
    }
    
  }
  
  class DataDumper extends AbstractSampleProcessor {
    public DataDumper() { /* TODO Auto-generated constructor stub */ }
    @Override
    public void processSample(SampleNode sn) throws IOException {
      loadData(sn);
      for (String param : d.getAllDisplayableNames(DATA_SET.ALL)) {
        double[] data = d.getData(param, true);
        double[] mm = Array.minMax(data);
        System.out.println(param + " min: " + mm[0] + "; max: " + mm[1]);
      }
    }
  }
  
  class LeafDataExporter extends AbstractSampleProcessor {
		public LeafDataExporter() {
			// TODO Auto-generated constructor stub
		}
    private static final int SAMPLE_SIZE = 1000;
    
    public void processSample(SampleNode sn) throws IOException {
      if (!Files.exists(sn.fcsFile)) return;
      loadPopsAndGates(sn);
      loadData(sn);
      
      ArrayList<String> params = new ArrayList<>();
      for (String s : EMInitializer.DATA_COLUMNS) {
        params.add(s);
      }
      HashSet<Gate> leafGates = sn.gating.getAllLeafGates();
//      leafGates.add(sn.gating.getRootGates().get(0));
      HashSet<Gate> map = new HashSet<>();
      for (Gate g : leafGates) {
        map.add(g);
      }
      String outputFileRoot = ext.verifyDirFormat(ext.parseDirectoryOfFile(sn.fcsFile)) + "sampling/" + ext.rootOf(sn.fcsFile, true) + "_";
      System.out.println("Exporting data for " + map.size() + " gates : " + params.toString());
      for (Gate g : map) {
        String outputFile = outputFileRoot + ext.replaceWithLinuxSafeCharacters(g.getID() + "_" + g.getXDimension().paramName + "_" + g.getYDimension().paramName + ".xln", false);
        System.out.println("... to file " + outputFile);
        PrintWriter writer = Files.getAppropriateWriter(outputFile);
        for (int i = 0; i < params.size(); i++) {
          writer.print(params.get(i));
          if (i < params.size() - 1) {
            writer.print("\t");
          }
        }
        writer.println();
        boolean[] incl = g.gate(d);
        int[] indices = Array.booleanArrayToIndices(incl);
        Random rand = new Random();
        for (int i = 0; i < SAMPLE_SIZE; i++) {
          int ind = indices[rand.nextInt(indices.length)];
          double[] line = d.getDataLine(params, ind);
          for (int l = 0; l < line.length; l++) {
            writer.print(line[l]);
            if (l < line.length - 1) {
              writer.print("\t");
            }
          }
          writer.println();
        }
        writer.flush();
        writer.close();
      }
      
      d.emptyAndReset();
      System.gc();
      
    }
    
  }
  
  class BasicSampleProcessor extends AbstractSampleProcessor {
    
    public void processSample(SampleNode sn) throws IOException {
      loadPopsAndGates(sn);
      loadData(sn);
      
      
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
        Element p = popMap.get(id);
        
        if (d != null) {
          boolean[] gating = gate.gate(d);
  //        gatingMap.put(gate.getName() == null || gate.getName().equals("") ? gate.getID() : gate.getName(), gating);
          p.setAttribute("count", "" + Array.booleanArraySum(gating));
        }
        
        Element g = gateMap.get(gate.getID());
  
        if (gate instanceof PolygonGate) {
          Element actGate = (Element) g.getElementsByTagName("gating:PolygonGate").item(0);
          NodeList nList2 = g.getElementsByTagName("data-type:fcs-dimension");
          int first = ((Element) popList.item(0)).getAttribute("data-type:name").equals(gate.getXDimension().paramName) ? 0 : 1;
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
          
          GateDimension gd = gate.getXDimension();
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
          if ((gd = gate.getYDimension()) != null) {
            dim1 = sn.doc.createElement("gating:dimension");
            
            if (gd instanceof RectangleGateDimension) {
              RectangleGateDimension rgd = (RectangleGateDimension) gd;
              if (Numbers.isFinite(rgd.getMax())) {
                dim1.setAttribute("gating:max", "" + Math.max(((RectangleGateDimension) gd).getMin(), ((RectangleGateDimension) gd).getMax()));
              }
              if (Numbers.isFinite(rgd.getMin())) {
                dim1.setAttribute("gating:min", "" + Math.min(((RectangleGateDimension) gd).getMin(), ((RectangleGateDimension) gd).getMax()));
              }
            }
            
            dim2 = sn.doc.createElement("data-type:fcs-dimension"); 
            dim2.setAttribute("data-type:name", gd.getParam());
            dim1.appendChild(dim2);
            actGate.appendChild(dim1);
          }
          
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
//      new WSPLoader().load("F:/Flow/713-715.wsp", AbstractSampleProcessor.class);
      new WSPLoader().load("F:/Flow/controlFCS/07-Nov-2016.wsp", LeafDataExporter.class);
//      new WSPLoader().load("F:/Flow/controlFCS/10-Aug-2016 ctl B.wsp", DataDumper.class);
      
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

  static HashMap<String, Gate> buildPopGraph(NodeList allGates, boolean flowJo) {
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

  static ArrayList<Gate> connectGates(HashMap<String, Gate> gateMap) {
    HashSet<Gate> rootGates = new HashSet<Gate>();
    for (Gate g : gateMap.values()) {
      if (null != g.parentID && !"".equals(g.parentID)) {
        g.parentGate = gateMap.get(g.parentID);
        g.parentGate.children.add(g);
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
    LIN_PARAMS.add("time");
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

