package org.genvisis.one.ben.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.GateFileUtils;
import org.genvisis.one.ben.fcs.gating.Gating;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class WSPLoader {

  private static final String[][] PANELS = {
    {"panel 1", "p1"},
    {"panel 2", "p2"},
  };
  private static final String WSP_EXT = ".wsp";
	private static final FilenameFilter WSP_FILTER = new FilenameFilter() {
		@Override
		public boolean accept(File dir, String name) {
			return name.endsWith(WSP_EXT);
		}
	};
  
  final HashMap<String, SampleNode> panel1Nodes = new HashMap<>();
  final HashMap<String, SampleNode> panel2Nodes = new HashMap<>();
  final ArrayList<SampleNode> allSamples = new ArrayList<>();
  final Logger log = new Logger();
  
  public boolean loadWorkspaces(String wspDir) {
  	String[] wspFiles = (new File(wspDir)).list(WSP_FILTER);
  	boolean allLoaded = true;
  	for (String f : wspFiles) {
  		if (!new File(wspDir + f).canRead()) {
  			System.err.println("Error - cannot access workspace file: " + wspDir + f);
  			continue;
  		}
  		try {
				loadSampleGating(wspDir + f);
			} catch (ParserConfigurationException | SAXException | IOException e) {
				log.reportException(e);
				allLoaded = false;
			}
  	}
  	return allLoaded;
  }
  
  private void loadSampleGating(String file) throws ParserConfigurationException, SAXException, IOException {
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
        Element e = (Element) n;
        String name = e.getAttributes().getNamedItem("name").getNodeValue();
        for (String chk : PANELS[0]) {
          if (name.equalsIgnoreCase(chk)) {
            panel1 = e;
            continue nodeLoop;
          }
        }
        for (String chk : PANELS[1]) {
          if (name.equalsIgnoreCase(chk)) {
            panel2 = e;
            continue nodeLoop;
          }
        }
      }
    }
    
    // TODO parse template gating
    
    ArrayList<String> panel1IDs;
    ArrayList<String> panel2IDs;
    
    panel1IDs = panel1 != null ? loadSampleList(panel1, file, 1) : new ArrayList<>();
    panel2IDs = panel2 != null ? loadSampleList(panel2, file, 2) : new ArrayList<>();
    
    NodeList samples = doc.getElementsByTagName("Sample");
    
    for (int i = 0, count = samples.getLength(); i < count; i++) {
      Element e = (Element) samples.item(i);
      String fcsFile = ((Element) e.getElementsByTagName("DataSet").item(0)).getAttribute("uri");
      
      Element transforms = (Element) e.getElementsByTagName("Transformations").item(0);
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
        log.reportError(e2.getMessage());
        sn.fcsFile = fcsFile;
      }
      sn.sampleNode = sampleNode;
      sn.doc = doc;
      Gating gs = new Gating();
      gs.setFile(file);
      NodeList nodes = sampleNode.getElementsByTagName("Population");
      gs.gateMap = GateFileUtils.buildPopGraph(nodes, false);
      gs.gateRoots = GateFileUtils.connectGates(gs.gateMap);
      gs.paramGateMap = GateFileUtils.parameterizeGates(gs.gateMap);
      for (Gate g : gs.gateMap.values()) {
          gs.getAllGateNames().add(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
      }
      sn.gating = gs;
      sn.savedTransforms = transformMap;
      
      if (panel1IDs.contains(id)) {
        panel1Nodes.put(ext.removeDirectoryInfo(sn.fcsFile), sn);
      }
      if (panel2IDs.contains(id)) {
        panel2Nodes.put(ext.removeDirectoryInfo(sn.fcsFile), sn);
      }
      allSamples.add(sn);
    }
    
  }
  
  private HashMap<String, AxisTransform> parseTransforms(Element transforms) {
  	HashMap<String, AxisTransform> map = new HashMap<>();
    NodeList nList = transforms.getChildNodes();
    for (int i = 0, count = nList.getLength(); i < count; i++) {
      Node n = nList.item(i);
      if (!n.getNodeName().startsWith("transforms:")) {
      	continue; 
      }
      
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
//        String len = ((Element) n).getAttribute("transforms:length");
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
  
  private ArrayList<String> loadSampleList(Element panelNode, String srcFile, int panel) {
  	Node node = panelNode.getElementsByTagName("SampleRefs").item(0);
  	if (node == null) {
  		System.err.println("No sample list tag for Panel " + panel + " in file: " + srcFile);
  		return new ArrayList<>();
  	}
    NodeList nList = ((Element) node).getElementsByTagName("SampleRef");
    ArrayList<String> sampleRef = new ArrayList<>();
    for (int i = 0; i < nList.getLength(); i++) {
      sampleRef.add(((Element) nList.item(i)).getAttribute("sampleID"));
    }
    return sampleRef;
  }
  
}

