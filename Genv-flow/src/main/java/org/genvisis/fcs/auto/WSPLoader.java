package org.genvisis.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.fcs.AbstractPanel2.AxisTransform;
import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.GateFileUtils;
import org.genvisis.fcs.gating.Gating;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class WSPLoader {

  private static final String WSP_EXT = ".wsp";
  private static final FilenameFilter WSP_FILTER = new FilenameFilter() {

    @Override
    public boolean accept(File dir, String name) {
      return name.endsWith(WSP_EXT);
    }
  };

  private final Set<Panel> panels = new HashSet<>();

  private final Map<Panel, Map<String, SampleNode>> panelNodeMap = new HashMap<>();

  final ArrayList<SampleNode> allSamples = new ArrayList<>();
  final Logger log = new Logger();

  public WSPLoader(Panel... panels) {
    for (Panel p : panels) {
      this.panels.add(p);
      panelNodeMap.put(p, new HashMap<>());
    }
  }

  public int loadWorkspaces(String wspD) {
    String wspDir = ext.verifyDirFormat(wspD);
    File dir = new File(wspDir);
    if (!dir.canRead()) {
      log.reportError("Cannot read workspace files in directory " + wspDir);
      return 0;
    }
    String[] wspFiles = dir.list(WSP_FILTER);
    log.reportTime("Processing " + wspFiles.length + " wsp files");
    int numLoaded = 0;
    for (String f : wspFiles) {
      File sub = new File(wspDir + f);
      if (!sub.canRead()) {
        System.err.println("Error - cannot access workspace file: " + wspDir + f);
        continue;
      }
      if (!sub.isDirectory()) {
        try {
          loadSampleGating(wspDir + f);
          numLoaded++;
        } catch (ParserConfigurationException | SAXException | IOException e) {
          log.reportException(e);
          numLoaded = -1;
        }
      }
    }
    for (File f : dir.listFiles()) {
      if (f.isDirectory()) {
        int subLoaded = loadWorkspaces(f.getAbsolutePath());
        if (subLoaded == -1) {
          numLoaded = -1;
        } else {
          numLoaded += subLoaded;
        }
      }
    }
    return numLoaded;
  }

  private void loadSampleGating(String file) throws ParserConfigurationException, SAXException,
                                             IOException {
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    Document doc = builder.parse(new File(file));
    doc.getDocumentElement().normalize();

    NodeList nList = doc.getElementsByTagName("GroupNode");
    Map<Panel, Element> panelElements = new HashMap<>();

    nodeLoop: for (int i = 0; i < nList.getLength(); i++) {
      Node n = nList.item(i);
      if (n.getNodeType() == Node.ELEMENT_NODE) {
        Element e = (Element) n;
        String name = e.getAttributes().getNamedItem("name").getNodeValue();
        for (Panel p : panels) {
          for (String chk : p.aliases) {
            if (name.equalsIgnoreCase(chk)) {
              panelElements.put(p, e);
              continue nodeLoop;
            }
          }
        }
      }
    }

    // TODO parse template gating

    Map<Panel, List<String>> panelIDs = new HashMap<>();
    for (Panel p : panels) {
      panelIDs.put(p, panelElements.containsKey(p) ? loadSampleList(panelElements.get(p), file, p)
                                                   : new ArrayList<>());
    }

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
      // try {
      // sn.fcsFile = URLDecoder.decode(fcsFile, "utf-8");
      // } catch (UnsupportedEncodingException e2) {
      // log.reportError(e2.getMessage());
      sn.fcsFile = fcsFile;
      if (sn.fcsFile.contains("%20")) {
        sn.fcsFile = sn.fcsFile.replaceAll("%20", " ");
      }
      // }
      sn.sampleNode = sampleNode;
      sn.doc = doc;
      Gating gs = new Gating();
      gs.setFile(file);
      NodeList nodes = sampleNode.getElementsByTagName("Population");
      gs.gateMap = GateFileUtils.buildPopGraph(nodes, false);
      gs.gateRoots = GateFileUtils.connectGates(gs.gateMap);
      gs.paramGateMap = GateFileUtils.parameterizeGates(gs.gateMap);
      for (Gate g : gs.gateMap.values()) {
        gs.getAllGateNames()
          .add(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());
      }
      sn.gating = gs;
      sn.savedTransforms = transformMap;
      sn.wspFile = file;

      boolean foundID = false;
      for (Panel p : panels) {
        if (panelIDs.get(p).contains(id)) {
          panelNodeMap.get(p).put(ext.removeDirectoryInfo(sn.fcsFile), sn);
          foundID = true;
        }
      }
      if (!foundID) {
        for (Panel p : panels) {
          for (String a : p.aliases) {
            if (sn.fcsFile.toLowerCase().contains(a)) {
              panelNodeMap.get(p).put(ext.removeDirectoryInfo(sn.fcsFile), sn);
            }
          }
        }
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

      String param = ((Element) ((Element) n).getElementsByTagName("data-type:parameter")
                                             .item(0)).getAttribute("data-type:name");
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
        // String len = ((Element) n).getAttribute("transforms:length");
        String rng = ((Element) n).getAttribute("transforms:maxRange");
        String neg = ((Element) n).getAttribute("transforms:neg");
        String wid = ((Element) n).getAttribute("transforms:width");
        String pos = ((Element) n).getAttribute("transforms:pos");
        double T = Double.parseDouble(rng);
        double W = Math.log10(Math.abs(Double.parseDouble(wid)));
        double M = Double.parseDouble(pos);
        double A = Double.parseDouble(neg);
        if (2 * W > M) {
          M = Math.ceil(2 * W);
        }
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

  private ArrayList<String> loadSampleList(Element panelNode, String srcFile, Panel panel) {
    Node node = panelNode.getElementsByTagName("SampleRefs").item(0);
    if (node == null) {
      System.err.println("No sample list tag for Panel " + panel.getName() + " in file: "
                         + srcFile);
      return new ArrayList<>();
    }
    NodeList nList = ((Element) node).getElementsByTagName("SampleRef");
    ArrayList<String> sampleRef = new ArrayList<>();
    for (int i = 0; i < nList.getLength(); i++) {
      sampleRef.add(((Element) nList.item(i)).getAttribute("sampleID"));
    }
    return sampleRef;
  }

  public Map<String, SampleNode> getPanelNodes(Panel p) {
    return panelNodeMap.containsKey(p) ? panelNodeMap.get(p) : new HashMap<>();
  }

  public static List<Panel> loadPanelsFromFile(InputStream panelDefFile) {
    List<Panel> panelsFound = new ArrayList<>();
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    try {
      DocumentBuilder builder = factory.newDocumentBuilder();
      Document doc = builder.parse(panelDefFile);
      doc.getDocumentElement().normalize();

      NodeList panels = doc.getElementsByTagName("panel");
      for (int i = 0, count = panels.getLength(); i < count; i++) {
        Element panelNode = (Element) panels.item(i);
        String name = panelNode.getElementsByTagName("name").item(0).getTextContent();
        NodeList aliasNodes = panelNode.getElementsByTagName("alias");
        String[] aliases = new String[aliasNodes.getLength()];
        for (int n = 0, countA = aliasNodes.getLength(); n < countA; n++) {
          aliases[n] = aliasNodes.item(n).getTextContent();
        }
        String[][] gateTree = null;
        Map<String, List<String>> specials = null;

        NodeList gateTreeNodes = panelNode.getElementsByTagName("gateTree");
        if (gateTreeNodes.getLength() > 0) {
          if (gateTreeNodes.getLength() > 1) {
            // TODO error
          }
          List<String[]> treeNodeList = new ArrayList<>();
          Element gateTreeNode = (Element) gateTreeNodes.item(0);
          NodeList treeNodes = gateTreeNode.getElementsByTagName("node");
          for (int iT = 0, countT = treeNodes.getLength(); iT < countT; iT++) {
            Element node = (Element) treeNodes.item(iT);
            String n = node.getAttribute("name");
            String p = node.getAttribute("parent");
            String[] value = p.isEmpty() ? new String[] {n} : new String[] {n, p};
            treeNodeList.add(value);
          }
          gateTree = treeNodeList.toArray(new String[treeNodeList.size()][]);
        }

        NodeList specialNodes = panelNode.getElementsByTagName("specials");
        if (specialNodes.getLength() > 0) {
          if (specialNodes.getLength() > 1) {
            // TODO error
          }
          specials = new HashMap<>();
          Element specialsNode = (Element) specialNodes.item(0);
          NodeList specialNodeList = specialsNode.getElementsByTagName("special");
          for (int iS = 0, countS = specialNodeList.getLength(); iS < countS; iS++) {
            Element specialNode = (Element) specialNodeList.item(iS);
            String key = specialNode.getElementsByTagName("key").item(0).getTextContent();
            List<String> values = new ArrayList<>();
            NodeList valueNodes = specialNode.getElementsByTagName("value");
            for (int iV = 0, countV = valueNodes.getLength(); iV < countV; iV++) {
              values.add(valueNodes.item(iV).getTextContent());
            }
            specials.put(key, values);
          }
        }

        panelsFound.add(new Panel(name, gateTree, specials, aliases));
      }
    } catch (ParserConfigurationException e) {
      e.printStackTrace();
    } catch (SAXException | IOException e) {
      e.printStackTrace();
    }
    return panelsFound;
  }

}
