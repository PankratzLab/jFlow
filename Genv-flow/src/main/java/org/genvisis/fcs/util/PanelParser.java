package org.genvisis.fcs.util;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.genvisis.fcs.auto.WSPLoader;
import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.Gating;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import com.google.common.base.Predicates;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

public class PanelParser {

  private Multimap<String, Integer> panels(String wspFile)
      throws SAXException, IOException, ParserConfigurationException {
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    Document doc = builder.parse(new File(wspFile));
    doc.getDocumentElement().normalize();

    NodeList nList = doc.getElementsByTagName("GroupNode");

    Multimap<String, Integer> panelNames = HashMultimap.create();

    for (int i = 0; i < nList.getLength(); i++) {
      Node n = nList.item(i);
      if (n.getNodeType() == Node.ELEMENT_NODE) {
        Element e = (Element) n;
        String name = e.getAttributes().getNamedItem("name").getNodeValue();

        e = ((Element) e.getElementsByTagName("Group").item(0));
        e = ((Element) e.getElementsByTagName("SampleRefs").item(0));
        if (e == null) {
          continue;
        }
        NodeList list = e.getElementsByTagName("SampleRef");
        int refs = list.getLength();
        for (int s = 0; s < refs; s++) {
          Element e1 = (Element) list.item(s);
          Integer id = Integer.parseInt(e1.getAttribute("sampleID"));
          panelNames.put(name, id);
        }
      }
    }

    Collection<Integer> c = panelNames.get("All Samples");
    Set<Integer> removeFromAll =
        panelNames
            .entries()
            .stream()
            .filter(
                Predicates.and(
                    e -> c.contains(e.getValue()), e -> !e.getKey().equals("All Samples")))
            .map(Entry::getValue)
            .collect(Collectors.toSet());

    for (Integer e : removeFromAll) {
      panelNames.remove("All Samples", e);
    }

    return panelNames;
  }

  private void parse(String wspFile, String outputFile)
      throws ParserConfigurationException, SAXException, IOException {

    Multimap<String, Integer> panels = panels(wspFile);

    WSPLoader loader = new WSPLoader();
    loader.loadWorkspaceGating(wspFile);

    List<SampleNode> allSamples = loader.getAllSamples();

    Document dom;
    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

    try {
      DocumentBuilder db = dbf.newDocumentBuilder();
      dom = db.newDocument();
      dom.setXmlStandalone(true);

      Element rootEle = dom.createElement("panels");
      dom.appendChild(rootEle);

      Element e = null;
      Element e1 = null;
      for (String panelStr : panels.keySet()) {
        e = dom.createElement("panel");
        rootEle.appendChild(e);

        e1 = dom.createElement("name");
        e1.appendChild(dom.createTextNode(panelStr));
        e.appendChild(e1);

        e1 = dom.createElement("alias");
        e1.appendChild(dom.createTextNode(panelStr));
        e.appendChild(e1);

        e1 = dom.createElement("gateTree");
        e.appendChild(e1);
        e = e1;

        Set<GateTreeNode> tree = new LinkedHashSet<>();
        for (SampleNode sn : allSamples) {
          if (!panels.get(panelStr).contains(Integer.parseInt(sn.id))) {
            continue;
          }
          Gating g = sn.getGating();
          List<Gate> roots = g.getRootGates();
          for (Gate root : roots) {
            tree.add(new GateTreeNode(root.getName(), null));
            recurse(tree, root);
          }
        }

        for (GateTreeNode node : tree) {
          e1 = dom.createElement("node");
          e1.setAttribute("name", node.gate);
          if (node.parent != null) {
            e1.setAttribute("parent", node.parent);
          }
          e.appendChild(e1);
        }
      }

      Transformer tr = TransformerFactory.newInstance().newTransformer();
      tr.setOutputProperty(OutputKeys.INDENT, "yes");
      tr.setOutputProperty(OutputKeys.METHOD, "xml");
      tr.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
      tr.setOutputProperty(OutputKeys.DOCTYPE_PUBLIC, "yes");
      tr.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

      tr.transform(new DOMSource(dom), new StreamResult(new FileOutputStream(outputFile)));
    } catch (Exception e) {
      e.printStackTrace();
    }

    System.out.println("Done!");
  }

  static class GateTreeNode {
    String gate;
    String parent;

    public GateTreeNode(String g, String p) {
      this.gate = g;
      this.parent = p;
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((gate == null) ? 0 : gate.hashCode());
      result = prime * result + ((parent == null) ? 0 : parent.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      GateTreeNode other = (GateTreeNode) obj;
      if (gate == null) {
        if (other.gate != null) return false;
      } else if (!gate.equals(other.gate)) return false;
      if (parent == null) {
        if (other.parent != null) return false;
      } else if (!parent.equals(other.parent)) return false;
      return true;
    }
  }

  private void recurse(Set<GateTreeNode> tree, Gate gate) {
    List<Gate> children = gate.getChildGates();
    for (Gate g : children) {
      tree.add(new GateTreeNode(g.getName(), gate.getName()));
    }
    for (Gate g : children) {
      recurse(tree, g);
    }
  }

  public static void main(String[] args)
      throws ParserConfigurationException, SAXException, IOException {
    CLI cli = new CLI(PanelParser.class);

    cli.addArg("wsp", "Example Workspace (.wsp) file", true);
    cli.addArg(
        "out",
        "Output file name (path optional; if not present, the output file will be created in the same directory as the input workspace file)",
        false);

    cli.parseWithExit(args);

    String wsp = cli.get("wsp");
    String out =
        cli.has("out")
            ? cli.get("out")
            : (Files.getNextAvailableFilename(
                ext.parseDirectoryOfFile(wsp) + "custom_panel_#.xml", "#"));
    new PanelParser().parse(wsp, out);
  }
}
