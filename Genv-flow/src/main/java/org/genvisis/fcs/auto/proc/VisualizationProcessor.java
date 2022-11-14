package org.genvisis.fcs.auto.proc;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import javax.swing.SwingConstants;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.fcs.FCSDataLoader;
import org.genvisis.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.fcs.FCSDataLoader.LOAD_STATE;
import org.genvisis.fcs.JFlow;
import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

public class VisualizationProcessor extends AbstractSampleProcessor {

  final String ovvrDir;
  final String ovvrSfx;
  final String ovvrMatch;
  final String outDir;
  final String clustDir;
  final String clustSfx;

  private final Multimap<String, AddlImage> addlImgs = HashMultimap.create();
  private final Map<String, ClusterOverride> clusterOverrides = new HashMap<>();

  public VisualizationProcessor(String outDir, String overrideDir, String overrideSuffix,
                                String overrideMatch, String clusterDir, String clusterSuffix,
                                String addlImgsFile, String clusterOverrideFile,
                                String dimOverrideFile) {
    super();
    this.outDir = outDir;
    this.ovvrDir = overrideDir;
    this.ovvrSfx = overrideSuffix;
    this.ovvrMatch = overrideMatch;
    this.clustDir = clusterDir;
    this.clustSfx = clusterSuffix;
    if (dimOverrideFile != null) {
      loadDimOverrides(dimOverrideFile);
    }
    if (addlImgsFile != null) {
      loadAddlImages(addlImgsFile);
    }
    if (clusterOverrideFile != null) {
      loadClusterOverrides(clusterOverrideFile);
    }
  }

  private void loadAddlImages(String addlImgsFile) {
    if (!Files.exists(addlImgsFile)) {
      System.err.println("Error - Couldn't find specified additional images file: " + addlImgsFile);
      return;
    }
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    try {
      DocumentBuilder builder = factory.newDocumentBuilder();
      Document doc = builder.parse(new File(addlImgsFile));
      doc.getDocumentElement().normalize();

      NodeList images = doc.getElementsByTagName("image");
      for (int i = 0, count = images.getLength(); i < count; i++) {
        Element imageNode = (Element) images.item(i);
        String parent = imageNode.getElementsByTagName("parent").item(0).getFirstChild()
                                 .getTextContent().trim();
        String name = imageNode.getElementsByTagName("name").item(0).getFirstChild()
                               .getTextContent().trim();
        String x = imageNode.getElementsByTagName("x").item(0).getFirstChild().getTextContent()
                            .trim();
        String y = imageNode.getElementsByTagName("y").item(0).getFirstChild().getTextContent()
                            .trim();
        addlImgs.put(parent, new AddlImage(x, y, parent, name));
      }

    } catch (ParserConfigurationException e) {
      e.printStackTrace();
    } catch (SAXException | IOException e) {
      e.printStackTrace();
    }
  }

  private void loadClusterOverrides(String file) {
    if (!Files.exists(file)) {
      System.err.println("Error - Couldn't find specified cluster override file: " + file);
      return;
    }
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    try {
      DocumentBuilder builder = factory.newDocumentBuilder();
      Document doc = builder.parse(new File(file));
      doc.getDocumentElement().normalize();

      NodeList clusters = doc.getElementsByTagName("cluster");
      for (int i = 0, count = clusters.getLength(); i < count; i++) {
        Element clusterNode = (Element) clusters.item(i);
        // load keys
        NodeList keys = clusterNode.getElementsByTagName("key");
        List<String> keyList = new ArrayList<>();
        for (int k = 0, countK = keys.getLength(); k < countK; k++) {
          String key = ((Element) keys.item(k)).getTextContent();
          keyList.add(key);
        }
        // load color keys
        NodeList colorKeys = clusterNode.getElementsByTagName("colorKey");
        String[] colorKey = new String[colorKeys.getLength()];
        for (int k = 0, countK = colorKeys.getLength(); k < countK; k++) {
          colorKey[k] = ((Element) colorKeys.item(k)).getTextContent();
        }
        // load which gate to use
        String type = clusterNode.getElementsByTagName("gate").item(0).getTextContent();
        boolean parent;
        switch (type.toLowerCase()) {
          case "parent":
            parent = true;
            break;
          case "self":
            parent = false;
            break;
          default:
            parent = true;
            break;
        }
        for (String k : keyList) {
          clusterOverrides.put(k, new ClusterOverride(colorKey, parent));
        }
      }

    } catch (ParserConfigurationException e) {
      e.printStackTrace();
    } catch (SAXException | IOException e) {
      e.printStackTrace();
    }
  }

  static class ClusterOverride {
    String[] overrideNames;
    boolean gateParent;

    public ClusterOverride(String[] clust, boolean gateParent) {
      this.overrideNames = clust;
      this.gateParent = gateParent;
    }
  }

  static class AddlImage {

    String xDim;
    String yDim;
    String parentName;
    String name;

    public AddlImage(String x, String y, String p, String n) {
      this.xDim = x;
      this.yDim = y;
      this.parentName = p;
      this.name = n;
    }
  }

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    System.gc();
    final JFlow fcp = new JFlow();
    fcp.getPanel().getColorScheme()[2] = new Color(128, 128, 128, 64);
    fcp.getPanel().allowSkip = false;

    long time1 = System.nanoTime();
    fcp.loadWorkspaceFile(sn.wspFile);
    Set<String> sampeIds = fcp.getWorkbench().getAllSamples();
    String id = sampeIds.toArray(new String[1])[0];
    final String fcsFile = fcp.getWorkbench().getSample(id).fcsFile;
    if (fcsFile.equals(sn.fcsFile)) {
      fcp.setCurrentSampleInWSP(sn.fcsFile);
    } else if (fcsFile.equals(ext.removeDirectoryInfo(sn.fcsFile))) {
      fcp.setCurrentSampleInWSP(ext.removeDirectoryInfo(sn.fcsFile));
    } else {
      fcp.setCurrentSampleInWSP(id);
    }

    long time2 = System.nanoTime();
    fcp.loadFile(sn.fcsFile, true);
    FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
    while (loader == null || loader.getLoadState() != LOAD_STATE.LOADED) {
      Thread.yield();
    }
    if (ovvrDir != null) {
      String kmeans = ovvrDir + ext.removeDirectoryInfo(sn.fcsFile) + ovvrSfx;
      if (!Files.exists(kmeans)) {
        log.reportError("K-Means boolean file missing!  Couldn't find " + kmeans);
      } else {
        loader.loadGateOverrides(kmeans, ovvrMatch);
      }
    }
    int rowCnt = loader.getCount();

    if (clustDir != null) {
      String clust = clustDir + ext.removeDirectoryInfo(sn.fcsFile) + clustSfx;
      if (!Files.exists(clust)) {
        log.reportError("Phenograph cluster file missing!  Couldn't find " + clust);
      } else {
        fcp.loadClusterAssignments(clust);
      }
    }

    long time3 = System.nanoTime();

    fcp.setSize(1000, 800);
    fcp.getPanel().setSize(800, 600);
    fcp.setSDVisible(false, false);
    fcp.setSDVisible(false, true);
    fcp.setMedianVisible(false, false);
    fcp.setMedianVisible(false, true);
    fcp.getPanel().setLayersInBase(new byte[] {0, 1, 99});
    fcp.setPlotType(PLOT_TYPE.DOT_PLOT);

    long time4 = System.nanoTime();
    ArrayList<String> gateNames = new ArrayList<>();
    ArrayList<long[]> gateTimes = new ArrayList<>();

    for (String s : fcp.getGatingStrategy().getAllGateNames()) {
      fcp.getPanel().forceSkip = true;
      long[] times = new long[4];

      Gate g = fcp.getGatingStrategy().gateMap.get(s);

      String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
      String outFile = outDir + cleanedName + "/" + cleanedName + "."
                       + ext.replaceWithLinuxSafeCharacters(g.getName());

      if (Files.exists(outFile + ".png")) {
        if (addlImgs.containsKey(g.getName())) {
          boolean all = true;
          for (AddlImage addl : addlImgs.get(g.getName())) {
            String outFile2 = outDir + cleanedName + "/" + cleanedName + "."
                              + ext.replaceWithLinuxSafeCharacters(addl.name);
            if (!Files.exists(outFile2 + ".png")) {
              all = false;
              break;
            }
          }
          if (all) {
            continue;
          }
        } else {
          continue;
        }
      }

      gateNames.add(s);

      fcp.clearClusterAssigns();
      if (clusterOverrides.containsKey(g.getName())) {
        ClusterOverride co = clusterOverrides.get(g.getName());
        if (co.gateParent) {
          fcp.gateSelected(g.getParentGate(), false);
        } else {
          fcp.gateSelected(g, false);
        }
        times[0] = System.nanoTime();

        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, co.overrideNames);
      } else {
        fcp.gateSelected(g.getParentGate(), false);
        times[0] = System.nanoTime();
        if (fcp.getClusterAssignments() != null) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
        } else {
          fcp.setPlotType(PLOT_TYPE.HEATMAP);
          g.setFillGate(true);
        }
      }

      g.setColor(1);
      g.setDisplayName(false);
      if (g.getXDimension() == null && g.getYDimension() != null) {
        // correct for swapped histogram
        g.setXDimension(g.getYDimension());
        g.setYDimension(null);
      }

      String x = g.getXDimension().getParam();
      while (!loader.getAllDisplayableNames(DATA_SET.ALL).contains(x)) {
        x = replaceName(x);
      }
      fcp.setXDataName(x);

      if (g.getYDimension() != null) {
        String y = g.getXDimension().getParam();
        while (!loader.getAllDisplayableNames(DATA_SET.ALL).contains(y)) {
          y = replaceName(y);
        }
        fcp.setYDataName(y);
      } else {
        fcp.setPlotType(PLOT_TYPE.HISTOGRAM);
      }
      fcp.getPanel().setTitle(g.getName());
      fcp.getPanel().setTitleLocation(SwingConstants.NORTH);
      fcp.getPanel().setDisplayTitle(true);
      fcp.getPanel().setTitleFontSize(20f);

      // fcp.setClassifierGate(g.getName());

      fcp.getPanel().classifierPrev = false;
      times[1] = System.nanoTime();

      fcp.getPanel().forceSkip = false;
      fcp.getPanel().setForceGatesChanged();
      fcp.getPanel().createImage();
      Thread.yield();
      times[2] = System.nanoTime();

      fcp.screencap(outFile + ".png");
      fcp.clearClusterAssigns();

      times[3] = System.nanoTime();

      g.setFillGate(false);

      gateTimes.add(times);

      createAddlImgs(fcp, loader, g, g.getName(), cleanedName);
    }

    loader.emptyAndReset();
    System.gc();

    long time5 = System.nanoTime();

    StringBuilder sb1 = new StringBuilder("TIMING-HDR").append("\t");
    sb1.append("FCS_FILE").append("\t");
    sb1.append("ROWS").append("\t");
    sb1.append("INIT").append("\t");
    sb1.append("WSP_LOAD").append("\t");
    sb1.append("LOAD").append("\t");
    sb1.append("CONF").append("\t");
    for (String s : gateNames) {
      sb1.append(s).append("\t");
      sb1.append("CONF").append("\t");
      sb1.append("CREATE").append("\t");
      sb1.append("SCREEN").append("\t");
    }
    sb1.append("CLEANUP");

    StringBuilder sb = new StringBuilder("TIMING").append("\t");
    sb.append(sn.fcsFile).append("\t");
    sb.append(rowCnt).append("\t");
    sb.append(TimeUnit.NANOSECONDS.toSeconds(time1)).append("\t");
    sb.append(TimeUnit.NANOSECONDS.toSeconds(time2)).append("\t");
    sb.append(TimeUnit.NANOSECONDS.toSeconds(time3)).append("\t");
    sb.append(TimeUnit.NANOSECONDS.toSeconds(time4)).append("\t");
    for (long[] timing : gateTimes) {
      for (long l : timing) {
        sb.append(TimeUnit.NANOSECONDS.toSeconds(l)).append("\t");
      }
    }
    sb.append(TimeUnit.NANOSECONDS.toSeconds(time5));

    System.out.println(sb1.toString());
    System.out.println(sb.toString());
  }

  private void createAddlImgs(final JFlow fcp, FCSDataLoader loader, Gate g, String gateName,
                              String cleanedName) {
    if (addlImgs.containsKey(gateName)) {
      fcp.getPanel().forceSkip = true;
      fcp.gateSelected(g, false);

      for (AddlImage addl : addlImgs.get(gateName)) {
        String outFile2 = outDir + cleanedName + "/" + cleanedName + "."
                          + ext.replaceWithLinuxSafeCharacters(addl.name);
        if (Files.exists(outFile2 + ".png")) {
          fcp.getPanel().forceSkip = false;
          continue;
        }
        fcp.getPanel().setTitle(addl.name);

        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        if (clusterOverrides.containsKey(g.getName())) {
          fcp.loadOverridesAsClusterColors(loader, clusterOverrides.get(g.getName()).overrideNames);
        } else {
          fcp.setPlotType(PLOT_TYPE.HEATMAP);
          g.setFillGate(true);
        }

        String x = addl.xDim;
        while (!loader.getAllDisplayableNames(DATA_SET.ALL).contains(x)) {
          x = replaceName(x);
        }
        fcp.setXDataName(x);

        String y = addl.yDim;
        while (!loader.getAllDisplayableNames(DATA_SET.ALL).contains(y)) {
          y = replaceName(y);
        }
        fcp.setYDataName(y);

        fcp.getPanel().forceSkip = false;
        fcp.getPanel().setForceGatesChanged();
        fcp.getPanel().createImage();
        Thread.yield();
        fcp.screencap(outFile2 + ".png");

        if (addlImgs.containsKey(addl.name)) {
          createAddlImgs(fcp, loader, g, addl.name, cleanedName);
        }
      }
    }
  }
}
