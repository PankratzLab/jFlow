package org.genvisis.one.ben.fcs.auto.proc;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import javax.swing.SwingConstants;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.genvisis.one.ben.flowannot.GateTree;

public class VisualizationProcessor implements SampleProcessor {

  final String autoDir;
  final String ovvrDir;
  final String ovvrSfx;
  final String ovvrMatch;
  final String outDir;

  public VisualizationProcessor(String a, String o, String m, String mS, String mM) {
    this.autoDir = a;
    this.outDir = o;
    this.ovvrDir = m;
    this.ovvrSfx = mS;
    this.ovvrMatch = mM;
  }

  private static String[][] hardcodedAddlImages = {{"effector helper Tcells (CCR7/CD45RA)",
                                                    "effector helper Tcells (CCR7- CD45RA+)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"naive helper Tcells (CCR7/CD45RA)",
                                                    "naive helper Tcells (CCR7+ CD45RA+)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"effector memory helper Tcells (CCR7/CD45RA)",
                                                    "effector memory helper Tcells (CCR7- CD45RA-)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"central memory helper Tcells (CCR7/CD45RA)",
                                                    "central memory helper Tcells (CCR7+ CD45RA-)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"effector cytotoxic Tcells (CCR7/CD45RA)",
                                                    "effector cytotoxic Tcells  (CCR7-  CD45RA+)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"naive cytotoxic Tcells (CCR7/CD45RA)",
                                                    "naive cytotoxic Tcells (CCR7+ , CD45RA+)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"effector memory cytotoxic Tcells (CCR7/CD45RA)",
                                                    "effector memory cytotoxic Tcells (CCR7- , CD45RA-)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"central memory cytotoxic Tcells (CCR7/CD45RA)",
                                                    "central memory cytotoxic Tcells (CCR7+ , CD45RA-)",
                                                    "Comp-BV 421-A (CCR7)",
                                                    "Comp-BV 711-A (CD45RA)"},
                                                   {"Effector_Cytotoxic Tcells (CD28/CD27)",
                                                    "effector cytotoxic Tcells (CCR7/CD45RA)",
                                                    "Comp-BV 510-A (CD28)", "Comp-BB 515-A (CD27)"},
                                                   {"Effector_Memory_Cytotoxic Tcells (CD28/CD27)",
                                                    "effector memory cytotoxic Tcells (CCR7/CD45RA)",
                                                    "Comp-BV 510-A (CD28)",
                                                    "Comp-BB 515-A (CD27)"},};

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

  static final Map<String, List<AddlImage>> addlImgs = new HashMap<>();
  static final Map<String, String> dimSwitch = new HashMap<>();
  {
    for (String[] addlImg : hardcodedAddlImages) {
      String parent;
      parent = addlImg[1];
      addlImgs.put(parent, new ArrayList<AddlImage>());
      addlImgs.get(parent).add(new AddlImage(addlImg[2], addlImg[3], parent, addlImg[0]));
    }
    dimSwitch.put("Comp-BV 605-A (CD95)", "Comp-BV605-A (CD95)");
    dimSwitch.put("Comp-BV 510-A (CD28)", "Comp-BV510-A (CD28)");
    dimSwitch.put("Comp-BB 515-A (CD27)", "Comp-BB515-A (CD27)");
    dimSwitch.put("Comp-BB515-A (CD27)", "Comp-FITC-A (CD27)");
    dimSwitch.put("Comp-BV 421-A (CCR7)", "Comp-BV421-A (CCR7)");
    dimSwitch.put("Comp-BV 711-A (CD45RA)", "Comp-BV711-A (CD45RA)");

  }

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    System.gc();
    final FCSPlot fcp = new FCSPlot();
    fcp.getPanel().getColorScheme()[2] = new Color(128, 128, 128, 64);
    fcp.getPanel().allowSkip = false;

    long time1 = System.nanoTime();
    fcp.loadWorkspaceFile(sn.wspFile);
    Set<String> sampeIds = fcp.getWorkbench().getAllSamples();
    String id = sampeIds.toArray(new String[1])[0];
    if (fcp.getWorkbench().getSample(id).fcsFile.equals(sn.fcsFile)) {
      fcp.setCurrentSampleInWSP(sn.fcsFile);
    } else if (fcp.getWorkbench().getSample(id).fcsFile.equals(ext.rootOf(sn.fcsFile, true))) {
      fcp.setCurrentSampleInWSP(ext.rootOf(sn.fcsFile, true));
    }
    fcp.setCurrentSampleInWSP(id);

    //    boolean hasAll = true;
    //
    //    for (String s : fcp.getGatingStrategy().getAllGateNames()) {
    //      Gate g = fcp.getGatingStrategy().gateMap.get(s);
    //
    //      String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
    //      String outFile = outDir + cleanedName + "/" + cleanedName + "."
    //                       + ext.replaceWithLinuxSafeCharacters(g.getName());
    //
    //      if (Files.exists(outFile + ".png")) {
    //        if (addlImgs.containsKey(g.getName())) {
    //          boolean all = true;
    //          for (AddlImage addl : addlImgs.get(g.getName())) {
    //            String outFile2 = outDir + cleanedName + "/" + cleanedName + "."
    //                              + ext.replaceWithLinuxSafeCharacters(addl.name);
    //            if (!Files.exists(outFile2 + ".png")) {
    //              all = false;
    //              break;
    //            }
    //          }
    //          if (all) {
    //            continue;
    //          }
    //        } else {
    //          continue;
    //        }
    //      }
    //
    //      hasAll = false;
    //      break;
    //    }
    //    if (hasAll) {
    //      log.report("All screenshots found for " + fcp.getGatingStrategy().getAllGateNames().size()
    //                 + " gates; Skipping FCS file: " + sn.fcsFile);
    //      return;
    //    }

    long time2 = System.nanoTime();
    fcp.loadFile(sn.fcsFile, true);
    FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
    while (loader == null || loader.getLoadState() != LOAD_STATE.LOADED) {
      Thread.yield();
    }
    if (ovvrDir != null) {
      loader.loadGateOverrides(ovvrDir + ext.removeDirectoryInfo(sn.fcsFile) + ovvrSfx, ovvrMatch);
    }
    int rowCnt = loader.getCount();

    long time3 = System.nanoTime();

    // String fNum = fcp.discoverFNumFile(autoDir);
    // if (fNum == null)
    // return;
    // fcp.loadAutoValues(fNum);

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
      fcp.gateSelected(g.getParentGate(), false);
      times[0] = System.nanoTime();

      if (ext.indexOfStr(g.getName(), GateTree.HELPER_SUB_PARENTS) >= 0) {
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.HELPER_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.CYTOTOXIC_SUB_PARENTS) >= 0) {
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.CYTOTOXIC_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.HELPER_SUB_IMGS) >= 0) {
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.HELPER_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.CYTOTOXIC_SUB_IMGS) >= 0) {
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.CYTOTOXIC_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.EFFECTOR_SUB_PARENTS) >= 0) {
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.EFFECTOR_MEM_SUB_PARENTS) >= 0) {
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_MEM_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.EFFECTOR_SUB_IMGS) >= 0) {
        fcp.gateSelected(g, false);
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_SUB_IMGS);

      } else if (ext.indexOfStr(g.getName(), GateTree.EFFECTOR_MEM_SUB_IMGS) >= 0) {
        fcp.gateSelected(g, false);
        fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
        g.setFillGate(false);
        fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_MEM_SUB_IMGS);

      } else {
        fcp.setPlotType(PLOT_TYPE.HEATMAP);
        g.setFillGate(true);
      }

      g.setColor(1);
      g.setDisplayName(false);
      if (g.getXDimension() == null && g.getYDimension() != null) {
        // correct for swapped histogram
        g.setXDimension(g.getYDimension());
        g.setYDimension(null);
      }
      fcp.setXDataName(g.getXDimension().getParam());
      if (g.getYDimension() != null) {
        fcp.setYDataName(g.getYDimension().getParam());
      } else {
        fcp.setPlotType(PLOT_TYPE.HISTOGRAM);
      }
      fcp.getPanel().setTitle(g.getName());
      fcp.getPanel().setTitleLocation(SwingConstants.NORTH);
      fcp.getPanel().setDisplayTitle(true);
      fcp.getPanel().setTitleFontSize(20f);

      //       fcp.setClassifierGate(g.getName());

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
    // fcp = null;
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

  private void createAddlImgs(final FCSPlot fcp, FCSDataLoader loader, Gate g, String gateName,
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

        if (ext.indexOfStr(addl.name, GateTree.HELPER_SUB_PARENTS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.HELPER_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.CYTOTOXIC_SUB_PARENTS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.CYTOTOXIC_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.HELPER_SUB_IMGS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.HELPER_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.CYTOTOXIC_SUB_IMGS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.CYTOTOXIC_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.EFFECTOR_SUB_PARENTS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.EFFECTOR_MEM_SUB_PARENTS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_MEM_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.EFFECTOR_SUB_IMGS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_SUB_IMGS);

        } else if (ext.indexOfStr(addl.name, GateTree.EFFECTOR_MEM_SUB_IMGS) >= 0) {
          fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
          g.setFillGate(false);
          fcp.loadOverridesAsClusterColors(loader, GateTree.EFFECTOR_MEM_SUB_IMGS);

        } else {
          fcp.setPlotType(PLOT_TYPE.HEATMAP);
          g.setFillGate(true);
        }

        String x = addl.xDim;
        while (!loader.getAllDisplayableNames(DATA_SET.ALL).contains(x)) {
          x = dimSwitch.get(x);
        }
        fcp.setXDataName(x);

        String y = addl.yDim;
        while (!loader.getAllDisplayableNames(DATA_SET.ALL).contains(y)) {
          y = dimSwitch.get(y);
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
