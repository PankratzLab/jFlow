package org.genvisis.fcs.auto.proc;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.concurrent.ConcurrentHashMap;

import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;

public class ConcordanceProcessor extends AbstractLoadingSampleProcessor {

  private static final String AUTO_FILE_SUFF = "def.txt.gz";

  private String autoDir;
  private String[] filesInAutoDir;
  String[][] autoData;
  private ConcurrentHashMap<String, ConcurrentHashMap<String, String>> results;
  private ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultsTree;

  public static final String[][] MATCH = {{"lymph", "Lymphocytes (SSC-A v FSC-A)"},
                                          {"Singlets", "Single Cells (FSC-H v FSC-W)"},
                                          {"PE.A", "Live cells (PE-)"},
                                          {"CD3.", "Tcells (CD3+ CD19-)"},};

  public ConcordanceProcessor(String autoDir,
                              ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap,
                              ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap2) {
    super();
    this.autoDir = autoDir;
    filesInAutoDir = (new File(autoDir)).list(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(AUTO_FILE_SUFF);
      }
    });
    this.results = resultMap;
    this.resultsTree = resultMap2;
  }

  @Override
  public void processSample(SampleNode sn, Logger log) throws IOException {
    loadPopsAndGates(sn);
    loadData(sn);

    loadAutoValues(sn);
    if (autoData == null) {
      System.err.println("Error - problem occured with auto-data.  Check log and try again for sample: "
                         + sn.fcsFile);
      return;
    }

    boolean[][] hand = new boolean[autoData[0].length][d.getCount()];
    String[] gateNames = new String[hand.length];
    for (int i = 0; i < hand.length; i++) {
      Gate gate = sn.gating.gateMap.get(ConcordanceProcessor.MATCH[i][1]);
      hand[i] = gate.gate(d);
      gateNames[i] = gate.getName();
    }

    boolean[][] auto = new boolean[autoData[0].length][];
    for (int i = 0; i < auto.length; i++) {
      String[] data = Matrix.extractColumn(autoData, i);
      auto[i] = new boolean[data.length];
      for (int k = 0; k < data.length; k++) {
        auto[i][k] = Boolean.valueOf(data[k]);
      }
    }

    ConcurrentHashMap<String, String> fileResults = new ConcurrentHashMap<>();
    ConcurrentHashMap<String, String> treeResults = new ConcurrentHashMap<>();
    for (int i = 0; i < auto.length; i++) {
      boolean[] prevHand = i == 0 ? ArrayUtils.booleanArray(hand[i].length, true) : hand[i - 1];
      boolean[] prevAuto = i == 0 ? ArrayUtils.booleanArray(auto[i].length, true) : auto[i - 1];
      boolean[] handV = hand[i];
      boolean[] autoV = auto[i];
      int[] conc = concordance(handV, autoV);
      int[] prevConc = prevConcordance(prevHand, handV, prevAuto, autoV);
      fileResults.put(gateNames[i], ArrayUtils.toStr(conc, "\t"));
      treeResults.put(gateNames[i], ArrayUtils.toStr(prevConc, "\t"));
    }
    results.put(sn.fcsFile, fileResults);
    resultsTree.put(sn.fcsFile, treeResults);

    System.out.println("Done with " + sn.fcsFile);

    d = null;
    System.gc();
  }

  private int[] prevConcordance(boolean[] prevH, boolean[] h, boolean[] prevA, boolean[] a) {
    int tp = 0, tn = 0, fp = 0, fn = 0;
    for (int i = 0; i < h.length; i++) {
      if (prevH[i] && prevA[i]) {
        if (h[i] && a[i]) {
          tp++;
        } else if (h[i] && !a[i]) {
          fn++;
        } else if (!h[i] && a[i]) {
          fp++;
        } else if (!h[i] && !a[i]) {
          tn++;
        }
      }
    }
    return new int[] {tp, tn, fp, fn};
  }

  /**
   * @param hand
   * @param auto
   * @return TP, TN, FP, FN
   */
  private int[] concordance(boolean[] hand, boolean[] auto) {
    int tp = 0, tn = 0, fp = 0, fn = 0;
    for (int i = 0; i < hand.length; i++) {
      if (hand[i] && auto[i]) {
        tp++;
      } else if (hand[i] && !auto[i]) {
        fn++;
      } else if (!hand[i] && auto[i]) {
        fp++;
      } else if (!hand[i] && !auto[i]) {
        tn++;
      }
    }
    return new int[] {tp, tn, fp, fn};
  }

  private void loadAutoValues(SampleNode sn) {
    String[] pts = sn.fcsFile.split("_");
    String fNum = null;
    for (String p : pts) {
      if (p.startsWith("F")) {
        fNum = p;
        break;
      }
    }
    if (fNum == null) {
      System.err.println("Error - couldn't parse F-number for fcs file: " + sn.fcsFile);
      return;
    }
    String file = null;
    for (String s : filesInAutoDir) {
      pts = s.split("_");
      for (String p : pts) {
        if (fNum.equals(p)) {
          file = s;
          break;
        }
      }
    }
    if (file == null) {
      System.err.println("Error - Couldn't match F-number of fcs file with auto-gated files: "
                         + sn.fcsFile + " | " + fNum);
      return;
    }
    autoData = HashVec.loadFileToStringMatrix(this.autoDir + file, true, null);
  }

}
