package org.genvisis.one.ben.fcs.auto;

import java.io.IOException;
import org.genvisis.one.ben.fcs.auto.proc.InclusionProcessor;
import org.genvisis.one.ben.fcs.auto.proc.PercentageAndCountWriterFactory;
import org.genvisis.one.ben.fcs.auto.proc.ProcessorFactory;
import org.genvisis.one.ben.fcs.auto.proc.SampleProcessor;
import org.genvisis.one.ben.fcs.auto.proc.VisualizationProcessor;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.ext;

public class FCSProcessingPipeline {

  public static final String[][] MATCH = {{"lymph", "Lymphocytes (SSC-A v FSC-A)"},
                                          {"Singlets", "Single Cells (FSC-H v FSC-W)"},
                                          {"PE.A", "Live cells (PE-)"},
                                          {"CD3.", "Tcells (CD3+ CD19-)"},};

  public FCSProcessingPipeline(String fcs, String wsp, String auto, String out, String highP,
                               String lowP, String ovvrDir, String ovvrSuff, String ovvrMatch,
                               String clusterDir, String clusterSuffix) {
    fcsDir = fcs;
    wspDir = wsp;
    autoDir = auto;
    outDir = out;
    highPrioFile = highP;
    lowPrioFile = lowP;
    this.ovvrDir = ovvrDir;
    this.ovvrSuff = ovvrSuff;
    this.ovvrMatch = ovvrMatch;
    this.clustDir = clusterDir;
    this.clustSfx = clusterSuffix;
  }

  private String fcsDir, wspDir, autoDir, outDir, highPrioFile, lowPrioFile;
  private String ovvrDir, ovvrSuff, ovvrMatch;
  private String clustDir, clustSfx;

  private void run(PIPELINE pipeToRun, int panel) throws IOException {

    ProcessorFactory<? extends SampleProcessor> pf = null;

    // pf = new GateAssignmentFactory();
    // pf = new PercentageWriterFactory();
    // pf = new LeafDataSamplerFactory("/scratch.global/cole0482/FCS/", new Logger());
    // pf = new ProcessorFactory<ConcordanceProcessor>() {
    //
    // ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap = new
    // ConcurrentHashMap<>();
    // ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMapTree = new
    // ConcurrentHashMap<>();
    // public void cleanup(Object owner) {
    // PrintWriter writer = Files.getAppropriateWriter(OUT + "concordance.xln");
    //
    // Set<String> files = resultMap.keySet();
    // for (String s : files) {
    // for (int i = 0; i < MATCH.length; i++) {
    // writer.println(s + "\t" + MATCH[i][1] + "\t" + resultMap.get(s).get(MATCH[i][1]));
    // }
    // }
    // writer.flush();
    // writer.close();
    //
    // writer = Files.getAppropriateWriter(OUT + "concordance_tree.xln");
    //
    // files = resultMapTree.keySet();
    // for (String s : files) {
    // for (int i = 0; i < MATCH.length; i++) {
    // writer.println(s + "\t" + MATCH[i][1] + "\t" + resultMapTree.get(s).get(MATCH[i][1]));
    // }
    // }
    // writer.flush();
    // writer.close();
    // };
    //
    // @Override
    // public ConcordanceProcessor createProcessor(Object owner, int index) {
    // return new ConcordanceProcessor(resultMap, resultMapTree);
    // }
    // };
    //
    // SamplingPipeline sp = new SamplingPipeline(1, null, WSP, FCS, null, OUT, pf);
    switch (pipeToRun) {
      case BOOL:
        pf = new ProcessorFactory<SampleProcessor>() {

          @Override
          public void cleanup(Object owner) {}

          @Override
          public SampleProcessor createProcessor(Object owner, int index) {
            return new InclusionProcessor(outDir, ovvrDir, ovvrSuff, ovvrMatch);
          }
        };
        break;
      case VIZ:
        pf = new ProcessorFactory<SampleProcessor>() {

          @Override
          public void cleanup(Object owner) {}

          @Override
          public SampleProcessor createProcessor(Object owner, int index) {
            return new VisualizationProcessor(autoDir, outDir, ovvrDir, ovvrSuff, ovvrMatch,
                                              clustDir, clustSfx);
          }
        };
        break;
      case PCTS_CNTS:
        pf = new PercentageAndCountWriterFactory(outDir, panel, ovvrDir, ovvrSuff, ovvrMatch);
        break;
      default:
        System.err.println("Error - pipeline option " + pipeToRun
                           + " not recognized!  Options are: "
                           + ArrayUtils.toStr(PIPELINE.values(), ", "));
        return;
    }

    SamplingPipeline sp = new SamplingPipeline(1, null, wspDir, fcsDir, null, outDir, panel,
                                               new String[] {highPrioFile, lowPrioFile}, pf);

    sp.run();

    while (!sp.checkFinished()) {
      Thread.yield();
    }

  }

  private static enum PIPELINE {
    BOOL, VIZ, PCTS_CNTS;
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;

    // private static final String FCS = "/home/pankrat2/shared/flow/fcs2/";
    // private static final String FCS = "/scratch.global/cole0482/FCS/src/";
    // private static final String FCS =
    // "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS CTL/";
    // private static final String WSP = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
    // private static final String OUT = "/scratch.global/cole0482/FCS/testConcordance/";

    // String fcs = "/home/pankrat2/shared/flow/fcs2/";
    // String wsp = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
    // String auto = "/home/pankrat2/shared/flow/testAutoGate/test2/gates2/";
    // String out = "/scratch.global/cole0482/FCS/testConcordance/";
    String fcs = "F:/Flow/boolGating/";
    String wsp = fcs;
    String auto = null;
    String out = fcs;
    String highPriorityFile = null;
    String lowPriorityFile = null;
    String gateOverrideDir = null;
    String gateOverrideMatchFile = null;
    String gateOverrideFileSuffix = ".boolMatrix.txt.gz";
    String clusterDir = null;
    String clusterSuffix = "_subFirst_TRUE_normalize_FALSE.IntMatrix.txt.gz";
    int panel = -1;
    PIPELINE pipe = PIPELINE.VIZ;
    //    boolean test = Files.isWindows();
    //    if (test) {
    //      String dir = "C:\\Users\\Ben\\Desktop\\dev\\work\\flow\\phenograph\\";
    //      fcs = dir + "src\\";
    //      wsp = dir + "src\\";
    //      out = dir + "out\\";
    //      gateOverrideDir = dir + "src\\";
    //      gateOverrideMatchFile = dir + "src\\ovvrMatch.txt";
    //      new FCSProcessingPipeline(fcs, wsp, auto, out, highPriorityFile, lowPriorityFile,
    //                                gateOverrideDir, gateOverrideFileSuffix, gateOverrideMatchFile).run(
    //                                                                                                    PIPELINE.VIZ,
    //                                                                                                    -1);
    //      return;
    //    }

    for (String arg : args) {
      if (arg.startsWith("fcs=")) {
        fcs = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("wsp=")) {
        wsp = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("auto=")) {
        auto = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("out=")) {
        out = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("panel=")) {
        panel = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("priority=")) {
        highPriorityFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("lowPriority=")) {
        lowPriorityFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("pipe=")) {
        pipe = PIPELINE.valueOf(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("gateOverride=")) {
        gateOverrideDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("gateOverrideMatchup=")) {
        gateOverrideMatchFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("clusterDir=")) {
        clusterDir = arg.split("=")[1];
        numArgs--;
      }
    }

    new FCSProcessingPipeline(fcs, wsp, auto, out, highPriorityFile, lowPriorityFile,
                              gateOverrideDir, gateOverrideFileSuffix, gateOverrideMatchFile,
                              clusterDir, clusterSuffix).run(pipe, panel);
  }
}
