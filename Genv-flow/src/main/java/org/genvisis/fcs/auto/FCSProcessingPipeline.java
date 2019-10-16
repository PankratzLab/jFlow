package org.genvisis.fcs.auto;

import java.io.IOException;

import org.genvisis.fcs.auto.proc.InclusionProcessor;
import org.genvisis.fcs.auto.proc.PercentageAndCountWriterFactory;
import org.genvisis.fcs.auto.proc.ProcessorFactory;
import org.genvisis.fcs.auto.proc.SampleProcessor;
import org.genvisis.fcs.auto.proc.VisualizationProcessor;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.ext;

public class FCSProcessingPipeline {

  public FCSProcessingPipeline(String fcs, String wsp, String auto, String out, String highP,
                               String lowP, String ovvrDir, String ovvrSuff, String ovvrMatch,
                               String clusterDir, String clusterSuffix, String addlImgs,
                               String clustersFile) {
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
    this.addlImgsFile = addlImgs;
    this.clusterOverrideFile = clustersFile;
  }

  private String fcsDir, wspDir, autoDir, outDir, highPrioFile, lowPrioFile;
  private String ovvrDir, ovvrSuff, ovvrMatch;
  private String clustDir, clustSfx;
  private String addlImgsFile, clusterOverrideFile;

  private void run(PIPELINE pipeToRun, int panel) throws IOException {

    ProcessorFactory<? extends SampleProcessor> pf = null;

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
                                              clustDir, clustSfx, addlImgsFile,
                                              clusterOverrideFile);
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
    String addlImgs = null;
    String clustersFiles = null;
    String clusterSuffix = "_subFirst_TRUE_normalize_FALSE.IntMatrix.txt.gz";
    int panel = -1;
    PIPELINE pipe = PIPELINE.VIZ;

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
      } else if (arg.startsWith("addlImgs=")) {
        addlImgs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("clusterOverrides=")) {
        clustersFiles = arg.split("=")[1];
        numArgs--;
      }
    }

    new FCSProcessingPipeline(fcs, wsp, auto, out, highPriorityFile, lowPriorityFile,
                              gateOverrideDir, gateOverrideFileSuffix, gateOverrideMatchFile,
                              clusterDir, clusterSuffix, addlImgs, clustersFiles).run(pipe, panel);
  }
}
