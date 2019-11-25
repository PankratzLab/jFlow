package org.genvisis.fcs.auto;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.genvisis.fcs.auto.proc.InclusionProcessor;
import org.genvisis.fcs.auto.proc.PercentageAndCountWriterFactory;
import org.genvisis.fcs.auto.proc.ProcessorFactory;
import org.genvisis.fcs.auto.proc.SampleProcessor;
import org.genvisis.fcs.auto.proc.VisualizationProcessor;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;

public class FCSProcessingPipeline {

  public FCSProcessingPipeline(String fcs, String wsp, String out, String highP, String lowP,
                               String ovvrDir, String ovvrSuff, String ovvrMatch, String clusterDir,
                               String clusterSuffix, String addlImgs, String clustersFile,
                               String dimensionOverrideFile) {
    fcsDir = fcs;
    wspDir = wsp;
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
    this.dimensionOverrideFile = dimensionOverrideFile;
  }

  private String fcsDir, wspDir, outDir, highPrioFile, lowPrioFile;
  private String ovvrDir, ovvrSuff, ovvrMatch;
  private String clustDir, clustSfx;
  private String addlImgsFile, clusterOverrideFile, dimensionOverrideFile;

  private void run(PIPELINE pipeToRun, List<Panel> panels) throws IOException {

    ProcessorFactory<? extends SampleProcessor> pf = null;

    switch (pipeToRun) {
      case BOOL:
        pf = new ProcessorFactory<SampleProcessor>() {

          @Override
          public void cleanup(Object owner) {}

          @Override
          public SampleProcessor createProcessor(Object owner, int index) {
            return new InclusionProcessor(outDir, ovvrDir, ovvrSuff, ovvrMatch,
                                          dimensionOverrideFile);
          }
        };
        break;
      case VIZ:
        pf = new ProcessorFactory<SampleProcessor>() {

          @Override
          public void cleanup(Object owner) {}

          @Override
          public SampleProcessor createProcessor(Object owner, int index) {
            return new VisualizationProcessor(outDir, ovvrDir, ovvrSuff, ovvrMatch, clustDir,
                                              clustSfx, addlImgsFile, clusterOverrideFile,
                                              dimensionOverrideFile);
          }
        };
        break;
      case PCTS_CNTS:
        pf = new PercentageAndCountWriterFactory(outDir, ovvrDir, ovvrSuff, ovvrMatch,
                                                 dimensionOverrideFile);
        break;
      default:
        System.err.println("Error - pipeline option " + pipeToRun
                           + " not recognized!  Options are: "
                           + ArrayUtils.toStr(PIPELINE.values(), ", "));
        return;
    }

    SamplingPipeline sp = new SamplingPipeline(1, null, wspDir, fcsDir, null, outDir, panels,
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
    CLI cli = new CLI(FCSProcessingPipeline.class);

    cli.addArg("fcs", "Directory with FCS files (finds all FCS files, including in subdirectories)",
               true);
    cli.addArg("wsp", "Directory with WSP files (finds all WSP files, including in subdirectories)",
               true);
    cli.addArg("out", "Output directory", true);
    cli.addArg("pipe",
               "Pipeline to run, one of " + Arrays.stream(PIPELINE.values()).map(PIPELINE::name)
                                                  .collect(Collectors.joining(", ")),
               true);

    cli.addArg("panelDefs", "XML file with panel definitions", true);

    cli.addArg("addlImgs", "XML file with additional image definitions", false);
    cli.addArg("dimensionOverrides", "XML file with fluorochrome name overrides", false);

    cli.addArg("clusterDir",
               "Directory with cluster assignment files - single-column files with gate name on the first line and a line for each data point following, with an integer on each line indicating which cluster each data point belongs to.",
               false);
    cli.addArg("clusterSfx", "Suffix of cluster assignment files", false);
    cli.addArg("clusterOverrides",
               "XML file defining custom overrides/assignments for existing gate cluster assignments",
               false);

    cli.addArg("gateOverrideDir",
               "Directory with files (one per FCS file) defining custom gate inclusions. Files are tab-delimited, with custom names on the first line and a line for each data point following with TRUE/FALSE to indicate if that datapoint is included in the gating.",
               false);
    cli.addArg("gateOverrideSfx", "Suffix of gate override files.", false);
    cli.addArg("gateOverrideMatchup",
               "File defining which existing gates should be overridden by the custom gates found in the gate override files.",
               false);

    cli.addArg("priority",
               "File containing fcs files (one per line, with spaces and non-path-safe characters replaced with underscores) that will be processed first",
               false);
    cli.addArg("lowPriority",
               "File containing fcs files (one per line, with spaces and non-path-safe characters replaced with underscores) that will be processed last",
               false);

    cli.parseWithExit(args);

    String panelDefFile = cli.get("panelDefs");
    String fcs = cli.get("fcs");
    String wsp = cli.get("wsp");
    String out = cli.get("out");
    PIPELINE pipe = PIPELINE.valueOf(cli.get("pipe"));

    String highPriorityFile = cli.has("priority") ? cli.get("priority") : null;
    String lowPriorityFile = cli.has("lowPriority") ? cli.get("lowPriority") : null;

    String addlImgs = cli.has("addlImgs") ? cli.get("addlImgs") : null;
    String clustersFile = cli.has("clusterOverrides") ? cli.get("clusterOverrides") : null;
    String dimensionOverrideFile = cli.has("dimensionOverrides") ? cli.get("dimensionOverrides")
                                                                 : null;
    String gateOverrideDir = cli.has("gateOverrideDir") ? cli.get("gateOverrideDir") : null;
    String gateOverrideSfx = cli.has("gateOverrideSfx") ? cli.get("gateOverrideSfx") : null;
    String gateOverrideMatch = cli.has("gateOverrideMatchup") ? cli.get("gateOverrideMatchup")
                                                              : null;

    String clusterDir = cli.has("clusterDir") ? cli.get("clusterDir") : null;
    String clusterSfx = cli.has("clusterSfx") ? cli.get("clusterSfx") : null;

    List<Panel> panels = WSPLoader.loadPanelsFromFile(new FileInputStream(panelDefFile));

    new FCSProcessingPipeline(fcs, wsp, out, highPriorityFile, lowPriorityFile, gateOverrideDir,
                              gateOverrideSfx, gateOverrideMatch, clusterDir, clusterSfx, addlImgs,
                              clustersFile, dimensionOverrideFile).run(pipe, panels);
  }

}
