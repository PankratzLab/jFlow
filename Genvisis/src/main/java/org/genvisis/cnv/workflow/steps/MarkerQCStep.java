package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.Collection;
import java.util.EnumSet;
import java.util.Set;
import java.util.StringJoiner;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.ListSelectionRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

public class MarkerQCStep extends Step {

  public static final String NAME = "Run Marker QC Metrics";
  public static final String DESC = "";

  public static MarkerQCStep create(Project proj, Step parseSamplesStep,
                                    Requirement<Integer> numThreadsReq) {
    String[] tgtMkrFiles = proj.TARGET_MARKERS_FILENAMES.getValue();
    final Requirement<File> targetMarkersReq = new Requirement.FileRequirement("targetMarkers",
                                                                               "A targetMarkers files listing the markers to QC.",
                                                                               tgtMkrFiles != null && tgtMkrFiles.length >= 1 ? new File(tgtMkrFiles[0])
                                                                                                                              : new File(""));
    final Set<String> sampleDataHeaders;
    if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue()) && proj.getSampleData(false) != null) {
      sampleDataHeaders = proj.getSampleData(false).getMetaHeaders();
    } else {
      sampleDataHeaders = Sets.newHashSet();
    }
    final Set<String> defaultBatchHeaders = Sets.intersection(sampleDataHeaders,
                                                              MarkerMetrics.DEFAULT_SAMPLE_DATA_BATCH_HEADERS);
    final ListSelectionRequirement batchHeadersReq = new ListSelectionRequirement("sampleDataColumns",
                                                                                  "SampleData column headers to use as batches for batch effects calculations",
                                                                                  sampleDataHeaders,
                                                                                  defaultBatchHeaders,
                                                                                  true);
    final Requirement<Boolean> exportAllReq = new Requirement.OptionalBoolRequirement("exportAll",
                                                                                      "Export all markers in project.",
                                                                                      true);
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(exportAllReq)
                                                                                 .add(targetMarkersReq))
                                                       .add(batchHeadersReq).add(numThreadsReq);
    return new MarkerQCStep(proj, reqSet, exportAllReq, targetMarkersReq, batchHeadersReq,
                            numThreadsReq);
  }

  final Project proj;
  Requirement<Boolean> exportAllReq;
  Requirement<File> targetMarkersReq;
  Requirement<Integer> numThreadsReq;
  ListSelectionRequirement batchHeadersReq;

  private MarkerQCStep(Project proj, RequirementSet reqSet, Requirement<Boolean> exportAllReq,
                       Requirement<File> targetMarkersReq, ListSelectionRequirement batchHeadersReq,
                       Requirement<Integer> numThreadsReq) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED));
    this.proj = proj;
    this.exportAllReq = exportAllReq;
    this.targetMarkersReq = targetMarkersReq;
    this.batchHeadersReq = batchHeadersReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Nothing to do here
  }

  @Override
  public void run(Variables variables) {
    boolean allMarkers = variables.get(exportAllReq);
    String tgtFile = allMarkers ? null : variables.get(targetMarkersReq).getAbsolutePath();
    boolean[] samplesToExclude = proj.getSamplesToExclude();
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    Set<String> batchHeaders = ImmutableSet.copyOf(variables.get(batchHeadersReq));
    MarkerMetrics.fullQC(proj, samplesToExclude, tgtFile, true, batchHeaders, numThreads);
    MarkerMetrics.filterMetrics(proj);
  }

  @Override
  public String getCommandLine(Variables variables) {
    boolean allMarkers = variables.get(exportAllReq);
    File tgtFile = variables.get(targetMarkersReq);
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    Collection<String> batchHeaders = variables.get(batchHeadersReq);
    StringJoiner args = new StringJoiner(" ");
    args.add(Files.getRunString());
    args.add(MarkerMetrics.class.getCanonicalName());
    args.add("-fullQC");
    args.add("proj=" + proj.getPropertyFilename());
    if (!allMarkers && tgtFile != null) args.add("markers=" + tgtFile.getAbsolutePath());
    String batchHeadersArg = String.join(",", batchHeaders);
    args.add(MarkerMetrics.BATCH_HEADERS_ARG + "=" + batchHeadersArg);
    args.add(PSF.Ext.NUM_THREADS_COMMAND + numThreads);
    return args.toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
    return Files.exists(markerMetricsFile);
  }

}
