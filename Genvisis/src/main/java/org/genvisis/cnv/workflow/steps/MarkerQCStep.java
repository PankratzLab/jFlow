package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.List;
import java.util.Map;
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
import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

public class MarkerQCStep extends Step {

  public static final String NAME = "Run Marker QC Metrics";
  public static final String DESC = "";

  private static final Requirement exportAllReq = new Requirement.OptionalBoolRequirement("Export all markers in project.",
                                                                                          true);

  public static MarkerQCStep create(Project proj, Step parseSamplesStep, Requirement numThreadsReq,
                                    double priority) {
    String[] tgtMkrFiles = proj.TARGET_MARKERS_FILENAMES.getValue();
    final Requirement targetMarkersReq = new Requirement.FileRequirement("A targetMarkers files listing the markers to QC.",
                                                                         tgtMkrFiles != null && tgtMkrFiles.length >= 1 ? tgtMkrFiles[0]
                                                                                                                        : "");
    final Set<String> sampleDataHeaders;
    if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue()) && proj.getSampleData(false) != null) {
      sampleDataHeaders = proj.getSampleData(false).getMetaHeaders();
    } else {
      sampleDataHeaders = Sets.newHashSet();
    }
    final Set<String> defaultBatchHeaders = Sets.intersection(sampleDataHeaders,
                                                              MarkerMetrics.DEFAULT_SAMPLE_DATA_BATCH_HEADERS);
    final ListSelectionRequirement batchHeadersReq = new ListSelectionRequirement("SampleData column headers to use as batches for batch effects calculations",
                                                                                  sampleDataHeaders,
                                                                                  defaultBatchHeaders,
                                                                                  true);

    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(exportAllReq)
                                                                                 .add(targetMarkersReq))
                                                       .add(batchHeadersReq).add(numThreadsReq);
    return new MarkerQCStep(reqSet, targetMarkersReq, batchHeadersReq, numThreadsReq, priority);
  }

  Requirement targetMarkersReq;
  Requirement numThreadsReq;
  ListSelectionRequirement batchHeadersReq;

  private MarkerQCStep(RequirementSet reqSet, Requirement targetMarkersReq,
                       ListSelectionRequirement batchHeadersReq, Requirement numThreadsReq,
                       double priority) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED), priority);
    this.targetMarkersReq = targetMarkersReq;
    this.batchHeadersReq = batchHeadersReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // Nothing to do here
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
    String tgtFile = allMarkers ? null : variables.get(this).get(targetMarkersReq);
    boolean[] samplesToExclude = proj.getSamplesToExclude();
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    Set<String> batchHeaders = ImmutableSet.copyOf(Requirement.ListSelectionRequirement.parseArgValString(variables.get(this)
                                                                                                                   .get(batchHeadersReq)));
    MarkerMetrics.fullQC(proj, samplesToExclude, tgtFile, true, batchHeaders, numThreads);
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
    String tgtFile = variables.get(this).get(targetMarkersReq);
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    List<String> batchHeaders = Requirement.ListSelectionRequirement.parseArgValString(variables.get(this)
                                                                                                .get(batchHeadersReq));
    StringJoiner args = new StringJoiner(" ");
    args.add(Files.getRunString());
    args.add(MarkerMetrics.class.getCanonicalName());
    args.add("-fullQC");
    args.add("proj=" + proj.getPropertyFilename());
    if (!allMarkers) args.add("markers=" + tgtFile);
    String batchHeadersArg = String.join(",", batchHeaders);
    args.add(MarkerMetrics.BATCH_HEADERS_ARG + "=" + batchHeadersArg);
    args.add(PSF.Ext.NUM_THREADS_COMMAND + numThreads);
    return args.toString();
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
    return Files.exists(markerMetricsFile);
  }

}
