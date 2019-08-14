package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Set;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gwas.utils.FurtherAnalysisQc;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.OptionalFileRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.utils.gwas.Qc;

public class IdentifyProblemMarkersStep extends Step {

  public static final String NAME = "Identify Problematic Markers";
  public static final String DESC = "";

  public static IdentifyProblemMarkersStep create(Project proj, final MarkerQCStep markerQCStep,
                                                  final FurtherAnalysisQCStep faqcStep) {
    final Requirement<Step> markerQCStepReq = new Requirement.StepRequirement(markerQCStep);
    final Requirement<Step> faqcStepReq = new Requirement.StepRequirement(faqcStep);
    final OptionalFileRequirement optionalDropFile = new OptionalFileRequirement("optionalDrops",
                                                                                 "A single-column file with no header of additional markers to exclude as problematic.",
                                                                                 new File(""));

    return new IdentifyProblemMarkersStep(proj, markerQCStepReq, faqcStepReq, optionalDropFile);
  }

  final Project proj;
  final Requirement<File> optionalDropsReq;

  public IdentifyProblemMarkersStep(Project proj, Requirement<Step> markerQCStepReq,
                                    Requirement<Step> faqcStepReq, Requirement<File> optDropReq) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(markerQCStepReq).add(faqcStepReq).add(optDropReq),
          EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.optionalDropsReq = optDropReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // no op
  }

  @Override
  public void run(Variables variables) {
    // multi-column, with header
    String file1 = proj.RESULTS_DIRECTORY.getValue(false, true)
                   + MarkerMetrics.MARKERS_TO_EXCLUDE_FILENAME;
    // single column, no header
    String file2 = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                   + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR + FurtherAnalysisQc.MARKER_QC_DROPS;
    // expecting single column, no header
    File file3 = variables.get(optionalDropsReq);

    Set<String> mkrs1 = HashVec.loadFileToHashSet(file1, new int[] {0}, "", true);
    Set<String> mkrs2 = HashVec.loadFileToHashSet(file2, false);
    Set<String> mkrs3 = file3.getPath().equals("") ? new HashSet<>()
                                                   : HashVec.loadFileToHashSet(file2, false);

    Set<String> combined = new HashSet<>(mkrs3);
    combined.addAll(mkrs2);
    combined.addAll(mkrs1);

    Files.writeIterable(combined, proj.FILTERED_MARKERS_FILENAME.getValue());
  }

  @Override
  public String getCommandLine(Variables variables) {
    return getStepCommandLine(proj, variables);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.FILTERED_MARKERS_FILENAME.getValue());
  }

  public static void main(String[] args) {
    Project proj = Step.parseProject(args);
    StepBuilder sb = new StepBuilder(proj);
    Step samplesStep = sb.generateSamplesParsingStep();
    MarkerQCStep markerQCStep = sb.generateMarkerQCStep(samplesStep);
    PlinkExportStep plinkStep = sb.generatePlinkExportStep(samplesStep);
    GwasQCStep qcStep = sb.generateGwasQCStep(plinkStep);
    AncestryStep ancStep = sb.generateAncestryStep(qcStep);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(plinkStep, qcStep, ancStep);
    IdentifyProblemMarkersStep step = IdentifyProblemMarkersStep.create(proj, markerQCStep,
                                                                        faqcStep);
    Variables variables = step.parseArguments(args);
    Step.run(proj, step, variables);
  }

}
