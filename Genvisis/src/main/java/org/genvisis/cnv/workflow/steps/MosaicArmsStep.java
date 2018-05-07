package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.common.Files;

public class MosaicArmsStep extends Step {

  private final Requirement numThreadsReq;

  public static MosaicArmsStep create(Project proj, final Step parseSamplesStep,
                                      Requirement numThreadsReq, double priority) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(numThreadsReq);
    return new MosaicArmsStep(numThreadsReq, reqSet, priority);
  }

  public static final String NAME = "Create Mosaic Arms File";
  public static final String DESC = "";

  private MosaicArmsStep(Requirement numThreadsReq, RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED), priority);
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {

    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    Mosaicism.findOutliers(proj);
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String kvCmd = "";

    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    if (numThreads != proj.NUM_THREADS.getValue()) {
      kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
    }

    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    if (kvCmd.length() > 0) {
      cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile)
         .append(kvCmd).append("\n");
    }
    return cmd.append(Files.getRunString())
              .append(" cnv.analysis.Mosaicism proj=" + proj.getPropertyFilename()).toString();
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    return Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
  }
}
