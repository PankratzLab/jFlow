package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;

public class MosaicArmsStep extends Step {

  private final Requirement<Integer> numThreadsReq;

  public static MosaicArmsStep create(Project proj, final Step parseSamplesStep,
                                      Requirement<Integer> numThreadsReq) {
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(numThreadsReq);
    return new MosaicArmsStep(numThreadsReq, reqSet);
  }

  public static final String NAME = "Create Mosaic Arms File";
  public static final String DESC = "";

  private MosaicArmsStep(Requirement<Integer> numThreadsReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED));
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj, Variables variables) {

    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
  }

  @Override
  public void run(Project proj, Variables variables) {
    Mosaicism.findOutliers(proj);
  }

  @Override
  public String getCommandLine(Project proj, Variables variables) {
    String kvCmd = "";

    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
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
  public boolean checkIfOutputExists(Project proj, Variables variables) {
    return Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
  }
}
